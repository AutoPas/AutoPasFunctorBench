#include <iostream>

// clang-format off
#ifdef __AVX__
#include <autopas/molecularDynamics/LJFunctorAVX.h>
#elif __ARM_FEATURE_SVE
#include <autopas/molecularDynamics/LJFunctorSVE.h>
#else
#error "Platform supports neither AVX nor ARM!"
#endif
// clang-format on

#include <autopas/molecularDynamics/MoleculeLJ.h>
#include <autopas/cells/FullParticleCell.h>
#include <autopas/utils/Timer.h>
#include <fstream>

// type aliases for ease of use
using Particle = autopas::MoleculeLJ<double>;
using Cell = autopas::FullParticleCell<autopas::MoleculeLJ<double>>;
// some constants that define the benchmark
constexpr bool shift{false};
constexpr bool mixing{false};
constexpr autopas::FunctorN3Modes functorN3Modes{autopas::FunctorN3Modes::Both};
constexpr bool newton3{true};
constexpr bool globals{false};

#ifdef __AVX__
using Functor = autopas::LJFunctorAVX<Particle, shift, mixing, functorN3Modes, globals>;
#elif __ARM_FEATURE_SVE
using Functor = autopas::LJFunctorSVE<Particle, shift, mixing, functorN3Modes, globals> ;
#endif

void checkFunctorType(const Functor &fun) {
    int identificationHits = 0;
#ifdef __AVX__
    if (dynamic_cast<const autopas::LJFunctorAVX<Particle, shift, mixing, functorN3Modes, globals> *>(&fun)) {
        std::cout << "Using AVX Functor" << std::endl;
        ++identificationHits;
    }
#endif
#ifdef __ARM_FEATURE_SVE
    if (dynamic_cast<const autopas::LJFunctorSVE<Particle, shift, mixing, functorN3Modes, globals> *>(&fun)) {
        std::cout << "Using SVE Functor" << std::endl;
        ++identificationHits;
    }
#endif
    if (identificationHits != 1) {
        throw std::runtime_error(
                "checkFunctorType matched "
                + std::to_string(identificationHits)
                + " types! There should only be one match.");
    }
}

double distSquared(std::array<double, 3> a, std::array<double, 3> b) {
    using autopas::utils::ArrayMath::sub;
    using autopas::utils::ArrayMath::dot;
    const auto c = sub(a, b);   // 3 FLOPS
    return dot(c, c);           // 3+2=5 FLOPs
}

std::map<std::string, autopas::utils::Timer> timer{
        {"Initialization",     autopas::utils::Timer()},
        {"Functor",            autopas::utils::Timer()},
        {"Output",             autopas::utils::Timer()},
        {"InteractionCounter", autopas::utils::Timer()},
};

void printTimer() {
    for (const auto &[name, t]: timer) {
        std::cout
                << std::setw(18)
                << std::left
                << name
                << " : "
                << std::setprecision(3)
                << std::setw(8)
                << static_cast<double>(timer[name].getTotalTime()) * 10e-9
                << " [s]\n";
    }
}

void initialization(Functor &functor, std::vector<Cell> &cells, const std::vector<size_t> &numParticlesPerCell,
                    double cutoff) {
    // initialize cells with randomly distributed particles
    timer.at("Initialization").start();
    for (size_t cellId = 0; cellId < numParticlesPerCell.size(); ++cellId) {
        for (size_t particleId = 0; particleId < numParticlesPerCell[cellId]; ++particleId) {
            Particle p{
                    {
                            // particles are next to each other in X direction
                            rand() / static_cast<double>(RAND_MAX) * cutoff + cutoff * cellId,
                            rand() / static_cast<double>(RAND_MAX) * cutoff,
                            rand() / static_cast<double>(RAND_MAX) * cutoff,
                    },
                    {0., 0., 0.,},
                    // every cell gets its own id space
                    particleId + ((std::numeric_limits<size_t>::max() / numParticlesPerCell.size()) * cellId),
                    0};
            cells[cellId].addParticle(p);
        }
        functor.SoALoader(cells[cellId], cells[cellId]._particleSoABuffer, 0);
    }
    timer.at("Initialization").stop();
}

void applyFunctor(Functor &functor, std::vector<Cell> &cells) {
    timer.at("Functor").start();
    functor.SoAFunctorPair(cells[0]._particleSoABuffer, cells[1]._particleSoABuffer, newton3);
    timer.at("Functor").stop();
}

void csvOutput(Functor &functor, std::vector<Cell> &cells) {
    timer.at("Output").start();
    std::ofstream csvFile("particles.csv");
    if (not csvFile.is_open()) {
        throw std::runtime_error("FILE NOT OPEN!");
    }
    csvFile << "CellId,ParticleId,rX,rY,rZ,fX,fY,fZ\n";
    for (size_t cellId = 0; cellId < cells.size(); ++cellId) {
        functor.SoAExtractor(cells[cellId], cells[cellId]._particleSoABuffer, 0);
        for (size_t particleId = 0; particleId < cells[cellId].numParticles(); ++particleId) {
            const auto &p = cells[cellId][particleId];
            using autopas::utils::ArrayUtils::to_string;
            csvFile << cellId << ","
                    << p.getID() << ","
                    << to_string(p.getR(), ",", {"", ""}) << ","
                    << to_string(p.getF(), ",", {"", ""})
                    << "\n";
        }
    }
    csvFile.close();
    timer.at("Output").stop();
}

std::tuple<size_t, size_t> countInteractions(std::vector<Cell> &cells, double cutoff) {
    timer.at("InteractionCounter").start();
    size_t calcsDist{0};
    size_t calcsForce{0};
    const auto cutoffSquared{cutoff * cutoff};
    for (const auto &p0: cells[0]) {
        for (const auto &p1: cells[1]) {
            ++calcsDist;
            if (distSquared(p0.getR(), p1.getR()) <= cutoffSquared) {
                ++calcsForce;
            }
        }
    }
    timer.at("InteractionCounter").stop();
    return {calcsDist, calcsForce};
}

/**
 * Mini benchmark tool to estimate the inner most kernel performance of AutoPas
 * @return
 */
int main() {

    constexpr double cutoff{3.}; // is also the cell size

    // choose functor based on available architecture
    Functor functor{cutoff};
    checkFunctorType(functor);

    constexpr double epsilon24{24.};
    constexpr double sigmaSquare{1.};
    functor.setParticleProperties(epsilon24, sigmaSquare);

    // define scenario
    const std::vector<size_t> numParticlesPerCell{1000, 1000};
    constexpr size_t iterations{100};
    size_t calcsDistTotal{0};
    size_t calcsForceTotal{0};
    // repeat the whole experiment multiple times and average results
    for (size_t iteration = 0; iteration < iterations; ++iteration) {
        std::vector<Cell> cells{2};

        initialization(functor, cells, numParticlesPerCell, cutoff);

        // TODO offer option to also test FunctorSingle
        // actual benchmark
        applyFunctor(functor, cells);

        // print particles to CSV for checking and prevent compiler from optimizing everything away.
        csvOutput(functor, cells);

        // gather data for analysis
        const auto [calcsDist, calcsForce] = countInteractions(cells, cutoff);
        calcsDistTotal += calcsDist;
        calcsForceTotal += calcsForce;
    }

    // print timer and statistics
    const auto gflops =
            static_cast<double>(calcsDistTotal * 8 + calcsForceTotal * functor.getNumFlopsPerKernelCall()) * 10e-9;

    using autopas::utils::ArrayUtils::operator<<;

    std::cout
            << "Iterations         : " << iterations << "\n"
            << "Particels per cell : " << numParticlesPerCell << "\n"
            << "Avgerage hit rate  : " << (static_cast<double>(calcsForceTotal) / calcsDistTotal) << "\n"
            << "GFLOPs             : " << gflops << "\n"
            << "GFLOPs/sec         : " << (gflops / (timer.at("Functor").getTotalTime() * 10e-9)) << "\n";

    printTimer();
}
