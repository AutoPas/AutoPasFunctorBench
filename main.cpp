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

/**
 * Mini benchmark tool to estimate the inner most kernel performance of AutoPas
 * @return
 */
int main() {

    // type aliases for ease of use
    using Particle = autopas::MoleculeLJ<double>;
    using Cell = autopas::FullParticleCell<autopas::MoleculeLJ<double>>;

    // some constants that define the benchmark
    constexpr double cutoff{3.}; // is also the cell size
    constexpr bool shift{false};
    constexpr bool mixing{false};
    constexpr autopas::FunctorN3Modes functorN3Modes{autopas::FunctorN3Modes::Both};
    constexpr bool newton3{true};
    constexpr bool globals{false};
    // choose functor based on available architecture
#ifdef __AVX__
    autopas::LJFunctorAVX<Particle, shift, mixing, functorN3Modes, globals> functor{cutoff};
#elif __ARM_FEATURE_SVE
    autopas::LJFunctorSVE<Particle, shift, mixing, functorN3Modes, globals> functor{cutoff};
#endif
    constexpr double epsilon24{24.};
    constexpr double sigmaSquare{1.};
    functor.setParticleProperties(epsilon24, sigmaSquare);

    // define scenario
    constexpr size_t numCells{2};
    std::array<Cell, numCells> cells;
    constexpr std::array<size_t, numCells> numParticlesPerCell{100, 50};

    // initialize cells with randomly distributed particles
    timer.at("Initialization").start();
    for (size_t cellId = 0; cellId < numCells; ++cellId) {
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
                    particleId + ((std::numeric_limits<size_t>::max() / numCells) * cellId),
                    0};
            cells[cellId].addParticle(p);
        }
        functor.SoALoader(cells[cellId], cells[cellId]._particleSoABuffer, 0);
    }
    timer.at("Initialization").stop();

    // actual benchmark
    timer.at("Functor").start();
    // TODO offer option to also test FunctorSingle
    functor.SoAFunctorPair(cells[0]._particleSoABuffer, cells[1]._particleSoABuffer, newton3);
    timer.at("Functor").stop();

    // print particles to CSV for checking and prevent compiler from optimizing everything away.
    timer.at("Output").start();
    std::ofstream csvFile("particles.csv");
    if (not csvFile.is_open()) {
        throw std::runtime_error("FILE NOT OPEN!");
    }
    csvFile << "CellId,ParticleId,rX,rY,rZ,fX,fY,fZ\n";
    for (size_t cellId = 0; cellId < numCells; ++cellId) {
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

    // count interactions
    timer.at("InteractionCounter").start();
    int calcsDist{0};
    int calcsForce{0};
    const auto cutoffSquared{cutoff * cutoff};
    for (const auto p0: cells[0]) {
        for (const auto p1: cells[1]) {
            ++calcsDist;
            if (distSquared(p0.getR(), p1.getR()) <= cutoffSquared) {
                ++calcsForce;
            }
        }
    }
    timer.at("InteractionCounter").stop();

    // print timer and statistics

    const auto gflops =
            static_cast<double>(calcsDist * 8 + calcsForce * functor.getNumFlopsPerKernelCall()) * 10e-9;

    using autopas::utils::ArrayUtils::operator<<;

    std::cout << "Particels per cell : " << numParticlesPerCell << "\n"
              << "Hit rate           : " << (static_cast<double>(calcsForce) / calcsDist) << "\n"
              << "GFLOPs             : " << gflops << "\n"
              << "GFLOPs/sec         : " << (gflops / (timer.at("Functor").getTotalTime() * 10e-9)) << "\n";

    printTimer();
}
