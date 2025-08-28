#include <iostream>

#include <molecularDynamicsLibrary/interpolationKernels/KryptonKernel.h>
#include <molecularDynamicsLibrary/PairwiseInterpolantFunctor.h>
#include <molecularDynamicsLibrary/AbInitioKryptonPairFunctor.h>
#include <molecularDynamicsLibrary/MoleculeLJ.h>

#include <autopas/cells/FullParticleCell.h>
#include <autopas/utils/Timer.h>
#include <fstream>

// type aliases for ease of use
using Particle = mdLib::MoleculeLJ;
using Cell = autopas::FullParticleCell<mdLib::MoleculeLJ>;

// some constants that define the benchmark
constexpr autopas::FunctorN3Modes functorN3Modes{autopas::FunctorN3Modes::Both};
constexpr bool newton3{true};
constexpr bool globals{false};

using InterpolantFunctor = mdLib::PairwiseInterpolantFunctor<mdLib::KryptonKernel, Particle, functorN3Modes, globals>;
using ReferenceFunctor = mdLib::AbInitioKryptonPairFunctor<Particle, functorN3Modes, globals>;

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
                << static_cast<double>(timer[name].getTotalTime()) * 1e-9
                << " [s]\n";
    }
}

void initialization(std::vector<Cell> &cells, const std::vector<size_t> &numParticlesPerCell,
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
    }
    timer.at("Initialization").stop();
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

    constexpr double cutoff{1.5}; // is also the cell size

    // choose functor based on available architecture

    mdLib::KryptonKernel kernel = mdLib::KryptonKernel{};

    const double a = 0.35;
    const double b = cutoff;

    std::vector<size_t> nodes {16};
    std::vector<double> splits {};

    InterpolantFunctor functor {kernel, cutoff, a, nodes, splits};
    ReferenceFunctor refereceFunctor {cutoff};

    // define scenario
    const std::vector<size_t> numParticlesPerCell{2000, 2000};
    constexpr size_t iterations{1000};
    size_t calcsDistTotal{0};
    size_t calcsForceTotal{0};
    // repeat the whole experiment multiple times and average results
    for (size_t iteration = 0; iteration < iterations; ++iteration) {
        std::vector<Cell> cells{2};

        initialization(cells, numParticlesPerCell, cutoff);

        timer.at("Functor").start();
        for (auto& p1 : cells[0]) {
            for (auto& p2 : cells[1]) {
                functor.AoSFunctor(p1, p2, newton3);
            }
        }
        timer.at("Functor").stop();
        
        // gather data for analysis
        const auto [calcsDist, calcsForce] = countInteractions(cells, cutoff);
        calcsDistTotal += calcsDist;
        calcsForceTotal += calcsForce;
    }
    using autopas::utils::ArrayUtils::operator<<;

    std::cout
            << "Iterations         : " << iterations << "\n"
            << "Particels per cell : " << numParticlesPerCell << "\n"
            << "Avgerage hit rate  : " << (static_cast<double>(calcsForceTotal) / calcsDistTotal) << "\n";

    printTimer();
}
