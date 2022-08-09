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
    autopas::utils::Timer timerInitialization;
    timerInitialization.start();
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
    timerInitialization.stop();

    // actual benchmark
    autopas::utils::Timer timerFunctor;
    timerFunctor.start();
    // TODO offer option to also test FunctorSingle
    functor.SoAFunctorPair(cells[0]._particleSoABuffer, cells[1]._particleSoABuffer, newton3);
    timerFunctor.stop();

    // print particles to CSV for checking and prevent compiler from optimizing everything away.
    autopas::utils::Timer timerOutput;
    timerOutput.start();
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
    timerOutput.stop();

    // print timer and statistics

    // Macro because we use argument as part of variable names and string literals
#define printTimer(name) {                                      \
    std::cout                                                   \
    << std::setw(14)                                            \
    << std::left                                                \
    << #name                                                    \
    << " : "                                                    \
    << std::setprecision(3)                                     \
    << std::setw(8)                                             \
    << static_cast<double>(timer##name.getTotalTime()) * 10e-9  \
    << " [s]\n";                                                \
}

    using autopas::utils::ArrayUtils::operator<<;

    std::cout << "Particels per cell: " << numParticlesPerCell << "\n";
    printTimer(Initialization);
    printTimer(Functor);
    printTimer(Output);

// Todo: Flop counting, GFLOPs/sec
}
