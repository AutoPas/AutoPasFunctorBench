#include <iostream>

#ifdef ENABLE_FAPP
#include "fj_tool/fapp.h"
#endif

// clang-format off

#if defined ENABLE_HWY
#include <molecularDynamicsLibrary/LJFunctorHWY.h>
#elif defined __AVX__
#include <molecularDynamicsLibrary/LJFunctorAVX.h>
#elif __ARM_FEATURE_SVE
#include <molecularDynamicsLibrary/LJFunctorSVE.h>
#else
#error "Platform supports neither AVX nor ARM!"
#endif
// clang-format on

#include <molecularDynamicsLibrary/MoleculeLJ.h>
#include <autopas/cells/FullParticleCell.h>
#include <autopas/utils/Timer.h>
#include <fstream>

// type aliases for ease of use
using Particle = mdLib::MoleculeLJ;
using Cell = autopas::FullParticleCell<mdLib::MoleculeLJ>;
// some constants that define the benchmark
constexpr bool shift{false};
constexpr bool mixing{false};
constexpr autopas::FunctorN3Modes functorN3Modes{autopas::FunctorN3Modes::Both};
constexpr bool newton3{true};
constexpr bool globals{false};

#if defined ENABLE_HWY
using Functor = mdLib::LJFunctorHWY<Particle, shift, mixing, functorN3Modes, globals>;
#elif __AVX__
using Functor = mdLib::LJFunctorAVX<Particle, shift, mixing, functorN3Modes, globals>;
#elif __ARM_FEATURE_SVE
using Functor = mdLib::LJFunctorSVE<Particle, shift, mixing, functorN3Modes, globals> ;
#endif

enum FunctorType {
    pair,
    single,
    verlet
};

void checkFunctorType(const Functor &fun) {
    int identificationHits = 0;

#if defined ENABLE_HWY
    if (dynamic_cast<const mdLib::LJFunctorHWY<Particle, shift, mixing, functorN3Modes, globals> *>(&fun)) {
        std::cout << "Using HWY Functor" << std::endl;
        ++identificationHits;
    }
#endif

#if defined __AVX__ && not defined ENABLE_HWY
    if (dynamic_cast<const mdLib::LJFunctorAVX<Particle, shift, mixing, functorN3Modes, globals> *>(&fun)) {
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
                << static_cast<double>(timer[name].getTotalTime()) * 1e-9
                << " [s]\n";
    }
}

void initialization(Functor &functor, std::vector<Cell> &cells, const std::vector<size_t> &numParticlesPerCell,
                    double cutoff, double hitRate) {
    // initialize cells with randomly distributed particles

    // TODO : consider hitRate

    timer.at("Initialization").start();
    cells[0].reserve(numParticlesPerCell[0]);
    cells[1].reserve(numParticlesPerCell[1]);
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
                    particleId % 5};
            cells[cellId].addParticle(p);
        }
        functor.SoALoader(cells[cellId], cells[cellId]._particleSoABuffer, 0, false);
    }
    timer.at("Initialization").stop();
}

void applyFunctorPair(Functor &functor, std::vector<Cell> &cells) {
    timer.at("Functor").start();
#ifdef ENABLE_FAPP
    fapp_start("SoAFunctorPair", 1, 0);
#endif
    functor.SoAFunctorPair(cells[0]._particleSoABuffer, cells[1]._particleSoABuffer, newton3);
#ifdef ENABLE_FAPP
    fapp_stop("SoAFunctorPair", 1, 0);
#endif
    timer.at("Functor").stop();
}

void applyFunctorSingle(Functor &functor, std::vector<Cell> &cells) {
        timer.at("Functor").start();
#ifdef ENABLE_FAPP
    fapp_start("SoAFunctorSingle", 1, 0);
#endif
    functor.SoAFunctorSingle(cells[0]._particleSoABuffer, newton3);
#ifdef ENABLE_FAPP
    fapp_stop("SoAFunctorSingle", 1, 0);
#endif
    timer.at("Functor").stop();
}

void applyFunctorVerlet(Functor &functor, std::vector<Cell> &cells) {
        timer.at("Functor").start();
#ifdef ENABLE_FAPP
    fapp_start("SoAFunctorVerlet", 1, 0);
#endif
    // TODO : prepare neighbor list
#ifdef ENABLE_FAPP
    fapp_stop("SoAFunctorVerlet", 1, 0);
#endif
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
        for (size_t particleId = 0; particleId < cells[cellId].getNumberOfParticles(autopas::IteratorBehavior::owned); ++particleId) {
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

void printHelp() {
    std::cout << "Usage: \n"
        << "./AutoPasFunctorBench <functor type> <iterations> <numParticles> <hit rate> \n"
        << "Possible values:\n"
        << "functor type: \"single\", \"pair\", \"verlet\"\n"
        << "iterations, numParticles: int\n"
        << "hit rate: double"
        << std::endl;
}

std::tuple<FunctorType, size_t, size_t, double> readCliInput(int argc, char* argv[]) {

    FunctorType type = FunctorType::single;
    size_t iterations {1000};
    size_t numParticles {1000};
    double hitRate {0.5};

    if (argc < 5) {
        printHelp();
        exit(1);
    }

    std::vector<std::string> types {"single", "pair", "verlet"};

    if (types[0].compare(argv[1]) == 0) {
        type = single;
    }
    else if (types[1].compare(argv[1]) == 0) {
        type = pair;
    }
    else if (types[2].compare(argv[1]) == 0) {
        type = verlet;
    }
    else {
        std::cout << "Unkown functor type: " << argv[1] << std::endl;
        printHelp();
        exit(1);
    }

    try {
        iterations = std::stoi(argv[2]);
        numParticles = std::stoi(argv[3]);
        hitRate = std::stod(argv[4]);
    }
    catch (std::exception &e) {
        std::cout << "Could not parse iterations, numParticles or hitRate. Please make sure that the types do match." << std::endl;
        printHelp();
        exit(1);
    }

    return {type, iterations, numParticles, hitRate};
}

/**
 * Mini benchmark tool to estimate the inner most kernel performance of AutoPas
 * @return
 */
int main(int argc, char* argv[]) {

    auto [functorType, iterations, numParticles, hitRate] = readCliInput(argc, argv);

    constexpr double cutoff{3.}; // is also the cell size

    // choose functor based on available architecture
    // todo this is now hard-coded to have mixing - this should be somewhat more flexible
    ParticlePropertiesLibrary<double, size_t> PPL{cutoff};
    Functor functor{cutoff/*, PPL*/};

    checkFunctorType(functor);

    // 5 site types to provide some variation (requiring gathering that is somewhat similar to a realistic scenario)
    if constexpr (mixing) {
        PPL.addSiteType(0,1.,1.,1.);
        PPL.addSiteType(1,1.,1.,1.);
        PPL.addSiteType(2,1.,1.,1.);
        PPL.addSiteType(3,1.,1.,1.);
        PPL.addSiteType(4,1.,1.,1.);
        PPL.calculateMixingCoefficients();
    } else {
        constexpr double epsilon24{24.};
        constexpr double sigmaSquare{1.};
        functor.setParticleProperties(epsilon24, sigmaSquare);
    }



    // define scenario
    const std::vector<size_t> numParticlesPerCell{numParticles, numParticles};
    size_t calcsDistTotal{0};
    size_t calcsForceTotal{0};
    // repeat the whole experiment multiple times and average results
    std::vector<Cell> cells{2};

    initialization(functor, cells, numParticlesPerCell, cutoff, hitRate);

    switch (functorType)
    {
    case FunctorType::pair:
        for (size_t iteration = 0; iteration < iterations; ++iteration) {
            // actual benchmark
            applyFunctorPair(functor, cells);
        }
        break;
    case FunctorType::single:
        for (size_t iteration = 0; iteration < iterations; ++iteration) {
            // actual benchmark
            applyFunctorSingle(functor, cells);
        }
        break;
    case FunctorType::verlet:
        for (size_t iteration = 0; iteration < iterations; ++iteration) {
            // actual benchmark
            applyFunctorVerlet(functor, cells); // TODO : actual implementation
        }
        break;
    default:
        throw std::runtime_error("No functor type matched");
    }
    
    // print particles to CSV for checking and prevent compiler from optimizing everything away.
    csvOutput(functor, cells);

    // gather data for analysis
    const auto [calcsDist, calcsForce] = countInteractions(cells, cutoff);
    calcsDistTotal += calcsDist;
    calcsForceTotal += calcsForce;

    // print timer and statistics
    const auto gflops =
            static_cast<double>(calcsDistTotal * 8 + calcsForceTotal * (newton3 ? 18 : 15)) * 1e-9;
//    const auto gflops =
//            static_cast<double>(calcsDistTotal * 8 + calcsForceTotal * functor.getNumFlopsPerKernelCall()) * 1e-9;

    using autopas::utils::ArrayUtils::operator<<;

    std::cout
            << "Iterations         : " << iterations << "\n"
            << "Particels per cell : " << numParticlesPerCell << "\n"
            << "Avgerage hit rate  : " << (static_cast<double>(calcsForceTotal) / calcsDistTotal) << "\n"
            << "GFLOPs             : " << gflops << "\n"
            << "GFLOPs/sec         : " << (gflops / (timer.at("Functor").getTotalTime() * 1e-9)) << "\n";

    printTimer();
}
