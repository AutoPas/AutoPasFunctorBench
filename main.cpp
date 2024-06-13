#include <iostream>


#include <molecularDynamicsLibrary/KryptonPairFunctor.h>

#include <molecularDynamicsLibrary/MoleculeLJ.h>
#include <autopas/cells/FullParticleCell.h>
#include <autopas/utils/Timer.h>
#include <autopas/utils/ArrayMath.h>
#include <fstream>

// type aliases for ease of use
using Particle = mdLib::MoleculeLJ;
using Cell = autopas::FullParticleCell<mdLib::MoleculeLJ>;
// some constants that define the benchmark
constexpr bool shift{false};
constexpr bool mixing{false};
constexpr autopas::FunctorN3Modes functorN3Modes{autopas::FunctorN3Modes::Both};
constexpr bool newton3{true};
constexpr bool globals{true};

enum FunctorType {
    aos,
    pair,
    single,
    verlet
};

using Functor = mdLib::KryptonPairFunctor<Particle, functorN3Modes, globals>;

void checkFunctorType(const Functor &fun) {
    int identificationHits = 0;
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

void resetTimer() {
    for (auto &[_, t]: timer) {
        t.reset();
    }
}

void initialization(Functor &functor, FunctorType type, std::vector<Cell> &cells, std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>>& neighborLists,
                    const std::vector<size_t> &numParticlesPerCell, double cutoff, double interactionLengthSquare, double hitRate) {
    // initialize cells with randomly distributed particles

    timer.at("Initialization").start();

    double cellLength = 0.;

    switch (type)
    {
    case aos: {
        // this is a formula determined by regression (based on a mapping from hitrate to div with random sample values)
        double div = 2.86*(hitRate*hitRate*hitRate)-4.13*(hitRate*hitRate)+2.81*hitRate+0.42;
        cellLength = cutoff / div;
        break;
    }
    case pair: {
        // this is a formula determined by regression (based on a mapping from hitrate to div with random sample values)
        double div = 2.86*(hitRate*hitRate*hitRate)-4.13*(hitRate*hitRate)+2.81*hitRate+0.42;
        cellLength = cutoff / div;
        break;
    }
    case single: {
        double div = 2.72*(hitRate*hitRate*hitRate)-4.02*(hitRate*hitRate)+2.47*hitRate+0.09;
        cellLength = cutoff / div;
        break;
    }
    default:
        cellLength = cutoff * 10;
        break;
    }

    cells[0].reserve(numParticlesPerCell[0]);
    cells[1].reserve(numParticlesPerCell[1]);
    for (size_t cellId = 0; cellId < numParticlesPerCell.size(); ++cellId) {
        for (size_t particleId = 0; particleId < numParticlesPerCell[cellId]; ++particleId) {
            Particle p{
                    {
                            // particles are next to each other in X direction
                            rand() / static_cast<double>(RAND_MAX) * cellLength + cellLength * cellId,
                            rand() / static_cast<double>(RAND_MAX) * cellLength,
                            rand() / static_cast<double>(RAND_MAX) * cellLength,
                    },
                    {0., 0., 0.,},
                    // every cell gets its own id space
                    particleId + ((std::numeric_limits<size_t>::max() / numParticlesPerCell.size()) * cellId),
                    particleId % 5};
            cells[cellId].addParticle(p);
        }
    }

    // for verlet lists, only consider first cell
    for (size_t i = 0; i < numParticlesPerCell[0]; ++i) {
        for (size_t j = newton3 ? i + 1 : 0; j < numParticlesPerCell[0]; ++j) {
            if (i == j) {
                continue;
            }

            auto dr = autopas::utils::ArrayMath::sub(cells[0][i].getR(), cells[0][j].getR());
            double dr2 = autopas::utils::ArrayMath::dot(dr, dr);
            if (dr2 <= interactionLengthSquare) {
                neighborLists.at(i).push_back(j);
            }
        }
    }

    timer.at("Initialization").stop();
}

void applyFunctorAoS(Functor &functor, std::vector<Cell> &cells) {
    timer.at("Functor").start();
    for (auto &p1 : cells[0]) {
        for (auto &p2 : cells[1]) {
            functor.AoSFunctor(p1, p2, newton3);
        }
    }
    timer.at("Functor").stop();
}
void applyFunctorPair(Functor &functor, std::vector<Cell> &cells) {
    timer.at("Functor").start();
//    functor.SoAFunctorPair(cells[0]._particleSoABuffer, cells[1]._particleSoABuffer, newton3);
    timer.at("Functor").stop();
}

void applyFunctorSingle(Functor &functor, std::vector<Cell> &cells) {
    timer.at("Functor").start();
//    functor.SoAFunctorSingle(cells[0]._particleSoABuffer, newton3);
    timer.at("Functor").stop();
}

void applyFunctorVerlet(Functor &functor, std::vector<Cell> &cells, const std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> &neighborLists) {
    timer.at("Functor").start();
    for (size_t i = 0; i < neighborLists.size(); ++i) {
        functor.SoAFunctorVerlet(cells[0]._particleSoABuffer, i, neighborLists[i], newton3);
    }
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
//        functor.SoAExtractor(cells[cellId], cells[cellId]._particleSoABuffer, 0);
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

std::tuple<size_t, size_t> countInteractions(std::vector<Cell> &cells, std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> &neighborLists, FunctorType type, double cutoff) {
    timer.at("InteractionCounter").start();
    size_t calcsDist{0};
    size_t calcsForce{0};
    const auto cutoffSquared{cutoff * cutoff};

    switch (type)
    {
    case aos:
        for (const auto &p0: cells[0]) {
            for (const auto &p1: cells[1]) {
                ++calcsDist;
                if (distSquared(p0.getR(), p1.getR()) <= cutoffSquared) {
                    ++calcsForce;
                }
            }
        }
        break;

    case pair:
        for (const auto &p0: cells[0]) {
            for (const auto &p1: cells[1]) {
                ++calcsDist;
                if (distSquared(p0.getR(), p1.getR()) <= cutoffSquared) {
                    ++calcsForce;
                }
            }
        }
        break;

    case single:

        for (const auto &p0: cells[0]) {
            for (const auto &p1: cells[0]) {
                ++calcsDist;
                if (distSquared(p0.getR(), p1.getR()) <= cutoffSquared) {
                    ++calcsForce;
                }
            }
        }
        break;

    case verlet:

        for (int i = 0 ; i < neighborLists.size(); ++i) {
            const auto& p0 = cells[0][i];

            for (const auto& p1 : neighborLists[i]) {
                ++calcsDist;
                if (distSquared(p0.getR(), cells[0][p1].getR()) <= cutoffSquared) {
                    ++calcsForce;
                }
            }
        }

        break;
    
    default:
        break;
    }

    
    timer.at("InteractionCounter").stop();
    return {calcsDist, calcsForce};
}

void printHelp() {
    std::cout << "Usage: \n"
        << "./AutoPasFunctorBench <functor type> <repetitions> <iterations> <numParticles> <hit rate> <outfile>\n"
        << "Possible values:\n"
        << "functor type: \"aos\", \"single\", \"pair\", \"verlet\"\n"
        << "repetitions, iterations, numParticles: int\n"
        << "hit rate: double\n"
        << "outfile: std::string"
        << std::endl;
}

template <typename T>
void writeListToJson(const std::vector<T>& data, const std::string fileName) {
    std::ofstream outFile(fileName);
    if (outFile.is_open()) {
        outFile << "{\n  \"times\": [";

        for (size_t i = 0; i < data.size(); ++i) {
            outFile << data[i];
            if (i < data.size() - 1) {
                outFile << ", ";
            }
        }

        outFile << "],\n";
        outFile << "  \"avg\": ";

        T sum = std::accumulate(data.begin(), data.end(), static_cast<T>(0));
        outFile << (static_cast<double>(sum) / static_cast<double>(data.size()));
        outFile << "\n}";

        outFile.close();
        std::cout << "Data successfully saved" << std::endl;
    } else {
        std::cerr << "Could not open the file for writing!" << std::endl;
        exit(1);
    }
}

std::tuple<FunctorType, size_t, size_t, size_t, double, std::string> readCliInput(int argc, char* argv[]) {

    FunctorType type = FunctorType::aos;
    size_t repetitions {1};
    size_t iterations {1000};
    size_t numParticles {1000};
    double hitRate {0.5};
    std::string outfile {"benchmark.json"};

    if (argc < 7) {
        printHelp();
        exit(1);
    }

    std::vector<std::string> types {"aos", "single", "pair", "verlet"};
    if (types[0].compare(argv[1]) == 0) {
        type = aos;
    }
    else if (types[1].compare(argv[1]) == 0) {
        type = single;
    }
    else if (types[2].compare(argv[1]) == 0) {
        type = pair;
    }
    else if (types[3].compare(argv[1]) == 0) {
        type = verlet;
    }
    else {
        std::cout << "Unkown functor type: " << argv[1] << std::endl;
        printHelp();
        exit(1);
    }

    try {
        repetitions = std::stoi(argv[2]);
        iterations = std::stoi(argv[3]);
        numParticles = std::stoi(argv[4]);
        hitRate = std::stod(argv[5]);
        outfile = argv[6];
    }
    catch (std::exception &e) {
        std::cout << "Could not parse repetitions, iterations, numParticles or hitRate. Please make sure that the types do match." << std::endl;
        printHelp();
        exit(1);
    }

    return {type, repetitions, iterations, numParticles, hitRate, outfile};
}

/**
 * Mini benchmark tool to estimate the inner most kernel performance of AutoPas
 * @return
 */
int main(int argc, char* argv[]) {

    auto [functorType, repetitions, iterations, numParticles, hitRate, outfile] = readCliInput(argc, argv);

    const double cutoff{3.}; // is also the cell size
    const double skin{std::pow(1./hitRate, 1./3.)};
    const double interactionLengthSquare{(cutoff * skin) * (cutoff * skin)};

    std::vector<uint64_t> times {};
    times.reserve(repetitions);

    for (int n = 0; n < repetitions; ++n) {

        // choose functor
        Functor functor{cutoff};
        checkFunctorType(functor);

        // define scenario
        const std::vector<size_t> numParticlesPerCell{numParticles, numParticles};
        size_t calcsDistTotal{0};
        size_t calcsForceTotal{0};
        // repeat the whole experiment multiple times and average results
        std::vector<Cell> cells{2};
        std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> neighborLists (numParticles);

        initialization(functor, functorType, cells, neighborLists, numParticlesPerCell, cutoff, interactionLengthSquare, hitRate);

        switch (functorType)
        {
        case FunctorType::aos:
            for (size_t iteration = 0; iteration < iterations; ++iteration) {
                // actual benchmark
                applyFunctorAoS(functor, cells);
            }
            break;
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
                applyFunctorVerlet(functor, cells, neighborLists);
            }
            break;
        default:
            throw std::runtime_error("No functor type matched");
        }
        
        // print particles to CSV for checking and prevent compiler from optimizing everything away.
        csvOutput(functor, cells);

        // gather data for analysis
        const auto [calcsDist, calcsForce] = countInteractions(cells, neighborLists, functorType, cutoff);
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

        times.push_back(timer["Functor"].getTotalTime());
        std::cout << "-----------------------------------" << "\n";
        resetTimer();
    }

    writeListToJson<uint64_t>(times, outfile);
}
