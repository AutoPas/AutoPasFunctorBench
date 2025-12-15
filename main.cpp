#include <iostream>

#include <molecularDynamicsLibrary/MoleculeLJ.h>
#include <molecularDynamicsLibrary/AxilrodTellerMutoFunctor.h>

#include <autopas/cells/FullParticleCell.h>
#include <autopas/utils/Timer.h>
#include <fstream>
#include <random>

// type aliases for ease of use
using Particle = mdLib::MoleculeLJ;
using Cell = autopas::FullParticleCell<Particle>;

// some constants that define the benchmark
constexpr bool mixing{false};
constexpr autopas::FunctorN3Modes functorN3Modes{autopas::FunctorN3Modes::Both};
constexpr bool newton3{true};
constexpr bool globals{false};

// using ATM = mdLib::AxilrodTellerFunctor<Particle, false, functorN3Modes, globals>;
using ATM = mdLib::AxilrodTellerMutoFunctor<Particle, false, functorN3Modes, globals>;

enum FunctorMode {
    AOS,
    SOASINGLE,
    SOAPAIR,
    SOATRIPLE
};

constexpr std::array functorsToTest = {AOS, SOASINGLE, SOAPAIR, SOATRIPLE};

double distSquared(std::array<double, 3> a, std::array<double, 3> b) {
    using autopas::utils::ArrayMath::sub;
    using autopas::utils::ArrayMath::dot;
    const auto c = sub(a, b);   // 3 FLOPS
    return dot(c, c);           // 3+2=5 FLOPs
}

std::map<std::string, autopas::utils::Timer> timer{
        {"Initialization",          autopas::utils::Timer()},
        {"AoSFunctor Loop",      autopas::utils::Timer()},
        {"AoSFunctor on Particles", autopas::utils::Timer()},
        {"SoAFunctor Single", autopas::utils::Timer()},
        {"SoAFunctor Pair", autopas::utils::Timer()},
        {"SoAFunctor Triple", autopas::utils::Timer()},
        {"Output",                  autopas::utils::Timer()},
        {"InteractionCounter",      autopas::utils::Timer()},
};

void printTimers(const FunctorMode mode) {
    auto printTimer = [&] (const auto &name) {
        std::cout
                << std::setw(20)
                << std::left
                << name
                << " : "
                << std::setprecision(4)
                << std::setw(8)
                << static_cast<double>(timer[name].getTotalTime()) * 1e-6
                << " [ms]\n";
    };

    // printTimer("Initialization");
    switch (mode) {
        case AOS:
            printTimer("AoSFunctor Loop");
            printTimer("AoSFunctor on Particles");
            break;
        case SOASINGLE:
            printTimer("SoAFunctor Single");
            break;
        case SOAPAIR:
            printTimer("SoAFunctor Pair");
            break;
        case SOATRIPLE:
            printTimer("SoAFunctor Triple");
            break;
    }
    // printTimer("Output");
    printTimer("InteractionCounter");
}

void generateParticles(ATM &functor, std::vector<Cell> &cells, const size_t numberOfParticles,
    const double cutoff, const FunctorMode mode){
    timer.at("Initialization").start();
    // generate randomly distributed particles
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, cutoff);

    auto fillCellWithParticles = [&] (Cell& cell, double xShift, double yShift, double zShift) {
        for (size_t particleId = 0; particleId < numberOfParticles; ++particleId) {
            Particle p{
                { dis(gen) + xShift, dis(gen) + yShift, dis(gen)},
                {0., 0., 0.,},
                particleId,
                0};
            cell.addParticle(p); // Add particle to the current cell
        }
    };

    switch (mode) {
        case AOS:
        fillCellWithParticles(cells[0], 0., 0., 0.);
            break;
        case SOATRIPLE:
            fillCellWithParticles(cells[2], 0., cutoff, 0.);
            functor.SoALoader(cells[2], cells[2]._particleSoABuffer, 0, false);
        case SOAPAIR:
            fillCellWithParticles(cells[1], cutoff, 0., 0.);
            functor.SoALoader(cells[1], cells[1]._particleSoABuffer, 0, false);
        case SOASINGLE:
            fillCellWithParticles(cells[0], 0., 0., 0.);
            functor.SoALoader(cells[0], cells[0]._particleSoABuffer, 0, false);
            break;
    }
    timer.at("Initialization").stop();
}

void applyAoSFunctor(ATM &functor, Cell &cell) {
    timer.at("AoSFunctor Loop").start();
    for(std::size_t i = 0; i < cell.size(); ++i) {
        for (std::size_t j = i + 1; j < cell.size(); ++j) {
            for (std::size_t k = j + 1; k < cell.size(); ++k) {
                timer.at("AoSFunctor on Particles").start();
                functor.AoSFunctor(cell[i], cell[j], cell[k], newton3);
                timer.at("AoSFunctor on Particles").stop();
            }
        }
    }
    timer.at("AoSFunctor Loop").stop();
}

void applySoAFunctorSingle(ATM &functor, Cell &cell) {
    timer.at("SoAFunctor Single").start();
    functor.SoAFunctorSingle(cell._particleSoABuffer, newton3);
    timer.at("SoAFunctor Single").stop();
}

void applySoAFunctorPair(ATM &functor, Cell &cell1, Cell &cell2) {
    timer.at("SoAFunctor Pair").start();
    functor.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, newton3);
    timer.at("SoAFunctor Pair").stop();
}

void applySoAFunctorTriple(ATM &functor, Cell &cell1, Cell &cell2, Cell &cell3) {
    timer.at("SoAFunctor Triple").start();
    functor.SoAFunctorTriple(cell1._particleSoABuffer, cell2._particleSoABuffer, cell3._particleSoABuffer, newton3);
    timer.at("SoAFunctor Triple").stop();
}

void applyFunctorOnParticles(ATM &functor, std::vector<Cell> &cells, const FunctorMode mode) {
    switch (mode) {
        case AOS:
            applyAoSFunctor(functor, cells[0]);
            break;
        case SOASINGLE:
            applySoAFunctorSingle(functor, cells[0]);
            break;
        case SOAPAIR:
            applySoAFunctorPair(functor, cells[0], cells[1]);
            break;
        case SOATRIPLE:
            applySoAFunctorTriple(functor, cells[0], cells[1], cells[2]);
            break;
    }
}

// void csvOutput(ATM &functor, std::vector<Cell> &cells) {
//     timer.at("Output").start();
//     std::ofstream csvFile("particles.csv");
//     if (not csvFile.is_open()) {
//         throw std::runtime_error("FILE NOT OPEN!");
//     }
//     for (size_t cellId = 0; cellId < cells.size(); ++cellId) {
//         functor.SoAExtractor(cells[cellId], cells[cellId]._particleSoABuffer, 0);
//         csvFile << "ParticleId,rX,rY,rZ,fX,fY,fZ\n";
//         for (auto p : cells[0]) {
//             using autopas::utils::ArrayUtils::to_string;
//             csvFile << p.getID() << ","
//                     << to_string(p.getR(), ",", {"", ""}) << ","
//                     << to_string(p.getF(), ",", {"", ""})
//                     << "\n";
//         }
//         csvFile.close();
//         timer.at("Output").stop();
//     }
// }

std::tuple<size_t, size_t> countInteractions(std::vector<Cell> &cells, const double cutoff, const FunctorMode mode) {
    timer.at("InteractionCounter").start();
    size_t calcsDist{0};
    size_t calcsForce{0};
    const auto cutoffSquared{cutoff * cutoff};

    switch (mode) {
        case AOS:
        case SOASINGLE:
            for (size_t i = 0; i < cells[0].size(); i++) {
                for (size_t j = i + 1; j < cells[0].size(); j++) {
                    if (distSquared(cells[0][i].getR(), cells[0][j].getR()) > cutoffSquared) continue;
                    for (size_t k = j + 1; k < cells[0].size(); k++) {
                        if (distSquared(cells[0][i].getR(), cells[0][k].getR()) > cutoffSquared) continue;
                        if (distSquared(cells[0][j].getR(), cells[0][k].getR()) > cutoffSquared) continue;
                        ++calcsForce;
                    }
                }
            }
            calcsDist = cells[0].size() * (cells[0].size() - 1) * (cells[0].size() - 2) / 6;
            break;
        case SOAPAIR:
            for (size_t i = 0; i < cells[0].size(); i++) {
                for (size_t j = i + 1; j < cells[0].size(); j++) {
                    if (distSquared(cells[0][i].getR(), cells[0][j].getR()) > cutoffSquared) continue;
                    for (size_t k = 0; k < cells[1].size(); k++) {
                        if (distSquared(cells[0][i].getR(), cells[1][k].getR()) > cutoffSquared) continue;
                        if (distSquared(cells[0][j].getR(), cells[1][k].getR()) > cutoffSquared) continue;
                        ++calcsForce;
                    }
                }
                for (size_t j = 0; j < cells[1].size(); j++) {
                    if (distSquared(cells[0][i].getR(), cells[1][j].getR()) > cutoffSquared) continue;
                    for (size_t k = j + 1; k < cells[1].size(); k++) {
                        if (distSquared(cells[0][i].getR(), cells[1][k].getR()) > cutoffSquared) continue;
                        if (distSquared(cells[1][j].getR(), cells[1][k].getR()) > cutoffSquared) continue;
                        ++calcsForce;
                    }
                }
            }
            calcsDist = cells[0].size() * cells[1].size() * (cells[0].size() + cells[1].size() - 2) / 2;
            break;
        case SOATRIPLE:
            for (size_t i = 0; i < cells[0].size(); i++) {
                for (size_t j = 0; j < cells[1].size(); j++) {
                    if (distSquared(cells[0][i].getR(), cells[1][j].getR()) > cutoffSquared) continue;
                    for (size_t k = 0; k < cells[2].size(); k++) {
                        if (distSquared(cells[0][i].getR(), cells[2][k].getR()) > cutoffSquared) continue;
                        if (distSquared(cells[1][j].getR(), cells[2][k].getR()) > cutoffSquared) continue;
                        ++calcsForce;
                    }
                }
            }
            calcsDist = cells[0].size() * cells[1].size() * cells[2].size();
            break;
        default:
            break;
    }
    timer.at("InteractionCounter").stop();
    return {calcsDist, calcsForce};
}

/**
 * Mini benchmark tool to estimate the inner most kernel performance of AutoPas
 * @return
 */
int main() {
    using autopas::utils::ArrayUtils::operator<<;

    constexpr double cutoff{3.};
    constexpr double nu{1.0};

    ATM functor{cutoff};
    functor.setParticleProperties(nu);

    // define scenario
    constexpr size_t numParticles{150};
    constexpr size_t iterations{100};

    size_t calcsDistTotal{0};
    size_t calcsForceTotal{0};

    std::cout << functor.getName() << " Benchmark: " <<
        "\nParticles per Cell: " << numParticles <<
            "\nIteration Average: " << iterations << "\n\n";

    // repeat the whole experiment multiple times and average results
    for (const auto &functorMode : functorsToTest) {
        for (size_t iteration = 0; iteration < iterations; ++iteration) {

            std::vector<Cell> cells{3};
            generateParticles(functor, cells, numParticles, cutoff, functorMode);

            // actual benchmark
            applyFunctorOnParticles(functor, cells, functorMode);

            // print particles to CSV
            // csvOutput(functor, cells);

            // gather data for analysis
            const auto [calcsDist, calcsForce] = countInteractions(cells, cutoff, functorMode);
            calcsDistTotal += calcsDist;
            calcsForceTotal += calcsForce;
        }

        std::string functorModeName;
        switch (functorMode) {
            case AOS: functorModeName = "AoSFunctor"; break;
            case SOASINGLE: functorModeName = "SoAFunctorSingle"; break;
            case SOAPAIR: functorModeName = "SoAFunctorPair"; break;
            case SOATRIPLE: functorModeName = "SoAFunctorTriple"; break;
            default: functorModeName = "unknown"; break;
        }
        std::cout << "\n------\nStatistics for " << functorModeName << "\n";

        printTimers(functorMode);

        std::cout << "Average hit rate     : " << static_cast<double>(calcsForceTotal) / static_cast<double>(calcsDistTotal) * 100 << " %\n"
                  << "Interactions         : " << calcsForceTotal / iterations << "\n\n";
    }
}