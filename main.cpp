#include <iostream>

#ifdef ENABLE_FAPP
#include "fj_tool/fapp.h"
#endif

#include <molecularDynamicsLibrary/MoleculeLJ.h>
#include <molecularDynamicsLibrary/AxilrodTellerFunctor.h>
#include <molecularDynamicsLibrary/TriwiseLUT.h>

#include <autopas/cells/FullParticleCell.h>
#include <autopas/utils/Timer.h>
#include <fstream>

// type aliases for ease of use
using Particle = mdLib::MoleculeLJ;
// some constants that define the benchmark
constexpr bool shift{false};
constexpr bool mixing{false};
constexpr autopas::FunctorN3Modes functorN3Modes{autopas::FunctorN3Modes::Both};
constexpr bool newton3{true};
constexpr bool globals{false};

// using ATM = mdLib::AxilrodTellerFunctor<Particle, false, false, functorN3Modes, globals>;
using ATM = mdLib::AxilrodTellerFunctor<Particle, false, true, functorN3Modes, globals>;

double distSquared(std::array<double, 3> a, std::array<double, 3> b) {
    using autopas::utils::ArrayMath::sub;
    using autopas::utils::ArrayMath::dot;
    const auto c = sub(a, b);   // 3 FLOPS
    return dot(c, c);           // 3+2=5 FLOPs
}

std::map<std::string, autopas::utils::Timer> timer{
        {"Initialization",          autopas::utils::Timer()},
        {"Functor on Triplet",      autopas::utils::Timer()},
        {"Functor on Particle Set", autopas::utils::Timer()},
        {"Output",                  autopas::utils::Timer()},
        {"InteractionCounter",      autopas::utils::Timer()},
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

std::vector<Particle> generateParticles(const size_t numberOfParticles, double cutoff){
    // generate randomly distributed particles
    timer.at("Initialization").start();
    std::vector<Particle> particles;
    for (size_t particleId = 0; particleId < numberOfParticles; ++particleId) {
        Particle p{
                {
                        rand() / static_cast<double>(RAND_MAX) * cutoff,
                        rand() / static_cast<double>(RAND_MAX) * cutoff,
                        rand() / static_cast<double>(RAND_MAX) * cutoff,
                },
                {0., 0., 0.,},
                // every cell gets its own id space
                particleId,
                0};
        particles.push_back(p);
    }
    timer.at("Initialization").stop();
    return particles;
}

void applyFunctorOnTriplet(ATM &functor, Particle &i, Particle &j, Particle &k) {
    timer.at("Functor on Triplet").start();
    functor.AoSFunctor(i, j, k, newton3);
    timer.at("Functor on Triplet").stop();
}

void applyFunctorOnParticleSet(ATM &functor, std::vector<Particle> &particles) {
    timer.at("Functor on Particle Set").start();
    for(std::size_t i = 0; i < particles.size(); ++i) {
        for (std::size_t j = i + 1; j < particles.size(); ++j) {
            for (std::size_t k = j + 1; k < particles.size(); ++k) {
                applyFunctorOnTriplet(functor, particles[i], particles[j], particles[k]);
            }
        }
    }
    timer.at("Functor on Particle Set").stop();
}

void csvOutput(ATM &functor, std::vector<Particle> &particles) {
    timer.at("Output").start();
    std::ofstream csvFile("particles.csv");
    if (not csvFile.is_open()) {
        throw std::runtime_error("FILE NOT OPEN!");
    }
    csvFile << "ParticleId,rX,rY,rZ,fX,fY,fZ\n";
    for (auto p : particles) {
        using autopas::utils::ArrayUtils::to_string;
        csvFile << p.getID() << ","
                << to_string(p.getR(), ",", {"", ""}) << ","
                << to_string(p.getF(), ",", {"", ""})
                << "\n";
    }
    csvFile.close();
    timer.at("Output").stop();
}

std::tuple<size_t, size_t> countInteractions(std::vector<Particle> &particles, double cutoff) {
    timer.at("InteractionCounter").start();
    size_t calcsDist{0};
    size_t calcsForce{0};
    const auto cutoffSquared{cutoff * cutoff};
    for (const auto &p0: particles) {
        for (const auto &p1: particles) {
            for (const auto &p2: particles) {
                if (p0 != p1 && p0 != p2 && p1 != p2) {
                    ++calcsDist;
                    if (distSquared(p0.getR(), p1.getR()) <= cutoffSquared) {
                        ++calcsForce;
                    }
                }
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

    constexpr double cutoff{3.};
    // choose functor based on available architecture

    mdLib::TriwiseLUT lut(500);
    ATM functor{cutoff, &lut};
    functor.setParticleProperties(1.0);

    // define scenario
    constexpr size_t numParticles{150};
    constexpr size_t iterations{100};
    size_t calcsDistTotal{0};
    size_t calcsForceTotal{0};
    // repeat the whole experiment multiple times and average results
    for (size_t iteration = 0; iteration < iterations; ++iteration) {
        std::vector<Particle> particles{generateParticles(numParticles, cutoff)};

        // actual benchmark
        applyFunctorOnParticleSet(functor, particles);

        // print particles to CSV
        csvOutput(functor, particles);

        // gather data for analysis
        const auto [calcsDist, calcsForce] = countInteractions(particles, cutoff);
        calcsDistTotal += calcsDist;
        calcsForceTotal += calcsForce;
    }

    // print timer and statistics
    // const auto gflops =
    //         static_cast<double>(calcsDistTotal * 8 + calcsForceTotal * functor.getNumFlopsPerKernelCall(0, 0, 0, true)) * 10e-9;

    using autopas::utils::ArrayUtils::operator<<;

    std::cout
            << "Iterations         : " << iterations << "\n"
            << "Particels per cell : " << numParticles << "\n"
            << "Avgerage hit rate  : " << (static_cast<double>(calcsForceTotal) / calcsDistTotal) << "\n"
            << "Interaction Counter: " << calcsForceTotal / 1000000 << " M\n\n";
            // << "GFLOPs             : " << gflops << "\n"
            // << "GFLOPs/sec         : " << (gflops / (timer.at("Functor on Triplet").getTotalTime() * 10e-9)) << "\n";

    printTimer();
}