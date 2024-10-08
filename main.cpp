#include <iostream>

#ifdef ENABLE_FAPP
#include "fj_tool/fapp.h"
#endif

#include <molecularDynamicsLibrary/MoleculeLJ.h>
#include <molecularDynamicsLibrary/ArgonFunctor.h>

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

using ArgonFunctor = mdLib::ArgonFunctor<Particle, functorN3Modes, globals>;

double distSquared(std::array<double, 3> a, std::array<double, 3> b) {
    using autopas::utils::ArrayMath::sub;
    using autopas::utils::ArrayMath::dot;
    const auto c = sub(a, b);   // 3 FLOPS
    return dot(c, c);           // 3+2=5 FLOPs
}

std::map<std::string, autopas::utils::Timer> timer{
        {"Initialization",          autopas::utils::Timer()},
        {"Functor",                 autopas::utils::Timer()},
        {"Output",                  autopas::utils::Timer()},
        {"InteractionCounter",      autopas::utils::Timer()},
};

void resetTimer() {
    for (const auto &[name, t]: timer) {
        timer[name].reset();
    }
}

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

void applyFunctorOnTriplet(ArgonFunctor &functor, Particle &i, Particle &j, Particle &k) {
    timer.at("Functor").start();
    functor.AoSFunctor(i, j, k, newton3);
    timer.at("Functor").stop();
}

void applyFunctorOnParticleSet(ArgonFunctor &functor, std::vector<Particle> &particles) {
    for(std::size_t i = 0; i < particles.size(); ++i) {
        for (std::size_t j = i + 1; j < particles.size(); ++j) {
            for (std::size_t k = j + 1; k < particles.size(); ++k) {
                applyFunctorOnTriplet(functor, particles[i], particles[j], particles[k]);
            }
        }
    }
}

void csvOutput(ArgonFunctor &functor, std::vector<Particle> &particles) {
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

    // Open the file
    std::ofstream file("../Benchmark.csv");

    // Check if the file is opened successfully
    if (!file.is_open()) {
        std::cout << "Error opening file!" << std::endl;
        return 1;
    }

    // Write headers
    file << "Radius,Particles,FunctionCalls,Interactions,HitRate,FunctorRuntime" << std::endl;

    constexpr double cutoff{3.};
    // choose functor based on available architecture
    ArgonFunctor functor{cutoff};

    // define scenario
    constexpr size_t numParticles{100};

    const std::vector<double> radiuses{{1, 1.5, 2, 2.5}};//, 3, 3.5, 4, 4.5, 5, 5.5}};
    for (auto i=0 ; i<100; i++) {
        for (auto radius : radiuses) {
            // generate particles at random positions
            std::vector<Particle> particles{generateParticles(numParticles, radius)};
            // actual benchmark
            applyFunctorOnParticleSet(functor, particles);
            // gather data for analysis
            const auto [calcsDistTotal, calcsForceTotal] = countInteractions(particles, cutoff);

            //using autopas::utils::ArrayUtils::operator<<;

            file
                    << radius << ","
                    << numParticles << ","
                    << calcsDistTotal << ","
                    << calcsForceTotal << ","
                    << (static_cast<double>(calcsForceTotal) / calcsDistTotal) << ","
                    << static_cast<double>(timer["Functor"].getTotalTime()) * 1e-9 << std::endl;

            resetTimer();
        }
    }

    file.close();

}