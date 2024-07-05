#include <iostream>

// Include the functor(s) to test
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
constexpr autopas::FunctorN3Modes functorN3Modes{autopas::FunctorN3Modes::Both};
constexpr bool newton3{true};
constexpr bool globals{true};

// Chose the Functor
using Functor = mdLib::KryptonPairFunctor<Particle, false, functorN3Modes, globals>;

/**
 * Mini validation setup that writes functor results to a CSV
 * @return
 */
int main(int argc, char* argv[]) {

    const double cutoff{300.};
    Functor functor{cutoff};

    auto particle1 = Particle({0., 0., 0.}, {0., 0., 0.}, 0, 0);
    auto particle2 = Particle({0., 0., 0.}, {0., 0., 0.}, 0, 0);

    // Open the file
    std::ofstream file("kryptonPair2.csv");

    // Check if the file is opened successfully
    if (!file.is_open()) {
        std::cout << "Error opening file!" << std::endl;
        return 1;
    }

    // Write headers
    file << "Distance [nm],Force [K/nm],Energy [K]" << std::endl;

    double epot{};
    double virial{};
    size_t resolution = 1000;
    double increment = 0.7 / resolution;
    double distance = 1.10798137169225 * 0.35;
    for (auto i = 0; i < resolution; i++) {
        // move second particle along x-Axis
        particle2.setF({0.,0.,0.});
        particle2.setR({distance, 0., 0.});
        // Compute the interactions
        functor.initTraversal();
        functor.AoSFunctor(particle1, particle2, newton3);
        functor.endTraversal(newton3);
        // Results
        auto force = particle1.getF();
        epot = functor.getPotentialEnergy();
        virial = functor.getVirial();

        // std::cout << "Distance,   " << distance << std::endl;
        // std::cout << "Force,      " << force[0] << std::endl;
        // std::cout << "Potential Energy:" << epot << std::endl;

        // Write to CSV
        file << std::fixed << std::setprecision(15) << distance << "," << force[0] << "," << epot << "," << virial << std::endl;
        distance += increment;
    }

    file.close();
}
