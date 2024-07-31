#include <iostream>

// Include the functor(s) to test
#include <molecularDynamicsLibrary/AxilrodTellerFunctor.h>

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
using Functor = mdLib::AxilrodTellerFunctor<Particle, false, functorN3Modes, globals>;

/**
 * Mini validation setup that writes functor results to a CSV
 * @return
 */
int main(int argc, char* argv[]) {

    const double cutoff{300.};
    Functor functor{cutoff};
    functor.setParticleProperties(1.6152500E-3);

    auto particle1 = Particle({0., 0., 0.}, {0., 0., 0.}, 0, 0);
    auto particle2 = Particle({0., 0., 0.}, {0., 0., 0.}, 0, 0);
    auto particle3 = Particle({0., 0., 0.}, {0., 0., 0.}, 0, 0);

    // Open the file
    std::ofstream file("../AxilrodTeller/EquilateralGeometry.csv");

    // Check if the file is opened successfully
    if (!file.is_open()) {
        std::cout << "Error opening file!" << std::endl;
        return 1;
    }

    // Write headers
    file << "Distance,Energy,Force_x,Force_y,Force_z" << std::endl;

    double epot{};
    size_t resolution = 1000;
    double increment = 0.001;
    double distance = 0.3;
    for (auto i = 0; i < resolution; i++) {
        // move second particle
        particle2.setF({0.,0.,0.});
        particle2.setR({0., distance, 0.});
        // move third particle
        particle3.setF({0.,0.,0.});
        particle3.setR({std::sqrt(3)/2 * distance, 0.5 * distance, 0.});
        // Compute the interactions
        functor.initTraversal();
        functor.AoSFunctor(particle1, particle2, particle3, newton3);
        functor.endTraversal(newton3);
        // Results
        auto force = particle1.getF();
        epot = functor.getPotentialEnergy();

        // std::cout << "Distance,   " << distance << std::endl;
        // std::cout << "Force,      " << force[0] << std::endl;
        // std::cout << "Potential Energy:" << epot << std::endl;

        // Write to CSV
        file << std::fixed << std::setprecision(15) << distance << "," << epot << "," << force[0] << "," << force[1] << "," << force[2]  << std::endl;
        distance += increment;
    }

    file.close();
}
