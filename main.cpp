#include <iostream>
#include <unistd.h> // for getopt
#include <fstream>
#include <variant>

// Include the functor(s) to test
#include <molecularDynamicsLibrary/AxilrodTellerMutoMultisiteFunctor.h>
#include <molecularDynamicsLibrary/MethaneMultisitePairwiseFunctor.h>
#include <molecularDynamicsLibrary/LJFunctor.h>

#include <molecularDynamicsLibrary/MultisiteMoleculeLJ.h>
#include <autopas/cells/FullParticleCell.h>
#include <autopas/utils/Timer.h>
#include <autopas/utils/ArrayMath.h>

// type aliases for ease of use
using Particle = mdLib::MultisiteMoleculeLJ;
using Cell = autopas::FullParticleCell<mdLib::MultisiteMoleculeLJ>;
// some constants that define the experiments
constexpr autopas::FunctorN3Modes functorN3Modes{autopas::FunctorN3Modes::Both};
constexpr bool newton3{true};
constexpr bool globals{true};

// Chose the Functors
using ATMFunctor = mdLib::AxilrodTellerMutoMultisiteFunctor<Particle, false, true, functorN3Modes, globals>;
using MethanePairFunctor = mdLib::MethaneMultisitePairwiseFunctor<Particle, false, functorN3Modes, globals>;
using LJFunctor = mdLib::LJFunctor<mdLib::MoleculeLJ, true, false, functorN3Modes, true, true>;

enum class FunctorOption {
    LJFunctor,
    MethanePairFunctor,
    ATMFunctor,
};

enum class GeometryOption {
    linear,
    equilateral,
    isosceles
};

std::string getCSVFileName(FunctorOption functorOption, GeometryOption geometryOption) {
    std::stringstream filenameStream;
    switch (functorOption) {
        case FunctorOption::LJFunctor:
            std::cout << "Using the Lennard-Jones Functor\n";
            filenameStream << "2b-lj";
        break;
        case FunctorOption::MethanePairFunctor:
            std::cout << "Using the pairwise Methane Functor\n";
            filenameStream << "methane-2b";
        break;
        case FunctorOption::ATMFunctor:
            std::cout << "Using the Axilrod-Teller-Muto Functor\n";
            filenameStream << "3b-atm";
        break;
        default:
            return "";
    }

    switch (geometryOption) {
        case GeometryOption::linear:
            std::cout << "Using a linear geometry\n";
            filenameStream << "-linear";
        break;
        case GeometryOption::equilateral:
            std::cout << "Using a equilateral triangle geometry\n";
            filenameStream << "-equilateral";
        break;
        case GeometryOption::isosceles:
            std::cout << "Using a isosceles triangle geometry\n";
            filenameStream << "-isosceles";
        break;
        default:
            return "";
    }
    filenameStream << ".csv";
    return filenameStream.str();
}

void setEquilateral(Particle &p2, Particle &p3, double distance) {
    p2.setR({0., distance, 0.});
    p3.setR({std::sqrt(3)/2 * distance, 0.5 * distance, 0.});
}

void setLinear(Particle &p1, Particle &p2, Particle &p3, double distance) {
    p1.setR({10.0, 0., 0.});
    p2.setR({10.0 + distance, 0., 0.});
    p3.setR({9.8, 0., 0.});
}

void setIsosceles(Particle& p1, Particle &p2, Particle &p3, double distance) {
    p1.setR({distance, 2.0, 0.0});
    p2.setR({0.0, 0.0, 0.0});
    p3.setR({0.0, 4.0, 0.0});
}

std::vector<std::array<double, 3>> createTetrahedralSites(double factor) {
    constexpr double distCH = 1.099;
    const auto siteDist = factor * distCH;
    std::vector<std::array<double, 3>> sitePositions{{-siteDist, 0.0, 0.0}};
    sitePositions.push_back({siteDist / 3., siteDist * (2.* std::sqrt(2.) / 3.), 0.0});
    sitePositions.push_back({siteDist / 3., - siteDist * std::sqrt(2.) / 3., siteDist * std::sqrt(6.) / 3.});
    sitePositions.push_back({siteDist / 3., - siteDist * std::sqrt(2.) / 3., -siteDist * std::sqrt(6.) / 3.});
    return sitePositions;
}


/**
 * Mini validation setup that writes functor results to a CSV
 * @return
 */
int main(int argc, char* argv[]) {

    FunctorOption functorOption{};
    GeometryOption geometryOption{};

    if (argc < 2) {
        std::cerr << "Usage: -f <functor> [-g <geometry>]" << std::endl;
        std::cerr << "functor: 0 : Lennard-Jones (default); 1 : Methane Multisite Pair; 2 : Axilrod-Teller\n";
        std::cerr << "Geometry: 0 : linear (default); 1 : equilateral; 2 : isosceles\n";
        return 1;
    }

    int option;
    while ((option = getopt(argc, argv, "f:g:")) != -1) {
        switch (option) {
            case 'f':
                functorOption = static_cast<FunctorOption>(std::stoi(optarg));
            break;
            case 'g':
                geometryOption = static_cast<GeometryOption>(std::stoi(optarg));
            break;
            default:
                std::cerr << "Usage: -f <functor> [-g <geometry>]" << std::endl;
                std::cerr << "functor: 0 : Lennard-Jones (default); 1 : Methane Multisite Pair; 2 : Axilrod-Teller\n";
                std::cerr << "Geometry: 0 : linear (default); 1 : equilateral; 2 : isosceles\n";
            return 1;
        }
    }

    constexpr double cutoff{300.};
    constexpr double epsilon = 200.8753;
    constexpr double r_eps = 4.015802;
    const double sigma = r_eps / (std::pow(2.0, 1.0 / 6.0));
    constexpr double nu = 1.61525e6;
    std::cout << "Sigma: " << sigma << std::endl;
    auto particle1 = Particle({0., 0., 0.}, {0., 0., 0.}, {1., 0., 0., 0.}, {0., 0., 0.}, 0, 0);
    auto particle2 = Particle({0., 0., 0.}, {0., 0., 0.}, {0., 0., 0., 1.}, {0., 0., 0.}, 0, 0);
    auto particle3 = Particle({0., 0., 0.}, {0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0.}, 0, 0);

    std::string filename = getCSVFileName(functorOption, geometryOption);

    ParticlePropertiesLibrary ppLibrary{cutoff};
    ppLibrary.addSiteType(0, 16.04);
    ppLibrary.addLJParametersToSite(0, epsilon, sigma);
    ppLibrary.addSiteType(1, 0);
    ppLibrary.addLJParametersToSite(1, epsilon, sigma);
    ppLibrary.addSiteType(2, 0);
    ppLibrary.addLJParametersToSite(2, epsilon, sigma);

    std::vector<std::array<double, 3>> positions{{0., 0., 0.}};
    auto positionsH088 = createTetrahedralSites(0.88);
    auto positionsH066 = createTetrahedralSites(-0.66);
    positions.insert(positions.end(), positionsH088.begin(), positionsH088.end());
    positions.insert(positions.end(), positionsH066.begin(), positionsH066.end());
    std::vector<size_t> siteIDs{0, 1, 1, 1, 1, 2, 2, 2, 2};
    ppLibrary.addMolType(0, siteIDs, positions, {1., 1., 1.});

    // Create the selected functor as a std::variant
    auto createFunctor = [&] () -> std::variant<LJFunctor, MethanePairFunctor, ATMFunctor> {
        switch (functorOption) {
            case FunctorOption::LJFunctor: {
                LJFunctor ljfunctor{cutoff};
                ljfunctor.setParticleProperties(24 * epsilon, sigma * sigma);
                return ljfunctor;
            }
            case FunctorOption::MethanePairFunctor:
                return MethanePairFunctor{cutoff, ppLibrary};
            case FunctorOption::ATMFunctor: {
                ATMFunctor atmfunctor{cutoff};
                atmfunctor.setParticleProperties(nu);
                return atmfunctor;
            }
        }
        throw std::invalid_argument("Unknown functor option!");
    };

    auto functor = createFunctor();

    // Open the file
    std::ofstream file(filename);

    // Check if the file is opened successfully
    if (!file.is_open()) {
        std::cout << "Error opening file!" << std::endl;
        return 1;
    }

    // Write headers
    file << "Distance [A],Energy [K],Force_x [K/A],Force_y [K/A],Force_z [K/A]" << std::endl;

    double distance = 2.5;
    constexpr double maxDistance = 8.5;
    constexpr double increment = 0.01;
    double epot{};

    while (distance <= maxDistance) {
        // Reset forces
        particle1.setF({0.,0.,0.});
        particle2.setF({0.,0.,0.});
        particle3.setF({0.,0.,0.});

        // Move particles
        switch (geometryOption) {
            case GeometryOption::linear: {
                setLinear(particle1, particle2, particle3, distance);
                break;
            }
            case GeometryOption::equilateral: {
                setEquilateral(particle2, particle3, distance);
                break;
            }
            case GeometryOption::isosceles: {
                setIsosceles(particle1, particle2, particle3, distance);
                break;
            }
        }

        // Do one functor interaction (AoS)
        std::visit([&](auto&& f) {
            f.initTraversal();

            using T = std::decay_t<decltype(f)>;
            if constexpr (std::is_same_v<T, LJFunctor> || std::is_same_v<T, MethanePairFunctor>) {
                // Call AoSFunctor with two particles for LJFunctor and KrPairFunctor
                f.AoSFunctor(particle1, particle2, newton3);
            } else if constexpr (std::is_same_v<T, ATMFunctor>) {
                // Call AoSFunctor with three particles for ATMFunctor and KrEATMFunctor
                f.AoSFunctor(particle1, particle2, particle3, newton3);
            }

            f.endTraversal(newton3);
            epot = f.getPotentialEnergy();
        }, functor);

        // Results
        auto force = particle1.getF();

        // Write to CSV
        file << distance << std::fixed << std::setprecision(15) << "," << epot << "," << force[0] << "," << force[1] << "," << force[2]  << std::endl;
        distance += increment;
    }

    file.close();
    return 0;
}
