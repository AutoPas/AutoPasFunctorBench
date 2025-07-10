#include <iostream>

#include <filesystem>
// Include the functor(s) to test
#include <molecularDynamicsLibrary/AxilrodTellerFunctor.h>
#include <molecularDynamicsLibrary/LJFunctor.h>
//#include <molecularDynamicsLibrary/ArgonFunctor.h>
#include <molecularDynamicsLibrary/KryptonExtendedATMFunctor.h>
#include <molecularDynamicsLibrary/LUT3B.h>

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
constexpr bool flops{true};
constexpr bool useLUTglob{false};
namespace mdFlexibleTypeDefs {
/**
 * If AutoPas is compiled with FLOP logging enabled, use functors with FLOP counting enabled.
 */
    constexpr bool countFLOPs =

            false;


/**
 * If md-flexible is compiled with globals calculations enabled, use functors which calculate globals.
 */
    constexpr bool calcGlobals =

            true;

}  // namespace mdFlexibleTypeDefs
// Chose the Functor
//using Functor = mdLib::AxilrodTellerFunctor<Particle, false, functorN3Modes, globals>;
//using Functor = mdLib::AxilrodTellerFunctor<Particle, false,true, functorN3Modes,globals, flops>


//WORKS FOR feat/3xa/lut
//using Functor = mdLib::AxilrodTellerFunctor<mdLib::MoleculeLJ, false, false, autopas::FunctorN3Modes::Both,mdFlexibleTypeDefs::calcGlobals, mdFlexibleTypeDefs::countFLOPs>;


//for auto-personal
//using Functor = mdLib::AxilrodTellerFunctor<mdLib::MoleculeLJ, false, true,false, autopas::FunctorN3Modes::Both,mdFlexibleTypeDefs::calcGlobals, mdFlexibleTypeDefs::countFLOPs>;
//using Functor = mdLib::LJFunctor<mdLib::MoleculeLJ, false, false,true, useLUTglob, autopas::FunctorN3Modes::Both, mdFlexibleTypeDefs::calcGlobals, mdFlexibleTypeDefs::countFLOPs, true>;


//for auto-personal argon
//using Functor  =  mdLib::ArgonFunctor<mdLib::MoleculeLJ, autopas::FunctorN3Modes::Both,true, false,mdFlexibleTypeDefs::calcGlobals >;


//for auto-personal krypton
using Functor =  mdLib::KryptonExtendedATMFunctor<mdLib::MoleculeLJ, autopas::FunctorN3Modes::Both, true, false, mdFlexibleTypeDefs::calcGlobals, mdFlexibleTypeDefs::countFLOPs>;

// noble gas krypton
//using Functor =  mdLib::KryptonExtendedATMFunctor<mdLib::MoleculeLJ, autopas::FunctorN3Modes::Both, mdFlexibleTypeDefs::calcGlobals, mdFlexibleTypeDefs::countFLOPs>;



//for Force_LUT
//using Functor = mdLib::AxilrodTellerFunctor<Particle, true,false, autopas::FunctorN3Modes::Both, false>
//using Functor = mdLib::AxilrodTellerFunctor<Particle, false, true, autopas::FunctorN3Modes::Both, true>;
void setEquilateral(Particle &p2, Particle &p3, double distance) {
    p2.
    setR({0., distance, 0.});
    p3.setR({std::sqrt(3)/2 * distance, 0.5 * distance, 0.});
}

void setLinear(Particle &p1, Particle &p2, Particle &p3, double distance) {
    p1.setR({2.0 + distance, 0., 0.});
    p2.setR({2.0, 0., 0.});
    p3.setR({1, 0., 0.});
}

void setIsosceles(Particle& p1, Particle &p2, Particle &p3, double distance) {
    p1.setR({distance, 2.0, 0.0});
    p2.setR({0.0, 0.0, 0.0});
    p3.setR({0.0, 4.0, 0.0});}

/**
 * Mini validation setup that writes functor results to a CSV
 * @return
 */
int main(int argc, char* argv[]) {
//    auto lut = mdLib::TriwiseLUT(1000);
auto lut3B = mdLib::LUT3B(0,8);
//lut3B.setNu(1.6152500E-3);
//lut3B.setResolution(200);
//    lut3B.fill_plain(9);
lut3B.setResolution(100);
    lut3B.fill_plain_krypton(8);
//auto lut3B = mdLib::LUT3B(800,9);
//auto lut2B = mdLib::LUT2B(100,6.25);
//auto lut3B = mdLib::LUT3B(100,6.25);

//    const double cutoff{6.25};
//   Functor functor{3, &lut3B};
//   Functor  functor(cutoff, &lut2B);
//    Functor functor{cutoff, &lut};
    // set nu
//    functor.setParticleProperties(1.6152500E-3);
//    functor.setParticleProperties(24,1 );




//argon
//lut3B.fill_plain(9);
Functor functor{8, &lut3B};
//Functor functor{6.25, &lut2B};
//    functor.setParticleProperties(1,1 );
//Functor functor{25};
//functor.setParticleProperties();
//lut3B.fill_plain_krypton(1.44);



//    auto particle1 = Particle({0., 0., 0.}, {0., 0., 0.}, 0, 0);
//    auto particle3 = Particle({2,2, 0.}, {0., 0., 0.}, 0, 0);
//    auto particle2 = Particle({0.50, 0.50, 0.}, {0., 0., 0.}, 0, 0);

    auto particle1 = Particle({0., 0., 0.}, {0., 0., 0.}, 0, 0);
    auto particle2 = Particle({2.5,5, 0.}, {0., 0., 0.}, 0, 0);
    auto particle3 = Particle({1, 0.50, 0.}, {0., 0., 0.}, 0, 0);

   auto res1 = lut3B.getLUTValuesKrypton(3,2,1);
   auto cas1 = lut3B.getLUTValuesKrypton(2,0.5,0.1);
   auto tas1 = lut3B.getLUTValuesKrypton(100,50,10);
  auto res2 =  lut3B.getLUTValuesKrypton(3,1,2);
    auto cas2 = lut3B.getLUTValuesKrypton(2,0.1,0.5);
    auto tas2 = lut3B.getLUTValuesKrypton(100,10,50);


    auto res3 = lut3B.getLUTValuesKrypton(2,1,3);
    auto res4 =  lut3B.getLUTValuesKrypton(2,3,1);
    auto cas4 = lut3B.getLUTValuesKrypton(2,0.5, 0.1);
    auto tas4 = lut3B.getLUTValuesKrypton(100,50,10);

    auto res5 = lut3B.getLUTValuesKrypton(1,2,3);
    auto res6 =  lut3B.getLUTValuesKrypton(1,3,2);


//    setEquilateral(particle3, particle2, 5);
//    setEquilateral(particle2, particle3, 5);
setLinear(particle1, particle2,particle3,1);

//    setIsosceles( particle2,particle3,particle1, 5);
//    setIsosceles( particle2,particle1,particle3, 5);

//    auto particle2 = Particle({0., 0., 0.}, {0., 0., 0.}, 0, 0);
//    auto particle1 = Particle({0, 5, 0.}, {0., 0., 0.}, 0, 0);
//    auto particle3 = Particle({4.337, 2.5, 0.}, {0., 0., 0.}, 0, 0);

//    setEquilateral(particle1, particle3, 5);


////krypton
//    auto particle1 = Particle({0.0, 0.0, 0.0}, {0., 0., 0.}, 0, 0);
//    auto particle2 = Particle({0.4, 0.0, 0.0}, {0., 0., 0.}, 0, 0);
//    auto particle3 = Particle({0.2, 0.346, 0.0}, {0., 0., 0.}, 0, 0);
//    auto particle1= Particle({-0.1, 0.0, 0.0}, {0., 0., 0.}, 0, 0);
//    auto particle3 = Particle({0.4, 0.0, 0.0}, {0., 0., 0.}, 0, 0);
//    auto particle2 = Particle({0.2, 0.346, 0.0}, {0., 0., 0.}, 0, 0);



//    auto


    // Open the file
    std::filesystem::path cwd = std::filesystem::current_path();
    std::cout << "Current working directory: " << cwd << std::endl;
//    std::ofstream file("../AxilrodTeller/Particle1AlongX.csv");
//    std::ofstream file("3xa_global_lut_true_newton_false.csv");
//    std::ofstream file("personal-AT-global-lut-true-NN-newton-false.csv");
//    std::ofstream file("personal-AT-global-lut-true-NN-newton-true.csv");
//    std::ofstream file("personal-LJ-global-lut-true-LL-newton-true.csv");
//    std::ofstream file("personal-AT-globalnoLUT-lut-true-NN-newton-true.csv");
//    std::ofstream file("personal_AT_Lut_true_NN_globs_testingCOS_lut.csv");
//    std::ofstream file("LJ_LuT1000_true_globs_lut_rearranged.csv");
//    std::ofstream file("krypton_yesLUT_newtfalse_local_cutsqrt9.csv");
//    std::ofstream file("krypton_yes_LUT_smaller_distance_correction.csv");
//    std::ofstream file("krypton_yes_LUT_smaller_distance_4array_yesNewton.csv");
//    std::ofstream file("Krypton_IsoTest_minus_in_CW_LUT.csv");
//    std::ofstream file("newKrypton_ISOtest_noLUT.csv");
//    std::ofstream file("newKrypton_Isotest_p2p1p3_withLUT.csv");
//    std::ofstream file("standardTriangle_rotate_withLUT.csv");
    std::ofstream file("emptyTest.csv");
//    std::ofstream file("equiTest23.csv");
//    std::ofstream file("linear123.csv");

    // Check if the file is opened successfully
    if (!file.is_open()) {
        std::cout << "Error opening file!" << std::endl;
        return 1;
    }

    // Write headers
    file << "Xposition,Energy,Force_x,Force_y,Force_z" << std::endl;










    double epot{};
    size_t resolution = 50;
    double increment = 0.001; //original
//    double increment = 0.0001;
    double distance = 0.;
    for (auto i = 0; i < resolution; i++) {
        distance += increment;
        //move first particle along x axis
        particle1.setR({distance, 0., 0});
        // set Forces to 0
        particle1.setF({0.,0.,0.});
        particle2.setF({0.,0.,0.});
        particle3.setF({0.,0.,0.});
        // Compute the interactions
        functor.initTraversal();
        if(i == 487){
            std::cout << "here";
        }

        if(i == 555){
            std::cout << "here";
        }
        functor.AoSFunctor(particle1, particle2, particle3, newton3);
//        functor.AoSFunctor(particle1, particle2, newton3);
        functor.endTraversal(newton3);
        // Results
        auto force1 = particle1.getF();
        auto force2 = particle2.getF();
        auto force3 = particle3.getF();
        epot = functor.getPotentialEnergy();

        auto vir = functor.getVirial();

        // Write to CSV
        file << std::fixed << std::setprecision(15) << distance<< "," << vir << "," << epot << "," << force1[0] << "," << force2[0] << "," << force3[0]  << std::endl;
//        file << std::fixed << std::setprecision(15) << distance<< "," << vir << "," << epot << "," << force1[0] << "," << force1[1] <<"," << force1[2] <<  "," << force2[0] << "," << force2[1] << "," << force2[2] << "," << force3[0]   << "," << force3[1] << "," << force3[2] <<std::endl;
//        file << std::fixed << std::setprecision(20) << distance<< "," << vir << "," << epot << "," << force1[0] << "," << force2[0]   << std::endl;
//        file << std::fixed << std::setprecision(350) <<  epot  << std::endl;
    }

    file.close();

    std::cout<< "DONE";
}
