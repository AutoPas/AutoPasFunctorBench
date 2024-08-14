#include <iostream>

#include <molecularDynamicsLibrary/ArgonInclude/RepulsiveTerm.h>
#include <molecularDynamicsLibrary/ArgonInclude/DispersionTerm.h>
#include <molecularDynamicsLibrary/ArgonInclude/CosineHandle.h>

#include <autopas/cells/FullParticleCell.h>
#include <autopas/utils/ArrayMath.h>
#include <fstream>

//SI units (nm and K)
static constexpr std::array<double, 23> A{
        {-5.39365445993314E+05,  // 000
         -1.00043526760807E+07,  // 001
         -1.80479894697093E+10,  // 011
         2.68023739515408E+07,   // 111
         5.17630401857978E+11,   // 002
         1.20250233629457E+07,   // 012
         -6.86507513446023E+07,  // 112
         7.73641060191982E+07,   // 022
         4.07116202374599E+07,   // 122
         2.57546504143754E+07,   // 222
         1.29463884038186E+07,   // 003
         -3.08989961490876E+11,  // 013
         3.29616043775900E+11,   // 113
         -1.21201021419532E+07,  // 023
         4.54508019194995E+07,   // 123
         3.22601026022283E+08,   // 033
         -1.79863484497154E+07,  // 004
         5.63204555102674E+08,   // 014
         7.64813924806795E+06,   // 114
         -8.82961781148373E+05,  // 024
         -1.02496007862500E+07,  // 005
         -3.04174890291515E+04,  // 015
         -2.83863618111236E+10}  // 006
};

static constexpr std::array<double, 23> alpha{
        {8.09052299484753E+00,  // 000
         9.52298731190775E+00,  // 001
         1.97867044131258E+01,  // 011
         8.63168953894591E+00,  // 111
         2.47643526123088E+01,  // 002
         8.39137345537347E+00,  // 012
         9.07955833453440E+00,  // 112
         1.39334614374608E+01,  // 022
         9.37640048180289E+00,  // 122
         8.01934231300047E+00,  // 222
         8.10590814604500E+00,  // 003
         2.22948530135448E+01,  // 013
         2.25887370431209E+01,  // 113
         7.87549358334693E+00,  // 023
         9.58307979519088E+00,  // 123
         1.44441527110870E+01,  // 033
         7.98634790509653E+00,  // 004
         1.43154883935442E+01,  // 014
         9.12235520967071E+00,  // 114
         7.92438461086463E+00,  // 024
         1.20141476840267E+01,  // 005
         7.09781720339141E+00,  // 015
         2.46296194255218E+01}  // 006
};

static constexpr std::array<double, 5> Z{
        {2.81011266190959E-04,   // 111
         -6.14241347348619E-05,  // 112
         8.72021611550457E-06,   // 122
         -4.82191783511564E-08,  // 222
         2.93702828075611E-06}   // 113
};

static constexpr std::array<double, 5> beta{
        {3.99870891182023E+02,  // 111
         2.82746852049202E+01,  // 112
         2.49749116804324E+01,  // 122
         3.93440001759947E+01,  // 222
         3.39906094408459E+01}  // 113
};

// type aliases for ease of use
    using DisplacementHandle = autopas::utils::ArrayMath::Argon::DisplacementHandle;

    namespace autopas::utils::ArrayMath {
        std::array<double, 3> rescale_array(std::array<double, 3> array, double s) {
            return array * s;
        }

        std::array<double, 3> sum_arrays(std::array<double, 3> array1, std::array<double, 3> array2) {
            return array1 + array2;
        }
    }

/**
 * Mini validation setup that writes functor results to a CSV
 * @return*
 **/
int main(int argc, char *argv[]) {

    enum ID {
        I, J, K
    };

    // particle 1 moving along x
    std::array<double, 3> positionI{};
    std::array<double, 3> positionJ{{1, 1, 0}};
    std::array<double, 3> positionK{{1, -1, 0}};

    // particle 1 moving along y
    /*std::array<double, 3> positionI{};
    std::array<double, 3> positionJ{{-1, 1, 0}};
    std::array<double, 3> positionK{{1, 1, 0}};*/

    // particle 1 moving along z
    /*std::array<double, 3> positionI{};
    std::array<double, 3> positionJ{{0, 1, 1}};
    std::array<double, 3> positionK{{0, -1, 1}};*/

    // equilateral triangle
    /*std::array<double, 3> positionI_unit{{0, 0, 0}};
    std::array<double, 3> positionJ_unit{{0, 1, 0}};
    std::array<double, 3> positionK_unit{{std::sqrt(3) / 2, 0.5, 0}};*/

    // symmetric linear geometry
    /*std::array<double, 3> positionI_unit{{0, 0, 0}};
    std::array<double, 3> positionJ_unit{{1, 0, 0}};
    std::array<double, 3> positionK_unit{{2, 0, 0}};*/

    // Open the file
    std::ofstream file("/Users/irene/TUM/thesis/AutoPas_Thesis/AutoPasFunctorBench/Argon/Particle1AlongX.csv");

    // Check if the file is opened successfully
    if (!file.is_open()) {
        std::cout << "Error opening file!" << std::endl;
        return 1;
    }

    // Write headers
    file << "Position,Energy,Force_x,Force_y,Force_z" << std::endl;


    size_t resolution = 2000;
    double increment = 0.001;
    double distance = 0.;
    for (auto i = 0; i < resolution; i++) {
        distance += increment;
        positionI = std::array<double, 3>{{distance, 0, 0}};

        // resize the triangle
        /*auto positionI = autopas::utils::ArrayMath::rescale_array(positionI_unit, r);
        auto positionJ = autopas::utils::ArrayMath::rescale_array(positionJ_unit, r);
        auto positionK = autopas::utils::ArrayMath::rescale_array(positionK_unit, r);*/

        const auto displacementHandleIJ = DisplacementHandle(positionI, positionJ, I, J);
        const auto displacementHandleJK = DisplacementHandle(positionJ, positionK, J, K);
        const auto displacementHandleKI = DisplacementHandle(positionK, positionI, K, I);

        const auto dispersionPotential = U_dispersive(Z, beta, displacementHandleIJ, displacementHandleJK,
                                                      displacementHandleKI);
        const auto repulsivePotential = U_repulsive(A, alpha, displacementHandleIJ, displacementHandleJK,
                                                    displacementHandleKI);
        const auto totalPotential = dispersionPotential + repulsivePotential;

        // Force on particle I
        const auto dispersiveForceI = F_dispersive<I>(Z, beta, displacementHandleIJ, displacementHandleJK,
                                                      displacementHandleKI);
        const auto repulsiveForceI = F_repulsive<I>(A, alpha, displacementHandleIJ, displacementHandleJK,
                                                    displacementHandleKI);
        const auto totalForceI = autopas::utils::ArrayMath::sum_arrays(dispersiveForceI, repulsiveForceI);

        const auto f = totalPotential;
        const auto nabla_f = totalForceI;
        // Write to CSV
        file << std::fixed << std::setprecision(15) << distance << "," << f << ","
        << nabla_f[0] << "," << nabla_f[1] << "," << nabla_f[2] << ","
        << std::endl;
    }

    file.close();
    return 0;
}

