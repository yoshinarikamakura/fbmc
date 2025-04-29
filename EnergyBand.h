#pragma once

#include "Definitions.h"
#include <array>
#include <random>
#include <string>
#include <vector>

class EnergyBand {

public:

/* Constructor
 *
 * === Input ===
 *  std::string band_filename ....Input file name
 */
    EnergyBand(std::string band_filename, std::mt19937& rng);

/* Energy of carrier
 *
 * === Input ===
 * State ....Wave vector (2*pi/a) and Band index
 *
 * === Output ===
 * double getEnergy ....Carrier energy (eV)
 */
    double getEnergy(State s);

/* Group vlocity of carrier
 *
 * === Input ===
 * State ....Wave vector (2*pi/a) and Band index
 *
 * === Output ===
 * Vector3 getVelocity ....Group Velocity (m/s)
*/
    Vector3 getVelocity(State s);

    inline double getNumberOfBands() const { return NB; }
    inline double getNumberOfTetrahedra() const { return NT; }
    inline double getUnitOfWaveVector() const { return K_UNIT; }
    inline double getFactorForDensityOfStates() const { return FACTOR_DOS; }
    inline double getMaximumEnergy() const { return emax; }
    inline std::vector<double> getMaximumEnergyPerBand() const { return emax_per_band; }
    inline Vector3 getMaximumWaveVector() const { return vkmax; }

    State getStateInTetrahedron(const double energy, const int nt, const int nb);
    double getTetrahedronDOS(const int nt, const int n, const double e);
    std::array<double, 4> getTetrahedronVertexEnergies(const int nt, const int nb);
    Vector3 getWaveVectorDifference(const State state, const int it);

private:
    std::mt19937& mt;
//    const std::string TABLE_DIRECTORY_NAME = "./TABLE/";
    double K_UNIT;
    double emax;
    std::vector<double> emax_per_band;
    Vector3 vkmax;
    double FACTOR_DOS;
    Vector3 FACTOR_GROUP_VELOCITY;

// Lattice Constant (m)
    double LATTICE_CONSTANT; 

// Number of devisions from Gamma-point to X-point in k-space
    int NK;

// Number of bands
    int NB;

// Domain of k-space
    Vector3 KMIN, KMAX;

// Number of triangles
    int NT;

// Number of grid points
    int NG;

// Table of grid points
    std::vector<int> grid_x;
    std::vector<int> grid_y;
    std::vector<int> grid_z;

// Table of grid number
    std::vector<std::vector<std::vector<int>>> grid_number;

// Table of tetrahedron vertex
    std::vector<int> tetrahedron_vertex1;
    std::vector<int> tetrahedron_vertex2;
    std::vector<int> tetrahedron_vertex3;
    std::vector<int> tetrahedron_vertex4;

// Table of cube
    std::vector<int> cube_x;
    std::vector<int> cube_y;
    std::vector<int> cube_z;

// Energy (eV)
    std::vector<std::vector<std::vector<std::vector<double>>>> grid_energy;

    double getDOS(const double energy);
    double loadBandFile(std::string band_filename);
};
