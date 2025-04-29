#pragma once
#include <list>
#include "Definitions.h"
#include "EnergyBand.h"

class FBMC {

public:
// ====== Public interface ======
 
//
// Constructor
//
// === Input ===
//  char* filename ....Input file name
//
    FBMC(std::string param_filename,
          std::string band_filename,
	  std::string phonon_filename,
	  double temperature,
	  std::mt19937& rng);

//
// Energy from wave vector and band index
// 
// === Input ===
// State ....Wave vector (2*pi/a) and Band index
//
// === Output ===
// double getEnergy ....Carrier energy (eV)
//
    inline double getEnergy(const State s) { return carrier_.getEnergy(s); }

//
// Group vlocity from wave vector and band index
//
// === Input ===
// State ....Wave vector (2*pi/a) and Band index
//
// === Output ===
// Vector3 getVelocity ....Group Velocity (m/s)
//
    inline Vector3 getVelocity(const State s) { return carrier_.getVelocity(s); }

//
// Select a carrier state randomly on an iso-energy surface
//
// === Input ===
// double energy ....Carrier energy (eV)
//
//  === Output ===
// State selectStateOnIsoEnergySurface ....Carrier state
//
    State selectStateOnIsoEnergySurface(const double energy);

//
// Carrier free flight
//
// === Input ===
// double dt ....Time step (s)
// Vector3 field ....Electric field (V/m)
// State initial_state ....Carrier state before free flight
//
//  === Output ===
// State flightFree ....Carrier state after free flight
//
    State flightFree(const double dt, const Vector3 field, const State initial_state);

//
// Carrier scattering
//
// === Input ===
// double dt ....Time step (s)
// State initial_state ....Carrier state before scattering
//
//  === Output ===
// State scatter ....Carrier state after scattering
// int scattering_mechanism ....Scattering mechanism
//
    State scatter(const double dt, const State initial_state, int& scattering_mechanism);



private:
// ====== Internal Data (Preset parameters) ======

// Directory name of input files
    const std::string TABLE_DIRECTORY_NAME_ = "./inputs/";

// Energy interval for DOSMAX table (eV)
    static constexpr double DELTAE_DOSMAX_ = 0.001;

// (eV)
    static constexpr double EPS_ESCAN_ = 0.01;

// Energy interval for MINMAX table (eV)
    static constexpr double DELTAE_MINMAX_ = 0.01;

// Threshold of infinite roop
    static constexpr int INFINITE_LOOP_THRESHOLD_ = 1000000;

// Minimum energy for scattering rate table (eV);
    static constexpr double MINIMUM_ENERGY_FOR_SCATTERING_RATE_TABLE_ = 0.001;

// Size of scattering rate table
    static constexpr int SIZE_OF_SCATTERING_RATE_TABLE_ = 1000;

// ====== Internal Data (Material parameters from input file) ======
 
// Crystal density (kg/m/m/m)
    double DENSITY_;

// Deformation potential parameter for acoustic phonon (?)
    double ADAC_;

// Deformation potential parameter for acoustic phonon (?)
    double BDAC_;

// Deformation potential parameter for optical phonon (?)
    double ADKOP_;

// Deformation potential parameter for optical phonon (?)
    double BDKOP_;

// Threshold energy of impact ionization (eV)
    double ETHII_;

// Prefactor of impact ionization rate (1/s)
    double SII_;

// Power in Keldysh-like formula of impact ionization ()
    double POWII_;

// Parameter to calculate secondary carrier energy
    double AII1_;

// Parameter to calculate secondary carrier energy
    double BII1_;

// Parameter to calculate secondary carrier energy
    double AII2_;

// Parameter to calculate secondary carrier energy
    double BII2_;

// ====== Internal Data (Objects) ======

    std::mt19937& mt_;
    EnergyBand carrier_;
    EnergyBand phonon_;

// ====== Internal Data ======

// Table of DOSmax
    std::vector<double> dosmax_;

// Size of minmax table
    int NE1_MINMAX_, NE41_MINMAX_;

// Minmax table
    std::vector<std::vector<std::list<std::pair<int, int>>>> minmax_;

// Factor used in carrier free flight
    double FACTOR_FREE_FLIGHT_;

// ====== Internal methods ======
    void constructDOSTables(std::string band_filename);
    void loadDOSmax(std::string band_filename);
    void makeDOSmaxTable(void);

    void loadMaterialParameter(std::string param_filename);

    void constructPhononScatteringTables(std::string param_filename, std::string band_filename, double temperature);

public:
    double getPhononScatteringRate(const int eta, const bool isAbsorption, const State initial_state);
private:

    State selectStateAfterPhononScattering(const int eta, const bool isAbsorption, const State initial_state);

    void loadPhononScatteringRate(std::string band_filename);
    void makePhononScatteringRateTable(void);
    double calculatePhononScatteringRate(const int eta, const bool isAbsorption, const State initial_state);
    double getToyMatrixElementForSilicon(State s, const bool isAbsorption, const double ei);
    std::vector<std::vector<double>> getMaxToyMatrixElement(const double eimax);

    double getImpactIonizationRate(const State initial_state);
    State selectStateAfterImpactIonization(const State initial_state);
};
