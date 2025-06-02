#include "FBMC.h"
#include <iostream>
using namespace std;

State FBMC::selectStateAfterImpactIonization(const State initial_state) {
    const double initial_energy = carrier_.getEnergy(initial_state); // Initial State Energy (eV)
    double secondary_carrier_energy;

    if (initial_energy > ETHII_) {
        secondary_carrier_energy = AII1_ * initial_energy + BII1_;
    }
    else {
        cerr << "# Error in selectStateAfterImpactIonization(),\n";
        cerr << "#   ===> Initial energy is less than ionization threshold.\n";
        cerr << "# initial_energy = " << initial_energy << endl;
        exit(EXIT_FAILURE);
    }
    State final_state = selectStateOnIsoEnergySurface(secondary_carrier_energy);

    final_state.r = initial_state.r;

    return final_state;
}

double FBMC::getImpactIonizationRate(const State initial_state) {
    const double initial_energy = carrier_.getEnergy(initial_state);

    double rate = 0.0;
    if (initial_energy > ETHII_) {
        rate = SII_ * pow(initial_energy - ETHII_, POWII_);
    }

    return rate;
}
