#include "FBMC.h"
#include <iostream>
using namespace std;

State FBMC::selectStateAfterImpactIonization(const State initial_state, array<double, 2>& eh_pair_energies) {
    const double initial_energy = carrier_.getEnergy(initial_state);
    double secondary_carrier_energy1, secondary_carrier_energy2, secondary_carrier_energy3;

    if (initial_energy > ETHII_) {
        // Secondary carrier with the same polarity
        secondary_carrier_energy1 = AII1_ * initial_energy + BII1_;

        // Secondary carrier with the opposite polarity
        secondary_carrier_energy2 = AII2_ * initial_energy + BII2_;

        // Secondary carrier with the same polarity
        secondary_carrier_energy3 = initial_energy - secondary_carrier_energy1 - secondary_carrier_energy2 - ETHII_;

	if (secondary_carrier_energy3 < 0.0) {
//            cerr << "# Error in selectStateAfterImpactIonization(),\n";
//            cerr << "#   ===> Energy concervation is violated.\n";
//            cerr << "# charge = " << initial_state.charge << endl;
//            cerr << "# initial_energy = " << initial_energy << endl;
//            exit(EXIT_FAILURE);
            secondary_carrier_energy3 = THERMAL_ENERGY_;
	}
    }
    else {
        cerr << "# Error in selectStateAfterImpactIonization(),\n";
        cerr << "#   ===> Initial energy is less than ionization threshold.\n";
        cerr << "# initial_energy = " << initial_energy << endl;
        exit(EXIT_FAILURE);
    }

    State final_state = selectStateOnIsoEnergySurface(secondary_carrier_energy1, initial_state.r, initial_state.charge);
    eh_pair_energies = {secondary_carrier_energy2, secondary_carrier_energy3};

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



