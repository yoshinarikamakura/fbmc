// Copyright (c) 2025 yoshinarikamakura
#include "FBMC.h"
#include <iostream>
using std::cerr;
using std::endl;

State FBMC::selectStateAfterImpactIonization(
            const State initial_state,
            const std::vector<double> dos_another_band,
            std::array<double, 2>& eh_pair_energies) {
    const double initial_energy = carrier_.getEnergy(initial_state);
    if (initial_energy < ETHII_) {
        cerr << "# Error in selectStateAfterImpactIonization(),\n";
        cerr << "#   ===> Initial energy is less than ionization threshold.\n";
        cerr << "# initial_energy = " << initial_energy << endl;
        exit(EXIT_FAILURE);
    }

    // 0 interacts with 2, then 1, 2, 3 are generated 
    double secondary_carrier_energy1; // Final state in the same band
    double secondary_carrier_energy2; // Final state in another band
    double secondary_carrier_energy3; // Final state in the same band

    if (dos_another_band.empty() == false) {
        // Calculate using random-k approximation
        std::uniform_real_distribution<double> urand(0.0, 1.0);
        const int SIZE_OF_DOS_TABLE_ANOTHER_BAND = dos_another_band.size();
        std::vector<double> cumulative2;
        cumulative2.resize(SIZE_OF_DOS_TABLE_ANOTHER_BAND);

        double previous2 = 0.0;
        for (int ie2 = 0; ie2 < SIZE_OF_DOS_TABLE_ANOTHER_BAND; ++ie2) {
            const double e2 = ie2 * DELTAE_DOS_;
            double integral = 0.0;
/*
            for (int ie1 = 0; ie1 < SIZE_OF_DOS_TABLE_; ++ie1) {
                const double e1 = ie1 * DELTAE_DOS_;
                const double e3 = initial_energy - e1 - e2 - EGAP_;
                const int ie3 = e3 / DELTAE_DOS_;
                if (ie3 >= 0 && ie3 < SIZE_OF_DOS_TABLE_) {
                    integral += dos_[ie3] * dos_[ie1];
                }
            }
*/
            for (int ie1 = 0; ie1 < SIZE_OF_DOS_TABLE_ - 1; ++ie1) {
                const double e1 = (ie1 + 0.5) * DELTAE_DOS_;
                const double e3 = initial_energy - e1 - e2 - EGAP_;
                const int ie3 = e3 / DELTAE_DOS_;
                if (ie3 >= 0 && ie3 < SIZE_OF_DOS_TABLE_) {
                    integral += dos_[ie3] * (dos_[ie1] + dos_[ie1 + 1])
                                * 0.5 * DELTAE_DOS_;
                }
            }
            cumulative2[ie2] = previous2 + dos_another_band[ie2] * integral;
            previous2 = cumulative2[ie2];
        } 

        double r2 = urand(mt_) * cumulative2[SIZE_OF_DOS_TABLE_ANOTHER_BAND - 1];
        int ie2;
        for (ie2 = 0; cumulative2[ie2] < r2; ie2++);
        secondary_carrier_energy2 = ie2 * DELTAE_DOS_;

        std::vector<double> cumulative1;
        cumulative1.resize(SIZE_OF_DOS_TABLE_);

        double previous1 = 0.0;
        for (int ie1 = 0; ie1 < SIZE_OF_DOS_TABLE_; ++ie1) {
            const double e1 = ie1 * DELTAE_DOS_;
            const double e3 = initial_energy - e1
                              - secondary_carrier_energy2 - EGAP_;
            const int ie3 = e3 / DELTAE_DOS_;
            if (ie3 >= 0 && ie3 < SIZE_OF_DOS_TABLE_) {
                cumulative1[ie1] = previous1 + dos_[ie3] * dos_[ie1]; 
            }
            else {
                cumulative1[ie1] = previous1; 
            }
            previous1 = cumulative1[ie1];
        } 

        double r1 = urand(mt_) * cumulative1[SIZE_OF_DOS_TABLE_ - 1];
        int ie1;
        for (ie1 = 0; cumulative1[ie1] < r1; ie1++);
        secondary_carrier_energy1 = ie1 * DELTAE_DOS_;
    }
    else {
        // Calculate deterministically
        // Secondary carrier with the same polarity
        secondary_carrier_energy1 = AII1_ * initial_energy + BII1_;

        // Secondary carrier with the opposite polarity
        secondary_carrier_energy2 = AII2_ * initial_energy + BII2_;
    }

    secondary_carrier_energy3 = initial_energy
                                - secondary_carrier_energy1
                                - secondary_carrier_energy2
                                - EGAP_;

    if (secondary_carrier_energy3 < 0.0) {
        cerr << "# Warning in selectStateAfterImpactIonization(),\n";
        cerr << "#   ===> Energy concervation is violated.\n";
        cerr << "# charge = " << initial_state.charge << endl;
        cerr << "# initial_energy = " << initial_energy << endl;
        // exit(EXIT_FAILURE);
        secondary_carrier_energy3 = THERMAL_ENERGY_;
    }

    State final_state =
        selectStateOnIsoEnergySurface(secondary_carrier_energy1,
                                      initial_state.r,
                                      initial_state.charge);
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



