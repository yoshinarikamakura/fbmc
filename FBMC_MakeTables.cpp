#include "FBMC.h"
#include <fstream>
#include <iostream>
using namespace std;

static int NE;

void FBMC::constructDOSTables(string band_filename) {
// Make a table of maximum DOS at a given energy
    const double emax = carrier_.getMaximumEnergy();
    NE = static_cast<int>(emax / DELTAE_DOSMAX_);
    cout << "# Size of dosmax table = " << NE << endl;
    dosmax_.resize(NE);
    loadDOSmax(band_filename);
    cout << "# dosmax_ table has been loaded.\n";

// Make a list of tetrahedron numbers containing a given energy
    const int NB = carrier_.getNumberOfBands();
    const int NT = carrier_.getNumberOfTetrahedra();
    int ie1_max = 0;
    int ie41_max = 0;
    for (int ib = 0; ib < NB; ++ib) {
    for (int it = 0; it < NT; ++it) {
	array<double, 4> e = carrier_.getTetrahedronVertexEnergies(it, ib);
        const int ie1 = static_cast<int>(e[0] / DELTAE_MINMAX_);
        const int ie41 = static_cast<int>(e[3] / DELTAE_MINMAX_) - ie1 + 1;
        if (ie1 > ie1_max) ie1_max = ie1;
        if (ie41 > ie41_max) ie41_max = ie41;
    }
    }

// Size of minmax table
    NE1_MINMAX_ = ie1_max + 1;
    NE41_MINMAX_ = ie41_max + 1;

// Dynamic memory allocation of 2D-array [NE1_MINMAX][NE41_MINMAX] of list
    minmax_.resize(NE1_MINMAX_);
    for (int ie1 = 0; ie1 < NE1_MINMAX_; ++ie1) {
        minmax_[ie1].resize(NE41_MINMAX_);
    }

    for (int ib = 0; ib < NB; ++ib) {
    for (int it = 0; it < NT; ++it) {
	array<double, 4> e = carrier_.getTetrahedronVertexEnergies(it, ib);
        const int ie1 = static_cast<int>(e[0] / DELTAE_MINMAX_);
        const int ie41 = static_cast<int>(e[3] / DELTAE_MINMAX_) - ie1 + 1;
        if (ie1 >= NE1_MINMAX_) {
            cerr << "# Error in constructTable(),\n";
            cerr << "#   ===> Use larger NE1_MINMAX.\n";
            exit(EXIT_FAILURE);
        }
        if (ie41 >= NE41_MINMAX_) {
            cerr << "# Error in constructTable(),\n";
            cerr << "#   ===> Use larger NE41_MINMAX.\n";
            exit(EXIT_FAILURE);
        }
        minmax_[ie1][ie41].push_back(make_pair(it, ib));
    }
    } 
    cout << "# minimax table has been created.\n";
}

void FBMC::loadDOSmax(std::string band_filename) {
    string dosmax_filename = TABLE_DIRECTORY_NAME_ + "DOSmax_" + band_filename;
    ifstream fin_dosmax(dosmax_filename, ios::in);

    if (fin_dosmax.is_open()) {
        cout << "# File: " + dosmax_filename + " is found.\n";
        for (int ie = 0; ie < NE; ++ie) {
            int ie_scan;
	    double tmp;
            fin_dosmax >> ie_scan >> tmp;
	    if (ie_scan == ie) {
                dosmax_[ie] = tmp;
	    }
	    else {
                cout << "ERROR in making dosmax table.\n";
                exit(EXIT_FAILURE);
	    }
        }
        fin_dosmax.close();
    }
    else {
        cout << "# Not found dosmax table.\n";
        cout << "# Now, making dosmax table. Please wait for a while.\n";

        makeDOSmaxTable();
        ofstream fout_dosmax(dosmax_filename, ios::out);
        for (int ie = 0; ie < NE; ++ie) {
            fout_dosmax << ie << ' ' << dosmax_[ie] << endl;
        }
        fout_dosmax.close();
    }
}

void FBMC::makeDOSmaxTable(void) {

    const int NB = carrier_.getNumberOfBands();
    const int NT = carrier_.getNumberOfTetrahedra();

    for (int ie = 0; ie < NE; ++ie) {
        const double energy_min = ie * DELTAE_DOSMAX_;
        const double energy_max = energy_min + DELTAE_DOSMAX_;

        double dosmax_tmp = 0.0;
        for (int ib = 0; ib < NB; ++ib) {
        for (int it = 0; it < NT; ++it) {
            array<double, 4> e = carrier_.getTetrahedronVertexEnergies(it, ib);

            if (e[0] < energy_max && e[3] > energy_min) {
                const double e_scan_min = (energy_min > e[0]) ? energy_min : e[0];
                const double e_scan_max = (energy_max < e[3]) ? energy_max : e[3];
                if (e_scan_max - e_scan_min < 1.0e-9) {
                    // Energy is constant in this tetrahedron, and thus dosmax = 0
                    break;
		}

		const double delta_e = (e_scan_max - e_scan_min) * EPS_ESCAN_;
                for (double e_scan = e_scan_min; e_scan < e_scan_max; e_scan += delta_e) {
                     const double dosmax_interval = carrier_.getTetrahedronDOS(it, ib, e_scan);
                    // Update dosmax_tmp
                    if (dosmax_interval > dosmax_tmp) {
                            dosmax_tmp = dosmax_interval;
                    }
                }
            }
        }
        }
        dosmax_[ie] = dosmax_tmp;
    }
}
