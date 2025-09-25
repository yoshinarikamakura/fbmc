// Copyright (c) 2025 yoshinarikamakura
#include "FBMC.h"
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
using std::cerr;
using std::cout;
using std::endl;

void FBMC::constructDOSTables(std::string band_filename) {
    if (carrier_.getMinimumEnergy() < -DELTAE_DOSMAX_) {
        cerr << "# Error in constructDOSTables(),\n";
        cerr << "#   ===> Negative band energy is in: "
             << band_filename << endl;
        exit(EXIT_FAILURE);
    }

    const double emax = carrier_.getMaximumEnergy();

// Make a table of DOS
    SIZE_OF_DOS_TABLE_ = static_cast<int>(emax / DELTAE_DOS_);
    cout << "# Size of dos table = " << SIZE_OF_DOS_TABLE_ << endl;
    dos_.resize(SIZE_OF_DOS_TABLE_);
    loadDOS(band_filename);
    cout << "# dos_ table has been loaded.\n";

// Make a table of maximum DOS at a given energy
    SIZE_OF_DOSMAX_TABLE_ = static_cast<int>(emax / DELTAE_DOSMAX_);
    cout << "# Size of dosmax table = " << SIZE_OF_DOSMAX_TABLE_ << endl;
    dosmax_.resize(SIZE_OF_DOSMAX_TABLE_);
    loadDOSmax(band_filename);
    cout << "# dosmax_ table has been loaded.\n";

// Make a list of tetrahedron numbers containing a given energy
    const int NB = carrier_.getNumberOfBands();
    const int NT = carrier_.getNumberOfTetrahedra();
    int ie1_max = 0;
    int ie41_max = 0;
    for (int ib = 0; ib < NB; ++ib) {
    for (int it = 0; it < NT; ++it) {
        std::array<double, 4> e = carrier_.getTetrahedronVertexEnergies(it, ib);
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
        std::array<double, 4> e = carrier_.getTetrahedronVertexEnergies(it, ib);
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
        minmax_[ie1][ie41].push_back(std::make_pair(it, ib));
    }
    }
    cout << "# minimax table has been created.\n";
}

void FBMC::loadDOS(std::string band_filename) {
    std::string dos_filename = TABLE_DIRECTORY_NAME_
                                  + "DOS_" + band_filename;
    std::ifstream fin_dos(dos_filename, std::ios::in);

    if (fin_dos.is_open()) {
        cout << "# File: " + dos_filename + " is found.\n";
        for (int ie = 0; ie < SIZE_OF_DOS_TABLE_; ++ie) {
            int ie_scan;
            double tmp;
            fin_dos >> ie_scan >> tmp;
            if (ie_scan == ie) {
                dos_[ie] = tmp;
            } else {
                cout << "ERROR in making dos table.\n";
                exit(EXIT_FAILURE);
            }
        }
        fin_dos.close();
    } else {
        cout << "# Not found dos table.\n";
        cout << "# Now, making dos table. Please wait for a while.\n";

        makeDOSTable();
        std::ofstream fout_dos(dos_filename, std::ios::out);
        for (int ie = 0; ie < SIZE_OF_DOS_TABLE_; ++ie) {
            fout_dos << ie << ' ' << dos_[ie] << endl;
        }
        fout_dos.close();
    }
}

void FBMC::loadDOSmax(std::string band_filename) {
    std::string dosmax_filename = TABLE_DIRECTORY_NAME_
                                  + "DOSmax_" + band_filename;
    std::ifstream fin_dosmax(dosmax_filename, std::ios::in);

    if (fin_dosmax.is_open()) {
        cout << "# File: " + dosmax_filename + " is found.\n";
        for (int ie = 0; ie < SIZE_OF_DOSMAX_TABLE_; ++ie) {
            int ie_scan;
            double tmp;
            fin_dosmax >> ie_scan >> tmp;
            if (ie_scan == ie) {
                dosmax_[ie] = tmp;
            } else {
                cout << "ERROR in making dosmax table.\n";
                exit(EXIT_FAILURE);
            }
        }
        fin_dosmax.close();
    } else {
        cout << "# Not found dosmax table.\n";
        cout << "# Now, making dosmax table. Please wait for a while.\n";

        makeDOSmaxTable();
        std::ofstream fout_dosmax(dosmax_filename, std::ios::out);
        for (int ie = 0; ie < SIZE_OF_DOSMAX_TABLE_; ++ie) {
            fout_dosmax << ie << ' ' << dosmax_[ie] << endl;
        }
        fout_dosmax.close();
    }
}

void FBMC::makeDOSTable(void) {
    for (int ie = 0; ie < SIZE_OF_DOS_TABLE_; ++ie) {
        const double energy = ie * DELTAE_DOS_;
        dos_[ie] = carrier_.getDOS(energy);
    }
}

void FBMC::makeDOSmaxTable(void) {
    const int NB = carrier_.getNumberOfBands();
    const int NT = carrier_.getNumberOfTetrahedra();

    for (int ie = 0; ie < SIZE_OF_DOSMAX_TABLE_; ++ie) {
        const double energy_min = ie * DELTAE_DOSMAX_;
        const double energy_max = energy_min + DELTAE_DOSMAX_;

        double dosmax_tmp = 0.0;
        for (int ib = 0; ib < NB; ++ib) {
        for (int it = 0; it < NT; ++it) {
            std::array<double, 4> e
                 = carrier_.getTetrahedronVertexEnergies(it, ib);

            if (e[0] < energy_max && e[3] > energy_min) {
                const double e_scan_min =
                     (energy_min > e[0]) ? energy_min : e[0];
                const double e_scan_max =
                    (energy_max < e[3]) ? energy_max : e[3];
                if (e_scan_max - e_scan_min < 1.0e-9) {
                    // Energy is almost constant in this tetrahedron,
                    // and thus dosmax = 0
                    break;
                }

                const double delta_e = (e_scan_max - e_scan_min) * EPS_ESCAN_;
                for (double e_scan = e_scan_min;
                     e_scan < e_scan_max;
                     e_scan += delta_e) {
                     const double dosmax_interval =
                         carrier_.getTetrahedronDOS(it, ib, e_scan);
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
