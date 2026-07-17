#include "crpropa/ParticleMass.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Common.h"
#include "crpropa/Units.h"

#include "kiss/convert.h"
#include "kiss/logger.h"

#include <vector>
#include <fstream>
#include <stdexcept>
#include <limits>

namespace crpropa {
struct PDGMassTable {
	bool initialized;
	std::vector<std::pair<int, double>> table;  // (abs(PDG), mass_kg)

	PDGMassTable() : initialized(false) {}

	void init() {
		std::string filename = getDataPath("particle_mass.txt");
		std::ifstream infile(filename.c_str());
		if (!infile.good())
			throw std::runtime_error("crpropa: could not open file " + filename);

		// Format per line: <pdgId> <mass_kg>
		// Comment lines start with '#'
		while (infile.good()) {
			if (infile.peek() != '#') {
				int pid;
				double mass_kg;
				if (infile >> pid >> mass_kg) {
					table.push_back(std::make_pair(std::abs(pid), mass_kg));
				}
			}
			infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}
		infile.close();
		initialized = true;
	}

	bool tryGetMass(int absPid, double &mass_kg) {
		if (!initialized) {
#pragma omp critical(init_pdgmass)
			{
				if (!initialized)
					init();
			}
		}
		for (size_t i = 0; i < table.size(); ++i) {
			if (table[i].first == absPid) {
				mass_kg = table[i].second;
				return true;
			}
		}
		return false;
	}
};

static PDGMassTable pdgMassTable;

struct NuclearMassTable {
	bool initialized;
	std::vector<double> table;

	NuclearMassTable() {
		initialized = false;
	}

	void init() {
		std::string filename = getDataPath("nuclear_mass.txt");
		std::ifstream infile(filename.c_str());

		if (!infile.good())
			throw std::runtime_error("crpropa: could not open file " + filename);

		table.assign((NUCLEAR_ZMAX + 1) * NUCLEAR_NSTRIDE, 0.0);

		int Z, N;
		double mass;
		while (infile.good()) {
			if (infile.peek() != '#') {
				infile >> Z >> N >> mass;
				if (Z <= NUCLEAR_ZMAX && N <= NUCLEAR_NMAX)
					table[Z * NUCLEAR_NSTRIDE + N] = mass;
			}
			infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}

		infile.close();
		initialized = true;
	}

	double getMass(std::size_t idx) {
		if (!initialized) {
#pragma omp critical(init)
			init();
		}
		if (table[idx] == 0.0)
			return 0.0; // triggers approximation in nuclearMass()
		return table[idx];
	}
};

static NuclearMassTable nuclearMassTable;

double particleMass(int id) {
	if (isNucleus(id))
		return nuclearMass(id);

	// prefer data-file lookup for everything non-nuclear
	double mass;
    const int aId = std::abs(id);
    if (pdgMassTable.tryGetMass(aId, mass))
        return mass;

	if (abs(id) == 11)
		return mass_electron;
	return 0.0;
}

double nuclearMass(int id) {
	int A = massNumber(id);
	int Z = chargeNumber(id);
	return nuclearMass(A, Z);
}

double nuclearMass(int A, int Z) {
	int N = A - Z;
	if ((A < 1) or (Z < 0) or (Z > A) or (Z > NUCLEAR_ZMAX) or (N > NUCLEAR_NMAX)) {
		KISS_LOG_WARNING <<
		"nuclearMass: nuclear mass not found in the mass table for " <<
	        "A = " << A << ", Z = " << Z << ". " <<
		"Approximated value used A * amu - Z * m_e instead.";
		return A * amu - Z * mass_electron;
	}
	double m = nuclearMassTable.getMass(Z * NUCLEAR_NSTRIDE + N);
	if (m == 0.0) {
		KISS_LOG_WARNING <<
		"nuclearMass: nuclear mass not found in the mass table for " <<
	        "A = " << A << ", Z = " << Z << ". " <<
		"Approximated value used A * amu - Z * m_e instead.";
		return A * amu - Z * mass_electron;
	}
	return m;
}

} // namespace crpropa
