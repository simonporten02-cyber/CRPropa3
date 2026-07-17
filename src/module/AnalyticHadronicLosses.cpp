#include "crpropa/module/AnalyticHadronicLosses.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"
#include <fstream>
#include <limits>
#include <stdexcept>
#include <cmath>

namespace crpropa {

    AnalyticHadronicLosses::AnalyticHadronicLosses(ref_ptr<Density> dens, double limit) :
        Module(), density(dens), limit(limit) {
        setDescription("AnalyticHadronicLosses");
    }

    AnalyticHadronicLosses::Table AnalyticHadronicLosses::loadTableFile(const std::string& absPath){
        Table T;
        std::ifstream in(absPath.c_str());
        if (!in.good()) throw std::runtime_error("HadronicLosses: could not open "+absPath);
        while (in.good()){
            if (in.peek()!='#'){
                double lgE, sig;
                in >> lgE >> sig;
                if (in){
                    T.logEGeV.push_back(lgE);
                    T.sigma_mb.push_back(sig);
                }
            }
            in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        in.close();
        if (!T.valid()) throw std::runtime_error("HadronicLosses: bad/short table "+absPath);
        return T;
    }

    void AnalyticHadronicLosses::setSigmaTableForPDG(int pdgAbs, const std::string& path){
        tableByPDG[std::abs(pdgAbs)] = loadTableFile(getDataPath(path));
    }

    void AnalyticHadronicLosses::process(Candidate* c) const {
        int a = std::abs(c->current.getId());
        if (!isPionOrKaon(a)) return;

        const Table* T = resolveTable(a);
        if (!T) return; // no pion/kaon table loaded

        // local frame step (same convention as Synchrotron)
        double z = c->getRedshift();
        double stepLocal = c->getCurrentStep() / (1 + z);
        if (stepLocal <= 0) return;

        // target density (m^-3) from your Density module
        Vector3d pos = c -> current.getPosition();
        const double nTarget = getDensityAtPosition(pos);
        if (nTarget <= 0) return;

        // interpolate sigma(E) and kappa(E) from tables (x = log10(E/GeV))
        double E  = c->current.getEnergy(); // J
        double lgE = std::log10(std::max(E / GeV, 1e-30));
        double sigma = interpolate(lgE, T->logEGeV, T->sigma_mb) * milli * barn; // m^2
        double kappa = 0.5; // normal assumption for high-energy collisions

        if (sigma <= 0.0) return;

        // continuous loss: dE/dx = n * sigma_inel * κ * E
        double dEdx = nTarget * sigma * kappa * E; // J/m
        double dE   = dEdx * stepLocal;

        double Enew = E - dE;
        c->current.setEnergy(Enew);

        if (dEdx > 0.0)
            c->limitNextStep(limit * E / dEdx);
    }


    const AnalyticHadronicLosses::Table* AnalyticHadronicLosses::resolveTable(int a) const{
        auto it = tableByPDG.find(a);
        if (it!=tableByPDG.end() && it->second.valid()) 
            return &it->second;
        return nullptr;
    }

    bool AnalyticHadronicLosses::isPionOrKaon(int a){
        return (a==211 || a==321 || a==311 || a==130 || a==310);
    }

    void AnalyticHadronicLosses::setLimit(double v){ 
        limit = v; 
    }
    double AnalyticHadronicLosses::getLimit() const { 
        return limit; 
    }

    std::string AnalyticHadronicLosses::getDescription() const {
        std::stringstream s;
        s << "π/K hadronic cooling (loss-only), density=";
        if (density) s << density->getDescription();
        else          s << "<none>";
        s << ", component=";
        s << ", step limiter=" << limit;
        return s.str();
    }

    double AnalyticHadronicLosses::getDensityAtPosition(Vector3d &pos) const {
        return density -> getDensity(pos);
    }


} // namespace crpropa