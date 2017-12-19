#ifndef SIMPLEGENERATOR_HH_
#define SIMPLEGENERATOR_HH_

struct SimpleEvent;

// Simple generator of an arbitrary number of independent jets
// with flat spectra in Pt, eta, and phi
class SimpleGenerator
{
public:
    // The constructor arguments are as follows:
    //
    //  njets          -- The number of jets to generate
    //
    //  partonCode     -- PYTHIA particle code of the parton
    //                    originating the jets. All jets
    //                    will be made with this parton type.
    //
    //  ptMin, ptMax   -- The minimum and the maximum transverse
    //                    momentum of the partons. The transverse
    //                    momentum spectrum of the partons
    //                    will be flat between these two values
    //
    //  etaMin, etaMax -- The minimum and the maximum pseudorapidity
    //                    of the partons
    //
    //  phiMin, phiMax -- The minimum and the maximum azimuthal angle
    //
    //  useOriginalSingleParticleGun  -- If "true", the code will use
    //                                   the original PYTHIA single jet
    //                                   gun to make jets. If "false",
    //                                   the jets will be made using
    //                                   one half of a dijet event.
    //
    // generateFlatInLogPt -- If "true", log(Pt) will have flat spectrum
    //                        between log(ptMin) and log(ptMax); otherwise
    //                        Pt will have flat spectrum between ptMin and
    //                        ptMax
    //
    SimpleGenerator(unsigned njets, unsigned partonCode,
                    double ptMin, double ptMax,
                    double etaMin, double etaMax,
                    double phiMin, double phiMax,
                    bool useOriginalSingleParticleGun=true,
                    bool generateFlatInLogPt=false);
    inline ~SimpleGenerator() {}

    // The following function populates the SimpleEvent with jets
    void generate(SimpleEvent* ev) const;

private:
    SimpleGenerator();

    unsigned njets_;
    unsigned partonCode_;
    double ptMin_;
    double ptMax_;
    double etaMin_;
    double etaMax_;
    double phiMin_;
    double phiMax_;
    double logMin_;
    double logMax_;
    bool single_;
    bool flatlog_;
};

#endif // SIMPLEGENERATOR_HH_
