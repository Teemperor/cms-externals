#include "SimpleEvent.hh"
#include "discretizeMCEvent.hh"

void discretizeMCEvent(const SimpleEvent& mcEvent, fftjet::Grid2d<Real>* calo,
                       const double magneticFieldB, const double caloRadius,
                       const double caloNoise, const double caloThreshold)
{
    calo->reset();

    // Cycle over jets within the MC event
    const unsigned njets = mcEvent.size();
    for (unsigned jetNumber=0; jetNumber<njets; ++jetNumber)
    {
        const SimulatedJet& jet(mcEvent[jetNumber]);

        // Cycle over particles within the jet
        const unsigned nParticles = jet.jetParticles.size();
        for (unsigned i=0; i<nParticles; ++i)
        {
            const PythiaJetGun::Particle& particle(jet.jetParticles[i]);

            // Check whether this particle reaches the calorimeter
            // and find the phi angle where it hits
            double caloPhi;
            if (phiAtRadius(particle, magneticFieldB, caloRadius, &caloPhi))
                // This particle did indeed reach the calorimeter.
                // Deposit its energy.
                calo->fillFast(particle.p4().eta(), caloPhi, particle.p4().e());
        }
    }

    // Add the calorimeter noise
    const unsigned nEta = calo->nEta();
    const unsigned nPhi = calo->nPhi();
    for (unsigned ieta=0; ieta<nEta; ++ieta)
        for (unsigned iphi=0; iphi<nPhi; ++iphi)
        {
            const double noiseE = PythiaJetGun::gaussRandom(0.0, caloNoise);
            calo->uncheckedFillBin(ieta, iphi, noiseE);
        }

    // Apply the calorimeter threshold. This is what is normally
    // done by every experiment to suppress noise contributions
    // and to reduce the readout size.
    calo->applyThreshold(caloThreshold);

    // Convert the energies in the grid into the transverse energies.
    // We will use the cell center to find the pseudorapidity, and
    // will assume that the mass associated with a single cell energy
    // deposit is 0. This procedure effectively sets the directions
    // of all particles which have reached a particular cell to the
    // direction towards the center of that cell. Something similar
    // happens in real calorimeters.
    for (unsigned ieta=0; ieta<nEta; ++ieta)
    {
        // Note that p = pt*cosh(eta)
        const double coshCaloEta = cosh(calo->etaBinCenter(ieta));
        for (unsigned iphi=0; iphi<nPhi; ++iphi)
        {
            const double e = calo->uncheckedAt(ieta, iphi);
            if (e > 0.0)
                calo->uncheckedSetBin(ieta, iphi, e/coshCaloEta);
        }
    }
}
