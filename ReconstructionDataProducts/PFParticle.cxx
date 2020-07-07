#include "ReconstructionDataProducts/PFParticle.h"

namespace gar {
    namespace rec {

        //--------------------------------------------------------------------------
        //Default constructor
        PFParticle::PFParticle(){
            // The default constructor is used e.g. by art::DataViewImpl::getByLabel
            // Make sure all Cluster objects are numbered, lest art deep-copy uninitialized
            // instances and then operator==() evaluates meaninglessly true.
            IDNumberGen::create(FirstNumber);
            fIDnumero = IDNumberGen::create()->getNewOne();
            return;
        }

        PFParticle::PFParticle(int type, float energy, float pos[3], float mom[3], float charge, std::vector<float> cov, int pdg, float goodness, size_t parent, std::vector<size_t> daughters)
        : fType(type),
        fEnergy(energy),
        fCharge(charge),
        fCov(cov),
        fParent(parent),
        fDaughters(daughters)
        {
            fPos[0] = pos[0];
            fPos[1] = pos[1];
            fPos[2] = pos[2];

            fMom[0] = mom[0];
            fMom[1] = mom[1];
            fMom[2] = mom[2];

            fPID.fPdg = pdg;
            fPID.fGoodness = goodness;
        }

        //--------------------------------------------------------------------------
        void PFParticle::setType(int type) {
            fType = type;
        }

        //--------------------------------------------------------------------------
        void PFParticle::setEnergy(float energy) {
            fEnergy = energy;
        }

        //--------------------------------------------------------------------------
        void PFParticle::setPosition(const float pos[3]) {
            fPos[0] = pos[0];
            fPos[1] = pos[1];
            fPos[2] = pos[2];
        }

        //--------------------------------------------------------------------------
        void PFParticle::setMomentum(const float mom[3]) {
            fMom[0] = mom[0];
            fMom[1] = mom[1];
            fMom[2] = mom[2];
        }

        //--------------------------------------------------------------------------
        void PFParticle::setMass(float mass) {
            fMass = mass;
        }

        //--------------------------------------------------------------------------
        void PFParticle::setCharge(float charge) {
            fCharge = charge;
        }

        //--------------------------------------------------------------------------
        void PFParticle::setCovMatrix(std::vector<float> cov) {
            fCov = cov;
        }

        //--------------------------------------------------------------------------
        void PFParticle::setParticleID(int pdg, float goodness) {
            fPID.fPdg = pdg;
            fPID.fGoodness = goodness;
        }

        //--------------------------------------------------------------------------
        void PFParticle::setParent(size_t parent) {
            fParent = parent;
        }

        //--------------------------------------------------------------------------
        void PFParticle::setDaughters(std::vector<size_t> daughters) {
            fDaughters = daughters;
        }

        //--------------------------------------------------------------------------
        std::ostream& operator<< (std::ostream& o, gar::rec::PFParticle const& h)
        {
            o << "Reconstructed Particle "
            << "\n\tType = "
            << h.Type()
            << "\n\tEnergy = "
            << h.Energy()
            << "\n\tMomentum = "
            << h.Momentum()[0] << ", " <<  h.Momentum()[1] << ", " << h.Momentum()[2] << ")"
            << "\n\tPID = "
            << h.GetParticleID().fPdg << ", goodness " << h.GetParticleID().fGoodness
            << "\n\tMass = "
            << h.Mass()
            << "\n\tID number = "
            << h.getIDNumber()
            << "\n\tPosition = "
            << h.Position()[0] << ", " <<  h.Position()[1] << ", " << h.Position()[2] << ")";

            return o;
        }



        //--------------------------------------------------------------------------
        // ID number methods
        bool PFParticle::operator==(const PFParticle& rhs) const {
            return (this->fIDnumero == rhs.fIDnumero);
        }

        bool PFParticle::operator!=(const PFParticle& rhs) const {
            return (this->fIDnumero != rhs.fIDnumero);
        }

        gar::rec::IDNumber PFParticle::getIDNumber() const { return fIDnumero; }
    }
}
