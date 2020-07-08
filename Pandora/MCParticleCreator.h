#ifndef MCPARTICLECREATOR_H
#define MCPARTICLECREATOR_H 1

#include "art/Framework/Principal/Event.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "Api/PandoraApi.h"

#include "CaloHitCreator.h"
#include "TrackCreator.h"

namespace gar {
    namespace gar_pandora {

        typedef std::vector< art::Ptr<simb::MCTruth> >      MCTruthVector;
        typedef std::vector< art::Ptr<simb::MCParticle> >   MCParticleVector;
        typedef std::vector< simb::MCParticle>              RawMCParticleVector;
        typedef std::map< art::Ptr<simb::MCTruth>,     MCParticleVector >             MCTruthToMCParticles;
        typedef std::map< art::Ptr<simb::MCParticle>,  art::Ptr<simb::MCTruth> >      MCParticlesToMCTruth;
        typedef std::map< int, art::Ptr<simb::MCParticle> >   MCParticleMap;

        class MCParticleCreator
        {
        public:

            class Settings
            {
            public:
                Settings();

                std::string    m_geantModuleLabel;             ///< The geant4 label
                std::string    m_generatorModuleLabel;         ///< The generator label
            };

            MCParticleCreator(const Settings &settings, const pandora::Pandora *const pPandora);

            ~MCParticleCreator();

            pandora::StatusCode CollectMCParticles(const art::Event &pEvent);
            pandora::StatusCode CreateMCParticles() const;
            pandora::StatusCode CreateTrackToMCParticleRelationships(const TrackVector &trackVector) const;
            pandora::StatusCode CreateCaloHitToMCParticleRelationships(const CalorimeterHitVector &calorimeterHitVector) const;

            void Reset();

        protected:

            pandora::StatusCode CollectMCParticles(const art::Event &pEvent, const std::string &label, MCParticleVector &particleVector);

            pandora::StatusCode CollectGeneratorMCParticles(const art::Event &pEvent, const std::string &label, RawMCParticleVector &particleVector);

            pandora::StatusCode CollectMCParticles(const art::Event &pEvent, const std::string &label, MCTruthToMCParticles &truthToParticles, MCParticlesToMCTruth &particlesToTruth);

        private:
            const Settings          m_settings;        ///< The mc particle creator settings
            const pandora::Pandora &m_pandora;         ///< Reference to the pandora object to create the mc particles
            float                   m_bField;          ///< The bfield

            MCParticleVector artMCParticleVector;
            RawMCParticleVector generatorArtMCParticleVector;
            MCTruthToMCParticles artMCTruthToMCParticles;
            MCParticlesToMCTruth artMCParticlesToMCTruth;
        };

        inline void MCParticleCreator::Reset()
        {
            artMCParticleVector.clear();
            generatorArtMCParticleVector.clear();
            artMCTruthToMCParticles.clear();
            artMCParticlesToMCTruth.clear();
        }
    }
}

#endif // #ifndef MCPARTICLECREATOR_H
