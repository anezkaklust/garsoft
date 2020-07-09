#include "art/Framework/Principal/Handle.h"
#include "nutools/MagneticField/MagneticField.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "CaloHitCreator.h"
#include "MCParticleCreator.h"
#include "TrackCreator.h"

#include "MCCheater/BackTracker.h"

#include <cmath>
#include <limits>

namespace gar {
    namespace gar_pandora {

        MCParticleCreator::MCParticleCreator(const Settings &settings, const pandora::Pandora *const pPandora)
        : m_settings(settings),
        m_pandora(*pPandora)
        {
            art::ServiceHandle<mag::MagneticField> magFieldService;
            G4ThreeVector zerovec(0, 0, 0);
            G4ThreeVector magfield = magFieldService->FieldAtPoint(zerovec);
            m_bField = magfield[0]; //x component at (0, 0, 0)
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        MCParticleCreator::~MCParticleCreator()
        {
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode MCParticleCreator::CollectMCParticles(const art::Event &pEvent)
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CollectMCParticles(pEvent, m_settings.m_geantModuleLabel, artMCParticleVector));
            if (!m_settings.m_generatorModuleLabel.empty())
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CollectGeneratorMCParticles(pEvent, m_settings.m_generatorModuleLabel, generatorArtMCParticleVector));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CollectMCParticles(pEvent, m_settings.m_geantModuleLabel, artMCTruthToMCParticles, artMCParticlesToMCTruth));

            return pandora::STATUS_CODE_SUCCESS;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode MCParticleCreator::CreateMCParticles() const
        {
            MCParticleMap particleMap;
            for (MCParticlesToMCTruth::const_iterator iter = artMCParticlesToMCTruth.begin(), iterEnd = artMCParticlesToMCTruth.end(); iter != iterEnd; ++iter)
            {
                const art::Ptr<simb::MCParticle> particle = iter->first;
                particleMap[particle->TrackId()] = particle;
            }

            int neutrinoCounter(0);
            for (MCTruthToMCParticles::const_iterator iter1 = artMCTruthToMCParticles.begin(), iterEnd1 = artMCTruthToMCParticles.end(); iter1 != iterEnd1; ++iter1)
            {
                const art::Ptr<simb::MCTruth> truth = iter1->first;
                if (truth->NeutrinoSet())
                {
                    const simb::MCNeutrino neutrino(truth->GetNeutrino());
                    ++neutrinoCounter;
                    const int neutrinoID(neutrinoCounter);

                    PandoraApi::MCParticle::Parameters mcParticleParameters;
                    mcParticleParameters.m_energy = neutrino.Nu().E();
                    mcParticleParameters.m_momentum = pandora::CartesianVector(neutrino.Nu().Pz(), neutrino.Nu().Py(), -neutrino.Nu().Px());
                    mcParticleParameters.m_vertex = pandora::CartesianVector(neutrino.Nu().Vz(), neutrino.Nu().Vy(), -neutrino.Nu().Vx()) * CLHEP::cm;
                    mcParticleParameters.m_endpoint = pandora::CartesianVector(neutrino.Nu().Vz(), neutrino.Nu().Vy(), -neutrino.Nu().Vx()) * CLHEP::cm;
                    mcParticleParameters.m_particleId = neutrino.Nu().PdgCode();
                    mcParticleParameters.m_mcParticleType = pandora::MC_3D;
                    mcParticleParameters.m_pParentAddress = (void*)((intptr_t)neutrinoID);

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(m_pandora, mcParticleParameters));

                    // Loop over associated particles
                    const MCParticleVector &particleVector = iter1->second;
                    for (MCParticleVector::const_iterator iter2 = particleVector.begin(), iterEnd2 = particleVector.end(); iter2 != iterEnd2; ++iter2)
                    {
                        const art::Ptr<simb::MCParticle> particle = *iter2;
                        const int trackID(particle->TrackId());
                        // Mother/Daughter Links
                        if (particle->Mother() == 0)
                        {
                            try
                            {
                                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(m_pandora, (void*)((intptr_t)neutrinoID), (void*)((intptr_t)trackID)));
                            }
                            catch (const pandora::StatusCodeException &)
                            {
                                LOG_WARNING("MCParticleCreator") << "CreatePandoraMCParticles - unable to create mc particle relationship, invalid information supplied " << std::endl;
                                continue;
                            }
                        }
                    }
                }
            }

            LOG_DEBUG("MCParticleCreator") << " Number of Pandora neutrinos: " << neutrinoCounter << std::endl;

            for (MCParticleMap::const_iterator iterI = particleMap.begin(), iterEndI = particleMap.end(); iterI != iterEndI; ++iterI)
            {
                const art::Ptr<simb::MCParticle> pMcParticle = iterI->second;

                // Lookup position and kinematics at start and end points
                const float vtxX(pMcParticle->Vx());
                const float vtxY(pMcParticle->Vy());
                const float vtxZ(pMcParticle->Vz());
                const float endX(pMcParticle->EndX());
                const float endY(pMcParticle->EndY());
                const float endZ(pMcParticle->EndZ());
                const float pX(pMcParticle->Px());
                const float pY(pMcParticle->Py());
                const float pZ(pMcParticle->Pz());
                const float E(pMcParticle->E());

                PandoraApi::MCParticle::Parameters mcParticleParameters;
                mcParticleParameters.m_energy = E;
                mcParticleParameters.m_particleId = pMcParticle->PdgCode();
                mcParticleParameters.m_mcParticleType = pandora::MC_3D;
                mcParticleParameters.m_pParentAddress = (void*)((intptr_t)pMcParticle->TrackId());
                mcParticleParameters.m_momentum = pandora::CartesianVector(pZ, pY, -pX);
                mcParticleParameters.m_vertex = pandora::CartesianVector(vtxZ, vtxY, -vtxX) * CLHEP::cm;
                mcParticleParameters.m_endpoint = pandora::CartesianVector(endZ, endY, -endX) * CLHEP::cm;

                LOG_DEBUG("MCParticleCreator") << " Adding MC Particle with parameters "
                << " mcParticleParameters.m_energy = " << mcParticleParameters.m_energy.Get()
                << " mcParticleParameters.m_particleId = " << mcParticleParameters.m_particleId.Get()
                << " mcParticleParameters.m_mcParticleType = " << mcParticleParameters.m_mcParticleType.Get()
                << " mcParticleParameters.m_pParentAddress = " << mcParticleParameters.m_pParentAddress.Get()
                << " mcParticleParameters.m_momentum = " << mcParticleParameters.m_momentum.Get()
                << " mcParticleParameters.m_vertex = " << mcParticleParameters.m_vertex.Get()
                << " mcParticleParameters.m_endpoint = " << mcParticleParameters.m_endpoint.Get();

                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(m_pandora, mcParticleParameters));

                // Create parent-daughter relationships
                const int id_mother(pMcParticle->Mother());
                MCParticleMap::const_iterator iterJ = particleMap.find(id_mother);
                if (iterJ != particleMap.end())
                {
                    try
                    {
                        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(m_pandora, (void*)((intptr_t)id_mother), (void*)((intptr_t)pMcParticle->TrackId())));
                    }
                    catch (const pandora::StatusCodeException &)
                    {
                        LOG_WARNING("MCParticleCreator") << "CreatePandoraMCParticles - Unable to create mc particle relationship, invalid information supplied " << std::endl;
                        continue;
                    }
                }
            }

            return pandora::STATUS_CODE_SUCCESS;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode MCParticleCreator::CreateTrackToMCParticleRelationships(const TrackVector &trackVector) const
        {
            return pandora::STATUS_CODE_SUCCESS;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode MCParticleCreator::CreateCaloHitToMCParticleRelationships(const CalorimeterHitVector &calorimeterHitVector) const
        {
            cheat::BackTrackerCore const* bt = gar::providerFrom<cheat::BackTracker>();
            // loop over all hits and fill in the map
            for( auto const& itr : calorimeterHitVector )
            {
                std::vector<gar::cheat::CalIDE> eveides = bt->CaloHitToCalIDEs(itr);
                // loop over all eveides for this hit
                for(size_t ieve = 0; ieve < eveides.size(); ieve++)
                {
                    // don't worry about eve particles that contribute less than 10% of the
                    // energy in the current hit
                    if( eveides[ieve].energyFrac < 0.1) continue;
                    const simb::MCParticle *part = bt->TrackIDToParticle(eveides[ieve].trackID);

                    if(!part) {
                        LOG_WARNING("MCParticleCreator")
                        << "Cannot find MCParticle for eveid: " << eveides[ieve].trackID;
                        continue;
                    }

                    const int hitID = eveides[ieve].trackID;
                    const int trackID = part->TrackId();
                    const float energyFrac(eveides[ieve].energyFrac);

                    try
                    {
                        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetCaloHitToMCParticleRelationship(m_pandora, (void*)((intptr_t)hitID), (void*)((intptr_t)trackID), energyFrac));
                    }
                    catch (const pandora::StatusCodeException &)
                    {
                        LOG_WARNING("MCParticleCreator") << "CreateCaloHitToMCParticleRelationships - unable to create calo hit to mc particle relationship, invalid information supplied " << std::endl;
                        continue;
                    }

                } // end loop over eve IDs for this hit
            }// end loop over hits

            return pandora::STATUS_CODE_SUCCESS;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode MCParticleCreator::CollectMCParticles(const art::Event &pEvent, const std::string &label, MCParticleVector &particleVector)
        {
            art::Handle< RawMCParticleVector > theParticles;
            pEvent.getByLabel(label, theParticles);

            if (!theParticles.isValid())
            {
                LOG_WARNING("MCParticleCreator") << "  Failed to find MC particles for label " << label << std::endl;
                return pandora::STATUS_CODE_NOT_FOUND;
            }

            LOG_DEBUG("MCParticleCreator") << "  Found: " << theParticles->size() << " MC particles " << std::endl;

            for (unsigned int i = 0; i < theParticles->size(); ++i)
            {
                const art::Ptr<simb::MCParticle> particle(theParticles, i);
                particleVector.push_back(particle);
            }

            return pandora::STATUS_CODE_SUCCESS;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode MCParticleCreator::CollectGeneratorMCParticles(const art::Event &pEvent, const std::string &label, RawMCParticleVector &particleVector)
        {
            art::Handle< std::vector<simb::MCTruth> > mcTruthBlocks;
            pEvent.getByLabel(label, mcTruthBlocks);

            if (!mcTruthBlocks.isValid())
            {
                LOG_WARNING("MCParticleCreator") << "  Failed to find MC Truth for generator " << label << std::endl;
                return pandora::STATUS_CODE_NOT_FOUND;
            }

            LOG_DEBUG("MCParticleCreator") << "  Found: " << mcTruthBlocks->size() << " MC truth blocks " << std::endl;

            if (mcTruthBlocks->size() != 1)
            throw cet::exception("MCParticleCreator") << " PandoraCollector::CollectGeneratorMCParticles --- Unexpected number of MC truth blocks ";

            const art::Ptr<simb::MCTruth> mcTruth(mcTruthBlocks, 0);
            for (int i = 0; i < mcTruth->NParticles(); ++i)
            {
                particleVector.push_back(mcTruth->GetParticle(i));
            }

            return pandora::STATUS_CODE_SUCCESS;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode MCParticleCreator::CollectMCParticles(const art::Event &pEvent, const std::string &label, MCTruthToMCParticles &truthToParticles, MCParticlesToMCTruth &particlesToTruth)
        {
            art::Handle< RawMCParticleVector > theParticles;
            pEvent.getByLabel(label, theParticles);

            if (!theParticles.isValid())
            {
                LOG_WARNING("MCParticleCreator") << "  Failed to find MC particles for label " << label << std::endl;
                return pandora::STATUS_CODE_NOT_FOUND;
            }

            LOG_DEBUG("MCParticleCreator") << "  Found: " << theParticles->size() << " MC particles " << std::endl;

            art::FindOneP<simb::MCTruth> theTruthAssns(theParticles, pEvent, label);

            for (unsigned int i = 0, iEnd = theParticles->size(); i < iEnd; ++i)
            {
                const art::Ptr<simb::MCParticle> particle(theParticles, i);
                const art::Ptr<simb::MCTruth> truth(theTruthAssns.at(i));
                truthToParticles[truth].push_back(particle);
                particlesToTruth[particle] = truth;
            }

            return pandora::STATUS_CODE_SUCCESS;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        MCParticleCreator::Settings::Settings()
        : m_geantModuleLabel( "" ),
        m_generatorModuleLabel( "" )
        {
        }
    }
}
