#include "CaloHitCreator.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "CoreUtils/ServiceUtil.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace gar {
    namespace gar_pandora {

        CaloHitCreator::CaloHitCreator(const Settings &settings, const pandora::Pandora *const pPandora) :
        m_settings(settings),
        m_pandora(*pPandora),
        m_eCalBarrelLayerThickness(0.8),
        m_eCalEndCapLayerThickness(0.8),
        m_calorimeterHitVector(0)
        {
            fGeo = gar::providerFrom<geo::Geometry>();
            std::string fEncoding = fGeo->GetECALCellIDEncoding();
            m_fieldDecoder = new gar::geo::BitFieldCoder( fEncoding );

            if ((m_eCalEndCapLayerThickness < std::numeric_limits<float>::epsilon()))
                throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        CaloHitCreator::~CaloHitCreator()
        {
            if ( m_fieldDecoder ) delete m_fieldDecoder;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode CaloHitCreator::CreateCaloHits(const art::Event *const pEvent)
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateECalCaloHits(pEvent));

            return pandora::STATUS_CODE_SUCCESS;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode CaloHitCreator::CreateECalCaloHits(const art::Event *const pEvent)
        {
            std::vector<art::Ptr<gar::rec::CaloHit>> hitVector;
            art::Handle< std::vector<gar::rec::CaloHit> > theHits;
            pEvent->getByLabel(m_settings.m_CaloHitCollection, theHits);

            if (!theHits.isValid())
            {
                LOG_DEBUG("CaloHitCreator") << "  Failed to find hits... " << std::endl;
                return pandora::STATUS_CODE_FAILURE;
            }
            else
            {
                LOG_DEBUG("CaloHitCreator") << "  Found: " << theHits->size() << " Hits " << std::endl;
            }

            for (unsigned int i = 0; i < theHits->size(); ++i)
            {
                const art::Ptr<gar::rec::CaloHit> hit(theHits, i);
                hitVector.push_back(hit);
            }

            for (std::vector<art::Ptr<gar::rec::CaloHit>>::const_iterator iter = hitVector.begin(), iterEnd = hitVector.end(); iter != iterEnd; ++iter)
            {
                art::Ptr<gar::rec::CaloHit> artPtrCaloHit = *iter;
                const gar::rec::CaloHit *pCaloHit = artPtrCaloHit.get();

                PandoraApi::CaloHit::Parameters caloHitParameters;
                caloHitParameters.m_hitType = pandora::ECAL;
                caloHitParameters.m_isDigital = false;
                caloHitParameters.m_layer = m_fieldDecoder->get(pCaloHit->CellID(), "layer");
                caloHitParameters.m_isInOuterSamplingLayer = false;
                this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                float absorberCorrection(1.);

                if (std::fabs(pCaloHit->Position()[2]) < m_settings.m_eCalEndCapInnerX)
                {
                    LOG_DEBUG("CaloHitCreator") << "IDS " << iter
                    << std::setw(15) << pCaloHit->CellID()
                    << std::setw(15) << pCaloHit->Position()[0]
                    << std::setw(15) << pCaloHit->Position()[1]
                    << std::setw(15) << pCaloHit->Position()[2]
                    << std::setw(5) << m_fieldDecoder->get(pCaloHit->CellID(), "system")
                    << std::setw(5) << m_fieldDecoder->get(pCaloHit->CellID(), "module")
                    << std::setw(5) << m_fieldDecoder->get(pCaloHit->CellID(), "stave")
                    << std::setw(5) << m_fieldDecoder->get(pCaloHit->CellID(), "layer")
                    << std::setw(5) << m_fieldDecoder->get(pCaloHit->CellID(), "cellX")
                    << std::setw(5) << m_fieldDecoder->get(pCaloHit->CellID(), "cellY")
                    << std::endl;

                    // this->GetBarrelCaloHitProperties(pCaloHit, barrelLayers, m_settings.m_eCalBarrelInnerSymmetry, caloHitParameters, m_settings.m_eCalBarrelNormalVector, absorberCorrection);
                    caloHitParameters.m_hadronicEnergy = m_settings.m_eCalToHadGeVBarrel * pCaloHit->Energy();
                }
                else
                {
                    // this->GetEndCapCaloHitProperties(pCaloHit, endcapLayers, caloHitParameters, absorberCorrection);
                    caloHitParameters.m_hadronicEnergy = m_settings.m_eCalToHadGeVEndCap * pCaloHit->Energy();
                }

                caloHitParameters.m_mipEquivalentEnergy = pCaloHit->Energy() * m_settings.m_eCalToMip * absorberCorrection;

                if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_eCalMipThreshold)
                    continue;

                caloHitParameters.m_electromagneticEnergy = m_settings.m_eCalToEMGeV * pCaloHit->Energy();

                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(m_pandora, caloHitParameters));
                m_calorimeterHitVector.push_back(pCaloHit);
            }

            return pandora::STATUS_CODE_SUCCESS;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        void CaloHitCreator::GetCommonCaloHitProperties(const gar::rec::CaloHit *const pCaloHit, PandoraApi::CaloHit::Parameters &caloHitParameters) const
        {
            const float *pCaloHitPosition(pCaloHit->Position());
            const pandora::CartesianVector positionVector(pCaloHitPosition[0], pCaloHitPosition[1], pCaloHitPosition[2]);

            caloHitParameters.m_cellGeometry = pandora::RECTANGULAR;
            caloHitParameters.m_positionVector = positionVector;
            caloHitParameters.m_expectedDirection = positionVector.GetUnitVector();
            caloHitParameters.m_pParentAddress = (void*)pCaloHit;
            caloHitParameters.m_inputEnergy = pCaloHit->Energy();
            caloHitParameters.m_time = pCaloHit->Time().first;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------


        //------------------------------------------------------------------------------------------------------------------------------------------

        int CaloHitCreator::GetNLayersFromEdge(const gar::rec::CaloHit *const pCaloHit) const
        {
            return 0;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        float CaloHitCreator::GetMaximumRadius(const gar::rec::CaloHit *const pCaloHit, const unsigned int symmetryOrder, const float phi0) const
        {
            return 0.;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        CaloHitCreator::Settings::Settings()
        : m_CaloHitCollection( "" ),
        m_eCalToMip(1.f),
        m_eCalMipThreshold(0.f),
        m_eCalToEMGeV(1.f),
        m_eCalToHadGeVBarrel(1.f),
        m_eCalToHadGeVEndCap(1.f),
        m_maxECalHitHadronicEnergy(10000.f),
        m_nOuterSamplingLayers(3),
        m_layersFromEdgeMaxRearDistance(250.f),
        m_eCalBarrelInnerPhi0(0.f),
        m_eCalBarrelInnerSymmetry(0.f),
        m_eCalEndCapOuterR(0.f),
        m_eCalEndCapOuterX(0.f),
        m_eCalBarrelOuterR(0.f),
        m_eCalBarrelNormalVector({0.0, 0.0, 1.0})
        {
        }
    }
}
