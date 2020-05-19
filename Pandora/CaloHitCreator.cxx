#include "CaloHitCreator.h"

#include <algorithm>
#include <cmath>
#include <limits>


namespace gar {
    namespace gar_pandora {

        CaloHitCreator::CaloHitCreator(const Settings &settings, const pandora::Pandora *const pPandora) :
        m_settings(settings),
        m_pandora(*pPandora),
        m_eCalBarrelLayerThickness(0.f),
        m_eCalEndCapLayerThickness(0.f),
        m_calorimeterHitVector(0)
        {
            m_eCalEndCapLayerThickness = 20.;

            if ((m_eCalEndCapLayerThickness < std::numeric_limits<float>::epsilon()))
                throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        CaloHitCreator::~CaloHitCreator()
        {
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