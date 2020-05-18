#include "CaloHitCreator.h"

#include <algorithm>
#include <cmath>
#include <limits>


namespace gar {
    namespace gar_pandora {

        CaloHitCreator::CaloHitCreator(const Settings &settings, const pandora::Pandora *const pPandora) :
        m_settings(settings),
        m_pandora(*pPandora),
        m_hCalBarrelLayerThickness(0.f),
        m_hCalEndCapLayerThickness(0.f),
        m_calorimeterHitVector(0)
        {
            m_hCalEndCapLayerThickness = 20.;

            if ((m_hCalEndCapLayerThickness < std::numeric_limits<float>::epsilon()))
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
            // Calo hit coordinate calculations
            const float barrelMaximumRadius(this->GetMaximumRadius(pCaloHit, m_settings.m_hCalBarrelOuterSymmetry, m_settings.m_hCalBarrelOuterPhi0));
            const float endCapMaximumRadius(this->GetMaximumRadius(pCaloHit, m_settings.m_hCalEndCapInnerSymmetryOrder, m_settings.m_hCalEndCapInnerPhiCoordinate));
            const float caloHitAbsZ(std::fabs(pCaloHit->Position()[2]));

            // Distance from radial outer
            float radialDistanceToEdge(std::numeric_limits<float>::max());

            if (caloHitAbsZ < m_settings.m_eCalBarrelOuterZ)
            {
                radialDistanceToEdge = (m_settings.m_hCalBarrelOuterR - barrelMaximumRadius) / m_hCalBarrelLayerThickness;
            }
            else
            {
                radialDistanceToEdge = (m_settings.m_hCalEndCapOuterR - endCapMaximumRadius) / m_hCalEndCapLayerThickness;
            }

            // Distance from rear of endcap outer
            float rearDistanceToEdge(std::numeric_limits<float>::max());

            if (caloHitAbsZ >= m_settings.m_eCalBarrelOuterZ)
            {
                rearDistanceToEdge = (m_settings.m_hCalEndCapOuterZ - caloHitAbsZ) / m_hCalEndCapLayerThickness;
            }
            else
            {
                const float rearDistance((m_settings.m_eCalBarrelOuterZ - caloHitAbsZ) / m_hCalBarrelLayerThickness);

                if (rearDistance < m_settings.m_layersFromEdgeMaxRearDistance)
                {
                    const float overlapDistance((m_settings.m_hCalEndCapOuterR - endCapMaximumRadius) / m_hCalEndCapLayerThickness);
                    rearDistanceToEdge = std::max(rearDistance, overlapDistance);
                }
            }

            return static_cast<int>(std::min(radialDistanceToEdge, rearDistanceToEdge));
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        float CaloHitCreator::GetMaximumRadius(const gar::rec::CaloHit *const pCaloHit, const unsigned int symmetryOrder, const float phi0) const
        {
            const float *pCaloHitPosition(pCaloHit->Position());

            if (symmetryOrder <= 2)
                return std::sqrt((pCaloHitPosition[0] * pCaloHitPosition[0]) + (pCaloHitPosition[1] * pCaloHitPosition[1]));

            float maximumRadius(0.f);
            const float twoPi(2.f * M_PI);

            for (unsigned int i = 0; i < symmetryOrder; ++i)
            {
                const float phi = phi0 + i * twoPi / static_cast<float>(symmetryOrder);
                float radius = pCaloHitPosition[0] * std::cos(phi) + pCaloHitPosition[1] * std::sin(phi);

                if (radius > maximumRadius)
                    maximumRadius = radius;
            }

            return maximumRadius;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        CaloHitCreator::Settings::Settings()
        : m_CaloHitCollections( StringVector() ),
        m_MuIDCaloHitCollections( StringVector() ),
        m_eCalToMip(1.f),
        m_hCalToMip(1.f),
        m_muonToMip(1.f),
        m_eCalMipThreshold(0.f),
        m_hCalMipThreshold(0.f),
        m_muonMipThreshold(0.f),
        m_eCalToEMGeV(1.f),
        m_eCalToHadGeVBarrel(1.f),
        m_eCalToHadGeVEndCap(1.f),
        m_hCalToEMGeV(1.f),
        m_hCalToHadGeV(1.f),
        m_maxHCalHitHadronicEnergy(10000.f),
        m_nOuterSamplingLayers(3),
        m_layersFromEdgeMaxRearDistance(250.f),
        m_hCalEndCapInnerSymmetryOrder(0),
        m_hCalEndCapInnerPhiCoordinate(0.f),
        m_stripSplittingOn(0),
        m_eCalBarrelOuterZ(0.f),
        m_hCalBarrelOuterZ(0.f),
        m_muonBarrelOuterZ(0.f),
        m_coilOuterR(0.f),
        m_eCalBarrelInnerPhi0(0.f),
        m_eCalBarrelInnerSymmetry(0.f),
        m_hCalBarrelInnerPhi0(0.f),
        m_hCalBarrelInnerSymmetry(0.f),
        m_muonBarrelInnerPhi0(0.f),
        m_muonBarrelInnerSymmetry(0.f),
        m_hCalEndCapOuterR(0.f),
        m_hCalEndCapOuterZ(0.f),
        m_hCalBarrelOuterR(0.f),
        m_hCalBarrelOuterPhi0(0.f),
        m_hCalBarrelOuterSymmetry(0.f),
        m_eCalBarrelNormalVector({0.0, 0.0, 1.0}),
        m_hCalBarrelNormalVector({0.0, 0.0, 1.0}),
        m_muonBarrelNormalVector({0.0, 0.0, 1.0})
        {
        }
    }
}