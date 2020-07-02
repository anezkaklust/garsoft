#include "CaloHitCreator.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "CoreUtils/ServiceUtil.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace gar {
    namespace gar_pandora {

        CaloHitCreator::CaloHitCreator(const Settings &settings, const pandora::Pandora *const pPandora)
        : m_settings(settings),
        m_pandora(*pPandora),
        m_eCalBarrelLayerThickness(0.f),
        m_eCalEndCapLayerThickness(0.f),
        m_calorimeterHitVector(0)
        {
            fGeo = gar::providerFrom<geo::Geometry>();
            std::string fEncoding = fGeo->GetECALCellIDEncoding();
            m_fieldDecoder = new gar::geo::BitFieldCoder( fEncoding );

            const std::vector<gar::geo::LayeredCalorimeterStruct::Layer> &barrelLayers = (fGeo->GetECALLayeredCalorimeterData()[gar::geo::LayeredCalorimeterData::BarrelLayout].get())->layers;

            const std::vector<gar::geo::LayeredCalorimeterStruct::Layer> &endcapLayers = (fGeo->GetECALLayeredCalorimeterData()[gar::geo::LayeredCalorimeterData::EndcapLayout].get())->layers;

            ///Take thicknesses from last layer (was like that before with gear)
            m_eCalEndCapLayerThickness = (endcapLayers.back().inner_thickness + endcapLayers.back().outer_thickness);
            m_eCalBarrelLayerThickness = (barrelLayers.back().inner_thickness + barrelLayers.back().outer_thickness);

            if ( (m_eCalEndCapLayerThickness < std::numeric_limits<float>::epsilon()) || (m_eCalBarrelLayerThickness < std::numeric_limits<float>::epsilon()) )
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

            const std::vector<gar::geo::LayeredCalorimeterStruct::Layer> &barrelLayers = (fGeo->GetECALLayeredCalorimeterData()[gar::geo::LayeredCalorimeterData::BarrelLayout].get())->layers;

            const std::vector<gar::geo::LayeredCalorimeterStruct::Layer> &endcapLayers = (fGeo->GetECALLayeredCalorimeterData()[gar::geo::LayeredCalorimeterData::EndcapLayout].get())->layers;

            for (std::vector<art::Ptr<gar::rec::CaloHit>>::const_iterator iter = hitVector.begin(), iterEnd = hitVector.end(); iter != iterEnd; ++iter)
            {
                art::Ptr<gar::rec::CaloHit> artPtrCaloHit = *iter;
                const gar::rec::CaloHit *pCaloHit = artPtrCaloHit.get();

                PandoraApi::CaloHit::Parameters caloHitParameters;
                caloHitParameters.m_hitType = pandora::ECAL;
                caloHitParameters.m_isDigital = false;
                caloHitParameters.m_layer = m_fieldDecoder->get(pCaloHit->CellID(), "layer");
                caloHitParameters.m_isInOuterSamplingLayer = (this->GetNLayersFromEdge(pCaloHit) <= m_settings.m_nOuterSamplingLayers);
                this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                float absorberCorrection(1.);

                if (std::fabs(pCaloHit->Position()[2]) < m_settings.m_eCalBarrelOuterZ)
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

                    this->GetBarrelCaloHitProperties(pCaloHit, barrelLayers, m_settings.m_eCalBarrelInnerSymmetry, caloHitParameters, m_settings.m_eCalBarrelNormalVector, absorberCorrection);
                    caloHitParameters.m_hadronicEnergy = m_settings.m_eCalToHadGeVBarrel * pCaloHit->Energy();
                }
                else
                {
                    this->GetEndCapCaloHitProperties(pCaloHit, endcapLayers, caloHitParameters, absorberCorrection);
                    caloHitParameters.m_hadronicEnergy = m_settings.m_eCalToHadGeVEndCap * pCaloHit->Energy();
                }

                caloHitParameters.m_mipEquivalentEnergy = pCaloHit->Energy() * m_settings.m_eCalToMip * absorberCorrection;

                if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_eCalMipThreshold)
                continue;

                caloHitParameters.m_electromagneticEnergy = m_settings.m_eCalToEMGeV * pCaloHit->Energy();

                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(m_pandora, caloHitParameters));
                m_calorimeterHitVector.push_back(artPtrCaloHit);
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

        void CaloHitCreator::GetEndCapCaloHitProperties(const gar::rec::CaloHit *const pCaloHit, const std::vector<gar::geo::LayeredCalorimeterStruct::Layer> &layers, PandoraApi::CaloHit::Parameters &caloHitParameters, float &absorberCorrection) const
        {
            caloHitParameters.m_hitRegion = pandora::ENDCAP;

            const int physicalLayer(std::min(static_cast<int>(caloHitParameters.m_layer.Get()), static_cast<int>(layers.size()-1)));
            caloHitParameters.m_cellSize0 = layers[physicalLayer].cellSize0;
            caloHitParameters.m_cellSize1 = layers[physicalLayer].cellSize1;

            double thickness = (layers[physicalLayer].inner_thickness+layers[physicalLayer].sensitive_thickness/2.0);
            double nRadLengths = layers[physicalLayer].inner_nRadiationLengths;
            double nIntLengths = layers[physicalLayer].inner_nInteractionLengths;
            double layerAbsorberThickness = (layers[physicalLayer].inner_thickness-layers[physicalLayer].sensitive_thickness/2.0);

            if(physicalLayer>0){
                thickness += (layers[physicalLayer-1].outer_thickness -layers[physicalLayer].sensitive_thickness/2.0);
                nRadLengths += layers[physicalLayer-1].outer_nRadiationLengths;
                nIntLengths += layers[physicalLayer-1].outer_nInteractionLengths;
                layerAbsorberThickness += (layers[physicalLayer-1].outer_thickness -layers[physicalLayer].sensitive_thickness/2.0);

            }

            caloHitParameters.m_cellThickness = thickness;
            caloHitParameters.m_nCellRadiationLengths = nRadLengths;
            caloHitParameters.m_nCellInteractionLengths = nIntLengths;

            if (caloHitParameters.m_nCellRadiationLengths.Get() < std::numeric_limits<float>::epsilon() || caloHitParameters.m_nCellInteractionLengths.Get() < std::numeric_limits<float>::epsilon())
            {
                LOG_ERROR("CaloHitCreator::GetEndcapCaloHitProperties")
                << "Calo hit has 0 radiation length or interaction length: not creating a Pandora calo hit.";
                throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
            }

            absorberCorrection = 1.;
            for (unsigned int i = 0, iMax = layers.size(); i < iMax; ++i)
            {
                float absorberThickness((layers[i].inner_thickness - layers[i].sensitive_thickness/2.0 ));

                if (i>0)
                absorberThickness += (layers[i-1].outer_thickness - layers[i-1].sensitive_thickness/2.0);

                if (absorberThickness < std::numeric_limits<float>::epsilon())
                continue;

                if (layerAbsorberThickness > std::numeric_limits<float>::epsilon())
                absorberCorrection = absorberThickness / layerAbsorberThickness;

                break;
            }

            caloHitParameters.m_cellNormalVector = (pCaloHit->Position()[2] > 0) ? pandora::CartesianVector(0, 0, 1) : pandora::CartesianVector(0, 0, -1);
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        void CaloHitCreator::GetBarrelCaloHitProperties( const gar::rec::CaloHit *const pCaloHit,
        const std::vector<gar::geo::LayeredCalorimeterStruct::Layer> &layers, unsigned int barrelSymmetryOrder, PandoraApi::CaloHit::Parameters &caloHitParameters, FloatVector const& normalVector, float &absorberCorrection ) const
        {
            caloHitParameters.m_hitRegion = pandora::BARREL;

            const int physicalLayer(std::min(static_cast<int>(caloHitParameters.m_layer.Get()), static_cast<int>(layers.size()-1)));
            caloHitParameters.m_cellSize0 = layers[physicalLayer].cellSize0;
            caloHitParameters.m_cellSize1 = layers[physicalLayer].cellSize1;

            double thickness = (layers[physicalLayer].inner_thickness+layers[physicalLayer].sensitive_thickness/2.0);
            double nRadLengths = layers[physicalLayer].inner_nRadiationLengths;
            double nIntLengths = layers[physicalLayer].inner_nInteractionLengths;

            double layerAbsorberThickness = (layers[physicalLayer].inner_thickness-layers[physicalLayer].sensitive_thickness/2.0);
            if(physicalLayer>0){
                thickness += (layers[physicalLayer-1].outer_thickness -layers[physicalLayer].sensitive_thickness/2.0);
                nRadLengths += layers[physicalLayer-1].outer_nRadiationLengths;
                nIntLengths += layers[physicalLayer-1].outer_nInteractionLengths;
                layerAbsorberThickness += (layers[physicalLayer-1].outer_thickness -layers[physicalLayer].sensitive_thickness/2.0);
            }

            caloHitParameters.m_cellThickness = thickness;
            caloHitParameters.m_nCellRadiationLengths = nRadLengths;
            caloHitParameters.m_nCellInteractionLengths = nIntLengths;

            if (caloHitParameters.m_nCellRadiationLengths.Get() < std::numeric_limits<float>::epsilon() || caloHitParameters.m_nCellInteractionLengths.Get() < std::numeric_limits<float>::epsilon())
            {
                LOG_ERROR("CaloHitCreator::GetBarrelCaloHitProperties")
                << "Calo hit has 0 radiation length or interaction length: not creating a Pandora calo hit.";
                throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
            }

            absorberCorrection = 1.;
            for (unsigned int i = 0, iMax = layers.size(); i < iMax; ++i)
            {
                float absorberThickness((layers[i].inner_thickness - layers[i].sensitive_thickness/2.0 ));

                if (i>0)
                absorberThickness += (layers[i-1].outer_thickness - layers[i-1].sensitive_thickness/2.0);

                if (absorberThickness < std::numeric_limits<float>::epsilon())
                continue;

                if (layerAbsorberThickness > std::numeric_limits<float>::epsilon())
                absorberCorrection = absorberThickness / layerAbsorberThickness;

                break;
            }

            if (barrelSymmetryOrder > 2)
            {
                const double phi = atan2( pCaloHit->Position()[1], pCaloHit->Position()[0] );

                LOG_DEBUG( "CaloHitCreator::GetBarrelCaloHitProperties" )
                << "This hit does not have any cellIDs set, will use phi-direction for normal vector "
                << " phi:" << std::setw(15) << phi*180/M_PI;

                caloHitParameters.m_cellNormalVector = pandora::CartesianVector( std::cos(phi), std::sin(phi) , 0.0 );
            }
            else
            {
                const float *pCaloHitPosition( pCaloHit->Position() );
                const float phi = std::atan2( pCaloHitPosition[1], pCaloHitPosition[0] );
                caloHitParameters.m_cellNormalVector = pandora::CartesianVector(std::cos(phi), std::sin(phi), 0);
            }
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        int CaloHitCreator::GetNLayersFromEdge(const gar::rec::CaloHit *const pCaloHit) const
        {
            // Calo hit coordinate calculations
            const float barrelMaximumRadius(this->GetMaximumRadius(pCaloHit, m_settings.m_eCalBarrelOuterSymmetry, m_settings.m_eCalBarrelOuterPhi0));
            const float endCapMaximumRadius(this->GetMaximumRadius(pCaloHit, m_settings.m_eCalEndCapInnerSymmetryOrder, m_settings.m_eCalEndCapInnerPhiCoordinate));
            const float caloHitAbsZ(std::fabs(pCaloHit->Position()[2]));

            // Distance from radial outer
            float radialDistanceToEdge(std::numeric_limits<float>::max());

            if (caloHitAbsZ < m_settings.m_eCalBarrelOuterZ)
            {
                radialDistanceToEdge = (m_settings.m_eCalBarrelOuterR - barrelMaximumRadius) / m_eCalBarrelLayerThickness;
            }
            else
            {
                radialDistanceToEdge = (m_settings.m_eCalEndCapOuterR - endCapMaximumRadius) / m_eCalEndCapLayerThickness;
            }

            // Distance from rear of endcap outer
            float rearDistanceToEdge(std::numeric_limits<float>::max());

            if (caloHitAbsZ >= m_settings.m_eCalBarrelOuterZ)
            {
                rearDistanceToEdge = (m_settings.m_eCalEndCapOuterZ - caloHitAbsZ) / m_eCalEndCapLayerThickness;
            }
            else
            {
                const float rearDistance((m_settings.m_eCalBarrelOuterZ - caloHitAbsZ) / m_eCalBarrelLayerThickness);

                if (rearDistance < m_settings.m_layersFromEdgeMaxRearDistance)
                {
                    const float overlapDistance((m_settings.m_eCalEndCapOuterR - endCapMaximumRadius) / m_eCalEndCapLayerThickness);
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
        : m_CaloHitCollection( "" ),
        m_eCalToMip(1.f),
        m_eCalMipThreshold(0.f),
        m_eCalToEMGeV(1.f),
        m_eCalToHadGeVBarrel(1.f),
        m_eCalToHadGeVEndCap(1.f),
        m_maxECalHitHadronicEnergy(10000.f),
        m_nOuterSamplingLayers(3),
        m_layersFromEdgeMaxRearDistance(250.f),
        m_eCalBarrelOuterR(0.f),
        m_eCalBarrelOuterZ(0.f),
        m_eCalBarrelInnerPhi0(0.f),
        m_eCalBarrelOuterPhi0(0.f),
        m_eCalBarrelInnerSymmetry(0.f),
        m_eCalBarrelOuterSymmetry(8),
        m_eCalBarrelNormalVector({0.0, 0.0, 1.0}),
        m_eCalEndCapOuterR(0.f),
        m_eCalEndCapOuterZ(0.f),
        m_eCalEndCapInnerSymmetryOrder(0),
        m_eCalEndCapInnerPhiCoordinate(0.f)
        {
        }
    }
}
