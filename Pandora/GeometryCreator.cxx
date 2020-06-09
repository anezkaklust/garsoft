#include "GeometryCreator.h"

#include "CoreUtils/ServiceUtil.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <exception>

namespace gar {
    namespace gar_pandora {

        GeometryCreator::GeometryCreator(const Settings &settings, const pandora::Pandora *const pPandora) :
        m_settings(settings),
        m_pPandora(*pPandora)
        {
            fGeo = gar::providerFrom<geo::Geometry>();
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        GeometryCreator::~GeometryCreator()
        {
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode GeometryCreator::CreateGeometry() const
        {
            try
            {
                SubDetectorTypeMap subDetectorTypeMap;
                this->SetMandatorySubDetectorParameters(subDetectorTypeMap);

                SubDetectorNameMap subDetectorNameMap;
                this->SetAdditionalSubDetectorParameters(subDetectorNameMap);

                std::string detectorName = "NDGAr";

                LOG_DEBUG( "GeometryCreator::CreateGeometry()" ) << "Creating geometry for detector " << detectorName << std::endl;

                for (SubDetectorTypeMap::const_iterator iter = subDetectorTypeMap.begin(), iterEnd = subDetectorTypeMap.end(); iter != iterEnd; ++iter)
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::SubDetector::Create(m_pPandora, iter->second));

                for (SubDetectorNameMap::const_iterator iter = subDetectorNameMap.begin(), iterEnd = subDetectorNameMap.end(); iter != iterEnd; ++iter){
                    LOG_DEBUG( "GeometryCreator::CreateGeometry()" ) << "Creating geometry for additional subdetector " << iter->first << std::endl;
                    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::SubDetector::Create(m_pPandora, iter->second));
                }
            }
            catch (std::exception &exception)
            {
                LOG_ERROR( "GeometryCreator::CreateGeometry()" ) << "Failure in pandora geometry creator, exception: " << exception.what() << std::endl;
                throw exception;
            }

            return pandora::STATUS_CODE_SUCCESS;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        void GeometryCreator::SetMandatorySubDetectorParameters(SubDetectorTypeMap &subDetectorTypeMap) const
        {
            PandoraApi::Geometry::SubDetector::Parameters eCalBarrelParameters, eCalEndCapParameters,
            muonBarrelParameters, muonEndCapParameters;

            this->SetDefaultSubDetectorParameters("ECalBarrel", pandora::ECAL_BARREL, eCalBarrelParameters);
            this->SetDefaultSubDetectorParameters("ECalEndcap", pandora::ECAL_BARREL, eCalBarrelParameters);

            subDetectorTypeMap[pandora::ECAL_BARREL] = eCalBarrelParameters;
            subDetectorTypeMap[pandora::ECAL_ENDCAP] = eCalEndCapParameters;

            PandoraApi::Geometry::SubDetector::Parameters trackerParameters;

            trackerParameters.m_subDetectorName = "Tracker";
            trackerParameters.m_subDetectorType = pandora::INNER_TRACKER;
            trackerParameters.m_innerRCoordinate = 0.f; // inner R of TPC
            trackerParameters.m_innerXCoordinate = 0.f;
            trackerParameters.m_innerPhiCoordinate = 0.f;
            trackerParameters.m_innerSymmetryOrder = 0;
            trackerParameters.m_outerRCoordinate = fGeo->TPCRadius(); // outer R of TPC
            trackerParameters.m_outerXCoordinate = fGeo->TPCLength() / 2.; // outer X of TPC
            trackerParameters.m_outerPhiCoordinate = 0.f;
            trackerParameters.m_outerSymmetryOrder = 0;
            trackerParameters.m_isMirroredInX = true;
            trackerParameters.m_nLayers = 0;
            subDetectorTypeMap[pandora::INNER_TRACKER] = trackerParameters;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        void GeometryCreator::SetAdditionalSubDetectorParameters(SubDetectorNameMap &subDetectorNameMap) const
        {

        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        void GeometryCreator::SetDefaultSubDetectorParameters(const std::string &subDetectorName, const pandora::SubDetectorType subDetectorType, PandoraApi::Geometry::SubDetector::Parameters &parameters) const
        {
            parameters.m_subDetectorName = subDetectorName;
            parameters.m_subDetectorType = subDetectorType;
            parameters.m_innerRCoordinate = inputParameters.extent[0]/dd4hep::mm;
            parameters.m_innerXCoordinate = inputParameters.extent[2]/dd4hep::mm;
            parameters.m_innerPhiCoordinate = inputParameters.inner_phi0;
            parameters.m_innerSymmetryOrder = inputParameters.inner_symmetry;
            parameters.m_outerRCoordinate = inputParameters.extent[1]/dd4hep::mm;
            parameters.m_outerXCoordinate = inputParameters.extent[3]/dd4hep::mm;
            parameters.m_outerPhiCoordinate = inputParameters.outer_phi0;
            parameters.m_outerSymmetryOrder = inputParameters.outer_symmetry;
            parameters.m_isMirroredInX = true;
            parameters.m_nLayers = layers.size();

            for (size_t i = 0; i< layers.size(); i++)
            {
                const dd4hep::rec::LayeredCalorimeterStruct::Layer & theLayer = layers.at(i);

                PandoraApi::Geometry::LayerParameters layerParameters;

                double totalNumberOfRadLengths = theLayer.inner_nRadiationLengths;
                double totalNumberOfIntLengths = theLayer.inner_nInteractionLengths;

                if(i>0){
                    //Add the numbers from previous layer's outer side
                    totalNumberOfRadLengths += layers.at(i-1).outer_nRadiationLengths;
                    totalNumberOfIntLengths += layers.at(i-1).outer_nInteractionLengths;
                }

                layerParameters.m_closestDistanceToIp = (theLayer.distance+theLayer.inner_thickness)/dd4hep::mm; //Distance to center of sensitive element
                layerParameters.m_nRadiationLengths = totalNumberOfRadLengths;
                layerParameters.m_nInteractionLengths = totalNumberOfIntLengths;

                parameters.m_layerParametersVector.push_back(layerParameters);
            }
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        GeometryCreator::Settings::Settings()
        {
        }
    }
}
