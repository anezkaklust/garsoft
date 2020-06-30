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
            trackerParameters.m_innerZCoordinate = 0.f;
            trackerParameters.m_innerPhiCoordinate = 0.f;
            trackerParameters.m_innerSymmetryOrder = 0;
            trackerParameters.m_outerRCoordinate = fGeo->TPCRadius(); // outer R of TPC
            trackerParameters.m_outerZCoordinate = fGeo->TPCLength() / 2.; // outer X of TPC
            trackerParameters.m_outerPhiCoordinate = 0.f;
            trackerParameters.m_outerSymmetryOrder = 0;
            trackerParameters.m_isMirroredInZ = true;
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
            parameters.m_innerRCoordinate = 0.f;
            parameters.m_innerZCoordinate = 0.f;
            parameters.m_innerPhiCoordinate = 0.f;
            parameters.m_innerSymmetryOrder = 8;
            parameters.m_outerRCoordinate = fGeo->TPCRadius(); // outer R of TPC
            parameters.m_outerZCoordinate = fGeo->TPCLength() / 2.; // outer X of TPC
            parameters.m_outerPhiCoordinate = 0.f;
            parameters.m_outerSymmetryOrder = 0;
            parameters.m_isMirroredInZ = true;
            parameters.m_nLayers = 60;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        GeometryCreator::Settings::Settings()
        {
        }
    }
}
