#ifndef GEOMETRY_CREATOR_H
#define GEOMETRY_CREATOR_H 1

#include "Api/PandoraApi.h"

#include "Geometry/Geometry.h"

namespace gar {
    namespace gar_pandora {

        class GeometryCreator
        {
        public:

            class Settings
            {
            public:
                Settings();        
            };

            GeometryCreator(const Settings &settings, const pandora::Pandora *const pPandora);

            ~GeometryCreator();

            pandora::StatusCode CreateGeometry() const;

        private:
            typedef std::map<pandora::SubDetectorType, PandoraApi::Geometry::SubDetector::Parameters> SubDetectorTypeMap;
            typedef std::map<std::string, PandoraApi::Geometry::SubDetector::Parameters> SubDetectorNameMap;

            const Settings          m_settings;                     ///< The geometry creator settings
            const pandora::Pandora &m_pPandora;                     ///< Address of the pandora object to create the geometry
        };
    }
}

#endif // #ifndef GEOMETRY_CREATOR_H