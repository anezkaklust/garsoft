#ifndef BFIELD_PLUGIN_H
#define BFIELD_PLUGIN_H 1

#include "Plugins/BFieldPlugin.h"
#include "Objects/CartesianVector.h"
#include "nug4/MagneticField/MagneticField.h"

namespace gar {
    namespace gar_pandora {

        class BFieldPlugin : public pandora::BFieldPlugin
        {
        public:
            BFieldPlugin(art::ServiceHandle<mag::MagneticField> &magFieldService);
            float GetBField(const pandora::CartesianVector &positionVector) const;

        private:
            pandora::StatusCode Initialize();
            pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

            art::ServiceHandle<mag::MagneticField> fieldService;
        };
    }
}

#endif // #ifndef BFIELD_PLUGIN_H
