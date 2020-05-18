#include "BFieldPlugin.h"

#include "Helpers/XmlHelper.h"

namespace gar {
    namespace gar_pandora {

        BFieldPlugin::BFieldPlugin()
        {
            /* nop */
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        float BFieldPlugin::GetBField(const pandora::CartesianVector &positionVector) const
        {
            double bfield[3] = {0.5, 0., 0.};
            return bfield[0];
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode BFieldPlugin::Initialize()
        {
            /* nop */
            return pandora::STATUS_CODE_SUCCESS;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode BFieldPlugin::ReadSettings(const pandora::TiXmlHandle /*xmlHandle*/)
        {
            /* nop */
            return pandora::STATUS_CODE_SUCCESS;
        }
    }
}