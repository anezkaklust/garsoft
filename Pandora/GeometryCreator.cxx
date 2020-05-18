#include "GeometryCreator.h"

namespace gar {
    namespace gar_pandora {

        GeometryCreator::GeometryCreator(const Settings &settings, const pandora::Pandora *const pPandora) :
        m_settings(settings),
        m_pPandora(*pPandora)
        {
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        GeometryCreator::~GeometryCreator()
        {
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        pandora::StatusCode GeometryCreator::CreateGeometry() const
        {
            return pandora::STATUS_CODE_SUCCESS;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        GeometryCreator::Settings::Settings()
        {
        }
    }
}