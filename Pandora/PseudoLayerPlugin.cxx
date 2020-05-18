#include "PseudoLayerPlugin.h"

//Pandora
#include "Pandora/AlgorithmHeaders.h"
#include "Helpers/XmlHelper.h"

using namespace pandora;

namespace gar {
  namespace gar_pandora {

    PseudoLayerPlugin::PseudoLayerPlugin( const Settings &settings ) :
        m_EndCapInnerZ(settings.m_innerZ),
        m_EndCapLayerThickness(settings.m_layerThickness)
    {
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    unsigned int PseudoLayerPlugin::GetPseudoLayer(const pandora::CartesianVector &positionVector) const
    {
      //TODO
      return 0;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    unsigned int PseudoLayerPlugin::GetPseudoLayerAtIp() const
    {
      return 0;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    StatusCode PseudoLayerPlugin::ReadSettings(const TiXmlHandle xmlHandle)
    {
        return STATUS_CODE_SUCCESS;
    }

  }
}