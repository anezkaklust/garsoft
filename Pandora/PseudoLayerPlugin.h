#ifndef PSEUDOLAYER_PLUGIN_H
#define PSEUDOLAYER_PLUGIN_H 1

#include "Objects/CartesianVector.h"

//Pandora
#include "Plugins/PseudoLayerPlugin.h"

// Pandora SDK
#include "Api/PandoraApi.h"

namespace gar {
  namespace gar_pandora {

    class PseudoLayerPlugin : public pandora::PseudoLayerPlugin
    {
    public:

      class Settings 
      {
      public:
        float    m_innerZ {} ;
        float    m_layerThickness {} ;
      };

      PseudoLayerPlugin( const Settings &settings );

    private:

      void InitializeGeometry() {}

      unsigned int GetPseudoLayer(const pandora::CartesianVector &positionVector) const override ;
      unsigned int GetPseudoLayerAtIp() const override ;

      pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

      float         m_EndCapInnerZ; //Z coordinate of start, current scenario = 0.0
      float         m_EndCapLayerThickness; // Layer Thickness of calo

    };
  }
}

#endif // #ifndef PSEUDOLAYER_PLUGIN_H