#ifndef CALO_HIT_CREATOR_H
#define CALO_HIT_CREATOR_H 1

#include "art/Framework/Principal/Event.h"

#include "Geometry/Geometry.h"
#include "ReconstructionDataProducts/CaloHit.h"

#include "Api/PandoraApi.h"


namespace gar {
    namespace gar_pandora {

        typedef std::vector<gar::rec::CaloHit *> CalorimeterHitVector;

        class CaloHitCreator
        {
        public:
            typedef std::vector<std::string> StringVector;
            typedef std::vector<float> FloatVector;

            class Settings
            {
            public:
                Settings();

                StringVector    m_CaloHitCollections;               ///< The calorimeter hit collections
                StringVector    m_MuIDCaloHitCollections;           ///< The muon id calorimeter hit collections
                
                float           m_eCalToMip;                            ///< The calibration from deposited ECal energy to mip
                float           m_hCalToMip;                            ///< The calibration from deposited HCal energy to mip
                float           m_muonToMip;                            ///< The calibration from deposited Muon energy to mip
                float           m_eCalMipThreshold;                     ///< Threshold for creating calo hits in the ECal, units mip
                float           m_hCalMipThreshold;                     ///< Threshold for creating calo hits in the HCal, units mip
                float           m_muonMipThreshold;                     ///< Threshold for creating calo hits in the HCal, units mip

                float           m_eCalToEMGeV;                          ///< The calibration from deposited ECal energy to EM energy
                float           m_eCalToHadGeVBarrel;                   ///< The calibration from deposited ECal barrel energy to hadronic energy
                float           m_eCalToHadGeVEndCap;                   ///< The calibration from deposited ECal endcap energy to hadronic energy
                float           m_hCalToEMGeV;                          ///< The calibration from deposited HCal energy to EM energy
                float           m_hCalToHadGeV;                         ///< The calibration from deposited HCal energy to hadronic energy

                float           m_maxHCalHitHadronicEnergy;             ///< The maximum hadronic energy allowed for a single hcal hit
                int             m_nOuterSamplingLayers;                 ///< Number of layers from edge for hit to be flagged as an outer layer hit
                float           m_layersFromEdgeMaxRearDistance;        ///< Maximum number of layers from candidate outer layer hit to rear of detector

                int             m_hCalEndCapInnerSymmetryOrder;         ///< HCal end cap inner symmetry order
                float           m_hCalEndCapInnerPhiCoordinate;         ///< HCal end cap inner phi coordinate

                // For Strip Splitting method.
                int             m_stripSplittingOn;                     ///< To use SSA, this should be true (default is false)
                
                float                         m_eCalBarrelOuterZ;                 ///< ECal barrel outer z coordinate
                float                         m_hCalBarrelOuterZ;                 ///< HCal barrel outer z coordinate
                float                         m_muonBarrelOuterZ;                 ///< Muon barrel outer z coordinate
                float                         m_coilOuterR;                       ///< Coil outer r coordinate

                float                         m_eCalBarrelInnerPhi0;              ///< ECal barrel inner phi0 coordinate
                unsigned int                  m_eCalBarrelInnerSymmetry;          ///< ECal barrel inner symmetry order
                float                         m_hCalBarrelInnerPhi0;              ///< HCal barrel inner phi0 coordinate
                unsigned int                  m_hCalBarrelInnerSymmetry;          ///< HCal barrel inner symmetry order
                float                         m_muonBarrelInnerPhi0;              ///< Muon barrel inner phi0 coordinate
                unsigned int                  m_muonBarrelInnerSymmetry;          ///< Muon barrel inner symmetry order

                float                         m_hCalEndCapOuterR;                 ///< HCal endcap outer r coordinate
                float                         m_hCalEndCapOuterZ;                 ///< HCal endcap outer z coordinate
                float                         m_hCalBarrelOuterR;                 ///< HCal barrel outer r coordinate
                float                         m_hCalBarrelOuterPhi0;              ///< HCal barrel outer phi0 coordinate
                unsigned int                  m_hCalBarrelOuterSymmetry;          ///< HCal barrel outer symmetry order

                FloatVector m_eCalBarrelNormalVector;
                FloatVector m_hCalBarrelNormalVector;
                FloatVector m_muonBarrelNormalVector;
            };

            CaloHitCreator(const Settings &settings, const pandora::Pandora *const pPandora);

            ~CaloHitCreator();

            pandora::StatusCode CreateCaloHits(const art::Event *const pEvent);

            const CalorimeterHitVector &GetCalorimeterHitVector() const;

            void Reset();

        private:

            pandora::StatusCode CreateECalCaloHits(const art::Event *const pEvent);

            void GetCommonCaloHitProperties(const gar::rec::CaloHit *const pCaloHit, PandoraApi::CaloHit::Parameters &caloHitParameters) const;

            int GetNLayersFromEdge(const  gar::rec::CaloHit *const pCaloHit) const;

            float GetMaximumRadius(const  gar::rec::CaloHit *const pCaloHit, const unsigned int symmetryOrder, const float phi0) const;

            const Settings                      m_settings;                         ///< The calo hit creator settings

            const pandora::Pandora &            m_pandora;                          ///< Reference to the pandora object to create calo hits

            float                               m_hCalBarrelLayerThickness;         ///< HCal barrel layer thickness
            float                               m_hCalEndCapLayerThickness;         ///< HCal endcap layer thickness

            CalorimeterHitVector                m_calorimeterHitVector;             ///< The calorimeter hit vector

            const geo::GeometryCore*            fGeo; //Geometry Manager
        };

        //------------------------------------------------------------------------------------------------------------------------------------------

        inline const CalorimeterHitVector &CaloHitCreator::GetCalorimeterHitVector() const
        {
            return m_calorimeterHitVector;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        inline void CaloHitCreator::Reset()
        {
            m_calorimeterHitVector.clear();
        }

    }
}

#endif // #ifndef CALO_HIT_CREATOR_H