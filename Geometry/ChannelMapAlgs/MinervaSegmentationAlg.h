#ifndef MINERVASEGMENTATIONALG_H
#define MINERVASEGMENTATIONALG_H

#include "Geometry/ChannelMapAlgs/SegmentationAlg.h"

#include "SimulationDataProducts/CaloDeposit.h"
#include "RawDataProducts/CaloRawDigit.h"

#include <string>
#include <vector>
#include <array>

namespace fhicl{
    class ParameterSet;
}

namespace gar {
    namespace geo {
        namespace seg {

            class MinervaSegmentationAlg: public SegmentationAlg {

            public:
                MinervaSegmentationAlg(fhicl::ParameterSet const& pset);

                MinervaSegmentationAlg(const BitFieldCoder* decoder, fhicl::ParameterSet const& pset);

                ~MinervaSegmentationAlg();

                void reconfigure(fhicl::ParameterSet const& pset) override;

                void Initialize(const gar::geo::GeometryCore& geo) override;

                std::array<double, 3> GetPosition(const gar::geo::GeometryCore& geo, const gar::raw::CellID_t& cID) const override;

                gar::raw::CellID_t GetCellID(const gar::geo::GeometryCore& geo, const unsigned int& det_id, const unsigned int& stave, const unsigned int& module, const unsigned int& layer, const unsigned int& slice, const std::array<double, 3>& localPosition) const override;

                const double& stripSizeX() const override { return _stripSizeX; }

                const double& stripSizeY() const { return _stripSizeY; }

                const double& layerDimX() const { return _layer_dim_X; }

                const double& layerDimY() const { return _layer_dim_Y; }

                const std::string& fieldNameX() const { return _xId; }

                const std::string& fieldNameY() const { return _yId; }

                const std::string& fieldNameZ() const { return _zId; }

                const std::string& fieldNameLayer() const { return _layerId; }

                const std::string& fieldNameSlice() const { return _sliceId; }

                const unsigned int& nPlanes() const override { return _nPlanes; }

                void setStripSizeX(double stripSize) { _stripSizeX = stripSize; }

                void setStripSizeY(double stripSize) { _stripSizeY = stripSize; }

                void setFieldNameX(const std::string& fieldName) { _xId = fieldName; }

                void setFieldNameY(const std::string& fieldName) { _yId = fieldName; }

                void setFieldNameLayer(const std::string& fieldName) { _layerId = fieldName; }

                void setFieldNameSlice(const std::string& fieldName) { _sliceId = fieldName; }

                bool isTile(const gar::raw::CellID_t& cID) const override { /* no op */ return false; }

                bool isBarrel(const gar::raw::CellID_t& cID) const override { /* no op */ return true; }

                void setLayerDimXY(const double& dimX, const double& dimY) const override { _layer_dim_X = dimX; _layer_dim_Y = dimY; }

                void setVariables(const double& innerangle, const double &endcapsidelength) const override { /* no op */ }

                void AddHitsMinerva(std::map< gar::raw::CellID_t, std::vector<gar::sdp::CaloDeposit> > &m_Deposits, std::vector<gar::sdp::CaloDeposit> &fDeposits) const;

            protected:

                gar::raw::CellID_t GetComplementaryCellID(gar::raw::CellID_t cellID, unsigned int comp) const;

                void PrintParameters() const override;

                /// the field name used for X
                std::string _xId;
                /// the field name used for Y
                std::string _yId;
                /// the field name used for Z
                std::string _zId;
                /// the field name used for layer
                std::string _layerId;
                /// the field name used for slice
                std::string _sliceId;
                /// the encoding string
                std::string _encoding;
                /// the strip size in X
                double _stripSizeX;
                /// the strip size in Y
                double _stripSizeY;
                /// fraction of tiles to remove at the edge
                double _frac;
                /// number of planes
                unsigned int _nPlanes;
                /// layer dimension in X
                mutable double _layer_dim_X;
                /// layer dimension in Y
                mutable double _layer_dim_Y;
            };
        }
    }
}

#endif
