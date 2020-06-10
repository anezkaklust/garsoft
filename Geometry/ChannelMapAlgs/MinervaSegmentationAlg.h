#ifndef MINERVASEGMENTATIONALG_H
#define MINERVASEGMENTATIONALG_H

#include "Geometry/ChannelMapAlgs/ECALSegmentationAlg.h"

#include <string>
#include <vector>
#include <array>

namespace fhicl{
    class ParameterSet;
}

namespace gar {
    namespace geo {
        namespace seg {

            class MinervaSegmentationAlg: public ECALSegmentationAlg {

            public:
                MinervaSegmentationAlg(fhicl::ParameterSet const& pset);

                MinervaSegmentationAlg(const BitFieldCoder* decoder, fhicl::ParameterSet const& pset);

                ~MinervaSegmentationAlg();

                void reconfigure(fhicl::ParameterSet const& pset) override;

                void Initialize(const gar::geo::GeometryCore& geo) override;

                std::array<double, 3> GetPosition(const gar::geo::GeometryCore& geo, const gar::raw::CellID_t& cID) const override;

                gar::raw::CellID_t GetCellID(const gar::geo::GeometryCore& geo, const unsigned int& det_id, const unsigned int& stave, const unsigned int& module, const unsigned int& layer, const unsigned int& slice, const std::array<double, 3>& localPosition) const override;

                const double& stripSizeX() const { return _stripSizeX; }

                const double& stripSizeY() const { return _stripSizeY; }

                const double& layerDimX() const { return _layer_dim_X; }

                const double& layerDimY() const { return _layer_dim_Y; }

                const double& layerDimZ() const { return _layer_dim_Z; }

                const std::string& fieldNameX() const { return _xId; }

                const std::string& fieldNameY() const { return _yId; }

                const std::string& fieldNameLayer() const { return _layerId; }

                const std::string& fieldNameSlice() const { return _sliceId; }

                void setStripSizeX(double stripSize) { _stripSizeX = stripSize; }

                void setStripSizeY(double stripSize) { _stripSizeY = stripSize; }

                void setFieldNameX(const std::string& fieldName) { _xId = fieldName; }

                void setFieldNameY(const std::string& fieldName) { _yId = fieldName; }

                void setFieldNameLayer(const std::string& fieldName) { _layerId = fieldName; }

                void setFieldNameSlice(const std::string& fieldName) { _sliceId = fieldName; }

                bool isTile(const gar::raw::CellID_t& cID) const override { /* no op */ return false; }

                bool isBarrel(const gar::raw::CellID_t& cID) const override { /* no op */ return true; }

                void setLayerDimXY(const double& dimX, const double& dimY) const override { /* no op */  return; }

            protected:

                void PrintParameters() const override;

                /// the field name used for X
                std::string _xId;
                /// the field name used for Y
                std::string _yId;
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
                /// layer dimension in X
                double _layer_dim_X;
                /// layer dimension in Y
                double _layer_dim_Y;
                /// layer dimension in Z
                double _layer_dim_Z;
            };
        }
    }
}

#endif
