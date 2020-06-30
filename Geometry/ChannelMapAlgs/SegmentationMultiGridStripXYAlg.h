#ifndef SEGMENTATIONMULTIGRIDSTRIPXYALG_H
#define SEGMENTATIONMULTIGRIDSTRIPXYALG_H

#include "Geometry/ChannelMapAlgs/SegmentationAlg.h"

#include <string>
#include <vector>
#include <array>

namespace fhicl{
    class ParameterSet;
}

namespace gar {
    namespace geo {
        namespace seg {

            class SegmentationMultiGridStripXYAlg: public SegmentationAlg {

            public:
                SegmentationMultiGridStripXYAlg(fhicl::ParameterSet const& pset);

                SegmentationMultiGridStripXYAlg(const BitFieldCoder* decoder, fhicl::ParameterSet const& pset);

                ~SegmentationMultiGridStripXYAlg();

                void reconfigure(fhicl::ParameterSet const& pset) override;

                void Initialize(const gar::geo::GeometryCore& geo) override;

                std::array<double, 3> GetPosition(const gar::geo::GeometryCore& geo, const gar::raw::CellID_t& cID) const override;

                gar::raw::CellID_t GetCellID(const gar::geo::GeometryCore& geo, const unsigned int& det_id, const unsigned int& stave, const unsigned int& module, const unsigned int& layer, const unsigned int& slice, const std::array<double, 3>& localPosition) const override;

                bool isTile(const gar::raw::CellID_t& cID) const override;

                bool isBarrel(const gar::raw::CellID_t& cID) const override;

                double getStripLength(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const gar::raw::CellID_t& cID) const override;

                std::pair<TVector3, TVector3> getStripEnds(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const gar::raw::CellID_t& cID) const override;

                std::pair<float, float> CalculateLightPropagation(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const gar::raw::CellID_t& cID) const override;

                std::array<double, 3> ReconstructStripHitPosition(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const float &xlocal, const gar::raw::CellID_t& cID) const override;

                const double& gridSizeX() const { return _gridSizeX; }

                const double& gridSizeY() const { return _gridSizeY; }

                const double& stripSizeX() const { return _stripSizeX; }

                const double& stripSizeY() const { return _stripSizeY; }

                const double& offsetX() const { return _offsetX; }

                const double& offsetY() const { return _offsetY; }

                const double& layerDimX() const { return _layer_dim_X; }

                const double& layerDimY() const { return _layer_dim_Y; }

                const std::string& fieldNameX() const { return _xId; }

                const std::string& fieldNameY() const { return _yId; }

                const std::string& fieldNameLayer() const { return _layerId; }

                const std::string& fieldNameSlice() const { return _sliceId; }

                const unsigned int& nLayers() const { return _nLayers; }

                void setGridSizeX(double cellSize) { _gridSizeX = cellSize; }

                void setGridSizeY(double cellSize) { _gridSizeY = cellSize; }

                void setStripSizeX(double stripSize) { _stripSizeX = stripSize; }

                void setStripSizeY(double stripSize) { _stripSizeY = stripSize; }

                void setOffsetX(double offset) { _offsetX = offset; }

                void setOffsetY(double offset) { _offsetY = offset; }

                void setFieldNameX(const std::string& fieldName) { _xId = fieldName; }

                void setFieldNameY(const std::string& fieldName) { _yId = fieldName; }

                void setFieldNameLayer(const std::string& fieldName) { _layerId = fieldName; }

                void setFieldNameSlice(const std::string& fieldName) { _sliceId = fieldName; }

                void setLayerDimXY(const double& dimX, const double& dimY) const { _layer_dim_X = dimX; _layer_dim_Y = dimY; }

            protected:

                void PrintParameters() const override;

                std::array<std::vector<unsigned int>, 2> TokenizeLayerVectors(std::vector<std::string> grid) const;

                /// the grid size in X
                double _gridSizeX;
                /// the coordinate offset in X
                double _offsetX;
                /// the grid size in Y
                double _gridSizeY;
                /// the coordinate offset in Y
                double _offsetY;
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
                /// number of layers
                unsigned int _nLayers;
                /// active slice number
                unsigned int _active_slice;
                /// the layers (start::end) for grid in Barrel
                std::vector<std::string> _gridBarrelLayers;
                /// the layers (start::end) for strips in Barrel
                std::vector<std::string> _stripBarrelLayers;
                /// the layers (start::end) for grid in Endcap
                std::vector<std::string> _gridEndcapLayers;
                /// the layers (start::end) for strips in Endcap
                std::vector<std::string> _stripEndcapLayers;
                /// layer dimension in X
                mutable double _layer_dim_X;
                /// layer dimension in Y
                mutable double _layer_dim_Y;
                /// are the strips on the same layer?
                bool _OnSameLayer;
            };
        }
    }
}

#endif
