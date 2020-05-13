#ifndef ECALSEGMENTATIONGRIDXYALG_H
#define ECALSEGMENTATIONGRIDXYALG_H

#include "Geometry/ChannelMapAlgs/ECALSegmentationAlg.h"
#include "Geometry/GeometryCore.h"

#include <string>
#include <vector>

namespace fhicl{
    class ParameterSet;
}

namespace gar {
    namespace geo {
        namespace seg {

            class ECALSegmentationGridXYAlg: public ECALSegmentationAlg {

            public:
                ECALSegmentationGridXYAlg(fhicl::ParameterSet const& pset);

                ECALSegmentationGridXYAlg(const BitFieldCoder* decoder, fhicl::ParameterSet const& pset);

                ~ECALSegmentationGridXYAlg();

                void reconfigure(fhicl::ParameterSet const& pset) override;

                void Initialize(const gar::geo::GeometryCore& geo) override;

                std::array<double, 3> GetPosition(const gar::geo::GeometryCore& geo, const gar::raw::CellID_t& cID) const override;

                gar::raw::CellID_t GetCellID(const gar::geo::GeometryCore& geo, const unsigned int& det_id, const unsigned int& stave, const unsigned int& module, const unsigned int& layer, const unsigned int& slice, const std::array<double, 3>& localPosition) const override;

                int getIDbyCellID(const gar::raw::CellID_t& cID, const char* id) const override;

                bool isTile(const gar::raw::CellID_t& cID) const override { return true; }

                bool isBarrel(const gar::raw::CellID_t& cID) const override;

                const double& gridSizeX() const { return _gridSizeX; }

                const double& gridSizeY() const { return _gridSizeY; }

                const double& offsetX() const { return _offsetX; }

                const double& offsetY() const { return _offsetY; }

                const std::string& fieldNameX() const { return _xId; }

                const std::string& fieldNameY() const { return _yId; }

                const unsigned int& nLayers() const { return _nLayers; }

                void setGridSizeX(double cellSize) { _gridSizeX = cellSize; }

                void setGridSizeY(double cellSize) { _gridSizeY = cellSize; }

                void setOffsetX(double offset) { _offsetX = offset; }

                void setOffsetY(double offset) { _offsetY = offset; }

                void setFieldNameX(const std::string& fieldName) { _xId = fieldName; }

                void setFieldNameY(const std::string& fieldName) { _yId = fieldName; }

                void setLayerDimXY(const double& dimX, const double& dimY) const { _layer_dim_X = dimX; _layer_dim_Y = dimY; }

            protected:

                void PrintParameters() const override;

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
                /// number of layers
                unsigned int _nLayers;

                /// layer dimension in X
                mutable double _layer_dim_X;
                /// layer dimension in Y
                mutable double _layer_dim_Y;
            };
        }
    }
}

#endif
