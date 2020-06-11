#include "Geometry/ChannelMapAlgs/SegmentationStripXAlg.h"

#include "fhiclcpp/ParameterSet.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geometry/GeometryCore.h"

#include <algorithm>
#include <sstream>
#include <iostream>
#include <cmath>

namespace gar {
    namespace geo {
        namespace seg {

            //----------------------------------------------------------------------------
            /// Default constructor used by derived classes passing the encoding string
            SegmentationStripXAlg::SegmentationStripXAlg(fhicl::ParameterSet const& pset)
            : SegmentationAlg(pset)
            {
                _type = "StripX";
                _description = "Cartesian segmentation in the local XY-plane with strip along the Y direction";

                std::cout << "######### gar::geo::seg::SegmentationStripXAlg() " << std::endl ;

                this->reconfigure(pset);
            }

            //----------------------------------------------------------------------------
            /// Default constructor used by derived classes passing an existing decoder
            SegmentationStripXAlg::SegmentationStripXAlg(const BitFieldCoder* decode, fhicl::ParameterSet const& pset)
            : SegmentationAlg(decode, pset)
            {
                _type = "StripX";
                _description = "Cartesian segmentation in the local XY-plane with strip along the Y direction";

                std::cout << "######### gar::geo::seg::SegmentationStripXAlg() " << std::endl ;

                this->reconfigure(pset);
            }

            //----------------------------------------------------------------------------
            SegmentationStripXAlg::~SegmentationStripXAlg()
            {
            }

            //----------------------------------------------------------------------------
            void SegmentationStripXAlg::reconfigure(fhicl::ParameterSet const& pset)
            {
                _stripSizeX = pset.get<double>("strip_size_x");
                _encoding = pset.get<std::string>("cellEncoding");

                _xId = pset.get<std::string>("identifier_x");
                _yId = pset.get<std::string>("identifier_y");
                _layerId = pset.get<std::string>("identifier_layer");
                _sliceId = pset.get<std::string>("identifier_slice");

                _nLayers = pset.get<unsigned int>("nlayers");

                _frac = 1./3.;

                this->PrintParameters();

                return;
            }

            //----------------------------------------------------------------------------
            void SegmentationStripXAlg::Initialize(const gar::geo::GeometryCore& geo)
            {

            }

            //----------------------------------------------------------------------------
            std::array<double, 3> SegmentationStripXAlg::GetPosition(const gar::geo::GeometryCore& geo, const gar::raw::CellID_t& cID) const
            {
                /* no op */
                return std::array<double, 3>{ 0., 0., 0. };
            }

            //----------------------------------------------------------------------------
            /// determine the cell ID based on the position
            gar::raw::CellID_t SegmentationStripXAlg::GetCellID(const gar::geo::GeometryCore& geo, const unsigned int& det_id, const unsigned int& stave, const unsigned int& module, const unsigned int& layer, const unsigned int& slice, const std::array<double, 3>& localPosition) const
            {
                gar::raw::CellID_t cID = 0;

                _decoder->set(cID, "system", det_id);
                _decoder->set(cID, "module", module);
                _decoder->set(cID, "stave", stave);
                _decoder->set(cID, "layer", layer);
                _decoder->set(cID, "slice", slice);

                double localX = localPosition[0];
                double localY = localPosition[1];

                //Segmentation in X
                int nCellsX = int(_layer_dim_X / _stripSizeX);
                int nCellsY = 1;

                int _cellIndexX = int ( localX / ( _layer_dim_X / nCellsX ) );
                int _cellIndexY = int ( localY / ( _layer_dim_Y / nCellsY ) );

                _decoder->set(cID, _xId, _cellIndexX);
                _decoder->set(cID, _yId, _cellIndexY);

                return cID;
            }

            //----------------------------------------------------------------------------
            void SegmentationStripXAlg::PrintParameters() const
            {
                std::cout << "cell encoding: " << _encoding << std::endl;
                std::cout << "identifier_x: " << _xId << std::endl;
                std::cout << "identifier_y: " << _yId << std::endl;
                std::cout << "strip_size_x: " << _stripSizeX << " cm" << std::endl;

                std::cout << "identifier_layer: " << _layerId << std::endl;
                std::cout << "identifier_slice: " << _sliceId << std::endl;
            }

            //----------------------------------------------------------------------------
            bool SegmentationStripXAlg::isTile(const gar::raw::CellID_t& cID) const
            {
                return false;
            }

            //----------------------------------------------------------------------------
            bool SegmentationStripXAlg::isBarrel(const gar::raw::CellID_t& cID) const
            {
                return true;
            }

            //----------------------------------------------------------------------------
            double SegmentationStripXAlg::getStripLength(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const gar::raw::CellID_t& cID) const
            {
                /* no op */
                return 0.;
            }

            //----------------------------------------------------------------------------
            std::pair<TVector3, TVector3> SegmentationStripXAlg::getStripEnds(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const gar::raw::CellID_t& cID) const
            {
                /* no op */
                return std::make_pair(TVector3(), TVector3());
            }

            //----------------------------------------------------------------------------
            std::pair<float, float> SegmentationStripXAlg::CalculateLightPropagation(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const gar::raw::CellID_t& cID) const
            {
                /* no op */
                return std::make_pair(0., 0.);
            }

            //----------------------------------------------------------------------------
            std::array<double, 3> SegmentationStripXAlg::ReconstructStripHitPosition(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const float &xlocal, const gar::raw::CellID_t& cID) const
            {
                /* no op */
                return std::array<double, 3>{ 0., 0., 0.};
            }

        }//seg
    } // geo
} //gar
