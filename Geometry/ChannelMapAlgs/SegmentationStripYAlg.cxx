#include "Geometry/ChannelMapAlgs/SegmentationStripYAlg.h"

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
            SegmentationStripYAlg::SegmentationStripYAlg(fhicl::ParameterSet const& pset)
            : SegmentationAlg(pset)
            {
                _type = "StripY";
                _description = "Cartesian segmentation in the local XY-plane with strip along the Y direction";

                std::cout << "######### gar::geo::seg::SegmentationStripYAlg() " << std::endl ;

                this->reconfigure(pset);
            }

            //----------------------------------------------------------------------------
            /// Default constructor used by derived classes passing an existing decoder
            SegmentationStripYAlg::SegmentationStripYAlg(const BitFieldCoder* decode, fhicl::ParameterSet const& pset)
            : SegmentationAlg(decode, pset)
            {
                _type = "StripY";
                _description = "Cartesian segmentation in the local XY-plane with strip along the Y direction";

                std::cout << "######### gar::geo::seg::SegmentationStripYAlg() " << std::endl ;

                this->reconfigure(pset);
            }

            //----------------------------------------------------------------------------
            SegmentationStripYAlg::~SegmentationStripYAlg()
            {
            }

            //----------------------------------------------------------------------------
            void SegmentationStripYAlg::reconfigure(fhicl::ParameterSet const& pset)
            {
                _stripSizeY = pset.get<double>("strip_size_y");
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
            void SegmentationStripYAlg::Initialize(const gar::geo::GeometryCore& geo)
            {

            }

            //----------------------------------------------------------------------------
            std::array<double, 3> SegmentationStripYAlg::GetPosition(const gar::geo::GeometryCore& geo, const gar::raw::CellID_t& cID) const
            {
                std::array<double, 3> cellPosition;
                
                int cellIndexX = _decoder->get(cID,_xId);
                int cellIndexY = _decoder->get(cID,_yId);

                //Segmentation in Y
                int nCellsX = 1;
                // int nCellsY = int(_layer_dim_Y / _stripSizeY);

                cellPosition[0] = ( cellIndexX + 0.5 ) * (_layer_dim_X / nCellsX ) - (_layer_dim_X / 2.);
                cellPosition[1] = ( cellIndexY + 0.5 ) * _stripSizeY;
                cellPosition[2] = 0.;

                return cellPosition;
            }

            //----------------------------------------------------------------------------
            /// determine the cell ID based on the position
            gar::raw::CellID_t SegmentationStripYAlg::GetCellID(const gar::geo::GeometryCore& geo, const unsigned int& det_id, const unsigned int& stave, const unsigned int& module, const unsigned int& layer, const unsigned int& slice, const std::array<double, 3>& localPosition) const
            {
                gar::raw::CellID_t cID = 0;

                _decoder->set(cID, "system", det_id);
                _decoder->set(cID, "module", module);
                _decoder->set(cID, "stave", stave);
                _decoder->set(cID, "layer", layer);
                _decoder->set(cID, "slice", slice);

                double localX = localPosition[0];
                double localY = localPosition[1];

                //Segmentation in Y
                int nCellsX = 1;
                int nCellsY = int(_layer_dim_Y / _stripSizeY);

                int _cellIndexX = int ( localX / ( _layer_dim_X / nCellsX ) );
                int _cellIndexY = int ( localY / ( _layer_dim_Y / nCellsY ) );

                _decoder->set(cID, _xId, _cellIndexX);
                _decoder->set(cID, _yId, _cellIndexY);

                return cID;
            }

            //----------------------------------------------------------------------------
            void SegmentationStripYAlg::PrintParameters() const
            {
                std::cout << "cell encoding: " << _encoding << std::endl;
                std::cout << "identifier_x: " << _xId << std::endl;
                std::cout << "identifier_y: " << _yId << std::endl;
                std::cout << "strip_size_y: " << _stripSizeY << " cm" << std::endl;

                std::cout << "identifier_layer: " << _layerId << std::endl;
                std::cout << "identifier_slice: " << _sliceId << std::endl;
            }

            //----------------------------------------------------------------------------
            bool SegmentationStripYAlg::isTile(const gar::raw::CellID_t& cID) const
            {
                return false;
            }

            //----------------------------------------------------------------------------
            bool SegmentationStripYAlg::isBarrel(const gar::raw::CellID_t& cID) const
            {
                return true;
            }

            //----------------------------------------------------------------------------
            double SegmentationStripYAlg::getStripLength(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const gar::raw::CellID_t& cID) const
            {
                /* no op */
                return 0.;
            }

            //----------------------------------------------------------------------------
            std::pair<TVector3, TVector3> SegmentationStripYAlg::getStripEnds(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const gar::raw::CellID_t& cID) const
            {
                /* no op */
                return std::make_pair(TVector3(), TVector3());
            }

            //----------------------------------------------------------------------------
            std::pair<float, float> SegmentationStripYAlg::CalculateLightPropagation(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const gar::raw::CellID_t& cID) const
            {
                /* no op */
                return std::make_pair(0., 0.);
            }

            //----------------------------------------------------------------------------
            std::array<double, 3> SegmentationStripYAlg::ReconstructStripHitPosition(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const float &xlocal, const gar::raw::CellID_t& cID) const
            {
                /* no op */
                return std::array<double, 3>{ 0., 0., 0.};
            }

        }//seg
    } // geo
} //gar
