#include "Geometry/ChannelMapAlgs/MinervaSegmentationAlg.h"

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
            MinervaSegmentationAlg::MinervaSegmentationAlg(fhicl::ParameterSet const& pset)
            : SegmentationAlg(pset)
            {
                _type = "StripXY";
                _description = "Cartesian segmentation in the local XY-plane, containing integer number of strips";

                std::cout << "######### gar::geo::seg::MinervaSegmentationAlg() " << std::endl ;

                this->reconfigure(pset);
            }

            //----------------------------------------------------------------------------
            /// Default constructor used by derived classes passing an existing decoder
            MinervaSegmentationAlg::MinervaSegmentationAlg(const BitFieldCoder* decode, fhicl::ParameterSet const& pset)
            : SegmentationAlg(decode, pset)
            {
                _type = "StripXY";
                _description = "Cartesian segmentation in the local XY-plane, containing integer number of strips";

                std::cout << "######### gar::geo::seg::MinervaSegmentationAlg() " << std::endl ;

                this->reconfigure(pset);
            }

            //----------------------------------------------------------------------------
            MinervaSegmentationAlg::~MinervaSegmentationAlg()
            {
            }

            //----------------------------------------------------------------------------
            void MinervaSegmentationAlg::reconfigure(fhicl::ParameterSet const& pset)
            {
                _stripSizeX = pset.get<double>("strip_size_x");
                _stripSizeY = pset.get<double>("strip_size_y");
                _encoding = pset.get<std::string>("cellEncoding");

                _xId = pset.get<std::string>("identifier_x");
                _yId = pset.get<std::string>("identifier_y");
                _zId = pset.get<std::string>("identifier_z");
                _layerId = pset.get<std::string>("identifier_layer");
                _sliceId = pset.get<std::string>("identifier_slice");

                _frac = 1./3.;

                this->PrintParameters();

                return;
            }

            //----------------------------------------------------------------------------
            void MinervaSegmentationAlg::Initialize(const gar::geo::GeometryCore& geo)
            {
                
            }

            //----------------------------------------------------------------------------
            std::array<double, 3> MinervaSegmentationAlg::GetPosition(const gar::geo::GeometryCore& geo, const gar::raw::CellID_t& cID) const
            {
                //Local origin for the Barrel in the middle of the layer
                //Local origin for the Endcal at the corner of the full stave
                std::array<double, 3> cellPosition;

                /* NO OP */

                return cellPosition;
            }

            //----------------------------------------------------------------------------
            /// determine the cell ID based on the position
            gar::raw::CellID_t MinervaSegmentationAlg::GetCellID(const gar::geo::GeometryCore& geo, const unsigned int& det_id, const unsigned int& stave, const unsigned int& module, const unsigned int& layer, const unsigned int& slice, const std::array<double, 3>& localPosition) const
            {
                gar::raw::CellID_t cID = 0;

                _decoder->set(cID, "system", det_id);
                _decoder->set(cID, "layer", layer);
                _decoder->set(cID, "slice", slice);

                double localX = localPosition[0];
                double localY = localPosition[1];
                double localZ = localPosition[2];

                if( localZ < 0 )
                {
                    //Segmentation in Y
                    int nCellsX = 1;
                    int nCellsY = int(_layer_dim_Y / _stripSizeY);

                    int _cellIndexX = int ( localX / ( _layer_dim_X / nCellsX ) );
                    int _cellIndexY = int ( localY / ( _layer_dim_Y / nCellsY ) );

                    _decoder->set(cID, _xId, _cellIndexX);
                    _decoder->set(cID, _yId, _cellIndexY);
                    _decoder->set(cID, _zId, 0);
                }

                if( localZ >= 0 )
                {
                    //Segmentation in X
                    int nCellsX = int(_layer_dim_X / _stripSizeX);
                    int nCellsY = 1;

                    int _cellIndexX = int ( localX / ( _layer_dim_X / nCellsX ) );
                    int _cellIndexY = int ( localY / ( _layer_dim_Y / nCellsY ) );

                    _decoder->set(cID, _xId, _cellIndexX);
                    _decoder->set(cID, _yId, _cellIndexY);
                    _decoder->set(cID, _zId, 1);
                }

                return cID;
            }

            //----------------------------------------------------------------------------
            void MinervaSegmentationAlg::PrintParameters() const
            {
                std::cout << "cell encoding: " << _encoding << std::endl;
                std::cout << "identifier_x: " << _xId << std::endl;
                std::cout << "identifier_y: " << _yId << std::endl;
                std::cout << "identifier_z: " << _zId << std::endl;
                std::cout << "strip_size_x: " << _stripSizeX << " cm" << std::endl;
                std::cout << "strip_size_y: " << _stripSizeY << " cm" << std::endl;

                std::cout << "identifier_layer: " << _layerId << std::endl;
                std::cout << "identifier_slice: " << _sliceId << std::endl;
            }

        }//seg
    } // geo
} //gar
