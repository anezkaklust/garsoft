#include "Geometry/ChannelMapAlgs/SegmentationMuIDAlg.h"

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
            SegmentationMuIDAlg::SegmentationMuIDAlg(fhicl::ParameterSet const& pset)
            : SegmentationAlg(pset)
            {
                _type = "MuID";
                _description = "Muon system segmentation with cross-strips depending on the slice";

                std::cout << "######### gar::geo::seg::SegmentationMuIDAlg() " << std::endl ;

                this->reconfigure(pset);
            }

            //----------------------------------------------------------------------------
            /// Default constructor used by derived classes passing an existing decoder
            SegmentationMuIDAlg::SegmentationMuIDAlg(const BitFieldCoder* decode, fhicl::ParameterSet const& pset)
            : SegmentationAlg(decode, pset)
            {
                _type = "MuID";
                _description = "Muon system segmentation with cross-strips depending on the slice";

                std::cout << "######### gar::geo::seg::SegmentationMuIDAlg() " << std::endl ;

                this->reconfigure(pset);
            }

            //----------------------------------------------------------------------------
            SegmentationMuIDAlg::~SegmentationMuIDAlg()
            {
            }

            //----------------------------------------------------------------------------
            void SegmentationMuIDAlg::reconfigure(fhicl::ParameterSet const& pset)
            {
                _stripSizeY = pset.get<double>("strip_size_y");
                _stripSizeX = _stripSizeY;
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
            void SegmentationMuIDAlg::Initialize(const gar::geo::GeometryCore& geo)
            {

            }

            //----------------------------------------------------------------------------
            std::array<double, 3> SegmentationMuIDAlg::GetPosition(const gar::geo::GeometryCore& geo, const gar::raw::CellID_t& cID) const
            {
                std::array<double, 3> cellPosition;

                int cellIndexX = _decoder->get(cID, _xId);
                int cellIndexY = _decoder->get(cID, _yId);
                int layer = _decoder->get(cID, "layer");

                if (layer % 2 == 0)
                {
                    //Segmentation in Y
                    int nCellsX = 1;
                    // int nCellsY = int(_layer_dim_Y / _stripSizeY);

                    cellPosition[0] = ( cellIndexX + 0.5 ) * (_layer_dim_X / nCellsX ) - (_layer_dim_X / 2.);
                    cellPosition[1] = ( cellIndexY + 0.5 ) * _stripSizeY;
                    cellPosition[2] = 0.;
                }

                if(layer % 2 != 0) {
                    //Segmentation in X
                    int nCellsY = 1;
                    // int nCellsX = int(_layer_dim_Y / _stripSizeY);

                    cellPosition[0] = ( cellIndexX + 0.5 ) * _stripSizeX;
                    cellPosition[1] = ( cellIndexY + 0.5 ) * (_layer_dim_Y / nCellsY ) - (_layer_dim_Y / 2.);
                    cellPosition[2] = 0.;
                }

                return cellPosition;
            }

            //----------------------------------------------------------------------------
            /// determine the cell ID based on the position
            gar::raw::CellID_t SegmentationMuIDAlg::GetCellID(const gar::geo::GeometryCore& geo, const unsigned int& det_id, const unsigned int& stave, const unsigned int& module, const unsigned int& layer, const unsigned int& slice, const std::array<double, 3>& localPosition) const
            {
                gar::raw::CellID_t cID = 0;

                _decoder->set(cID, "system", det_id);
                _decoder->set(cID, "module", module);
                _decoder->set(cID, "stave", stave);
                _decoder->set(cID, "layer", layer);
                _decoder->set(cID, "slice", slice);

                double localX = localPosition[0];
                double localY = localPosition[1];

                //Segmentation depending on the slice in the layer (1 layer with Sc/Fe/Sc/Fe/Sc)
                //Slices 1/3/5 active
                int nCellsX = 1;
                int nCellsY = 1;

                if(layer % 2 == 0) {
                    //Segmentation in Y
                    nCellsX = 1;
                    nCellsY = int(_layer_dim_Y / _stripSizeY);

                }

                if(layer % 2 != 0) {
                    //Segmentation in X
                    nCellsY = 1;
                    nCellsX = int(_layer_dim_X / _stripSizeX);
                }

                int _cellIndexX = int ( localX / ( _layer_dim_X / nCellsX ) );
                int _cellIndexY = int ( localY / ( _layer_dim_Y / nCellsY ) );

                _decoder->set(cID, _xId, _cellIndexX);
                _decoder->set(cID, _yId, _cellIndexY);

                return cID;
            }

            //----------------------------------------------------------------------------
            void SegmentationMuIDAlg::PrintParameters() const
            {
                std::cout << "cell encoding: " << _encoding << std::endl;
                std::cout << "identifier_x: " << _xId << std::endl;
                std::cout << "identifier_y: " << _yId << std::endl;
                std::cout << "strip_size_x: " << _stripSizeX << " cm" << std::endl;
                std::cout << "strip_size_y: " << _stripSizeY << " cm" << std::endl;

                std::cout << "identifier_layer: " << _layerId << std::endl;
                std::cout << "identifier_slice: " << _sliceId << std::endl;
                std::cout << "nLayers: " << _nLayers << std::endl;
            }

            //----------------------------------------------------------------------------
            bool SegmentationMuIDAlg::isTile(const gar::raw::CellID_t& cID) const
            {
                return false;
            }

            //----------------------------------------------------------------------------
            bool SegmentationMuIDAlg::isBarrel(const gar::raw::CellID_t& cID) const
            {
                return true;
            }

            //----------------------------------------------------------------------------
            double SegmentationMuIDAlg::getStripLength(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const gar::raw::CellID_t& cID) const
            {
                double stripLength = 0.;
                int layer = _decoder->get(cID, "layer");

                if(layer == 1 || layer == 3) stripLength = _layer_dim_X;
                if(layer == 2) stripLength = _layer_dim_Y;

                return stripLength;
            }

            //----------------------------------------------------------------------------
            std::pair<TVector3, TVector3> SegmentationMuIDAlg::getStripEnds(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const gar::raw::CellID_t& cID) const
            {
                int layer = _decoder->get(cID, "layer");

                TVector3 stripEnd1(0., 0., 0.);
                TVector3 stripEnd2(0., 0., 0.);

                if(layer == 1 || layer == 3) {
                    stripEnd1.SetX(-_layer_dim_X/2.);
                    stripEnd2.SetX(_layer_dim_X/2.);
                    stripEnd1.SetY(local[1]);
                    stripEnd2.SetY(local[1]);
                    stripEnd1.SetZ(local[2]);
                    stripEnd2.SetZ(local[2]);
                }

                if(layer == 2) {
                    stripEnd1.SetY(-_layer_dim_Y/2.);
                    stripEnd2.SetY(_layer_dim_Y/2.);
                    stripEnd1.SetX(local[0]);
                    stripEnd2.SetX(local[0]);
                    stripEnd1.SetZ(local[2]);
                    stripEnd2.SetZ(local[2]);
                }


                return std::make_pair(stripEnd1, stripEnd2);
            }

            //----------------------------------------------------------------------------
            std::pair<float, float> SegmentationMuIDAlg::CalculateLightPropagation(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const gar::raw::CellID_t& cID) const
            {
                float c = (CLHEP::c_light * CLHEP::mm / CLHEP::ns) / CLHEP::cm; // in cm/ns
                //time1 is left SiPM
                float time1 = 0.;
                //time2 is right SiPM
                float time2 = 0.;

                int layer = _decoder->get(cID, "layer");

                if(layer == 1 || layer == 3) {
                    time1 = ( _layer_dim_X / 2 + local[0] ) / c;
                    time2 = ( _layer_dim_X / 2 - local[0] ) / c;
                }

                if(layer == 2) {
                    time1 = ( _layer_dim_Y / 2 + local[1] ) / c;
                    time2 = ( _layer_dim_Y / 2 - local[1] ) / c;
                }

                return std::make_pair(time1, time2);
            }

            //----------------------------------------------------------------------------
            std::array<double, 3> SegmentationMuIDAlg::ReconstructStripHitPosition(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const float &xlocal, const gar::raw::CellID_t& cID) const
            {
                int layer = _decoder->get(cID, "layer");
                float pos = xlocal;
                std::array<double, 3> newlocal;

                //Need to check if the local position is bigger than the strip length!
                // bool isBarrel = this->isBarrel(cID);
                float stripLength = this->getStripLength(geo, local, cID);

                if( std::abs(pos) > stripLength / 2. ) pos = (pos > 0) ? stripLength / 2. : -stripLength / 2.;

                if(layer == 1 || layer == 3) {
                    newlocal = {pos, local[1], local[2]};
                }

                if(layer == 2) {
                    newlocal = {local[0], pos, local[2]};
                }

                return newlocal;
            }

        }//seg
    } // geo
} //gar
