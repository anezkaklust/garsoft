#include "Geometry/ChannelMapAlgs/ECALSegmentationMultiGridStripXYAlg.h"

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
            ECALSegmentationMultiGridStripXYAlg::ECALSegmentationMultiGridStripXYAlg(fhicl::ParameterSet const& pset)
            : ECALSegmentationAlg(pset)
            {
                _type = "MultiGridStripXY";
                _description = "Cartesian segmentation in the local XY-plane, containing integer number of tiles/strips/cells";

                std::cout << "######### gar::geo::seg::ECALSegmentationMultiGridStripXYAlg() " << std::endl ;

                this->reconfigure(pset);
            }

            //----------------------------------------------------------------------------
            /// Default constructor used by derived classes passing an existing decoder
            ECALSegmentationMultiGridStripXYAlg::ECALSegmentationMultiGridStripXYAlg(const BitFieldCoder* decode, fhicl::ParameterSet const& pset)
            : ECALSegmentationAlg(decode, pset)
            {
                _type = "MultiGridStripXY";
                _description = "Cartesian segmentation in the local XY-plane, containing integer number of tiles/strips/cells";

                std::cout << "######### gar::geo::seg::ECALSegmentationMultiGridStripXYAlg() " << std::endl ;

                this->reconfigure(pset);
            }

            //----------------------------------------------------------------------------
            ECALSegmentationMultiGridStripXYAlg::~ECALSegmentationMultiGridStripXYAlg()
            {
            }

            //----------------------------------------------------------------------------
            void ECALSegmentationMultiGridStripXYAlg::reconfigure(fhicl::ParameterSet const& pset)
            {
                _gridSizeX = pset.get<double>("grid_size_x");
                _gridSizeY = pset.get<double>("grid_size_y");

                _stripSizeX = pset.get<double>("strip_size_x");
                _stripSizeY = pset.get<double>("strip_size_y");
                _encoding = pset.get<std::string>("cellEncoding");

                _offsetX = pset.get<double>("offset_x");
                _offsetY = pset.get<double>("offset_y");

                _xId = pset.get<std::string>("identifier_x");
                _yId = pset.get<std::string>("identifier_y");
                _layerId = pset.get<std::string>("identifier_layer");
                _sliceId = pset.get<std::string>("identifier_slice");

                _gridBarrelLayers = pset.get< std::vector<std::string> >("grid_barrel_layers");
                _stripBarrelLayers = pset.get< std::vector<std::string> >("strip_barrel_layers");

                _gridEndcapLayers = pset.get< std::vector<std::string> >("grid_endcap_layers");
                _stripEndcapLayers = pset.get< std::vector<std::string> >("strip_endcap_layers");

                _nLayers = pset.get<unsigned int>("nlayers");

                _frac = 1./3.;

                _OnSameLayer = pset.get<bool>("strip_on_same_layer");

                this->PrintParameters();

                return;
            }

            //----------------------------------------------------------------------------
            void ECALSegmentationMultiGridStripXYAlg::Initialize(const gar::geo::GeometryCore& geo)
            {

            }

            //----------------------------------------------------------------------------
            std::array<double, 3> ECALSegmentationMultiGridStripXYAlg::GetPosition(const gar::geo::GeometryCore& geo, const gar::raw::CellID_t& cID) const
            {
                //Local origin for the Barrel in the middle of the layer
                //Local origin for the Endcal at the corner of the full stave
                std::array<double, 3> cellPosition;

                //Need to differentiate case tile and strips based on layer and slice
                bool isTile = this->isTile(cID);
                bool isBarrel = this->isBarrel(cID);

                if(isTile) {
                    //Need to check if the tile is at the edge of the layer!
                    cellPosition[0] = binToPosition(_decoder->get(cID, _xId), _gridSizeX, _offsetX);
                    cellPosition[1] = binToPosition(_decoder->get(cID, _yId), _gridSizeY, _offsetY);
                    cellPosition[2] = 0.;

                    if( isBarrel ) {

                        //Check if the position is outside of the layer size
                        if( std::fabs(cellPosition[0]) > _layer_dim_X/2. ) {

                            cellPosition[0] = cellPosition[0] / std::fabs(cellPosition[0]) * ( _layer_dim_X - _frac ) / 2.;
                        }

                        if(std::fabs(cellPosition[1]) > _layer_dim_Y/2. ) {

                            cellPosition[1] = cellPosition[1] / std::fabs(cellPosition[1]) * ( _layer_dim_Y - _frac ) / 2.;
                        }
                    }

                    //TODO
                    //Problem here is that the full layer is the endcap stave, no "layer shape"
                    //the problem will occur for hits that are at the limit
                    //checking if x or y are over the layer dim does not work (as the stave is not a square/rectangle)
                    if( not isBarrel ) {

                        float r_point = std::sqrt( cellPosition[0]*cellPosition[0] + cellPosition[1]*cellPosition[1] );

                        //Check if the point is outside.... layer_dim_x = layer_dim_y = apothem
                        if(r_point > _layer_dim_Y) {

                            if( std::fabs(cellPosition[0]) < geo.GetECALEndcapSideLength() / 2. || std::fabs(cellPosition[1]) < geo.GetECALEndcapSideLength() / 2. ) {

                                if( std::fabs(cellPosition[0]) < geo.GetECALEndcapSideLength() / 2. && std::fabs(cellPosition[1]) > _layer_dim_Y ) {

                                    cellPosition[1] = cellPosition[1] / std::fabs(cellPosition[1]) * ( _layer_dim_Y - _frac ) ;
                                }

                                if( std::fabs(cellPosition[1]) < geo.GetECALEndcapSideLength() / 2. && std::fabs(cellPosition[0]) > _layer_dim_X ) {

                                    cellPosition[0] = cellPosition[0] / std::fabs(cellPosition[0]) * ( _layer_dim_X - _frac ) ;
                                }

                            } else if ( std::fabs(cellPosition[0]) > geo.GetECALEndcapSideLength() / 2. && std::fabs(cellPosition[1]) > geo.GetECALEndcapSideLength() / 2. ) {
                                //Complicated part //do NOTHING the hit will be dropped
                            }
                        }
                    }
                }
                else
                {
                    int cellIndexX = _decoder->get(cID,_xId);
                    int cellIndexY = _decoder->get(cID,_yId);

                    if( (_OnSameLayer && _decoder->get(cID, _sliceId) == 2) || (not _OnSameLayer && _decoder->get(cID, _layerId)%2 == 0) )
                    {
                        //Segmentation in Y
                        int nCellsX = 1;
                        // int nCellsY = int(_layer_dim_Y / _stripSizeY);

                        if(isBarrel) {
                            cellPosition[0] = ( cellIndexX + 0.5 ) * (_layer_dim_X / nCellsX ) - (_layer_dim_X / 2.);
                            cellPosition[1] = ( cellIndexY + 0.5 ) * _stripSizeY;
                            cellPosition[2] = 0.;
                        }
                        else {
                            cellPosition[1] = ( cellIndexY + 0.5 ) * _stripSizeY;
                            cellPosition[0] = ( cellIndexX + 0.5 ) * ( _layer_dim_X / nCellsX );
                            if( cellPosition[1] > geo.GetECALEndcapSideLength() / 2.) {
                                //Correct for the endcap
                                float angle = geo.GetECALInnerAngle() - M_PI/2;
                                int n = std::round( ( cellPosition[1] - ( geo.GetECALEndcapSideLength() / 2. ) ) / _stripSizeX );
                                cellPosition[0] -= n * ( _stripSizeX / std::tan(angle) ) / 2.;
                            }
                            cellPosition[2] = 0.;
                        }
                    }

                    if( (_OnSameLayer && _decoder->get(cID, _sliceId) == 3) || (not _OnSameLayer && _decoder->get(cID, _layerId)%2 != 0) )
                    {
                        //Segmentation in X
                        // int nCellsX = int(_layer_dim_X / _stripSizeX);
                        int nCellsY = 1;

                        if(isBarrel) {
                            cellPosition[0] = ( cellIndexX + 0.5 ) * _stripSizeX;
                            cellPosition[1] = ( cellIndexY + 0.5 ) * (_layer_dim_Y / nCellsY ) - (_layer_dim_Y / 2.);
                            cellPosition[2] = 0.;
                        }
                        else {
                            cellPosition[0] = ( cellIndexX + 0.5 ) * _stripSizeX;
                            cellPosition[1] = ( cellIndexY + 0.5 ) * ( _layer_dim_Y / nCellsY );
                            if( cellPosition[0] > geo.GetECALEndcapSideLength() / 2.) {
                                //Correct for the endcap
                                float angle = geo.GetECALInnerAngle() - M_PI/2;
                                int n = std::round( ( cellPosition[0] - ( geo.GetECALEndcapSideLength() / 2. ) ) / _stripSizeY );
                                cellPosition[1] -= n * ( _stripSizeY / std::tan(angle) ) / 2.;
                            }
                            cellPosition[2] = 0.;
                        }
                    }
                }

                return cellPosition;
            }

            //----------------------------------------------------------------------------
            /// determine the cell ID based on the position
            gar::raw::CellID_t ECALSegmentationMultiGridStripXYAlg::GetCellID(const gar::geo::GeometryCore& geo, const unsigned int& det_id, const unsigned int& stave, const unsigned int& module, const unsigned int& layer, const unsigned int& slice, const std::array<double, 3>& localPosition) const
            {
                gar::raw::CellID_t cID = 0;

                _decoder->set(cID, "system", det_id);
                _decoder->set(cID, "module", module);
                _decoder->set(cID, "stave", stave);
                _decoder->set(cID, "layer", layer);
                _decoder->set(cID, "slice", slice);

                //Need to differentiate case tile and strips based on layer and slice
                bool isTile = this->isTile(cID);

                double localX = localPosition[0];
                double localY = localPosition[1];

                if(isTile)
                {
                    _decoder->set(cID, _xId, positionToBin(localX, _gridSizeX, _offsetX));
                    _decoder->set(cID, _yId, positionToBin(localY, _gridSizeY, _offsetY));
                }
                else
                {
                    if( (_OnSameLayer && slice == 2) || (not _OnSameLayer && layer%2 == 0) )
                    {
                        //Segmentation in Y
                        int nCellsX = 1;
                        int nCellsY = int(_layer_dim_Y / _stripSizeY);

                        int _cellIndexX = int ( localX / ( _layer_dim_X / nCellsX ) );
                        int _cellIndexY = int ( localY / ( _layer_dim_Y / nCellsY ) );

                        _decoder->set(cID, _xId, _cellIndexX);
                        _decoder->set(cID, _yId, _cellIndexY);
                    }

                    if( (_OnSameLayer && slice == 3) || (not _OnSameLayer && layer%2 != 0) )
                    {
                        //Segmentation in X
                        int nCellsX = int(_layer_dim_X / _stripSizeX);
                        int nCellsY = 1;

                        int _cellIndexX = int ( localX / ( _layer_dim_X / nCellsX ) );
                        int _cellIndexY = int ( localY / ( _layer_dim_Y / nCellsY ) );

                        _decoder->set(cID, _xId, _cellIndexX);
                        _decoder->set(cID, _yId, _cellIndexY);
                    }
                }

                return cID;
            }

            //----------------------------------------------------------------------------
            void ECALSegmentationMultiGridStripXYAlg::PrintParameters() const
            {
                std::cout << "cell encoding: " << _encoding << std::endl;
                std::cout << "identifier_x: " << _xId << std::endl;
                std::cout << "identifier_y: " << _yId << std::endl;
                std::cout << "grid_size_x: " << _gridSizeX << " cm" << std::endl;
                std::cout << "grid_size_y: " << _gridSizeY << " cm" << std::endl;
                std::cout << "strip_size_x: " << _stripSizeX << " cm" << std::endl;
                std::cout << "strip_size_y: " << _stripSizeY << " cm" << std::endl;
                std::cout << "offset_x: " << _offsetX << " cm" << std::endl;
                std::cout << "offset_y: " << _offsetY << " cm" << std::endl;

                std::cout << "identifier_layer: " << _layerId << std::endl;
                std::cout << "identifier_slice: " << _sliceId << std::endl;

                std::cout << "grid_barrel_layers: ";
                for(unsigned int i = 0; i < _gridBarrelLayers.size(); i++)
                    std::cout << _gridBarrelLayers.at(i) << " ";
                std::cout << std::endl;

                std::cout << "strip_barrel_layers: ";
                for(unsigned int i = 0; i < _stripBarrelLayers.size(); i++)
                    std::cout << _stripBarrelLayers.at(i) << " ";
                std::cout << std::endl;

                std::cout << "grid_endcap_layers: ";
                for(unsigned int i = 0; i < _gridEndcapLayers.size(); i++)
                    std::cout << _gridEndcapLayers.at(i) << " ";
                std::cout << std::endl;

                std::cout << "strip_endcap_layers: ";
                for(unsigned int i = 0; i < _stripEndcapLayers.size(); i++)
                    std::cout << _stripEndcapLayers.at(i) << " ";
                std::cout << std::endl;

                std::cout << "strip_on_same_layer: " << _OnSameLayer << std::endl;
            }

            //----------------------------------------------------------------------------
            std::array<std::vector<unsigned int>, 2> ECALSegmentationMultiGridStripXYAlg::TokenizeLayerVectors(std::vector<std::string> grid) const
            {
                std::vector<unsigned int> _gridFirst;
                std::vector<unsigned int> _gridLast;

                for(unsigned int i = 0; i < grid.size(); i++)
                {
                    std::vector<std::string> descriptors ;
                    Tokenizer t( descriptors ,':') ;

                    std::string value = grid.at(i);
                    std::for_each( value.begin(), value.end(), t ) ;

                    _gridFirst.push_back(std::atoi(descriptors.at(0).c_str()));
                    _gridLast.push_back(std::atoi(descriptors.at(1).c_str()));
                }

                std::array<std::vector<unsigned int>, 2> _list = {{_gridFirst, _gridLast}};
                return _list;
            }

            //----------------------------------------------------------------------------
            bool ECALSegmentationMultiGridStripXYAlg::isTile(const gar::raw::CellID_t& cID) const
            {
                bool isTile = false;

                unsigned int det_id = _decoder->get(cID, "system");
                unsigned int layer = _decoder->get(cID, _layerId);

                std::array<std::vector<unsigned int>, 2> _list;

                if(det_id == 1)
                    _list = TokenizeLayerVectors(_gridBarrelLayers);
                if(det_id == 2)
                    _list =  TokenizeLayerVectors(_gridEndcapLayers);

                //Check if it is tile configuration
                for(unsigned int i = 0; i < _list.at(0).size(); i++)
                {
                    if( layer >= _list.at(0).at(i) && layer <= _list.at(1).at(i) )
                    {
                        isTile = true;
                        break;
                    }
                    else{
                        continue;
                    }
                }

                return isTile;
            }

            //----------------------------------------------------------------------------
            bool ECALSegmentationMultiGridStripXYAlg::isBarrel(const gar::raw::CellID_t& cID) const
            {
                bool isBarrel = true;

                int det_id = _decoder->get(cID, "system");
                int module = _decoder->get(cID, "module");
                if( det_id == 2 && (module == 0 || module == 6) ) isBarrel = false;

                return isBarrel;
            }

            //----------------------------------------------------------------------------
            double ECALSegmentationMultiGridStripXYAlg::getStripLength(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const gar::raw::CellID_t& cID) const
            {
                double stripLength = 0.;
                bool isBarrel = this->isBarrel(cID);

                //Strip along X
                if( (_OnSameLayer && _decoder->get(cID, _sliceId) == 2) || (not _OnSameLayer && _decoder->get(cID, _layerId)%2 == 0) ) {
                    if(isBarrel) {
                        stripLength = _layer_dim_X;
                    } else {
                        float angle = geo.GetECALInnerAngle() - M_PI/2;
                        int n = std::round( (local[1] - ( geo.GetECALEndcapSideLength() / 2. )) / _stripSizeX );
                        if ( n < 0 ) n = 0;
                        stripLength = _layer_dim_X - ( n * ( _stripSizeX / std::tan(angle) ) );
                    }
                }
                if( (_OnSameLayer && _decoder->get(cID, _sliceId) == 3) || (not _OnSameLayer && _decoder->get(cID, _layerId)%2 != 0) ) {
                    //Strip along Y
                    if(isBarrel) {
                        stripLength = _layer_dim_Y;
                    } else {
                        float angle = geo.GetECALInnerAngle() - M_PI/2;
                        int n = std::round( (local[1] - ( geo.GetECALEndcapSideLength() / 2. )) / _stripSizeY );
                        if ( n < 0 ) n = 0;
                        stripLength = _layer_dim_Y - ( n * ( _stripSizeY / std::tan(angle) ) );
                    }
                }

                return stripLength;
            }

            //----------------------------------------------------------------------------
            std::pair<TVector3, TVector3> ECALSegmentationMultiGridStripXYAlg::getStripEnds(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const gar::raw::CellID_t& cID) const
            {
                //Local origin for the Barrel in the middle of the layer
                //Local origin for the Endcal at the corner of the full stave
                TVector3 stripEnd1(0., 0., 0.);
                TVector3 stripEnd2(0., 0., 0.);
                bool isBarrel = this->isBarrel(cID);

                //Strip along X (local)
                if( (_OnSameLayer && _decoder->get(cID, _sliceId) == 2) || (not _OnSameLayer && _decoder->get(cID, _layerId)%2 == 0) ) {

                    if(isBarrel) {
                        stripEnd1.SetX(-_layer_dim_X/2.);
                        stripEnd2.SetX(_layer_dim_X/2.);
                    } else {
                        stripEnd1.SetX(0.);
                        //not perfect at the edges...
                        float angle = geo.GetECALInnerAngle() - M_PI/2;
                        int n = std::round( (local[1] - ( geo.GetECALEndcapSideLength() / 2. )) / _stripSizeX );
                        if ( n < 0 ) n = 0;
                        stripEnd2.SetX( _layer_dim_X - ( n * ( _stripSizeX / std::tan(angle) ) ) );
                    }

                    stripEnd1.SetY(local[1]);
                    stripEnd2.SetY(local[1]);
                    stripEnd1.SetZ(local[2]);
                    stripEnd2.SetZ(local[2]);
                }
                //Strip along Y (local)
                if( (_OnSameLayer && _decoder->get(cID, _sliceId) == 3) || (not _OnSameLayer && _decoder->get(cID, _layerId)%2 != 0) ) {

                    if(isBarrel) {
                        stripEnd1.SetY(-_layer_dim_Y/2.);
                        stripEnd2.SetY(_layer_dim_Y/2.);
                    } else {
                        stripEnd1.SetY(0.);
                        float angle = geo.GetECALInnerAngle() - M_PI/2;
                        int n = std::round( (local[0] - (geo.GetECALEndcapSideLength() / 2.) ) / _stripSizeY );
                        if ( n < 0 ) n = 0;
                        stripEnd2.SetY( _layer_dim_Y - ( n * ( _stripSizeY / std::tan(angle) ) ) );//not perfect at the edges...
                    }

                    stripEnd1.SetX(local[0]);
                    stripEnd2.SetX(local[0]);
                    stripEnd1.SetZ(local[2]);
                    stripEnd2.SetZ(local[2]);
                }

                return std::make_pair(stripEnd1, stripEnd2);
            }

            //----------------------------------------------------------------------------
            std::pair<float, float> ECALSegmentationMultiGridStripXYAlg::CalculateLightPropagation(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const gar::raw::CellID_t& cID) const
            {
                //Propagate the light to the SiPM on the side
                //convert c to mm/ns
                // c   = 29.9792458 cm/ns
                float c = (CLHEP::c_light * CLHEP::mm / CLHEP::ns) / CLHEP::cm; // in cm/ns
                //time1 is left SiPM
                float time1 = 0.;
                //time2 is right SiPM
                float time2 = 0.;
                //check if it is barrel or not
                bool isBarrel = this->isBarrel(cID);

                //Strip along X (local)
                if( (_OnSameLayer && _decoder->get(cID, _sliceId) == 2) || (not _OnSameLayer && _decoder->get(cID, _layerId)%2 == 0) )
                {
                    if(isBarrel) {
                        time1 = ( _layer_dim_X / 2 + local[0] ) / c;
                        time2 = ( _layer_dim_X / 2 - local[0] ) / c;
                    } else {
                        float angle = geo.GetECALInnerAngle() - M_PI/2;
                        int n = std::round( (local[0] - ( geo.GetECALEndcapSideLength() / 2. )) / _stripSizeY );
                        if ( n < 0 ) n = 0;
                        float strip_dim = _layer_dim_X - ( n * ( _stripSizeY / std::tan(angle) ) );
                        time1 = ( strip_dim / 2. + local[0] ) / c;
                        time2 = ( strip_dim / 2. - local[0] ) / c;
                    }
                }

                //Strip along Y (local)
                if( (_OnSameLayer && _decoder->get(cID, _sliceId) == 3) || (not _OnSameLayer && _decoder->get(cID, _layerId)%2 != 0) )
                {
                    if(isBarrel) {
                        time1 = ( _layer_dim_Y / 2 + local[1] ) / c;
                        time2 = ( _layer_dim_Y / 2 - local[1] ) / c;
                    } else {
                        float angle = geo.GetECALInnerAngle() - M_PI/2;
                        int n = std::round( (local[1] - ( geo.GetECALEndcapSideLength() / 2. )) / _stripSizeX );
                        if ( n < 0 ) n = 0;
                        float strip_dim = _layer_dim_Y - ( n * ( _stripSizeX / std::tan(angle) ) );
                        time1 = ( strip_dim / 2. + local[1] ) / c;
                        time2 = ( strip_dim / 2. - local[1] ) / c;
                    }
                }

                return std::make_pair(time1, time2);
            }

            //----------------------------------------------------------------------------
            std::array<double, 3> ECALSegmentationMultiGridStripXYAlg::ReconstructStripHitPosition(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const float &xlocal, const gar::raw::CellID_t& cID) const
            {
                float pos = xlocal;
                //Need to check if the local position is bigger than the strip length!
                bool isBarrel = this->isBarrel(cID);
                float stripLength = this->getStripLength(geo, local, cID);

                if(isBarrel) {
                    if( std::abs(pos) > stripLength / 2. ) pos = (pos > 0) ? stripLength / 2. : -stripLength / 2.;
                } else {
                    if( std::abs(pos) > stripLength ) pos = (pos > 0) ? stripLength : 0.;
                }

                std::array<double, 3> newlocal;

                if( (_OnSameLayer && _decoder->get(cID, _sliceId) == 2) || (not _OnSameLayer && _decoder->get(cID, _layerId)%2 == 0) )
                    newlocal = {pos, local[1], local[2]};

                if( (_OnSameLayer && _decoder->get(cID, _sliceId) == 3) || (not _OnSameLayer && _decoder->get(cID, _layerId)%2 != 0) ) //Strip along Y
                    newlocal = {local[0], pos, local[2]};

                return newlocal;
            }

        }//seg
    } // geo
} //gar
