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
                    //Need to check in which half of the square is the hit --> {y, z} point is below or over the diagonal
                    //(y = z (odd) or y = -z (even))

                    int nCellsX = 1;
                    int nCellsY = int(_layer_dim_Y / (_stripSizeY * 2));

                    int _cellIndexX = int ( localX / ( _layer_dim_X / nCellsX ) );
                    int _cellIndexY = int ( localY / ( _layer_dim_Y / nCellsY ) );

                    //Transform the localX/Y to the local of this cell
                    // double cellOriginX = 0.;
                    double cellOriginY = ( _cellIndexY + 0.5 ) * (_stripSizeY * 2);
                    double cellOriginZ = - 1.;

                    double localy = localY - cellOriginY; //transform it
                    double localz = localZ - cellOriginZ; //transform it
                    bool above = false;
                    bool below = false;
                    if(localy < 0) {
                        //Need to check if the point is below or above y = z
                        if(localy > localz) {
                            above = true;
                        }
                        if(localy < localz) {
                            below = true;
                        }
                        if(above) {
                            _decoder->set(cID, _xId, _cellIndexX);
                            _decoder->set(cID, _yId, _cellIndexY);
                            _decoder->set(cID, _zId, 0);
                            _decoder->set(cID, "triangle", 1);
                        }
                        if(below) {
                            _decoder->set(cID, _xId, _cellIndexX);
                            _decoder->set(cID, _yId, _cellIndexY);
                            _decoder->set(cID, _zId, 0);
                            _decoder->set(cID, "triangle", 0);
                        }
                    }
                    else if( localy > 0 ) {
                        //Need to check if the point is below or above y = -z
                        if(localy > -localz) {
                            above = true;
                        }
                        if(localy < -localz) {
                            below = true;
                        }
                        if(above) {
                            _decoder->set(cID, _xId, _cellIndexX);
                            _decoder->set(cID, _yId, _cellIndexY);
                            _decoder->set(cID, _zId, 0);
                            _decoder->set(cID, "triangle", 3);
                        }
                        if(below) {
                            _decoder->set(cID, _xId, _cellIndexX);
                            _decoder->set(cID, _yId, _cellIndexY);
                            _decoder->set(cID, _zId, 0);
                            _decoder->set(cID, "triangle", 2);
                        }
                    }
                    else {
                        //exception
                        MF_LOG_ERROR("MinervaSegmentationAlg") << "Cannot determine if localX > | < 0";
                        return 0;
                    }
                } //end if z < 0

                if( localZ >= 0 )
                {
                    //Segmentation in X
                    int nCellsX = int(_layer_dim_X / _stripSizeX);
                    int nCellsY = 1;

                    int _cellIndexX = int ( localX / ( _layer_dim_X / nCellsX ) );
                    int _cellIndexY = int ( localY / ( _layer_dim_Y / nCellsY ) );

                    //Transform the localX/Y to the local of this cell
                    double cellOriginX = ( _cellIndexX + 0.5 ) * (_stripSizeX * 2);
                    // double cellOriginY = 0.;
                    double cellOriginZ = 1.;

                    double localx = localX - cellOriginX; //transform it
                    double localz = localZ - cellOriginZ; //transform it
                    bool above = false;
                    bool below = false;
                    if(localx < 0) {
                        //Need to check if the point is below or above y = z
                        if(localx > localz) {
                            above = true;
                        }
                        if(localx < localz) {
                            below = true;
                        }
                        if(above) {
                            _decoder->set(cID, _xId, _cellIndexX);
                            _decoder->set(cID, _yId, _cellIndexY);
                            _decoder->set(cID, _zId, 1);
                            _decoder->set(cID, "triangle", 1);
                        }
                        if(below) {
                            _decoder->set(cID, _xId, _cellIndexX);
                            _decoder->set(cID, _yId, _cellIndexY);
                            _decoder->set(cID, _zId, 1);
                            _decoder->set(cID, "triangle", 0);
                        }
                    }
                    else if( localx > 0 ) {
                        //Need to check if the point is below or above y = -z
                        if(localx > -localz) {
                            above = true;
                        }
                        if(localx < -localz) {
                            below = true;
                        }
                        if(above) {
                            _decoder->set(cID, _xId, _cellIndexX);
                            _decoder->set(cID, _yId, _cellIndexY);
                            _decoder->set(cID, _zId, 1);
                            _decoder->set(cID, "triangle", 3);
                        }
                        if(below) {
                            _decoder->set(cID, _xId, _cellIndexX);
                            _decoder->set(cID, _yId, _cellIndexY);
                            _decoder->set(cID, _zId, 1);
                            _decoder->set(cID, "triangle", 2);
                        }
                    }
                    else {
                        //exception
                        MF_LOG_ERROR("MinervaSegmentationAlg") << "Cannot determine if localX > | < 0";
                        return 0;
                    }
                } //end if z >= 0

                return cID;
            }

            //------------------------------------------------------------------------------
            gar::raw::CellID_t MinervaSegmentationAlg::GetComplementaryCellID(gar::raw::CellID_t cellID, unsigned int comp) const
            {
                gar::raw::CellID_t cID = 0;

                if(comp == 0) {
                    //Need to get the lower cellX/Y
                    if(_decoder->get(cellID, "cellX") == 0){
                        _decoder->set(cID, "system", _decoder->get(cellID, "system"));
                        _decoder->set(cID, "layer", _decoder->get(cellID, "layer"));
                        _decoder->set(cID, "slice", _decoder->get(cellID, "slice"));
                        _decoder->set(cID, "cellX", _decoder->get(cellID, "cellX"));
                        _decoder->set(cID, "cellY", _decoder->get(cellID, "cellY")+1);
                        _decoder->set(cID, "cellZ", _decoder->get(cellID, "cellZ"));
                        _decoder->set(cID, "triangle", comp);
                    }
                    if(_decoder->get(cellID, "cellY") == 0){
                        _decoder->set(cID, "system", _decoder->get(cellID, "system"));
                        _decoder->set(cID, "layer", _decoder->get(cellID, "layer"));
                        _decoder->set(cID, "slice", _decoder->get(cellID, "slice"));
                        _decoder->set(cID, "cellX", _decoder->get(cellID, "cellX")+1);
                        _decoder->set(cID, "cellY", _decoder->get(cellID, "cellY"));
                        _decoder->set(cID, "cellZ", _decoder->get(cellID, "cellZ"));
                        _decoder->set(cID, "triangle", comp);
                    }
                }

                if(comp == 1 || comp == 2) {
                    //Keep the same cellX/cellY
                    _decoder->set(cID, "system", _decoder->get(cellID, "system"));
                    _decoder->set(cID, "layer", _decoder->get(cellID, "layer"));
                    _decoder->set(cID, "slice", _decoder->get(cellID, "slice"));
                    _decoder->set(cID, "cellX", _decoder->get(cellID, "cellX"));
                    _decoder->set(cID, "cellY", _decoder->get(cellID, "cellY"));
                    _decoder->set(cID, "cellZ", _decoder->get(cellID, "cellZ"));
                    _decoder->set(cID, "triangle", comp);
                }

                if(comp == 3) {
                    //Need to get the upper cellX/Y
                    if(_decoder->get(cellID, "cellX") == 0){
                        _decoder->set(cID, "system", _decoder->get(cellID, "system"));
                        _decoder->set(cID, "layer", _decoder->get(cellID, "layer"));
                        _decoder->set(cID, "slice", _decoder->get(cellID, "slice"));
                        _decoder->set(cID, "cellX", _decoder->get(cellID, "cellX"));
                        _decoder->set(cID, "cellY", _decoder->get(cellID, "cellY")-1);
                        _decoder->set(cID, "cellZ", _decoder->get(cellID, "cellZ"));
                        _decoder->set(cID, "triangle", comp);
                    }
                    if(_decoder->get(cellID, "cellY") == 0){
                        _decoder->set(cID, "system", _decoder->get(cellID, "system"));
                        _decoder->set(cID, "layer", _decoder->get(cellID, "layer"));
                        _decoder->set(cID, "slice", _decoder->get(cellID, "slice"));
                        _decoder->set(cID, "cellX", _decoder->get(cellID, "cellX")-1);
                        _decoder->set(cID, "cellY", _decoder->get(cellID, "cellY"));
                        _decoder->set(cID, "cellZ", _decoder->get(cellID, "cellZ"));
                        _decoder->set(cID, "triangle", comp);
                    }
                }

                return cID;
            }

            //------------------------------------------------------------------------------
            void MinervaSegmentationAlg::AddHitsMinerva(std::map< gar::raw::CellID_t, std::vector<gar::sdp::CaloDeposit> > &m_Deposits, std::vector<gar::sdp::CaloDeposit> &fDeposits) const
            {
                //Loop over the hits in the map and add them together
                for(auto &it : m_Deposits)
                {
                    gar::raw::CellID_t cellID = it.first;

                    //Check this cellID
                    //if it is a triangle = 0 or 3
                    //if it is triangle = 1 or 2
                    //need to look for the complementary cellID, add them to this hit and remove from the deposits (to avoid double hit creation)

                    int triangleNb = _decoder->get(cellID, "triangle");

                    //Case bottom or upper triangle, look for complementary cellID
                    if(triangleNb == 0 || triangleNb == 3)
                    {
                        if(triangleNb == 0) {
                            std::vector<gar::sdp::CaloDeposit> vechit = it.second;
                            std::sort(vechit.begin(), vechit.end()); //sort per time

                            float esum = 0.;
                            float stepLength = 0.;
                            float time = vechit.at(0).Time();
                            int trackID = vechit.at(0).TrackID();
                            double pos[3] = { vechit.at(0).X(), vechit.at(0).Y(), vechit.at(0).Z() };

                            for(auto const &hit : vechit) {
                                esum += hit.Energy();
                                stepLength += hit.StepLength();
                            }

                            //need to check if cellX or cellY are - 1
                            gar::raw::CellID_t complementary_cellID = this->GetComplementaryCellID(cellID, 3);
                            auto find = m_Deposits.find( complementary_cellID );
                            if(find != m_Deposits.end()) {
                                std::vector<gar::sdp::CaloDeposit> vechit_comp = find->second;
                                std::sort(vechit_comp.begin(), vechit_comp.end()); //sort per time

                                for(auto const &hit : vechit_comp) {
                                    esum += hit.Energy();
                                    stepLength += hit.StepLength();
                                }

                                //remove the element from the map now to avoid double counting
                                m_Deposits.erase(find->first);
                            }

                            fDeposits.emplace_back( trackID, time, esum, pos, cellID, stepLength );
                        }

                        if(triangleNb == 3) {

                            std::vector<gar::sdp::CaloDeposit> vechit = it.second;
                            std::sort(vechit.begin(), vechit.end()); //sort per time

                            float esum = 0.;
                            float stepLength = 0.;
                            float time = vechit.at(0).Time();
                            int trackID = vechit.at(0).TrackID();
                            double pos[3] = { vechit.at(0).X(), vechit.at(0).Y(), vechit.at(0).Z() };

                            for(auto const &hit : vechit) {
                                esum += hit.Energy();
                                stepLength += hit.StepLength();
                            }

                            //need to check if cellX or cellY are + 1
                            gar::raw::CellID_t complementary_cellID = this->GetComplementaryCellID(cellID, 0);
                            auto find = m_Deposits.find( complementary_cellID );
                            if(find != m_Deposits.end()) {
                                std::vector<gar::sdp::CaloDeposit> vechit_comp = find->second;
                                std::sort(vechit_comp.begin(), vechit_comp.end()); //sort per time

                                for(auto const &hit : vechit_comp) {
                                    esum += hit.Energy();
                                    stepLength += hit.StepLength();
                                }

                                //remove the element from the map now to avoid double counting
                                m_Deposits.erase(find->first);
                            }

                            fDeposits.emplace_back( trackID, time, esum, pos, cellID, stepLength );
                        }
                    }

                    if(triangleNb == 1 || triangleNb == 2)
                    {
                        if(triangleNb == 1) {
                            std::vector<gar::sdp::CaloDeposit> vechit = it.second;
                            std::sort(vechit.begin(), vechit.end()); //sort per time

                            float esum = 0.;
                            float stepLength = 0.;
                            float time = vechit.at(0).Time();
                            int trackID = vechit.at(0).TrackID();
                            double pos[3] = { vechit.at(0).X(), vechit.at(0).Y(), vechit.at(0).Z() };

                            for(auto const &hit : vechit) {
                                esum += hit.Energy();
                                stepLength += hit.StepLength();
                            }

                            //need to check if cellX or cellY are the same
                            gar::raw::CellID_t complementary_cellID = this->GetComplementaryCellID(cellID, 2);
                            auto find = m_Deposits.find( complementary_cellID );
                            if(find != m_Deposits.end()) {
                                std::vector<gar::sdp::CaloDeposit> vechit_comp = find->second;
                                std::sort(vechit_comp.begin(), vechit_comp.end()); //sort per time

                                for(auto const &hit : vechit_comp) {
                                    esum += hit.Energy();
                                    stepLength += hit.StepLength();
                                }

                                //remove the element from the map now to avoid double counting
                                m_Deposits.erase(find->first);
                            }
                            fDeposits.emplace_back( trackID, time, esum, pos, cellID, stepLength );
                        }

                        if(triangleNb == 2) {
                            std::vector<gar::sdp::CaloDeposit> vechit = it.second;
                            std::sort(vechit.begin(), vechit.end()); //sort per time

                            float esum = 0.;
                            float stepLength = 0.;
                            float time = vechit.at(0).Time();
                            int trackID = vechit.at(0).TrackID();
                            double pos[3] = { vechit.at(0).X(), vechit.at(0).Y(), vechit.at(0).Z() };

                            for(auto const &hit : vechit) {
                                esum += hit.Energy();
                                stepLength += hit.StepLength();
                            }

                            //need to check if cellX or cellY are the same
                            gar::raw::CellID_t complementary_cellID = this->GetComplementaryCellID(cellID, 1);
                            auto find = m_Deposits.find( complementary_cellID );
                            if(find != m_Deposits.end()) {
                                std::vector<gar::sdp::CaloDeposit> vechit_comp = find->second;
                                std::sort(vechit_comp.begin(), vechit_comp.end()); //sort per time

                                for(auto const &hit : vechit_comp) {
                                    esum += hit.Energy();
                                    stepLength += hit.StepLength();
                                }

                                //remove the element from the map now to avoid double counting
                                m_Deposits.erase(find->first);
                            }
                            fDeposits.emplace_back( trackID, time, esum, pos, cellID, stepLength );
                        }
                    }

                    //remove the element from the map now
                    m_Deposits.erase(it.first);
                }
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
