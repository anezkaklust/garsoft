#include "Geometry/ChannelMapAlgs/ECALSegmentationAlg.h"

#include "fhiclcpp/ParameterSet.h"

#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iomanip>

namespace gar {
    namespace geo {
        namespace seg {

            //----------------------------------------------------------------------------
            ECALSegmentationAlg::ECALSegmentationAlg(fhicl::ParameterSet const& pset) :
                _name("Segmentation"), _type("Segmentation"), _decoder(new BitFieldCoder(pset.get<std::string>("cellEncoding"))) {
                }

                //----------------------------------------------------------------------------
                ECALSegmentationAlg::ECALSegmentationAlg(const BitFieldCoder* newDecoder, fhicl::ParameterSet const& pset) :
                    _name("Segmentation"), _type("Segmentation"), _decoder(newDecoder) {
                    }

                    //----------------------------------------------------------------------------
                    ECALSegmentationAlg::~ECALSegmentationAlg() {
                        if (_decoder != 0) {
                            delete _decoder;
                        }
                    }

                    //----------------------------------------------------------------------------
                    void ECALSegmentationAlg::setDecoder(const BitFieldCoder* newDecoder) {
                        if ( _decoder == newDecoder )
                        return; //self assignment
                        else {
                            delete _decoder;
                        }
                        _decoder = newDecoder;
                    }

                    //----------------------------------------------------------------------------
                    double ECALSegmentationAlg::binToPosition(gar::raw::CellID_t bin, double cellSize, double offset) {
                        return bin * cellSize + offset;
                    }

                    //----------------------------------------------------------------------------
                    int ECALSegmentationAlg::positionToBin(double position, double cellSize, double offset) {
                        return int(floor((position + 0.5 * cellSize - offset) / cellSize));
                    }

                    //----------------------------------------------------------------------------
                    double ECALSegmentationAlg::getStripLength(const gar::geo::GeometryCore& geo, const gar::raw::CellID_t& cID) const {
                        return 0.;
                    }

                    //----------------------------------------------------------------------------
                    std::pair<TVector3, TVector3> ECALSegmentationAlg::getStripEnds(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const gar::raw::CellID_t& cID) const {
                        return std::make_pair( TVector3(0, 0, 0), TVector3(0, 0, 0) );
                    }

                    //----------------------------------------------------------------------------
                    std::pair<float, float> ECALSegmentationAlg::CalculateLightPropagation(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const gar::raw::CellID_t& cID) const {
                        return std::make_pair( 0., 0. );
                    }

                    //----------------------------------------------------------------------------
                    std::array<double, 3> ECALSegmentationAlg::ReconstructStripHitPosition(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const float &xlocal, const gar::raw::CellID_t& cID) const {
                        return std::array<double, 3>{ {0., 0., 0.} };
                    }
                }//seg
            }//geo
        } //gar
