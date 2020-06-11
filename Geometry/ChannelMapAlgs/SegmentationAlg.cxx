#include "Geometry/ChannelMapAlgs/SegmentationAlg.h"

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
            SegmentationAlg::SegmentationAlg(fhicl::ParameterSet const& pset) :
                _name("Segmentation"), _type("Segmentation"), _encoding(pset.get<std::string>("cellEncoding")), _decoder(new BitFieldCoder(pset.get<std::string>("cellEncoding"))) {
                }

                //----------------------------------------------------------------------------
                SegmentationAlg::SegmentationAlg(const BitFieldCoder* newDecoder, fhicl::ParameterSet const& pset) :
                    _name("Segmentation"), _type("Segmentation"), _encoding(pset.get<std::string>("cellEncoding")), _decoder(newDecoder) {
                    }

                    //----------------------------------------------------------------------------
                    SegmentationAlg::~SegmentationAlg() {
                        if (_decoder != 0) {
                            delete _decoder;
                        }
                    }

                    //----------------------------------------------------------------------------
                    void SegmentationAlg::setDecoder(const BitFieldCoder* newDecoder) {
                        if ( _decoder == newDecoder )
                        return; //self assignment
                        else {
                            delete _decoder;
                        }
                        _decoder = newDecoder;
                    }

                    //----------------------------------------------------------------------------
                    double SegmentationAlg::binToPosition(gar::raw::CellID_t bin, double cellSize, double offset) {
                        return bin * cellSize + offset;
                    }

                    //----------------------------------------------------------------------------
                    int SegmentationAlg::positionToBin(double position, double cellSize, double offset) {
                        return int( floor((position + 0.5 * cellSize - offset) / cellSize) );
                    }

                    //----------------------------------------------------------------------------
                    double SegmentationAlg::getStripLength(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const gar::raw::CellID_t& cID) const {
                        return 0.;
                    }

                    //----------------------------------------------------------------------------
                    std::pair<TVector3, TVector3> SegmentationAlg::getStripEnds(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const gar::raw::CellID_t& cID) const {
                        return std::make_pair( TVector3(0, 0, 0), TVector3(0, 0, 0) );
                    }

                    //----------------------------------------------------------------------------
                    std::pair<float, float> SegmentationAlg::CalculateLightPropagation(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const gar::raw::CellID_t& cID) const {
                        return std::make_pair( 0., 0. );
                    }

                    //----------------------------------------------------------------------------
                    std::array<double, 3> SegmentationAlg::ReconstructStripHitPosition(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const float &xlocal, const gar::raw::CellID_t& cID) const {
                        return std::array<double, 3>{ {0., 0., 0.} };
                    }
                }//seg
            }//geo
        } //gar
