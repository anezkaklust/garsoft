#include "Geometry/ChannelMapAlgs/SegmentationGridXYAlg.h"

#include "fhiclcpp/ParameterSet.h"

#include <algorithm>
#include <sstream>
#include <iostream>
#include <cmath>

namespace gar {
    namespace geo {
        namespace seg {

            //----------------------------------------------------------------------------
            /// Default constructor used by derived classes passing the encoding string
            SegmentationGridXYAlg::SegmentationGridXYAlg(fhicl::ParameterSet const& pset)
            : SegmentationAlg(pset)
            {
                _type = "GridXY";
                _description = "Cartesian segmentation in the local XY-plane";

                std::cout << " ######### gar::geo::seg::SegmentationGridXYAlg() " << std::endl ;

                this->reconfigure(pset);
            }

            //----------------------------------------------------------------------------
            /// Default constructor used by derived classes passing an existing decoder
            SegmentationGridXYAlg::SegmentationGridXYAlg(const BitFieldCoder* decode, fhicl::ParameterSet const& pset)
            : SegmentationAlg(decode, pset)
            {
                _type = "GridXY";
                _description = "Cartesian segmentation in the local XY-plane";

                std::cout << " ######### gar::geo::seg::SegmentationGridXYAlg() " << std::endl ;

                this->reconfigure(pset);
            }

            //----------------------------------------------------------------------------
            SegmentationGridXYAlg::~SegmentationGridXYAlg()
            {
            }

            //----------------------------------------------------------------------------
            void SegmentationGridXYAlg::reconfigure(fhicl::ParameterSet const& pset)
            {
                _xId = pset.get<std::string>("identifier_x");
                _yId = pset.get<std::string>("identifier_y");
                _encoding = pset.get<std::string>("cellEncoding");

                _gridSizeX = pset.get<double>("grid_size_x");
                _gridSizeY = pset.get<double>("grid_size_y");

                _offsetX = pset.get<double>("offset_x");
                _offsetY = pset.get<double>("offset_y");

                _nLayers = pset.get<unsigned int>("nlayers");

                PrintParameters();

                return;
            }

            //----------------------------------------------------------------------------
            void SegmentationGridXYAlg::Initialize(const gar::geo::GeometryCore& )
            {

            }

            //----------------------------------------------------------------------------
            std::array<double, 3> SegmentationGridXYAlg::GetPosition(const gar::geo::GeometryCore& , const gar::raw::CellID_t& cID) const
            {
                std::array<double, 3> cellPosition;
                cellPosition[0] = binToPosition(_decoder->get(cID, _xId), _gridSizeX, _offsetX);
                cellPosition[1] = binToPosition(_decoder->get(cID, _yId), _gridSizeY, _offsetY);
                cellPosition[2] = 0.;

                return cellPosition;
            }

            //----------------------------------------------------------------------------
            /// determine the cell ID based on the position
            gar::raw::CellID_t SegmentationGridXYAlg::GetCellID(const gar::geo::GeometryCore& , const unsigned int& det_id, const unsigned int& stave, const unsigned int& module, const unsigned int& layer, const unsigned int& slice, const std::array<double, 3>& localPosition) const
            {
                gar::raw::CellID_t cID = 0;

                _decoder->set(cID, "system", det_id);
                _decoder->set(cID, "module", module);
                _decoder->set(cID, "stave", stave);
                _decoder->set(cID, "layer", layer);
                _decoder->set(cID, "slice", slice);

                _decoder->set(cID, _xId, positionToBin(localPosition[0], _gridSizeX, _offsetX));
                _decoder->set(cID, _yId, positionToBin(localPosition[1], _gridSizeY, _offsetY));

                return cID;
            }

            //----------------------------------------------------------------------------
            bool SegmentationGridXYAlg::isBarrel(const gar::raw::CellID_t& cID) const
            {
                bool isBarrel = true;
                // The octogonal ECAL geometries have a module 6.  Is this code only for
                // that particular case?
                int det_id = _decoder->get(cID, "system");
                int module = _decoder->get(cID, "module");
                if( det_id == 2 && (module == 0 || module == 6) ) isBarrel = false;

                return isBarrel;
            }

            //----------------------------------------------------------------------------
            void SegmentationGridXYAlg::PrintParameters() const
            {
                std::cout << "cell encoding: " << _encoding << std::endl;
                std::cout << "identifier_x: " << _xId << std::endl;
                std::cout << "identifier_y: " << _yId << std::endl;
                std::cout << "grid_size_x: " << _gridSizeX << " cm" << std::endl;
                std::cout << "grid_size_y: " << _gridSizeY << " cm" << std::endl;
                std::cout << "offset_x: " << _offsetX << " cm" << std::endl;
                std::cout << "offset_y: " << _offsetY << " cm" << std::endl;
            }

        }//seg
    } //geo
} //gar
