#ifndef ECALSEGMENTATIONALG_H
#define ECALSEGMENTATIONALG_H

#include <vector>
#include <map>

#include "Geometry/GeometryCore.h"
#include "Geometry/BitFieldCoder.h"
#include "RawDataProducts/CaloRawDigit.h"

#include "TVector3.h"

namespace fhicl{
    class ParameterSet;
}

namespace gar {
    namespace geo {

        class GeometryCore;

        typedef unsigned long long ulong64 ;

        class ECALSegmentationAlg {

        public:
            virtual ~ECALSegmentationAlg();

            virtual const std::string& name() const {
                return _name;
            }

            virtual void setName(const std::string& value) {
                _name = value;
            }

            virtual const std::string& type() const {
                return _type;
            }

            virtual const std::string& description() const {
                return _description;
            }

            virtual const BitFieldCoder* decoder()  const {
                return _decoder;
            }

            virtual const double& gridSizeX() const {
                return _gridSizeX;
            }

            virtual const double& stripSizeX() const {
                return _stripSizeX;
            }

            virtual void setDecoder(const BitFieldCoder* decoder);

            //Pure virtual member functions
            virtual void reconfigure(fhicl::ParameterSet const& pset) = 0;

            virtual void Initialize(const gar::geo::GeometryCore& geo) = 0;

            virtual std::array<double, 3> GetPosition(const gar::geo::GeometryCore& geo, const long64& cID) const = 0;

            virtual raw::CellID_t GetCellID(const gar::geo::GeometryCore& geo, const unsigned int& det_id, const unsigned int& stave, const unsigned int& module, const unsigned int& layer, const unsigned int& slice, const std::array<double, 3>& localPosition) const = 0;

            virtual int getIDbyCellID(const raw::CellID_t& cID, const char* id) const = 0;

            virtual void PrintParameters() const = 0;

            virtual bool isTile(const long long int& cID) const = 0;

            virtual void setLayerDimXY(const double& dimX, const double& dimY) const = 0;

            //Non-pure virtual member functions
            virtual double getStripLength(const gar::geo::GeometryCore& geo, const raw::CellID_t& cID) const;

            virtual std::pair<TVector3, TVector3> getStripEnds(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const raw::CellID_t& cID) const;

            virtual std::pair<float, float> CalculateLightPropagation(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const raw::CellID_t& cID) const;

            virtual std::array<double, 3> ReconstructStripHitPosition(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const float &xlocal, const raw::CellID_t& cID) const;

        protected:
            ECALSegmentationAlg(fhicl::ParameterSet const& pset);

            ECALSegmentationAlg(const BitFieldCoder* decoder, fhicl::ParameterSet const& pset);

            static double binToPosition(raw::CellID_t bin, double cellSize, double offset = 0);

            static int positionToBin(double position, double cellSize, double offset = 0);

            std::string _name;

            std::string _type;

            std::string _description;

            double _gridSizeX;

            double _stripSizeX;

            const BitFieldCoder* _decoder = 0;

        private:
            ECALSegmentationAlg(const ECALSegmentationAlg&);
        };
    }
}

#endif
