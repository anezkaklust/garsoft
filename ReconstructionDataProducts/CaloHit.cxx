//
//  Hit.cxx
//
//  Created by Brian Rebel on 10/6/16.
//
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CoreUtils/ServiceUtil.h"

#include "ReconstructionDataProducts/CaloHit.h"

namespace gar {
    namespace rec {

        //--------------------------------------------------------------------------
        CaloHit::CaloHit()
        : fEnergy(0.),
        fCellID(0.)
        {
            IDNumberGen::create(FirstNumber);
            fIDnumero = IDNumberGen::create()->getNewOne();

            fPosition[0] = 0.;
            fPosition[1] = 0.;
            fPosition[2] = 0.;

            fTime = std::make_pair(0., 0.);

            return;
        }

        bool CaloHit::operator==(const CaloHit& rhs) const {
            return (this->fIDnumero == rhs.fIDnumero);
        }

        bool CaloHit::operator!=(const CaloHit& rhs) const {
            return (this->fIDnumero != rhs.fIDnumero);
        }

        gar::rec::IDNumber CaloHit::getIDNumber() const {return fIDnumero;}


        //--------------------------------------------------------------------------
        CaloHit::CaloHit(float energy, float time, float *pos, raw::CellID_t cellID)
        : fEnergy  (energy  )
        , fCellID     (cellID  )
        {
            IDNumberGen::create(FirstNumber);
            fIDnumero = IDNumberGen::create()->getNewOne();

            fPosition[0] = pos[0];
            fPosition[1] = pos[1];
            fPosition[2] = pos[2];

            fTime = std::make_pair(time, 0.);

            return;
        }

        //--------------------------------------------------------------------------
        CaloHit::CaloHit(float energy, std::pair<float, float> time, float *pos, raw::CellID_t cellID)
        : fEnergy  (energy  )
        , fTime(time)
        , fCellID     (cellID  )
        {
            IDNumberGen::create(FirstNumber);
            fIDnumero = IDNumberGen::create()->getNewOne();

            fPosition[0] = pos[0];
            fPosition[1] = pos[1];
            fPosition[2] = pos[2];

            return;
        }

        //--------------------------------------------------------------------------
        bool CaloHit::operator< (const CaloHit &rhs) const
        {
            const TVector3 rhsPos(rhs.Position()[0], rhs.Position()[1], rhs.Position()[2]);
            const TVector3 thisPos(this->Position()[0], this->Position()[1], this->Position()[2]);
            const TVector3 deltaPosition(rhsPos - thisPos);

            if (std::fabs(deltaPosition.z()) > std::numeric_limits<float>::epsilon())
            return (deltaPosition.z() > std::numeric_limits<float>::epsilon());

            if (std::fabs(deltaPosition.x()) > std::numeric_limits<float>::epsilon())
            return (deltaPosition.x() > std::numeric_limits<float>::epsilon());

            if (std::fabs(deltaPosition.y()) > std::numeric_limits<float>::epsilon())
            return (deltaPosition.y() > std::numeric_limits<float>::epsilon());

            return (this->Energy() > rhs.Energy());
        }

        //--------------------------------------------------------------------------
        std::ostream& operator<< (std::ostream& o, gar::rec::CaloHit const& h)
        {

            o << "CaloHit "
            << "\n\tID number = "
            << h.getIDNumber()
            << "\n\tposition = ("
            << h.Position()[0]
            << ", "
            << h.Position()[1]
            << ", "
            << h.Position()[2]
            << ")"
            << "\n\tenergy = "
            << h.Energy()
            << "\n\t time: "
            << h.Time().first << " " << h.Time().second
            << " cellID: "
            << h.CellID();

            return o;
        }

        //--------------------------------------------------------------------------
        const unsigned int CaloHit::GetDetectorID() const
        {
            gar::geo::GeometryCore const* fGeo = gar::providerFrom<geo::Geometry>();
            return fGeo->getIDbyCellID(this->CellID(), "system");
        }

        //--------------------------------------------------------------------------
        const unsigned int CaloHit::GetLayer() const
        {
            gar::geo::GeometryCore const* fGeo = gar::providerFrom<geo::Geometry>();
            return fGeo->getIDbyCellID(this->CellID(), "layer");
        }

        //--------------------------------------------------------------------------
        const unsigned int CaloHit::GetModule() const
        {
            gar::geo::GeometryCore const* fGeo = gar::providerFrom<geo::Geometry>();
            return fGeo->getIDbyCellID(this->CellID(), "module");
        }

        //--------------------------------------------------------------------------
        const unsigned int CaloHit::GetStave() const
        {
            gar::geo::GeometryCore const* fGeo = gar::providerFrom<geo::Geometry>();
            return fGeo->getIDbyCellID(this->CellID(), "stave");
        }

        //--------------------------------------------------------------------------
        const unsigned int CaloHit::GetCellLengthScale() const
        {
            gar::geo::GeometryCore const* fGeo = gar::providerFrom<geo::Geometry>();

            if(fGeo->isTile(this->CellID()))
            {
                return std::sqrt( fGeo->getTileSize() * fGeo->getTileSize() );
            }
            else
            {
                std::array<double, 3> point = {this->Position()[0], this->Position()[1], this->Position()[2]};
                return std::sqrt( fGeo->getStripWidth() * fGeo->getStripLength(point, this->CellID()) );
            }
        }

    } // rec
} // gar
