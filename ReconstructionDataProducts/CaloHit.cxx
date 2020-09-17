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
        fCellID(0.),
        fLayer(0.)
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
        CaloHit::CaloHit(float energy, float time, float *pos, raw::CellID_t cellID, unsigned int layer)
        : fEnergy  (energy  )
        , fCellID  (cellID  )
        , fLayer   (layer   )
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
        CaloHit::CaloHit(float energy, std::pair<float, float> time, float *pos, raw::CellID_t cellID, unsigned int layer)
        : fEnergy (energy  )
        , fTime   (time    )
        , fCellID (cellID  )
        , fLayer  (layer   )
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
        const unsigned int CaloHit::GetCellLengthScale() const
        {
            std::array<double, 3> point = {this->Position()[0], this->Position()[1], this->Position()[2]};
            if(gar::providerFrom<geo::Geometry>()->isTile(point, this->CellID()))
            {
                return std::sqrt( gar::providerFrom<geo::Geometry>()->getTileSize(point) * gar::providerFrom<geo::Geometry>()->getTileSize(point) );
            }
            else
            {
                return std::sqrt( gar::providerFrom<geo::Geometry>()->getStripWidth(point) * gar::providerFrom<geo::Geometry>()->getStripLength(point, this->CellID()) );
            }
        }

    } // rec
} // gar
