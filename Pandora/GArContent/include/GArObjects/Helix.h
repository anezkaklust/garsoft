#ifndef GARCONTENT_HELIX_H
#define GARCONTENT_HELIX_H 1

#include "Objects/CartesianVector.h"

#include "Pandora/StatusCodes.h"

namespace gar_content
{

    /**
    *  @brief  Helix class
    */
    class Helix
    {
    public:
        /**
        *  @brief  Constructor using canonical (LEP-wise) parameterisation
        *
        *  @param  phi0 phi angle of momentum vector at the point of closest approach to IP in R-Phi plane
        *  @param  d0 signed distance of closest approach in R-Phi plane
        *  @param  x0 x coordinate of the point of closest approach to IP in R-Phi plane
        *  @param  omega signed curvature
        *  @param  tanLambda tangent of dip angle
        *  @param  bField magnetic field (in Tesla)
        */
        Helix(const float phi0, const float d0, const float x0, const float omega, const float tanlambda, const float bField);

        /**
        *  @brief  Constructor
        *
        *  @param  position position of the reference point
        *  @param  momentum momentum vector at the reference point
        *  @param  charge particle charge
        *  @param  bField magnetic field (in Tesla)
        */
        Helix(const pandora::CartesianVector &position, const pandora::CartesianVector &momentum, const float charge, const float bField);

        /**
        *  @brief  Get helix intersection point with a plane parallel to x axis. The plane is defined by two coordinates in the
        *          plane (z0,y0) and a normal vector (az,ay).
        *
        *  @param  z0 z coordinate in the specified plane
        *  @param  y0 y coordinate in the specified plane
        *  @param  az z component of vector normal to specified plane
        *  @param  ay y component of vector normal to specified plane
        *  @param  referencePoint the reference point of the helix
        *  @param  intersectionPoint to receive the coordinates of the intersection point
        */
        pandora::StatusCode GetPointInZY(const float z0, const float y0, const float az, const float ay, const pandora::CartesianVector &referencePoint,
        pandora::CartesianVector &intersectionPoint) const;

        /**
        *  @brief  Get helix intersection point with a plane parallel to x axis. The plane is defined by two coordinates in the
        *          plane (z0,y0) and a normal vector (az,ay).
        *
        *  @param  z0 z coordinate in the specified plane
        *  @param  y0 y coordinate in the specified plane
        *  @param  az z component of vector normal to specified plane
        *  @param  ay y component of vector normal to specified plane
        *  @param  referencePoint the reference point of the helix
        *  @param  intersectionPoint to receive the coordinates of the intersection point
        *  @param  genericTime to receive the generic time (helix length, from reference point to intersection, divided by particle momentum)
        */
        pandora::StatusCode GetPointInZY(const float z0, const float y0, const float az, const float ay, const pandora::CartesianVector &referencePoint,
        pandora::CartesianVector &intersectionPoint, float &genericTime) const;

        /**
        *  @brief  Get helix intersection point with a plane perpendicular to x axis.
        *
        *  @param  xPlane the x coordinate for the specified plane
        *  @param  referencePoint the reference point of the helix
        *  @param  intersectionPoint to receive the coordinates of the intersection point
        */
        pandora::StatusCode GetPointInX(const float xPlane, const pandora::CartesianVector &referencePoint, pandora::CartesianVector &intersectionPoint) const;

        /**
        *  @brief  Get helix intersection point with a plane perpendicular to x axis.
        *
        *  @param  xPlane the x coordinate for the specified plane
        *  @param  referencePoint the reference point of the helix
        *  @param  intersectionPoint to receive the coordinates of the intersection point
        *  @param  genericTime to receive the generic time (helix length, from reference point to intersection, divided by particle momentum)
        */
        pandora::StatusCode GetPointInX(const float xPlane, const pandora::CartesianVector &referencePoint, pandora::CartesianVector &intersectionPoint,
        float &genericTime) const;

        /**
        *  @brief  Get coordinates of helix intersection with cylinder, aligned along z-axis
        *
        *  @param  radius the radius of the cylinder
        *  @param  referencePoint the reference point of the helix
        *  @param  intersectionPoint to receive the coordinates of the intersection point
        */
        pandora::StatusCode GetPointOnCircle(const float radius, const pandora::CartesianVector &referencePoint, pandora::CartesianVector &intersectionPoint) const;

        /**
        *  @brief  Get coordinates of helix intersection with cylinder, aligned along z-axis
        *
        *  @param  radius the radius of the cylinder
        *  @param  referencePoint the reference point of the helix
        *  @param  intersectionPoint to receive the coordinates of the intersection point
        *  @param  genericTime to receive the generic time (helix length, from reference point to intersection, divided by particle momentum)
        */
        pandora::StatusCode GetPointOnCircle(const float radius, const pandora::CartesianVector &referencePoint, pandora::CartesianVector &intersectionPoint,
        float &genericTime) const;

        /**
        *  @brief  Get distance of the closest approach of helix to an arbitrary point in space
        *
        *  @param  point coordinates of the specified point
        *  @param  distance to receive a vector of distances from helix to point in the following projections:
        *          x component: distance in R-Phi plane
        *          y-component: distance along Z axis
        *          z-component: 3D distance magnitude
        */
        pandora::StatusCode GetDistanceToPoint(const pandora::CartesianVector &point, pandora::CartesianVector &distance) const;

        /**
        *  @brief  Get distance of the closest approach of helix to an arbitrary point in space
        *
        *  @param  point coordinates of the specified point
        *  @param  distance to receive a vector of distances from helix to point in the following projections:
        *          x component: distance in R-Phi plane
        *          y-component: distance along Z axis
        *          z-component: 3D distance magnitude
        *  @param  genericTime to receive the generic time (helix length, from reference point to intersection, divided by particle momentum)
        */
        pandora::StatusCode GetDistanceToPoint(const pandora::CartesianVector &point, pandora::CartesianVector &distance, float &genericTime) const;

        /**
        *  @brief  Get distance between two helices
        *
        *  @param  pHelix address of a second helix
        *  @param  positionOfClosestApproach to receive position of the point of closest approach
        *  @param  v0momentum to receive the v0 momentum
        *  @param  helixDistance to receive the distance between the two helices
        */
        pandora::StatusCode GetDistanceToHelix(const Helix *const pHelix, pandora::CartesianVector &positionOfClosestApproach, pandora::CartesianVector &v0momentum,
        float &helixDistance) const;

        /**
        *  @param  Get extrapolated momentum at a specified position
        *
        *  @param  position the specified position
        *
        *  @return the extrapolated momentum
        */
        pandora::CartesianVector GetExtrapolatedMomentum(const pandora::CartesianVector &position) const;

        /**
        *  @brief  Get momentum of particle at the point of closest approach to IP
        *
        *  @return the momentum of particle
        */
        const pandora::CartesianVector &GetMomentum() const;

        /**
        *  @brief  Get reference point of track
        *
        *  @return the reference point of track
        */
        const pandora::CartesianVector &GetReferencePoint() const;

        /**
        *  @brief  Get phi angle of the momentum vector at the point of closest approach to IP
        *
        *  @return the phi angle of the momentum vector
        */
        float GetPhi0() const;

        /**
        *  @brief  Get x signed distance of closest approach to IP in the R-Phi plane
        *
        *  @return the signed distance of closest approach
        */
        float GetD0() const;

        /**
        *  @brief  Get x coordinate of the point of closest approach to IP in the R-Phi plane
        *
        *  @return the x coordinate of the point of closest approach
        */
        float GetX0() const;

        /**
        *  @brief  Get signed curvature of the track
        *
        *  @return the signed curvature of the track
        */
        float GetOmega() const;

        /**
        *  @brief  Get tangent of dip angle of the track
        *
        *  @return the tangent of dip angle of the track
        */
        float GetTanLambda() const;

        /**
        *  @brief  Get transverse momentum of the track
        *
        *  @return the transverse momentum of the track
        */
        float GetPzy() const;

        /**
        *  @brief  Get charge
        *
        *  @return the charge
        */
        float GetCharge() const;

        /**
        *  @brief  Get z coordinate of circumference
        *
        *  @return the z coordinate of circumference
        */
        float GetZCentre() const;

        /**
        *  @brief  Get y coordinate of circumference
        *
        *  @return the y coordinate of circumference
        */
        float GetYCentre() const;

        /**
        *  @brief  Get radius of circumference
        *
        *  @return the radius of circumference
        */
        float GetRadius() const;

    private:
        static const float FCT;
        static const float TWO_PI;
        static const float HALF_PI;

        pandora::CartesianVector     m_referencePoint;       ///< The coordinates of the reference point
        pandora::CartesianVector     m_momentum;             ///< The momentum vector at reference point

        float               m_phi0;                 ///< phi0 in canonical parameterization
        float               m_d0;                   ///< d0 in canonical parameterisation
        float               m_x0;                   ///< x0 in canonical parameterisation
        float               m_omega;                ///< signed curvature in canonical parameterisation
        float               m_tanLambda;            ///< tanLambda
        float               m_pzy;                  ///< The transverse momentum
        float               m_charge;               ///< The particle charge
        float               m_zCentre;              ///< The circle centre z coordinate
        float               m_yCentre;              ///< The circle centre y coordinate
        float               m_radius;               ///< The radius of circle in ZY plane
        float               m_phiRefPoint;          ///< Phi w.r.t. (Z0, Y0) of circle at reference point
        float               m_phiAtPCA;             ///< Phi w.r.t. (Z0, Y0) of circle at point of closest approach
        float               m_zAtPCA;               ///< z coordinate at point of closest approach
        float               m_yAtPCA;               ///< y coordinate at point of closest approach
        float               m_pzAtPCA;              ///< Momentum z component at point of closest approach
        float               m_pyAtPCA;              ///< Momentum y component at point of closest approach
        float               m_phiMomRefPoint;       ///< Phi of Momentum vector at reference point
    };

    //------------------------------------------------------------------------------------------------------------------------------------------

    inline pandora::StatusCode Helix::GetPointInZY(const float z0, const float y0, const float az, const float ay, const pandora::CartesianVector &referencePoint,
    pandora::CartesianVector &intersectionPoint) const
    {
        float genericTime;
        return this->GetPointInZY(z0, y0, az, ay, referencePoint, intersectionPoint, genericTime);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    inline pandora::StatusCode Helix::GetPointInX(const float xPlane, const pandora::CartesianVector &referencePoint, pandora::CartesianVector &intersectionPoint) const
    {
        float genericTime;
        return this->GetPointInX(xPlane, referencePoint, intersectionPoint, genericTime);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    inline pandora::StatusCode Helix::GetPointOnCircle(const float radius, const pandora::CartesianVector &referencePoint, pandora::CartesianVector &intersectionPoint) const
    {
        float genericTime;
        return this->GetPointOnCircle(radius, referencePoint, intersectionPoint, genericTime);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    inline pandora::StatusCode Helix::GetDistanceToPoint(const pandora::CartesianVector &point, pandora::CartesianVector &distance) const
    {
        float genericTime;
        return this->GetDistanceToPoint(point, distance, genericTime);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    inline const pandora::CartesianVector &Helix::GetMomentum() const
    {
        return m_momentum;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    inline const pandora::CartesianVector &Helix::GetReferencePoint() const
    {
        return m_referencePoint;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    inline float Helix::GetPhi0() const
    {
        if (m_phi0 < 0.)
        return m_phi0 + TWO_PI;

        return m_phi0;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    inline float Helix::GetD0() const
    {
        return m_d0;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    inline float Helix::GetX0() const
    {
        return m_x0;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    inline float Helix::GetOmega() const
    {
        return m_omega;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    inline float Helix::GetTanLambda() const
    {
        return m_tanLambda;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    inline float Helix::GetPzy() const
    {
        return m_pzy;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    inline float Helix::GetCharge() const
    {
        return m_charge;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    inline float Helix::GetZCentre() const
    {
        return m_zCentre;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    inline float Helix::GetYCentre() const
    {
        return m_yCentre;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    inline float Helix::GetRadius() const
    {
        return m_radius;
    }

} // namespace gar_content

#endif // #ifndef GARCONTENT_HELIX_H
