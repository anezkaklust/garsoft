#include "Utils.h"

Utils::Utils()
: _seed(0), _rando(new TRandom3(0))
{
    _origin[0] = 0.;
    _origin[1] = 0.;
    _origin[2] = 0.;
}

Utils::~Utils()
{
    delete _rando;
}

void Utils::SetSeed(int seed)
{
    _seed = seed;
    _rando->SetSeed(_seed);
}

void Utils::SetOrigin(double *origin)
{
    _origin[0] = origin[0];
    _origin[1] = origin[1];
    _origin[2] = origin[2];
}

bool Utils::PointInFiducial(TVector3 point)
{
    //TPC Fiducial volume defined as
    //R < 260 cm
    //|X| < 250 cm
    bool isInFiducial = true;

    float r_point = std::sqrt( point.Y()*point.Y() + point.Z()*point.Z() );
    if( r_point > _TPCFidRadius ) isInFiducial = false;
    if( r_point < _TPCFidRadius && std::abs(point.X()) > _TPCFidLength ) isInFiducial = false;

    return isInFiducial;
}

bool Utils::PointInTPC(TVector3 point)
{
    //TPC volume defined as
    //R < 260 cm
    //|X| < 250 cm
    if(PointInFiducial(point)) return true;
    bool isInTPC = true;

    float r_point = std::sqrt( point.Y()*point.Y() + point.Z()*point.Z() );
    if( r_point > _TPCRadius ) isInTPC = false;
    if( r_point < _TPCRadius && std::abs(point.X()) > _TPCLength ) isInTPC = false;

    return isInTPC;
}

bool Utils::PointInCalo(TVector3 point)
{
    //Barrel Radius 278 cm
    //Endcap starts at 364 cm
    bool isInCalo = false;
    float r_point = std::sqrt( point.Y()*point.Y() + point.Z()*point.Z() );
    //in the Barrel
    if( r_point > _ECALInnerRadius && r_point < _ECALOuterRadius && std::abs(point.X()) < _ECALStartX ) isInCalo = true;
    //in the Endcap
    if( r_point < _ECALInnerRadius && std::abs(point.X()) > _ECALStartX && std::abs(point.X()) < _ECALEndX ) isInCalo = true;

    return isInCalo;
}

bool Utils::isThroughCalo(TVector3 point)
{
    return !PointInTPC(point) && !PointInCalo(point);
}

bool Utils::hasDecayedInCalo(TVector3 epoint)
{
    return PointInCalo(epoint);
}

bool Utils::isBackscatter(TVector3 spoint, TVector3 epoint)
{
    bool isBackscatter = false;

    //check if started in the ECAL but made it to the tracker
    float r_spoint = std::sqrt( spoint.Y()*spoint.Y() + spoint.Z()*spoint.Z() );
    float r_epoint = std::sqrt( epoint.Y()*epoint.Y() + epoint.Z()*epoint.Z() );

    //in the Barrel
    if( r_spoint > _ECALInnerRadius && r_epoint < _TPCRadius ) isBackscatter = true;
    //in the Endcap
    if( (r_spoint < _ECALInnerRadius && r_epoint < _TPCRadius) && ( std::abs(spoint.X()) > _ECALStartX && std::abs(epoint.X()) < _TPCLength ) ) isBackscatter = true;

    return isBackscatter;
}

bool Utils::isBremsstrahlung(TVector3 spoint, int pdg, int motherpdg)
{
    bool isBremsstrahlung = false;

    //Check if it has origin in the tracker and that the pdg is photon and mother is electron/positron (most probable)
    if(PointInFiducial(spoint) && pdg == 22 && std::abs(motherpdg) == 11) isBremsstrahlung = true;

    return isBremsstrahlung;
}
