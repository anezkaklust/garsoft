/// \file    GeometryDrawer.cxx
/// \brief   Class to aid in the rendering of Geometry objects
/// \author  messier@indiana.edu
/// \version $Id: GeometryDrawer.cxx,v 1.2 2010/11/11 18:11:22 p-novaart Exp $
#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TLine.h"
#include "TBox.h"
#include "TText.h"
#include "TMath.h"

#include "Geometry/Geometry.h"
#include "nutools/EventDisplayBase/View2D.h"
#include "nutools/EventDisplayBase/View3D.h"
#include "EventDisplay/GeometryDrawer.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace gar {
namespace evd{

  //......................................................................
  GeometryDrawer::GeometryDrawer()
  {
  }

  //......................................................................
  GeometryDrawer::~GeometryDrawer()
  {
  }

  //......................................................................
  void GeometryDrawer::DetOutline3D(evdb::View3D*        view)
  {
    art::ServiceHandle<geo::Geometry> geo;

    double xlo =  -geo->DetLength()/2;
    double xhi =  geo->DetLength()/2;
    double ylo = -geo->DetHalfHeight();
    double yhi =  geo->DetHalfHeight();
    //double zlo =  -geo->DetHalfWidth();
    double zhi =  geo->DetHalfWidth();
    double r=yhi;
    double fracinner = 788.0/2580.;

    double ang = TMath::Pi()*2.0/18.0;

    int c = kGray;
    int s = 1;
    int w = 1;
    TPolyLine3D& spos = view->AddPolyLine3D(19, c, w, s);
    TPolyLine3D& sposi = view->AddPolyLine3D(19, c, w, s);
    for (int i=0;i<19;++i)
      {
	spos.SetPoint(i,xhi,r*TMath::Cos(i*ang),r*TMath::Sin(i*ang));
	sposi.SetPoint(i,xhi,fracinner*r*TMath::Cos(i*ang),fracinner*r*TMath::Sin(i*ang));
      }
    TPolyLine3D& sneg = view->AddPolyLine3D(19, c, w, s);
    TPolyLine3D& snegi = view->AddPolyLine3D(19, c, w, s);
    for (int i=0;i<19;++i)
      {
	sneg.SetPoint(i,xlo,r*TMath::Cos(i*ang),r*TMath::Sin(i*ang));
	snegi.SetPoint(i,xlo,fracinner*r*TMath::Cos(i*ang),fracinner*r*TMath::Sin(i*ang));
      }

    c = kGray+2;
    s = 1;
    w = 1;
    for (int i=0;i<18;++i)
      {
         TPolyLine3D& gridt = view->AddPolyLine3D(4, c, s, w);
	 gridt.SetPoint(0,xlo,fracinner*r*TMath::Cos(i*ang),fracinner*r*TMath::Sin(i*ang));
	 gridt.SetPoint(1,xlo,r*TMath::Cos(i*ang),r*TMath::Sin(i*ang));
	 gridt.SetPoint(2,xhi,r*TMath::Cos(i*ang),r*TMath::Sin(i*ang));
	 gridt.SetPoint(3,xhi,fracinner*r*TMath::Cos(i*ang),fracinner*r*TMath::Sin(i*ang));
      }

    // Indicate coordinate system
    double x0 = -0.20;     // Center location of the key
    double y0 =  1.10*ylo; // Center location of the key
    double z0 = -0.10*zhi; // Center location of the key  
    double sz =  0.20*zhi; // Scale size of the key in z direction

    c = kBlue;

    TPolyLine3D& xaxis = view->AddPolyLine3D(2, c, s, w);
    TPolyLine3D& yaxis = view->AddPolyLine3D(2, c, s, w);
    TPolyLine3D& zaxis = view->AddPolyLine3D(2, c, s, w);
    xaxis.SetPoint(0, x0,    y0, z0);
    xaxis.SetPoint(1, sz+x0, y0, z0);

    yaxis.SetPoint(0, x0, y0,     z0);
    yaxis.SetPoint(1, x0, y0+sz,  z0);

    zaxis.SetPoint(0, x0, y0, z0);
    zaxis.SetPoint(1, x0, y0, z0+sz);

    TPolyLine3D& xpoint = view->AddPolyLine3D(3, c, s, w);
    TPolyLine3D& ypoint = view->AddPolyLine3D(3, c, s, w);
    TPolyLine3D& zpoint = view->AddPolyLine3D(3, c, s, w);
  
    xpoint.SetPoint(0, 0.95*sz+x0, y0, z0-0.05*sz);
    xpoint.SetPoint(1, 1.00*sz+x0, y0, z0);
    xpoint.SetPoint(2, 0.95*sz+x0, y0, z0+0.05*sz);

    ypoint.SetPoint(0, x0, 0.95*sz+y0, z0-0.05*sz);
    ypoint.SetPoint(1, x0, 1.00*sz+y0, z0);
    ypoint.SetPoint(2, x0, 0.95*sz+y0, z0+0.05*sz);

    zpoint.SetPoint(0, x0-0.05*sz, y0, 0.95*sz+z0);
    zpoint.SetPoint(1, x0+0.00*sz, y0, 1.00*sz+z0);
    zpoint.SetPoint(2, x0+0.05*sz, y0, 0.95*sz+z0);

    TPolyLine3D& zleg = view->AddPolyLine3D(4, c, s, w);  
    zleg.SetPoint(0,  x0-0.05*sz, y0+0.05*sz, z0+1.05*sz);
    zleg.SetPoint(1,  x0+0.05*sz, y0+0.05*sz, z0+1.05*sz);
    zleg.SetPoint(2,  x0-0.05*sz, y0-0.05*sz, z0+1.05*sz);
    zleg.SetPoint(3,  x0+0.05*sz, y0-0.05*sz, z0+1.05*sz);

    TPolyLine3D& yleg = view->AddPolyLine3D(5, c, s, w);  
    yleg.SetPoint(0,  x0-0.05*sz, y0+1.15*sz, z0);
    yleg.SetPoint(1,  x0+0.00*sz, y0+1.10*sz, z0);
    yleg.SetPoint(2,  x0+0.00*sz, y0+1.05*sz, z0);
    yleg.SetPoint(3,  x0+0.00*sz, y0+1.10*sz, z0);
    yleg.SetPoint(4,  x0+0.05*sz, y0+1.15*sz, z0);

    TPolyLine3D& xleg = view->AddPolyLine3D(7, c, s, w);  
    xleg.SetPoint(0,  x0+1.05*sz, y0+0.05*sz, z0-0.05*sz);
    xleg.SetPoint(1,  x0+1.05*sz, y0+0.00*sz, z0-0.00*sz);
    xleg.SetPoint(2,  x0+1.05*sz, y0+0.05*sz, z0+0.05*sz);
    xleg.SetPoint(3,  x0+1.05*sz, y0+0.00*sz, z0-0.00*sz);
    xleg.SetPoint(4,  x0+1.05*sz, y0-0.05*sz, z0-0.05*sz);
    xleg.SetPoint(5,  x0+1.05*sz, y0+0.00*sz, z0-0.00*sz);
    xleg.SetPoint(6,  x0+1.05*sz, y0-0.05*sz, z0+0.05*sz);

  }

}
}// namespace
////////////////////////////////////////////////////////////////////////
