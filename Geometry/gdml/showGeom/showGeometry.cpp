#include "TGeoManager.h"
#include "TGLViewer.h"
#include "TPad.h"
#include "TApplication.h"

#include <iostream>
#include <stdlib.h>

void showGeometry(std::string file, int vislevel, bool checkOverlaps)
{
    TGeoManager *geo = new TGeoManager();
    geo->Import(file.c_str());
    geo->DefaultColors();

    if(checkOverlaps){
        geo->CheckOverlaps(1e-5,"d100000000");
        geo->PrintOverlaps();
    }

    geo->SetVisOption(1);
    geo->SetVisLevel(vislevel);

    geo->GetTopVolume()->Draw("ogl");

    TGLViewer * v = (TGLViewer *)gPad->GetViewer3D();
    v->SetStyle(TGLRnrCtx::kOutline);
    v->SetSmoothPoints(kTRUE);
    v->SetLineScale(0.5);
    v->UpdateScene();
}

int main(int argc, char **argv)
{
    if(argc < 4)
    {
        std::cout << "./showGeometry <file> <vislevel: 1 to 7> <checkOverlaps: 0 or 1>" << std::endl;
        return -1;
    }

    bool checkOverlaps = false;
    std::string file = argv[1];
    std::cout << "Geometry file " << file << std::endl;
    int vislevel = std::atoi(argv[2]);
    if(vislevel < 1 || vislevel > 7){
        std::cout << "Set Visibility level to default" << std::endl;
        vislevel = 5;
    }
    std::cout << "VisLevel set to " << vislevel << std::endl;
    std::string checkOverlaps_str = argv[3];
    if( (checkOverlaps_str != "0" && checkOverlaps_str != "1") || checkOverlaps_str == "0" )
    std::cout << "Check for Overlaps set to false" << std::endl;
    else{
        std::cout << "Check for Overlaps set to true" << std::endl;
        checkOverlaps = true;
    }

    TApplication theApp("theApp", &argc, argv) ;
    showGeometry(file, vislevel, checkOverlaps);
    theApp.Run();

    return 0;
}
