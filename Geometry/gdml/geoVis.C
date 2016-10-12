//
//  geom.C
//  garsoft-mrb
//
//  Created by Brian Rebel on 10/12/16.
//  Copyright © 2016 Brian Rebel. All rights reserved.
//
#include "TGeoManager.h"
#include "TFile.h"
#include "TSystem.h"

typedef struct _drawopt {
  const char* volume;
  int         color;
} drawopt;

void geoVis(TString volName="volCryostat"){
  
  gSystem->Load("libGeom");
  gSystem->Load("libGdml");
  
  TGeoManager::Import("generic_detector.gdml");
  
  drawopt opt[] = {
      //   {"volWorld",                 0},
      //   {"volDetEnclosure",          kWhite},
      //   {"volCryostat",              kOrange},
      //   {"volTPCWirePlaneLengthSide", kCyan+3},
      //   {"volTPCWirePlaneWidthSide", kRed},
    {"volTPCWidthFace",    kRed},
    {"volTPCLengthFace",   kCyan+5},
    {"volTPCBottomFace",   kOrange},
    {"volTubBottom",       kOrange+7},
    {"voltheX",            kOrange+7},
    {"volArgon_solid_L",   kRed},
    {"volArgon_cap_L",     kOrange},
    {"volArgon_cap_front", kOrange},
    {"volTPC",             kOrange},
    {"volDetEnclosure",    kBlue},
    {"volTPCActive",       kGreen},
    {0, 0}
  };
  
  /*for (int i=0;; ++i) {
   if (opt[i].volume==0) break;
   gGeoManager->FindVolumeFast(opt[i].volume)->SetLineColor(opt[i].color);
   }*/
  
//  TList* mat = gGeoManager->GetListOfMaterials();
//  TIter next(mat);
//  TObject *obj;
  // while (obj = next()) {
  //  obj->Print();
  // }

  gGeoManager->CheckOverlaps(0.01);
  gGeoManager->PrintOverlaps();
  gGeoManager->SetMaxVisNodes(70000);

  //gGeoManager->GetTopVolume()->Draw();
  //gGeoManager->FindVolumeFast(volName)->Draw();

  TFile *tf = new TFile("generic_detector.root", "RECREATE");
  gGeoManager->Write();
  tf->Close();
}
