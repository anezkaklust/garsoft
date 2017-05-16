/**
 * @file    RawDataDrawer.cxx
 * @brief   Class to aid in the rendering of RawData objects
 * @author  messier@indiana.edu
 * 
 * This class prepares the rendering of the raw digits content for the 2D view
 * of the display. In particular, it fills the 2D view of detected charge
 * vs. time and wire number for a single TPC.
 * The code is not ready to support an arbitrary TPC from the detector;
 * it can be fixed to support that, but a good deal of the calling code also
 * has to change.
 * 
 * Notes from Gianluca Petrillo (petrillo@fnal.gov) on August 19, 2015
 * --------------------------------------------------------------------
 * 
 * As of August 2015, this class performs some preprocessing for the data to
 * be displayed.
 * The main argument here is that the display assigns to each plane a space
 * of 600x250 pixels, while the amound of data is typically 500x5000 elements.
 * These numbers vary wildly: while the first number is close to the current
 * default size for a new window, that window can be resized and windows many
 * times larger can be obtained; for the data content, a MicroBooNE plane
 * contains 2400x9600 elements (that's roughly 25 millions, and there is three
 * of them), each of the two LArIAT planes contains 240x3200 (less than 1
 * million), and a single TPC of DUNE's 35t prototype about 350x4800 (but it
 * is larger when extended readout window modes are used).
 * The code produces TBox'es to be rendered by the nutools event display
 * infrastructure. The code before August 2015 would produce one box for each
 * of the aggregated charge; charge could be aggregated in time by FHiCL
 * configuration to merge TDC ticks. This typically bloats rendering time,
 * since rendering of all boxes is attempted.
 * The new code performs dynamic aggregation after discovering the actual size
 * of the graphical viewport, and it submits at most one TBox per pixel.
 * Additional improvement is caching of the uncompressed raw data, so that
 * following zooming is faster, and especially a way to bypass the decompression
 * when the original data is not compressed in the first place, that saves
 * a bit of time and quite some memory.
 * 
 * @todo There are probably a number of glitches and shortcuts in the current
 * preprocessing implementation. If they become a problem, they can probably be
 * fixed on demand. Examples include:
 * - no alignment of the boxes with the wire and tick numbers
 * - possible border effects
 * - the fact that only the maximum charge is displayed on each box , and no
 *   dynamic charge range detection is in place
 * - the first drawing is performed with a grid that is not the final one
 *   (because for some reason the frame starts larger than it should and is
 *   resized later)
 * - the drawing honours the zoom region the user selected, even if the viewport
 *   ends up being larger (that is, if upstream decides that the zoom region is
 *   too small and a larger area will be drawn, the additional area will be
 *   blank)
 */
#include <cmath>       // std::abs(), ...
#include <utility>     // std::pair<>, std::move()
#include <memory>      // std::unique_ptr()
#include <tuple>
#include <limits>      // std::numeric_limits<>
#include <type_traits> // std::add_const_t<>, ...
#include <typeinfo>    // to use typeid()
#include <algorithm>   // std::fill(), std::find_if(), ...
#include <cstddef>     // std::ptrdiff_t

#include "TH1F.h"
#include "TMarker3DBox.h"
#include "TFrame.h"
#include "TVirtualPad.h"

#include "EventDisplay/ChangeTrackers.h"
#include "EventDisplay/RawDataDrawer.h"
#include "nutools/EventDisplayBase/View3D.h"
#include "nutools/EventDisplayBase/EventHolder.h"
#include "EventDisplay/ColorDrawingOptions.h"
#include "EventDisplay/RawDrawingOptions.h"
#include "RawDataProducts/RawDigit.h"
#include "Geometry/Geometry.h"
#include "DetectorInfo/GArPropertiesService.h"
#include "DetectorInfo/DetectorPropertiesService.h"


#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/demangle.h"

namespace {
    // Utility function to make uniform error messages.
  void writeErrMsg(const char* fcn,
                   cet::exception const& e)
  {
    mf::LogWarning("RecoBaseDrawer") << "RecoBaseDrawer::" << fcn
    << " failed with message:\n"
    << e;
  }
}


namespace gar {
  namespace evd {
    
    //......................................................................
    RawDataDrawer::RawDataDrawer()
    {
    }
    
    //......................................................................
    RawDataDrawer::~RawDataDrawer()
    {
    }
    
    //......................................................................
    void RawDataDrawer::FillQHisto(gar::raw::RawDigit const& dig,
                                   TH1F*                     histo)
    {
      for(size_t t = 0; t < dig.NADC(); ++t)
        histo->Fill(t, 1. * dig.ADC(t) - 1. * dig.Pedestal());

      return;
    }//end loop over raw hits

    //......................................................................
    void RawDataDrawer::RawDigit3D(const art::Event& evt,
                                   evdb::View3D*     view)
    {
      // Check if we're supposed to draw raw hits at all
      art::ServiceHandle<evd::RawDrawingOptions>   rawopt;
      art::ServiceHandle<evd::ColorDrawingOptions> cdopt;
      art::ServiceHandle<geo::Geometry>            geom;
      
      if(rawopt->fDrawRawOrReco > 0) return;
      
      std::vector<const raw::RawDigit*> rawhits;
      this->GetRawDigits(evt, rawhits);
      
      fEventTQHist->Reset();
      fDigitTQHist->Reset();
      
      auto const& colorSet = cdopt->RawQ();
      
      for(auto dig : rawhits){

        this->FillQHisto(*dig, fEventTQHist);
        if(dig->Channel() == rawopt->fChannel)
          this->FillQHisto(*dig, fDigitTQHist);
        
        // draw the digit as time, y, z with the latter coordinates
        // given by the channel map.
        float xyz[3];
        
        geom->ChannelToPosition((unsigned int)dig->Channel(),
                                xyz);
        
        // loop over the ADC values for this digit and draw the value if it
        // is above the threshold
        TMarker3DBox *box = nullptr;
        for(size_t t = 0; t < dig->NADC(); ++t){
          if(dig->ADC(t) < rawopt->fMinSignal) continue;
          
          xyz[0] = 1. * t;
          box = &(view->AddMarker3DBox(xyz[0],
                                          xyz[1],
                                          xyz[2],
                                          0.5,    // the extent is 1/2 tick
                                          0.5 * geom->ChannelPitch(),
                                          0.5 * geom->ChannelPitch()));
          
          box->SetFillStyle(1001);
          box->SetFillColor(colorSet.GetColor(dig->ADC(t)));
          box->SetBit(kCannotPick);
        } // end loop over ADC values for the digit
      }//end loop over raw digits
      
      return;
    }

    //......................................................................
    void RawDataDrawer::GetRawDigits(art::Event                        const& evt,
                                     std::vector<const raw::RawDigit*>      & digits)
    {
      digits.clear();
      
      art::ServiceHandle<evd::RawDrawingOptions> rawopt;

      std::vector<const raw::RawDigit*> temp;
      
      try{
        evt.getView(rawopt->fRawDataLabel, temp);
        for(size_t t = 0; t < temp.size(); ++t){
          digits.push_back(temp[t]);
        }
      }
      catch(cet::exception& e){
        writeErrMsg("GetDigits", e);
      }
      
      return;
    } // RawDataDrawer::GetRawDigits()
    
  } // namespace evd
}

////////////////////////////////////////////////////////////////////////
