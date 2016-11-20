////////////////////////////////////////////////////////////////////////
/// $Id: SimChannel.h,v 1.3 2010/03/26 20:08:36 brebel Exp $
///
/// \file  SimChannel.h
///
/// \brief object containing MC truth information necessary for making RawDigits 
/// and doing back tracking
///
/// \author  seligman@nevis.columbia.edu
///
////////////////////////////////////////////////////////////////////////

#ifndef SIM_SIMCHANNEL_H
#define SIM_SIMCHANNEL_H

#include <string>
#include <vector>
#include <set>
#include <stdint.h>

namespace gar {
  namespace sdp {
    
    
    struct TrackIDE{
      int   trackID;    ///< Geant4 supplied trackID
      float energyFrac; ///< fraction of hit energy from the particle with this trackID
      float energy;     ///< energy from the particle with this trackID
      TrackIDE() {}
      TrackIDE(int id, float ef, float e) : trackID(id), energyFrac(ef), energy (e) {}
    };
    
    class IDE{
    public:
      
      IDE();
      
#ifndef __GCCXML__
      //constructor for IDEs applying G4 offset...
      IDE(IDE const&, int);
      IDE(int   id,
          float numE,
          float xpos,
          float ypos,
          float zpos,
          float e)
      : trackID     (id)
      , numElectrons(numE)
      , energy      (e)
      , x           (xpos)
      , y           (ypos)
      , z           (zpos)
      {}
      
#endif
      
      int   trackID;      ///< Geant4 supplied track ID
      float numElectrons; ///< total number of electrons for this track ID and time
      float energy;       ///< total energy deposited for this track ID and time
      float x;            ///< x position of ionization
      float y;            ///< y position of ionization
      float z;            ///< z position of ionization
    };
    
    struct TDCIDE
    {
      unsigned short        fTDC;  ///< TDC tick
      std::vector<sdp::IDE> fIDEs; ///< IDEs for this TDC
      
      bool operator<(TDCIDE const& other) const;
    };
    
    class SimChannel
    {
    public:
      
      // Default constructor
      SimChannel();
      
    private:
      
      unsigned int          fChannel; ///< electronics channel associated with these gar::sdp::Electrons
      std::set<sdp::TDCIDE> fTDCIDEs; ///< set of TDCIDE structs
      
      
#ifndef __GCCXML__
      
      TDCIDE const& FindTDCIDE(unsigned short tdc) const;
      
    public:
      
      explicit SimChannel(unsigned int channel);
      
      // method to add ionization electrons and energy to this channel
      void AddIonizationElectrons(int trackID,
                                  unsigned int tdc,
                                  double numberElectrons,
                                  double *xyz,
                                  double energy);
      
      
      
      unsigned int Channel() const;
      
      // method to return a collection of IDE structs for all geant4
      // track ids represented between startTDC and endTDC
      std::vector<sdp::IDE> TrackIDsAndEnergies(unsigned int startTDC,
                                                unsigned int endTDC) const;
      
      const std::set<sdp::TDCIDE>& TDCIDEs() const;
      
      // The number of ionization electrons associated with this channel for the
      // specified TDC.
      double Charge(unsigned int tdc) const;
      double Energy(unsigned int tdc) const;
      
      // A vector of TrackIDEs for a range of TDCs
      std::vector<sdp::TrackIDE> TrackIDEs(unsigned int startTDC,
                                           unsigned int endTDC) const;
      
      bool operator<  (const SimChannel& other) const;
      bool operator== (const SimChannel& other) const;
      
      std::pair<int,int> MergeSimChannel(const SimChannel&, int);
      
      //@{
      /**
       * @brief Dumps the full content of the SimChannel into a stream
       * @param OSTREAM an ostream-line stream object
       * @param out the stream to send the information into
       * @param indent indentation of the lines
       * @param indent_first indentation for the first line (default: as indent)
       */
      template <class OSTREAM>
      void Dump(OSTREAM          & out,
                std::string const& indent,
                std::string const& first_indent) const;
      template <class OSTREAM>
      void Dump(OSTREAM          & out,
                std::string const& indent = "") const { Dump(out, indent, indent); }
      //@}
      
#endif
      
    };
    
  } // namespace sdp
} // gar

#ifndef __GCCXML__

inline bool                              gar::sdp::SimChannel::operator<  (const sdp::SimChannel& other) const { return fChannel < other.Channel();  }
inline bool                              gar::sdp::SimChannel::operator== (const sdp::SimChannel& other) const { return fChannel == other.Channel(); }
inline std::set<gar::sdp::TDCIDE> const& gar::sdp::SimChannel::TDCIDEs()                                 const { return fTDCIDEs; }
inline unsigned int                      gar::sdp::SimChannel::Channel()                                 const { return fChannel; }


// -----------------------------------------------------------------------------
// ---  template implementation
// ---
template <class OSTREAM>
void gar::sdp::SimChannel::Dump(OSTREAM          & out,
                                std::string const& indent,
                                std::string const& first_indent) const
{
  out
  << first_indent
  << "channel #"
  << Channel()
  << " read "
  << fTDCIDEs.size()
  << " TDCs:\n";
  
  double         channel_energy = 0.;
  double         channel_charge = 0.;
  unsigned short tdc            = 0;
  for(auto const& TDCinfo: fTDCIDEs) {
    tdc = TDCinfo.fTDC;
    
    out
    << indent
    << "  TDC #"
    << tdc
    << " with "
    << TDCinfo.fIDEs.size()
    << " IDEs\n";
    
    double tdc_energy = 0., tdc_charge = 0.;
    for (auto const& ide: TDCinfo.fIDEs) {
      out
      << indent
      << "    ("
      << ide.x
      << ", "
      << ide.y
      << ", "
      << ide.z << ") "
      << ide.numElectrons
      << " electrons, "
      << ide.energy
      << " MeV (trkID="
      << ide.trackID
      << ")\n";
      
      tdc_energy += ide.energy;
      tdc_charge += ide.numElectrons;
    } // for IDEs
    
    out
    << indent
    << "    => TDC #"
    << tdc
    << " CH #"
    << Channel()
    << " collected "
    << tdc_energy
    << " electrons and "
    << tdc_energy
    << " MeV\n";
    
    channel_energy += tdc_energy;
    channel_charge += tdc_charge;
  } // for TDCs
  
  out
  << indent
  << "  => channel #"
  << Channel()
  << " collected "
  << channel_charge
  << " electrons and "
  << channel_energy
  << " MeV\n";
} // gar::sdp::SimChannel::Dump<>()

#endif

#endif // sdp_SIMCHANNEL_H

////////////////////////////////////////////////////////////////////////
