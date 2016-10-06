////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapAlg.h
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_CHANNELMAPALG_H
#define GEO_CHANNELMAPALG_H

// Framework libraries
#include "cetlib/exception.h"

// ROOT libraries
#include "TVector3.h"

// C/C++ standard libraries
#include <vector>
#include <map>
#include <set>

namespace gar{
  namespace geo{
    
    /// Exception thrown on invalid wire number (e.g. NearestWireID())
    class InvalidChannelIDError: public cet::exception {
    public:
      InvalidChannelIDError(std::string const& cat): cet::exception(cat) {}
      
      InvalidChannelIDError(std::string const& cat,
                            int                bad_wire,
                            int                better_wire = -1)
      : cet::exception    (cat)
      , wire_number       (bad_wire)
      , better_wire_number(better_wire)
      {}
      
      int wire_number        = std::numeric_limits<int>::min(); ///< the invalid wire number
      int better_wire_number = std::numeric_limits<int>::min(); ///< a suggestion for a good wire number

    }; // class InvalidWireIDError
    
    
    class ChannelMapAlg{
      
    public:
      
      virtual ~ChannelMapAlg() = default;
      
      virtual void          Initialize()         = 0;
      virtual void          Uninitialize()       = 0;
      virtual unsigned int  Nchannels()    const = 0;
      
    protected:
      
    };
  }
} // namespace gar

#endif // GEO_CHANNELMAPALG_H

