#include "garsoft/StandardRecord/SRNeutrino.h"

#include <limits>

// This value will cause an error as soon as it is used, which should help find
// uninitialized variables.
const float kNaN = std::numeric_limits<float>::signaling_NaN();

namespace caf
{
  SRNeutrino::SRNeutrino() :
    pdg(0), ccnc(-1), mode(-1), interactionType(-1),
    Q2(kNaN), W(kNaN), X(kNaN), Y(kNaN), theta(kNaN), T(kNaN),
    gint(-1), tgtPDG(-1), weight(kNaN), gT(kNaN)
  {
  }
}
