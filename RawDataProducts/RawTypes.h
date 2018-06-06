#ifndef GAR_RAWTYPES_H
#define GAR_RAWTYPES_H

#include <vector>
#include <cstdint>
#include <cstddef>

namespace gar
{
  namespace raw
  {
    typedef short ADC_t;
    typedef std::vector<ADC_t> ADCvector_t;

    typedef enum _compress {
      kNone,       ///< no compression 
      kHuffman,    ///< Huffman Encoding
      kZeroSuppression,  ///< Zero Suppression algorithm
      kZeroHuffman,  ///< Zero Suppression followed by Huffman Encoding
      kDynamicDec  ///< Dynamic decimation
    } Compress_t;
  }
}
#endif // GAR_RAWTYPES_H
