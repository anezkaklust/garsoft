/// \file    raw.cxx
/// \brief   raw data utilities
/// \author  trj@fnal.gov
/// with thanks to Brian Rebel and Jonathan Insler

#include "RawDataProducts/raw.h"
#include "RawDataProducts/RawTypes.h"

#include <iostream>
#include <string>
#include <bitset>
#include <numeric> // std::adjacent_difference()
#include <iterator> // std::back_inserter()

#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace gar {
  namespace raw {

    //----------------------------------------------------------
    // no arguments -- can only run Huffman for now (or other future parameterless compression methods)
    void Compress(gar::raw::ADCvector_t    &adc, 
		  gar::raw::Compress_t     compress)
    {
      if(compress == gar::raw::kHuffman) CompressHuffman(adc);
      else if(compress == raw::kZeroHuffman){
	throw cet::exception("gar::raw")
	  << "Compress method called for kZeroHuffman but no threshold or pedestal value";
      }
      else { 
	throw cet::exception("gar::raw")
	  << "Compress method called with compression type: " << compress << " but no threshold or pedestal value";
      }
      return;
    }

    //----------------------------------------------------------
    int Compress(gar::raw::ADCvector_t &adc, 
		  gar::raw::Compress_t  compress, 
		  gar::raw::ADC_t       zerothreshold,
		  size_t                ticksbefore,
		  size_t                ticksafter)
    {
      int retval = 1;
      if(compress == raw::kHuffman) 
	{
	  CompressHuffman(adc);
	}
      else if(compress == raw::kZeroSuppression) 
	{
	  retval = ZeroSuppression(adc,zerothreshold,ticksbefore,ticksafter);
	}
      else if(compress == raw::kZeroHuffman)
	{
	  retval = ZeroSuppression(adc,zerothreshold,ticksbefore,ticksafter);
	  CompressHuffman(adc);
        }    

      return retval;
    }


    //----------------------------------------------------------
    // Zero suppression function
    int ZeroSuppression(gar::raw::ADCvector_t &adc, 
			 gar::raw::ADC_t       zerothreshold,
			 size_t                ticksbefore_in,
			 size_t                ticksafter_in)
    {
      const size_t adcsize = adc.size();

      size_t ticksbefore = ticksbefore_in;
      if (adcsize < ticksbefore_in) ticksbefore = adcsize;
      size_t ticksafter = ticksafter_in;
      if (adcsize < ticksafter_in) ticksafter = adcsize;

      gar::raw::ADCvector_t zerosuppressed(adcsize);
      size_t maxblocks = adcsize/2 + 1;
      gar::raw::ADCvector_t blockbegin(maxblocks);  // use the ADC data structure to save TOC info
      gar::raw::ADCvector_t blocksize(maxblocks);

      size_t nblocks = 0;
      size_t zerosuppressedsize = 0;
 
      bool inablock = false;
      bool taflag = false;   // an ADC over threshold is coming up within ticksbefore of the current point
      size_t nta = 0;        // count of how many ticks have gone by since seeing one over threshold ticksbefore in the future
      bool tbflag = false;
      size_t ntb = 0;

      // get us started -- see if there are any adc's over threshold in the first ticksbefore+1 samples
      for (size_t i=0; i<ticksbefore; ++i)
	{
	  if (adc[i] > zerothreshold)
	    {
	      taflag = true;
	      break;
	    }
	}

      for(size_t i = 0; i < adcsize; ++i){
	gar::raw::ADC_t adc_current_value = adc[i];
    
	// determine if there is an ADC over threshold within nticksbefore of the current point in the future

	if ( i < adcsize-ticksbefore )
	  {
	    if (adc[i+ticksbefore] > zerothreshold) 
	      {
		taflag = true;
	      }
	    else
	      {
		++nta;
	      }
	  }
	if (nta > ticksbefore) 
	  {
	    taflag = false;
	    nta = 0;
	  }

	if ( adc_current_value > zerothreshold )
	  {
	    ntb = ticksafter + 1;  // need to count the tick we're on
	  }
	else
	  {
	    if (ntb>0) --ntb;
	  }
	tbflag = (ntb>0);

	if( taflag || tbflag )
	  {
	    if(!inablock){
	      blockbegin[nblocks] = i;
	      blocksize[nblocks] = 0;
	      inablock = true;
	    }
      	    zerosuppressed[zerosuppressedsize] = adc[i];
	    zerosuppressedsize++;
	    blocksize[nblocks]++;		
	  }
	else
	  {
	    if (inablock)
	      {
		nblocks++;  
		inablock = false;
	      }
	  } 
      }
      if (inablock)
	{
	  nblocks++;  
	}

      adc.resize(2+nblocks+nblocks+zerosuppressedsize);

      adc[0] = adcsize; //fill first entry in adc with length of uncompressed vector
      adc[1] = nblocks;

    
      for(unsigned int i = 0; i < nblocks; ++i)
	adc[i+2] = blockbegin[i];

      for(unsigned int i = 0; i < nblocks; ++i)
	adc[i+nblocks+2] = blocksize[i];

      for(unsigned int i = 0; i < zerosuppressedsize; ++i)
	adc[i+nblocks+nblocks+2] = zerosuppressed[i];
 
      return nblocks;  // for use in discarding rawdigit in case it's all zero
    }


    //----------------------------------------------------------
    // Reverse zero suppression function, filling in suppressed ADC values with pedestal

    void ZeroUnsuppression(const gar::raw::ADCvector_t   &adc, 
			   gar::raw::ADCvector_t         &uncompressed,
			   gar::raw::ADC_t               pedestal)
    {
      const int lengthofadc = adc[0];
      const int nblocks = adc[1];

      uncompressed.resize(lengthofadc);
      for (int i = 0;i < lengthofadc; ++i){
	uncompressed[i] = pedestal;
      }
    
      int zerosuppressedindex = nblocks*2 + 2;

      for(int i = 0; i < nblocks; ++i){ //loop over each nonzero block of the compressed vector
      
	for(int j = 0; j < adc[2+nblocks+i]; ++j){//loop over each block size

	  //set uncompressed value
	  uncompressed[adc[2+i]+j] = adc[zerosuppressedindex];
	  zerosuppressedindex++;

	}
      }

      return;
    }


    //----------------------------------------------------------
    // if the compression type is kNone, copy the adc vector into the uncompressed vector
    void Uncompress(const gar::raw::ADCvector_t& adc, 
		    gar::raw::ADCvector_t      &uncompressed, 
		    gar::raw::Compress_t          compress)
    {
      if(compress == raw::kHuffman) UncompressHuffman(adc, uncompressed);
      else if(compress == raw::kNone){
	for(unsigned int i = 0; i < adc.size(); ++i) uncompressed[i] = adc[i];
      }
      else {
	throw cet::exception("raw")
	  << "raw::Uncompress() does not support compression #"
	  << ((int) compress) << " without a pedestal specified";
      }
      return;

    }
  
    //----------------------------------------------------------
    // if the compression type is kNone, copy the adc vector into the uncompressed vector
    void Uncompress(const gar::raw::ADCvector_t& adc, 
		    gar::raw::ADCvector_t      &uncompressed, 
		    ADC_t               pedestal,
		    gar::raw::Compress_t          compress)
    {
      if(compress == raw::kHuffman) UncompressHuffman(adc, uncompressed);
      else if(compress == raw::kZeroSuppression){
	ZeroUnsuppression(adc, uncompressed, pedestal);
      }
      else if(compress == raw::kZeroHuffman){
	gar::raw::ADCvector_t tmp(2*adc[0]);
	UncompressHuffman(adc, tmp);
	ZeroUnsuppression(tmp, uncompressed, pedestal);
      }
      else if(compress == raw::kNone){
	for(unsigned int i = 0; i < adc.size(); ++i) uncompressed[i] = adc[i];
      }
      else {
	throw cet::exception("raw")
	  << "raw::Uncompress() does not support compression #"
	  << ((int) compress);
      }
      return;
    }
  

    // the current Huffman Coding scheme used by uBooNE is
    // based on differences between adc values in adjacent time bins
    // the code is 
    // no change for 4 ticks --> 1
    // no change for 1 tick  --> 01
    // +1 change             --> 001
    // -1 change             --> 0001
    // +2 change             --> 00001
    // -2 change             --> 000001
    // +3 change             --> 0000001
    // -3 change             --> 00000001
    // abs(change) > 3       --> write actual value to short
    // use 15th bit to set whether a block is encoded or raw value
    // 1 --> Huffman coded, 0 --> raw
    // pad out the lowest bits in a word with 0's
    void CompressHuffman(gar::raw::ADCvector_t &adc)
    {
      gar::raw::ADCvector_t const orig_adc(std::move(adc));
    
      // diffs contains the difference between an element of adc and the previous
      // one; the first entry is never used.
      std::vector<short> diffs;
      diffs.reserve(orig_adc.size());
      std::adjacent_difference
	(orig_adc.begin(), orig_adc.end(), std::back_inserter(diffs));
    
      // prepare adc for the new data; we kind-of-expect the size,
      // so we pre-allocate it; we might want to shrink-to-fit at the end
      adc.clear();
      adc.reserve(orig_adc.size());
      // now loop over the diffs and do the Huffman encoding
      adc.push_back(orig_adc.front());
      unsigned int curb = 15U;

      std::bitset<16> bset;
      bset.set(15);

      for(size_t i = 1U; i < diffs.size(); ++i){

	switch (diffs[i]) {
	  // if the difference is 0, check to see what the next 3 differences are
        case 0 : {
	  if(i < diffs.size() - 3){
	    // if next 3 are also 0, set the next bit to be 1
	    if(diffs[i+1] == 0 && diffs[i+2] == 0 && diffs[i+3] == 0){
	      if(curb > 0){
		--curb;
		bset.set(curb);
		i += 3;
		continue;
	      }
	      else{	    
		adc.push_back(bset.to_ulong());
		
		// reset the bitset to be ready for the next word
		bset.reset();
		bset.set(15);
		bset.set(14); // account for the fact that this is a zero diff
		curb = 14;
		i += 3; 
		continue;
	      } // end if curb is not big enough to put current difference in bset	  
	    } // end if next 3 are also zero
	    else{
	      // 0 diff is encoded as 01, so move the current bit one to the right
	      if(curb > 1){
		curb -= 2;
		bset.set(curb);
		continue;
	      } // end if the current bit is large enough to set this one
	      else{	    
		adc.push_back(bset.to_ulong());
		// reset the bitset to be ready for the next word
		bset.reset();
		bset.set(15);
		bset.set(13); // account for the fact that this is a zero diff
		curb = 13;
		continue;
	      } // end if curb is not big enough to put current difference in bset	  	    
	    } // end if next 3 are not also 0
	  }// end if able to check next 3
	  else{
	    // 0 diff is encoded as 01, so move the current bit one to the right
	    if(curb > 1){
	      curb -= 2;
	      bset.set(curb);
	      continue;
	    } // end if the current bit is large enough to set this one
	    else{	    
	      adc.push_back(bset.to_ulong());
	      // reset the bitset to be ready for the next word
	      bset.reset();
	      bset.set(15);
	      bset.set(13); // account for the fact that this is a zero diff
	      curb = 13;
	      continue;
	    } // end if curb is not big enough to put current difference in bset	  
	  }// end if not able to check the next 3
	  break;
	}// end if current difference is zero
	case 1: {
	  if(curb > 2){
	    curb -= 3;
	    bset.set(curb);
	  }
	  else{
	    adc.push_back(bset.to_ulong());
	    // reset the bitset to be ready for the next word
	    bset.reset();
	    bset.set(15);
	    bset.set(12); // account for the fact that this is a +1 diff
	    curb = 12;
	  } // end if curb is not big enough to put current difference in bset	  
	  break;
	} // end if difference = 1
        case -1: {
	  if(curb > 3){
	    curb -= 4;
	    bset.set(curb);
	  }
	  else{
	    adc.push_back(bset.to_ulong());
	    // reset the bitset to be ready for the next word
	    bset.reset();
	    bset.set(15);
	    bset.set(11); // account for the fact that this is a -1 diff
	    curb = 11;
	  } // end if curb is not big enough to put current difference in bset	  
	  break;
	}// end if difference = -1
        case 2: {
	  if(curb > 4){
	    curb -= 5;
	    bset.set(curb);
	  }
	  else{
	    adc.push_back(bset.to_ulong());
	    // reset the bitset to be ready for the next word
	    bset.reset();
	    bset.set(15);
	    bset.set(10); // account for the fact that this is a +2 diff
	    curb = 10;
	  } // end if curb is not big enough to put current difference in bset	  
	  break;
	}// end if difference = 2
        case -2: {
	  if(curb > 5){
	    curb -= 6;
	    bset.set(curb);
	  }
	  else{
	    adc.push_back(bset.to_ulong());
	    // reset the bitset to be ready for the next word
	    bset.reset();
	    bset.set(15);
	    bset.set(9); // account for the fact that this is a -2 diff
	    curb = 9;
	  } // end if curb is not big enough to put current difference in bset	  
	  break;
	}// end if difference = -2
        case 3: {
	  if(curb > 6){
	    curb -= 7;
	    bset.set(curb);
	  }
	  else{
	    adc.push_back(bset.to_ulong());
	    // reset the bitset to be ready for the next word
	    bset.reset();
	    bset.set(15);
	    bset.set(8); // account for the fact that this is a +3 diff
	    curb = 8;
	  } // end if curb is not big enough to put current difference in bset	  
	  break;
	}// end if difference = 3
        case -3: {
	  if(curb > 7){
	    curb -= 8;
	    bset.set(curb);
	  }
	  else{
	    adc.push_back(bset.to_ulong());
	    // reset the bitset to be ready for the next word
	    bset.reset();
	    bset.set(15);
	    bset.set(7); // account for the fact that this is a -3 diff
	    curb = 7;
	  } // end if curb is not big enough to put current difference in bset	  
	  break;
	}// end if difference = -3
        default: {
	  // if the difference is too large that we have to put the entire adc value in:
	  // put the current value into the adc vec unless the current bit is 15, then there 
	  // were multiple large difference values in a row
	  if(curb != 15){
	    adc.push_back(bset.to_ulong());
	  }
          
	  bset.reset();
	  bset.set(15);
	  curb = 15;
          
	  // put the current adcvalue in adc, with its bit 15 set to 0
	  if(orig_adc[i] > 0) adc.push_back(orig_adc[i]);
	  else{
	    std::bitset<16> tbit(-orig_adc[i]);
	    tbit.set(14);
	    adc.push_back(tbit.to_ulong());
	  } 
	  break;
        } // if |difference| > 3
	}// switch diff[i]
      }// end loop over differences

      //write out the last bitset
      adc.push_back(bset.to_ulong());
    
      // this would reduce global memory usage,
      // at the cost of a new allocation and copy
      //  adc.shrink_to_fit();

    } // CompressHuffman()
    //--------------------------------------------------------
    // need to decrement the bit you are looking at to determine the deltas as that is how
    // the bits are set
    void UncompressHuffman(const gar::raw::ADCvector_t& adc, 
			   gar::raw::ADCvector_t      &uncompressed)
    {
    
      //the first entry in adc is a data value by construction
      uncompressed[0] = adc[0];

      unsigned int curu = 1;
      short curADC = uncompressed[0];

      // loop over the entries in adc and uncompress them according to the
      // encoding scheme above the CompressHuffman method
      for(unsigned int i = 1; i < adc.size() && curu < uncompressed.size(); ++i){

	std::bitset<16> bset(adc[i]);

	int numu = 0;

	//check the 15 bit to see if this entry is a full data value or not
	if( !bset.test(15) ){
	  curADC = adc[i];
	  if(bset.test(14)){
	    bset.set(14, false);
	    curADC = -1*bset.to_ulong();
	  }
	  uncompressed[curu] = curADC;

	  ++curu;
	}
	else{

	  int  b       = 14;
	  int  lowestb = 0;

	  // ignore any padding with zeros in the lower order bits
	  while( !bset.test(lowestb) && lowestb < 15) ++lowestb;

	  if(lowestb > 14){
	    mf::LogWarning("raw.cxx") << "encoded entry has no set bits!!! " 
				      << i << " "
				      << bset.to_string< char,std::char_traits<char>,std::allocator<char> >(); 
	    continue;
	  }

	  while( b >= lowestb){ 

	    // count the zeros between the current bit and the next on bit
	    int zerocnt = 0;
	    while( !bset.test(b-zerocnt) && b-zerocnt > lowestb) ++zerocnt;

	    b -= zerocnt;

	    if(zerocnt == 0){
	      for(int s = 0; s < 4; ++s){
		uncompressed[curu] = curADC;
		++curu;
		++numu;
		if(curu > uncompressed.size()-1) break;
	      }
	      --b;
	    }
	    else if(zerocnt == 1){
	      uncompressed[curu] = curADC;
	      ++curu;
	      ++numu;
	      --b;
	    }
	    else if(zerocnt == 2){
	      curADC += 1;
	      uncompressed[curu] = curADC;
	      ++curu;
	      ++numu;
	      --b;
	    }
	    else if(zerocnt == 3){
	      curADC -= 1;
	      uncompressed[curu] = curADC;
	      ++curu;
	      ++numu;
	      --b;
	    }
	    else if(zerocnt == 4){
	      curADC += 2;
	      uncompressed[curu] = curADC;
	      ++curu;
	      ++numu;
	      --b;
	    }
	    else if(zerocnt == 5){
	      curADC -= 2;
	      uncompressed[curu] = curADC;
	      ++curu;
	      ++numu;
	      --b;
	    }
	    else if(zerocnt == 6){
	      curADC += 3;
	      uncompressed[curu] = curADC;
	      ++curu;
	      ++numu;
	      --b;
	    }
	    else if(zerocnt == 7){
	      curADC -= 3;
	      uncompressed[curu] = curADC;
	      ++curu;
	      ++numu;
	      --b;
	    }

	    if(curu > uncompressed.size() - 1) break;

	  }// end loop over bits
 
	  if(curu > uncompressed.size() - 1) break;

	}// end if this entry in the vector is encoded

      }// end loop over entries in adc

      return;
    }

  } // namespace raw
} // namespace gar
