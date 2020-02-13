#include "CAF.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <string>

////////////////////////////////////////////////////////////////////////
//// Class:       cafanatree
//// File Name:   cafanatree_module.cc
////
//// Authors: Tanaz Mohayai and Eldwan Brianne
//// To run this module:
//// 1) cd Build
//// 2) rm -rf *
//// 3) cmake ${module home direcory}
//// 4) make install
//// 5) cd ${module home direcory}
//// 6) bin/cafanatree_module --infile ${name of the anatree file that is output from anatree module, not to be confused with edepsim file} --outfile ${a name of your choosing for the output file}
//////////////////////////////////////////////////////////////////////////

void ShowHelp()
{
    std::cout << "./cafanatree_module --infile <inputfile> --outfile <outputfile> --correct4origin <0/1>" << std::endl;
}

int main( int argc, char const *argv[] )
{
    if( argc == 1 || ((argc == 2) && ((std::string("--help") == argv[1]) || (std::string("-h") == argv[1]))) || argc < 7 ){
        ShowHelp();
        return 2;
    }

    if( argv[1] != std::string("--infile") || argv[3] != std::string("--outfile") || argv[5] != std::string("--correct4origin") ) {
        ShowHelp();
        return -2;
    }

    // get command line options
    std::string outfile = "";
    std::string infile = "";
    std::string correct4origin = "";
    int p = 0;
    while( p < argc )
    {
        if( argv[p] == std::string("--infile") ){
            infile = argv[p+1];
            p++;
        }
        else if( argv[p] == std::string("--outfile") ){
            outfile = argv[p+1];
            p++;
        }
        else if( argv[p] == std::string("--correct4origin") ){
            correct4origin = argv[p+1];
            p++;
        }
        else{
            p++;
        }
    }

    if(correct4origin != "0" && correct4origin != "1")
    {
        ShowHelp();
        return -2;
    }

    printf( "Making CAF from tree dump: %s\n", infile.c_str() );
    printf( "Output CAF file: %s\n", outfile.c_str() );
    printf( "Correct for Origin: %s\n", correct4origin.c_str() );

    CAF *caf = new CAF(infile, outfile, std::atoi(correct4origin.c_str()));
    if(not caf->BookTFile()) return -1;

    caf->loop();

    caf->WriteTTree();
    printf( "-30-\n" );

    return 0;
}
