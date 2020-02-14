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
    std::cout << "./cafanatree_module --infile <inputfile> --outfile <outputfile> --correct4origin <0/1> --originTPC <x> <y> <z> (in cm)" << std::endl;
}

int main( int argc, char const *argv[] )
{
    if( argc == 1 || ((argc == 2) && ((std::string("--help") == argv[1]) || (std::string("-h") == argv[1]))) || argc < 8 ){
        ShowHelp();
        return 2;
    }

    if( argv[1] != std::string("--infile") || argv[3] != std::string("--outfile") || argv[5] != std::string("--correct4origin") || argv[7] != std::string("--originTPC") ) {
        ShowHelp();
        return -2;
    }

    // get command line options
    std::string outfile = "";
    std::string infile = "";
    std::string correct4origin = "";
    std::string x, y, z = "";
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
        else if( argv[p] == std::string("--originTPC") ){
            if( argc > 8 && argc == 11 )
            {
                x = argv[p+1];
                y = argv[p+2];
                z = argv[p+3];
                p += 3;
            }
            else{
                std::cout << "Missing an origin coordinate!" << std::endl;
                ShowHelp();
                return -2;
            }
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

    if( x == "" && y == "" && z == "" )
    {
        printf("No TPC offset given, defaulting to (0, 0, 0)!!\n");
        x = y = z = "0";
    }

    printf( "Making CAF from tree dump: %s\n", infile.c_str() );
    printf( "Output CAF file: %s\n", outfile.c_str() );
    printf( "Correct for Origin: %s\n", correct4origin.c_str() );
    printf( "TPC offset: (%s, %s, %s) cm\n", x.c_str(), y.c_str(), z.c_str() );

    double originTPC[3] = {std::atof(x.c_str()), std::atof(y.c_str()), std::atof(z.c_str())};

    CAF *caf = new CAF(infile, outfile, std::atoi(correct4origin.c_str()), &originTPC[0]);
    if(not caf->BookTFile()) return -1;

    caf->loop();

    caf->WriteTTree();
    printf( "-30-\n" );

    return 0;
}
