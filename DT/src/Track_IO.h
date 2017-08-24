#ifndef Track_IO_h
#define Track_IO_h

/* system headers */
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <sys/times.h>
#include <sys/types.h>
#include <time.h> 
#include <iostream>
#include <iomanip>
#include <fstream>
#include "TFile.h"
#include "TString.h"

using namespace std;

// /////////////////////////////////////////////////////////////////////
//                                                                    //
// Silvia Pesente Februray 2010                                       //
//                                                                    //
// PATTERN RECOGNITION:                                               //
// Class for all IO operations for PATTERN RECOGNITION                //
//                                                                    //
// /////////////////////////////////////////////////////////////////////

class Track_IO {
 public:
  // default constructor
  Track_IO();
  ~Track_IO();
  
  
  // operations
	// raw file
  //  void readRawFile(char* rawFileName);
  FILE * openINFile(TString fin, int ID);
  FILE * openTXTFile(int runN, int maxEvent);
  ofstream * openOUTHBFile(int runN, int maxEvent);
  TFile * openOUTRootFile(int runN, int maxEvent, int chside=0);
  TFile * openOUTHistoFile(int runN, int maxEvent);

  private:
  FILE *infile;
  FILE *txtfile;
  TFile *fo_tree;
  TFile *fo_histo;
  ofstream *HBfile;
 
};
#endif /*Track_IO_h*/

