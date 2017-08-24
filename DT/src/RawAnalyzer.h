// //////////////////////////////////////////////////////////////////////
//                                                                     //
// Silvia Pesente January 2010                                         //
//                                                                     //
// PATTERN RECOGNITION:                                                //
// Read rawfile,                                                       //
// Do Pattern Recognition,                                             // 
// Analyze track                                                       //
// Create and save tree and store histos into files                    //
//                                                                     //
// 100202 : SV TOMIO merged here                                       //
// //////////////////////////////////////////////////////////////////////


#include <iostream>
#include <algorithm>

#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <map>
#include <math.h>
#include <stdlib.h>
#include <string.h>
using std::string;  
//using namespace std;

#include "TFile.h"
/* #include "TROOT.h" */

#include "Geom.h"
#include "Save_HistosAndTree.h"
#include "Track_IO.h"
#include "TTrigCalibration.h"

// SV 100203-05 include TOMTOOL classes
#include "RawHistos.h"

// SV 100210 include EM classes
/* #include "../EM/ImgAnalyzer.h" */

class RawAnalyzer {
  
 public:
  RawAnalyzer();
  ~RawAnalyzer();
  
  // set and get RawHistos pointer
  inline void setRawHistosPtr(RawHistos * _ptr) { _rawHistos = _ptr; return; };
  inline bool RawHistosPtr() { return _rawHistos; }

  // set and get TTrigCalibration pointer
  inline void setTTrigCalibPtr(TTrigCalibration * _ptr) { _ttrigCalib = _ptr; return; };
  inline bool TTrigCalibPtr() { return _ttrigCalib; }

  /*  // set and get ImgAnalyzer pointer
  inline void setImgAnalyzerPtr(ImgAnalyzer * _ptr) { _imgAnalyzer = _ptr; return; };
  inline bool ImgAnalyzerPtr() { return _imgAnalyzer; }
  */

  // member functions
  void goAnalysis(TString fin, int maxEvent=100, int runN=0, int runTrig=0);
  void goAnalysis(TString fin, int maxEvent=100, int runN=0, int runTrig=0, bool ttrig=0);
  void readWord(bool &break_flag, FILE *infile);
  void collectHit();
  void fillMap();
  void checkMap();
  void fillMap_t0();
  void checkMap_t0();
  void fillMap_ttrig(int runTrig);
  void checkMap_ttrig(int runTrig);
  int getTDCid( int ros, int rob, int tdc, int cha );
  int getTubeId( int se, int sl, int lay, int tube );
  int getTube(int idTube);
  int getLay(int idTube);
  int getSL(int idTube);
  int getSe(int idTube);  
  int getSLId( int se, int sl);

 private:
  map<int, int> chmap;
  map<int, float> t0map;
  map<int, float> ttrigmap;
  map<int,int>::iterator iter;
  map<int,float>::iterator iter_t0;
  map<int,float>::iterator iter_ttrig;
  
  Track *track;
  HITCollection *hits;
  Geom *geo;
  TimeCorr *corr;
  Save_HistosAndTree *dump;
  Track_IO *inout;
  
  FILE *rawfile;
  FILE *fo_txt;
  TFile *fo_Tree;
  TTree *tree;
  TFile *fo_Histo;
  int maxwords;
  ofstream *HBFile;
  
  float dtime;
  //variable declaration: the value is filled at the first occurance
  long Event_Id, Bunch_Id, TDC_Id, channel, rawtime, ROB_Id;
  //long int is 32-bit=4-byte word
  long evbu, type, TTCcount;
  
  int numEvent;
  int numEventDAQ;
  bool check_header;
  int nHEADER,ngheader,nTRAILER,ngtrailer,ngroup,nERROR,nDEBUG;

  
  // TOMTOOL classes
  RawHistos        * _rawHistos;
  TTrigCalibration * _ttrigCalib;
  /*   ImgAnalyzer      * _imgAnalyzer;   */
  
  //OBJECTS USED TO ORGANIZE AND PRESENT DATA
  TCanvas               *fC1;
  
};

