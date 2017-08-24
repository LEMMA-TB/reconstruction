#ifndef READERROS8_CUT_H
#define READERROS8_CUT_H

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Silvia Pesente January 2010
/// 20170609 Sara Vanini adapted to read data from ROS 8 LEMMA testbeam
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

#include "Geom.h"
#include "Save_HistosAndTree.h"
#include "Track_IO.h"
#include "TTrigCalibration.h"
// SV 100203-05 include TOMTOOL classes
#include "RawHistos.h"
#include "Occupancy.h" //ALTEA

class ReaderROS8_CUT
{
public:
    ReaderROS8_CUT();
    ~ReaderROS8_CUT();

    void setDebug(bool debug) { m_debug=debug; return;}

    int readEvent(FILE *infile, HITCollection *hits);
    void readROS(FILE * infile, unsigned int & wordCount, unsigned int &rosOffset, HITCollection *hits);
    void readTDCGroup(int finalROSWordCount, unsigned int tdcData, FILE * infile, unsigned int & wordCount, unsigned int rosChID, HITCollection *hits);
    void readPU(FILE * infile, unsigned int & wordCount, unsigned int &puOffset);
    unsigned int readWord(FILE *infile, unsigned int &wordCount);



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
 void goAnalysis(TString fin, int maxEvent=100, int runN=0, int runTrig=0, bool ttrigs=0, bool n2chambers=0, int chside=0);

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
 map<int,float>::iterator iter_t0;
 map<int,float>::iterator iter_ttrig;

 Geom *geo;
 TimeCorr *corr;
 Save_HistosAndTree *dump;
 Track_IO *inout;

 FILE *fo_txt;
 TFile *fo_Tree;
 TTree *tree;
 TFile *fo_Histo;

 int m_chside;


 // TOMTOOL classes
 RawHistos        * _rawHistos;
 TTrigCalibration * _ttrigCalib;

 // debug glag
 bool m_debug;
 bool m_integrity;

 //ALTEA
 Occupancy	  *_occupancy;
 ofstream fout;

};

#endif // READERROS8_CUT_H
