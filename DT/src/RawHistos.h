// 2008/01/28
// Sara Vanini : histograms for muon radiography data acquisition monitor
// SV 100118 put in class for TOM 

#ifndef RAWHISTOS_h
#define RAWHISTOS_h

#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <map>
#include <fstream>
#include <algorithm>

// ROOT utilities
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGCanvas.h"
#include "TPaveText.h"
#include "TPad.h"
#include "TFrame.h"
#include "TString.h"
#include "TFile.h"

// TOM classes
#include "HITCollection.h"


class RawHistos {

public:
  RawHistos();
  ~RawHistos();

  // member functions
  void initHistos();

  // build canvas
  void buildCanvasTDC(TCanvas * c=0);
  void buildCanvasCH1(TCanvas * c=0);
  void buildCanvasCH2(TCanvas * c=0);
  void buildCanvasSLS(TCanvas * c=0);
  void buildCanvasHITS(TCanvas * c=0);

  void updateCanvas(int flag);
  void fillHistos(HITCollection * hits);

  // flag for histos filling and RawHisto pointer
  inline void setFlagFillHistos(bool flag) { _flagFillHistos = flag; return; };
  inline bool flagFillHistos() { return _flagFillHistos; }
  inline bool flagCanvasTDC() { return _flagCanvasTDC; }
  inline bool flagCanvasCH1() { return _flagCanvasCH1; }
  inline bool flagCanvasCH2() { return _flagCanvasCH2; }
  inline bool flagCanvasSLS() { return _flagCanvasSLS; }
  inline bool flagCanvasHITS() { return _flagCanvasHITS; }

  inline void setCanvas(int flag) { _nCanvas = flag; return; };
  inline int getCanvas() { return _nCanvas; }
 
  void dumpHistos(char * file);


private:
  //for debugging
  bool DEBUG_FLAG;

  //utilities
  bool _flagFillHistos;
  bool _flagCanvasTDC, _flagCanvasCH1, _flagCanvasCH2, _flagCanvasSLS, _flagCanvasHITS;
  int _nCanvas;

  //occupancy histo cuts
  int minTime;
  int maxTime;

  //hit difference histo cuts
  int minTimeHit;
  int maxTimeHit;

  //define map array and function to fill it
  int nros;
  int nrob;
  int ntdc;
  int ncha;

  //graphic tools
  TCanvas 	*c1, *ch1, *ch2, *sls, *hits;
  TPad 		*pad1, *pad2, *pad3;
  TPad		*pad_hits_ch1, *pad_time_ch1;
  TPad		*pad_hits_ch2, *pad_time_ch2;
  TPad 		*pad_hits_sls, *pad_time_sls;
  TPad 		*pad_hits_hits;

  //histos
  TH1F 		* htime_sl[8];
  TH1F 		* hoccSL1_lay[4]; 
  TH1F 		* hoccSL2_lay[4]; 
  TH1F 		* hoccCH1_lay[12]; 
  TH1F 		* hoccCH2_lay[12]; 
  TH1F 		* hhit_diff[4];
  TH2F 		* hocc, * htimech;
  TH1F 		* htime;
};
#endif /*RAWHISTOS_h*/

