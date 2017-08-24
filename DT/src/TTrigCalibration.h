// /////////////////////////////////////////////////////////////////////
//                                                                    //
// Silvia Pesente February 2010                                       //
//                                                                    //
// From CMSSW (DTTTrigCalibration class):                             //
// Analyzer class which fills time box plots                          //
// with SL granularity for TTrig computation,                         //
// fits the rising edge of the time box with                          //
// the integral of a gaussian returning the                           //
// mean value and the sigma,                                          //
// and write results to a txt file.                                   //
// The time boxes are written to a root file.                         //
//                                                                    //
// /////////////////////////////////////////////////////////////////////


#ifndef TTrigCalibration_h
#define TTrigCalibration_h

#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <map>
#include <fstream>
#include <algorithm>

// ROOT utilities
#include "TROOT.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGCanvas.h"
#include "TPaveText.h"
#include "TPad.h"
#include "TFrame.h"
#include "TString.h"
#include "TFile.h"
#include "HITCollection.h"


class TTrigCalibration {

public:
  TTrigCalibration( int runTrig ); //ALTEA
  ~TTrigCalibration();

  // member functions
  void initHistos();
  void initVariables();

  // build canvas
  void buildCanvasTROB(TCanvas * c=0);
  void buildCanvasTSL(TCanvas * c=0);

  void computeTTrig(int flag);
  void fillHistos(HITCollection * hits);

  // flag for histos filling and RawHisto pointer
  inline void setFlagFillHistos(bool flag) { _flagFillHistos = flag; return; };
  // inline bool setFlagFillHistos(bool flag) { cout<<"flag "<<flag<<endl;_flagFillHistos = flag; return _flagFillHistos; };
  inline bool flagFillHistos() { return _flagFillHistos; }
  inline bool flagCanvasTROB() { return _flagCanvasTROB; }
  inline bool flagCanvasTSL() { return _flagCanvasTSL; }

  inline void setCanvas(int flag) { _nCanvas = flag; return; };
  inline int getCanvas() { return _nCanvas; }

  
  void dumpHistos(char * file);
  void saveTTrigFile(int runN);
  void dumpTTrigs(char * calibFileName);

  void fitTimeBox(TH1F *hTimeBox, double& TTRIG_ROB, double& RMS_ROB);
  // Automatically compute the seeds the range to be used for time box fit
  void getFitSeeds(TH1F *hTBox, double& mean, double& sigma, double& tBoxMax, double& xFitMin,double& xFitMax);
/*   // Define the integral of the gaussian to be used in the fit */
/*   double intGauss(double *x, double *par); */


private:
  //for debugging
  bool DEBUG_TTRIG;
  bool DEBUG_TTRIGBIN;

  //utilities
  bool _flagFillHistos;
  bool _flagCanvasTROB, _flagCanvasTSL;
  int _nCanvas;


  float TTrig_ROB_mean[6];
  float TTrig_ROB_RMS[6];
  float TTrig_mean[3];
  float TTrig_RMS[3]; 
  
  //time boxes (ns) for each ROB
  int hROBWidth;
  int hROBEdgeL;
  int hROBEdgeH;
  //time boxes (ns) for each SL
  int hSLWidth;
  int hSLEdgeL;
  int hSLEdgeH;

  //cut on raw time
  float rtime_min;
  float rtime_max;

  //graphic tools
  TCanvas 	*TBOX_ROB;
  TCanvas       *TBOX_SL;
  TPad 		*pad_time_ROB;
  TPad		*pad_time_SL;


  //histos
  TH1F          * htbox_ROB[8];
  TH1F          * htbox_SL[3];

  TFile *hDebugFile;

};

// Define the integral of the gaussian to be used in the fit
double intGauss(double *x, double *par);

#endif /*TTrigCalibration_h*/

