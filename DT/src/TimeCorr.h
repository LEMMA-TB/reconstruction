#ifndef TimeCorr_h
#define TimeCorr_h

#include <iostream>
#include <vector>
#include <math.h>
#include "TDirectory.h"
#include "TSpline.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "HITCollection.h"

// /////////////////////////////////////////////////////////////////////
//                                                                    //
// Silvia Pesente March 2009                                          //
//                                                                    //
// PATTERN RECOGNITION:                                               //
// Class to compute all TIME CORRECTIONS:                             //
//                                                                    //
// /////////////////////////////////////////////////////////////////////


class TimeCorr {
 public:
  // default constructor
  TimeCorr(){};
  
  // destructor
  ~TimeCorr(); 
  
  void InitSpline();
  /*   TSpline3 * GetSpline(int i); */
  float Get_Time_SlopeCorr(float time, float slope_raw);
  float Get_Time_DriftCorr(float time, int SL, int mean_wire);
  float Get_SigmaAngCorr(float slope_raw);
  float Get_LinearityCorr(int ch, int sl, float time, float slope_raw);
  
  void print();


 private:
  TFile *fsp;
  TSpline3 * spline[8];
  //TF1 * spline[8];
  
};
#endif /*TimeCorr_h*/
