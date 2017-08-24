#ifndef Segment_h
#define Segment_h

#include <iostream>
#include <vector>
#include <math.h>
#include <utility>
using namespace std;


#include "HITColl_Layer.h"
#include "Geom.h"
#include "FIT.h"
#include "TimeCorr.h"
#include "TMath.h"


// /////////////////////////////////////////////////////////////////////
//                                                                    //
// Silvia Pesente March 2009                                          //
//                                                                    //
// PATTERN RECOGNITION:                                               //
// Class to SOLVE LEFT/RIGHT AMBIGUITY of each segment                //
// and compute SEGMENT PARAMETERS:                                    //
// SLOPE, X0, T0, CHI2, NPT                                           //
// from the fit_t0 of each segment (PHI or THETA)                     //
//                                                                    //
// /////////////////////////////////////////////////////////////////////


class Segment {
 public:

  // default constructor
  Segment();

  // default destructor
  ~Segment();

  HITColl_Seg * Get_PHISegment(HITColl_Seg *hit_seg, float w_xcorr, TimeCorr *corr); 
  HITColl_Seg * Get_SLSegment(HITColl_Seg *hit_seg, float w_xcorr, TimeCorr *corr, float T0phi); 
  HITColl_Seg * RejectDoubleHits(HITColl_Seg *hit_seg); 
  vector <pair<float,float> > Xmean_SL(HITColl_Seg *hit_seg, int SL,int & size ); 
  pair<float,float> Get_Xmean_SL(HITColl_Seg *hit_seg,int CH,int SL,float mean); 
  pair<float,float> Get_SlopePHI(pair<float,float> PHI1, pair<float,float> PHI2); 
  void Set_HITTime_Drift(HITColl_Seg *hit_seg, float w_xcorr, TimeCorr *corr);
  void Set_HITTime_Slope(HITColl_Seg *hit_seg, pair<float,float> seg_par, TimeCorr *corr);
  void Set_HITTime_Linearity(HITColl_Seg *hit_seg, pair<float,float> seg_par, TimeCorr *corr);
  void Set_HITTime_T0(HITColl_Seg *hit_seg, float T0);
  HITColl_Seg * SelectHITsforFIT(HITColl_Seg *hit_seg, pair<float,float> seg_par); 
  HITColl_Seg * SelectHITsforSLFIT(HITColl_Seg *hit_seg, float xmean); 
  int Get_NSolvedHIT(HITColl_Seg *hit_seg);
  int Get_NnotSolvedHIT(HITColl_Seg *hit_seg);
  bool HITisDouble(HITColl_Seg *hit_seg,int i);
  void Set_DoubleHITs(HITColl_Seg *hit_seg);
  void SolveAmbiguity(double sigmaSlope,double sigmaTime,HITColl_Seg *hit_seg,int nrVarv,bool change_code);  
  void Set_HITCode(int *cc, int *hh);
  HITColl_Seg * RejectHits0(HITColl_Seg *hit_seg);
  int * Get_DoubleHITxL(HITColl_Seg *hit_seg,int nL);
  int * Get_NX_xLs(HITColl_Seg *hit_seg,int nL);
  double Set_Slope(pair<double,double> slope);
  double Get_Slope();
  double Set_X0(pair<double,double> X0);
  double Get_X0();
  double Set_T0(pair<double,double> T0);
  double Get_T0();
  double Set_T0_fin(pair<double,double> T0_fin);
  double Get_T0_fin();
  int   Set_NPT(int npt);
  int   Get_NPT();
  bool  Set_FitOK(bool okFit);
  bool  Get_FitOK();
  double Set_Chi2(double chi2);
  double Get_Chi2();
  bool Set_IsDouble(bool isdouble);
  bool Set_IsGood(bool isgood);
  bool IsDouble();
  bool IsGood();
  
  void printSegment();
  
  
 private:
  HITColl_Seg *hit_seg;
  FIT *fit;
  HITColl_Layer *layer[8]; 
  int *L_2hit;
  int *NX_xLs;
  int ii, ii1, ii2; 
  double seg_Slope, seg_X0, seg_T0, seg_T0_fin, seg_chi2;
  int seg_NPT;
  bool onlyoneseg;
  bool seg_isdouble, seg_isgood, seg_okFit;  

};

#endif /*Segment_h*/

