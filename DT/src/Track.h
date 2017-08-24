#ifndef Track_h
#define Track_h

#include <iostream>
#include <vector>
#include <math.h>
#include <utility>
using namespace std;


#include "Segment.h"


// /////////////////////////////////////////////////////////////////////
//                                                                    //
// Silvia Pesente April 2009                                          //
//                                                                    //
// PATTERN RECOGNITION:                                               //
// Class to compute TRACK PARAMETERS:                                 //
// SLOPE, X0, T0, CHI2, NPT                                           //
// from the fit_global of each chamber (PHI + THETA)                  //
//                                                                    //
// /////////////////////////////////////////////////////////////////////

class Track {
 public:

  // default constructor
  Track();

  // default destructor
  ~Track();

  void SelectTrack(HITCollection *hits, TimeCorr *corr, bool n2chambers=1);
  HITColl_Seg * SelectSeg(HITCollection *hits, int i, int CH, bool phi);
  float Get_MeanWire(HITColl_Seg *hit_seg);
  bool Track_IsGood();
  bool Set_IsGood(int s, bool isgood);
  bool Get_IsGood(int s);
  double Get_Slope(int s);
  double Set_Slope(int s,double slope);
  double Get_X0(int s);
  double Set_X0(int s,double x0);
  int Get_NPT(int s);
  int Set_NPT(int s,int npt);
  double Get_Chi2(int s);
  double Set_Chi2(int s,double chi2);
  double Get_T0(int s);
  double Set_T0(int s,double t0);
  double Get_T0_fin(int s);
  double Set_T0_fin(int s,double t0);
  double Get_Slope_glo(int s);
  double Set_Slope_glo(int s,double slope_glo);
  double Get_erSlope_glo(int s);
  double Set_erSlope_glo(int s,double erslope_glo);
  double Get_X0_glo(int s);
  double Set_X0_glo(int s,double x0_glo);
  double Get_erX0_glo(int s);
  double Set_erX0_glo(int s,double erx0_glo);
  int Get_NPT_glo(int s);
  int Set_NPT_glo(int s,int npt_glo);
  double Get_Chi2_glo(int s);
  double Set_Chi2_glo(int s,double chi2_glo);
  double Get_T0_glo(int s);
  double Set_T0_glo(int s,double t0_glo);
  double Get_erT0_glo(int s);
  double Set_erT0_glo(int s,double ert0_glo);
  double Get_T0_glo_fin(int s);
  double Set_T0_glo_fin(int s,double t0_glo_fin);
  bool Set_IsGood_glo(int s, bool isgood_glo);
  bool Get_IsGood_glo(int s);
  bool Get_MP_CH1_IsOk_glo();
  bool Get_MP_CH2_IsOk_glo();
  bool Get_SegIsAbsorbed(int s);
  void printTrack();
  
 private:
  FIT *fitg;
  HITColl_Seg *hit_seg[6];
  HITColl_Seg *hit_seg_glo[6];
  Segment *seg[6]; 
  bool isGood;
  bool track_IsGood[6];
  double track_Slope[6], track_X0[6], track_T0[6], track_Chi2[6], track_T0_fin[6];
  int track_NPT[6];
  bool track_IsGood_glo[4];
  double track_Slope_glo[4], track_X0_glo[4], track_T0_glo[4], track_T0_glo_fin[4], track_Chi2_glo[4];
  double track_erSlope_glo[4], track_erX0_glo[4], track_erT0_glo[4];
  int track_NPT_glo[4];
  bool MP_ch1_isok, MP_ch2_isok;
  bool SegIsAbsorbed[6];
  
};

#endif /*Track_h*/

