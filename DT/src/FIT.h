#ifndef FIT_h
#define FIT_h

#include <iostream>
#include <vector>

#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "TDecompSVD.h"
#include <math.h>

// /////////////////////////////////////////////////////////////////////
//                                                                    //
// Silvia Pesente March 2009                                          //
//                                                                    //
// PATTERN RECOGNITION:                                               //
// Class to perform FIT                                               //
//                                                                    //
// /////////////////////////////////////////////////////////////////////


class FIT {
 public:
  
  ~FIT(); 
  void FIT_simple(int nrPnts, double *ax, double *ay, float &m, float &q, float &chi2);  
  void FIT_t0(double sigmaPhi,double sigmaTimes,int nrPnts,int nrMinPnts,int nrVarv, double *ax, double *ay, int *ac, std::pair<double,double> &m, std::pair<double,double> &q, std::pair<double,double> &t0, double &chi2, int &NPT, bool &okFit);
  void FIT_global(double *sigmaPhi,int *nrPnts,int *nrMinPnts,int *nrVarv, std::vector<std::vector<double> > axf, std::vector<std::vector<double> > at, std::vector<std::vector<double> > ay, std::vector<std::vector<int> > ac, double *m, double *sm, double *q, double *sq, double &t0, double &st0, double &vdrift, double *chi2, int *NPT, bool &okFit, int &fit_CH);
  
  
};
#endif /*FIT_h*/
