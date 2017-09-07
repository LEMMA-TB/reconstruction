#ifndef Save_HistosAndTree_h
#define Save_HistosAndTree_h

/* system headers */
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <sys/times.h>
#include <sys/types.h>
#include <time.h> 
#include <iostream>
#include <iomanip>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TGCanvas.h"
#include "TPaveText.h"
#include "TPad.h"
#include "TFrame.h"

#include "Track.h"



// //////////////////////////////////////////////////////////////////////
//                                                                     //
// Silvia Pesente January 2010                                         //
//                                                                     //
// PATTERN RECOGNITION:                                                //
// Class for storing histograms and trees into root files              //
//                              			   	       //
// //////////////////////////////////////////////////////////////////////

class Save_HistosAndTree {
 public:
  // default constructor
  Save_HistosAndTree();
  ~Save_HistosAndTree();
  
  void set1chamber(){m_nmaxseg=2; m_nmaxseg_glo=0; return; }

  // operations
  void initHB(int runID, int numEvent, ofstream *HBFile);
  void dumpHB(Track *track, HITCollection *hits,int numEvent,ofstream *HBFile);

  void initTree();
  void bookTree(TTree*, bool n2chambers=1);
  void cleanTree();
  void init_Statistics();
  void dumpTree(Track *track, HITCollection *hits,int numEvent,TTree *tree);
  void dumpTree_Track(Track *track, HITCollection *hits,int numEvent,TTree *tree);
  void dumpTree_Hits(HITCollection *hits,int numEvent,TTree *tree);
  void dumpHisto(Track *track, HITCollection *hits,int numEvent);
  void initHistos();
  void writeHistos();
  void resetHistos();
  void deleteHistos();

  void fillVar_glo(int npc, int onseg, double m, double erm, double a, double era, double t0, double chi2, double * res, int p);
  void fillVar(int npc, int onseg, double m, double a, double t0, double chi2, double * res, int p, double *xh);

  void compute_Statistics(Track *track, HITCollection *hits,int numEvent);
  void write_Statistics(FILE *fo_txt);

  private:
  // variables for output tree
  int m_nmaxseg;
  int m_nmaxseg_glo;
  int onevent;
  int onseg, onseg_glo;
  int onhit;
  int   * ohlay, * ohwire, * ohtime_in, * ohtime, * ohtrig, * ohtime_tube, *oNhit_tube; 
  float * osegS, * osegX, * osegK, * osegT0;
  int   * osegN;
  float * osl1r, * osl2r, * osl3r, * osl4r, * osl5r, * osl6r, * osl7r, * osl8r;
  float * osxh1,* osxh2,* osxh3,* osxh4,* osxh5,* osxh6,* osxh7,* osxh8;
  float * osegS_glo, * osegX_glo, * osegK_glo, * osegT0_glo;
  float * osegerS_glo, * osegerX_glo;
  int   * osegN_glo;
  float * osl1r_glo, * osl2r_glo, * osl3r_glo, * osl4r_glo, * osl5r_glo, * osl6r_glo, * osl7r_glo, * osl8r_glo;
  

  /* struct definition */
  typedef struct {
    char bor[4]; // Begin of Run flag
    char type[5]; // Type of Run (MC - RAW)
    int number;   // Run number
    char time[30]; // Run starting time
    int  nChamber;   // Number of chambers
  } runHEADER;
  
  typedef struct {
    char type[5]; // Type of Event
    int number;   // Event number
    int nHits[4];    // Number of hits (per chamber)
  } evtHEADER;

  typedef struct {
    int nTracks;  // number of tracks, 0 if any
  } trkHEADER;

  typedef struct {
    int chamber;   // Chamber
    int layer;      //Layer number
    int XorZ;      //Layer type (0=X, 2=Z)
    double coord;    // Measured coordinate (X if XorZ=0, Z if XorZ=2)
    double dist;    // distance from wire (X if XorZ=0, Z if XorZ=2)
    double coordMC;    // MC coordinate
    double quota;   // Hit y coordinate (same as layer)
  } HIT_HB;

  typedef struct {
    int chamber;   // Chamber
    int nPoints;   // Number of points of the fit
    int XorZ;      // Layer type (0=X (phi), 2=Z (theta))
    double slope;   // Track slope (tan(angle)) - from the fit
    double inter;   // Intercept - from the fit
    double erSlope; // Error on the slope
    double erInter; // Error on the intercept
    double erCorr;  // Correlation error
    double tZero;   // t0 (in MC always = -1000);
    double chi2;    // Fit chi2
    int   lay[8];  // Layer number for each point of the track (8 maximum)
    double res[8];  // Residuals for each point of the track (8 maximum)
  } TRACK_HB;


  //histos
  static const int ev_graph=10;
  int count_graph;
  TGraph        *g1_hit[ev_graph], *g1_hit_in[ev_graph];
  TGraph        *g2_hit[ev_graph], *g2_hit_in[ev_graph];
  TGraph        *g3_hit[ev_graph], *g3_hit_in[ev_graph];
  TGraph        *g4_hit[ev_graph], *g4_hit_in[ev_graph];
  TGraph        *g5_hit[ev_graph];
  TGraph        *g6_hit[ev_graph];
  TH1F          *h_Nhit;
  TH1F          *h_tempo[6];
  TH1F          *h_Nhit_xCH[4];
  TH1F          *h_Nhit_xSL[8];
  TH1F          *h_Nhit_xL[4];
  TH1F          *h_Nhit_xSL13;
  TH1F          *h_SLOPE[6];
  TH1F          *h_dphi;
  TH1F          *h_dphi_glo;
  TH1F          *h_dthe;
  TH1F          *h_dthe_glo;
  TH1F          *h_dT0;
  TH1F          *h_dT0_glo;
  TH2F          *h_dT0_slope_glo[4];
  TH2F          *h_dT0_NPT_glo[4];
  TH1F          *h_dT0_fin;
  TH1F          *h_dT0_fin2;
  TH2F          *h_T0ch1_T0ch2;
  TH2F          *h_T0ch1_T0ch2_glo;
  TH1F          *h_X0[6];
  TH1F          *h_T0[6];
  TH1F          *h_CHI2[6];
  TH1F          *h_NPT[6];
  TH1F          *h_SLOPE_glo[4];
  TH1F          *h_X0_glo[4];
  TH1F          *h_T0_glo[2];
  TH1F          *h_erSLOPE_glo[4];
  TH1F          *h_erX0_glo[4];
  TH1F          *h_erT0_glo[2];
  TH1F          *h_MP_CH_glo[2];
  TH1F          *h_resCH1Phi_glo;
  TH1F          *h_resCH1The_glo;
  TH1F          *h_resCH2Phi_glo;
  TH1F          *h_resCH2The_glo;
  TH1F          *h_resSLup;
  TH1F          *h_resSLdown;
  TH1F          *h_resCH1phi_glo[8];
  TH1F          *h_resCH2phi_glo[8];
  TH1F          *h_resCH1the_glo[4];
  TH1F          *h_resCH2the_glo[4];
  TH1F          *h_resSL1[4];
  TH1F          *h_resSL2[4];
  TH2F          *h_lincorr[6][8];
  TH2F          *h_lincorr8[8];

  // variables for statistics_file
  int ev_analyzed;
  int ev_CH1_fitOk, ev_CH2_fitOk, ev_SL1_fitOk, ev_SL2_fitOk;
  int ev_CH1_fitgloOk, ev_CH2_fitgloOk;
  int ev_CH1_fitgloOk_NPT8, ev_CH1_fitgloOk_NPT7, ev_CH1_fitgloOk_NPT6;
  int ev_CH2_fitgloOk_NPT8, ev_CH2_fitgloOk_NPT7, ev_CH2_fitgloOk_NPT6;
  int ev_CH1_fitgloOk_MP, ev_CH2_fitgloOk_MP;
  int ev_CH1_fitgloOk_MP_NPT8, ev_CH1_fitgloOk_MP_NPT7, ev_CH1_fitgloOk_MP_NPT6;
  int ev_CH2_fitgloOk_MP_NPT8, ev_CH2_fitgloOk_MP_NPT7, ev_CH2_fitgloOk_MP_NPT6;
  float res_CH1_phi_mean, res_CH1_the_mean, res_CH1_phi_sigma, res_CH1_the_sigma; 
  float res_CH2_phi_mean, res_CH2_the_mean, res_CH2_phi_sigma, res_CH2_the_sigma; 
  float res_SL1_mean, res_SL2_mean, res_SL1_sigma, res_SL2_sigma; 


};
#endif /*Save_HistosAndTree_h*/

