#include "TFile.h"
#include "TTree.h"

#include <vector>
#include <iostream>

using namespace std;

void seg.C( string filename1, string filename2, int runNum ) {
    
    TFile *file1 = new TFile( filename1 );
    TFile *file2 = new TFile( filename2 );
    
    TTree *t1 = ( TTree* )file1->Get( "RADMU" );
    TTree *t2 = ( TTree* )file2->Get( "RADMU" );
    
    // set branch address tree1
    
    int event1, Nhits1, trtime1, rawT_in1[120], rawT1[120], rawT_tube1[120], Nhit_tube1[120];
    int lay1[120], tube1[120];
    t1->SetBranchAddress( "EVENT", &event1 );
    t1->SetBranchAddress( "DTBX_nhit", &Nhits1);
    t1->SetBranchAddress( "TRTIME", &trtime1 );
    t1->SetBranchAddress( "DTBX_lay", lay1 );
    t1->SetBranchAddress( "DTBX_tube", tube1 );
    t1->SetBranchAddress( "DTBX_time_in", rawT_in1);
    t1->SetBranchAddress( "DTBX_time", rawT1);
    t1->SetBranchAddress( "DTBX_time_tube", rawT_tube1 );
    t1->SetBranchAddress( "DTBX_Nhit_tube", Nhit_tube1 );
    
    Int_t iseg1;
    Int_t segN1[50];
    Float_t segX1[50], segS1[50], segK1[50], segT1[50];
    Float_t r11[50], r21[50], r31[50], r41[50], r51[50], r61[50], r71[50], r81[50];
    Float_t xh11[50], xh21[50], xh31[50], xh41[50], xh51[50], xh61[50], xh71[50], xh81[50];
    
    t1->SetBranchAddress( "SEG_ns", &iseg1 );
    t1->SetBranchAddress( "SEG_sx",  segX1 );
    t1->SetBranchAddress( "SEG_ss",  segS1 );
    t1->SetBranchAddress( "SEG_sk",  segK1 );
    t1->SetBranchAddress( "SEG_sn",  segN1 );
    t1->SetBranchAddress( "SEG_t0",  segT1 );

    t1->SetBranchAddress( "SEG_1r",  r11 );
    t1->SetBranchAddress( "SEG_2r",  r21 );
    t1->SetBranchAddress( "SEG_3r",  r31 );
    t1->SetBranchAddress( "SEG_4r",  r41 );
    t1->SetBranchAddress( "SEG_5r",  r51 );
    t1->SetBranchAddress( "SEG_6r",  r61 );
    t1->SetBranchAddress( "SEG_7r",  r71 );
    t1->SetBranchAddress( "SEG_8r",  r81 );
    
    t1->SetBranchAddress( "SEG_xh1", xh11 );
    t1->SetBranchAddress( "SEG_xh2", xh21 );
    t1->SetBranchAddress( "SEG_xh3", xh31 );
    t1->SetBranchAddress( "SEG_xh4", xh41 );
    t1->SetBranchAddress( "SEG_xh5", xh51 );
    t1->SetBranchAddress( "SEG_xh6", xh61 );
    t1->SetBranchAddress( "SEG_xh7", xh71 );
    t1->SetBranchAddress( "SEG_xh8", xh81 );
    
    //set branch address tree2
    
    int event2, Nhits2, trtime2, rawT_in2[120], rawT2[120], rawT_tube2[120], Nhit_tube2[120];
    int lay2[120], tube2[120];
    t2->SetBranchAddress( "EVENT", &event2 );
    t2->SetBranchAddress( "DTBX_nhit", &Nhits2);
    t2->SetBranchAddress( "TRTIME", &trtime2 );
    t2->SetBranchAddress( "DTBX_lay", lay2 );
    t2->SetBranchAddress( "DTBX_tube", tube2 );
    t2->SetBranchAddress( "DTBX_time_in", rawT_in2 );
    t2->SetBranchAddress( "DTBX_time", rawT2 );
    t2->SetBranchAddress( "DTBX_time_tube", rawT_tube2 );
    t2->SetBranchAddress( "DTBX_Nhit_tube", Nhit_tube2 );

    Int_t iseg2;
    Int_t segN2[50];
    Float_t segX2[50], segS2[50], segK2[50], segT2[50];
    Float_t r12[50], r22[50], r32[50], r42[50], r52[50], r62[50], r72[50], r82[50];
    Float_t xh12[50], xh22[50], xh32[50], xh42[50], xh52[50], xh62[50], xh72[50], xh82[50];

    t2->SetBranchAddress( "SEG_ns", &iseg2 );
    t2->SetBranchAddress( "SEG_sx",  segX2 );
    t2->SetBranchAddress( "SEG_ss",  segS2 );
    t2->SetBranchAddress( "SEG_sk",  segK2 );
    t2->SetBranchAddress( "SEG_sn",  segN2 );
    t2->SetBranchAddress( "SEG_t0",  segT2 );

    t2->SetBranchAddress( "SEG_1r",  r12 );
    t2->SetBranchAddress( "SEG_2r",  r22 );
    t2->SetBranchAddress( "SEG_3r",  r32 );
    t2->SetBranchAddress( "SEG_4r",  r42 );
    t2->SetBranchAddress( "SEG_5r",  r52 );
    t2->SetBranchAddress( "SEG_6r",  r62 );
    t2->SetBranchAddress( "SEG_7r",  r72 );
    t2->SetBranchAddress( "SEG_8r",  r82 );
    
    t2->SetBranchAddress( "SEG_xh1", xh12 );
    t2->SetBranchAddress( "SEG_xh2", xh22 );
    t2->SetBranchAddress( "SEG_xh3", xh32 );
    t2->SetBranchAddress( "SEG_xh4", xh42 );
    t2->SetBranchAddress( "SEG_xh5", xh52 );
    t2->SetBranchAddress( "SEG_xh6", xh62 );
    t2->SetBranchAddress( "SEG_xh7", xh72 );
    t2->SetBranchAddress( "SEG_xh8", xh82 );
    
    // set variables for new tree
    int m_nmaxseg = 6;
    int nmaxhit = 100;
    int onevent = 0;
    int onseg = m_nmaxseg;
    int onhit = ;
    int * ohlay = new int[nmaxhit];
    int * ohwire = new int[nmaxhit];
    int * ohtime_in = new int[nmaxhit];
    int * ohtime = new int[nmaxhit];
    int * ohtrig = new int[nmaxhit];
    int * ohtime_tube = new int[108];
    int * oNhit_tube = new int[18]; 
    float * osegS = new float[m_nmaxseg];
    float * osegX = new float[m_nmaxseg];
    float * osegK = new float[m_nmaxseg];
    float * osegT0 = new float[m_nmaxseg];
    int   * osegN = new float[m_nmaxseg];
    float * osl1r = new float[m_nmaxseg];
    float * osl2r = new float[m_nmaxseg];
    float * osl3r = new float[m_nmaxseg];
    float * osl4r = new float[m_nmaxseg];
    float * osl5r = new float[m_nmaxseg];
    float * osl6r = new float[m_nmaxseg];
    float * osl7r = new float[m_nmaxseg];
    float * osl8r = new float[m_nmaxseg];
    float * osxh1 = new float[m_nmaxseg];
    float * osxh2 = new float[m_nmaxseg];
    float * osxh3 = new float[m_nmaxseg];
    float * osxh4 = new float[m_nmaxseg];
    float * osxh5 = new float[m_nmaxseg];
    float * osxh6 = new float[m_nmaxseg];
    float * osxh7 = new float[m_nmaxseg];
    float * osxh8 = new float[m_nmaxseg];
    
    // new tree: tree3
    
    char filename3[200];
    sprintf( filename3, "./output/Run_%d_merged.root", runNum );
    
    TFile *file3 = new TFile( filename3 );
    TTree *t3 = new TTree( "RADMU" );
    
    t3->Branch( "EVENT", &onevent,  "EVENT/I" );
    t3->Branch( "DTBX_nhit", &onhit, "Nhits/I" );
    t3->Branch( "TRTIME",     ohtrig,"trigT[Nhits]/I" );
    t3->Branch( "DTBX_lay",   ohlay, "hlay[Nhits]/I" );
    t3->Branch( "DTBX_tube",  ohwire,"htube[Nhits]/I" );
    t3->Branch( "DTBX_time_in",ohtime_in,"htime_in[Nhits]/I" );
    t3->Branch( "DTBX_time",  ohtime,"htime[Nhits]/I" );
    t3->Branch( "DTBX_time_tube",ohtime_tube,"rtime[108]/I" );
    t3->Branch( "DTBX_Nhit_tube",oNhit_tube,"Nhit[18]/I" );
    
    t3->Branch( "SEG_ns", &onseg,  "ntes/I" );
    t3->Branch( "SEG_sx",  osegX,  "X[ntes]/F" );
    t3->Branch( "SEG_ss",  osegS,  "SLOPE[ntes]/F" );
    t3->Branch( "SEG_sk",  osegK,  "CHI2[ntes]/F" );
    t3->Branch( "SEG_sn",  osegN,  "NPT[ntes]/I" );
    t3->Branch( "SEG_t0",  osegT0, "T0[ntes]/F" );
    t3->Branch( "SEG_1r",  osl1r,   "l1r[ntes]/F" );
    t3->Branch( "SEG_2r",  osl2r,   "l2r[ntes]/F" );
    t3->Branch( "SEG_3r",  osl3r,   "l3r[ntes]/F" );
    t3->Branch( "SEG_4r",  osl4r,   "l4r[ntes]/F" );
    t3->Branch( "SEG_5r",  osl5r,   "l5r[ntes]/F" );
    t3->Branch( "SEG_6r",  osl6r,   "l6r[ntes]/F" );
    t3->Branch( "SEG_7r",  osl7r,   "l7r[ntes]/F" );
    t3->Branch( "SEG_8r",  osl8r,   "l8r[ntes]/F" );

    t3->Branch( "SEG_xh1",  osxh1,   "xh1[ntes]/F" );
    t3->Branch( "SEG_xh2",  osxh2,   "xh2[ntes]/F" );
    t3->Branch( "SEG_xh3",  osxh3,   "xh3[ntes]/F" );
    t3->Branch( "SEG_xh4",  osxh4,   "xh4[ntes]/F" );
    t3->Branch( "SEG_xh5",  osxh5,   "xh5[ntes]/F" );
    t3->Branch( "SEG_xh6",  osxh6,   "xh6[ntes]/F" );
    t3->Branch( "SEG_xh7",  osxh7,   "xh7[ntes]/F" );
    t3->Branch( "SEG_xh8",  osxh8,   "xh8[ntes]/F" );
    
    
    
}
