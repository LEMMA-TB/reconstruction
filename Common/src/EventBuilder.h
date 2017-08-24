#ifndef EVENTBUILDER_H
#define EVENTBUILDER_H


// NEW VERSION MODIFIED BY A. LORENZON 23.08.2017

//---- STL
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

#include "TFile.h"
#include "TTree.h"

using namespace std;

class EventBuilder
{
public:
    
    EventBuilder();
    ~EventBuilder(){}
    
    void setDebug(bool debug);
    void openAlignments(std::string alignment_file);
    void setAlignment(bool doAli);
    void calibrate();
    void openDataFiles(std::string inputDTFile, std::string inputSiFile, std::string outputFile);
    void matchEvents(int nEvents);
    int  getDtEvent();
    int  getSiEvent();
    void dumpGlobalEvent();
    void dumpOutput();
    
private:
    
    // --- global variables
    static const Int_t nMaxHits           = 100;
    static const Int_t nSiLayers          = 8; 
    static const Int_t nMaxClusters       = 5;
    static const Int_t nDevaChannels      = 6;
    static const Int_t nCherenkovChannels = 4;
   
    // --- si data file
    TFile *m_siFile;
    TTree *t2;
    
    // --- Declaration of input si tree leaf types
    Int_t           si_iev;
    Int_t           si_nhits;
    Int_t           si_ievmu;
    Int_t           si_subdet[nMaxHits];
    Float_t         si_xh[nMaxHits];
    Float_t         si_xh_PulseHeight[nMaxHits];
    Float_t         si_yh[nMaxHits];
    Float_t         si_yh_PulseHeight[nMaxHits];
    Float_t         si_zh[nMaxHits];
    Int_t           si_itrack[nMaxHits];
    Int_t           si_S2_PulseHeight;
    Int_t           si_S2_TimeSpectrum;
    Int_t           si_LemmaTrigger_PulseHeight;
    Int_t           si_LemmaTrigger_TimeSpectrum;
    Int_t           si_Deva_PulseHeight[nDevaChannels];
    Int_t           si_Deva_TimeSpectrum[nDevaChannels];
    Int_t           si_Cherenkov_PulseHeight[nCherenkovChannels];
    Int_t           si_Cherenkov_TimeSpectrum[nCherenkovChannels];
    Int_t           si_S3_PulseHeight;
    Int_t           si_S3_TimeSpectrum;
    Int_t           si_PbGlass_PulseHeight;
    Int_t           si_PbGlass_TimeSpectrum;
    
    Int_t           m_isi;
    
    // --- Declaration of si tree branches
    TBranch        *b_iev;   
    TBranch        *b_ievmu; 
    TBranch        *b_nhits; 
    TBranch        *b_subdet;
    TBranch        *b_xh;
    TBranch        *b_xh_PulseHeight;
    TBranch        *b_yh;
    TBranch        *b_yh_PulseHeight;
    TBranch        *b_zh;   
    TBranch        *b_itrack;
    TBranch        *b_S2_PulseHeight;
    TBranch        *b_S2_TimeSpectrum;
    TBranch        *b_LemmaTrigger_PulseHeight;
    TBranch        *b_LemmaTrigger_TimeSpectrum;
    TBranch        *b_Deva_PulseHeight;
    TBranch        *b_Deva_TimeSpectrum;
    TBranch        *b_Cherenkov_PulseHeight;
    TBranch        *b_Cherenkov_TimeSpectrum;
    TBranch        *b_S3_PulseHeight;
    TBranch        *b_S3_TimeSpectrum;
    TBranch        *b_PbGlass_PulseHeight;
    TBranch        *b_PbGlass_TimeSpectrum;
    
    // --- dt data file
    TFile * m_dtFile;
    TTree * m_tree;
    
    // --- Declaration of input mu tree leaf type
    Int_t           dtEvent;
    Int_t           nseg;
    Int_t           segN[2];
    Float_t         s1p[2];   //[ntes]
    Float_t         s2p[2];   //[ntes]
    Float_t         s3p[2];   //[ntes]
    Float_t         s4p[2];   //[ntes]
    Float_t         s5p[2];   //[ntes]
    Float_t         s6p[2];   //[ntes]
    Float_t         s7p[2];   //[ntes]
    Float_t         s8p[2];   //[ntes]
    
    Int_t m_idt;
    
    // --- Declaration of mu tree branches
    TBranch        *b_EVENT;   
    TBranch        *b_ntes;   
    TBranch        *b_SEG_sn; 
    TBranch        *b_SEG_xh1;
    TBranch        *b_SEG_xh2;
    TBranch        *b_SEG_xh3;
    TBranch        *b_SEG_xh4;
    TBranch        *b_SEG_xh5;
    TBranch        *b_SEG_xh6;
    TBranch        *b_SEG_xh7;
    TBranch        *b_SEG_xh8;
    
    // --- output file
    TFile *         m_outFile;
    TTree *         m_outTree;
    
    // --- Declaration of output leaf types
    Int_t           iev;
    Int_t           siiev;
    Int_t           nhits;
    Int_t           subdet[nMaxHits];
    Float_t         xh[nMaxHits];
    Float_t         xh_PulseHeight[nMaxHits];
    Float_t         yh[nMaxHits];
    Float_t         yh_PulseHeight[nMaxHits];
    Float_t         zh[nMaxHits];
    Int_t           itrack[nMaxHits];
    Int_t           S2_PulseHeight;
    Int_t           S2_TimeSpectrum;
    Int_t           LemmaTrigger_PulseHeight;
    Int_t           LemmaTrigger_TimeSpectrum;
    Int_t           Deva_PulseHeight[nDevaChannels];
    Int_t           Deva_TimeSpectrum[nDevaChannels];
    Int_t           Cherenkov_PulseHeight[nCherenkovChannels];
    Int_t           Cherenkov_TimeSpectrum[nCherenkovChannels];
    Int_t           S3_PulseHeight;
    Int_t           S3_TimeSpectrum;
    Int_t           PbGlass_PulseHeight;
    Int_t           PbGlass_TimeSpectrum;
    
    // --- Alignment variables
    bool  doAlignment;
    float x0SiAndMu[nSiLayers+1][3];
    float phi0Theta0Mu          [2];
    
    // --- Detector ID variables
    map<Int_t,Int_t> map_detID;
    const Int_t detIDs[nSiLayers+1]={10,20,30,40,51,50,56,55,70};

    ifstream myfile;
    bool m_debug;
};

#endif // EVENTBUILDER_H
