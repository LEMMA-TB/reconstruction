#ifndef OCCUPANCY_H
#define OCCUPANCY_H

#include <iostream>
#include "HITCollection.h"
#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TString.h"

class Occupancy
{
 public:
    Occupancy();
    ~Occupancy();
  
    // member functions
    void initVariables();
    void initHistos();

    void fillOccHistos( HITCollection *hits );
    void saveOccHistos();

 private:
    
    // histos
    TH1F *occLay[12];
    TH1F *Number_of_hits_per_event;

    // histos' variables
    int occLayWidth;
    int occLayEdgeL;
    int occLayEdgeH;

    //
    TFile *outputFile;
};

#endif
