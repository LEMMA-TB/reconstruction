#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <stdio.h>

using namespace std;

void seg( TString filename1, TString filename2 ) {
    
    TFile *file1 = new TFile( filename1 );
    TFile *file2 = new TFile( filename2 );
    
    TTree *t1 = ( TTree* )file1->Get( "RADMU" );
    TTree *t2 = ( TTree* )file2->Get( "RADMU" );
    
    int event1;
    int event2;
    
    t1->SetBranchAddress( "EVENT", &event1 );
    
    t2->SetBranchAddress( "EVENT", &event2 );
    
    int evnum1 = t1->GetEntries();
    int evnum2 = t2->GetEntries();
    
    for( int ientry = 0; ientry < TMath::min( evnum1, evnum2 ); ientry++ ) {
        
        t1->GetEntry( ientry );
        t2->GetEntry( ientry );
        
        
        
    }
    
}
