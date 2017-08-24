#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TAxis.h"

void efficiency(  Int_t RunNumber ) {
    
    char line[50];
    
    sprintf( line, "output/Run_%d_SiDT_neg.root", RunNumber );
    
    TFile *file = new TFile( line );
    
    TTree *tree = ( TTree* )file->Get( "LEMMA" );
    
    const Int_t nMaxHits=100;
    
    Int_t nhits;            tree->SetBranchAddress( "nhits",  &nhits );
    Float_t xh[nMaxHits];   tree->SetBranchAddress( "xh",     xh     );
    Float_t yh[nMaxHits];   tree->SetBranchAddress( "yh",     yh     );
    Int_t subdet[nMaxHits]; tree->SetBranchAddress( "subdet", subdet );
    Int_t itrack[nMaxHits]; tree->SetBranchAddress( "itrack", itrack );

    Int_t Counter_All     = 0;
    Int_t Counter_Partial = 0;
    
    Int_t EvNum = tree->GetEntries();
    
    for( Int_t ientry=0; ientry<EvNum; ++ientry ) {
        
        Int_t Nhits_10=0, Nhits_20=0, Nhits_30=0, Nhits_40=0, Nhits_50=0, Nhits_55=0, Nhits_70=0; 
        
        tree->GetEntry( ientry );
        
        for ( Int_t i=0; i<nhits; ++i ) {
            
            if ( subdet[i]==10 ) {
                if( xh[i]!=-999 ) Nhits_10 = 1;
            } else if ( subdet[i]==20 ) {
                if( xh[i]!=-999 ) Nhits_20 = 1;
            } else if ( subdet[i]==30 ) {
                if( xh[i]!=-999 ) Nhits_30 = 1;
            } else if ( subdet[i]==40 ) {
                if( xh[i]!=-999 ) Nhits_40 = 1;
            } else if ( subdet[i]==50 ) {
                if( xh[i]!=-999 ) Nhits_50 = 1;
            } else if ( subdet[i]==55 ) {
                if( xh[i]!=-999 ) Nhits_55 = 1;
            } else if ( subdet[i]==70 ) {
                if( xh[i]!=-999 ) Nhits_70 = 1;
            }
        }
        
//         cout << Nhits_10 << " " << Nhits_20 << " " <<  Nhits_30 << " " <<  Nhits_40 << " " <<  Nhits_50 << " " <<  Nhits_55 << " " <<  Nhits_70 << endl; 
        
        if( Nhits_10 && Nhits_20 && Nhits_30 && Nhits_40 && Nhits_50 && Nhits_55 ) {
            Counter_Partial += 1;
//             Counter_Partial += (Nhits_10 + Nhits_20 + Nhits_30 + Nhits_40 + Nhits_50 + Nhits_55 );
//             cout << Counter_Partial << endl;
        }
        
        if( Nhits_10 && Nhits_20 && Nhits_30 && Nhits_40 && Nhits_50 && Nhits_55 && Nhits_70 ) {
            Counter_All += 1;
//             Counter_All += (Nhits_10 + Nhits_20 + Nhits_30 + Nhits_40 + Nhits_50 + Nhits_55 + Nhits_70 );
//             cout << Counter_All << endl;
        }

    }
    cout << "EvNum = " << EvNum << endl;
    cout << "Counter_Partial = " << Counter_Partial << endl;
    cout << "Counter_All = " << Counter_All << endl;
    cout << "Efficiency = " << Counter_All*1.0/Counter_Partial << endl;
    
}
