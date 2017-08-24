#include <iostream>
#include "TFile.h"
#include "TTree.h"

using namespace std;

void FindEvent( int RunNumber ) {
    
    TH1F *xhit_double_mu = new TH1F( "xhit_double_mu", "X hit of double tracks after event building", 240, -1200, +1200 );
    
    char line[50];
    sprintf( line, "Run_%d_SibigDT_merged.root", RunNumber );
    
    TFile *file = new TFile( line );
    
    TTree *tree = ( TTree* )file->Get( "LEMMA" );
    
    Int_t   iev;         tree->SetBranchAddress( "iev",    &iev    );
    Int_t   nhits;       tree->SetBranchAddress( "nhits",  &nhits  );
    Int_t   subdet[100]; tree->SetBranchAddress( "subdet", subdet  );
    Float_t xh[100];     tree->SetBranchAddress( "xh",     xh      );
    Float_t yh[100];     tree->SetBranchAddress( "yh",     yh      );
    Int_t   itrack[100]; tree->SetBranchAddress( "itrack", itrack  );
    
    Float_t xhit_pos[8];
    Float_t xhit_neg[8];
    Float_t yhit_pos[4];
    Float_t yhit_neg[4];
        
    Int_t double_mu_counter=0;
    
    Int_t numEvent = tree->GetEntries();

    for( Int_t ientry=0; ientry<numEvent; ++ientry ) {
        
        tree->GetEntry( ientry );
        
        Int_t xhit_pos_counter=0;
        Int_t xhit_neg_counter=0;
        Int_t yhit_pos_counter=0;
        Int_t yhit_neg_counter=0;
        
        for( Int_t i=0; i<8; ++i ) {
            xhit_pos[i]=-999;
            xhit_neg[i]=-999;
        }
        
        for( Int_t i=0; i<4; ++i ) {
            yhit_pos[i]=-999;
            yhit_neg[i]=-999;
        }
        
        for( Int_t isubdet=0; isubdet<nhits; ++isubdet ) {
            
            if( subdet[isubdet] == 70 ) {

                if ( itrack[isubdet] == -1 ) {
                    
                    if ( xh[isubdet]>+37.5 ) {
                        xhit_neg[xhit_neg_counter]=xh[isubdet];
                        xhit_neg_counter++;
                    }
                    
                    if ( yh[isubdet]!=-999 ) {
                        yhit_neg[yhit_neg_counter]=yh[isubdet];
                        yhit_neg_counter++;
                    }
                    
                } // end if itrack == -1
                
                if ( itrack[isubdet] == +1 ) {
                    
                    if ( xh[isubdet]!=-999 && xh[isubdet]<-37.5 ) {
                        xhit_pos[xhit_pos_counter]=xh[isubdet];
                        xhit_pos_counter++;
                    }
                    
                    if ( yh[isubdet]!=-999 ) {
                        yhit_pos[yhit_pos_counter]=yh[isubdet];
                        yhit_pos_counter++;
                    }   
                    
                } // end if itrack == +1

            } // end if subdet == 70
            
        } // end loop on hits

        if( xhit_pos_counter>=6 && yhit_pos_counter>=3 && xhit_neg_counter>=6 && yhit_neg_counter>=3 ) {
            ++double_mu_counter;
            
            /*cout << "******* EVENT " << iev << " ******* " << endl;
            for( Int_t i=0; i<8; ++i ) { 
                
                cout << "Neg xh = " << xhit_neg[i] << "\t Pos xh = " << xhit_pos[i] << endl;
                
            }
            
            cout << endl;
            
            for( Int_t i=0; i<4; ++i) {
                 
                 cout << "Neg yh = " << yhit_neg[i] << "\t Pos yh = " << yhit_pos[i] << endl;
                 
            }
            
            cout << endl;*/
            
            xhit_double_mu->Fill( xhit_neg[0]*10 );
            xhit_double_mu->Fill( xhit_pos[0]*10 );
            
        }        
        
    } // end loop on entries

    cout << "Pairs of tracks = " << double_mu_counter << endl;
    
    TCanvas *c1 = new TCanvas();
    TPaveStats *ptstats =  new TPaveStats(0.6885855,0.7175825,0.8888661,0.8783054,"brNDC");
    ptstats->SetBorderSize(1);
    ptstats->SetFillColor(0);
    ptstats->SetLineColor(0);
    ptstats->SetTextAlign(12);
    ptstats->SetTextFont(42);
    ptstats->SetOptStat(1110);
    ptstats->SetOptFit(0);
    xhit_double_mu->GetListOfFunctions()->Add( ptstats );
    ptstats->SetParent( xhit_double_mu );
    TAxis *Xaxis = xhit_double_mu->GetXaxis();
    Xaxis->SetTitle( "X hit [mm]" );
    TAxis *Yaxis = xhit_double_mu->GetYaxis();
    Yaxis->SetTitle( "Counts" );
    xhit_double_mu->SetLineWidth(2);
    xhit_double_mu->Draw();
    
}
