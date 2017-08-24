#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TAxis.h"

void looper( TTree *tree, TH1F *HitsPerSegment_phi, TH1F *HitsPerSegment_theta, TH1F *HitsPerTrack );

void histos( Int_t RunNumber ) {
    
    ////////////////////////// open root file and initialize histograms //////////////////////////////////////
    
    char line[50];
//     char bigline[1000];
    
    sprintf( line, "output/Run_%4d_DT_pos.root", RunNumber );    
    TFile *file_pos = new TFile( line );
    
    sprintf( line, "output/Run_%4d_DT_neg.root", RunNumber ); 
    TFile *file_neg = new TFile( line );
    
    TTree *tree_pos = ( TTree * )file_pos->Get( "RADMU" );
    TTree *tree_neg = ( TTree * )file_neg->Get( "RADMU" );
    
    TH1F *HitsPerSegment_theta_pos = new TH1F( "HitsPerSegment_theta_pos", "Number of hits / segment in SL THETA", 11,  -0.5, 10.5 );
    TH1F *HitsPerSegment_phi_pos   = new TH1F( "HitsPerSegment_phi_pos",   "Number of hits / segment in SL PHI",   21,  -0.5, 20.5 );
    TH1F *HitsPerTrack_pos         = new TH1F( "HitsPerTrack_pos",         "Number of hits / track",               21,  -0.5, 20.5 );
    
    TH1F *HitsPerSegment_theta_neg = new TH1F( "HitsPerSegment_theta_neg", "Number of hits / segment in SL THETA", 11,  -0.5, 10.5 );
    TH1F *HitsPerSegment_phi_neg   = new TH1F( "HitsPerSegment_phi_neg",   "Number of hits / segment in SL PHI",   21,  -0.5, 20.5 );
    TH1F *HitsPerTrack_neg         = new TH1F( "HitsPerTrack_neg",         "Number of hits / track",               21,  -0.5, 20.5 );
    
    TH1F *Xhit_Track_posANDneg     = new TH1F( "Xhit_Track_posANDneg",     "X hit of reconstructed tracks",               100, -150, 150 ); //1bin = 3cm
    TH1F *Xhit_Track_posPLUSneg    = new TH1F( "Xhit_Track_posPLUSneg",    "Sum of X hits of reconstructed tracks",        50, -100, 100 );
    TH1F *Xhit_Track_posMINUSneg   = new TH1F( "Xhit_Track_posMINUSneg",   "Difference of X hits of reconstructed tracks", 50, -250, -50 );
    
    TH2F *Xhit_Track_posVSneg      = new TH2F( "Xhit_Track_posVSneg",      "X hit of pos track VS X hit of neg track", 50, 0, 150, 50, -150, 0 );
    
    TH1F *HitsPerSegmentMatching_theta_pos = new TH1F( "HitsPerSegmentMatching_theta_pos", "Number of hits in matching tracks / segment in SL THETA", 11,  -0.5, 10.5 );
    TH1F *HitsPerSegmentMatching_phi_pos   = new TH1F( "HitsPerSegmentMatching_phi_pos",   "Number of hits in matching tracks / segment in SL PHI",   21,  -0.5, 20.5 );
    TH1F *HitsPerTrackMatching_pos         = new TH1F( "HitsPerTrackMatching_pos",         "Number of hits in matching tracks / track",               21,  -0.5, 20.5 );
    TH1F *HitsPerSegmentMatching_theta_neg = new TH1F( "HitsPerSegmentMatching_theta_neg", "Number of hits in matching tracks / segment in SL THETA", 11,  -0.5, 10.5 );
    TH1F *HitsPerSegmentMatching_phi_neg   = new TH1F( "HitsPerSegmentMatching_phi_neg",   "Number of hits in matching tracks / segment in SL PHI",   21,  -0.5, 20.5 );
    TH1F *HitsPerTrackMatching_neg         = new TH1F( "HitsPerTrackMatching_neg",         "Number of hits in matching tracks / track",               21,  -0.5, 20.5 );
    
    ////////////////////////// fill histograms //////////////////////////////////////////////////////////////
    
    looper( tree_pos, HitsPerSegment_phi_pos, HitsPerSegment_theta_pos, HitsPerTrack_pos );
    looper( tree_neg, HitsPerSegment_phi_neg, HitsPerSegment_theta_neg, HitsPerTrack_neg );
    
    Int_t NumberOfSegments_pos;               tree_pos->SetBranchAddress( "SEG_ns", &NumberOfSegments_pos );
    Int_t NumberOfSegments_neg;               tree_neg->SetBranchAddress( "SEG_ns", &NumberOfSegments_neg );
    Int_t HitsPerSegmentMatching_pos[50];     tree_pos->SetBranchAddress( "SEG_sn",    HitsPerSegmentMatching_pos );
    Int_t HitsPerSegmentMatching_neg[50];     tree_neg->SetBranchAddress( "SEG_sn",    HitsPerSegmentMatching_neg );
    
    Float_t xhit1_pos[2];         tree_pos->SetBranchAddress( "SEG_xh1",            xhit1_pos );
    Float_t xhit1_neg[2];         tree_neg->SetBranchAddress( "SEG_xh1",            xhit1_neg );
    Float_t xhit2_pos[2];         tree_pos->SetBranchAddress( "SEG_xh2",            xhit2_pos );
    Float_t xhit2_neg[2];         tree_neg->SetBranchAddress( "SEG_xh2",            xhit2_neg );        
    Float_t xhit3_pos[2];         tree_pos->SetBranchAddress( "SEG_xh3",            xhit3_pos );
    Float_t xhit3_neg[2];         tree_neg->SetBranchAddress( "SEG_xh3",            xhit3_neg );
    Float_t xhit4_pos[2];         tree_pos->SetBranchAddress( "SEG_xh4",            xhit4_pos );
    Float_t xhit4_neg[2];         tree_neg->SetBranchAddress( "SEG_xh4",            xhit4_neg );
    Float_t xhit5_pos[2];         tree_pos->SetBranchAddress( "SEG_xh5",            xhit5_pos );
    Float_t xhit5_neg[2];         tree_neg->SetBranchAddress( "SEG_xh5",            xhit5_neg );
    Float_t xhit6_pos[2];         tree_pos->SetBranchAddress( "SEG_xh6",            xhit6_pos );
    Float_t xhit6_neg[2];         tree_neg->SetBranchAddress( "SEG_xh6",            xhit6_neg );
    Float_t xhit7_pos[2];         tree_pos->SetBranchAddress( "SEG_xh7",            xhit7_pos );
    Float_t xhit7_neg[2];         tree_neg->SetBranchAddress( "SEG_xh7",            xhit7_neg );
    Float_t xhit8_pos[2];         tree_pos->SetBranchAddress( "SEG_xh8",            xhit8_pos );
    Float_t xhit8_neg[2];         tree_neg->SetBranchAddress( "SEG_xh8",            xhit8_neg );
    
    Float_t* xhit[16];
    
    xhit[0]  = xhit1_pos; 
    xhit[1]  = xhit2_pos; 
    xhit[2]  = xhit3_pos; 
    xhit[3]  = xhit4_pos; 
    xhit[4]  = xhit5_pos; 
    xhit[5]  = xhit6_pos; 
    xhit[6]  = xhit7_pos; 
    xhit[7]  = xhit8_pos; 
    xhit[8]  = xhit1_neg;
    xhit[9]  = xhit2_neg;
    xhit[10] = xhit3_neg;
    xhit[11] = xhit4_neg;
    xhit[12] = xhit5_neg;
    xhit[13] = xhit6_neg;
    xhit[14] = xhit7_pos;
    xhit[15] = xhit8_neg;
    
    Int_t EvNum = tree_pos->GetEntries(); //same of tree_neg->GetEntries()
    
    for( Int_t ientry=0; ientry<EvNum; ++ientry ) {
        
        tree_pos->GetEntry( ientry );
        tree_neg->GetEntry( ientry );

        if ( NumberOfSegments_pos!=0 && NumberOfSegments_neg!=0 ) {
            
            Int_t npt_pos_phi;
            Int_t npt_neg_phi;
            Int_t npt_sum=0;
                        
            for( Int_t isegment=0; isegment<NumberOfSegments_pos; ++isegment ) {
                
                Int_t ch = static_cast<Int_t>( HitsPerSegmentMatching_pos[isegment]/1000 ); // <0 THETA 
                                                                                            // >0 PHI
                Int_t npt  = HitsPerSegmentMatching_pos[isegment] - 100*static_cast<Int_t>( HitsPerSegmentMatching_pos[isegment]/100 );
            
                if( ch<0 ){

                    npt = -npt;
                    HitsPerSegmentMatching_theta_pos->Fill( npt );
                    npt_sum += npt;
                
                } else if( ch>0 ) {

                    HitsPerSegmentMatching_phi_pos->Fill( npt );
                    npt_sum += npt;
                    npt_pos_phi = npt;
                
                } else {
                
                    std::cout << "ch = 0" << endl;
                
                }
                
                if( npt_sum ) HitsPerTrackMatching_pos->Fill( npt_sum );
            
            } // end loop on pos segments
            
            npt_sum=0;
            
            for( Int_t isegment=0; isegment<NumberOfSegments_neg; ++isegment ) {
                
                Int_t ch = static_cast<Int_t>( HitsPerSegmentMatching_neg[isegment]/1000 ); // <0 THETA 
                                                                                            // >0 PHI
                Int_t npt  = HitsPerSegmentMatching_neg[isegment] - 100*static_cast<Int_t>( HitsPerSegmentMatching_neg[isegment]/100 );
            
                if( ch<0 ){

                    npt = -npt;
                    HitsPerSegmentMatching_theta_neg->Fill( npt );
                    npt_sum += npt;
                
                } else if( ch>0 ) {

                    HitsPerSegmentMatching_phi_neg->Fill( npt );
                    npt_sum += npt;
                    npt_neg_phi = npt;
                
                } else {
                
                    std::cout << "ch = 0" << endl;
                
                }
                
                if( npt_sum ) HitsPerTrackMatching_neg->Fill( npt_sum );
            
            } //end loop on neg segments
            
//             if( npt_neg_phi>=7 && npt_pos_phi>=7 ) {
                
                Int_t i=0;
                while( ( xhit[i][0]==-999 || xhit[i+8][0]==-999 ) && i<8 ) {
                    ++i;
                }
            
//             sprintf( bigline, "%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n %1d \n",
//                      xhit1_pos[0], xhit2_pos[0], xhit3_pos[0], xhit4_pos[0], xhit5_pos[0], xhit6_pos[0], xhit7_pos[0], xhit8_pos[0],
//                      xhit1_neg[0], xhit2_neg[0], xhit3_neg[0], xhit4_neg[0], xhit5_neg[0], xhit6_neg[0], xhit7_neg[0], xhit8_neg[0],
//                      i );
//             
//             cout << bigline;

                Xhit_Track_posANDneg->Fill( -xhit[i][0] );
                Xhit_Track_posANDneg->Fill( -xhit[i+8][0] );
                Xhit_Track_posPLUSneg->Fill( -xhit[i][0]-xhit[i+8][0] );
                Xhit_Track_posMINUSneg->Fill( -xhit[i][0]+xhit[i+8][0] );
                Xhit_Track_posVSneg->Fill( -xhit[i+8][0], -xhit[i][0] );
                
//             } // end if npt
            
        } //end if number of segments
        
    } // end loop on entries
    
    
    ////////////////////////// draw histograms //////////////////////////////////////////////////////////////
    
    TPaveStats *ptstats;
    TAxis *Xaxis = new TAxis();
    TAxis *Yaxis = new TAxis();
    TCanvas *c1  = new TCanvas( "c1", "c1", 1500, 750 );
    
    c1->Divide(2, 2);
    c1->cd(1);

    ptstats = new TPaveStats(0.6885855,0.7175825,0.8888661,0.8783054,"brNDC");
    ptstats->SetBorderSize(1);
    ptstats->SetFillColor(0);
    ptstats->SetLineColor(0);
    ptstats->SetTextAlign(12);
    ptstats->SetTextFont(42);
    ptstats->SetOptStat(1110);
    ptstats->SetOptFit(0);
    HitsPerSegment_phi_pos->GetListOfFunctions()->Add( ptstats );
    ptstats->SetParent( HitsPerSegment_phi_pos );
    Xaxis = HitsPerSegment_phi_pos->GetXaxis();
    Xaxis->SetTitle( "Number of hits" );
    Yaxis = HitsPerSegment_phi_pos->GetYaxis();
    Yaxis->SetTitle( "Counts" );
    HitsPerSegment_phi_pos->SetLineWidth(2);
    HitsPerSegment_phi_pos->Draw();

    c1->cd(2);
    
    ptstats = new TPaveStats(0.6885855,0.7175825,0.8888661,0.8783054,"brNDC");
    ptstats->SetBorderSize(1);
    ptstats->SetFillColor(0);
    ptstats->SetLineColor(0);
    ptstats->SetTextAlign(12);
    ptstats->SetTextFont(42);
    ptstats->SetOptStat(1110);
    ptstats->SetOptFit(0);
    HitsPerSegment_phi_neg->GetListOfFunctions()->Add( ptstats );
    ptstats->SetParent( HitsPerSegment_phi_neg );
    Xaxis = HitsPerSegment_phi_neg->GetXaxis();
    Xaxis->SetTitle( "Number of hits" );
    Yaxis = HitsPerSegment_phi_neg->GetYaxis();
    Yaxis->SetTitle( "Counts" );
    HitsPerSegment_phi_neg->SetLineWidth(2);
    HitsPerSegment_phi_neg->SetLineColor(2);
    HitsPerSegment_phi_neg->Draw();
    
    c1->cd(3);
    
    ptstats = new TPaveStats(0.6885855,0.7175825,0.8888661,0.8783054,"brNDC");
    ptstats->SetBorderSize(1);
    ptstats->SetFillColor(0);
    ptstats->SetLineColor(0);
    ptstats->SetTextAlign(12);
    ptstats->SetTextFont(42);
    ptstats->SetOptStat(1110);
    ptstats->SetOptFit(0);
    HitsPerSegment_theta_pos->GetListOfFunctions()->Add( ptstats );
    ptstats->SetParent( HitsPerSegment_theta_pos );
    Xaxis = HitsPerSegment_theta_pos->GetXaxis();
    Xaxis->SetTitle( "Number of hits" );
    Yaxis = HitsPerSegment_theta_pos->GetYaxis();
    Yaxis->SetTitle( "Counts" );
    HitsPerSegment_theta_pos->SetLineWidth(2);
    HitsPerSegment_theta_pos->Draw();
    
    c1->cd(4);
    ptstats = new TPaveStats(0.6885855,0.7175825,0.8888661,0.8783054,"brNDC");
    ptstats->SetBorderSize(1);
    ptstats->SetFillColor(0);
    ptstats->SetLineColor(0);
    ptstats->SetTextAlign(12);
    ptstats->SetTextFont(42);
    ptstats->SetOptStat(1110);
    ptstats->SetOptFit(0);
    HitsPerSegment_theta_neg->GetListOfFunctions()->Add( ptstats );
    ptstats->SetParent( HitsPerSegment_theta_neg );
    Xaxis = HitsPerSegment_theta_neg->GetXaxis();
    Xaxis->SetTitle( "Number of hits" );
    Yaxis = HitsPerSegment_theta_neg->GetYaxis();
    Yaxis->SetTitle( "Counts" );
    HitsPerSegment_theta_neg->SetLineWidth(2);
    HitsPerSegment_theta_neg->SetLineColor(2);
    HitsPerSegment_theta_neg->Draw();
    
    TCanvas *c2 = new TCanvas( "c2", "c2", 1500, 750 );
    c2->Divide(2, 2);
    
    c2->cd(1);
    
    ptstats = new TPaveStats(0.6885855,0.7175825,0.8888661,0.8783054,"brNDC");
    ptstats->SetBorderSize(1);
    ptstats->SetFillColor(0);
    ptstats->SetLineColor(0);
    ptstats->SetTextAlign(12);
    ptstats->SetTextFont(42);
    ptstats->SetOptStat(1110);
    ptstats->SetOptFit(0);
    HitsPerTrack_pos->GetListOfFunctions()->Add( ptstats );
    ptstats->SetParent( HitsPerTrack_pos );
    Xaxis = HitsPerTrack_pos->GetXaxis();
    Xaxis->SetTitle( "Number of hits" );
    Yaxis = HitsPerTrack_pos->GetYaxis();
    Yaxis->SetTitle( "Counts" );
    HitsPerTrack_pos->SetLineWidth(2);
    HitsPerTrack_pos->Draw();
    
    c2->cd(2);
    
    ptstats = new TPaveStats(0.6885855,0.7175825,0.8888661,0.8783054,"brNDC");
    ptstats->SetBorderSize(1);
    ptstats->SetFillColor(0);
    ptstats->SetLineColor(0);
    ptstats->SetTextAlign(12);
    ptstats->SetTextFont(42);
    ptstats->SetOptStat(1110);
    ptstats->SetOptFit(0);
    HitsPerTrack_neg->GetListOfFunctions()->Add( ptstats );
    ptstats->SetParent( HitsPerTrack_neg );
    Xaxis = HitsPerTrack_neg->GetXaxis();
    Xaxis->SetTitle( "Number of hits" );
    Yaxis = HitsPerTrack_neg->GetYaxis();
    Yaxis->SetTitle( "Counts" );
    HitsPerTrack_neg->SetLineWidth(2);
    HitsPerTrack_neg->SetLineColor(2);
    HitsPerTrack_neg->Draw();
    
    TCanvas *c3 = new TCanvas( "c3", "c3");

    ptstats = new TPaveStats(0.6885855,0.7175825,0.8888661,0.8783054,"brNDC");
    ptstats->SetBorderSize(1);
    ptstats->SetFillColor(0);
    ptstats->SetLineColor(0);
    ptstats->SetTextAlign(12);
    ptstats->SetTextFont(42);
    ptstats->SetOptStat(1110);
    ptstats->SetOptFit(0);
    Xhit_Track_posANDneg->GetListOfFunctions()->Add( ptstats );
    ptstats->SetParent( Xhit_Track_posANDneg );
    Xaxis = Xhit_Track_posANDneg->GetXaxis();
    Xaxis->SetTitle( "X hit [cm]" );
    Yaxis = Xhit_Track_posANDneg->GetYaxis();
    Yaxis->SetTitle( "Counts" );
    Xhit_Track_posANDneg->SetLineWidth(2);
    Xhit_Track_posANDneg->Draw();
    
    TCanvas *c4 = new TCanvas( "c4", "c4" );
    
    ptstats = new TPaveStats(0.6885855,0.7175825,0.8888661,0.8783054,"brNDC");
    ptstats->SetBorderSize(1);
    ptstats->SetFillColor(0);
    ptstats->SetLineColor(0);
    ptstats->SetTextAlign(12);
    ptstats->SetTextFont(42);
    ptstats->SetOptStat(1110);
    ptstats->SetOptFit(0);
    Xhit_Track_posVSneg->GetListOfFunctions()->Add( ptstats );
    ptstats->SetParent( Xhit_Track_posVSneg );
    Xaxis = Xhit_Track_posVSneg->GetXaxis();
    Xaxis->SetTitle( "X hit neg [cm]" );
    Yaxis = Xhit_Track_posVSneg->GetYaxis();
    Yaxis->SetTitle( "X hit pos [cm]" );
    Xhit_Track_posVSneg->Draw( "colz" );
    
    TCanvas *c5 = new TCanvas( "c5", "c5" );
    
    ptstats = new TPaveStats(0.6885855,0.7175825,0.8888661,0.8783054,"brNDC");
    ptstats->SetBorderSize(1);
    ptstats->SetFillColor(0);
    ptstats->SetLineColor(0);
    ptstats->SetTextAlign(12);
    ptstats->SetTextFont(42);
    ptstats->SetOptStat(1110);
    ptstats->SetOptFit(0);
    Xhit_Track_posPLUSneg->GetListOfFunctions()->Add( ptstats );
    ptstats->SetParent( Xhit_Track_posPLUSneg );
    Xaxis = Xhit_Track_posPLUSneg->GetXaxis();
    Xaxis->SetTitle( "X hit pos + X hit neg [cm]" );
    Yaxis = Xhit_Track_posPLUSneg->GetYaxis();
    Yaxis->SetTitle( "Counts" );
    Xhit_Track_posPLUSneg->SetLineWidth(2);
    Xhit_Track_posPLUSneg->Draw();
    
    TCanvas *c6 = new TCanvas( "c6", "c6" );
    
    ptstats = new TPaveStats(0.6885855,0.7175825,0.8888661,0.8783054,"brNDC");
    ptstats->SetBorderSize(1);
    ptstats->SetFillColor(0);
    ptstats->SetLineColor(0);
    ptstats->SetTextAlign(12);
    ptstats->SetTextFont(42);
    ptstats->SetOptStat(1110);
    ptstats->SetOptFit(0);
    Xhit_Track_posMINUSneg->GetListOfFunctions()->Add( ptstats );
    ptstats->SetParent( Xhit_Track_posMINUSneg );
    Xaxis = Xhit_Track_posMINUSneg->GetXaxis();
    Xaxis->SetTitle( "X hit pos - X hit neg [cm]" );
    Yaxis = Xhit_Track_posMINUSneg->GetYaxis();
    Yaxis->SetTitle( "Counts" );
    Xhit_Track_posMINUSneg->SetLineWidth(2);
    Xhit_Track_posMINUSneg->Draw();
    
    TCanvas *c7 = new TCanvas( "c7", "c7", 1500, 750 );
    c7->Divide(2, 2);
    c7->cd(1);

    ptstats = new TPaveStats(0.6885855,0.7175825,0.8888661,0.8783054,"brNDC");
    ptstats->SetBorderSize(1);
    ptstats->SetFillColor(0);
    ptstats->SetLineColor(0);
    ptstats->SetTextAlign(12);
    ptstats->SetTextFont(42);
    ptstats->SetOptStat(1110);
    ptstats->SetOptFit(0);
    HitsPerSegmentMatching_phi_pos->GetListOfFunctions()->Add( ptstats );
    ptstats->SetParent( HitsPerSegmentMatching_phi_pos );
    Xaxis = HitsPerSegmentMatching_phi_pos->GetXaxis();
    Xaxis->SetTitle( "Number of hits" );
    Yaxis = HitsPerSegmentMatching_phi_pos->GetYaxis();
    Yaxis->SetTitle( "Counts" );
    HitsPerSegmentMatching_phi_pos->SetLineWidth(2);
    HitsPerSegmentMatching_phi_pos->Draw();

    c7->cd(2);
    
    ptstats = new TPaveStats(0.6885855,0.7175825,0.8888661,0.8783054,"brNDC");
    ptstats->SetBorderSize(1);
    ptstats->SetFillColor(0);
    ptstats->SetLineColor(0);
    ptstats->SetTextAlign(12);
    ptstats->SetTextFont(42);
    ptstats->SetOptStat(1110);
    ptstats->SetOptFit(0);
    HitsPerSegmentMatching_phi_neg->GetListOfFunctions()->Add( ptstats );
    ptstats->SetParent( HitsPerSegmentMatching_phi_neg );
    Xaxis = HitsPerSegmentMatching_phi_neg->GetXaxis();
    Xaxis->SetTitle( "Number of hits" );
    Yaxis = HitsPerSegmentMatching_phi_neg->GetYaxis();
    Yaxis->SetTitle( "Counts" );
    HitsPerSegmentMatching_phi_neg->SetLineWidth(2);
    HitsPerSegmentMatching_phi_neg->SetLineColor(2);
    HitsPerSegmentMatching_phi_neg->Draw();
    
    c7->cd(3);
    
    ptstats = new TPaveStats(0.6885855,0.7175825,0.8888661,0.8783054,"brNDC");
    ptstats->SetBorderSize(1);
    ptstats->SetFillColor(0);
    ptstats->SetLineColor(0);
    ptstats->SetTextAlign(12);
    ptstats->SetTextFont(42);
    ptstats->SetOptStat(1110);
    ptstats->SetOptFit(0);
    HitsPerSegmentMatching_theta_pos->GetListOfFunctions()->Add( ptstats );
    ptstats->SetParent( HitsPerSegmentMatching_theta_pos );
    Xaxis = HitsPerSegmentMatching_theta_pos->GetXaxis();
    Xaxis->SetTitle( "Number of hits" );
    Yaxis = HitsPerSegmentMatching_theta_pos->GetYaxis();
    Yaxis->SetTitle( "Counts" );
    HitsPerSegmentMatching_theta_pos->SetLineWidth(2);
    HitsPerSegmentMatching_theta_pos->Draw();
    
    c7->cd(4);
    ptstats = new TPaveStats(0.6885855,0.7175825,0.8888661,0.8783054,"brNDC");
    ptstats->SetBorderSize(1);
    ptstats->SetFillColor(0);
    ptstats->SetLineColor(0);
    ptstats->SetTextAlign(12);
    ptstats->SetTextFont(42);
    ptstats->SetOptStat(1110);
    ptstats->SetOptFit(0);
    HitsPerSegmentMatching_theta_neg->GetListOfFunctions()->Add( ptstats );
    ptstats->SetParent( HitsPerSegmentMatching_theta_neg );
    Xaxis = HitsPerSegmentMatching_theta_neg->GetXaxis();
    Xaxis->SetTitle( "Number of hits" );
    Yaxis = HitsPerSegmentMatching_theta_neg->GetYaxis();
    Yaxis->SetTitle( "Counts" );
    HitsPerSegmentMatching_theta_neg->SetLineWidth(2);
    HitsPerSegmentMatching_theta_neg->SetLineColor(2);
    HitsPerSegmentMatching_theta_neg->Draw();
    
    TCanvas *c8 = new TCanvas( "c8", "c8", 1500, 750 );
    c8->Divide(2, 2);
    
    c8->cd(1);
    
    ptstats = new TPaveStats(0.6885855,0.7175825,0.8888661,0.8783054,"brNDC");
    ptstats->SetBorderSize(1);
    ptstats->SetFillColor(0);
    ptstats->SetLineColor(0);
    ptstats->SetTextAlign(12);
    ptstats->SetTextFont(42);
    ptstats->SetOptStat(1110);
    ptstats->SetOptFit(0);
    HitsPerTrackMatching_pos->GetListOfFunctions()->Add( ptstats );
    ptstats->SetParent( HitsPerTrackMatching_pos );
    Xaxis = HitsPerTrackMatching_pos->GetXaxis();
    Xaxis->SetTitle( "Number of hits" );
    Yaxis = HitsPerTrackMatching_pos->GetYaxis();
    Yaxis->SetTitle( "Counts" );
    HitsPerTrackMatching_pos->SetLineWidth(2);
    HitsPerTrackMatching_pos->Draw();
    
    c8->cd(2);
    
    ptstats = new TPaveStats(0.6885855,0.7175825,0.8888661,0.8783054,"brNDC");
    ptstats->SetBorderSize(1);
    ptstats->SetFillColor(0);
    ptstats->SetLineColor(0);
    ptstats->SetTextAlign(12);
    ptstats->SetTextFont(42);
    ptstats->SetOptStat(1110);
    ptstats->SetOptFit(0);
    HitsPerTrackMatching_neg->GetListOfFunctions()->Add( ptstats );
    ptstats->SetParent( HitsPerTrackMatching_neg );
    Xaxis = HitsPerTrackMatching_neg->GetXaxis();
    Xaxis->SetTitle( "Number of hits" );
    Yaxis = HitsPerTrackMatching_neg->GetYaxis();
    Yaxis->SetTitle( "Counts" );
    HitsPerTrackMatching_neg->SetLineWidth(2);
    HitsPerTrackMatching_neg->SetLineColor(2);
    HitsPerTrackMatching_neg->Draw();
    
    ////////////////////////// write histograms //////////////////////////////////////////////////////////////
    
    sprintf( line, "output/Run_%4d_DT_histos.root", RunNumber );
    TFile *fileoutput = new TFile( line, "RECREATE" );
    
    HitsPerSegment_phi_pos->Write();
    HitsPerSegment_phi_neg->Write();
    HitsPerSegment_theta_pos->Write();
    HitsPerSegment_theta_neg->Write();
    Xhit_Track_posANDneg->Write();
    Xhit_Track_posPLUSneg->Write();
    Xhit_Track_posMINUSneg->Write();
    Xhit_Track_posVSneg->Write();
    HitsPerSegmentMatching_phi_pos->Write();
    HitsPerSegmentMatching_phi_neg->Write();
    HitsPerSegmentMatching_theta_pos->Write();
    HitsPerSegmentMatching_theta_neg->Write();
    
}



void looper( TTree *tree, TH1F *HitsPerSegment_phi, TH1F *HitsPerSegment_theta, TH1F *HitsPerTrack ) {
    
    Int_t NumberOfSegments;    tree->SetBranchAddress( "SEG_ns", &NumberOfSegments );
    Int_t HitsPerSegment[50];  tree->SetBranchAddress( "SEG_sn",    HitsPerSegment );
    
    Int_t EvNum = tree->GetEntries();
    
    //loop on entries
    for( Int_t ientry=0; ientry<EvNum; ++ientry ) {
        
        tree->GetEntry( ientry );
        
        Int_t npt_sum = 0;
        
        //loop on segments
        for( Int_t isegment=0; isegment<NumberOfSegments; ++isegment ) {
            
            Int_t ch = static_cast<Int_t>( HitsPerSegment[isegment]/1000 ); // <0 THETA 
                                                                            // >0 PHI
            Int_t npt  = HitsPerSegment[isegment] - 100*static_cast<Int_t>( HitsPerSegment[isegment]/100 );
            
            if( ch<0 ){

                npt = -npt;
                HitsPerSegment_theta->Fill( npt );
                npt_sum += npt;
                
            } else if( ch>0 ) {

                HitsPerSegment_phi->Fill( npt );
                npt_sum += npt;
                
            } else {
                
                std::cout << "ch = 0" << endl;
                
            }
            
        } //end loop on segments
        
        if( npt_sum ) HitsPerTrack->Fill( npt_sum );
        
    } //end loop on entries
    
}
