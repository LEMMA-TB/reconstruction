#include <iostream>
#include <stdio.h>
#include <math.h>
#include <TObject.h>
#include <TROOT.h>
#include <TRint.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TFrame.h>
#include <TMath.h>
#include <TFormula.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCut.h>
#include <TObject.h>
#include <TString.h>
#include <TNtuple.h>
#include <TSpline.h>

TFile *fsp, *fin;

TCanvas *TCspline;
 
TH2F *h_lincorr8[8];
TH1F *h1_lincorr8[8];
TGraph *g_lincorr[8];

TSpline3 *sp[8];

bool G2=1;
bool G4=0;

//static const int nbin=78;
static const int nbin=39;
//static const int nbin=19;

float x[nbin];
float y[nbin];

void create_spline8(){

  TCspline = new TCanvas("TCspline","spline",480,10,800,800);
  TCspline->SetFillColor(10);
  TCspline->SetBorderSize(1);
  TCspline->Divide(3,3);
  
  if(G2)
    fsp = new TFile("spline8G2.root","recreate");
  else if(G4)
    fsp = new TFile("spline8G4.root","recreate");
  else  
    fsp = new TFile("spline8.root","recreate");

  fin = new TFile("/data/radmu/Patt_Rec/INFO/prova_HIT_r811_4000000ev_nolin.root");
  
  for(int slo=0; slo<8; slo++){
    TCspline->cd(slo+1);
    char fileName[300];
    sprintf(fileName,"h_lincorr_slope%d",slo);
      h_lincorr8[slo]=(TH2F*)gDirectory->Get(fileName);
      if(G2)
	h_lincorr8[slo]->FitSlicesY(0,0,-1,0,"G2");
      else 
	if(G4)
	  h_lincorr8[slo]->FitSlicesY(0,0,-1,0,"G4");
	else
	  h_lincorr8[slo]->FitSlicesY();
      char fileName1[300];
      sprintf(fileName1,"h_lincorr_slope%d_1",slo);
      h1_lincorr8[slo]=(TH1F*)gDirectory->Get(fileName1);
      for(int ii=0;ii<nbin;ii++){
	if(G2){
	  x[ii]= 2.5 + 10.*ii;
	  y[ii]=h1_lincorr8[slo]->GetBinContent(2*ii+1);   
	}
	else 
	  if(G4){
	    x[ii]= 7.5 + 20.*ii;
	    y[ii]=h1_lincorr8[slo]->GetBinContent(4*ii+2);   
	  }
	  else{
	    x[ii]= 2.5 + 5.*ii;
	    y[ii]=h1_lincorr8[slo]->GetBinContent(ii+1);   
	  }
      }
      g_lincorr[slo]=new TGraph(nbin,x,y);
      char fileName2[300];
      sprintf(fileName2,"g_lincorr_slope%d",slo);
      g_lincorr[slo]->SetName(fileName2);
      
      g_lincorr[slo]->Draw("ap");
      
      char fileName3[300];
      sprintf(fileName3,"sp_%d",slo);
      sp[slo]=new TSpline3(fileName3,g_lincorr[slo]);
      sp[slo]->SetName(fileName3);     
      //TCspline->cd(slo+1);
      sp[slo]->Draw("same");
      
      fsp->cd();
      sp[slo]->Write();
      fin->cd();
      
  }
  
  fsp->cd();
  fsp->Close();
  
  fin->cd();
  fin->Close();
  
  return;
}
