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

TH2F *h_lincorr[6][8];
TH1F *h1_lincorr[6][8];
TGraph *g_lincorr[6][8];

TSpline3 *sp[6][8];

bool G2=0;
bool G4=1;

//static const int nbin=78;
//static const int nbin=39;
static const int nbin=19;

float x[nbin];
float y[nbin];

void create_spline(){
  
  fsp = new TFile("spline.root","recreate");
  
  fin = new TFile("/data/radmu/Patt_Rec/INFO/prova_HIT_r811_4000000ev_nolin.root");
  
  for(int i=0; i<6; i++)  
    for(int slo=0; slo<8; slo++){
      char fileName[300];
      sprintf(fileName,"h_lincorr_seg%d_slope%d",i,slo);
      h_lincorr[i][slo]=(TH2F*)gDirectory->Get(fileName);
      if(G2)
	h_lincorr[i][slo]->FitSlicesY(0,0,-1,100,"G2");
      else 
	if(G4)
	  h_lincorr[i][slo]->FitSlicesY(0,0,-1,10,"G4");
	else
	  h_lincorr[i][slo]->FitSlicesY();
      char fileName1[300];
      sprintf(fileName1,"h_lincorr_seg%d_slope%d_1",i,slo);
      h1_lincorr[i][slo]=(TH1F*)gDirectory->Get(fileName1);
      for(int ii=0;ii<nbin;ii++){
	if(G2){
	  x[ii]= 2.5 + 10.*ii;
	  y[ii]=h1_lincorr[i][slo]->GetBinContent(2*ii+1);   
	}
	else 
	  if(G4){
	    x[ii]= 7.5 + 20.*ii;
	    y[ii]=h1_lincorr[i][slo]->GetBinContent(4*ii+2);   
	  }
	  else{
	    x[ii]= 2.5 + 5.*ii;
	    y[ii]=h1_lincorr[i][slo]->GetBinContent(ii+1);   
	  }
      }
      g_lincorr[i][slo]=new TGraph(nbin,x,y);
      char fileName2[300];
      sprintf(fileName2,"g_lincorr_seg%d_slope%d",i,slo);
      g_lincorr[i][slo]->SetName(fileName2);
      
      if(i==0 && slo==0)
	g_lincorr[i][slo]->Draw("ap");
      
      char fileName3[300];
      sprintf(fileName3,"sp_%d_%d",i,slo);
      sp[i][slo]=new TSpline3(fileName3,g_lincorr[i][slo]);
      sp[i][slo]->SetName(fileName3);     
      if(i==0 && slo==0)
	sp[i][slo]->Draw();
      fsp->cd();
      sp[i][slo]->Write();
      fin->cd();
      
    }
  
   fsp->cd();
   fsp->Close();

   fin->cd();
   fin->Close();
  
  return;
}
