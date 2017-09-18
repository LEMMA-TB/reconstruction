#include <iostream>
#include <iomanip>
#include <string>

#include <vector>
#include "TH1.h"
#include "TH2D.h"
#include "TH2.h"
#include "TH3F.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TChain.h"


void HistStyle(TH1 *hist);

bool debug = false;
bool cutEvents = false;
bool cutIron = false;

//set default plain style, for grace!
gROOT->Reset();
//gStyle->SetErrorX(0);
gROOT->SetStyle("Plain");
gStyle->SetOptStat(1);
gStyle->SetOptStat(1111);
gStyle->SetOptFit(1);

void dtplots(std::string file, int maxEvents)
{
  TChain * T = new TChain("RADMU");
//  T->Add("output/Radmufit_r4390_1000007ev_PR.root");
//  T->Add("output/Radmufit_r4390_3000001ev_PR.root");

  T->Add(file.data());

  //book histograms
  // DTBX HITS
  TH1F * hnhits = new TH1F("hnhits","hit number per event",30,0.,30);
  TH1F * hocc_sl1 = new TH1F("hocc_sl1","Occupancy SL 1",100,0,100);
  TH1F * hocc_sl2 = new TH1F("hocc_sl2","Occupancy SL 2",100,0,100);
  TH1F * hocc_sl3 = new TH1F("hocc_sl3","Occupancy SL 3",100,0,100);
  TH1F * htime_sl1 = new TH1F("htime_sl1","Drift time SL 1",600,-100,500);
  TH1F * htime_sl2 = new TH1F("htime_sl2","Drift time SL 2",600,-100,500);
  TH1F * htime_sl3 = new TH1F("htime_sl3","Drift time SL 3",600,-100,500);

  TH1F * hSignal = new TH1F("hSignal","Signal per sl",3,1,4);
  TH1F * hNoise = new TH1F("hNoise","Noise per sl",3,1,4);
  TH1F * hSoverN = new TH1F("hSoverN","Signal over Noise per sl",3,1,4);

  TH1F * hocc_lay[12];
  for (int il=0; il<12; il++ ){
      std::string name = "htime_lay";
      std::stringstream sstm;
      sstm << name << il;
      std::string result = sstm.str();
      hocc_lay[il] = new TH1F(result.c_str(),"Occupancy LAYER",70,0,70);
  }

  //T0 FIT
  TH1F * ht0_ch1 = new TH1F("ht0_ch1","chamber 1: event by event t0",100,-50.,50);

  //SEG
  TH1F * hnpt = new TH1F("hnpt","npt parameter",30,0.,30);
  TH1F * hslope_ch1_phi = new TH1F("hslope_ch1_phi","slope in chamber 1 phi view",400,-2.,2.);
  TH1F * hslope_ch1_theta = new TH1F("hslope_ch1_theta","slope in chamber 1 theta view",400,-2.,2.);
  TH1F * hX_ch1 = new TH1F("hX_ch1","intercept in chamber 1 phi view",150,-150.,150.);
  TH1F * hZ_ch1 = new TH1F("hZ_ch1","intercept in chamber 1 theta view",150,-150.,150.);

  TH2F * hslope1vsX1 = new TH2F("hslope1vsX1","slope ch 1 vs X ch 1",300,-150.,150.,100,-2.,2.);

  TH1F * hres_ch1 = new TH1F("hres_ch1","residuals chamber 1",100,0.,0.1);
  TH1F * hres_ch1_phi1 = new TH1F("hres_ch1_phi1","residuals chamber 1 - PHI 1",1000,-0.1,0.1);
  TH1F * hres_ch1_phi2 = new TH1F("hres_ch1_phi2","residuals chamber 1 - PHI 2",1000,-0.1,0.1);
  TH1F * hres_ch1_theta = new TH1F("hres_ch1_theta","residuals chamber 1 - THETA",1000,-0.1,0.1);

  //DTBX block
  int Nhits, rawT_in[120], rawT[120];
  int lay[120], tube[120];
  T->SetBranchAddress( "DTBX_nhit", &Nhits);
  T->SetBranchAddress( "DTBX_time_in", rawT_in);
  T->SetBranchAddress( "DTBX_time", rawT);
  T->SetBranchAddress( "DTBX_lay", lay );
  T->SetBranchAddress( "DTBX_tube", tube );

  //SEG block
  Int_t iseg;
  Int_t segN[50];
  Float_t segX[50], segS[50], segK[50], segT[50];
  Int_t l1c[50], l2c[50], l3c[50], l4c[50], l5c[50], l6c[50], l7c[50], l8c[50];
  Float_t r1[50], r2[50], r3[50], r4[50], r5[50], r6[50], r7[50], r8[50];

  // associate variables to branches
  T->SetBranchAddress( "SEG_ns", &iseg );
  T->SetBranchAddress( "SEG_sx",  segX );
  T->SetBranchAddress( "SEG_ss",  segS );
  T->SetBranchAddress( "SEG_sk",  segK );
  T->SetBranchAddress( "SEG_sn",  segN );
  T->SetBranchAddress( "SEG_t0",  segT );

  T->SetBranchAddress( "SEG_1r",  r1 );
  T->SetBranchAddress( "SEG_2r",  r2 );
  T->SetBranchAddress( "SEG_3r",  r3 );
  T->SetBranchAddress( "SEG_4r",  r4 );
  T->SetBranchAddress( "SEG_5r",  r5 );
  T->SetBranchAddress( "SEG_6r",  r6 );
  T->SetBranchAddress( "SEG_7r",  r7 );
  T->SetBranchAddress( "SEG_8r",  r8 );


  //get number of events
  Int_t evnum = T->GetEntries();
  //cout << " number of events " << evnum << endl;

  //loop over the events
  for (Int_t ientry=0 ; ientry < TMath::Min(maxEvents,evnum); ientry ++)
  {
    if(floor(ientry/1000)*1000==ientry)
                cout << "Event " << ientry << endl;

    //read tree
    T->GetEntry(ientry);

    // ** RAW HITS
    hnhits->Fill(Nhits);
    for(int ih=0; ih < Nhits; ih++){
        // layer
        int layer = lay[ih] - 1300;
        int sl;
        if(layer<=4)
            sl=1;
        else if(layer>=5 && layer<=8)
            sl=2;
        else if(layer>=8)
            sl=3;

        // fill signa to noise histograms
        float tdriftMin = 0.;
        float tdriftMax = 400.;

        if(rawT[ih]>tdriftMin && rawT[ih]<tdriftMax)
            hSignal->Fill(sl);
        else
            hNoise->Fill(sl);


        // fill occupancy
        hocc_lay[layer-1]->Fill(tube[ih]);

        if(layer <=4){
            hocc_sl1->Fill(tube[ih]);
            htime_sl1->Fill(rawT[ih]);
        }
        if(layer <=8 && layer >=5){
            hocc_sl2->Fill(tube[ih]);
            htime_sl2->Fill(rawT[ih]);
        }
        if(layer >=9){
            hocc_sl3->Fill(tube[ih]);
            htime_sl3->Fill(rawT[ih]);
        }
    }

    // ** SEG
    Float_t slope_ch1_phi = 9999.;
    Float_t slope_ch1_theta = 9999.;
    Float_t X_ch1 = 9999.;
    Float_t Z_ch1 = 9999.;
    Float_t T0_ch1 = 9999.;
    Float_t res_ch1[12] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}; //residuals in micron

    //cout << " number of segments " << iseg << endl;
    for ( Int_t is = 0; is < iseg; is++ ){
        //for LNL 2 ch settings
        Int_t ch = static_cast<Int_t>(segN[is]/1000);

        Int_t npt  = segN[is] - 100*static_cast<Int_t>(segN[is]/100);
        if(ch<0){
            ch = -ch;
            npt = -npt;
        }
        hnpt->Fill(npt);
        //cout << "segN[is] " << segN[is] << " ch " << ch << " npt " << npt << endl;

        //select 3,4-hit theta segment on ch 1
        if(segN[is]==-1304 || segN[is]==-1303)
        {
            slope_ch1_theta = segS[is];
            Z_ch1 = segX[is];

            res_ch1[4] = r1[is];
            res_ch1[5] = r2[is];
            res_ch1[6] = r3[is];
            res_ch1[7] = r4[is];
        }

        //select segments over the whole chamber
        if( npt<=8 && npt>=6 && segN[is]>0)
        {
            //select 8 hit fit
            //if(npt==8)
            {
                if(ch==1)
                {
                    slope_ch1_phi = segS[is];
                    X_ch1 = segX[is];
                    T0_ch1 = segT[is];

                    res_ch1[0] = r1[is];
                    res_ch1[1] = r2[is];
                    res_ch1[2] = r3[is];
                    res_ch1[3] = r4[is];
                    res_ch1[8] = r5[is];
                    res_ch1[9] = r6[is];
                    res_ch1[10] = r7[is];
                    res_ch1[11] = r8[is];

                    }//end ch=1
            }// end if npt=8
        } // end select chamber
    }// end segment loop

    if(slope_ch1_phi != 9999.)
    {
        // *** residuals analysis
        Int_t npt_ch1 = 0;
        Float_t sigma_fix = 250.; //micron
        Float_t res_ch1_sq_sum = 0.;
        Float_t sigma_ch1 = 0;
        Float_t sigmaTimes = 2.;

        for(int ir=0; ir<4; ir++)
            hres_ch1_phi1->Fill(res_ch1[ir]);

        for(int ir=8; ir<12; ir++)
            hres_ch1_phi2->Fill(res_ch1[ir]);

        // npt fit
        for(int l=0; l<12; l++)
        {
            if(res_ch1[l] != -999.)
                npt_ch1 += 1;

            //residual square sum
            res_ch1_sq_sum += res_ch1[l]*res_ch1[l];
        }

        // sigma
        sigma_ch1 = TMath::Sqrt( res_ch1_sq_sum / (npt_ch1-5) ); // micron

//cout << " sigma_ch1 " << sigma_ch1 << " sigma_ch2 " << sigma_ch2 << endl;
//cout << " npt_ch1 " << npt_ch1 << " npt_ch2 " << npt_ch2 << endl;
//cout << " res_ch1_sq_sum " << res_ch1_sq_sum << " res_ch2_sq_sum " << res_ch2_sq_sum << endl;

        hres_ch1->Fill(sigma_ch1);
    }

    if(slope_ch1_theta != 9999.)
    {
        hslope_ch1_theta->Fill(slope_ch1_theta);
        hZ_ch1->Fill(Z_ch1);

        for(int ir=4; ir<8; ir++)
            hres_ch1_theta->Fill(res_ch1[ir]);
    }

    if(slope_ch1_phi != 9999.)
    {
        hslope_ch1_phi->Fill(slope_ch1_phi);
        hX_ch1->Fill(X_ch1);
        hslope1vsX1->Fill(X_ch1,slope_ch1_phi);
        ht0_ch1->Fill(T0_ch1);
    }

    //end segment analysis
    }//end event loop



  //draw histograms
  TCanvas * csn = new TCanvas("csn","csn",1500,1000);
  csn->cd(1);
//  csn->Divide(2,2);
//  hSignal->DrawClone();
//  csn->cd(2);
//  hNoise->DrawClone();
//  csn->cd(3);
  hSignal->Divide(hNoise);
  hSignal->Draw();
  hSignal->SetTitle("Signal / Noise ");
  Double_t norm = hSignal->GetEntries();
  hSignal->Scale(1/norm);

  TCanvas * chit = new TCanvas("chit","chit",1500,1000);
  chit->Divide(2,2);
  chit->cd(1);
  hnhits->Draw();
  chit->cd(2);
  hocc_sl1->Draw();
  chit->cd(3);
  hocc_sl2->Draw();
  chit->cd(4);
  hocc_sl3->Draw();

  TCanvas * ctime = new TCanvas("ctime","ctime",1500,1000);
  ctime->Divide(3,1);
  ctime->cd(1);
  htime_sl1->Draw();
  ctime->cd(2);
  htime_sl2->Draw();
  ctime->cd(3);
  htime_sl3->Draw();

  TPaveText *tpt1 = new TPaveText(0.85,0.4,0.95,0.7,"brNDC");
  tpt1->SetFillColor(kWhite);
  tpt1->SetTextColor(kRed);
  tpt1->SetTextAlign(12);
  tpt1->SetBorderSize(1);
  tpt1->AddText("PHI 1");
  TPaveText *tpt2 = new TPaveText(0.85,0.4,0.95,0.7,"brNDC");
  tpt2->SetFillColor(kWhite);
  tpt2->SetTextColor(kRed);
  tpt2->SetTextAlign(12);
  tpt2->SetBorderSize(1);
  tpt2->AddText("THETA");
  TPaveText *tpt3 = new TPaveText(0.85,0.4,0.95,0.7,"brNDC");
  tpt3->SetFillColor(kWhite);
  tpt3->SetTextColor(kRed);
  tpt3->SetTextAlign(12);
  tpt3->SetBorderSize(1);
  tpt3->AddText("PHI 3");

  TCanvas * cocc_lay = new TCanvas("cocc_lay","cocc_lay",1600,1000);
  cocc_lay->Divide(3,4);
  for(int i=0; i<3; i++){
      cocc_lay->cd(1+i);
      hocc_lay[3+4*i]->Draw();
      cocc_lay->cd(4+i);
      hocc_lay[2+4*i]->Draw();
      cocc_lay->cd(7+i);
      hocc_lay[1+4*i]->Draw();
      cocc_lay->cd(10+i);
      hocc_lay[0+4*i]->Draw();
  }
  cocc_lay->cd(1);
  tpt1->Draw();
  cocc_lay->cd(2);
  tpt2->Draw();
  cocc_lay->cd(3);
  tpt3->Draw();

  TCanvas * cfit = new TCanvas("cfit","cfit",1600,1000);
  cfit->Divide(3,2);
  cfit->cd(1);
  hnpt->Draw();
  cfit->cd(2);
  hslope_ch1_phi->SetXTitle("tg(phi)");
  hslope_ch1_phi->Draw();
  cfit->cd(3);
  hX_ch1->SetXTitle("X (cm)");
  hX_ch1->Draw();
  cfit->cd(4);
  ht0_ch1->Draw();
  cfit->cd(5);
  hslope_ch1_theta->SetXTitle("tg(theta)");
  hslope_ch1_theta->Draw();
  cfit->cd(6);
  hZ_ch1->SetXTitle("Z (cm)");
  hZ_ch1->Draw();
  TPaveText *tpc = new TPaveText(0.65,0.55,0.98,0.65,"brNDC");
  tpc->SetFillColor(kWhite);
  tpc->SetTextAlign(12);
  tpc->AddText("CHAMBER 1");
  tpc->AddText("Slope and Intercept");
  tpc->Draw();

  TCanvas * cres = new TCanvas("cres","cres",1500,500);
  cres->Divide(3,1);
  cres->cd(1);
  hres_ch1_phi1->Draw();
  cres->cd(2);
  hres_ch1_phi2->Draw();
  cres->cd(3);
  hres_ch1_theta->Draw();
}

void HistStyle(TH1 *hist)
{
  //gROOT->Reset();
  //gStyle->SetErrorX(0);
  //gROOT->SetStyle("Plain");
  //gStyle->SetOptStat(1);
  //gStyle->SetOptTitle(kFALSE);

  //hist->SetLineWidth(2);
  hist->GetXaxis()->SetTitleOffset(1.);
  hist->GetYaxis()->SetTitleOffset(1.3);
  hist->GetZaxis()->SetTitleOffset(1.);
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetZaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleSize(0.04);
  hist->GetYaxis()->SetTitleSize(0.04);
  hist->GetZaxis()->SetTitleSize(0.04);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetZaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetLabelSize(0.04);
  hist->GetZaxis()->SetLabelSize(0.04);
  hist->SetStats(1);

  //hist->SetLineColor(kYellow);
  //hist->SetMarkerStyle(21);
  hist->SetMarkerSize(0.9);

}


void dtplots_newntuple(std::string file, int maxEvents)
{
    TChain * T = new TChain("LEMMA-DT");
    T->Add(file.data());

    //book histograms
    // DTBX HITS
    TH1F * hnhits = new TH1F("hnhits","hit number per event",30,0.,30);
    TH1F * hnhits_sl1 = new TH1F("hnhits_sl1","hit number per event SL 1",30,0.,30);
    TH1F * hnhits_sl2 = new TH1F("hnhits_sl2","hit number per event SL 2",30,0.,30);
    TH1F * hnhits_sl3 = new TH1F("hnhits_sl3","hit number per event SL 3",30,0.,30);
    TH1F * hocc_sl1 = new TH1F("hocc_sl1","Occupancy SL 1",100,0,100);
    TH1F * hocc_sl2 = new TH1F("hocc_sl2","Occupancy SL 2",100,0,100);
    TH1F * hocc_sl3 = new TH1F("hocc_sl3","Occupancy SL 3",100,0,100);
    TH1F * htime_sl1 = new TH1F("htime_sl1","Drift time SL 1",600,-100,500);
    TH1F * htime_sl2 = new TH1F("htime_sl2","Drift time SL 2",600,-100,500);
    TH1F * htime_sl3 = new TH1F("htime_sl3","Drift time SL 3",600,-100,500);

    TH1F * hSignal = new TH1F("hSignal","Signal per sl",3,1,4);
    TH1F * hNoise = new TH1F("hNoise","Noise per sl",3,1,4);
    TH1F * hSoverN = new TH1F("hSoverN","Signal over Noise per sl",3,1,4);

    TH1F * hocc_lay[12];
    for (int il=0; il<12; il++ ){
        std::string name = "htime_lay";
        std::stringstream sstm;
        sstm << name << il;
        std::string result = sstm.str();
        hocc_lay[il] = new TH1F(result.c_str(),"Occupancy LAYER",70,0,70);
    }

    //T0 FIT
    TH1F * ht0_ch1 = new TH1F("ht0_ch1","chamber 1: event by event t0",100,-50.,50);

    //SEG
    TH1F * hnpt = new TH1F("hnpt","npt parameter",30,0.,30);
    TH1F * hslope_ch1_phi = new TH1F("hslope_ch1_phi","slope in chamber 1 phi view",400,-2.,2.);
    TH1F * hslope_ch1_theta = new TH1F("hslope_ch1_theta","slope in chamber 1 theta view",400,-2.,2.);
    TH1F * hX_ch1 = new TH1F("hX_ch1","intercept in chamber 1 phi view",150,-150.,150.);
    TH1F * hZ_ch1 = new TH1F("hZ_ch1","intercept in chamber 1 theta view",150,-150.,150.);

    TH2F * hslope1vsX1 = new TH2F("hslope1vsX1","slope ch 1 vs X ch 1",300,-150.,150.,100,-2.,2.);

    TH1F * hres_ch1 = new TH1F("hres_ch1","residuals chamber 1",100,0.,0.1);
    TH1F * hres_ch1_phi1 = new TH1F("hres_ch1_phi1","residuals chamber 1 - PHI 1",1000,-0.1,0.1);
    TH1F * hres_ch1_phi2 = new TH1F("hres_ch1_phi2","residuals chamber 1 - PHI 2",1000,-0.1,0.1);
    TH1F * hres_ch1_theta = new TH1F("hres_ch1_theta","residuals chamber 1 - THETA",1000,-0.1,0.1);

    //DTBX block
    int Nhits;
    float rawT[5000];
    int lay[5000], tube[5000];
    T->SetBranchAddress( "DT_NHits", &Nhits);
    T->SetBranchAddress( "DT_lay", lay );
    T->SetBranchAddress( "DT_tube", tube );
    T->SetBranchAddress( "DT_dtime", rawT);

    /*
    //SEG block
    Int_t iseg;
    Int_t segN[50];
    Float_t segX[50], segS[50], segK[50], segT[50];
    Int_t l1c[50], l2c[50], l3c[50], l4c[50], l5c[50], l6c[50], l7c[50], l8c[50];
    Float_t r1[50], r2[50], r3[50], r4[50], r5[50], r6[50], r7[50], r8[50];

    // associate variables to branches
    T->SetBranchAddress( "SEG_ns", &iseg );
    T->SetBranchAddress( "SEG_sx",  segX );
    T->SetBranchAddress( "SEG_ss",  segS );
    T->SetBranchAddress( "SEG_sk",  segK );
    T->SetBranchAddress( "SEG_sn",  segN );
    T->SetBranchAddress( "SEG_t0",  segT );

    T->SetBranchAddress( "SEG_1r",  r1 );
    T->SetBranchAddress( "SEG_2r",  r2 );
    T->SetBranchAddress( "SEG_3r",  r3 );
    T->SetBranchAddress( "SEG_4r",  r4 );
    T->SetBranchAddress( "SEG_5r",  r5 );
    T->SetBranchAddress( "SEG_6r",  r6 );
    T->SetBranchAddress( "SEG_7r",  r7 );
    T->SetBranchAddress( "SEG_8r",  r8 );
    */

    //get number of events
    Int_t evnum = T->GetEntries();
    //cout << " number of events " << evnum << endl;

    //loop over the events
    for (Int_t ientry=0 ; ientry < TMath::Min(maxEvents,evnum); ientry ++)
    {
        if(floor(ientry/1000)*1000==ientry)
            cout << "Event " << ientry << endl;

        //read tree
        T->GetEntry(ientry);

        // ** RAW HITS
        hnhits->Fill(Nhits);

        int NHits_sl1 = 0;
        int NHits_sl2 = 0;
        int NHits_sl3 = 0;

        for(int ih=0; ih < Nhits; ih++){
            // layer
            int layer = lay[ih];
            int sl;
            if(layer<=4)
                sl=1;
            else if(layer>=5 && layer<=8)
                sl=2;
            else if(layer>=8)
                sl=3;

            // fill signa to noise histograms
            float tdriftMin = 0.;
            float tdriftMax = 400.;

            if(rawT[ih]>tdriftMin && rawT[ih]<tdriftMax)
                hSignal->Fill(sl);
            else
                hNoise->Fill(sl);

            // fill occupancy
            hocc_lay[layer-1]->Fill(tube[ih]);

            if(layer <=4){
                hocc_sl1->Fill(tube[ih]);
                htime_sl1->Fill(rawT[ih]);
                NHits_sl1 += 1;
            }
            if(layer <=8 && layer >=5){
                hocc_sl2->Fill(tube[ih]);
                htime_sl2->Fill(rawT[ih]);
                NHits_sl2 += 1;
            }
            if(layer >=9){
                hocc_sl3->Fill(tube[ih]);
                htime_sl3->Fill(rawT[ih]);
                NHits_sl3 += 1;
            }
        } // end hit loop

        hnhits_sl1->Fill(NHits_sl1);
        hnhits_sl2->Fill(NHits_sl2);
        hnhits_sl3->Fill(NHits_sl3);

        /*
        // ** SEG
        Float_t slope_ch1_phi = 9999.;
        Float_t slope_ch1_theta = 9999.;
        Float_t X_ch1 = 9999.;
        Float_t Z_ch1 = 9999.;
        Float_t T0_ch1 = 9999.;
        Float_t res_ch1[12] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}; //residuals in micron

        //cout << " number of segments " << iseg << endl;
        for ( Int_t is = 0; is < iseg; is++ ){
            //for LNL 2 ch settings
            Int_t ch = static_cast<Int_t>(segN[is]/1000);

            Int_t npt  = segN[is] - 100*static_cast<Int_t>(segN[is]/100);
            if(ch<0){
                ch = -ch;
                npt = -npt;
            }
            hnpt->Fill(npt);
            //cout << "segN[is] " << segN[is] << " ch " << ch << " npt " << npt << endl;

            //select 3,4-hit theta segment on ch 1
            if(segN[is]==-1304 || segN[is]==-1303)
            {
                slope_ch1_theta = segS[is];
                Z_ch1 = segX[is];

                res_ch1[4] = r1[is];
                res_ch1[5] = r2[is];
                res_ch1[6] = r3[is];
                res_ch1[7] = r4[is];
            }

            //select segments over the whole chamber
            if( npt<=8 && npt>=6 && segN[is]>0)
            {
                //select 8 hit fit
                //if(npt==8)
                {
                    if(ch==1)
                    {
                        slope_ch1_phi = segS[is];
                        X_ch1 = segX[is];
                        T0_ch1 = segT[is];

                        res_ch1[0] = r1[is];
                        res_ch1[1] = r2[is];
                        res_ch1[2] = r3[is];
                        res_ch1[3] = r4[is];
                        res_ch1[8] = r5[is];
                        res_ch1[9] = r6[is];
                        res_ch1[10] = r7[is];
                        res_ch1[11] = r8[is];

                    }//end ch=1
                }// end if npt=8
            } // end select chamber
        }// end segment loop

        if(slope_ch1_phi != 9999.)
        {
            // *** residuals analysis
            Int_t npt_ch1 = 0;
            Float_t sigma_fix = 250.; //micron
            Float_t res_ch1_sq_sum = 0.;
            Float_t sigma_ch1 = 0;
            Float_t sigmaTimes = 2.;

            for(int ir=0; ir<4; ir++)
                hres_ch1_phi1->Fill(res_ch1[ir]);

            for(int ir=8; ir<12; ir++)
                hres_ch1_phi2->Fill(res_ch1[ir]);

            // npt fit
            for(int l=0; l<12; l++)
            {
                if(res_ch1[l] != -999.)
                    npt_ch1 += 1;

                //residual square sum
                res_ch1_sq_sum += res_ch1[l]*res_ch1[l];
            }

            // sigma
            sigma_ch1 = TMath::Sqrt( res_ch1_sq_sum / (npt_ch1-5) ); // micron

            //cout << " sigma_ch1 " << sigma_ch1 << " sigma_ch2 " << sigma_ch2 << endl;
            //cout << " npt_ch1 " << npt_ch1 << " npt_ch2 " << npt_ch2 << endl;
            //cout << " res_ch1_sq_sum " << res_ch1_sq_sum << " res_ch2_sq_sum " << res_ch2_sq_sum << endl;

            hres_ch1->Fill(sigma_ch1);
        }

        if(slope_ch1_theta != 9999.)
        {
            hslope_ch1_theta->Fill(slope_ch1_theta);
            hZ_ch1->Fill(Z_ch1);

            for(int ir=4; ir<8; ir++)
                hres_ch1_theta->Fill(res_ch1[ir]);
        }

        if(slope_ch1_phi != 9999.)
        {
            hslope_ch1_phi->Fill(slope_ch1_phi);
            hX_ch1->Fill(X_ch1);
            hslope1vsX1->Fill(X_ch1,slope_ch1_phi);
            ht0_ch1->Fill(T0_ch1);
        }

        //end segment analysis
        */
    }//end event loop



    //draw histograms
    TCanvas * csn = new TCanvas("csn","csn",1500,1000);
    csn->cd(1);
    //  csn->Divide(2,2);
    //  hSignal->DrawClone();
    //  csn->cd(2);
    //  hNoise->DrawClone();
    //  csn->cd(3);
    hSignal->Divide(hNoise);
    hSignal->Draw();
    hSignal->SetTitle("Signal / Noise ");
    Double_t norm = hSignal->GetEntries();
    hSignal->Scale(1/norm);

    TCanvas * chit = new TCanvas("chit","chit",1500,1000);
    chit->Divide(3,2);
//    chit->cd(1);
//    hnhits->Draw();
    chit->cd(1);
    hnhits_sl1->Draw();
    chit->cd(2);
    hnhits_sl2->Draw();
    chit->cd(3);
    hnhits_sl3->Draw();
    chit->cd(4);
    hocc_sl1->Draw();
    chit->cd(5);
    hocc_sl2->Draw();
    chit->cd(6);
    hocc_sl3->Draw();

    TCanvas * ctime = new TCanvas("ctime","ctime",1500,1000);
    ctime->Divide(3,1);
    ctime->cd(1);
    htime_sl1->Draw();
    ctime->cd(2);
    htime_sl2->Draw();
    ctime->cd(3);
    htime_sl3->Draw();

    TPaveText *tpt1 = new TPaveText(0.85,0.4,0.95,0.7,"brNDC");
    tpt1->SetFillColor(kWhite);
    tpt1->SetTextColor(kRed);
    tpt1->SetTextAlign(12);
    tpt1->SetBorderSize(1);
    tpt1->AddText("PHI 1");
    TPaveText *tpt2 = new TPaveText(0.85,0.4,0.95,0.7,"brNDC");
    tpt2->SetFillColor(kWhite);
    tpt2->SetTextColor(kRed);
    tpt2->SetTextAlign(12);
    tpt2->SetBorderSize(1);
    tpt2->AddText("THETA");
    TPaveText *tpt3 = new TPaveText(0.85,0.4,0.95,0.7,"brNDC");
    tpt3->SetFillColor(kWhite);
    tpt3->SetTextColor(kRed);
    tpt3->SetTextAlign(12);
    tpt3->SetBorderSize(1);
    tpt3->AddText("PHI 3");

    TCanvas * cocc_lay = new TCanvas("cocc_lay","cocc_lay",1600,1000);
    cocc_lay->Divide(3,4);
    for(int i=0; i<3; i++){
        cocc_lay->cd(1+i);
        hocc_lay[3+4*i]->Draw();
        cocc_lay->cd(4+i);
        hocc_lay[2+4*i]->Draw();
        cocc_lay->cd(7+i);
        hocc_lay[1+4*i]->Draw();
        cocc_lay->cd(10+i);
        hocc_lay[0+4*i]->Draw();
    }
    cocc_lay->cd(1);
    tpt1->Draw();
    cocc_lay->cd(2);
    tpt2->Draw();
    cocc_lay->cd(3);
    tpt3->Draw();

    TCanvas * cfit = new TCanvas("cfit","cfit",1600,1000);
    cfit->Divide(3,2);
    cfit->cd(1);
    hnpt->Draw();
    cfit->cd(2);
    hslope_ch1_phi->SetXTitle("tg(phi)");
    hslope_ch1_phi->Draw();
    cfit->cd(3);
    hX_ch1->SetXTitle("X (cm)");
    hX_ch1->Draw();
    cfit->cd(4);
    ht0_ch1->Draw();
    cfit->cd(5);
    hslope_ch1_theta->SetXTitle("tg(theta)");
    hslope_ch1_theta->Draw();
    cfit->cd(6);
    hZ_ch1->SetXTitle("Z (cm)");
    hZ_ch1->Draw();
    TPaveText *tpc = new TPaveText(0.65,0.55,0.98,0.65,"brNDC");
    tpc->SetFillColor(kWhite);
    tpc->SetTextAlign(12);
    tpc->AddText("CHAMBER 1");
    tpc->AddText("Slope and Intercept");
    tpc->Draw();

    TCanvas * cres = new TCanvas("cres","cres",1500,500);
    cres->Divide(3,1);
    cres->cd(1);
    hres_ch1_phi1->Draw();
    cres->cd(2);
    hres_ch1_phi2->Draw();
    cres->cd(3);
    hres_ch1_theta->Draw();
}

void HistStyle(TH1 *hist)
{
    //gROOT->Reset();
    //gStyle->SetErrorX(0);
    //gROOT->SetStyle("Plain");
    //gStyle->SetOptStat(1);
    //gStyle->SetOptTitle(kFALSE);

    //hist->SetLineWidth(2);
    hist->GetXaxis()->SetTitleOffset(1.);
    hist->GetYaxis()->SetTitleOffset(1.3);
    hist->GetZaxis()->SetTitleOffset(1.);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetZaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetTitleSize(0.04);
    hist->GetYaxis()->SetTitleSize(0.04);
    hist->GetZaxis()->SetTitleSize(0.04);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetZaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetLabelSize(0.04);
    hist->GetZaxis()->SetLabelSize(0.04);
    hist->SetStats(1);

    //hist->SetLineColor(kYellow);
    //hist->SetMarkerStyle(21);
    hist->SetMarkerSize(0.9);

}


