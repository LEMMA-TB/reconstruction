#include "RawHistos.h"
#include "HITCollection.h"


RawHistos::RawHistos() {

  // init useful variables
  //for debugging
  DEBUG_FLAG = 0;

  //occupancy histo cuts --> FIX tdc units or nseconds?
  minTime = 1200;
  maxTime = 2200;

  //hit difference histo cuts
  minTimeHit = 1400;
  maxTimeHit = 2000;

  //define map array and function to fill it
  nros =  1;
  nrob =  19;
  ntdc =  4;
  ncha = 32;

  // init histograms
  initHistos();

  // reset flags
  _flagCanvasTDC=0;
  _flagCanvasCH1=0;
  _flagCanvasCH2=0;
  _flagCanvasSLS=0;
  _flagCanvasHITS=0;

  return;
};

RawHistos::~RawHistos() {

  for(int i=0; i<4; i++)
  {
	delete hoccSL1_lay[   i];
	delete hoccSL2_lay[   i];
        delete hhit_diff[i];
  }

  for(int i=0; i<12; i++)
  {
	delete hoccCH1_lay[   i];
	delete hoccCH2_lay[   i];
  }

  for(int j=0; j<8; j++)
  	delete htime_sl[  j];

  delete hocc;
  delete htime; 
  delete htimech;
 
  if(_flagCanvasTDC){
    delete pad1;
    delete pad2;
    delete pad3;
  }

  if(_flagCanvasCH1){
    delete pad_hits_ch1;
    delete pad_time_ch1;
  }

  if(_flagCanvasCH2){
    delete pad_hits_ch2;
    delete pad_time_ch2;
  }

  if(_flagCanvasSLS){
    delete pad_hits_sls;
    delete pad_time_sls;
  }

  if(_flagCanvasHITS){
    delete pad_hits_hits;
  }

  if(c1)
    delete c1;

  if(ch1)
    delete ch1;

  if(ch2)
    delete ch2;

  if(sls)
    delete sls;

  if(hits)
    delete hits;


  return;
};

void RawHistos::initHistos(){
  
  //occupancy histograms: all layers from bottom to top 4+4+12+12
  int nCell = 72;
  int i = 0;

  for(i=0; i<4; i++)
  {
  	TString nameLayer( "hocc" );
	int layNum = i+1;	
	nameLayer += "_sl1_lay_";
	nameLayer += layNum;

	hoccSL1_lay[   i] = new TH1F( nameLayer,
                             nameLayer, nCell, 0.5, nCell + 0.5 );
  }

  for(i=0; i<4; i++)
  {
  	TString nameLayer( "hocc" );
	int layNum = i+1;	
	nameLayer += "_sl2_lay_";
	nameLayer += layNum;

	hoccSL2_lay[   i] = new TH1F( nameLayer,
                             nameLayer, nCell, 0.5, nCell + 0.5 );
  }

  for(i=0; i<12; i++)
  {
  	TString nameLayer( "hocc" );
	int layNum = i+1;	
	nameLayer += "_ch1_lay_";
	nameLayer += layNum;

	hoccCH1_lay[   i] = new TH1F( nameLayer,
                             nameLayer, nCell, 0.5, nCell + 0.5 );
  }

  for(i=0; i<12; i++)
  {
  	TString nameLayer( "hocc" );
	int layNum = i+1;	
	nameLayer += "_ch2_lay_";
	nameLayer += layNum;

	hoccCH2_lay[   i] = new TH1F( nameLayer,
                             nameLayer, nCell, 0.5, nCell + 0.5 );
  }

  //time boxes in SLs histograms from bottom to top 1+1+3+3
  int histoWidth = 2000;
  int histoEdgeL = 1000;
  int histoEdgeH = 3000;

  for(int j=0; j<8; j++)
  {
	TString tname( "htime" );
  	tname += "_sl_";
  	tname += j + 1;

    	htime_sl[  j] = new TH1F( tname,
                                tname,
                                histoWidth, histoEdgeL, histoEdgeH );
  }

  //second - first hit time histogram 
    for(i=0; i<4; i++)
  {
  	TString nameCh( "hhit_diff_ch_" );
	int chNum = i+8;	
	nameCh += chNum;

	hhit_diff[i] = new TH1F(nameCh,nameCh,1250,0.,1250.);
  }

  //TDC histograms
  hocc  = new TH2F("hocc","Channel vs TDC",200,0,99,66,-1,32);
  htime = new TH1F("htime","TimeBox",2000,1000,3000); 
  htimech = new TH2F("htimech","Channel vs time",2000,1000,3000,8000,0,4000); 

  return;
};

void RawHistos::buildCanvasTDC(TCanvas * c1) {

  if(DEBUG_FLAG)
    cout << "RawHistos::buildCanvasTDC" << endl; 
  //build the MINICRATE canvas with pads
  if(c1==0)
    c1 = new TCanvas("c1","RADMU monitor: TDC",1,10,800,800);
  c1->SetFillColor(10);
  //gStyle->SetFrameFillColor(18); 
  c1->SetBorderSize(2);
  TPaveText * title = new TPaveText(.2,0.96,.8,.995);
  title->SetFillColor(33);
  title->AddText("RADMU MONITOR: TDC hits");
  title->Draw();
  pad1 = new TPad("pad1","The pad with hits",0.05,0.50,0.95,0.92,19,1,1);
  pad1->Draw();
  pad2 = new TPad("pad2","The pad with times",0.05,0.05,0.50,0.45,19,1,1);
  pad2->Draw();
  pad3 = new TPad("pad3","The pad with times",0.50,0.05,0.95,0.45,19,1,1);
  pad3->Draw();

  // flag
  _flagCanvasTDC = true;

  return;
}

void RawHistos::buildCanvasCH2(TCanvas * ch2) {

  //Chamber 2 canvas
  if(ch2==0)
    ch2 = new TCanvas("ch2","RADMU monitor:CH 2",160,10,800,800);
  ch2->SetFillColor(10);
  ch2->SetBorderSize(1);
  
  TPaveText * title_ch2 = new TPaveText(.2,0.96,.8,.995);
  title_ch2->SetFillColor(33);
  title_ch2->AddText("RADMU MONITOR: chamber 2");
  title_ch2->Draw();

  pad_hits_ch2 = new TPad("pad_hits_ch2","The pad with hits",0.01,0.01,0.5,0.95,19,1,1);
  pad_hits_ch2->Divide(1,12);  
  pad_hits_ch2->Draw();

  pad_time_ch2 = new TPad("pad_time_ch2","The pad with times",0.5,0.01,0.99,0.95,19,1,1);
  pad_time_ch2->Divide(1,3);
  pad_time_ch2->Draw();

  // flag
  _flagCanvasCH2 = true;

  return;
}

void RawHistos::buildCanvasCH1(TCanvas * ch1) {

  //Chamber 1 canvas
  if(ch1==0)
    ch1 = new TCanvas("ch1","RADMU monitor:CH 1",320,10,800,800);
  ch1->SetFillColor(10);
  ch1->SetBorderSize(1);
 
  TPaveText * title_ch1 = new TPaveText(.2,0.96,.8,.995);
  title_ch1->SetFillColor(33);
  title_ch1->AddText("RADMU MONITOR: chamber 1");
  title_ch1->Draw();
  
  pad_hits_ch1 = new TPad("pad_hits_ch1","The pad with hits",0.01,0.01,0.5,0.95,19,1,1);
  pad_hits_ch1->Divide(1,12);  
  pad_hits_ch1->Draw();

  pad_time_ch1 = new TPad("pad_time_ch1","The pad with times",0.5,0.01,0.99,0.95,19,1,1);
  pad_time_ch1->Divide(1,3);
  pad_time_ch1->Draw();

  // flag
  _flagCanvasCH1 = true;

  return;
}

void RawHistos::buildCanvasSLS(TCanvas * sls) {

  //SLs canvas
  if(sls==0)
    sls = new TCanvas("sls","RADMU monitor:SLs",480,10,800,800);
  sls->SetFillColor(10);
  sls->SetBorderSize(1);
 
  TPaveText * title_sls = new TPaveText(.2,0.96,.8,.995);
  title_sls->SetFillColor(33);
  title_sls->AddText("RADMU MONITOR: SLs");
  title_sls->Draw();
  
  pad_hits_sls = new TPad("pad_hits_sls","The pad with hits",0.01,0.01,0.5,0.95,19,1,1);
  pad_hits_sls->Divide(1,12);  
  pad_hits_sls->Draw();

  pad_time_sls = new TPad("pad_time_sls","The pad with times",0.5,0.01,0.99,0.95,19,1,1);
  pad_time_sls->Divide(1,3);
  pad_time_sls->Draw();

  // flag
  _flagCanvasSLS = true;

  return;
}

void RawHistos::buildCanvasHITS(TCanvas * hits) {

  //hit difference canvas
  if(hits==0)
    hits = new TCanvas("hits","RADMU monitor: Hit Difference",640,10,800,800);
  hits->SetFillColor(10);
  hits->SetBorderSize(1);

  TPaveText * title_hits = new TPaveText(.2,0.96,.8,.995);
  title_hits->SetFillColor(33);
  title_hits->AddText("RADMU MONITOR: TDC-time Hit difference");
  title_hits->Draw();

  pad_hits_hits = new TPad("pad_hits_hits","The pad with hits time difference",0.01,0.01,0.95,0.95,19,1,1);
  pad_hits_hits->Divide(2,2);  
  pad_hits_hits->Draw();

  // flag
  _flagCanvasHITS = true;

  return;

};

void RawHistos::fillHistos(HITCollection * hits){

  if(DEBUG_FLAG)
    cout << "RawHistos::fillHistos : number of hits= " << hits->Get_NumberHITS() << endl;

  // loop on hits
  for(int i=0; i<hits->Get_NumberHITS(); i++) {

	// get hit and fill histos
   	HIT *hit = hits->hit(i);
	int se = hit->CH_ID();
        int sl = hit->SL_ID();
	int lay = hit->L_ID();
	int wire = hit->wire_ID();
	int TDC_Id = hit->TDC_ID();
	int ROB_Id = hit->ROB_ID();
	int channel = hit->channel_ID();

	float time = hit->rtime();
	//float dtime = hit->dtime();

	htime->Fill(time);

	//fill tdc histograms
	// ROB: 0..24;  TDC:0..3;  ch:0..31
	int tdcflag = TDC_Id + 4*ROB_Id;
	int chflag = channel + 32*tdcflag;
	if(time > minTime && time < maxTime)
		hocc->Fill(tdcflag,channel);
	htimech->Fill(time,chflag);

	// fill occupancy
	if(time > minTime && time < maxTime) {
		if(se == 8)
			hoccSL1_lay[lay-1]->Fill(wire);
		if(se == 9)
			hoccSL2_lay[lay-1]->Fill(wire);
		if(se == 10)
			hoccCH1_lay[lay-1 + (sl-1)*4]->Fill(wire);
		if(se == 11)
			hoccCH2_lay[lay-1 + (sl-1)*4]->Fill(wire);
	}
							
	if(DEBUG_FLAG)
		cout << "Filling : lay " << lay << "  sl " << sl << 
			" tube " << wire << endl;				

	htime_sl[sl]->Fill(time);

	if(DEBUG_FLAG)
		cout 	<< " Leading measurment " << " TDC " << setbase(10) << TDC_Id 
			<< " ch " << channel << " time " << time << endl;

  } //end hit loop

  return;

}

void RawHistos::updateCanvas(int flag) {
  // NB flag: 0=TDC, 1=ch1, 2=ch2, 3=sls, 4=hits

  if(DEBUG_FLAG)
    cout << "RawHistos::updateCanvas() " << ", _flagCanvasTDC " << _flagCanvasTDC
		<< ", _flagCanvasCH2 " << _flagCanvasCH2 << ", _flagCanvasCH1 " << _flagCanvasCH1 
		<< ", _flagCanvasSLS " << _flagCanvasSLS << ", _flagCanvasHITS " << _flagCanvasHITS 
		<< " *** flag " << flag << endl; 

  //TDC
  if(flag==0){
    pad1->cd();
    hocc->Draw();	
    pad1->Update();
    pad2->cd();
    htimech->Draw();
    pad2->Update(); 
    pad3->cd();
    htime->Draw();
    pad3->Update();
  }

  //chamber 2 plots
  if(flag==2){
    //ch2->cd();
    pad_hits_ch2->cd(1);
    hoccCH2_lay[11]->Draw();
    pad_hits_ch2->cd(2);
    hoccCH2_lay[10]->Draw();
    pad_hits_ch2->cd(3);
    hoccCH2_lay[9]->Draw();
    pad_hits_ch2->cd(4);
    hoccCH2_lay[8]->Draw();
    pad_hits_ch2->cd(5);
    hoccCH2_lay[7]->Draw();
    pad_hits_ch2->cd(6);
    hoccCH2_lay[6]->Draw();
    pad_hits_ch2->cd(7);
    hoccCH2_lay[5]->Draw();
    pad_hits_ch2->cd(8);
    hoccCH2_lay[4]->Draw();
    pad_hits_ch2->cd(9);
    hoccCH2_lay[3]->Draw();
    pad_hits_ch2->cd(10);
    hoccCH2_lay[2]->Draw();
    pad_hits_ch2->cd(11);
    hoccCH2_lay[1]->Draw();
    pad_hits_ch2->cd(12);
    hoccCH2_lay[0]->Draw();

    pad_time_ch2->cd(1);
    htime_sl[7]->Draw();
    pad_time_ch2->cd(2);
    htime_sl[6]->Draw();
    pad_time_ch2->cd(3);
    htime_sl[5]->Draw();

    //ch2->Update();	
  }

  //chamber 1 plots
  if(flag==1){
    pad_hits_ch1->cd(1);
    hoccCH1_lay[11]->Draw();
    pad_hits_ch1->cd(2);
    hoccCH1_lay[10]->Draw();
    pad_hits_ch1->cd(3);
    hoccCH1_lay[9]->Draw();
    pad_hits_ch1->cd(4);
    hoccCH1_lay[8]->Draw();
    pad_hits_ch1->cd(5);
    hoccCH1_lay[7]->Draw();
    pad_hits_ch1->cd(6);
    hoccCH1_lay[6]->Draw();
    pad_hits_ch1->cd(7);
    hoccCH1_lay[5]->Draw();
    pad_hits_ch1->cd(8);
    hoccCH1_lay[4]->Draw();
    pad_hits_ch1->cd(9);
    hoccCH1_lay[3]->Draw();
    pad_hits_ch1->cd(10);
    hoccCH1_lay[2]->Draw();
    pad_hits_ch1->cd(11);
    hoccCH1_lay[1]->Draw();
    pad_hits_ch1->cd(12);
    hoccCH1_lay[0]->Draw();
    pad_hits_ch1->Update();
    pad_time_ch1->cd(1);
    htime_sl[4]->Draw();
    pad_time_ch1->cd(2);
    htime_sl[3]->Draw();
    pad_time_ch1->cd(3);
    htime_sl[2]->Draw();
    //ch1->Update();
  }

  //SLs plots
  if(flag==3){
    //sls->cd();
    pad_hits_sls->cd(1);
    hoccSL2_lay[3]->Draw();
    pad_hits_sls->cd(2);
    hoccSL2_lay[2]->Draw();
    pad_hits_sls->cd(3);
    hoccSL2_lay[1]->Draw();
    pad_hits_sls->cd(4);
    hoccSL2_lay[0]->Draw();
    pad_hits_sls->cd(9);
    hoccSL1_lay[3]->Draw();
    pad_hits_sls->cd(10);
    hoccSL1_lay[2]->Draw();
    pad_hits_sls->cd(11);
    hoccSL1_lay[1]->Draw();
    pad_hits_sls->cd(12);
    hoccSL1_lay[0]->Draw();
    pad_hits_sls->Update();
    pad_time_sls->cd(1);
    htime_sl[1]->Draw();
    pad_time_sls->cd(3);
    htime_sl[0]->Draw();		
    //sls->Update();
  }

  // hit difference plot
  if(flag==4){
    pad_hits_hits->cd(1);
    hhit_diff[3]->Draw();
    pad_hits_hits->cd(2);
    hhit_diff[2]->Draw();
    pad_hits_hits->cd(3);
    hhit_diff[1]->Draw();
    pad_hits_hits->cd(4);
    hhit_diff[0]->Draw();
    pad_hits_hits->Update();
    //hits->Update();
  }

  if(DEBUG_FLAG)
    cout << "RawHistos::updateCanvas() : END" << endl;

  return;
}; 

void RawHistos::dumpHistos(char* fileName){

  TFile * file = new TFile(fileName,"RECREATE");
  file->cd();

  //TDC plots    
  hocc->Write();	
  htimech->Write();
  htime->Write();

  //chamber 2 plots
  hoccCH2_lay[11]->Write();
  hoccCH2_lay[10]->Write();
  hoccCH2_lay[9]->Write();
  hoccCH2_lay[8]->Write();
  hoccCH2_lay[7]->Write();
  hoccCH2_lay[6]->Write();
  hoccCH2_lay[5]->Write();
  hoccCH2_lay[4]->Write();
  hoccCH2_lay[3]->Write();
  hoccCH2_lay[2]->Write();
  hoccCH2_lay[1]->Write();
  hoccCH2_lay[0]->Write();

  htime_sl[7]->Write();
  htime_sl[6]->Write();
  htime_sl[5]->Write();

  //chamber 1 plots
  hoccCH1_lay[11]->Write();
  hoccCH1_lay[10]->Write();
  hoccCH1_lay[9]->Write();
  hoccCH1_lay[8]->Write();
  hoccCH1_lay[7]->Write();
  hoccCH1_lay[6]->Write();
  hoccCH1_lay[5]->Write();
  hoccCH1_lay[4]->Write();
  hoccCH1_lay[3]->Write();
  hoccCH1_lay[2]->Write();
  hoccCH1_lay[1]->Write();
  hoccCH1_lay[0]->Write();
  htime_sl[4]->Write();
  htime_sl[3]->Write();
  htime_sl[2]->Write();

  //SLs plots
  hoccSL2_lay[3]->Write();
  hoccSL2_lay[2]->Write();
  hoccSL2_lay[1]->Write();
  hoccSL2_lay[0]->Write();
  hoccSL1_lay[3]->Write();
  hoccSL1_lay[2]->Write();
  hoccSL1_lay[1]->Write();
  hoccSL1_lay[0]->Write();
  htime_sl[1]->Write();
  htime_sl[0]->Write();		

  // hit difference plot
  hhit_diff[3]->Write();
  hhit_diff[2]->Write();
  hhit_diff[1]->Write();
  hhit_diff[0]->Write();

  file->Close();
  delete file;

  return;
};

