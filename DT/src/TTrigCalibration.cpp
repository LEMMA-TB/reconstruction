#include "TTrigCalibration.h"

//ALTEA  
//int N_ROB[19]={0, 1, 2, 6, 3, 4, 5, 12, 13, 14, 24, 15, 16, 17, 18, 19, 21, 22, 23};
int N_ROB[6]={0, 1, 2, 4, 5, 6};

TTrigCalibration::TTrigCalibration( int runTrig ) {

  // init useful variables
  //for debugging
  DEBUG_TTRIG = 0;
  DEBUG_TTRIGBIN = 0;

  TString tfilename = "./output/DTTimeBoxFitter_";
  tfilename += runTrig;
  tfilename += ".root"; //ALTEA
  
  hDebugFile = new TFile( tfilename, "RECREATE"); //ALTEA

  // init histograms and variables
  initVariables();
  initHistos();

  // reset flags
  //  _flagFillHistos=0;
  _flagCanvasTROB=0;
  _flagCanvasTSL=0;

  return;
}

TTrigCalibration::~TTrigCalibration() {
  
  for(int i=0; i<8; i++)
    {
      delete htbox_ROB[i];
    }
  
  for(int i=0; i<3; i++)
    {
      delete htbox_SL[i];
    }
  
  if(_flagCanvasTROB){
    delete pad_time_ROB;
  }
  
  if(_flagCanvasTSL){
    delete pad_time_SL;
  }
  
  if(TBOX_ROB)
    delete TBOX_ROB;

  if(TBOX_SL)
    delete TBOX_SL;
  
  
  return;
}

void TTrigCalibration::initHistos(){
  
  if(DEBUG_TTRIG)
    cout << "TTrigCalibration::initHistos" << endl; 

  for(int j=0; j<8; j++)
    {
      TString tBname( "htbox" );
      tBname += "_ROB_";
      tBname += j ;
      
      htbox_ROB[  j] = new TH1F( tBname,
				 tBname,
				 hROBWidth, hROBEdgeL, hROBEdgeH );
    }

  for(int j=0; j<3; j++)
    {
      TString tBname( "htbox" );
      tBname += "_SL_";
      tBname += j+1 ;
      
      htbox_SL[  j] = new TH1F( tBname,
				tBname,
				hSLWidth, hSLEdgeL, hSLEdgeH );
    }


  if(DEBUG_TTRIG)
    cout << "TTrigCalibration::initHistos DONE" << endl; 

  return;

}


void TTrigCalibration::initVariables(){
  
  if(DEBUG_TTRIG)
    cout << "TTrigCalibration::initVariables" << endl; 

  for(int i=0;i<6;i++) TTrig_ROB_mean[i]=0.;
  for(int i=0;i<6;i++) TTrig_ROB_RMS[i]=0.;
  for(int i=0;i<3;i++) TTrig_mean[i]=0.;
  for(int i=0;i<3;i++) TTrig_RMS[i]=0.;
 
  //ALTEA
  //time boxes (ns) for each ROB
  hROBWidth = 4000;
  hROBEdgeL = 500;
  hROBEdgeH = 4500;

  //time boxes (ns) for each SL
  hSLWidth = 4000;
  hSLEdgeL = 500;
  hSLEdgeH = 4500;

  //cut on raw time
  rtime_min=3000;
  rtime_max=4000;

  if(DEBUG_TTRIG)
    cout << "TTrigCalibration::initVariables DONE" << endl; 

  return;

}

void TTrigCalibration::buildCanvasTROB(TCanvas * TBOX_ROB) {

  if(DEBUG_TTRIG)
    cout << "TTrigCalibration::buildCanvasTROB" << endl; 

  if(TBOX_ROB==0)
    TBOX_ROB = new TCanvas("TBOX_ROB","RADMU monitor: TBOX ROB",640,10,800,800);
  TBOX_ROB->SetFillColor(10);
  TBOX_ROB->SetBorderSize(2);

  TPaveText * title_TBOX_ROB = new TPaveText(.2,0.96,.8,.995);
  title_TBOX_ROB->SetFillColor(33);
  title_TBOX_ROB->AddText("RADMU MONITOR: TBOX ROB");
  title_TBOX_ROB->Draw();
  pad_time_ROB = new TPad("pad_time_ROB","The pad with time box for each ROB",0.01,0.01,0.95,0.95,19,1,1);
  pad_time_ROB->Divide(5,5);  
  pad_time_ROB->Draw();

  // flag
  _flagCanvasTROB = true;

  return;
}

void TTrigCalibration::buildCanvasTSL(TCanvas * TBOX_SL) {

  if(DEBUG_TTRIG)
    cout << "TTrigCalibration::buildCanvasTSL" << endl; 

  if(TBOX_SL==0)
    TBOX_SL = new TCanvas("TBOX_SL","RADMU monitor: TBOX SL",640,10,800,800);
  TBOX_SL->SetFillColor(10);
  TBOX_SL->SetBorderSize(2);
  
  TPaveText * title_TBOX_SL = new TPaveText(.2,0.96,.8,.995);
  title_TBOX_SL->SetFillColor(33);
  title_TBOX_SL->AddText("RADMU MONITOR: TBOX SL");
  title_TBOX_SL->Draw();

  pad_time_SL = new TPad("pad_time_SL","The pad with time box for each SL",0.01,0.01,0.95,0.95,19,1,1);
  pad_time_SL->Divide(3,3);  
  pad_time_SL->Draw();


  // flag
  _flagCanvasTSL = true;

  return;
}



void TTrigCalibration::fillHistos(HITCollection * hits){
  
  if(DEBUG_TTRIG)
    cout << "TTrigCalibration::fillHistos : number of hits= " << hits->Get_NumberHITS() << endl;
  
  // loop on hits
  for(int i=0; i<hits->Get_NumberHITS(); i++) {
    
    // get hit and fill histos
    HIT *hit = hits->hit(i);
    if(hit->rtime()>rtime_min && hit->rtime()<rtime_max){

       float dtime = hit->rtime() - hit->t0();
       //float dtime = hit->rtime();

      int ROB_Id = hit->ROB_ID();
      htbox_ROB[ROB_Id]->Fill(dtime);
      
      //int sl=0;
      //if(hit->CH_ID()==8) sl=0;
      //if(hit->CH_ID()==9) sl=1;
      //if(hit->CH_ID()==10){
	//if(hit->SL_ID()==1)sl=2;
	//if(hit->SL_ID()==2)sl=3;
	//if(hit->SL_ID()==3)sl=4;
      //}
      //if(hit->CH_ID()==11){
	//if(hit->SL_ID()==1)sl=5;
	//if(hit->SL_ID()==2)sl=6;
	//if(hit->SL_ID()==3)sl=7;
      //}
	  
      //htbox_SL[sl]->Fill(dtime);
      
      //ALTEA
      htbox_SL[hit->SL_ID()-1]->Fill(dtime);
      
    } //close cut on time
  } //end hit loop
  
  return;
  
}



void TTrigCalibration::computeTTrig(int flag) {
  // NB flag: 0=TBOX_ROB, 1=TBOX_SL
  
  if(DEBUG_TTRIG)
    cout << "TTrigCalibration::computeTTrig() " << ", _flagCanvasTROB " << _flagCanvasTROB
	 << ", _flagCanvasTSL " << _flagCanvasTSL << " *** flag " << flag << endl; 
  
  //TROB
  if(flag==0){
    //    TBOX_ROB->cd();
    double TTRIG_ROB_ORD[6];
    double RMS_ROB_ORD[6];
    
    for(int i=0;i<8;i++){
      pad_time_ROB->cd(i+1);
      htbox_ROB[i]->Draw();	
      
      if( htbox_ROB[i]->GetEntries()>0 ) {
	
	double TTRIG_ROB,RMS_ROB;
	fitTimeBox(htbox_ROB[i],TTRIG_ROB,RMS_ROB);
	
	if(DEBUG_TTRIG)
	  cout<<" ROB n."<<i<<"  TTRIG (ns) = "<<TTRIG_ROB
	      <<"  RMS (ns) = "<<RMS_ROB<<endl;
	
	if(i<=2) {
	  N_ROB[i]=i;
	  TTRIG_ROB_ORD[i]=TTRIG_ROB;
	  RMS_ROB_ORD[i]=RMS_ROB;
	}
	if(i==6) {
	  N_ROB[3]=i;
	  TTRIG_ROB_ORD[3]=TTRIG_ROB;
	  RMS_ROB_ORD[3]=RMS_ROB;
	}
	if(i>=3 && i<=5) {
	  N_ROB[i+1]=i;
	  TTRIG_ROB_ORD[i+1]=TTRIG_ROB;
	  RMS_ROB_ORD[i+1]=RMS_ROB;
	}
	if(i>=12 && i<=14) {
	  N_ROB[i-5]=i;
	  TTRIG_ROB_ORD[i-5]=TTRIG_ROB;
	  RMS_ROB_ORD[i-5]=RMS_ROB;
	}
	if(i==24) {
	  N_ROB[10]=i;
	  TTRIG_ROB_ORD[10]=TTRIG_ROB;
	  RMS_ROB_ORD[10]=RMS_ROB;
	}
	if(i>=15 && i<=19) {
	  N_ROB[i-4]=i;
	  TTRIG_ROB_ORD[i-4]=TTRIG_ROB;
	  RMS_ROB_ORD[i-4]=RMS_ROB;
	}
	if(i>=21 && i<=23) {
	  N_ROB[i-5]=i;
	  TTRIG_ROB_ORD[i-5]=TTRIG_ROB;
	  RMS_ROB_ORD[i-5]=RMS_ROB;
	}
      }
    }
    for(int i=0;i<19;i++){
      cout<<" ROB n."<<N_ROB[i]<<"  TTRIG (ns) = "<<TTRIG_ROB_ORD[i]<<"  RMS (ns) = "<<RMS_ROB_ORD[i]<<endl;
      if(RMS_ROB_ORD[i]<20.){
	TTrig_ROB_mean[i]=TTRIG_ROB_ORD[i];
	TTrig_ROB_RMS[i]=RMS_ROB_ORD[i];
      }
    }
    
    //    TBOX_ROB->Update();
    
  }
  
  //ALTEA nel nostro caso flag==1

  if(flag==1){
    // TBOX_SL->cd();
    double TTRIG_SL_ORD[3];
    double RMS_SL_ORD[3];

    for(int i=0;i<3;i++){
      pad_time_SL->cd(i+1);
      htbox_SL[i]->Write();	
      
	if( htbox_SL[i]->GetEntries()>0 ) {
	double TTRIG_SL;
	double RMS_SL;
	fitTimeBox(htbox_SL[i],TTRIG_SL,RMS_SL);
    	TTRIG_SL_ORD[i]=TTRIG_SL;
	RMS_SL_ORD[i]=RMS_SL;
	
      }
    }
    for(int i=0;i<3;i++){
      cout<<" SL n."<<i+1<<"  TTRIG (ns) = "<<TTRIG_SL_ORD[i]<<"  RMS (ns) = "<<RMS_SL_ORD[i]<<endl;

      //ALTEA io non voglio che sia <20
      //if(RMS_SL_ORD[i]<20.){ 
      //TTrig_mean[i]=TTRIG_SL_ORD[i];
	//TTrig_RMS[i]=RMS_SL_ORD[i];
      //}
      
	  TTrig_mean[i]=TTRIG_SL_ORD[i];
	  TTrig_RMS[i]=RMS_SL_ORD[i];
    }
    
    // TBOX_SL->Update();
  }
  
  
  if(DEBUG_TTRIG)
    cout << "TTrigCalibration::updateCanvas() : END" << endl;
  
  return;
  
} 

void TTrigCalibration::dumpHistos(char* fileName){
  
  TFile * file = new TFile(fileName,"RECREATE");
  file->cd();
  
  //TDC plots    
  for(int i=0;i<8;i++) htbox_ROB[i]->Write();	
  for(int i=0;i<3;i++) htbox_SL[i]->Write();

  file->Close();
  delete file;
  
  //print desidered canvas
  TString nameOut("htbox_SL_");
  nameOut += fileName;
  nameOut += ".eps";
  
  TBOX_SL->Print(nameOut,"eps");
  
  return;
}


void TTrigCalibration::saveTTrigFile(int runN){
  
  for(int i=0;i<6;i++)
    cout<<" ROB n."<<N_ROB[i]<<"  TTRIG (ns) = "<<TTrig_ROB_mean[i]<<"  RMS (ns) = "<<TTrig_ROB_RMS[i]<<endl;
  //int CH[8]={8,9,10,10,10,11,11,11};
  //int SL[8]={1,3,1,2,3,1,2,3};
  for(int i=0;i<3;i++)
    cout << "SL n. " << i+1 << "	TTRIG (ns) = " << TTrig_mean[i] << "		RMS (ns) = " << TTrig_RMS[i] << endl;
    //cout<<" CH n."<<CH[i]<<" SL n."<<SL[i]<<"  TTRIG (ns) = "<<TTrig_mean[i]<<"  RMS (ns) = "<<TTrig_RMS[i]<<endl;
  
  
  FILE *fttrig;
  char fname[200];
  sprintf(fname,"./ttrig/ttrig_%d.txt",runN);
  fttrig=fopen(fname,"w");
  if(DEBUG_TTRIG)
    cout << "Opening file " << fname << endl;
  for(int i=0;i<3;i++)
    fprintf(fttrig,"%4d %4d %10.2f  %10.2f\n",/*CH[i],SL[i]*/11, i+1, TTrig_mean[i],TTrig_RMS[i] );
  fclose(fttrig);
  
  
  return;
}

void TTrigCalibration::dumpTTrigs(char * calibFileName){
  
  cout << "-----------  TTrigCalibration::dumpTTrigs   ----------" << endl;

  for(int i=0;i<6;i++)
    cout<<" ROB n."<<N_ROB[i]<<"  TTRIG (ns) = "<<TTrig_ROB_mean[i]<<"  RMS (ns) = "<<TTrig_ROB_RMS[i]<<endl;
  //int CH[8]={8,9,10,10,10,11,11,11};
  //int SL[8]={1,3,1,2,3,1,2,3};
  for(int i=0;i<3;i++)
    cout << "SL n. " << i+1 << "	TTRIG (ns) = " << TTrig_mean[i] << "		RMS (ns) = " << TTrig_RMS[i] << endl;
    //cout<<" CH n."<<CH[i]<<" SL n."<<SL[i]<<"  TTRIG (ns) = "<<TTrig_mean[i]<<"  RMS (ns) = "<<TTrig_RMS[i]<<endl;
  
  
  FILE *fttrig;
  fttrig=fopen(calibFileName,"w");
  if(DEBUG_TTRIG)
    cout << "Opening file " << calibFileName << endl;
  for(int i=0;i<3;i++)
    fprintf(fttrig,"%4d %4d %10.2f  %10.2f\n",/*CH[i],SL[i],*/11, i+1, TTrig_mean[i],TTrig_RMS[i] );
  fclose(fttrig);
  
  
  return;
}




void TTrigCalibration::fitTimeBox(TH1F *hTimeBox, double& mean, double& sigma) {
  
  mean=999;
  sigma=999;
  
  // Get seeds for the fit
  // The TimeBox range to be fitted (the rising edge)
  double xFitMin=0;     // Min value for the fit
  double xFitMax=0;     // Max value for the fit 
  double xValue=0;      // The mean value of the gaussian
  double xFitSigma=0;   // The sigma of the gaussian
  double tBoxMax=0;     // The max of the time box, it is used as seed for gaussian integral
  
  TH1F *hTimeBoxForSeed = (TH1F*) hTimeBox->Clone(); //FIXME: test
  
  getFitSeeds(hTimeBoxForSeed, xValue, xFitSigma, tBoxMax, xFitMin, xFitMax);
  
  // Define the fitting function and use fit seeds

  TF1 *fIntGaus = new TF1("IntGauss", intGauss, xFitMin, xFitMax, 3); 

  fIntGaus->SetParName(0, "Constant");
  fIntGaus->SetParameter(0, tBoxMax);
  fIntGaus->SetParName(1, "Mean");
  fIntGaus->SetParameter(1, xValue);
  fIntGaus->SetParName(2, "Sigma");
  fIntGaus->SetParameter(2, xFitSigma);
  fIntGaus->SetLineColor(kRed);
  
  // Fit the histo
  char *option = "Q";
  //   if(theVerbosityLevel >= 2)
  //     option = "";
  
  hTimeBox->Fit("IntGauss", option, "",xFitMin, xFitMax);
  
  // Get fitted parameters
  mean =  fIntGaus->GetParameter("Mean");
  sigma = fIntGaus->GetParameter("Sigma");
  //   double constant = fIntGaus->GetParameter("Constant");
  double chiSquare = fIntGaus->GetChisquare()/fIntGaus->GetNDF();
  
  if(DEBUG_TTRIG) {
    cout << " === Fit Results: " << endl;
    cout << "     Fitted mean = " << mean << endl;
    cout << "     Fitted sigma = " << sigma << endl;
    cout << "     Reduced Chi Square = " << chiSquare << endl;
  }
  //return make_pair(mean, sigma);
}


void TTrigCalibration::getFitSeeds(TH1F *hTBox, double& mean, double& sigma, double& tBoxMax,
		 double& xFitMin, double& xFitMax) {
  if(DEBUG_TTRIG)
    cout << " === Looking for fit seeds in Time Box:" << endl;
  
  
  // The approximate width of the time box
  static const int tBoxWidth = 400; //FIXME: tune it
  
  int nBins = hTBox->GetNbinsX();
  const int xMin = (int)hTBox->GetXaxis()->GetXmin();
  const int xMax = (int)hTBox->GetXaxis()->GetXmax();
  const int nEntries =  (int)hTBox->GetEntries();
  
  double binValue = (double)(xMax-xMin)/(double)nBins;
  
  // Compute a threshold for TimeBox discrimination
  const double threshold = binValue*nEntries/(double)(tBoxWidth*2.);
  if(DEBUG_TTRIG)
    cout << "   Threshold for logic time box is (# entries): " <<  threshold << endl;
  
  
  while(threshold > hTBox->GetMaximum()/2.) {
    cout << " Rebinning!" << endl;
    hTBox->Rebin(2);
    nBins = hTBox->GetNbinsX();
    binValue = (double)(xMax-xMin)/(double)nBins;
  }
  
  hDebugFile->cd();
  TString hLName = TString(hTBox->GetName())+"L";
  TH1F hLTB(hLName.Data(), "Logic Time Box", nBins, xMin, xMax);
  // Loop over all time box bins and discriminate them accordigly to the threshold
  for(int i = 1; i <= nBins; i++) {
    if(hTBox->GetBinContent(i) > threshold)
      hLTB.SetBinContent(i, 1);
  }
  hLTB.Write();
  
  // Look for the time box in the "logic histo" and save beginning and lenght of each plateau
  vector< pair<int, int> > startAndLenght;
  if(DEBUG_TTRIG)
    cout << "   Look for rising and folling edges of logic time box: " << endl;
  int start = -1;
  int lenght = -1;
  for(int j = 1; j < nBins;j++) {
    int diff = (int)hLTB.GetBinContent(j+1)-(int)hLTB.GetBinContent(j);
    if(diff == 1) { // This is a rising edge
      start = j;
      lenght = 1;
      if(DEBUG_TTRIGBIN) {
	cout << "     >>>----" << endl;
	cout << "      bin: " << j << " is a rising edge" << endl;
      }
    } else if(diff == -1) { // This is a falling edge
      startAndLenght.push_back(make_pair(start, lenght));
      if(DEBUG_TTRIGBIN) {
	cout << "      bin: " << j << " is a falling edge, lenght is: " << lenght << endl;
	cout << "     <<<----" << endl;
      }
      start = -1;
      lenght = -1;
    } else if(diff == 0 && (int)hLTB.GetBinContent(j) == 1) { // This is a bin within the plateau
      lenght ++;
    }
  }
  
  // Look for the plateau of the right lenght
  if(DEBUG_TTRIG)
    cout << "    Look for the best interval:" << endl;
  int delta = 999999;
  int beginning = -1;
  int tbWidth = -1;
  for(vector< pair<int, int> >::const_iterator stAndL = startAndLenght.begin();
      stAndL != startAndLenght.end();
      stAndL++) {
    if(abs((*stAndL).second - tBoxWidth) < delta) {
      delta = abs((*stAndL).second - tBoxWidth);
      beginning = (*stAndL).first;
      tbWidth = (*stAndL).second;
      if(DEBUG_TTRIG)
	cout << "   Candidate: First: " <<  beginning
	     << ", width: " << tbWidth
	     << ", delta: " << delta << endl;
    }
  }
  
  mean = xMin + beginning*binValue;
  sigma = 10; //FIXME: estimate it!
  
  tBoxMax = hTBox->GetMaximum();
  
  // Define the fit range
  xFitMin = mean-5.*sigma;
  xFitMax = mean+5.*sigma;
  
  if(DEBUG_TTRIG) {
    cout << "      Time Box Rising edge: " << mean << endl;
    cout << "      Time Box Width: " << tbWidth*binValue << endl;
    cout << "    = Seeds and range for fit:" << endl;
    cout << "       Seed mean = " << mean << endl;
    cout << "       Seed sigma = " << sigma << endl;
    cout << "       Fitting from = " << xFitMin << " to " << xFitMax << endl << endl;
  }
}


double intGauss(double *x, double *par) {
  double media = par[1];
  double sigma = par[2];
  double var = (x[0]-media)/(sigma*sqrt(2.));
  
  return 0.5*par[0]*(1+TMath::Erf(var));
  
}
