#include "TimeCorr.h"

bool debug_TimeCorr  = false;

TimeCorr::~TimeCorr()
{   
  return;
}

//void TimeCorr::Init(){
void TimeCorr::InitSpline(){
  
  fsp = new TFile("utils/splineN8new.root");
  for(int slo=0;slo<8;slo++){
    char fn[20];
    sprintf(fn,"sp%d",slo);
    spline[slo]=(TSpline3*)gDirectory->Get(fn);
    }
  /*
    spline[0]=(TSpline3*)gDirectory->Get("sp0") ;  
    spline[1]=(TSpline3*)gDirectory->Get("sp1") ;  
    spline[2]=(TSpline3*)gDirectory->Get("sp2") ;  
    spline[3]=(TSpline3*)gDirectory->Get("sp3") ;  
    spline[4]=(TSpline3*)gDirectory->Get("sp4") ;  
    spline[5]=(TSpline3*)gDirectory->Get("sp5") ;  
    spline[6]=(TSpline3*)gDirectory->Get("sp6") ;  
    spline[7]=(TSpline3*)gDirectory->Get("sp7") ;  
  */
  //   spline[0]=new TSpline3(*( (TSpline3*)gDirectory->Get("sp0") ));  
  //   spline[1]=new TSpline3(*( (TSpline3*)gDirectory->Get("sp1") ));  
  //   spline[2]=new TSpline3(*( (TSpline3*)gDirectory->Get("sp2") ));  
  //   spline[3]=new TSpline3(*( (TSpline3*)gDirectory->Get("sp3") ));  
  //   spline[4]=new TSpline3(*( (TSpline3*)gDirectory->Get("sp4") ));  
  //   spline[5]=new TSpline3(*( (TSpline3*)gDirectory->Get("sp5") ));  
  //   spline[6]=new TSpline3(*( (TSpline3*)gDirectory->Get("sp6") ));  
  //   spline[7]=new TSpline3(*( (TSpline3*)gDirectory->Get("sp7") ));  
  
  fsp->Close();
  delete fsp;
  fsp=NULL;

  return;  
}

float TimeCorr::Get_Time_SlopeCorr(float time, float slope_raw){
  
//   float dtime_slopecorr = time + convToNs*20.5 * pow(slope_raw,2);
  float dtime_slopecorr = time + 19.79 * pow(slope_raw,2);
  if(debug_TimeCorr)
    printf("Time after slope correction = %.1f\n",dtime_slopecorr);
  
  return dtime_slopecorr;
}

float TimeCorr::Get_Time_DriftCorr(float time, int SL, int mean_wire){
  
  float x_prop;
  if(SL==1 || SL==3)
    x_prop = mean_wire * 4.2;
    //x_prop = (59-mean_wire) * 4.2;
  else
    x_prop = (73 - mean_wire) * 4.2; 
  
  float dtime_driftcorr = time - (x_prop / velWireProp);
  if(debug_TimeCorr)
    printf("Time after drift correction = %.1f\n",dtime_driftcorr);
  
  return dtime_driftcorr;
  
}

float TimeCorr::Get_SigmaAngCorr(float slope_raw){ 
  
  float sigma_angle ;
  if(slope_raw!=0.)
    sigma_angle = 0.03 + 0.025 * (1-TMath::Cos(TMath::ATan(slope_raw)));
  else
    sigma_angle = 0.1;
  if(debug_TimeCorr)
    printf("Sigma angle correction = %.2f\n",sigma_angle);
  
  return sigma_angle;
  
}

float TimeCorr::Get_LinearityCorr(int ch, int sl, float time, float slope_raw){ 

//   int seg=-1;
//   if(ch==11) 
//     if(sl==1 || sl==3)
//       seg=0;
//     else
//       seg=1;
//   if(ch==10) 
//     if(sl==1 || sl==3)
//       seg=2;
//     else
//       seg=3;
//   if(ch==9) seg=4;
//   if(ch==8) seg=5;

//  float m[9]={0.,0.08749,0.1763269,0.267949,0.363970,0.466308,0.577350,0.700207,0.839099};
  float m[9]={0.,0.08749,0.1763269,0.267949,0.363970,0.466308,0.577350,0.700207,3.};
  //  float maxt[8]={381.,384.,392.,402.,402.,402.,397.,392.};
  float maxt=380.;
  
  float dtime_lincorr=0.; 
  //   float tshift=12.5;
  float tshift=0.;

  for(int i=0;i<8;i++){
    if( TMath::Abs(slope_raw)>=m[i] && TMath::Abs(slope_raw)<m[i+1] ){
      if((time+tshift)>2.5 && (time+tshift)<maxt)
	dtime_lincorr = time - (spline[i]->Eval(time+tshift));   
      else if((time+tshift)<=2.5)
	dtime_lincorr = time - (spline[i]->Eval(2.5));
      else if((time+tshift)>=maxt)
	dtime_lincorr = time - (spline[i]->Eval(maxt));
    }
  }
  if(debug_TimeCorr) printf("Time after linear correction %.1f\n",dtime_lincorr); 
  
  return dtime_lincorr;
  
}


