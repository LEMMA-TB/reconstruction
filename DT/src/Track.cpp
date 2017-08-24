#include "Track.h"

Track::Track() { 

  if(DEBUG_TRACK)
      cout << "Track constructor" << endl;

  for(int i=0;i<6;i++) {
    seg[i]=NULL;
    hit_seg[i]=NULL;
  }
}

Track::~Track()
{

  for(int i=0;i<6;i++) {
    if(seg[i]) delete seg[i];
    if(hit_seg[i]) delete hit_seg[i];
  }

  return;
}
 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 

HITColl_Seg * Track::SelectSeg(HITCollection *hits, int i, int CH, bool phi){
  
  hit_seg[i]->selectHIT(hits,CH,phi);
  
  return hit_seg[i];
}

void Track::SelectTrack(HITCollection *hits, TimeCorr *corr, bool n2chambers)
{
  if(DEBUG_TRACK)
        cout << "Track::SelectTrack" << endl;

  isGood=false;
  
  if(DEBUG_TRACK && hits->Get_NumberHITS()>0)
    printf("\n...Event %d \nInto Track::SelectTrack...\n",(hits->hit(0))->EvNumber());
  
  if( hits->Get_NumberHITS()>MaxNHit ){
    if(DEBUG_TRACK) printf("===> Track with more than %d HITS, track not analyzed!!\n", MaxNHit); 
    isGood=false;
  }
  else
    isGood=true;
  
  // TODO mettere controllo anche sulle hit max dei singoli segmenti!!!
  
  int nseg=6;
  int nseg_glo=4;

  // 20170720 for LEMMA only 1 chamber
  if(!n2chambers)
      nseg_glo=0;

  float w_xcorr[nseg]; 
  for(int i=0;i<nseg;i++) w_xcorr[i]=0.;
  
  int min_pt[6]={MinNHit_Phi,MinNHit_Theta,MinNHit_Phi,MinNHit_Theta,MinNHit_Theta,MinNHit_Theta};
  int max_pt[6]={MaxNHit_Phi,MaxNHit_Theta,MaxNHit_Phi,MaxNHit_Theta,MaxNHit_Theta,MaxNHit_Theta};
  int min_pt_abs[6]={2,1,2,1,1,1};
  
  for(int i=0;i<nseg;i++){
    hit_seg[i]=NULL;
    hit_seg[i]=new HITColl_Seg();
  }
  
  hit_seg[0]=SelectSeg(hits,0,11,true);
  hit_seg[1]=SelectSeg(hits,1,11,false);
  hit_seg[2]=SelectSeg(hits,2,10,true); //inutile
  hit_seg[3]=SelectSeg(hits,3,10,false); //inutile
  hit_seg[4]=SelectSeg(hits,4,9,true);  //inutile
  hit_seg[5]=SelectSeg(hits,5,8,true);  //inutile
  
  w_xcorr[0]=Get_MeanWire(hit_seg[1]);
  w_xcorr[1]=Get_MeanWire(hit_seg[0]);
  w_xcorr[2]=Get_MeanWire(hit_seg[3]);
  w_xcorr[3]=Get_MeanWire(hit_seg[2]);
  w_xcorr[4]=w_xcorr[2];
  w_xcorr[5]=w_xcorr[2];
  
  int nrVarv_ch1[4]={3,3,0,0};
  int nrVarv_ch2[4]={0,0,3,3};
  int nrVarv[4];
  double Slope[nseg],X0[nseg],T0[nseg],Chi2[nseg],T0_fin[nseg];
  double Slope_glo[nseg_glo],X0_glo[nseg_glo],T0_glo[nseg_glo],T0_glo_fin[nseg_glo],Chi2_glo[nseg_glo];
  double erSlope_glo[nseg_glo],erX0_glo[nseg_glo],erT0_glo[nseg_glo];
  int NPT[nseg];
  int NPT_glo[nseg_glo];
  bool ISGood[nseg];
  bool ISGood_glo[nseg_glo];
  MP_ch1_isok=false;
  MP_ch2_isok=false;

  for(int i=0;i<nseg;i++)
    {
      //ISGood[i]=true; //ALTEA  
      ISGood[i]=false; 
      Slope[i]=-999.; X0[i]=-999.; T0[i]=-999.; NPT[i]=-999; Chi2[i]=-999.; T0_fin[i]=-999.;
      SegIsAbsorbed[i]=false;

      Set_IsGood(i,ISGood[i]); 
      Set_Slope(i,Slope[i]);
      Set_X0(i,X0[i]);
      Set_T0(i,T0[i]);
      Set_T0_fin(i,T0_fin[i]);
      Set_Chi2(i,Chi2[i]);
      Set_NPT(i,NPT[i]);
    }
  for(int i=0;i<nseg_glo;i++)
    {
      ISGood_glo[i]=false; Slope_glo[i]=-999.; X0_glo[i]=-999.; T0_glo[i]=-999.; T0_glo_fin[i]=-999.; NPT_glo[i]=-999; Chi2_glo[i]=-999.;
      erSlope_glo[i]=-999.; erX0_glo[i]=-999.; erT0_glo[i]=-999.;

      Set_IsGood_glo(i,ISGood_glo[i]); 
      Set_Slope_glo(i,Slope_glo[i]);
      Set_X0_glo(i,X0_glo[i]);
      Set_T0_glo(i,T0_glo[i]);
      Set_T0_glo_fin(i,T0_glo_fin[i]);
      Set_Chi2_glo(i,Chi2_glo[i]);
      Set_NPT_glo(i,NPT_glo[i]);
      Set_erSlope_glo(i,erSlope_glo[i]);
      Set_erX0_glo(i,erX0_glo[i]);
      Set_erT0_glo(i,erT0_glo[i]);
    }
  
  double T0phi1=-999.;
  double T0phi2=-999.;
  
  for(int i=0;i<nseg;i++){
    seg[i]=NULL;
    seg[i]=new Segment();  
  }
  
  if(isGood){
    
    for(int i=0;i<nseg;i++){
      if( hit_seg[i]->Get_NumberHITS()<=min_pt_abs[i]){
          SegIsAbsorbed[i]=true;
          if(DEBUG_TRACK)
              printf("Event %d: Seg %d is absorbed, N. hits = %d <= %d\n",(hits->hit(0))->EvNumber(),i,hit_seg[i]->Get_NumberHITS(),min_pt_abs[i]);
      }
    }
    
    for(int i=0;i<2;i++){
      if( hit_seg[i]->Get_NumberHITS()<min_pt[i] ||  hit_seg[i]->Get_NumberHITS()>max_pt[i]){
          if(DEBUG_TRACK)
              printf("1: Seg.%d is not good ---> track rejected!\n",i);
          isGood=false;
          return;
      }
    }
    
    for(int i=0;i<nseg;i++){
      hit_seg_glo[i]=NULL;
    }
    
    for(int i=0;i<nseg;i++){
      if( hit_seg[i]->Get_NumberHITS()>=min_pt[i] &&  hit_seg[i]->Get_NumberHITS()<=max_pt[i]){
          if(i==0){
              hit_seg_glo[i]=seg[i]->Get_PHISegment(hit_seg[i],w_xcorr[i],corr);
            if(seg[i]->IsGood()){ //ALTEA
                 if( TMath::Abs(seg[i]->Get_T0()) < T0_max && TMath::Abs(seg[i]->Get_T0_fin()) < T0_max)
                     T0phi1=(seg[i]->Get_T0()+seg[i]->Get_T0_fin()) - (19.79 * pow(seg[i]->Get_Slope(),2));
                 else if(TMath::Abs(seg[i]->Get_T0_fin()) < T0_max)
                     T0phi1=(seg[i]->Get_T0_fin()) - (19.79 * pow(seg[i]->Get_Slope(),2));
                 else if(TMath::Abs(seg[i]->Get_T0()) < T0_max)
                     T0phi1=(seg[i]->Get_T0()) - (19.79 * pow(seg[i]->Get_Slope(),2));
                 else
                     T0phi1=0;
 	          } 
            else{ //ALTEA
 	             if(DEBUG_TRACK)
 	                 printf("2: Seg.%d is not good ---> track rejected!\n",i);
 	             isGood=false;
 	             return;
 	          } 
	      }
	      else if(i==1){
              hit_seg_glo[i]=seg[i]->Get_SLSegment(hit_seg[i],w_xcorr[i],corr,T0phi1);
              if(!seg[i]->IsGood()){
                if(DEBUG_TRACK) //ALTEA
                    printf("3: Seg.%d is not good ---> track rejected!\n",i);
                isGood=false;
                return;
              }
	     }
// 	     ALTEA2
// 	     else if(i==2){
//              cout << "Non ci sono mai" << endl;
//              hit_seg_glo[i]=seg[i]->Get_PHISegment(hit_seg[i],w_xcorr[i],corr);
//              if(seg[i]->IsGood())
//                  if( TMath::Abs(seg[i]->Get_T0()) < T0_max && TMath::Abs(seg[i]->Get_T0_fin()) < T0_max )
//                      T0phi2=(seg[i]->Get_T0()+seg[i]->Get_T0_fin()) - (19.79 * pow(seg[i]->Get_Slope(),2));
//                  else if( TMath::Abs(seg[i]->Get_T0()) < T0_max )
//                      T0phi2=(seg[i]->Get_T0()) - (19.79 * pow(seg[i]->Get_Slope(),2));
//                  else if( TMath::Abs(seg[i]->Get_T0_fin()) < T0_max )
//                      T0phi2=(seg[i]->Get_T0_fin()) - (19.79 * pow(seg[i]->Get_Slope(),2));
//                  else
//                      T0phi2=T0phi1; // sistemare e mettere tof corretto!!!
//              else
//                 T0phi2=T0phi1; // sistemare e mettere tof corretto!!!
// 	     }
	     else{
             hit_seg_glo[i]=seg[i]->Get_SLSegment(hit_seg[i],w_xcorr[i],corr,T0phi2);	
	     }
	
	
    if(seg[i]->IsGood()){ //ALTEA
 	  ISGood[i]=seg[i]->IsGood();
	  Slope[i]=seg[i]->Get_Slope();
	  X0[i]=seg[i]->Get_X0();
	  T0[i]=seg[i]->Get_T0();
	  T0_fin[i]=seg[i]->Get_T0_fin();
	  Chi2[i]=seg[i]->Get_Chi2();
	  NPT[i]=seg[i]->Get_NPT();
	  
	  if(i==0 || i==2){
	    if( TMath::Abs(seg[i]->Get_T0()) < T0_max && TMath::Abs(seg[i]->Get_T0_fin()) < T0_max )
	      T0_glo_fin[i]=seg[i]->Get_T0() + seg[i]->Get_T0_fin();
	    else if( TMath::Abs(seg[i]->Get_T0()) < T0_max )
	      T0_glo_fin[i]=seg[i]->Get_T0();
	    else if( TMath::Abs(seg[i]->Get_T0_fin()) < T0_max )
	      T0_glo_fin[i]=seg[i]->Get_T0_fin();
	  }
 	}
	
	if(DEBUG_TRACK) {
	  printf("Segment %d: \n  IsDouble: %d, IsGood: %d, OKFIT: %d\n",i,seg[i]->IsDouble(),seg[i]->IsGood(),seg[i]->Get_FitOK());
	  printf("Alla fine di tutto: slope %.3f, X0 %.1f, T0 %.1f, T0_fin %.1f, NPT %d, Chi2 %.3f\n",seg[i]->Get_Slope(),seg[i]->Get_X0(),seg[i]->Get_T0(),seg[i]->Get_T0_fin(),seg[i]->Get_NPT(),seg[i]->Get_Chi2());
	}
	
	if(DEBUG_TRACK) 
	  printf("T0phi1 = %.1f, T0phi2 = %.1f\n",T0phi1,T0phi2);
	
      }
      Set_IsGood(i,ISGood[i]); 
      Set_Slope(i,Slope[i]);
      Set_X0(i,X0[i]);
      Set_T0(i,T0[i]);
      Set_T0_fin(i,T0_fin[i]);
      Set_Chi2(i,Chi2[i]);
      Set_NPT(i,NPT[i]);
    }
    
    if(DEBUG_TRACK)
      for(int i=0;i<nseg_glo;i++) 
	cout<<"seg i="<<i<<" Get_IsGood "<<Get_IsGood(i)<<endl;
    
    //     if( Get_IsGood(0) && Get_IsGood(1) && Get_IsGood(2) && Get_IsGood(3) ){
    if( Get_IsGood(0) && Get_IsGood(1) ){
      for(int i=0;i<nseg_glo;i++) 
	if(DEBUG_TRACK && Get_IsGood(i)) 
	  printf("hit_seg_glo[%d] has %d hit\n",i,hit_seg_glo[i]->Get_NumberHITS());
      
      double sigmaPhi[4]={0.,0.,0.,0.}; 
      int nrMinPnts[4]={0,0,0,0}; 
      vector<double> axf_sl[4],ay_sl[4],at_sl[4];
      vector<int> ac_sl[4];
      vector<vector<double> > axf,ay,at;
      vector<vector<int> > ac;
      int nrPnts[4]={0,0,0,0};

      double m[4]={-999.,-999.,-999.,-999.};
      double q[4]={-999.,-999.,-999.,-999.};
      double t0=0.;
      double sm[4]={-999.,-999.,-999.,-999.};
      double sq[4]={-999.,-999.,-999.,-999.};
      double st0=0.;
      
      double vdrift=0;
      double chi2[4]={-999.,-999.,-999.,-999.};
      int NPT[4]={0,0,0,0};
      bool okFit=false;
// ALTEA2
//       for(int i=0;i<nseg_glo;i++)
// 	if(Get_IsGood(i))
// 	  {          
// 	    sigmaPhi[i]=corr->Get_SigmaAngCorr(seg[i]->Get_Slope());
// 	    nrPnts[i]=seg[i]->Get_NPT();
// 	    nrMinPnts[i]=min_pt[i];
// 	    for(int j=0;j<nrPnts[i];j++){
// 	      HIT *hit=hit_seg_glo[i]->hit(j);
// 	      axf_sl[i].push_back(hit->x_wire_ID());
// 	      ay_sl[i].push_back(hit->y_wire_ID());
// 	      at_sl[i].push_back(hit->dtime());
// 	      if(hit->code_ID().first)
// 		ac_sl[i].push_back(-1);
// 	      else
// 		ac_sl[i].push_back(1);
// 	    }
// 	    axf.push_back(axf_sl[i]);
// 	    ay.push_back(ay_sl[i]);
// 	    ac.push_back(ac_sl[i]);
// 	    at.push_back(at_sl[i]);
// 	  }
// 	else{
// 	  for(int j=0;j<nrPnts[i];j++){
// 	    axf_sl[i].push_back(-999.);
// 	    ay_sl[i].push_back(-999.);
// 	    at_sl[i].push_back(-999.);
// 	    ac_sl[i].push_back(-999);
// 	  }
// 	  axf.push_back(axf_sl[i]);
// 	  ay.push_back(ay_sl[i]);
// 	  ac.push_back(ac_sl[i]);
// 	  at.push_back(at_sl[i]);
// 	}
      
      
      for(int i=0;i<nseg_glo;i++) 
	nrVarv[i]=nrVarv_ch1[i];
      
    int fit_CH=1;
    fitg->FIT_global(sigmaPhi, nrPnts, nrMinPnts, nrVarv, axf, at, ay, ac, m, sm, q, sq, t0, st0, vdrift, chi2, NPT, okFit, fit_CH);  
      
      for(int i=0;i<2;i++){
	if(T0_glo_fin[i] != -999.)
	  T0_glo_fin[i] += t0;
	else
	  T0_glo_fin[i] = t0;
	T0_glo[i] = t0;
	Slope_glo[i]=m[i];
	X0_glo[i]=q[i];
	Chi2_glo[i]=chi2[i];
	NPT_glo[i]=NPT[i];
	ISGood_glo[i]=okFit;
	erT0_glo[i]=st0;
	erSlope_glo[i]=sm[i];
	erX0_glo[i]=sq[i];
      }
      // ***********************************************

      
      if( nrPnts[1]==4 ){
	double xth[4]={-999.,-999.,-999.,-999.};
	for(int pi=0;pi<4;pi++){
	  int pii=-999;
	  if(ay[1][pi] >= 7. && ay[1][pi] <= 8) pii=0;
	  if(ay[1][pi] >= 8. && ay[1][pi] <= 9) pii=1;
	  if(ay[1][pi] >=10. && ay[1][pi] <=11) pii=2;
	  if(ay[1][pi] >=11. && ay[1][pi] <=12) pii=3;
	  
	  xth[pii] = axf[1][pi]+ac[1][pi]*vdrift*at[1][pi] - (ac[1][pi]*vdrift*t0);
	  
	  if(DEBUG_TRACK)
	    printf("xth[%d]=%.1f, ay[%d]=%.3f\n",pii,xth[pii],pi,ay[1][pi]);
	}
	
	float MX_CH1_123 = ((xth[0]+xth[2])/2. - xth[1])/(sigmaPhi[1]*1.22);
	float MX_CH1_234 = ((xth[1]+xth[3])/2. - xth[2])/(sigmaPhi[1]*1.22);
	float MX_CH1_134 = ((xth[0]+2.*xth[3])/3. - xth[2])/(sigmaPhi[1]*1.247);
	float MX_CH1_124 = ((2.*xth[0]+xth[3])/3. - xth[1])/(sigmaPhi[1]*1.247);
	//printf("MX_CH1_123=%.2f, MX_CH1_234=%.2f, MX_CH1_134=%.2f, MX_CH1_124=%.2f\n",MX_CH1_123,MX_CH1_234,MX_CH1_134,MX_CH1_124);
	Double_t cut=2.5;
	
	if( MX_CH1_234>-cut && MX_CH1_234<cut && 
	    MX_CH1_123>-cut && MX_CH1_123<cut &&
	    MX_CH1_134>-cut && MX_CH1_134<cut &&
	    MX_CH1_124>-cut && MX_CH1_124<cut ) 
	  MP_ch1_isok=true;
	
      }
      
      // ***********************************************
      
      t0=0.;
      
//     ALTEA 2
//       if(Get_IsGood(2) && Get_IsGood(3)){
// 	for(int i=0;i<nseg_glo;i++) 
// 	  nrVarv[i]=nrVarv_ch2[i];
// 	fit_CH=2;
// 	fitg->FIT_global(sigmaPhi, nrPnts, nrMinPnts, nrVarv, axf, at, ay, ac, m, sm, q, sq, t0, st0, vdrift, chi2, NPT, okFit, fit_CH);  
// 	
// 	for(int i=2;i<4;i++){
// 	  if(T0_glo_fin[i] != -999.)
// 	    T0_glo_fin[i] += t0;
// 	  else
// 	    T0_glo_fin[i] = t0;
// 	  T0_glo[i] = t0;
// 	  Slope_glo[i]=m[i];
// 	  X0_glo[i]=q[i];
// 	  Chi2_glo[i]=chi2[i];
// 	  NPT_glo[i]=NPT[i];
// 	  ISGood_glo[i]=okFit;
// 	  erT0_glo[i]=st0;
// 	  erSlope_glo[i]=sm[i];
// 	  erX0_glo[i]=sq[i];
// 	}
// 	
// 	if( nrPnts[3]==4 ){
// 	  double xth[4]={-999.,-999.,-999.,-999.};
// 	  for(int pi=0;pi<4;pi++){
// 	    int pii=-999;
// 	    if(ay[3][pi] >=-12. && ay[3][pi] <=-11) pii=0;
// 	    if(ay[3][pi] >=-11. && ay[3][pi] <=-10.) pii=1;
// 	    if(ay[3][pi] >= -9. && ay[3][pi] <= -8.) pii=2;
// 	    if(ay[3][pi] >= -8. && ay[3][pi] <= -7.) pii=3;
// 	    xth[pii] = axf[3][pi]+ac[3][pi]*vdrift*at[3][pi] - (ac[3][pi]*vdrift*t0);
// 	    
// 	    if(DEBUG_TRACK)
// 	      printf("x2th[%d]=%.1f, ay[%d]=%.3f\n",pii,xth[pii],pi,ay[1][pi]);
// 	  }
// 	  
// 	  float MX_CH2_123 = ((xth[0]+xth[2])/2. - xth[1])/(sigmaPhi[3]*1.22);
// 	  float MX_CH2_234 = ((xth[1]+xth[3])/2. - xth[2])/(sigmaPhi[3]*1.22);
// 	  float MX_CH2_134 = ((xth[0]+2.*xth[3])/3. - xth[2])/(sigmaPhi[3]*1.247);
// 	  float MX_CH2_124 = ((2.*xth[0]+xth[3])/3. - xth[1])/(sigmaPhi[3]*1.247);
// 	  //printf("MX_CH2_123=%.2f, MX_CH2_234=%.2f, MX_CH2_134=%.2f, MX_CH2_124=%.2f\n",MX_CH2_123,MX_CH2_234,MX_CH2_134,MX_CH2_124);
// 	  Double_t cut=2.5;
// 	  
// 	  if( 
// 	     MX_CH2_234>-cut && MX_CH2_234<cut && 
// 	     MX_CH2_123>-cut && MX_CH2_123<cut && 
// 	     MX_CH2_134>-cut && MX_CH2_134<cut &&
// 	     MX_CH2_124>-cut && MX_CH2_124<cut ) 
// 	    MP_ch2_isok=true;
// 	  
// 	}
//       }
      
      if(DEBUG_TRACK){
	for(int i=0;i<nseg_glo;i++) 
	  printf("Alla fine del fit globale su seg %d: \nGet_IsGood=%d, NPT=%d, m=%.3f, q=%.0f, \nchi2=%.2f, t0=%.1f, t0_tot=%.1f\n",i,ISGood_glo[i],NPT_glo[i],Slope_glo[i],X0_glo[i],Chi2_glo[i],T0_glo[i],T0_glo_fin[i]);  
	printf("MP_CH1=%d, MP_CH2=%d\n",MP_ch1_isok,MP_ch2_isok);
      }
      
      for(int i=0;i<nseg_glo;i++){
	
 	Set_Slope_glo(i,Slope_glo[i]);
	Set_X0_glo(i,X0_glo[i]);
	Set_T0_glo(i,T0_glo[i]);
	Set_T0_glo_fin(i,T0_glo_fin[i]);
	Set_Chi2_glo(i,Chi2_glo[i]);
 	Set_NPT_glo(i,NPT_glo[i]);
	Set_IsGood_glo(i,ISGood_glo[i]);
 	Set_erSlope_glo(i,erSlope_glo[i]);
	Set_erX0_glo(i,erX0_glo[i]);
	Set_erT0_glo(i,erT0_glo[i]);
      }

    }
  } // close if(isGood)
  
  return;
}

 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 

float Track::Get_MeanWire(HITColl_Seg *hit_seg){
  float N=0.; int j=0;
  if(hit_seg->Get_NumberHITS()>0){
    for(int i=0;i<hit_seg->Get_NumberHITS();i++){
      HIT *hit=hit_seg->hit(i);
      N += hit->wire_ID();
      j++; 
    }
    N /= float(j);
  }
  
  if(DEBUG_TRACK)
    printf("Mean Wire for correction %.1f\n",N);

  return (N);
}
 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


bool Track::Track_IsGood(){
  return isGood;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool Track::Set_IsGood(int s, bool isgood){
  track_IsGood[s]=isgood;
  return track_IsGood[s];
}

bool Track::Get_IsGood(int s){
  return track_IsGood[s];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double Track::Set_Slope(int s, double slope){
  track_Slope[s]=slope;
  return track_Slope[s];
}

double Track::Get_Slope(int s){
  return track_Slope[s];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double Track::Set_X0(int s, double X0){
  track_X0[s]=X0;
  return track_X0[s];
}

double Track::Get_X0(int s){
  return track_X0[s];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double Track::Set_T0(int s, double T0){
  track_T0[s]=T0;
  return track_T0[s];
}

double Track::Get_T0(int s){
  return track_T0[s];
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double Track::Set_T0_fin(int s, double T0_fin){
  track_T0_fin[s]=T0_fin;
  return track_T0_fin[s];
}

double Track::Get_T0_fin(int s){
  return track_T0_fin[s];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double Track::Set_Chi2(int s, double Chi2){
  track_Chi2[s]=Chi2;
  return track_Chi2[s];
}

double Track::Get_Chi2(int s){
  return track_Chi2[s];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int Track::Set_NPT(int s, int NPT){
  track_NPT[s]=NPT;
  return track_NPT[s];
}

int Track::Get_NPT(int s){
  return track_NPT[s];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double Track::Set_Slope_glo(int s, double slope_glo){
  track_Slope_glo[s]=slope_glo;
  return track_Slope_glo[s];
}

double Track::Get_Slope_glo(int s){
  return track_Slope_glo[s];
}

double Track::Set_erSlope_glo(int s, double erslope_glo){
  track_erSlope_glo[s]=erslope_glo;
  return track_erSlope_glo[s];
}

double Track::Get_erSlope_glo(int s){
  return track_erSlope_glo[s];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double Track::Set_X0_glo(int s, double X0_glo){
  track_X0_glo[s]=X0_glo;
  return track_X0_glo[s];
}

double Track::Get_X0_glo(int s){
  return track_X0_glo[s];
}

double Track::Set_erX0_glo(int s, double erX0_glo){
  track_erX0_glo[s]=erX0_glo;
  return track_erX0_glo[s];
}

double Track::Get_erX0_glo(int s){
  return track_erX0_glo[s];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double Track::Set_T0_glo(int s, double T0_glo){
  track_T0_glo[s]=T0_glo;
  return track_T0_glo[s];
}

double Track::Get_T0_glo(int s){
  return track_T0_glo[s];
}

double Track::Set_erT0_glo(int s, double erT0_glo){
  track_erT0_glo[s]=erT0_glo;
  return track_erT0_glo[s];
}

double Track::Get_erT0_glo(int s){
  return track_erT0_glo[s];
}

double Track::Set_T0_glo_fin(int s, double T0_glo_fin){
  track_T0_glo_fin[s]=T0_glo_fin;
  return track_T0_glo_fin[s];
}

double Track::Get_T0_glo_fin(int s){
  return track_T0_glo_fin[s];
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double Track::Set_Chi2_glo(int s, double Chi2_glo){
  track_Chi2_glo[s]=Chi2_glo;
  return track_Chi2_glo[s];
}

double Track::Get_Chi2_glo(int s){
  return track_Chi2_glo[s];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int Track::Set_NPT_glo(int s, int NPT_glo){
  track_NPT_glo[s]=NPT_glo;
  return track_NPT_glo[s];
}

int Track::Get_NPT_glo(int s){
  return track_NPT_glo[s];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool Track::Set_IsGood_glo(int s, bool isgood_glo){
  track_IsGood_glo[s]=isgood_glo;
  return track_IsGood_glo[s];
}

bool Track::Get_IsGood_glo(int s){
  return track_IsGood_glo[s];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool Track::Get_MP_CH1_IsOk_glo(){
  return MP_ch1_isok;
}

bool Track::Get_MP_CH2_IsOk_glo(){
  return MP_ch2_isok;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool Track::Get_SegIsAbsorbed(int s){
  return SegIsAbsorbed[s];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Track::printTrack()
{
  
}
