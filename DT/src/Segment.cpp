#include "Segment.h"

bool drift_corr = true;
bool slope_corr = true;
bool t0_corr = true;
bool lin_corr = true;

static const float cut_SLs = 1.5; // in cm
static const float cut_SL = 4.; // in cm

Segment::Segment()
{
  return;
}

Segment::~Segment()
{
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


HITColl_Seg * Segment::Get_PHISegment(HITColl_Seg *hit_seg, float w_xcorr, TimeCorr *corr)
{
  if(DEBUG_SEG) printf("\n ...Start computing PHI SEGMENT...\n");
  
  seg_Slope=-999.; seg_X0=-999.; seg_T0=-999.; seg_T0_fin=-999.; seg_chi2=-999.; seg_NPT=-999; 
  seg_isdouble=false; 
  seg_isgood=true; 
  seg_okFit=false; 
  onlyoneseg=false;
  
  hit_seg = RejectDoubleHits(hit_seg);
  
  if(DEBUG_SEG) 
    printf("n.hit after RejectDoubleHits: %d\n",hit_seg->Get_NumberHITS());
  if(hit_seg->Get_NumberHITS()<MinNHit_Phi){
    Set_IsDouble(false);
    Set_IsGood(false);
    //count;
    return hit_seg;
  }
  
  // drift time correction
  if(drift_corr) Set_HITTime_Drift(hit_seg,w_xcorr,corr);
  
  vector< pair<float,float> > _mean1;
  vector< pair<float,float> > _mean3;
  int size1=0;
  int size3=0;
  _mean1=Xmean_SL(hit_seg, 1,size1);
  if(DEBUG_SEG) 
    printf("   N.hit tot=%d - if double: N.hit1=%d, N.hit2=%d\n",ii,ii1,ii2);
  _mean3=Xmean_SL(hit_seg, 3,size3);
  if(DEBUG_SEG) {
    printf("   N.hit tot=%d - if double: N.hit1=%d, N.hit2=%d\n",ii,ii1,ii2);
    printf("\nSize SL 1 %d, Size SL 3 %d\n",size1,size3);
  }
  
  if(size1==2 && size3==2) {
    Set_IsDouble(true);
    Set_IsGood(false);
    return hit_seg;
  }
  
  if(size1==0 && size3==0) {
    Set_IsDouble(false);
    Set_IsGood(false);
    return hit_seg;
  }
  
  pair<float,float> PHI1=make_pair(0,0);
  pair<float,float> PHI3=make_pair(0,0);
  pair<float,float> seg_par=make_pair(-999.,-999.);
  float slope0=1000000.;
  pair<float,float> seg_par0=make_pair(0,0);
  float xmean_SL=-999.;
  
  if(size1 != 0 && size3 != 0){
    for(int i=0;i<size1;i++)
      for(int j=0;j<size3;j++){
	PHI1=_mean1[i];
	PHI3=_mean3[j];
	seg_par0=Get_SlopePHI(PHI1,PHI3);
	if(fabs(seg_par0.first)<fabs(slope0)){
	  slope0 =seg_par0.first;
	  seg_par=make_pair(seg_par0.first,seg_par0.second);
	}
      }
  }
  else 
    if(size1 != 0 || size3 != 0){
      onlyoneseg=true;
      if(DEBUG_SEG) printf("Only one segment with N.hit %d: size1 = %d, size3 = %d\n",hit_seg->Get_NumberHITS(),size1,size3);
      if(hit_seg->Get_NumberHITS()<MinNHit_1SLPhi)
	  {
        Set_IsDouble(false);
        Set_IsGood(false);
      }
      else if( size1==2 || size3==2 )
	  {
        if(DEBUG_SEG)
            printf("Seg is double, segment rejected\n");
        Set_IsDouble(true);
        Set_IsGood(false);
        }
        else if(hit_seg->Get_NumberHITS()<= MaxNHit_1SLPhi) 
        {
            if(size1!=0) xmean_SL=_mean1[0].first;
            if(size3!=0) xmean_SL=_mean3[0].first;
            if(DEBUG_SEG) printf("xmean_SL = %.1f\n",xmean_SL);
            Set_IsGood(true);
        } 
        else{
            Set_IsDouble(false);
            Set_IsGood(false);
        }
    }
    else{
      Set_IsDouble(false);
      Set_IsGood(false);
    }
  
    if(!IsGood()) 
        return hit_seg;
  
    if(DEBUG_SEG && !onlyoneseg) printf("Slope after selection: %.1f\n",seg_par.first);
    double sigmaSlope=0.;
    double sigmaTime=0.;
    int nrVarv=3;
    bool change_code=false;
  
    if(onlyoneseg==false){
        // slope time correction
        if(slope_corr) Set_HITTime_Slope(hit_seg,seg_par,corr);
        hit_seg = SelectHITsforFIT(hit_seg,seg_par);
        sigmaSlope = corr->Get_SigmaAngCorr(seg_par.first);
        if(DEBUG_SEG) printf("Sigma SlopeCorr: %.2f\n",sigmaSlope);
    }
    else{
        hit_seg = SelectHITsforSLFIT(hit_seg,xmean_SL);      
        if(hit_seg->Get_NumberHITS()<MinNHit_1SLPhi){
            Set_IsGood(false);
            return hit_seg;
        }
    
        nrVarv=2; 
        change_code=false;
        sigmaTime=sigmaTimes_1SL;
    
        SolveAmbiguity(sigmaSlope,sigmaTime,hit_seg,nrVarv,change_code);
    
        if(DEBUG_SEG) 
            printf("Linear Fit: is good = %d, slope = %.3f, X0 = %.1f, t0 = %.1f, NPT = %d, chi2 = %.3f\n",Get_FitOK(),Get_Slope(),Get_X0(),Get_T0(),Get_NPT(),Get_Chi2()); 
    
        if(!Get_FitOK()){
            Set_IsGood(false); //qui
            return hit_seg;
        }
    
        seg_okFit=false; // lo riazzero!
        seg_par = make_pair(Get_Slope(),Get_X0());
        // slope time correction
        if(slope_corr) Set_HITTime_Slope(hit_seg,seg_par,corr);
        hit_seg = SelectHITsforFIT(hit_seg,seg_par);
        sigmaSlope = corr->Get_SigmaAngCorr(seg_par.first);
        if(DEBUG_SEG) printf("Sigma SlopeCorr: %.2f\n",sigmaSlope);
    }
  
    if( (onlyoneseg==false && hit_seg->Get_NumberHITS()<MinNHit_Phi)
        || (onlyoneseg==true && hit_seg->Get_NumberHITS()<MinNHit_1SLPhi) )
    {
      Set_IsGood(false);
      return hit_seg;
    }
  
  if(onlyoneseg==false){
    nrVarv=3; 
    change_code=false;
    sigmaTime=sigmaTimes_Phi;
    
    SolveAmbiguity(sigmaSlope,sigmaTime,hit_seg,nrVarv,change_code);
     
    if(!Get_FitOK()){
       Set_IsGood(false); //qui
       return hit_seg;
     }
    
    seg_okFit=false; // lo riazzero!
    // T0 time correction
    if(t0_corr) Set_HITTime_T0(hit_seg,Get_T0());
  }
  
  // Linearity correction
  if(lin_corr) Set_HITTime_Linearity(hit_seg,seg_par,corr);
  
  if(!onlyoneseg) 
    nrVarv=3; 
  else
    nrVarv=2;
  
  change_code=true;
  sigmaTime=sigmaTimes_Phi;
  
  SolveAmbiguity(sigmaSlope,sigmaTime,hit_seg,nrVarv,change_code);
  
  if(t0_corr) Set_HITTime_T0(hit_seg,Get_T0_fin());
  
   if(!Get_FitOK()){
     Set_IsGood(false);  //qui
     return hit_seg;
   }
  
  Set_IsDouble(false);
  Set_IsGood(true);	
  
  if(DEBUG_SEG) 
    printf("Linear Fit: is good = %d, slope = %.2f, X0 = %.1f, t0 = %.1f, NPT = %d, chi2 = %.3f\n",Get_FitOK(),Get_Slope(),Get_X0(),Get_T0(),Get_NPT(),Get_Chi2()); 
  
  return hit_seg;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


HITColl_Seg * Segment::Get_SLSegment(HITColl_Seg *hit_seg, float w_xcorr, TimeCorr *corr, float T0phi)
{
  if(DEBUG_SEG) printf("\n ...Start computing THETA SEGMENT...\n");
  
  HIT *hit=hit_seg->hit(0);
  int SL = hit->SL_ID();
  
  seg_Slope=-999.; seg_X0=-999.; seg_T0=-999.; seg_T0_fin=-999.; seg_NPT=-999; seg_chi2=-999.; 
  seg_isdouble=false; 
  seg_isgood=true; 
  seg_okFit=false; 
  onlyoneseg=true;
  
  hit_seg = RejectDoubleHits(hit_seg);
  if(DEBUG_SEG) printf("n.hit after RejectDoubleHits: %d\n",hit_seg->Get_NumberHITS());
  if(hit_seg->Get_NumberHITS()<MinNHit_Theta){
    if(DEBUG_SEG)
      printf("After RejectDoubleHits seg. has less then %d hits,\n I cannot say anything about SLOPE\n",MinNHit_Theta);
    Set_IsGood(false);
    return hit_seg;
  }  
  
  // drift time correction
  if(drift_corr) Set_HITTime_Drift(hit_seg,w_xcorr,corr);
  // T0 time correction
  if(t0_corr) Set_HITTime_T0(hit_seg,T0phi);
  
  pair<float,float> seg_par0=make_pair(0,0);
  vector< pair<float,float> > _mean2;
  int size2=0;
  _mean2=Xmean_SL(hit_seg, SL,size2);
  
  if(DEBUG_SEG){
    printf("Size SL 2 %d\n",size2);
    printf(" N.hit tot=%d - if double: N.hit1=%d, N.hit2=%d\n",ii,ii1,ii2);
  }
  
  if(size2==0){
    Set_IsGood(false);
    return hit_seg;
  }
  
  if(size2==2){
    if(DEBUG_SEG)
      printf("Seg is double, segment rejected\n");
    Set_IsDouble(true);
    Set_IsGood(false);
    return hit_seg;
  }
  
  pair<float,float> seg_par=make_pair(0,0);
  
  double sigmaSlope=0;
  double sigmaTime=0;
  sigmaSlope=corr->Get_SigmaAngCorr(seg_par.first);
  if(DEBUG_SEG) printf("Sigma SlopeCorr: %.2f\n",sigmaSlope);
  
  if(DEBUG_SEG)
    printf("\nHo una sola traccia---->OK! continuo con l'analisi\n");
  
  hit_seg = SelectHITsforSLFIT(hit_seg,_mean2[0].first);
  if(hit_seg->Get_NumberHITS()<MinNHit_Theta){
    if(DEBUG_SEG) 
      printf("Dopo il SelectHITforSLFIT la traccia ha meno di %d hit e non posso dire nulla sulla slope\n --> rigetto il segmento\n",MinNHit_Theta);
    Set_IsGood(false);
    return hit_seg;
  }
  
  int nrVarv=2;
  bool change_code=false;
  sigmaTime=sigmaTimes_1SL;
  
  SolveAmbiguity(sigmaSlope,sigmaTime,hit_seg,nrVarv,change_code);
  
  if(DEBUG_SEG) 
    printf("Linear Fit: is good = %d, slope = %.2f, X0 = %.1f, t0 = %.1f, NPT = %d, chi2 = %.3f\n",Get_FitOK(), Get_Slope(),Get_X0(),Get_T0(),Get_NPT(),Get_Chi2()); 
  
  if(!Get_FitOK()){
    Set_IsGood(false);
    return hit_seg;
  }
  
  seg_okFit=false; // lo riazzero!
  seg_par=make_pair(Get_Slope(),Get_X0());
  // slope time correction
  if(slope_corr) Set_HITTime_Slope(hit_seg,seg_par,corr);
  // linearity time correction
  if(lin_corr) Set_HITTime_Linearity(hit_seg,seg_par,corr);
  sigmaSlope=corr->Get_SigmaAngCorr(seg_par.first);
  if(DEBUG_SEG) printf("Sigma SlopeCorr: %.2f\n",sigmaSlope);
  
  hit_seg = SelectHITsforFIT(hit_seg,seg_par);
  
  if(hit_seg->Get_NumberHITS()<MinNHit_Theta){
    if(DEBUG_SEG) 
      printf("Dopo il SelectHITforSLFIT la traccia ha meno di 2 hit e non posso dire nulla sulla slope\n --> butto il segmento\n");
    Set_IsDouble(false);
    Set_IsGood(false);
    return hit_seg;
  }
  
  nrVarv=2;
  change_code=true;
  sigmaTime=sigmaTimes_1SL;
  
  SolveAmbiguity(sigmaSlope,sigmaTime,hit_seg,2,change_code);
  
  if(DEBUG_SEG) 
    printf("Linear Fit: slope = %.2f, X0 = %.1f, NPT = %d, chi2 = %.3f\n",Get_Slope(),Get_X0(),Get_NPT(),Get_Chi2()); 
  if(!Get_FitOK()){
    Set_IsGood(false);
    return hit_seg;
  }
  Set_IsDouble(false);
  Set_IsGood(true);
  
  return hit_seg;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// Check for more than 1 hit/layer;
// if (N.hit/layer) are more than 3, reject hits of that layer
// TODO: potrei fare questo controllo con le HITColl_Layer,
//       pero` dovrei legare il numero di hit dell HITCOll_Seg
//       con il numero della hit del HITColl_Layer
HITColl_Seg * Segment::RejectDoubleHits(HITColl_Seg *hit_seg)
{

  bool hit_deleted=0;
  if(DEBUG_SEG){
    printf(" Starting Segment::RejectDoubleHits \n");
    printf("hit_seg->Get_NumberHITS() = %d\n",hit_seg->Get_NumberHITS()); 
    
    for(int i=0;i<hit_seg->Get_NumberHITS();i++){
      
      HIT *hit=hit_seg->hit(i);
      cout<<"N.hit "<< i <<", Event "<<hit->EvNumber()<<", SL "<<hit->SL_ID()<<", L "<<hit->L_ID()<<", wire "<<hit->wire_ID()<<": l "<<hit->code_ID().first<<", r "<<hit->code_ID().second<< ", xl "<<hit->Get_XPair(hit->dtime()).first<<  ", xr "<<hit->Get_XPair(hit->dtime()).second << ", y "<<hit->y_wire_ID()<<endl;      
      
    }
  }
  
  int hit_xL[8]={0,0,0,0,0,0,0,0};
  for(int i=0;i<hit_seg->Get_NumberHITS();i++){
    HIT *hit=hit_seg->hit(i);
    if(hit->SL_ID()==1 || hit->SL_ID()==2 ){
      hit_xL[(hit->L_ID())-1]++;
    }
    else{
      hit_xL[(hit->L_ID())+3]++;
    }
  }
  
  
  // Check for more than 1 hit/layer;
  for(int j=0;j<8;j++){
    if(hit_xL[j]>1){
      
      int n_hit[ hit_xL[j] ];
      int w_hit[ hit_xL[j] ];
      
      for(int i=0;i<hit_xL[j];i++)
	{n_hit[i]=0;w_hit[i]=0;}
      int wi=0;
      if(DEBUG_SEG) 
	printf("L %d has n. hit %d\n",j+1,hit_xL[j]);
      
      for(int i=0;i<hit_seg->Get_NumberHITS();i++){
	HIT *hit=hit_seg->hit(i);
	if( ((hit->SL_ID()==1 || hit->SL_ID()==2 ) && hit->L_ID()==(j+1))
	    || ( hit->SL_ID()==3 && hit->L_ID()==(j-3))){
	  n_hit[wi]=i;
	  w_hit[wi]=hit->wire_ID();
	  wi++;
	}
      }
      
      if(DEBUG_SEG && wi>1) 
	printf("n.hit prima: %d\n",hit_seg->Get_NumberHITS());
      // if (N.hit/layer) are more than 3, reject hits of that layer
      if(wi>3) 
	for(int k=0;k<wi;k++) {
	  HIT *hit=hit_seg->hit(n_hit[wi-1-k]);
	  hit->Change_LCode(false);
	  hit->Change_RCode(false);
	  hit_seg->eraseHIT(n_hit[wi-1-k]);
	  hit_deleted=1;
	}
      else 
	if(wi==2){
	  if(w_hit[0]==w_hit[1]){
	    // 	  HIT *hit=hit_seg->hit(n_hit[1]);
	    // 	  hit->Change_LCode(false);
	    // 	  hit->Change_RCode(false);
	    //    hit_seg->eraseHIT(n_hit[1]);
	    //    hit_deleted=1;
	  }
	}
	else 
	  if(wi==3){
	    if(w_hit[0]==w_hit[1] && w_hit[0]==w_hit[2]){
	      HIT *hit=hit_seg->hit(n_hit[2]);
	      hit->Change_LCode(false);
	      hit->Change_RCode(false);
	      hit_seg->eraseHIT(n_hit[2]);
	      //  hit_seg->eraseHIT(n_hit[1]);
	      // hit_deleted=1;
	    }	    
	    else
	      if(w_hit[1]==w_hit[2]){
		//	hit_seg->eraseHIT(n_hit[2]);
		//	hit_deleted=1;
	      }
	      else if(w_hit[0]==w_hit[2]){
		//	hit_seg->eraseHIT(n_hit[2]);
		//	hit_deleted=1;
	      }
	      else if(w_hit[0]==w_hit[1]){
		//	hit_seg->eraseHIT(n_hit[1]);
		//	hit_deleted=1;
	      }	  
	      else if(w_hit[0]!=w_hit[1] && w_hit[0]!=w_hit[2] && w_hit[1]!=w_hit[2] ){
		for(int k=0;k<wi;k++) {
		  HIT *hit=hit_seg->hit(n_hit[wi-1-k]);
		  hit->Change_LCode(false);
		  hit->Change_RCode(false);
		  hit_seg->eraseHIT(n_hit[wi-1-k]);   // da sistemare!!!!
		  hit_deleted=1;
		}
	      }
	  }
      
      if(DEBUG_SEG && wi>1) printf("n.hit dopo: %d\n",hit_seg->Get_NumberHITS());
    } // close if(hit_xL[j]>1)
  } // close check layer with more than 1 hit
  
  if(DEBUG_SEG && hit_deleted==1)
    for(int i=0;i<hit_seg->Get_NumberHITS();i++){
      
      HIT *hit=hit_seg->hit(i);
      cout<<"N.hit "<< i <<", Event "<<hit->EvNumber()<<", SL "<<hit->SL_ID()<<", L "<<hit->L_ID()<<", wire "<<hit->wire_ID()<<": l "<<hit->code_ID().first<<", r "<<hit->code_ID().second<< ", xl "<<hit->Get_XPair(hit->dtime()).first<<  ", xr "<<hit->Get_XPair(hit->dtime()).second << ", y "<<hit->y_wire_ID()<<endl;      
      
    }
  
  return hit_seg;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// cerco le semicelle con piu' hit e determino X_medio e Y_medio
vector< pair<float,float> > Segment::Xmean_SL(HITColl_Seg *hit_seg, int SL, int & size ) 
{
  float mean0=0, mean1=0, mean2=0; 
  float mean[2]={0.,0.};
  ii=0, ii1=0, ii2=0;
  int w_min=100;
  int w_max=0;
  
  HIT *hit=hit_seg->hit(0);
  int CH=hit->CH_ID();
  
  for(int i=0;i<hit_seg->Get_NumberHITS();i++){
    HIT *hit=hit_seg->hit(i);
    if(hit->SL_ID()==SL){
      mean0 += hit->wire_ID();
      ii++;
      if(hit->wire_ID()<w_min) w_min=hit->wire_ID();
      if(hit->wire_ID()>w_max) w_max=hit->wire_ID();
    }
  }
  
  if(ii!=0)
    if( fabs(w_max-w_min)<=2 ){
      mean0 /= ii;     
      mean[0] = mean0;
      mean[1] = -999.;
      
      if(DEBUG_SEG){
	printf("\nfabs(w_max-w_min)<=2\n");
	printf("SL %d Filo medio = %.1f, Filo min = %d, Filo max = %d,\n",SL,mean[0],w_min,w_max);
      }
    }
    else{
      if(fabs(w_max-w_min)<=4 ){
	for(int i=0;i<hit_seg->Get_NumberHITS();i++){
	  HIT *hit=hit_seg->hit(i);
	  if(hit->SL_ID()==SL){
	    if(fabs(hit->wire_ID()-(w_min))<2){
	      mean1 += hit->wire_ID();
	      ii1++;
	    }
	    if(fabs(hit->wire_ID()-(w_max))<2){
	      mean2 += hit->wire_ID();
	      ii2++;
	    }
	  }
	}
	if(DEBUG_SEG){
	  printf("\nfabs(w_max-w_min)<=4\n");
	  printf("SL %d Filo medio = %.1f, Filo min = %d, Filo max = %d,\n",SL,mean[0],w_min,w_max);
	}
      }
      else{
	for(int i=0;i<hit_seg->Get_NumberHITS();i++){
	  HIT *hit=hit_seg->hit(i);
	  if(hit->SL_ID()==SL){
	    if(fabs(hit->wire_ID()-(w_min+1))<=2){
	      mean1 += hit->wire_ID();
	      ii1++;
	    }
	    if(fabs(hit->wire_ID()-(w_max-1))<=2){
	      mean2 += hit->wire_ID();
	      ii2++;
	    }
	  }
	}
	if(DEBUG_SEG){
	  printf("\nfabs(w_max-w_min)>4\n");
	  printf("SL %d Filo medio = %.1f, Filo min = %d, Filo max = %d,\n",SL,mean[0],w_min,w_max);
	}
      }
      
      if(ii1>=2 && ii2>=2){
	if(fabs(ii1-ii2)<=1){
	  mean[0]=mean1/ii1;
	  mean[1]=mean2/ii2;
	  if(DEBUG_SEG) 
	    printf("SL %d Filo min = %d, Filo max = %d,\n",SL,w_min,w_max);
	  if(DEBUG_SEG) 
	    printf("SL %d Filo mean 1 = %.1f, Filo mean 2 = %.1f\n",SL,mean[0],mean[1]);
	}
	else if(ii1>ii2){
	  mean[0]=mean1/ii1;
	  mean[1]=-999.;
	  if(DEBUG_SEG) 
	    printf("SL %d Filo min = %d, Filo max = %d,\n",SL,w_min,w_max);
	  if(DEBUG_SEG) 
	    printf("SL %d Filo mean 1 = %.1f, Filo mean 2 = %.1f\n",SL,mean[0],mean[1]);
	}
	else if(ii2>ii1){
	  mean[0]=mean2/ii2;
	  mean[1]=-999.;
	  if(DEBUG_SEG) 
	    printf("SL %d Filo min = %d, Filo max = %d,\n",SL,w_min,w_max);
	  if(DEBUG_SEG) 
	    printf("SL %d Filo mean 1 = %.1f, Filo mean 2 = %.1f\n",SL,mean[0],mean[1]);
	}
      }
      else if(ii1>=2 || ii2>=2){
	if(ii1>ii2){
	  mean[0]=mean1/ii1;
	  mean[1]=-999.;
	  if(DEBUG_SEG) 
	    printf("SL %d Filo min = %d, Filo max = %d,\n",SL,w_min,w_max);
	  if(DEBUG_SEG) 
	    printf("SL %d Filo mean 1 = %.1f, Filo mean 2 = %.1f\n",SL,mean[0],mean[1]);
	}
	else if(ii2>ii1){
	  mean[0]=mean2/ii2;
	  mean[1]=-999.;
	  if(DEBUG_SEG) 
	    printf("SL %d Filo min = %d, Filo max = %d,\n",SL,w_min,w_max);
	  if(DEBUG_SEG) 
	    printf("SL %d Filo mean 1 = %.1f, Filo mean 2 = %.1f\n",SL,mean[0],mean[1]);
	}
      }
    }
  
  vector< pair<float,float> > _mean;
  
  if(mean[0]!=0 && mean[1]!=0)
    if(mean[1]==-999.)
      _mean.push_back(Get_Xmean_SL(hit_seg,CH,SL,mean[0]));
    else{
      _mean.push_back(Get_Xmean_SL(hit_seg,CH,SL,mean[0]));
      _mean.push_back(Get_Xmean_SL(hit_seg,CH,SL,mean[1]));
    }
  
  size=_mean.size();
  return _mean;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


pair<float,float> Segment::Get_Xmean_SL(HITColl_Seg *hit_seg, int CH, int SL, float mean){
  
  Geom *geo=new Geom();
  
  float x_mean=0;
  x_mean=geo->get_x_wire(CH, SL, 1, int(mean));
  
  if(DEBUG_SEG) 
    printf("SL %d X_filo medio = %.1f\n",SL,x_mean);
  
  int n_sc=7; // numero di semicelle da controllare
  int check_nhit[n_sc]; for(int i=0;i<n_sc;i++) check_nhit[i]=0;
  float check_xhit[n_sc]; for(int i=0;i<n_sc;i++) check_xhit[i]=0;
  float check_yhit[n_sc]; for(int i=0;i<n_sc;i++) check_yhit[i]=0;
  
  for(int i=0;i<hit_seg->Get_NumberHITS();i++){
    
    HIT *hit=hit_seg->hit(i);
    if(hit->SL_ID()==SL){
      for(int j=0; j<n_sc; j++){

	if( (hit->Get_XPair(hit->dtime()).first) >= (x_mean+(-4+j)*2.1)
	    &&
	    (hit->Get_XPair(hit->dtime()).first) <= (x_mean+(-3+j)*2.1) ){
	  check_xhit[j] += (hit->Get_XPair(hit->dtime()).first);
	  check_yhit[j] += (hit->y_wire_ID());
	  check_nhit[j] ++;
	}
	if( (hit->Get_XPair(hit->dtime()).second) >= (x_mean+(-4+j)*2.1)
	    &&
	    (hit->Get_XPair(hit->dtime()).second) <= (x_mean+(-3+j)*2.1) ){
	  check_xhit[j] += (hit->Get_XPair(hit->dtime()).second);
	  check_yhit[j] += (hit->y_wire_ID());
	  check_nhit[j] ++;
	}	
      }
    }
    
  }
  for(int j=0; j<n_sc; j++){
    if(check_nhit[j]!=0){
      check_xhit[j] /= check_nhit[j];
      check_yhit[j] /= check_nhit[j];
    }
  }  
  
  float Xm=0; 
  float Ym=0; 
  int Nn=0;
  int count=0; 
  for(int j=0; j<n_sc; j++){
    if(check_nhit[j]>Nn){
      Xm=check_xhit[j];
      Ym=check_yhit[j];
      Nn=check_nhit[j];
      count=1;
    }
    else if(check_nhit[j]==Nn && check_nhit[j]!=0){
      Xm += check_xhit[j];
      Ym += check_yhit[j];
      count ++;
      Nn=check_nhit[j];
    }
  }
  if(count>1) {
    Xm /= count;
    Ym /= count;
  }
  
  if(DEBUG_SEG) {
    printf("N.hit into collection %d\n", hit_seg->Get_NumberHITS());
    for(int j=0; j<n_sc; j++) 
      printf("SL%d: N.hit x semicella %d, X_hits %.1f, Y_hits %.1f\n", SL, check_nhit[j], check_xhit[j], check_yhit[j]);
    printf("        Xm con piu' hit = %.1f\n", Xm);
    printf("        Ym con piu' hit = %.1f\n", Ym);
  }
  
  delete geo;

  return make_pair(Xm,Ym);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


pair<float,float> Segment::Get_SlopePHI(pair<float,float> PHI1, pair<float,float> PHI2)
{
  float Slope=-999.;
  float X0=-999.;
  Slope=(PHI1.first-PHI2.first)/(PHI1.second-PHI2.second);
  X0=-PHI1.second*Slope+PHI1.first;
  
  if(DEBUG_SEG){
    printf("...............Slope %.2f, X0 %.1f\n",Slope,X0);
    printf("...............Angle %.1f, X0 %.1f\n",atan(Slope)*180./3.14,X0);
  }
  
  return make_pair(Slope,X0);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void Segment::Set_HITTime_Drift(HITColl_Seg *hit_seg, float w_xcorr, TimeCorr *corr)
{
  for(int i=0;i<hit_seg->Get_NumberHITS();i++){
    HIT *hit=hit_seg->hit(i);
    
    float time=hit->dtime();
    if(DEBUG_SEG)
      printf("(time-t0-ttrig) = %.1f(ns) ",hit->dtime());
    time=corr->Get_Time_DriftCorr(time, hit->SL_ID(), w_xcorr);
    if(DEBUG_SEG)
      printf("----> time after DriftCorr = %.1f(ns)\n",time);
    
    hit->Set_CorrTime(time);
  }
  
  return;
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void Segment::Set_HITTime_Slope(HITColl_Seg *hit_seg, pair<float,float> seg_par, TimeCorr *corr)
{
  
  for(int i=0;i<hit_seg->Get_NumberHITS();i++){
    HIT *hit=hit_seg->hit(i);
    
    float time=hit->dtime();
    if(DEBUG_SEG)
      printf("time %.1f(ns) ",hit->dtime());
    time=corr->Get_Time_SlopeCorr(time, seg_par.first);
    if(DEBUG_SEG)
      printf("----> after SlopeCorr = %.1f(ns)\n",time);
    hit->Set_CorrTime(time);
  }
  
  return;
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void Segment::Set_HITTime_Linearity(HITColl_Seg *hit_seg, pair<float,float> seg_par, TimeCorr *corr)
{
  
  for(int i=0;i<hit_seg->Get_NumberHITS();i++){
    HIT *hit=hit_seg->hit(i);
    
    float time=hit->dtime();
    if(DEBUG_SEG)
      printf("time %.1f(ns) ",hit->dtime());
    time=corr->Get_LinearityCorr(hit->CH_ID(),hit->SL_ID(),time, seg_par.first);
    if(DEBUG_SEG)
      printf("----> after LinerityCorr = %.1f(ns)\n",time);
    hit->Set_CorrTime(time);
  }
  
  return;
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void Segment::Set_HITTime_T0(HITColl_Seg *hit_seg, float T0)
{
  //  if(T0==-999.) 
  if(TMath::Abs(T0)>T0_max) 
    return;
  
  for(int i=0;i<hit_seg->Get_NumberHITS();i++){
    HIT *hit=hit_seg->hit(i);
    
    if(DEBUG_SEG)
      printf("time %.1f(ns) ",hit->dtime());
    
    hit->Set_CorrTime(hit->dtime()-T0);
    
    if(DEBUG_SEG)
      printf("----> after t0 subtraction = %.1f(ns)\n",hit->dtime());
  }
  
  return;
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// change hit->Code_ID of each hit of collection,
// if the distance between one of the hit_pair and the raw_linear_fit
// is more than 1 cm, the Code_ID of that hit_pair is set false
// If left and right solutions are false, the hit is rejected 
HITColl_Seg * Segment::SelectHITsforFIT(HITColl_Seg *hit_seg, pair<float,float> seg_par)
{
  if(DEBUG_SEG) 
    printf("\nStarting SelectHITsforFIT \nN.hit into SelectHITsforFIT: %d\n",hit_seg->Get_NumberHITS());
  
  float cutx = cut_SLs + fabs(seg_par.first);

  int Nhit=hit_seg->Get_NumberHITS();
  if(DEBUG_SEG) printf("N.hit into SelectHITsforFIT: %d\n",Nhit);
  
  int Ntoberejected=0;
  
  if(DEBUG_SEG)
    printf("  Check if the distance between hits and temptative track is less than %.1f\n",cutx);
  
  for(int i=0;i<hit_seg->Get_NumberHITS();i++){
    if(DEBUG_SEG) 
      printf("N.hit %d\n",i);
    
    HIT *hit=hit_seg->hit(i);
    
    if(DEBUG_SEG){
      printf("  Xl=%.1f, Xr=%.1f, Y=%.1f, \n",hit->Get_XPair(hit->dtime()).first, hit->Get_XPair(hit->dtime()).second, hit->y_wire_ID());
      printf("  distanza dalla retta l %.2f  r %.2f\n",(hit->Get_XPair(hit->dtime()).first) -( seg_par.first*hit->y_wire_ID() + seg_par.second),(hit->Get_XPair(hit->dtime()).second)-( seg_par.first*hit->y_wire_ID() + seg_par.second) );  
    }
    
    if( fabs((hit->Get_XPair(hit->dtime()).first) -( seg_par.first*hit->y_wire_ID() + seg_par.second))>cutx ) 
      hit->Change_LCode(false);
    if( fabs((hit->Get_XPair(hit->dtime()).second)-( seg_par.first*hit->y_wire_ID() + seg_par.second))>cutx ) 
      hit->Change_RCode(false);
    
    if(hit->code_ID().first==false && hit->code_ID().second==false)
      Ntoberejected++;
    
  } // close for()
  
  if(DEBUG_SEG)
    printf("N.hit to be rejected %d\n",Ntoberejected);
  
  for(int i=Nhit-1;i>=0;i--){
    HIT *hit=hit_seg->hit(i);
    if(hit->code_ID().first==false && hit->code_ID().second==false){
      hit_seg->eraseHIT(i); 
      if(DEBUG_SEG)
	printf("....deleting hit %d...\n",i);
    }
  }
  
  if(DEBUG_SEG)
  {
    printf("N.hit into collection at the end of SelectHITsforFIT: %d\n", hit_seg->Get_NumberHITS());
    for(int i=0;i<hit_seg->Get_NumberHITS();i++){
      HIT *hit=hit_seg->hit(i);
      cout<<"N.hit "<< i <<": left = "<<hit->code_ID().first<<", right = "<<hit->code_ID().second<<endl; 
    }   
  }

  return hit_seg;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// change hit->Code_ID of each hit of collection,
// if the distance between one of the hit_pair and the mean 
// is more than 6 cm???, the Code_ID of that hit_pair is set false
// If left and right solutions are false, the hit is rejected 
HITColl_Seg * Segment::SelectHITsforSLFIT(HITColl_Seg *hit_seg, float xmean)
{
  if(DEBUG_SEG) 
    printf("\nStarting SelectHITsforSLFIT \nN.hit into SelectHITsforSLFIT: %d\n",hit_seg->Get_NumberHITS());
  
  int Nhit=hit_seg->Get_NumberHITS();
  if(DEBUG_SEG) printf("N.hit into SelectHITsforFIT: %d\n",Nhit);
  
  int Ntoberejected=0;
  
  for(int i=0;i<hit_seg->Get_NumberHITS();i++){
    if(DEBUG_SEG) 
      printf("N.hit %d\n",i);
    
    HIT *hit=hit_seg->hit(i);
    
    if(DEBUG_SEG){
      printf("  Xl=%.1f, Xr=%.1f, Y=%.1f, \n",hit->Get_XPair(hit->dtime()).first, hit->Get_XPair(hit->dtime()).second, hit->y_wire_ID());
//       printf("  distanza dalla retta l %.2f  r %.2f\n",(hit->Get_XPair(hit->dtime()).first) -( hit->y_wire_ID() - seg_par.second)/seg_par.first,(hit->Get_XPair(hit->dtime()).second)-( hit->y_wire_ID() - seg_par.second)/seg_par.first);  
      printf("  distanza dalla media l %.2f  r %.2f\n",(hit->Get_XPair(hit->dtime()).first) - xmean,(hit->Get_XPair(hit->dtime()).second) - xmean );  
    }
    
    if( fabs((hit->Get_XPair(hit->dtime()).first) - xmean ) > cut_SL ) 
      hit->Change_LCode(false);
    if( fabs((hit->Get_XPair(hit->dtime()).second) - xmean ) > cut_SL ) 
      hit->Change_RCode(false);
    
    if(hit->code_ID().first==false && hit->code_ID().second==false)
      Ntoberejected++;
    
  } // close for()
  
  if(DEBUG_SEG)
    printf("N.hit to be rejected %d\n",Ntoberejected);
  
  for(int i=Nhit-1;i>=0;i--){
    HIT *hit=hit_seg->hit(i);
    if(hit->code_ID().first==false && hit->code_ID().second==false){
      hit_seg->eraseHIT(i); 
      if(DEBUG_SEG)
	printf("....deleting hit %d...\n",i);
    }
  }
  
  if(DEBUG_SEG)
  {
    printf("N.hit into collection at the end of SelectHITsforFIT: %d\n", hit_seg->Get_NumberHITS());
    for(int i=0;i<hit_seg->Get_NumberHITS();i++){
      HIT *hit=hit_seg->hit(i);
      cout<<"N.hit "<< i <<": left = "<<hit->code_ID().first<<", right = "<<hit->code_ID().second<<endl; 
    }   
  }
  
  return hit_seg;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


int Segment::Get_NnotSolvedHIT(HITColl_Seg *hit_seg){
  int N=0;
  for(int i=0;i<hit_seg->Get_NumberHITS();i++){
    HIT *hit=hit_seg->hit(i);
    if( (hit->code_ID().first==true) && (hit->code_ID().second==true) )
      N++;
  }  
  return N;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


int Segment::Get_NSolvedHIT(HITColl_Seg *hit_seg){
  int Ns=0;
  for(int i=0;i<hit_seg->Get_NumberHITS();i++){
    HIT *hit=hit_seg->hit(i);
    if( ((hit->code_ID().first==true) && (hit->code_ID().second==false))
	||
	((hit->code_ID().first==false) && (hit->code_ID().second==true)) )
      Ns++;
  }  
  return Ns;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


bool Segment::HITisDouble(HITColl_Seg *hit_seg,int k){
  bool isDouble=false;
  HIT *hit=hit_seg->hit(k);
  int SL=hit->SL_ID();
  int L=hit->L_ID();
  for(int i=0;i<hit_seg->Get_NumberHITS();i++){
    HIT *hit=hit_seg->hit(i);
    if(i!=k)
      if(hit->SL_ID()==SL && hit->L_ID()==L)
        isDouble=true;
  }
  
  return isDouble;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void Segment::Set_DoubleHITs(HITColl_Seg *hit_seg){
  
  for(int i=0;i<hit_seg->Get_NumberHITS();i++){
    HIT *hit=hit_seg->hit(i);
    bool itis=HITisDouble(hit_seg,i);
    hit->Set_IsDouble(itis);
  }
  return;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void Segment::SolveAmbiguity(double sigmaSlope,double sigmaTime,HITColl_Seg *hit_seg,int nrVarv,bool change_code){
  
    //ALTEA 2
    //TH1F *chi2distribution = new TH1F("chi2distribution", "chi2distribution", 100000, 0, 100000);
    Int_t Nhits = 0;
    int nrPnts = 0; // N.layer con N.hit>0
    Int_t NnotSolv = 0;
    Int_t NSolv = 0;

    Int_t nrMinPnts=0, nrMinPnts0=0;
    HIT *hit=hit_seg->hit(0);
    if( (hit->CH_ID()==11 || hit->CH_ID()==10) && (hit->SL_ID()!=2) )
        if(!onlyoneseg)
            nrMinPnts=MinNHit_Phi;
        else
            nrMinPnts=MinNHit_1SLPhi;
    else
        nrMinPnts=MinNHit_Theta;
  
    nrMinPnts0=nrMinPnts;
  
    Nhits = hit_seg->Get_NumberHITS();
    NnotSolv=Get_NnotSolvedHIT(hit_seg);
    NSolv=Get_NSolvedHIT(hit_seg);
    Set_DoubleHITs(hit_seg);
    if(DEBUG_SEG) 
        printf("N.Hit to fit %d, solved %d, not solved %d\n",Nhits,NSolv,NnotSolv);
    for(int i=0;i<hit_seg->Get_NumberHITS();i++){ 
        HIT *hit=hit_seg->hit(i);
        if(DEBUG_SEG) 
            if(hit->IsDouble()) printf("Hit.%d is double\n",i);  
    }
  
    int nL=0;
    //  HIT *hit=hit_seg->hit(0);
    if(hit->SL_ID()==2) nL=4;
    else nL=8;
  
  L_2hit=NULL;
  L_2hit=new int [nL];
  for(int i=0;i<nL;i++) L_2hit[i]=0;  
  Get_DoubleHITxL(hit_seg,nL); 
    
  nrPnts=Nhits;
  for(int i=0;i<nL;i++)
    if(L_2hit[i]==2) nrPnts -= 1;
    else if(L_2hit[i]==3) nrPnts -= 2;
  
  if(DEBUG_SEG){
    printf("N hits %d\n",Nhits);
    for(int i=0;i<nL;i++)printf(" L_2hit[%d]=%d\n",i,L_2hit[i]);
    printf("nrPnts %d\n",nrPnts);
  }
  
  delete [] L_2hit;
  L_2hit=NULL;
  

  if(nrPnts<=2){
    if(DEBUG_SEG)
      printf("Ho meno di 2 hit su due diversi Layer--->NON POSSO FARE IL FIT ed ESCO!!!\n");
//     delete [] NX_xLs;
//     NX_xLs=NULL;
    
    return;
  }
  
  // *************************************************************************************
    
  NX_xLs=NULL;
  NX_xLs=new int [nL];

  for(int i=0;i<nL;i++) NX_xLs[i]=0;  
  Get_NX_xLs(hit_seg,nL); 
  if(DEBUG_SEG) for(int i=0;i<nL;i++)cout<<" NX_xLs["<<i<<"] "<<NX_xLs[i]<<endl;
  
  if(DEBUG_SEG)
    for(int i=0;i<hit_seg->Get_NumberHITS();i++){
      HIT *hit=hit_seg->hit(i);
      printf("N.hit %d\n",i);
      hit->print();
    }

  // *************************************************************************************


  // TODO da sistemare;
  //   HITColl_Layer *layer[8]; 
  for(int i=0;i<8;i++){
    layer[i]= new HITColl_Layer;
    layer[i]->selectHIT(hit_seg, i+1);
  }
  
  int N[8];
  for(int i=0;i<8;i++) N[i]=layer[i]->Get_NumberHITS();
  
//   if(DEBUG_SEG)
//     for(int i=0;i<8;i++) printf("Numero di HIT x Layer %d = %d\n",i+1,N[i]);
  
  bool nohit[8];    for(int i=0;i<8;i++) nohit[i]=false;
  double xx[8];      for(int i=0;i<8;i++) xx[i]=0.;
  double yy[8];      for(int i=0;i<8;i++) yy[i]=0.;
  int cc[8];        for(int i=0;i<8;i++) cc[i]=0;
  int cc_fin[8];    for(int i=0;i<8;i++) cc_fin[i]=0;
  int hh_fin[8];    for(int i=0;i<8;i++) hh_fin[i]=0;
  int hh[8];        for(int i=0;i<8;i++) hh[i]=0;
  int c2h[8];       for(int i=0;i<8;i++) c2h[i]=1;
  for(int i=0;i<8;i++) 
    if(N[i]==0) { N[i]=1; nohit[i]=true; }
  
  int n_comb=1; for(int i=0;i<nL;i++) if(NX_xLs[i]>1) n_comb *=NX_xLs[i];
  if(DEBUG_SEG) printf("N. combinations %d\n\n",n_comb);  
  bool okFit = false;
  bool okFit0 = false;
  int comb=0;
  int k;   
  double Chi2[8]; for(int i=0;i<8;i++) Chi2[i]=100000.;  //CHI2 threshold
  double chi2=-999.;
  int NPT=0;  
  double chi20=0;
  int NPT0=0;

  pair<double,double>  m=make_pair(0,0); 
  pair<double,double>  q=make_pair(0,0);
  pair<double,double>  t0=make_pair(0,0);

  pair<double,double>  m0=make_pair(0,0); 
  pair<double,double>  q0=make_pair(0,0);
  pair<double,double>  t00=make_pair(0,0);

  nL=8;
  
  Double_t ax[nrPnts]; for(int i=0;i<nrPnts;i++) ax[i]=0.;
  Double_t ay[nrPnts]; for(int i=0;i<nrPnts;i++) ay[i]=0.;
  Int_t ac[nrPnts]; for(int i=0;i<nrPnts;i++) ac[i]=0;
  Int_t ah[nrPnts]; for(int i=0;i<nrPnts;i++) ah[i]=0;
  Int_t ac_fin[nrPnts]; for(int i=0;i<nrPnts;i++) ac_fin[i]=0;
//   Int_t ah_fin[nrPnts]; for(int i=0;i<nrPnts;i++) ah_fin[i]=0;
  
  for(int i0=0;i0<N[0];i0++){
    if(nohit[0]==true){xx[0]=-999.; hh[0]=0; cc[0]=0;}
    else{
      HIT * hit=layer[0]->hit(i0);
      yy[0]=hit->y_wire_ID();
      if(hit->IsSolved()){
	if(hit->code_ID().first==1){
	  xx[0]=hit->Get_XPair(hit->dtime()).first;
	  cc[0]=-1;
	  hh[0]=(i0+1);
	}
	if(hit->code_ID().second==1){
	  xx[0]=hit->Get_XPair(hit->dtime()).second;
	  cc[0]=1;
	  hh[0]=(i0+1);
	}
      }
      else 
	if(!(hit->IsSolved())){
	  if(c2h[0]==1){
	    xx[0]=hit->Get_XPair(hit->dtime()).first;
	    cc[0]=-1;
	    hh[0]=(i0+1);
	  }
	  else if(c2h[0]==2){
	    xx[0]=hit->Get_XPair(hit->dtime()).second;
	    cc[0]=1;
	    hh[0]=(i0+1);
	  }
	  c2h[0]++;
	  if(c2h[0]==3) c2h[0]=1;
	}
    }
    for(int i1=0;i1<N[1];i1++){
      if(nohit[1]==true){xx[1]=-999.; hh[1]=0; cc[1]=0;}
      else{
	HIT * hit=layer[1]->hit(i1);
	yy[1]=hit->y_wire_ID();
	if(hit->IsSolved()){
	  if(hit->code_ID().first==1){
	    xx[1]=hit->Get_XPair(hit->dtime()).first;
	    cc[1]=-1;
	    hh[1]=(i1+1);
	  }
	  if(hit->code_ID().second==1){
	    xx[1]=hit->Get_XPair(hit->dtime()).second;
	    cc[1]=1;
	    hh[1]=(i1+1);
	  }
	}
	else 
	  if(!(hit->IsSolved())){
	    if(c2h[1]==1){
	      xx[1]=hit->Get_XPair(hit->dtime()).first;
	      cc[1]=-1;
	      hh[1]=(i1+1);
	    }
	    else if(c2h[1]==2){
	      xx[1]=hit->Get_XPair(hit->dtime()).second;
	      cc[1]=1;
	      hh[1]=(i1+1);
	    }
	    c2h[1]++;
	    if(c2h[1]==3) c2h[1]=1;
	  }
      }
      for(int i2=0;i2<N[2];i2++){
	if(nohit[2]==true){xx[2]=-999; hh[2]=0; cc[2]=0;}
	else{
	  HIT * hit=layer[2]->hit(i2);
	  yy[2]=hit->y_wire_ID();
	  if(hit->IsSolved()){
	    if(hit->code_ID().first==1){
	      xx[2]=hit->Get_XPair(hit->dtime()).first;
	      cc[2]=-1;
	      hh[2]=(i2+1);
	    }
	    if(hit->code_ID().second==1){
	      xx[2]=hit->Get_XPair(hit->dtime()).second;
	      cc[2]=1;
	      hh[2]=(i2+1);
	    }
	  }
	  else 
	    if(!(hit->IsSolved())){
	      if(c2h[2]==1){
		xx[2]=hit->Get_XPair(hit->dtime()).first;
		cc[2]=-1;
		hh[2]=(i2+1);
	      }
	      else if(c2h[2]==2){
		xx[2]=hit->Get_XPair(hit->dtime()).second;
		cc[2]=1;
		hh[2]=(i2+1);
	      }
	      c2h[2]++;
	      if(c2h[2]==3) c2h[2]=1;
	    }
	}
	for(int i3=0;i3<N[3];i3++){
	  if(nohit[3]==true){xx[3]=-999; hh[3]=0; cc[3]=0;}
	  else{
	    HIT * hit=layer[3]->hit(i3);
	    yy[3]=hit->y_wire_ID();
	    if(hit->IsSolved()){
	      if(hit->code_ID().first==1){
		xx[3]=hit->Get_XPair(hit->dtime()).first;
		cc[3]=-1;
		hh[3]=(i3+1);
	      }
	      if(hit->code_ID().second==1){
		xx[3]=hit->Get_XPair(hit->dtime()).second;
		cc[3]=1;
		hh[3]=(i3+1);
	      }
	    }
	    else 
	      if(!(hit->IsSolved())){
		if(c2h[3]==1){
		  xx[3]=hit->Get_XPair(hit->dtime()).first;
		  cc[3]=-1;
		  hh[3]=(i3+1);
		}
		else if(c2h[3]==2){
		  xx[3]=hit->Get_XPair(hit->dtime()).second;
		  cc[3]=1;
		  hh[3]=(i3+1);
		}
		c2h[3]++;
		if(c2h[3]==3) c2h[3]=1;
	      }
	  }
	  for(int i4=0;i4<N[4];i4++){
	    if(nohit[4]==true){xx[4]=-999; hh[4]=0; cc[4]=0;}
	    else{
	      HIT * hit=layer[4]->hit(i4);
	      yy[4]=hit->y_wire_ID();
	      if(hit->IsSolved()){
		if(hit->code_ID().first==1){
		  xx[4]=hit->Get_XPair(hit->dtime()).first;
		  cc[4]=-1;
		  hh[4]=(i4+1);
		}
		if(hit->code_ID().second==1){
		  xx[4]=hit->Get_XPair(hit->dtime()).second;
		  cc[4]=1;
		  hh[4]=(i4+1);
		}
	      }
	      else 
		if(!(hit->IsSolved())){
		  if(c2h[4]==1){
		    xx[4]=hit->Get_XPair(hit->dtime()).first;
		    cc[4]=-1;
		    hh[4]=(i4+1);
		  }
		  else if(c2h[4]==2){
		    xx[4]=hit->Get_XPair(hit->dtime()).second;
		    cc[4]=1;
		    hh[4]=(i4+1);
		  }
		  c2h[4]++;
		  if(c2h[4]==3) c2h[4]=1;
		}
	    }
	    for(int i5=0;i5<N[5];i5++){
	      if(nohit[5]==true){xx[5]=-999; hh[5]=0; cc[5]=0;}
	      else{
		HIT * hit=layer[5]->hit(i5);
		yy[5]=hit->y_wire_ID();
		if(hit->IsSolved()){
		  if(hit->code_ID().first==1){
		    xx[5]=hit->Get_XPair(hit->dtime()).first;
		    cc[5]=-1;
		    hh[5]=(i5+1);
		  }
		  if(hit->code_ID().second==1){
		    xx[5]=hit->Get_XPair(hit->dtime()).second;
		    cc[5]=1;
		    hh[5]=(i5+1);
		  }
		}
		else 
		  if(!(hit->IsSolved())){
		    if(c2h[5]==1){
		      xx[5]=hit->Get_XPair(hit->dtime()).first;
		      cc[5]=-1;
		      hh[5]=(i5+1);
		    }
		    else if(c2h[5]==2){
		      xx[5]=hit->Get_XPair(hit->dtime()).second;
		      cc[5]=1;
		      hh[5]=(i5+1);
		    }
		    c2h[5]++;
		    if(c2h[5]==3) c2h[5]=1;
		  }
	      }
	      for(int i6=0;i6<N[6];i6++){
		if(nohit[6]==true){xx[6]=-999; hh[6]=0; cc[6]=0;}
		else{
		  HIT * hit=layer[6]->hit(i6);
		  yy[6]=hit->y_wire_ID();
		  if(hit->IsSolved()){
		    if(hit->code_ID().first==1){
		      xx[6]=hit->Get_XPair(hit->dtime()).first;
		      cc[6]=-1;
		      hh[6]=(i6+1);
		    }
		    if(hit->code_ID().second==1){
		      xx[6]=hit->Get_XPair(hit->dtime()).second;
		      cc[6]=1;
		      hh[6]=(i6+1);
		    }
		  }
		  else 
		    if(!(hit->IsSolved())){
		      if(c2h[6]==1){
			xx[6]=hit->Get_XPair(hit->dtime()).first;
			cc[6]=-1;
			hh[6]=(i6+1);
		      }
		      else if(c2h[6]==2){
			xx[6]=hit->Get_XPair(hit->dtime()).second;
			cc[6]=1;
			hh[6]=(i6+1);
		      }
		      c2h[6]++;
		      if(c2h[6]==3) c2h[6]=1;
		    }
		}
		for(int i7=0;i7<N[7];i7++){
		  if(nohit[7]==true){xx[7]=-999; hh[7]=0; cc[7]=0;}
		  else{
		    HIT * hit=layer[7]->hit(i7);
		    yy[7]=hit->y_wire_ID();
		    if(hit->IsSolved()){
		      if(hit->code_ID().first==1){
			xx[7]=hit->Get_XPair(hit->dtime()).first;
			cc[7]=-1;
			hh[7]=(i7+1);
		      }
		      if(hit->code_ID().second==1){
			xx[7]=hit->Get_XPair(hit->dtime()).second;
			cc[7]=1;
			hh[7]=(i7+1);
		      }
		    }
		    else 
		      if(!(hit->IsSolved())){
			if(c2h[7]==1){
			  xx[7]=hit->Get_XPair(hit->dtime()).first;
			  cc[7]=-1;
			  hh[7]=(i7+1);
			}
			else if(c2h[7]==2){
			  xx[7]=hit->Get_XPair(hit->dtime()).second;
			  cc[7]=1;
			  hh[7]=(i7+1);
			}
			c2h[7]++;
			if(c2h[7]==3) c2h[7]=1;
		      }
		  }
		  
		  k=0;
		  for(int n=0;n<8;n++){
		    if(xx[n]!=-999) {
		      ax[k]=xx[n];
		      ay[k]=yy[n];
		      ac[k]=cc[n];
		      ah[k]=hh[n];
		      k++;
		    }
		  }
		  if(DEBUG_SEG){
		    for(int i=0;i<nrPnts;i++) printf("ax[%d]=%.2f, ",i,ax[i]);
		    printf("\n");
		    for(int i=0;i<nrPnts;i++) printf("ay[%d]=%.2f, ",i,ay[i]);
		    printf("\n");
		    for(int i=0;i<nrPnts;i++) printf("ac[%d]=%d, ",i,ac[i]);
		    printf("\n");
		    for(int i=0;i<nrPnts;i++) printf("ah[%d]=%d, ",i,ah[i]);
		    printf("\n");
		  }
		  
		  chi20=0;
		  NPT0=0;
		  okFit0=false;
		  
		  fit->FIT_t0(sigmaSlope,sigmaTime,nrPnts,nrMinPnts0,nrVarv,ax,ay,ac,m0,q0,t00,chi20,NPT0,okFit0); //<----- FIT 
		  
		  if(DEBUG_SEG) printf("comb %d, okFit = %d, m = %.2f, q = %.0f, NPT = %d, chi2 = %.3f, T0 = %.1f\n\n", comb, okFit0, m0.first, q0.first, NPT0, chi20,t00.first);
		  if(okFit0)
		    if(NPT0>=NPT){ 
// 		      if(chi20<Chi2[NPT0-1]){
// 			  Chi2[NPT0-1]=chi20;
		      if((chi20+(TMath::Abs(t00.first)/25.))<Chi2[NPT0-1]){
                  
                  Chi2[NPT0-1]=chi20+(TMath::Abs(t00.first)/25.);
                  chi2=chi20;
                  NPT=NPT0;
                  for(int i=0;i<8;i++) cc_fin[i]=cc[i];
                  for(int i=0;i<nrPnts;i++) ac_fin[i]=ac[i];
                  //for(int i=0;i<nrPnts;i++) ah_fin[i]=ah[i];
                  for(int i=0;i<8;i++) hh_fin[i]=hh[i];
                  m.first=m0.first;
                  q.first=q0.first;
                  t0.first=t00.first;
                  m.second=m0.second;
                  q.second=q0.second;
                  t0.second=t00.second;
                  nrMinPnts0 = TMath::Max(nrMinPnts,(NPT-1));
                  okFit=okFit0;

                  //chi2distribution->Fill(chi2);
		      }
		    }
		  
		  comb++;
		  
		}		
	      }     
	    }	    
	  }	  
	}	
      }
    }    
  }
  
  if(change_code==true) {
    if(DEBUG_SEG){
      for(int i=0;i<8;i++) printf("cc[%d]=%d ",i,cc_fin[i]);  
      printf("\n");  
    }
    
    int jj=0;
    for(int i=0;i<8;i++){
      int ccold=cc_fin[i];
      if(ccold!=0){
	cc_fin[i]=ac_fin[jj];
	jj++;
      }
    }
    
    if(DEBUG_SEG){
      for(int i=0;i<8;i++) printf("cc[%d]=%d ",i,cc_fin[i]);  
      printf("\n");  
      for(int i=0;i<nrPnts;i++) printf("ac[%d]=%d ",i,ac_fin[i]);    
      printf("\n\n");
    }  
    
    Set_HITCode(cc_fin,hh_fin);
    RejectHits0(hit_seg);
    
  }
  
  if(okFit)
    {  // TODO: controllare!!!
      Set_Slope(m); 
      Set_X0(q); 
      
//       if( TMath::Abs(t0+T0_shift)>T0_max ){
// 	t0=0.;
//       }
      if(!change_code) 
	Set_T0(t0);
      else 
	Set_T0_fin(t0);
      
      Set_NPT(NPT);
      Set_Chi2(chi2);
      Set_FitOK(okFit);
    }
  
  if(DEBUG_SEG){
    double tt0=-999.;
    if(!change_code) tt0=Get_T0();
    else tt0=Get_T0_fin();
    printf("\nokFit = %d, Slope = %.3f, X0 = %.1f, t0 = %.1f, Chi2 = %.3f\n", Get_FitOK(), Get_Slope(),Get_X0(),tt0,Get_Chi2());
    for(int i=0;i<hit_seg->Get_NumberHITS();i++){
      HIT *hit=hit_seg->hit(i);
      cout<<"N.hit "<< i <<": left = "<<hit->code_ID().first<<", right = "<<hit->code_ID().second<<endl; 
    }   
  }
  
  for(int i=0;i<8;i++) {  
    delete layer[i]; layer[i]=NULL;
  }
  
  delete [] NX_xLs;
  NX_xLs=NULL;
  
  return;
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Set or REJECT (TODO) nel caso di fit con t0!!!
void Segment::Set_HITCode(int *cc, int *hh){
  for(int i=0;i<8;i++){
    for(int j=0;j<layer[i]->Get_NumberHITS();j++){
      HIT * hit=layer[i]->hit(j);
      {
	if(hh[i] == (j+1)){
	  if( cc[i] == -1 )
	    {
	      hit->Change_LCode(true);
	      hit->Change_RCode(false);
	    }
	  if( cc[i] == 1 )
	    {
	      hit->Change_LCode(false);
	      hit->Change_RCode(true);
	    }
	  if( cc[i] == 0 )
	    {
	      hit->Change_LCode(false);
	      hit->Change_RCode(false);
	    }
	}
      }
      if( hit->IsDouble() ){
	if(hh[i] != (j+1)){
	  hit->Change_LCode(false);
	  hit->Change_RCode(false);
	}
      }
    }
  }
  
  return;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


HITColl_Seg * Segment::RejectHits0(HITColl_Seg *hit_seg){
  
  if(DEBUG_SEG)
    printf(" Starting Segment::RejectHits0 \n");
  
  int N=hit_seg->Get_NumberHITS();
  for(int i=(N-1);i>=0;i--){
    HIT * hit=hit_seg->hit(i);
    if(hit->code_ID().first==0 && hit->code_ID().second==0){
      if(DEBUG_SEG)
      printf("deleting hit %d\n",i);
      hit_seg->eraseHIT(i);
    }
  }  
  
  return hit_seg;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


bool Segment::Set_IsDouble(bool isdouble){
  seg_isdouble=isdouble;
  return seg_isdouble;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


bool Segment::Set_IsGood(bool isgood){
  seg_isgood=isgood;
  return seg_isgood;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


bool Segment::IsGood(){
  return seg_isgood;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


bool Segment::IsDouble(){
  return seg_isdouble;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


double Segment::Set_Slope(pair<double,double> slope){
  seg_Slope=slope.first;
  return seg_Slope;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


double Segment::Get_Slope(){
  
  return seg_Slope;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


double Segment::Set_X0(pair<double,double> X0){
  seg_X0=X0.first;
  return seg_X0;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


double Segment::Get_X0(){
  
  return seg_X0;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


double Segment::Set_T0(pair<double,double> t0){
  seg_T0=t0.first;
  return seg_T0;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


double Segment::Get_T0(){
  
  return seg_T0;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


double Segment::Set_T0_fin(pair<double,double> t0_fin){
  seg_T0_fin=t0_fin.first;
  return seg_T0_fin;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


double Segment::Get_T0_fin(){
  
  return seg_T0_fin;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


int Segment::Set_NPT(int npt){
  seg_NPT=npt;
  return seg_NPT;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


int Segment::Get_NPT(){
  
  return seg_NPT;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


bool Segment::Set_FitOK(bool okFit){
  seg_okFit=okFit;
  return seg_okFit;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


bool Segment::Get_FitOK(){
  
  return seg_okFit;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


double Segment::Set_Chi2(double chi2){
  seg_chi2=chi2;
  return seg_chi2;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


double Segment::Get_Chi2(){
  
  return seg_chi2;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


int * Segment::Get_DoubleHITxL(HITColl_Seg *hit_seg,int nL){
  
  for(int i=0;i<hit_seg->Get_NumberHITS();i++){
    
    HIT *hit=hit_seg->hit(i);
    for(int jj=0;jj<4;jj++){
      if(hit->SL_ID()==1 || hit->SL_ID()==2 )    
	if(hit->L_ID()==(jj+1))
	  L_2hit[jj]++;
      if(hit->SL_ID()==3)    
	if(hit->L_ID()==(jj+1))
	  L_2hit[jj+4]++;
    }

  }
  
  if(DEBUG_SEG)  
    for(int i=0;i<nL;i++)cout<<"L_2hit["<<i<<"] "<<L_2hit[i]<<endl;
  
  return L_2hit;
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


int * Segment::Get_NX_xLs(HITColl_Seg *hit_seg,int nL){

  for(int i=0;i<hit_seg->Get_NumberHITS();i++){
    
    HIT *hit=hit_seg->hit(i);
    for(int jj=0;jj<4;jj++){
      if(hit->SL_ID()==1 || hit->SL_ID()==2 )    
	if(hit->L_ID()==(jj+1)){
	  if(hit->code_ID().first==true)	  
	    NX_xLs[jj]++;
	  if(hit->code_ID().second==true)	  
	    NX_xLs[jj]++;
	}
      if(hit->SL_ID()==3)    
	if(hit->L_ID()==(jj+1)){
	  if(hit->code_ID().first==true)	  
	    NX_xLs[jj+4]++;
	  if(hit->code_ID().second==true)	  
	    NX_xLs[jj+4]++;
	}
    }
  }
  
  if(DEBUG_SEG)  
    for(int i=0;i<nL;i++)cout<<"NX_xLs["<<i<<"] "<<NX_xLs[i]<<endl;
  
  return NX_xLs;
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void Segment::printSegment()
{
  
}
