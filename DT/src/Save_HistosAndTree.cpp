#include "Save_HistosAndTree.h"

static const bool calcola_spline=1;

Save_HistosAndTree::Save_HistosAndTree(){
  count_graph=0;
  m_nmaxseg = 6;
  m_nmaxseg_glo = 4;

  return;
}

Save_HistosAndTree::~Save_HistosAndTree(){
  return;
}

void Save_HistosAndTree::initHB(int runID, int numEvent, ofstream *HBFile){
  
  
  time_t in_time;
  time(&in_time);
  
  // create run header, and write header

  runHEADER runHeader;
  runHeader.number = runID;
  sprintf(runHeader.type,"DATA");
  sprintf(runHeader.time,ctime(&in_time));
  sprintf(runHeader.bor,"BOR"); 
  runHeader.nChamber = 4;
  
  HBFile->write((char *)&runHeader,sizeof(runHeader));
  if(DEBUG_HB){
    cout << " runHeader: size " << sizeof(runHeader) << endl;
    cout << "            number " << runHeader.number << endl;
    cout << "            type " << runHeader.type << endl;
    cout << "            bor " << runHeader.bor << endl;
    cout << "            nCH " << runHeader.nChamber << endl;
  }

  return;
  
}


void Save_HistosAndTree::dumpHB(Track *track, HITCollection *hits,int numEvent,ofstream *HBFile){
  
  if(DEBUG_HB) 
    printf("Into dumpHB...\n");
  
  if(track->Track_IsGood()){
    //     if(track->Get_IsGood_glo(0) && track->Get_IsGood_glo(1)
    //        && track->Get_IsGood_glo(2) && track->Get_IsGood_glo(3)){
    if(track->Get_IsGood_glo(0) && track->Get_IsGood_glo(1)
       && ( (track->Get_IsGood_glo(2) && track->Get_IsGood_glo(3)) 
	    || (track->Get_SegIsAbsorbed(2) && track->Get_SegIsAbsorbed(3)) ) ){
      
      if(DEBUG_HB) 
	printf("If track->Get_IsGood(0) && track->Get_IsGood(1)...\n");
      
      char boe[4];
      sprintf(boe,"BOE");
      HBFile->write((char *)&boe,sizeof(boe));      
      if(DEBUG_HB)
	cout << " Begin of event " << boe << endl;
      
      evtHEADER evtHeader;
      evtHeader.number = numEvent;
      sprintf(evtHeader.type,"DATA");
      for(int i=0;i<4;i++) 
	evtHeader.nHits[i]=0;
      //       if(track->Get_NPT_glo(0)!=-999 && track->Get_NPT_glo(1)!=-999)
      // 	evtHeader.nHits[0] = int(track->Get_NPT_glo(0)+track->Get_NPT_glo(1));
      //       if(track->Get_NPT_glo(2)!=-999 && track->Get_NPT_glo(3)!=-999)
      // 	evtHeader.nHits[1] = int(track->Get_NPT_glo(2)+track->Get_NPT_glo(3));
      //       if(track->Get_NPT(4)!=-999)
      // 	evtHeader.nHits[2] = int(track->Get_NPT(4));
      //       if(track->Get_NPT(5)!=-999)
      // 	evtHeader.nHits[3] = int(track->Get_NPT(5));
      HBFile->write((char *)&evtHeader,sizeof(evtHeader));
      if(DEBUG_HB){
	cout << " evtHeader: size " << sizeof(evtHeader) << endl;
	cout << "            number " << evtHeader.number << endl;
	cout << "            type " << evtHeader.type << endl;
	for(int i=0;i<4;i++) 
	  printf("            nHits[%d] %d\n",i,evtHeader.nHits[i]);
      }
      
      trkHEADER trkHeader; 
      trkHeader.nTracks = 2;
      if(track->Get_IsGood_glo(2) || track->Get_SegIsAbsorbed(2) )
	trkHeader.nTracks++;
      if(track->Get_IsGood_glo(3) || track->Get_SegIsAbsorbed(3) ) 
	trkHeader.nTracks++;
      if(track->Get_IsGood(4) || track->Get_SegIsAbsorbed(4) ) 
	trkHeader.nTracks++;
      if(track->Get_IsGood(5) || track->Get_SegIsAbsorbed(5) ) 
	trkHeader.nTracks++;

      // Writing tracks header (how many tracks do I expect when reading)
      HBFile->write((char *)&trkHeader,sizeof(trkHeader));
      if(DEBUG_HB){
	cout << " trkHeader: size " << sizeof(trkHeader) << endl;
	cout << "            ntracks " << trkHeader.nTracks << endl;
      }

      TRACK_HB track_HB;
      
      int nlay=8; 
      double res_glo[m_nmaxseg][nlay];
      int seg=0; 
      for(int jj=0;jj<m_nmaxseg;jj++)
	for(int j=0;j<nlay;j++){
	  res_glo[jj][j]=-999.;  
	}
      int jres=0;
      
      for(int ih=0;ih<hits->Get_NumberHITS();ih++){
	HIT *hit=hits->hit(ih);
	if(hit->code_ID().first==1 || hit->code_ID().second==1){
	  if( (hit->CH_ID()==11 || hit->CH_ID()==10) && hit->SL_ID()==3 )
	    jres=4+hit->L_ID()-1;
	  else 
	    jres=hit->L_ID()-1;
	  
	  if(hit->CH_ID()==11)
	    if(hit->SL_ID()!=2) 
	      seg=0;
	    else
	      seg=1;
	  if(hit->CH_ID()==10)
	    if(hit->SL_ID()!=2) 
	      seg=2;
	    else
	      seg=3;
	  if(hit->CH_ID()==9)
	    seg=4;
	  if(hit->CH_ID()==8)  
	    seg=5;
	  
	  if(seg<4){
	    if(hit->code_ID().first==1) 
	      res_glo[seg][jres]= hit->Get_XPair(hit->dtime()).first - (-1*0.00547*track->Get_T0_glo(seg)+track->Get_Slope_glo(seg)*(hit->y_wire_ID())+track->Get_X0_glo(seg));
	    if(hit->code_ID().second==1) 
	      res_glo[seg][jres]= hit->Get_XPair(hit->dtime()).second - (1*0.00547*track->Get_T0_glo(seg)+track->Get_Slope_glo(seg)*(hit->y_wire_ID())+track->Get_X0_glo(seg));
	  }
	  else{
	    if(hit->code_ID().first==1) 
	      res_glo[seg][jres]= hit->Get_XPair(hit->dtime()).first - (track->Get_Slope(seg)*(hit->y_wire_ID())+track->Get_X0(seg));
	    if(hit->code_ID().second==1) 
	      res_glo[seg][jres]= hit->Get_XPair(hit->dtime()).second - (track->Get_Slope(seg)*(hit->y_wire_ID())+track->Get_X0(seg));
	  }
	  
	}
	
      } // close loop on hits
      
      
      for(int i=0; i<m_nmaxseg; i++){
	track_HB.chamber=-999;
	track_HB.nPoints = -999;
	track_HB.XorZ=-999;
	track_HB.slope=-999.;
	track_HB.inter=-999.;
	track_HB.erSlope=-999.;
	track_HB.erInter=-999.;
	track_HB.erCorr=-999.;
	track_HB.tZero=-999.;
	track_HB.chi2=-999.;
	for(int j=0;j<nlay;j++) {
	  track_HB.res[j]=-999.;
	  track_HB.lay[j]=-999;
	}
	
	
	if(i<4 && ( track->Get_IsGood_glo(i) || track->Get_SegIsAbsorbed(i) ) ){
	  
	  if(i==0 || i==1) 
	    track_HB.chamber=0;
	  else 
	    track_HB.chamber=1;

	  if(track->Get_SegIsAbsorbed(i))
	    track_HB.nPoints = -1;
	  else
	    track_HB.nPoints = track->Get_NPT_glo(i);	    
	  
	  if(i==0 || i==2)
	    track_HB.XorZ=0;
	  else
	    track_HB.XorZ=2;
	  
	  track_HB.slope=track->Get_Slope_glo(i);
	  track_HB.inter=track->Get_X0_glo(i);
	  track_HB.erSlope=track->Get_erSlope_glo(i);
	  track_HB.erInter=track->Get_erX0_glo(i);
	  track_HB.tZero=track->Get_T0_glo_fin(i);
	  track_HB.chi2=track->Get_Chi2_glo(i);
	  for(int j=0;j<nlay;j++) {
	    track_HB.res[j]=res_glo[i][j];
	    track_HB.lay[j]=j;
	  }
	  
	  if(DEBUG_HB){
	    cout << "\n Track: chamber " << track_HB.chamber <<
	      " nPoints " << track_HB.nPoints << 
	      " XorZ " << track_HB.XorZ <<
	      " slope " << track_HB.slope <<
	      " inter " << track_HB.inter <<
	      " erSlope " << track_HB.erSlope << 
	      " erInter " << track_HB.erInter <<
	      " tZero " << track_HB.tZero <<
	      " chi2 " << track_HB.chi2 <<endl;
	    cout << "        res " << track_HB.res[0] << " " <<  track_HB.res[1] << " " <<  
	      track_HB.res[2] << " " <<  track_HB.res[3] << " " <<
	      track_HB.res[4] << " " <<  track_HB.res[5] << " " <<  
	      track_HB.res[6] << " " <<  track_HB.res[7] << " " << 
	      endl;	
	    cout << "        lay " << track_HB.lay[0] << " " <<  track_HB.lay[1] << " " <<  
	      track_HB.lay[2] << " " <<  track_HB.lay[3] << " " <<
	      track_HB.lay[4] << " " <<  track_HB.lay[5] << " " <<  
	      track_HB.lay[6] << " " <<  track_HB.lay[7] << " " << 
	      endl;	
	  
	  }
	  HBFile->write((char *)&track_HB,sizeof(track_HB));
	
	} // close if<4
	
	else if( track->Get_IsGood_glo(0) && track->Get_IsGood_glo(1) 
		 && ( (track->Get_IsGood_glo(2) && track->Get_IsGood_glo(3)) 
		      || (track->Get_SegIsAbsorbed(2) && track->Get_SegIsAbsorbed(3)) )
		 && (i>=4 && (track->Get_IsGood(i) || track->Get_SegIsAbsorbed(i)) ) ){
	  
 	  
	  if(i==4) 
	    track_HB.chamber=2;
	  else if(i==5) 
	    track_HB.chamber=3;
	  if(track->Get_SegIsAbsorbed(i))
	    track_HB.nPoints = -1;
	  else
	    track_HB.nPoints = track->Get_NPT(i);
	  
	  track_HB.XorZ=0;
	  track_HB.slope=track->Get_Slope(i);
	  track_HB.inter=track->Get_X0(i);
	  track_HB.erSlope=-999.;
	  track_HB.erInter=-999.;
	  track_HB.tZero=-999.;
	  track_HB.chi2=track->Get_Chi2(i);
	  for(int j=0;j<nlay;j++) {
	    track_HB.res[j]=res_glo[i][j];
	    track_HB.lay[j]=j;
	  }
	  
	  if(DEBUG_HB){
	    cout << "\n Track: chamber " << track_HB.chamber <<
	      " nPoints " << track_HB.nPoints << 
	      " XorZ " << track_HB.XorZ <<
	      " slope " << track_HB.slope <<
	      " inter " << track_HB.inter <<
	      " erSlope " << track_HB.erSlope << 
	      " erInter " << track_HB.erInter <<
	      " tZero " << track_HB.tZero <<
	      " chi2 " << track_HB.chi2 <<endl;
	    cout << " res " << track_HB.res[0] << " " <<  track_HB.res[1] << " " <<  
	      track_HB.res[2] << " " <<  track_HB.res[3] << " " <<
	      track_HB.res[4] << " " <<  track_HB.res[5] << " " <<  
	      track_HB.res[6] << " " <<  track_HB.res[7] << " " << 
	      endl;	
	    cout << " lay " << track_HB.lay[0] << " " <<  track_HB.lay[1] << " " <<  
	      track_HB.lay[2] << " " <<  track_HB.lay[3] << " " <<
	      track_HB.lay[4] << " " <<  track_HB.lay[5] << " " <<  
	      track_HB.lay[6] << " " <<  track_HB.lay[7] << " " << 
	      endl;	
	  }
	  
	  HBFile->write((char *)&track_HB,sizeof(track_HB));
	  
	}
	
	if(DEBUG_HB) 
	  printf("All variables filled...\ni=%d\n",i);
	
      }
      
      
      if(DEBUG_HB) 
	printf("All variables filled...\n");
      
      if(DEBUG_HB) 
	printf("Exiting dumpHB...\n");
      
    } // close if(track->Get_IsGood_glo(0) && track->Get_IsGood_glo(1))
    
  } // close if(Track_IsGood()...)
  
  
  if(DEBUG_HB) 
    printf("Exiting dumpHB...\n");
  
  return;
}


void Save_HistosAndTree::dumpTree(Track *track, HITCollection *hits,int numEvent,TTree *tree){
  
    if(DEBUG_TREE)
      printf("Start filling tree...\n");

    cleanTree();

    if(track->Track_IsGood()){
        if(track->Get_IsGood(0) && track->Get_IsGood(1)){
            dumpTree_Hits(hits,numEvent,tree);
            dumpTree_Track(track,hits,numEvent,tree);
        }
    }

    onevent = numEvent; //ALTEA
    tree->Fill();
  
    return;
}


void Save_HistosAndTree::dumpTree_Hits(HITCollection *hits,int numEvent,TTree *tree){

  if(DEBUG_STOREHIT) 
    cout<<"N.hit: "<<hits->Get_NumberHITS()<<endl;
  
  if(DEBUG_STOREHIT) 
    if(hits->Get_NumberHITS()<200 && numEvent<11 ){
      printf("\nEvent %d, N.hits %d\n", numEvent, hits->Get_NumberHITS());
      for(int i=0;i<hits->Get_NumberHITS();i++){
	HIT *hit=hits->hit(i);
	printf("Hit %d: CH=%d, SL=%d, L=%d, W=%d, dx=%.1f, xW=%.1f, yW=%.1f, t=%.0f, dt=%.0f\n",i,hit->CH_ID(),hit->SL_ID(),hit->L_ID(),hit->wire_ID(),(hit->dtime())*0.00547,hit->x_wire_ID(),hit->y_wire_ID(),hit->rtime(),hit->dtime());
      }
    }
  
  onhit=hits->Get_NumberHITS();
  if(DEBUG_STOREHIT)
    printf("onhits = %d\n",onhit);
    
  for(int i=0; i<onhit; i++){
    if(DEBUG_STOREHIT)
      printf("Loop on hit\n");
    HIT *hit=hits->hit(i);
    ohtrig[i] = int( hit->rtime() - hit->dtime_in() );
    int hlay=0;  //---->NTUPLA
    if(hit->CH_ID()==11){
      int lay = 0;
      if(hit->SL_ID()==1) lay = hit->L_ID();
      if(hit->SL_ID()==2) lay = 4 + hit->L_ID();
      if(hit->SL_ID()==3) lay = 8 + hit->L_ID();
      hlay = 1300 + lay;
    }
    if(hit->CH_ID()==10){
      int lay = 0;
      if(hit->SL_ID()==1) lay = hit->L_ID();
      if(hit->SL_ID()==2) lay = 4 + hit->L_ID();
      if(hit->SL_ID()==3) lay = 8 + hit->L_ID();
      hlay = 2300 + lay;
    }
    if(hit->CH_ID()==9){
      int lay = 0;
      if(hit->SL_ID()==3) lay = hit->L_ID();
      hlay = 3300 + lay;
    }
    if(hit->CH_ID()==8){
      int lay = 0;
      if(hit->SL_ID()==1) lay = hit->L_ID();
      hlay = 4300 + lay;
    }
    ohlay[i] = hlay;
    ohwire[i] = hit->wire_ID();
    ohtime_in[i] = hit->dtime_in();
    ohtime[i] = hit->dtime();
    
    if(hit->CH_ID()==9 && hit->SL_ID()==3 && (((hit->L_ID()==1||hit->L_ID()==3) && (hit->wire_ID()>71 && hit->wire_ID()<77)) || ((hit->L_ID()==2||hit->L_ID()==4) && (hit->wire_ID()>71 && hit->wire_ID()<76))) && hit->dtime_in()>=-10){
      int nfilo=0;
      int fili=18;
      if(hit->L_ID()==1){
	if(hit->wire_ID()==72) nfilo=0;
	if(hit->wire_ID()==73) nfilo=1;
	if(hit->wire_ID()==74) nfilo=2;
	if(hit->wire_ID()==75) nfilo=3;	
	if(hit->wire_ID()==76) nfilo=4;
      }
      if(hit->L_ID()==2){
	if(hit->wire_ID()==72) nfilo=5;
	if(hit->wire_ID()==73) nfilo=6;
	if(hit->wire_ID()==74) nfilo=7;
	if(hit->wire_ID()==75) nfilo=8;	
      }
      if(hit->L_ID()==3){
	if(hit->wire_ID()==72) nfilo=9;
	if(hit->wire_ID()==73) nfilo=10;
	if(hit->wire_ID()==74) nfilo=11;
	if(hit->wire_ID()==75) nfilo=12;	
	if(hit->wire_ID()==76) nfilo=13;
      }
      if(hit->L_ID()==4){
	if(hit->wire_ID()==72) nfilo=14;
	if(hit->wire_ID()==73) nfilo=15;
	if(hit->wire_ID()==74) nfilo=16;
	if(hit->wire_ID()==75) nfilo=17;	
      }
      
      
      if(ohtime_tube[nfilo]==-999){
	ohtime_tube[nfilo] = hit->rtime();
	oNhit_tube[nfilo]++;
	}
	else 
	  if(ohtime_tube[nfilo]!=-999 && ohtime_tube[nfilo+fili*1]==-999 && ohtime_tube[nfilo+fili*2]==-999 && ohtime_tube[nfilo+fili*3]==-999 && ohtime_tube[nfilo+fili*4]==-999 && ohtime_tube[nfilo+fili*5]==-999){
	    ohtime_tube[nfilo+fili*1] = hit->rtime(); 
	    oNhit_tube[nfilo]++;
	  }	    
	  else 
	    if(ohtime_tube[nfilo]!=-999 && ohtime_tube[nfilo+fili*1]!=-999 && ohtime_tube[nfilo+fili*2]==-999 && ohtime_tube[nfilo+fili*3]==-999 && ohtime_tube[nfilo+fili*4]==-999 && ohtime_tube[nfilo+fili*5]==-999){
	      ohtime_tube[nfilo+fili*2] = hit->rtime(); 
	      oNhit_tube[nfilo]++;
	    }
	    else 
	      if(ohtime_tube[nfilo]!=-999 && ohtime_tube[nfilo+fili*1]!=-999 && ohtime_tube[nfilo+fili*2]!=-999 && ohtime_tube[nfilo+fili*3]==-999 && ohtime_tube[nfilo+fili*4]==-999 && ohtime_tube[nfilo+fili*5]==-999){
		ohtime_tube[nfilo+fili*3] = hit->rtime(); 
		oNhit_tube[nfilo]++;
	      }
	      else
		if(ohtime_tube[nfilo]!=-999 && ohtime_tube[nfilo+fili*1]!=-999 && ohtime_tube[nfilo+fili*2]!=-999 && ohtime_tube[nfilo+fili*3]!=-999 && ohtime_tube[nfilo+fili*4]==-999 && ohtime_tube[nfilo+fili*5]==-999){
		  ohtime_tube[nfilo+fili*4] = hit->rtime(); 
		  oNhit_tube[nfilo]++;
		}
		else
		  if(ohtime_tube[nfilo]!=-999 && ohtime_tube[nfilo+fili*1]!=-999 && ohtime_tube[nfilo+fili*2]!=-999 && ohtime_tube[nfilo+fili*3]!=-999 && ohtime_tube[nfilo+fili*4]!=-999 && ohtime_tube[nfilo+fili*5]==-999){
		    ohtime_tube[nfilo+fili*5] = hit->rtime(); 
		    oNhit_tube[nfilo]++;
		  }
    }
  }
  
  
  
  if(DEBUG_STOREHIT)
    printf("onevent = %d\n",numEvent);
  onevent = numEvent;
  
  return;
  
}

void Save_HistosAndTree::dumpTree_Track(Track *track, HITCollection *hits,int numEvent,TTree *tree){
  
  if(DEBUG_STORETRACK) 
    printf("Into dumpTree_Track...\n");

  // check track
  if(!track->Track_IsGood()) return;

  // track loop
  for(int it=0; it<m_nmaxseg; it++){

      if(!track->Get_X0(it)) return;

      //old-nasty code...
      int k=-999;
      if(it==0 && track->Get_IsGood(0) ) k= 1300+track->Get_NPT(it);
      if(it==1 && track->Get_IsGood(1) ) k=-1300-track->Get_NPT(it);
      if(it==2 && track->Get_IsGood(2) ) k= 2300+track->Get_NPT(it);
      if(it==3 && track->Get_IsGood(3) ) k=-2300-track->Get_NPT(it);
      if(it==4 && track->Get_IsGood(4) ) k= 3300+track->Get_NPT(it);
      if(it==5 && track->Get_IsGood(5) ) k= 4300+track->Get_NPT(it);

      float slope = track->Get_Slope(it);
      float X0 = track->Get_X0(it);
      float T0 = track->Get_T0(it);
      float chi2 = track->Get_Chi2(it);

      // compute residuals
      double res[8]={-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999};
      double xhit[8]={-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999};
      for(int ih=0;ih<hits->Get_NumberHITS();ih++){
          HIT *hit=hits->hit(ih);
          int hlay = hit->L_ID();
          int hsl = hit->SL_ID();
          int hch = hit->CH_ID();

          // reject hits not on the segment
          bool rejectHit = true;
          if(it==0 && hch==11 && (hsl==1 || hsl==3)) rejectHit = false;
          if(it==1 && hch==11 && hsl==2) rejectHit = false;
          if(it==2 && hch==10 && (hsl==1 || hsl==3)) rejectHit = false;
          if(it==3 && hch==10 && hsl==2) rejectHit = false;
          if(it==4 && hch==9) rejectHit = false;
          if(it==5 && hch==8) rejectHit = false;
          if(hit->code_ID().first==0 && hit->code_ID().second==0) rejectHit = true;

          if(rejectHit) continue;

          float hX = -999.;
          if(hit->code_ID().first==1) hX = hit->Get_XPair(hit->dtime()).first;
          if(hit->code_ID().second==1) hX = hit->Get_XPair(hit->dtime()).second;

          int jres = (hsl==3)? hlay-1+4 : hlay-1;

          // SV 20170724 for LEMMA tb fill with hX
          res[jres]= hX - (slope*(hit->y_wire_ID())+X0);
          xhit[jres]=hX;
      }// end hit loop

      // fill....
      fillVar(k,it,slope,X0,T0,chi2,res,it,xhit);

//      /// debugging
//      std::cout << "*** RESIDUALS :";
//      for(int ir=0; ir<8; ir++)
//          std::cout << "   " << res[ir];
//      std::cout << std::endl;

//      std::cout << "slope " << slope << std::endl;
//      std::cout << "X " << X0 << std::endl;
  }// end track loop
	
  if(DEBUG_STORETRACK)
      printf("All variables filled...\n");

  onevent = numEvent;
  onseg = m_nmaxseg;

  if(DEBUG_STORETRACK)
    printf("Exiting dumpTree_Track...\n");
  
  return;
}

/// SV 20170725 this is the old BUGGY code....
//void Save_HistosAndTree::dumpTree_Track(Track *track, HITCollection *hits,int numEvent,TTree *tree){

//    if(DEBUG_STORETRACK)
//        printf("Into dumpTree_Track...\n");

//    if(track->Track_IsGood()){
//        if(track->Get_IsGood(0) && track->Get_IsGood(1)){

//            if(DEBUG_STORETRACK)
//                printf("If track->Get_IsGood(0) && track->Get_IsGood(1)...\n");

//            int inseg=0;
//            for(int i=0; i<m_nmaxseg; i++){

//                int nlay=0;
//                if(i==0 || i==2) nlay=8;
//                else nlay=4;
//                double res[nlay]; for(int j=0;j<nlay;j++) res[j]=-999.;
//                double res_glo[nlay]; for(int j=0;j<nlay;j++) res_glo[j]=-999.;
//                int jres=0;
//                for(int ih=0;ih<hits->Get_NumberHITS();ih++){
//                    HIT *hit=hits->hit(ih);
//                    if(hit->code_ID().first==1 || hit->code_ID().second==1){
//                        if( (hit->CH_ID()==11 || hit->CH_ID()==10) && hit->SL_ID()==3 )
//                            jres=4+hit->L_ID()-1;
//                        else
//                            jres=hit->L_ID()-1;

//                        if(track->Get_IsGood(i)){
//                            if(hit->code_ID().first==1)
//                                res[jres]= hit->Get_XPair(hit->dtime()).first - (track->Get_Slope(i)*(hit->y_wire_ID())+track->Get_X0(i));
//                            if(hit->code_ID().second==1)
//                                res[jres]= hit->Get_XPair(hit->dtime()).second - (track->Get_Slope(i)*(hit->y_wire_ID())+track->Get_X0(i));

//                            std::cout << "\n Track " << i << "---> Res[" << jres << "]=" << res[jres] << std::endl;
//                            std::cout << "X left " << hit->Get_XPair(hit->dtime()).first << ", X right " <<  hit->Get_XPair(hit->dtime()).second;
//                            std::cout << " ---> X fit " << (track->Get_Slope(i)*(hit->y_wire_ID())+track->Get_X0(i)) << std::endl;
//                            std::cout << "Segment " << i << ", slope " << track->Get_Slope(i) << ", X " << track->Get_X0(i) << std::endl;
//                        }// close good track

//                        if(track->Get_IsGood_glo(i)){
//                            if(hit->code_ID().first==1)
//                                res_glo[jres]= hit->Get_XPair(hit->dtime()).first - (-1*0.00547*track->Get_T0_glo(i)+track->Get_Slope_glo(i)*(hit->y_wire_ID())+track->Get_X0_glo(i));
//                            if(hit->code_ID().second==1)
//                                res_glo[jres]= hit->Get_XPair(hit->dtime()).second - (1*0.00547*track->Get_T0_glo(i)+track->Get_Slope_glo(i)*(hit->y_wire_ID())+track->Get_X0_glo(i));
//                        }
//                    }// close ID
//                } // close loop on hits


//                if(DEBUG_STORETRACK)
//                    printf("computing k (NPT)...\n");

//                int k=-999;
//                if(i==0 && track->Get_IsGood(0) ) k= 1300+track->Get_NPT(i);
//                if(i==1 && track->Get_IsGood(1) ) k=-1300-track->Get_NPT(i);
//                if(i==2 && track->Get_IsGood(2) ) k= 2300+track->Get_NPT(i);
//                if(i==3 && track->Get_IsGood(3) ) k=-2300-track->Get_NPT(i);
//                if(i==4 && track->Get_IsGood(4) ) k= 3300+track->Get_NPT(i);
//                if(i==5 && track->Get_IsGood(5) ) k= 4300+track->Get_NPT(i);

//                if(i<m_nmaxseg){
//                    fillVar(k,inseg,track->Get_Slope(i),track->Get_X0(i),track->Get_T0(i),track->Get_Chi2(i),res,i);

//                    /// debugging
//                    std::cout << "*** RESIDUALS :";
//                    for(int ir=0; ir<nlay; ir++)
//                        std::cout << "   " << res[ir];
//                    std::cout << std::endl;

//                    std::cout << "slope " << track->Get_Slope(i) << std::endl;
//                    std::cout << "X " << track->Get_X0(i) << std::endl;


//                }

//                if(DEBUG_STORETRACK)
//                    printf("computing k_glo (NPT_glo)...\n");

//                int k_glo=-999;
//                if(i==0 && track->Get_IsGood_glo(0) ) k_glo= 1300+track->Get_NPT_glo(i);
//                if(i==1 && track->Get_IsGood_glo(1) ) k_glo=-1300-track->Get_NPT_glo(i);
//                if(i==2 && track->Get_IsGood_glo(2) ) k_glo= 2300+track->Get_NPT_glo(i);
//                if(i==3 && track->Get_IsGood_glo(3) ) k_glo=-2300-track->Get_NPT_glo(i);

//                if(i<m_nmaxseg_glo){
//                    fillVar_glo(k_glo,inseg,track->Get_Slope_glo(i),track->Get_erSlope_glo(i),track->Get_X0_glo(i),track->Get_erX0_glo(i),track->Get_T0_glo_fin(i),track->Get_Chi2_glo(i),res_glo,i);
//                }

//                if(DEBUG_STORETRACK)
//                    printf("All variables filled...\ninseg=%d\n",inseg);

//                inseg++;
//            }

//            if(DEBUG_STORETRACK)
//                printf("All variables filled...\n");

//            onevent = numEvent;
//            onseg = inseg;
//            onseg_glo = inseg;

//            if(DEBUG_STORETRACK)
//                printf("Exiting dumpTree_Track...\n");

//        } // close if(track->Get_IsGood(0) && track->Get_IsGood(1))
//    } // close if(Track_IsGood()...)

//    if(DEBUG_STORETRACK)
//        printf("Exiting dumpTree_Track...\n");

//    return;
//}


void Save_HistosAndTree::dumpHisto(Track *track, HITCollection *hits,int numEvent){

  if(track->Track_IsGood()){
    if(track->Get_IsGood(0) && track->Get_IsGood(1)){
      if(track->Get_IsGood(2) && track->Get_IsGood(3)){
	h_dphi->Fill(TMath::ATan(track->Get_Slope(2))-TMath::ATan(track->Get_Slope(0)));
	h_dthe->Fill(TMath::ATan(track->Get_Slope(3))-TMath::ATan(track->Get_Slope(1)));
	h_dphi_glo->Fill(TMath::ATan(track->Get_Slope_glo(2))-TMath::ATan(track->Get_Slope_glo(0)));
	h_dthe_glo->Fill(TMath::ATan(track->Get_Slope_glo(3))-TMath::ATan(track->Get_Slope_glo(1)));
	if( track->Get_T0(0)!=0. && track->Get_T0(2)!=0. && 
	    (track->Get_T0(2)-track->Get_T0(0)) != 0.){
	  h_dT0->Fill(track->Get_T0(2)-track->Get_T0(0));
	  h_T0ch1_T0ch2->Fill(track->Get_T0(0),track->Get_T0(2));
	}
	if( track->Get_T0_fin(0)!=0. && track->Get_T0_fin(2)!=0. && 
	    (track->Get_T0_fin(2)-track->Get_T0_fin(0)) != 0.){
	  h_dT0_fin->Fill(track->Get_T0_fin(2)-track->Get_T0_fin(0));
	  if(TMath::Abs(track->Get_T0(2))<30. && TMath::Abs(track->Get_T0(0))<30.)
	    h_dT0_fin2->Fill(track->Get_T0(2)+track->Get_T0_fin(2)-track->Get_T0(0)-track->Get_T0_fin(0));
	  
	}
	if( track->Get_T0_glo_fin(0)!=0. && track->Get_T0_glo_fin(2)!=0. && 
	    (track->Get_T0_glo_fin(2)-track->Get_T0_glo_fin(0)) != 0.){
	  h_dT0_glo->Fill(track->Get_T0_glo_fin(2)-track->Get_T0_glo_fin(0));
	  for(int jj=0;jj<4;jj++){
	    h_dT0_slope_glo[jj]->Fill((track->Get_T0_glo_fin(2)-track->Get_T0_glo_fin(0)),track->Get_Slope_glo(jj));
	    h_dT0_NPT_glo[jj]->Fill((track->Get_T0_glo_fin(2)-track->Get_T0_glo_fin(0)),track->Get_NPT_glo(jj));
	  }
	  h_T0ch1_T0ch2_glo->Fill(track->Get_T0_glo_fin(0),track->Get_T0_glo_fin(2));
	}
      }
      for(int i=0;i<6;i++){
	if(track->Get_IsGood(i)){
	  h_SLOPE[i]->Fill(TMath::ATan(track->Get_Slope(i)));
	  h_X0[i]->Fill(track->Get_X0(i));
	  h_T0[i]->Fill(track->Get_T0(i));
	  h_NPT[i]->Fill(track->Get_NPT(i));
	  h_CHI2[i]->Fill(track->Get_Chi2(i));
	}
      }
    } // close if(track->Get_IsGood(0) && track->Get_IsGood(1))
    
    for(int i=0;i<4;i++){
      if(track->Get_IsGood_glo(i)){
	h_SLOPE_glo[i]-> Fill(TMath::ATan(track->Get_Slope_glo(i)));
	h_erSLOPE_glo[i]-> Fill(TMath::ATan(track->Get_erSlope_glo(i)));
	h_X0_glo[i]->Fill(track->Get_X0_glo(i));
	h_erX0_glo[i]->Fill(track->Get_erX0_glo(i));	
	if(i==0){
	  h_T0_glo[i]->Fill(track->Get_T0_glo_fin(i));
	  h_erT0_glo[i]->Fill(track->Get_erT0_glo(i));
	}
	if(i==2){
	  h_T0_glo[1]->Fill(track->Get_T0_glo_fin(i));
	  h_erT0_glo[1]->Fill(track->Get_erT0_glo(i));
	}
      }
    }
    
  } // close if(track->Track_IsGood())

  if(track->Track_IsGood()){
    if(track->Get_IsGood_glo(0) && track->Get_IsGood_glo(1)){
      
      float sigmaPhi[4]={0.,0.,0.,0.};
      for(int i=0;i<4;i++)
	sigmaPhi[i]=0.03 + 0.025 * (1-TMath::Cos(TMath::ATan(track->Get_Slope_glo(i))));
      
      if(track->Get_NPT_glo(1)==4){
	float xth[4]={-999.,-999.,-999.,-999.};
	for(int ih=0;ih<hits->Get_NumberHITS();ih++){
	  HIT *hit=hits->hit(ih);
	  if( hit->code_ID().first==1 || hit->code_ID().second==1 ){
	    if( hit->CH_ID()==11 && hit->SL_ID()==2 )
	      for(int i=0;i<4;i++){
		if( hit->L_ID()==(i+1) ){
		  if(hit->code_ID().first==1) 
		    xth[i]= hit->Get_XPair(hit->dtime()).first;
		  if(hit->code_ID().second==1) 
		    xth[i]= hit->Get_XPair(hit->dtime()).second;
		  if(DEBUG_STOREHIT)
		    printf("...xth[%d]=%.1f, ay[%d]=%.3f\n",i,xth[i],i,hit->y_wire_ID());
		}
	      }
	  }
	} // close loop on hits

	if(xth[0]>-999. && xth[1]>-999. && xth[2]>-999. && xth[3]>-999.){
	  h_MP_CH_glo[0]->Fill(((xth[0]+xth[2])/2. - xth[1])/(sigmaPhi[1]*1.22));
	  h_MP_CH_glo[0]->Fill(((xth[1]+xth[3])/2. - xth[2])/(sigmaPhi[1]*1.22));
	  h_MP_CH_glo[0]->Fill(((xth[0]+2.*xth[3])/3. - xth[2])/(sigmaPhi[1]*1.247));
	  h_MP_CH_glo[0]->Fill(((2.*xth[0]+xth[3])/3. - xth[1])/(sigmaPhi[1]*1.247));
	  if(DEBUG_STOREHIT)
	    {
	      printf("   MP_012=%.2f\n",((xth[0]+xth[2])/2. - xth[1])/(sigmaPhi[1]*1.22));
	      printf("   MP_123=%.2f\n",((xth[1]+xth[3])/2. - xth[2])/(sigmaPhi[1]*1.22));
	      printf("   MP_023=%.2f\n",((xth[0]+2.*xth[3])/3. - xth[2])/(sigmaPhi[1]*1.247));
	      printf("   MP_023=%.2f\n",((2.*xth[0]+xth[3])/3. - xth[1])/(sigmaPhi[1]*1.247));
	    }
	}
	
	
	if(track->Get_IsGood_glo(2) && track->Get_IsGood_glo(3)){
	  if(track->Get_NPT_glo(3)==4){
	    float x2th[4]={-999.,-999.,-999.,-999.};
	    for(int ih=0;ih<hits->Get_NumberHITS();ih++){
	      HIT *hit=hits->hit(ih);
	      if( hit->code_ID().first==1 || hit->code_ID().second==1 ){
		if( hit->CH_ID()==10 && hit->SL_ID()==2 )
		  for(int i=0;i<4;i++){
		    if( hit->L_ID()==(i+1) ){
		      if(hit->code_ID().first==1) 
			x2th[i]= hit->Get_XPair(hit->dtime()).first;
		      if(hit->code_ID().second==1) 
			x2th[i]= hit->Get_XPair(hit->dtime()).second;
		      if(DEBUG_STOREHIT)
			printf("...x2th[%d]=%.1f, ay[%d]=%.3f\n",i,x2th[i],i,hit->y_wire_ID());
		    }
		  }
	      }
	    } // close loop on hits
	    
	    if(x2th[0]>-999. && x2th[1]>-999. && x2th[2]>-999. && x2th[3]>-999.){
	      h_MP_CH_glo[1]->Fill(((x2th[0]+x2th[2])/2. - x2th[1])/(sigmaPhi[3]*1.22));
	      h_MP_CH_glo[1]->Fill(((x2th[1]+x2th[3])/2. - x2th[2])/(sigmaPhi[3]*1.22));
	      h_MP_CH_glo[1]->Fill(((x2th[0]+2.*x2th[3])/3. - x2th[2])/(sigmaPhi[3]*1.247));
	      h_MP_CH_glo[1]->Fill(((2.*x2th[0]+x2th[3])/3. - x2th[1])/(sigmaPhi[3]*1.247));
	      if(DEBUG_STOREHIT)
		{
		printf("   MP_012=%.2f\n",((x2th[0]+x2th[2])/2. - x2th[1])/(sigmaPhi[3]*1.22));
		printf("   MP_123=%.2f\n",((x2th[1]+x2th[3])/2. - x2th[2])/(sigmaPhi[3]*1.22));
		printf("   MP_023=%.2f\n",((x2th[0]+2.*x2th[3])/3. - x2th[2])/(sigmaPhi[3]*1.247));
		printf("   MP_023=%.2f\n",((2.*x2th[0]+x2th[3])/3. - x2th[1])/(sigmaPhi[3]*1.247));
	      }
	    }
	  }
	  
	} // close if(track->Get_IsGood(2) && track->Get_IsGood(3))
      } // close if(track->Get_NPT_glo(3)==4)
    } // close if(track->Get_IsGood(0) && track->Get_IsGood(1))
  } // close if(track->Track_IsGood())

  if(track->Track_IsGood()){
    if(track->Get_IsGood(0) && track->Get_IsGood(1)){
      
      if(DEBUG_STORETRACK) 
	printf("Filling Histo for res\n");
      
      int inseg=0;
      for(int i=0; i<m_nmaxseg; i++){
	
	int nlay=0; 
	if(i==0 || i==2) nlay=8;
	else nlay=4;
    double res[nlay];
    for(int j=0;j<nlay;j++) res[j]=-999.;
    double res_glo[nlay]; for(int j=0;j<nlay;j++) res_glo[j]=-999.;
	
	// distance track-wire ...to compute linear correction
	double dft[nlay]; 
	double tdrift[nlay]; 
	double xhit[nlay];
	double yhit[nlay];
	
	if(calcola_spline){
	  for(int j=0;j<nlay;j++) dft[j]=-999.; 
	  for(int j=0;j<nlay;j++) tdrift[j]=-999.; 
	  for(int j=0;j<nlay;j++) xhit[j]=-999.; 
	  for(int j=0;j<nlay;j++) yhit[j]=-999.; 
	}


    int jres=0;
    for(int ih=0;ih<hits->Get_NumberHITS();ih++){
      HIT *hit=hits->hit(ih);
      if(hit->code_ID().first==1 || hit->code_ID().second==1){
        if( (hit->CH_ID()==11 || hit->CH_ID()==10) && hit->SL_ID()==3 )
          jres=4+hit->L_ID()-1;
        else
          jres=hit->L_ID()-1;

        if(track->Get_IsGood(i)){
         if(hit->code_ID().first==1)
              res[jres]=hit->Get_XPair(hit->dtime()).first - (track->Get_Slope(i)*(hit->y_wire_ID())+track->Get_X0(i));

            if(hit->code_ID().second==1)
            res[jres]=hit->Get_XPair(hit->dtime()).second - (track->Get_Slope(i)*(hit->y_wire_ID())+track->Get_X0(i));
          if(calcola_spline){
        dft[jres]= hit->x_wire_ID() - (track->Get_Slope(i)*(hit->y_wire_ID())+track->Get_X0(i));
        if(hit->code_ID().first==1)
          tdrift[jres]= -(hit->dtime());
        if(hit->code_ID().second==1)
          tdrift[jres]= hit->dtime();
        if(hit->code_ID().first==1)
          xhit[jres]= hit->Get_XPair(hit->dtime()).first;
        if(hit->code_ID().second==1)
          xhit[jres]= hit->Get_XPair(hit->dtime()).second;
        yhit[jres]= hit->y_wire_ID();
          }
        }
	    
        if(track->Get_IsGood_glo(i)){
          if(hit->code_ID().first==1)
        res_glo[jres]= hit->Get_XPair(hit->dtime()).first - (-1*0.00547*track->Get_T0_glo(i)+track->Get_Slope_glo(i)*(hit->y_wire_ID())+track->Get_X0_glo(i));
          if(hit->code_ID().second==1)
        res_glo[jres]= hit->Get_XPair(hit->dtime()).second - (1*0.00547*track->Get_T0_glo(i)+track->Get_Slope_glo(i)*(hit->y_wire_ID())+track->Get_X0_glo(i));
        }
      }
    } // close loop on hits
	

	if(i==0)  
	  for(int ii=0;ii<8;ii++) {
	    h_resCH1phi_glo[ii]->Fill(10000.*res_glo[ii]);
	    h_resCH1Phi_glo->Fill(10000.*res_glo[ii]);
	  }
	if(i==1)
	  for(int ii=0;ii<4;ii++){
	    h_resCH1the_glo[ii]->Fill(10000.*res_glo[ii]);
	    h_resCH1The_glo->Fill(10000.*res_glo[ii]);
	  }
	if(i==2)  
	  for(int ii=0;ii<8;ii++){
	    h_resCH2phi_glo[ii]->Fill(10000.*res_glo[ii]);
	    h_resCH2Phi_glo->Fill(10000.*res_glo[ii]);
	  }
	if(i==3)  
	  for(int ii=0;ii<4;ii++){
	    h_resCH2the_glo[ii]->Fill(10000.*res_glo[ii]);
	    h_resCH2The_glo->Fill(10000.*res_glo[ii]);
	  }
	if(i==4)  
	  for(int ii=0;ii<4;ii++){
	    h_resSL1[ii]->Fill(10000.*res[ii]);
	    h_resSLup->Fill(10000.*res[ii]);
	  }
	if(i==5)  
	  for(int ii=0;ii<4;ii++){
	    h_resSL2[ii]->Fill(10000.*res[ii]);
 	    h_resSLdown->Fill(10000.*res[ii]);
	  }

	if(calcola_spline && 
	   ( (track->Get_NPT(i)==8 && (i==0 || i==2) ) 
	     || 
	     (track->Get_NPT(i)==4 && (i!=0 && i!=2) ) ) ){
	  
	  //float m[9]={0.,0.08749,0.1763269,0.267949,0.363970,0.466308,0.577350,0.700207,0.839099};
	  float m[9]={0.,0.08749,0.1763269,0.267949,0.363970,0.466308,0.577350,0.700207,3.};
	  double vdrift = 0.00547;
	  double sum_x=0.;
	  double sum_y=0.;
	  double sum_x2=0.;
	  double sum_y2=0.;
	  double sum_xy=0.;
	  int npti     =0;            /* Numero di punti per vista phi */
	  double  delta[8];                 /* Linear fit, y = a + b x               */
	  double  a[8];		      /* a[layer]    */
	  double  b[8];     /* b[layer]           */
	  double resid[8];	
	  

	  for(int slo=0;slo<8;slo++){ // cut on slopes (8 intervals)
	    
	    if( TMath::Abs(track->Get_Slope(i))>=m[slo] && TMath::Abs(track->Get_Slope(i))<m[slo+1] )
	      if(i==0 || i==2){
		npti=7;
		for(int lay=0; lay<8; lay++){
		  for(int w=0 ; w<8 ; w++ ){              /* Run sui fili */
		    if ( w != lay ){                /* Esclusione del filo j */
		      sum_x += yhit[w];
		      sum_y += xhit[w];
		      sum_x2 += yhit[w]*yhit[w];
		      sum_y2 += xhit[w]*xhit[w];
		      sum_xy += xhit[w]*yhit[w];
		    }
		  }
		  delta[lay] = npti* sum_x2 - sum_x * sum_x;
		  a[lay] = (sum_x2 * sum_y - sum_x * sum_xy) / delta[lay];
		  b[lay] = (npti * sum_xy - sum_x * sum_y) / delta[lay];
		  resid[lay] = xhit[lay]/vdrift - (a[lay] + b[lay] * yhit[lay])/vdrift;
		  
 		  if ( tdrift[lay] < 0 )
 		    resid[lay] = -resid[lay];
		  
		  h_lincorr[i][slo]->Fill(tdrift[lay],resid[lay]);
		  h_lincorr8[slo]->Fill(tdrift[lay],resid[lay]);
		}
	      } // close if(i==0 || i==2)
	      else{
		npti=3;
		for(int lay=0; lay<4; lay++){
		  for(int w=0 ; w<4 ; w++ ){              /* Run sui fili */
		    if ( w != lay ){                /* Esclusione del filo j */
		      sum_x += yhit[w];
		      sum_y += xhit[w];
		      sum_x2 += yhit[w]*yhit[w];
		      sum_y2 += xhit[w]*xhit[w];
		      sum_xy += xhit[w]*yhit[w];
		    }
		  }

		  delta[lay] = npti* sum_x2 - sum_x * sum_x;
		  a[lay] = (sum_x2 * sum_y - sum_x * sum_xy) / delta[lay];
		  b[lay] = (npti * sum_xy - sum_x * sum_y) / delta[lay];
		  resid[lay] = xhit[lay]/vdrift - (a[lay] + b[lay] * yhit[lay])/vdrift;
		  
 		  if ( tdrift[lay] < 0 )
 		    resid[lay] = -resid[lay];
		  
		  h_lincorr[i][slo]->Fill(tdrift[lay],resid[lay]);
		}
		
	      } // close else()
	  } // close loop on slopes
	} // close if(calcola_spline)
	
	
      } // loop on segments
      

    } //close if(track->Get_IsGood(0)...)
  } // close if(Track_IsGood())
  

  if(hits->Get_NumberHITS()<1000){
    if(hits->Get_NumberHITS() != 0)
      h_Nhit->Fill(hits->Get_NumberHITS());
    int hit_xCH[4]={0,0,0,0}; 
    int hit_xSL[8]={0,0,0,0,0,0,0,0}; 
    int hit_xL[4]={0,0,0,0}; 
    
    for(int i=0;i<hits->Get_NumberHITS();i++){
      HIT *hit=hits->hit(i);
      if(DEBUG_STOREHIT)
	hit->print();
      if(hit->code_ID().first==1 || hit->code_ID().second==1){
	if(hit->CH_ID()==11 && hit->SL_ID()!=2) 
	  h_tempo[0]->Fill(hit->dtime());		  
	else if(hit->CH_ID()==11) 
	  h_tempo[1]->Fill(hit->dtime());		  
	if(hit->CH_ID()==10  && hit->SL_ID()!=2) 
	  h_tempo[2]->Fill(hit->dtime());		  
	else if(hit->CH_ID()==10) 
	  h_tempo[3]->Fill(hit->dtime());		  
	if(hit->CH_ID()==9) h_tempo[4]->Fill(hit->dtime());		  
	if(hit->CH_ID()==8) h_tempo[5]->Fill(hit->dtime());
      }		  
      for(int i=0;i<2;i++)
	if( hit->CH_ID()==(11-i) ){
	  hit_xCH[i]++;
	  if(hit->SL_ID()==1) hit_xSL[0+3*i]++;
	  if(hit->SL_ID()==2) hit_xSL[1+3*i]++;
	  if(hit->SL_ID()==3) hit_xSL[2+3*i]++;
	}
      for(int i=0;i<2;i++)
	if( hit->CH_ID()==(9-i) ){
	  hit_xCH[i+2]++;
	  if(hit->SL_ID()==1) hit_xSL[7]++;
	  if(hit->SL_ID()==3) hit_xSL[6]++;
	}
      if(hit->CH_ID()==11)
	if(hit->SL_ID()==1)
	  for(int i=0;i<4;i++) 
	    if(hit->L_ID()==(i+1)) 
	      hit_xL[i]++; 
    }
    for(int i=0;i<4;i++) if(hit_xCH[i] != 0) h_Nhit_xCH[i]->Fill(hit_xCH[i]);
    for(int i=0;i<8;i++) if(hit_xSL[i] != 0) h_Nhit_xSL[i]->Fill(hit_xSL[i]);
    for(int i=0;i<4;i++) if(hit_xL[i] != 0) h_Nhit_xL[i]->Fill(hit_xL[i]);

    int ii=0;
    int ii2=0;
    int ii3=0;
    if(numEvent==113 || numEvent==123 ||numEvent==125 ||numEvent==129 ||numEvent==138 || numEvent==142 ||numEvent==154 ||numEvent==166 ||numEvent==168 ||numEvent==170){
      for(int i=0;i<hits->Get_NumberHITS();i++){
	HIT *hit=hits->hit(i);
	if(hit->CH_ID()==11){
	  g1_hit[count_graph]->SetPoint(ii,hit->x_wire_ID() + (hit->dtime())*0.00547,hit->y_wire_ID());
	  g2_hit[count_graph]->SetPoint(ii,hit->x_wire_ID() - (hit->dtime())*0.00547,hit->y_wire_ID());
	  g1_hit_in[count_graph]->SetPoint(ii,hit->x_wire_ID() + (hit->dtime_in())*0.00547,hit->y_wire_ID());
	  g2_hit_in[count_graph]->SetPoint(ii,hit->x_wire_ID() - (hit->dtime_in())*0.00547,hit->y_wire_ID());
	  ii++;
	}
	if(hit->CH_ID()==10){
	  g3_hit[count_graph]->SetPoint(ii2,hit->x_wire_ID() + (hit->dtime())*0.00547,hit->y_wire_ID());
	  g4_hit[count_graph]->SetPoint(ii2,hit->x_wire_ID() - (hit->dtime())*0.00547,hit->y_wire_ID());
	  g3_hit_in[count_graph]->SetPoint(ii2,hit->x_wire_ID() + (hit->dtime_in())*0.00547,hit->y_wire_ID());
	  g4_hit_in[count_graph]->SetPoint(ii2,hit->x_wire_ID() - (hit->dtime_in())*0.00547,hit->y_wire_ID());
	  ii2++;
	}
	if(hit->CH_ID()==9){
	  g5_hit[count_graph]->SetPoint(ii3,hit->x_wire_ID() + (hit->dtime())*0.00547,hit->y_wire_ID());
	  g6_hit[count_graph]->SetPoint(ii3,hit->x_wire_ID() - (hit->dtime())*0.00547,hit->y_wire_ID());
	  ii3++;
	}
      }
      count_graph++;
    } // close if(numEvent==***)

  } // close if(hits->Get_NumberHITS()<1000)
  

  
  return;
  
}

void Save_HistosAndTree::initHistos(){
  
 
  // *** init histograms for debugging
  for(int i=0;i<ev_graph;i++){
    g1_hit[i]=new TGraph();
    g2_hit[i]=new TGraph();
    g1_hit[i]->SetName(Form("g1_hit_ev%d",i+1));
    g2_hit[i]->SetName(Form("g2_hit_ev%d",i+1));
    g1_hit[i]->SetMarkerStyle(5);
    g2_hit[i]->SetMarkerStyle(5);
    g1_hit[i]->SetMarkerColor(2);
    g2_hit[i]->SetMarkerColor(9);
    g1_hit_in[i]=new TGraph();
    g2_hit_in[i]=new TGraph();
    g1_hit_in[i]->SetName(Form("g1_hit_ev_in%d",i+1));
    g2_hit_in[i]->SetName(Form("g2_hit_ev_in%d",i+1));
    g1_hit_in[i]->SetMarkerStyle(5);
    g2_hit_in[i]->SetMarkerStyle(5);
    g1_hit_in[i]->SetMarkerColor(1);
    g2_hit_in[i]->SetMarkerColor(1);
    g3_hit[i]=new TGraph();
    g4_hit[i]=new TGraph();
    g3_hit[i]->SetName(Form("g3_hit_ev%d",i+1));
    g4_hit[i]->SetName(Form("g4_hit_ev%d",i+1));
    g3_hit[i]->SetMarkerStyle(5);
    g4_hit[i]->SetMarkerStyle(5);
    g3_hit[i]->SetMarkerColor(2);
    g4_hit[i]->SetMarkerColor(9);
    g3_hit_in[i]=new TGraph();
    g4_hit_in[i]=new TGraph();
    g3_hit_in[i]->SetName(Form("g3_hit_ev_in%d",i+1));
    g4_hit_in[i]->SetName(Form("g4_hit_ev_in%d",i+1));
    g3_hit_in[i]->SetMarkerStyle(5);
    g4_hit_in[i]->SetMarkerStyle(5);
    g3_hit_in[i]->SetMarkerColor(1);
    g4_hit_in[i]->SetMarkerColor(1);
    g5_hit[i]=new TGraph();
    g6_hit[i]=new TGraph();
    g5_hit[i]->SetName(Form("g5_hit_ev%d",i+1));
    g6_hit[i]->SetName(Form("g6_hit_ev%d",i+1));
    g5_hit[i]->SetMarkerStyle(5);
    g6_hit[i]->SetMarkerStyle(5);
    g5_hit[i]->SetMarkerColor(2);
    g6_hit[i]->SetMarkerColor(9);
  }
  for(int i=0;i<6;i++)
    h_tempo[i]=new TH1F(Form("h_tempo_%d",i),Form("Time (for seg%d) (ns)",i),1000,-100,1900);
  h_Nhit=new TH1F("h_Nhit","h_Nhit",1000,0,1000);
  for(int i=0;i<4;i++)
    h_Nhit_xCH[i]=new TH1F(Form("h_Nhit_xCH%d",i),Form("Nhit x CH %d",i),100,0,100);
  for(int i=0;i<8;i++)
    h_Nhit_xSL[i]=new TH1F(Form("h_Nhit_xSL%d",i),Form("Nhit x SL %d",i),100,0,100);
  for(int i=0;i<4;i++)
    h_Nhit_xL[i]=new TH1F(Form("h_Nhit_xL%d",i+1),Form("Nhit x L %d, CH 11, SL 1",i+1),10,0,10);


  h_Nhit_xSL13=new TH1F("h_Nhit_xSL13","Nhit x SL 13 dopo la cura",30,0,30);
  for(int i=0;i<6;i++){
    h_SLOPE[i]=new TH1F(Form("h_SLOPE_%d",i),Form("ATan(Slope) (rad) seg %d",i),1000,-2.,2.);
    h_X0[i]=new TH1F(Form("h_X0_%d",i),Form("X0 seg %d",i),1000,-200.,200.);
    h_T0[i]=new TH1F(Form("h_T0_%d",i),Form("T0 seg %d",i),1000,-100.,100.);
    h_CHI2[i]=new TH1F(Form("h_CHI1_%d",i),Form("CHI2 seg %d",i),1000,-200.,200.);
    h_NPT[i]=new TH1F(Form("h_NPT_%d",i),Form("NPT seg %d",i),20,0.,20.);
  }
  h_dphi=new TH1F("h_dphi","d(ATan(Slope)_{#Phi} (rad)",1000,-0.2,0.2);
  h_dthe=new TH1F("h_dthe","d(ATan(Slope)_{#Theta} (rad)",1000,-0.2,0.2);
  h_dT0=new TH1F("h_dT0","(T0#phi_{CH2} - T0#phi_{CH1}) first fit",2000,-100.,100.);
  h_T0ch1_T0ch2=new TH2F("h_T0ch1_T0ch2","T0#phi_{CH2} VS T0#phi_{CH1} ",100,-50.,50.,100,-50.,50.);
  h_dphi_glo=new TH1F("h_dphi_glo","d(ATan(Slope)_{#Phi} (rad) GLOBAL_FIT",1000,-0.2,0.2);
  h_dthe_glo=new TH1F("h_dthe_glo","d(ATan(Slope)_{#Theta} (rad) GLOBAL_FIT",1000,-0.2,0.2);
  h_dT0_glo=new TH1F("h_dT0_glo","(T0#phi_{CH2} - T0#phi_{CH1}) GLOBAL_FIT + SINGLE FIT",2000,-100.,100.);
  for(int i=0;i<4;i++){
    h_dT0_slope_glo[i]=new TH2F(Form("h_dT0_slope_glo_%d",i),Form("Slope_%d VS T0#phi_{CH2} - T0#phi_{CH1} GLOBAL_FIT",i),400,-20.,20.,200,-1,1);
    h_dT0_NPT_glo[i]=new TH2F(Form("h_dT0_NPT_glo_%d",i),Form("NPT_%d VS T0#phi_{CH2} - T0#phi_{CH1} GLOBAL_FIT",i),400,-20.,20.,10,0,10);
  }
  h_dT0_fin=new TH1F("h_dT0_fin","(T0#phi_{CH2} - T0#phi_{CH1}) second fit",2000,-100.,100.);
  h_dT0_fin2=new TH1F("h_dT0_fin2","T0#phi_{CH2} - T0#phi_{CH1} 1st+2nd fit",2000,-100.,100.);
  h_T0ch1_T0ch2_glo=new TH2F("h_T0ch1_T0ch2_glo","T0#phi_{CH2} VS T0#phi_{CH1} GLOBAL_FIT",100,-50.,50.,100,-50.,50.);
  
  for(int i=0;i<4;i++){
    h_SLOPE_glo[i]=new TH1F(Form("h_SLOPE_glo_%d",i),Form("Glo: ATan(Slope) (rad) seg %d",i),4000,-2.,2.);
    h_erSLOPE_glo[i]=new TH1F(Form("h_erSLOPE_glo_%d",i),Form("Glo: sigma ATan(Slope) (rad) seg %d",i),1000,0.,0.1);
    h_X0_glo[i]=new TH1F(Form("h_X0_glo_%d",i),Form("Glo: X0 (cm) seg %d",i),3000,-200,200);
    h_erX0_glo[i]=new TH1F(Form("h_erX0_glo_%d",i),Form("Glo: sigma X0 (cm) seg %d",i),100,0.,1.);
  }
  
  for(int i=0;i<2;i++){
    h_T0_glo[i]=new TH1F(Form("h_T0_glo_%d",i),Form("Glo: T0 (ns) seg %d",i),100,-50.,50.);
    h_erT0_glo[i]=new TH1F(Form("h_erT0_glo_%d",i),Form("Glo: sigma T0 (ns) seg %d",i),100,-50.,50.);
    h_MP_CH_glo[i]=new TH1F(Form("h_MP_CH%d_glo",i+1),Form("Mean Position, vista #Theta CH%d",i+1),1000,-10.,10.);
  }

  for(int i=0;i<8;i++){
    h_resCH1phi_glo[i]=new TH1F(Form("h_resCH1phi_glo_%d",i),Form("Residui CH1 Phi layer %d (um)",i),1000,-2000.,2000.);
    h_resCH2phi_glo[i]=new TH1F(Form("h_resCH2phi_glo_%d",i),Form("Residui CH2 Phi layer %d (um)",i),1000,-2000.,2000.);
  }
  for(int i=0;i<4;i++){
    h_resCH1the_glo[i]=new TH1F(Form("h_resCH1the_glo_%d",i),Form("Residui CH1 Theta layer %d (um)",i),1000,-2000.,2000.);
    h_resCH2the_glo[i]=new TH1F(Form("h_resCH2the_glo_%d",i),Form("Residui CH2 Theta layer %d (um)",i),1000,-2000.,2000.);
    h_resSL1[i]=new TH1F(Form("h_resSL1_%d",i),Form("Residui SL1 layer %d (um)",i),1000,-2000.,2000.);
    h_resSL2[i]=new TH1F(Form("h_resSL2_%d",i),Form("Residui SL2 layer %d (um)",i),1000,-2000.,2000.);
  }
  
  h_resCH1Phi_glo=new TH1F("h_resCH1Phi_glo","Residui CH1 #Phi (um)",1000,-2000.,2000.);
  h_resCH1The_glo=new TH1F("h_resCH1The_glo","Residui CH1 #Theta (um)",1000,-2000.,2000.);
  h_resCH2Phi_glo=new TH1F("h_resCH2Phi_glo","Residui CH2 #Phi (um)",1000,-2000.,2000.);
  h_resCH2The_glo=new TH1F("h_resCH2The_glo","Residui CH2 #Theta (um)",1000,-2000.,2000.);
  h_resSLup=new TH1F("h_resSLup","Residui SL1 (um)",1000,-2000.,2000.);
  h_resSLdown=new TH1F("h_resSLdown","Residui SL2 (um)",1000,-2000.,2000.);
  
  if(calcola_spline)
    for(int i=0; i<m_nmaxseg; i++)
      for(int slo=0; slo<8; slo++)  
	h_lincorr[i][slo]=new TH2F(Form("h_lincorr_seg%d_slope%d",i,slo),Form("dist track-wire (ns) VS res (ns), seg%d slope%d",i,slo),90,0,450,100,-50,50);

  for(int slo=0; slo<8; slo++)  
    h_lincorr8[slo]=new TH2F(Form("h_lincorr_slope%d",slo),Form("dist track-wire (ns) VS res (ns), vista phi slope%d",slo),90,0,450,100,-50,50);
  return;
  
}

void Save_HistosAndTree::writeHistos(){
  for(int i=0;i<ev_graph;i++){
    g1_hit[i]->Write();  
    g2_hit[i]->Write();  
    g3_hit[i]->Write();  
    g4_hit[i]->Write();  
    g5_hit[i]->Write();  
    g6_hit[i]->Write();  
    g1_hit_in[i]->Write();  
    g2_hit_in[i]->Write();  
    g3_hit_in[i]->Write();  
    g4_hit_in[i]->Write();  
  }
  
  for(int i=0;i<6;i++)
    h_tempo[i]->Write();  
  h_Nhit->Write();  
  for(int i=0;i<4;i++)
    h_Nhit_xCH[i]->Write();
  for(int i=0;i<8;i++)
    h_Nhit_xSL[i]->Write();
  for(int i=0;i<4;i++)
    h_Nhit_xL[i]->Write();
  h_Nhit_xSL13->Write();
  for(int i=0;i<6;i++){
    h_SLOPE[i]->Write(); 
    h_X0[i]->Write(); 
    h_T0[i]->Write(); 
    h_NPT[i]->Write(); 
    h_CHI2[i]->Write(); 
  } 
  h_dphi->Write(); 
  h_dthe->Write(); 
  h_dT0->Write(); 
  h_T0ch1_T0ch2->Write(); 
  h_dphi_glo->Write(); 
  h_dthe_glo->Write(); 
  h_dT0_glo->Write(); 
  for(int i=0;i<4;i++){
    h_dT0_slope_glo[i]->Write(); 
    h_dT0_NPT_glo[i]->Write(); 
  }
  h_dT0_fin->Write(); 
  h_dT0_fin2->Write(); 
  h_T0ch1_T0ch2_glo->Write(); 
  
  for(int i=0;i<4;i++){
    h_SLOPE_glo[i]->Write();
    h_erSLOPE_glo[i]->Write();
    h_X0_glo[i]->Write();
    h_erX0_glo[i]->Write();
  }
  for(int i=0;i<2;i++){
    h_T0_glo[i]->Write();
    h_erT0_glo[i]->Write();
    h_MP_CH_glo[i]->Write();
  }
  
  for(int i=0;i<8;i++){
    h_resCH1phi_glo[i]->Write();
    h_resCH2phi_glo[i]->Write();
  }
  for(int i=0;i<4;i++){
    h_resCH1the_glo[i]->Write();
    h_resCH2the_glo[i]->Write();
    h_resSL1[i]->Write();
    h_resSL2[i]->Write();
  }

  h_resCH1Phi_glo->Write();
  h_resCH1The_glo->Write();
  h_resCH2Phi_glo->Write();
  h_resCH2The_glo->Write();
  h_resSLup->Write();
  h_resSLdown->Write();

  if(calcola_spline)
    for(int i=0; i<m_nmaxseg; i++)
      for(int slo=0; slo<8; slo++)  
	h_lincorr[i][slo]->Write();  
  if(calcola_spline)
    for(int slo=0; slo<8; slo++)  
      h_lincorr8[slo]->Write();  
  
  return;
}

void Save_HistosAndTree::resetHistos()
{
  // *** Reset histograms
  for(int i=0;i<6;i++)
    h_tempo[i]->Reset();  
  h_Nhit->Reset();  
  for(int i=0;i<4;i++)
    h_Nhit_xCH[i]->Reset();
  for(int i=0;i<8;i++)
    h_Nhit_xSL[i]->Reset();
  for(int i=0;i<4;i++)
    h_Nhit_xL[i]->Reset();
  h_Nhit_xSL13->Reset();
  for(int i=0;i<6;i++){
    h_SLOPE[i]->Reset(); 
    h_X0[i]->Reset(); 
    h_T0[i]->Reset(); 
    h_CHI2[i]->Reset(); 
    h_NPT[i]->Reset(); 
  } 
  h_dphi->Reset(); 
  h_dthe->Reset(); 
  h_dT0->Reset(); 
  h_T0ch1_T0ch2->Reset(); 
  h_dphi_glo->Reset(); 
  h_dthe_glo->Reset(); 
  h_dT0_glo->Reset(); 
  for(int i=0;i<4;i++){
  h_dT0_slope_glo[i]->Reset(); 
  h_dT0_NPT_glo[i]->Reset(); 
  }
  h_dT0_fin->Reset(); 
  h_dT0_fin2->Reset(); 
  h_T0ch1_T0ch2_glo->Reset(); 
  
  for(int i=0;i<4;i++){
    h_SLOPE_glo[i]->Reset();
    h_erSLOPE_glo[i]->Reset();
    h_X0_glo[i]->Reset();
    h_erX0_glo[i]->Reset();
  }
  for(int i=0;i<2;i++){
    h_T0_glo[i]->Reset();
    h_erT0_glo[i]->Reset();
    h_MP_CH_glo[i]->Reset();
  }

  for(int i=0;i<8;i++){
    h_resCH1phi_glo[i]->Reset();
    h_resCH2phi_glo[i]->Reset();
  }
  for(int i=0;i<4;i++){
    h_resCH1the_glo[i]->Reset();
    h_resCH2the_glo[i]->Reset();
    h_resSL1[i]->Reset();
    h_resSL2[i]->Reset();
  }
  
  h_resCH1Phi_glo->Reset();
  h_resCH1The_glo->Reset();
  h_resCH2Phi_glo->Reset();
  h_resCH2The_glo->Reset();
  h_resSLup->Reset();
  h_resSLdown->Reset();

  if(calcola_spline)
    for(int i=0; i<m_nmaxseg; i++)
      for(int slo=0; slo<8; slo++)  
	h_lincorr[i][slo]->Reset();  
  if(calcola_spline)
    for(int slo=0; slo<8; slo++)  
      h_lincorr8[slo]->Reset();  
  
  return;
}

void Save_HistosAndTree::deleteHistos()
{
  // *** Delete histograms 
  for(int i=0;i<6;i++)
    h_tempo[i]->Delete();  
  h_Nhit->Delete();  
  for(int i=0;i<4;i++)
    h_Nhit_xCH[i]->Delete();
  for(int i=0;i<8;i++)
    h_Nhit_xSL[i]->Delete();
  for(int i=0;i<4;i++)
    h_Nhit_xL[i]->Delete();
  h_Nhit_xSL13->Delete();
  for(int i=0;i<6;i++){
    h_SLOPE[i]->Delete();  
    h_X0[i]->Delete();  
    h_T0[i]->Delete();  
    h_CHI2[i]->Delete();  
    h_NPT[i]->Delete();  
  }
  h_dphi->Delete(); 
  h_dthe->Delete(); 
  h_dT0->Delete(); 
  h_T0ch1_T0ch2->Delete(); 
  h_dphi_glo->Delete(); 
  h_dthe_glo->Delete(); 
  h_dT0_glo->Delete(); 
  for(int i=0;i<4;i++){
  h_dT0_slope_glo[i]->Delete(); 
  h_dT0_NPT_glo[i]->Delete(); 
  }
  h_dT0_fin->Delete(); 
  h_dT0_fin2->Delete(); 
  h_T0ch1_T0ch2_glo->Delete(); 
  
  for(int i=0;i<4;i++){
    h_SLOPE_glo[i]->Delete();
    h_erSLOPE_glo[i]->Delete();
    h_X0_glo[i]->Delete();
    h_erX0_glo[i]->Delete();
  }
  for(int i=0;i<2;i++){
    h_T0_glo[i]->Delete();
    h_erT0_glo[i]->Delete();
    h_MP_CH_glo[i]->Delete();
  }
  
  for(int i=0;i<8;i++){
    h_resCH1phi_glo[i]->Delete();
    h_resCH2phi_glo[i]->Delete();
  }
  for(int i=0;i<4;i++){
    h_resCH1the_glo[i]->Delete();
    h_resCH2the_glo[i]->Delete();
    h_resSL1[i]->Delete();
    h_resSL2[i]->Delete();
  }

  h_resCH1Phi_glo->Delete();
  h_resCH1The_glo->Delete();
  h_resCH2Phi_glo->Delete();
  h_resCH2The_glo->Delete();
  h_resSLup->Delete();
  h_resSLdown->Delete();

  if(calcola_spline)
    for(int i=0; i<m_nmaxseg; i++)
      for(int slo=0; slo<8; slo++)  
	h_lincorr[i][slo]->Delete();  
  if(calcola_spline)
    for(int slo=0; slo<8; slo++)  
      h_lincorr8[slo]->Delete();  
  
  return;
}

void Save_HistosAndTree::initTree(){
  
  if ( DEBUG_TREE ) cout << "Init output tree" << endl;
  
  onevent = 0;

  int nmaxhit = 100;
  ohtrig = new int[nmaxhit];     // ttrig+t0
  ohlay = new int[nmaxhit];
  ohwire = new int[nmaxhit];
  ohtime_in = new int[nmaxhit];  // time dtime_in()
  ohtime = new int[nmaxhit];     // time dtime()
  ohtime_tube = new int[108];  // time new tube to be tested
  oNhit_tube = new int[18];  // Nhit new tube to be tested
  
  onseg = m_nmaxseg;
  onseg_glo = m_nmaxseg_glo;
  osegS = new float[m_nmaxseg];
  osegX = new float[m_nmaxseg];
  osegK = new float[m_nmaxseg];
  osegN = new int[m_nmaxseg];
  osegT0 = new float[m_nmaxseg];
  osl1r = new float[m_nmaxseg];
  osl2r = new float[m_nmaxseg];
  osl3r = new float[m_nmaxseg];
  osl4r = new float[m_nmaxseg];
  osl5r = new float[m_nmaxseg];
  osl6r = new float[m_nmaxseg];
  osl7r = new float[m_nmaxseg];
  osl8r = new float[m_nmaxseg];

  osxh1 = new float[m_nmaxseg];
  osxh2 = new float[m_nmaxseg];
  osxh3 = new float[m_nmaxseg];
  osxh4 = new float[m_nmaxseg];
  osxh5 = new float[m_nmaxseg];
  osxh6 = new float[m_nmaxseg];
  osxh7 = new float[m_nmaxseg];
  osxh8 = new float[m_nmaxseg];

  osegS_glo = new float[m_nmaxseg];
  osegerS_glo = new float[m_nmaxseg];
  osegX_glo = new float[m_nmaxseg];
  osegerX_glo = new float[m_nmaxseg];
  osegK_glo = new float[m_nmaxseg];
  osegN_glo = new int[m_nmaxseg];
  osegT0_glo = new float[m_nmaxseg];
  osl1r_glo = new float[m_nmaxseg];
  osl2r_glo = new float[m_nmaxseg];
  osl3r_glo = new float[m_nmaxseg];
  osl4r_glo = new float[m_nmaxseg];
  osl5r_glo = new float[m_nmaxseg];
  osl6r_glo = new float[m_nmaxseg];
  osl7r_glo = new float[m_nmaxseg];
  osl8r_glo = new float[m_nmaxseg];


  return; 
}


void Save_HistosAndTree::bookTree(TTree* tree,bool n2chambers)
{
  if ( DEBUG_TREE ) cout << "Book output tree" << endl;

  tree->Branch( "EVENT", &onevent,  "EVENT/I" );
  tree->Branch( "DTBX_nhit", &onhit, "Nhits/I" );
  tree->Branch( "TRTIME",     ohtrig,"trigT[Nhits]/I" );
  tree->Branch( "DTBX_lay",   ohlay, "hlay[Nhits]/I" );
  tree->Branch( "DTBX_tube",  ohwire,"htube[Nhits]/I" );
  tree->Branch( "DTBX_time_in",ohtime_in,"htime_in[Nhits]/I" );
  tree->Branch( "DTBX_time",  ohtime,"htime[Nhits]/I" );
  tree->Branch( "DTBX_time_tube",ohtime_tube,"rtime[108]/I" );
  tree->Branch( "DTBX_Nhit_tube",oNhit_tube,"Nhit[18]/I" );
  tree->Branch( "SEG_ns", &onseg,  "ntes/I" );
  tree->Branch( "SEG_sx",  osegX,  "X[ntes]/F" );
  tree->Branch( "SEG_ss",  osegS,  "SLOPE[ntes]/F" );
  tree->Branch( "SEG_sk",  osegK,  "CHI2[ntes]/F" );
  tree->Branch( "SEG_sn",  osegN,  "NPT[ntes]/I" );
  tree->Branch( "SEG_t0",  osegT0, "T0[ntes]/F" );
  tree->Branch( "SEG_1r",  osl1r,   "l1r[ntes]/F" );
  tree->Branch( "SEG_2r",  osl2r,   "l2r[ntes]/F" );
  tree->Branch( "SEG_3r",  osl3r,   "l3r[ntes]/F" );
  tree->Branch( "SEG_4r",  osl4r,   "l4r[ntes]/F" );
  tree->Branch( "SEG_5r",  osl5r,   "l5r[ntes]/F" );
  tree->Branch( "SEG_6r",  osl6r,   "l6r[ntes]/F" );
  tree->Branch( "SEG_7r",  osl7r,   "l7r[ntes]/F" );
  tree->Branch( "SEG_8r",  osl8r,   "l8r[ntes]/F" );

  tree->Branch( "SEG_xh1",  osxh1,   "xh1[ntes]/F" );
  tree->Branch( "SEG_xh2",  osxh2,   "xh2[ntes]/F" );
  tree->Branch( "SEG_xh3",  osxh3,   "xh3[ntes]/F" );
  tree->Branch( "SEG_xh4",  osxh4,   "xh4[ntes]/F" );
  tree->Branch( "SEG_xh5",  osxh5,   "xh5[ntes]/F" );
  tree->Branch( "SEG_xh6",  osxh6,   "xh6[ntes]/F" );
  tree->Branch( "SEG_xh7",  osxh7,   "xh7[ntes]/F" );
  tree->Branch( "SEG_xh8",  osxh8,   "xh8[ntes]/F" );

  if(n2chambers){
      tree->Branch( "SEG_ns_glo", &onseg_glo,  "ntes_glo/I" );
      tree->Branch( "SEG_sx_glo",  osegX_glo,  "X_glo[ntes_glo]/F" );
      tree->Branch( "SEG_ersx_glo",osegerX_glo,"Xer_glo[ntes_glo]/F" );
      tree->Branch( "SEG_ss_glo",  osegS_glo,  "SLOPE_glo[ntes_glo]/F" );
      tree->Branch( "SEG_erss_glo",osegerS_glo,"SLOPEer_glo[ntes_glo]/F" );
      tree->Branch( "SEG_sk_glo",  osegK_glo,  "CHI2_glo[ntes_glo]/F" );
      tree->Branch( "SEG_sn_glo",  osegN_glo,  "NPT_glo[ntes_glo]/I" );
      tree->Branch( "SEG_t0_glo",  osegT0_glo, "T0_glo[ntes_glo]/F" );
      tree->Branch( "SEG_1r_glo",  osl1r_glo,   "l1r_glo[ntes_glo]/F" );
      tree->Branch( "SEG_2r_glo",  osl2r_glo,   "l2r_glo[ntes_glo]/F" );
      tree->Branch( "SEG_3r_glo",  osl3r_glo,   "l3r_glo[ntes_glo]/F" );
      tree->Branch( "SEG_4r_glo",  osl4r_glo,   "l4r_glo[ntes_glo]/F" );
      tree->Branch( "SEG_5r_glo",  osl5r_glo,   "l5r_glo[ntes_glo]/F" );
      tree->Branch( "SEG_6r_glo",  osl6r_glo,   "l6r_glo[ntes_glo]/F" );
      tree->Branch( "SEG_7r_glo",  osl7r_glo,   "l7r_glo[ntes_glo]/F" );
      tree->Branch( "SEG_8r_glo",  osl8r_glo,   "l8r_glo[ntes_glo]/F" );
}
  if ( DEBUG_TREE ) cout << "End booking" << endl;
  return;
}

void Save_HistosAndTree::cleanTree() {
  if ( DEBUG_TREE ) cout << "Clean output tree" << endl;
  
  int i;
  for ( i = 0; i < onseg; i++ ){
    osegS[i] = -999.;
    osegX[i] = -999.;
    osegK[i] = -999.;
    osegN[i] = -999;
    osegT0[i] = -999.;
    osl1r[i] = -999.;
    osl2r[i] = -999.;
    osl3r[i] = -999.;
    osl4r[i] = -999.;
    osl5r[i] = -999.;
    osl6r[i] = -999.;
    osl7r[i] = -999.;
    osl8r[i] = -999.;

    osxh1[i] = -999.;
    osxh2[i] = -999.;
    osxh3[i] = -999.;
    osxh4[i] = -999.;
    osxh5[i] = -999.;
    osxh6[i] = -999.;
    osxh7[i] = -999.;
    osxh8[i] = -999.;
 }

  for ( i = 0; i < onseg; i++ ){
    osegS_glo[i] = -999.;
    osegerS_glo[i] = -999.;
    osegX_glo[i] = -999.;
    osegerX_glo[i] = -999.;
    osegK_glo[i] = -999.;
    osegN_glo[i] = -999;
    osegT0_glo[i] = -999.;
    osl1r_glo[i] = -999.;
    osl2r_glo[i] = -999.;
    osl3r_glo[i] = -999.;
    osl4r_glo[i] = -999.;
    osl5r_glo[i] = -999.;
    osl6r_glo[i] = -999.;
    osl7r_glo[i] = -999.;
    osl8r_glo[i] = -999.;
  }
  
  for ( i = 0; i < 100; i++ ){
    ohtrig[i] = -999;
    ohlay[i] = -999;
    ohwire[i] = -999;
    ohtime_in[i] = -999;
    ohtime[i] = -999;
  }
  for ( i = 0; i < 108; i++ )
    ohtime_tube[i] = -999;
  for ( i = 0; i < 18; i++ )
    oNhit_tube[i] = 0;
  
  onseg = 0;
  onseg_glo = 0;
  onevent = 0;
  onhit = 0;
  
  if ( DEBUG_TREE ) cout << "End clean" << endl;
  return;
}

void Save_HistosAndTree::fillVar(int npc, int onseg, double m, double a, double t0, double chi2, double * res, int p, double  * xh)
{
  
  if(DEBUG_STORETRACK)
    printf("Filling Variables...\n");

  osegN[onseg] = npc;
  osegS[onseg] = m;
  osegX[onseg] = a;
  osegT0[onseg] = t0;
  osegK[onseg] = chi2;

  osl1r[onseg] = res[0];
  osl2r[onseg] = res[1];
  osl3r[onseg] = res[2];
  osl4r[onseg] = res[3];

  osxh1[onseg] = xh[0];
  osxh2[onseg] = xh[1];
  osxh3[onseg] = xh[2];
  osxh4[onseg] = xh[3];


  if(p==0||p==2)
    {
      osl5r[onseg] = res[4];
      osl6r[onseg] = res[5];
      osl7r[onseg] = res[6];
      osl8r[onseg] = res[7];

      osxh5[onseg] = xh[4];
      osxh6[onseg] = xh[5];
      osxh7[onseg] = xh[6];
      osxh8[onseg] = xh[7];
    }

  if(DEBUG_STORETRACK) 
    printf("Variables filled...\n");

  return;

}

void Save_HistosAndTree::fillVar_glo(int npc, int onseg, double m, double erm, double a, double era, double t0, double chi2, double * res, int p)
{
  if(DEBUG_STORETRACK) 
    printf("Filling Variables global...\n");

  osegN_glo[onseg] = npc;
  osegS_glo[onseg] = m;
  osegerS_glo[onseg] = erm;
  osegX_glo[onseg] = a;
  osegerX_glo[onseg] = era;
  osegT0_glo[onseg] = t0;
  osegK_glo[onseg] = chi2;
  osl1r_glo[onseg] = res[0];
  osl2r_glo[onseg] = res[1];
  osl3r_glo[onseg] = res[2];
  osl4r_glo[onseg] = res[3];
  
  if(p==0||p==2)
    {
      osl5r_glo[onseg] = res[4];
      osl6r_glo[onseg] = res[5];
      osl7r_glo[onseg] = res[6];
      osl8r_glo[onseg] = res[7];
    }

  if(DEBUG_STORETRACK) 
    printf("Variables global filled...\n");

  return;


}


void Save_HistosAndTree::init_Statistics(){
  // initialize variables for statistics_file
  ev_analyzed=0;
  ev_CH1_fitOk=0; ev_CH2_fitOk=0; ev_SL1_fitOk=0; ev_SL2_fitOk=0;
  ev_CH1_fitgloOk=0; ev_CH2_fitgloOk=0;
  ev_CH1_fitgloOk_NPT8=0; ev_CH1_fitgloOk_NPT7=0; ev_CH1_fitgloOk_NPT6=0;
  ev_CH2_fitgloOk_NPT8=0; ev_CH2_fitgloOk_NPT7=0; ev_CH2_fitgloOk_NPT6=0;
  ev_CH1_fitgloOk_MP=0; ev_CH2_fitgloOk_MP=0;
  ev_CH1_fitgloOk_MP_NPT8=0; ev_CH1_fitgloOk_MP_NPT7=0; ev_CH1_fitgloOk_MP_NPT6=0;
  ev_CH2_fitgloOk_MP_NPT8=0; ev_CH2_fitgloOk_MP_NPT7=0; ev_CH2_fitgloOk_MP_NPT6=0;
  res_CH1_phi_mean=-999.; res_CH1_the_mean=-999.; res_CH1_phi_sigma=-999.; res_CH1_the_sigma=-999.; 
  res_CH2_phi_mean=-999.; res_CH2_the_mean=-999.; res_CH2_phi_sigma=-999.; res_CH2_the_sigma=-999.; 
  res_SL1_mean=-999.; res_SL2_mean=-999.; res_SL1_sigma=-999.; res_SL2_sigma=-999.; 

  return;
}

void Save_HistosAndTree::compute_Statistics(Track *track, HITCollection *hits,int numEvent){
  
  ev_analyzed=numEvent;
  if(track->Track_IsGood()){
    if(track->Get_IsGood(0) && track->Get_IsGood(1)){
      ev_CH1_fitOk++;

      if(track->Get_IsGood(2) && track->Get_IsGood(3))
	ev_CH2_fitOk++;
      if(track->Get_IsGood(4))
	ev_SL1_fitOk++;
      if(track->Get_IsGood(5))
	ev_SL2_fitOk++;
      
      if(track->Get_IsGood_glo(0) && track->Get_IsGood_glo(1)){
	ev_CH1_fitgloOk++;
	if(track->Get_NPT_glo(0)==8)
	  ev_CH1_fitgloOk_NPT8++;
	if(track->Get_NPT_glo(0)==7)
	  ev_CH1_fitgloOk_NPT7++;
	if(track->Get_NPT_glo(0)==6)
	  ev_CH1_fitgloOk_NPT6++;
	
	if(track->Get_MP_CH1_IsOk_glo()){
	  ev_CH1_fitgloOk_MP++;
	  if(track->Get_NPT_glo(0)==8)
	    ev_CH1_fitgloOk_MP_NPT8++;
	  if(track->Get_NPT_glo(0)==7)
	    ev_CH1_fitgloOk_MP_NPT7++;
	  if(track->Get_NPT_glo(0)==6)
	    ev_CH1_fitgloOk_MP_NPT6++;
	}
      }
      if(track->Get_IsGood_glo(2) && track->Get_IsGood_glo(3)){
	ev_CH2_fitgloOk++;
	if(track->Get_NPT_glo(2)==8)
	  ev_CH2_fitgloOk_NPT8++;
	if(track->Get_NPT_glo(2)==7)
	  ev_CH2_fitgloOk_NPT7++;
	if(track->Get_NPT_glo(2)==6)
	  ev_CH2_fitgloOk_NPT6++;
	
 	if(track->Get_MP_CH2_IsOk_glo()){
	  ev_CH2_fitgloOk_MP++;
	  if(track->Get_NPT_glo(2)==8)
	    ev_CH2_fitgloOk_MP_NPT8++;
	  if(track->Get_NPT_glo(2)==7)
	    ev_CH2_fitgloOk_MP_NPT7++;
	  if(track->Get_NPT_glo(2)==6)
	    ev_CH2_fitgloOk_MP_NPT6++;
	}
      }
    }
  }
  
  return;
}


void Save_HistosAndTree::write_Statistics(FILE *fo_txt){
  
  fprintf (fo_txt, "Events analyzed %d\n",ev_analyzed);
  fprintf (fo_txt, "---> FIT OK  (FIT CH1 has to be Ok) \n");
  fprintf (fo_txt, "     CH1 %d,    CH2 %d,    Sl1 %d,    Sl2 %d\n",ev_CH1_fitOk, ev_CH2_fitOk, ev_SL1_fitOk, ev_SL2_fitOk);
  fprintf (fo_txt, "---> GLOBALFIT OK  (FIT CH1 has to be Ok) \n");
  fprintf (fo_txt, "     CH1 %d:  NPT 8: %d,  NPT 7 %d,  NPT 6 %d\n",ev_CH1_fitgloOk,ev_CH1_fitgloOk_NPT8,ev_CH1_fitgloOk_NPT7,ev_CH1_fitgloOk_NPT6);
  fprintf (fo_txt, "     CH2 %d:  NPT 8: %d,  NPT 7 %d,  NPT 6 %d\n",ev_CH2_fitgloOk,ev_CH2_fitgloOk_NPT8,ev_CH2_fitgloOk_NPT7,ev_CH2_fitgloOk_NPT6);
  fprintf (fo_txt, "---> GLOBALFIT OK and Mean Position OK  (FIT CH1 has to be Ok) \n");
  fprintf (fo_txt, "     CH1 %d:  NPT 8: %d,  NPT 7 %d,  NPT 6 %d\n",ev_CH1_fitgloOk_MP,ev_CH1_fitgloOk_MP_NPT8,ev_CH1_fitgloOk_MP_NPT7,ev_CH1_fitgloOk_MP_NPT6);
  fprintf (fo_txt, "     CH2 %d:  NPT 8: %d,  NPT 7 %d,  NPT 6 %d\n",ev_CH2_fitgloOk_MP,ev_CH2_fitgloOk_MP_NPT8,ev_CH2_fitgloOk_MP_NPT7,ev_CH2_fitgloOk_MP_NPT6);
  fprintf (fo_txt, "---> MEAN and SIGMA of RESIDUALS (um) (FITGLO OK for CH1 and CH2, FIT OK for SL1 and SL2)\n");
  fprintf (fo_txt, "     CH1:  PHI:  X0 = %.1f  RMS = %.1f,  THETA:  X0 = %.1f  RMS = %.1f\n",
       h_resCH1Phi_glo->GetMean(),h_resCH1Phi_glo->GetRMS(),h_resCH1The_glo->GetMean(),h_resCH1The_glo->GetRMS());
  fprintf (fo_txt, "     CH2:  PHI:  X0 = %.1f  RMS = %.1f,  THETA:  X0 = %.1f  RMS = %.1f\n",
       h_resCH2Phi_glo->GetMean(),h_resCH2Phi_glo->GetRMS(),h_resCH2The_glo->GetMean(),h_resCH2The_glo->GetRMS());
  fprintf (fo_txt, "     SL1:        X0 = %.1f  RMS = %.1f\n",
       h_resSLup->GetMean(),h_resSLup->GetRMS());
  fprintf (fo_txt, "     SL2:        X0 = %.1f  RMS = %.1f\n",
       h_resSLdown->GetMean(),h_resSLdown->GetRMS());
  
  return;
  
}
