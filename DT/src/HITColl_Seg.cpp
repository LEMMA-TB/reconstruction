#include "HITColl_Seg.h"

static int CH_HVmod = 10;
static int SL_HVmod = 3;
static int LAY_HVmod = 4;

HITColl_Seg::~HITColl_Seg()
{
  int N=Get_NumberHITS();
  if(DEBUG_HITCOLL_S)
    printf("...Into ~HITColl_Seg() deleting N.%d hits...\n",N);
  for(int h=(N-1);h>=0;h--){
    eraseHIT(h);
  }
  
  return;
}

void HITColl_Seg::selectHIT(HITCollection * hits, int CH, bool phi) 
{
  if(DEBUG_HITCOLL_S)
    cout << "\n\n *** HITColl_Seg::selectHIT" << endl;

  int j=0;  
  for(int i=0;i<hits->Get_NumberHITS();i++){
    HIT *hit=hits->hit(i);
    
    int SL1=0;
    int SL2=0;
    if(CH==10||CH==11)
        if(phi==false){
            SL1=2;
            SL2=0;
        } else if(phi==true){
            SL1=1;
            SL2=3;
        }
        if(CH==9){
            SL1=3;
            SL2=0;
        }
        if(CH==8){
            SL1=1;
            SL2=0;
        }
    
//// 20170405 if hit is in the layer with HV modified, to be excluded from fit
//    if( hit->CH_ID()==CH_HVmod && hit->SL_ID()==SL_HVmod && hit->L_ID()==LAY_HVmod)
//       continue;

    if( hit->CH_ID()==CH && (hit->SL_ID()==SL1 || hit->SL_ID()==SL2) )
      if( hit->dtime_in()>-10. && hit->dtime_in()<410.5 )
	{
	  addHIT(hit);
	  if(DEBUG_HITCOLL_S)
	    printHIT(j);
	  j++;
	}
      else {
	hit->Change_LCode(false);
	hit->Change_RCode(false);
      }
    
  }
  return;
}

void HITColl_Seg::selectHIT(HITCollection * hits, int CH, int SL) 
{
  if(DEBUG_HITCOLL_S)
    cout << "\n\n *** HITColl_Seg::selectHIT" << endl;
  
  int j=0;
  for(int i=0;i<hits->Get_NumberHITS();i++){
    HIT *hit=hits->hit(i);
    
    if(hit->CH_ID()==CH && hit->SL_ID()==SL 
       &&
       hit->dtime_in()>-10. && hit->dtime_in()<410.5)
      {
	addHIT(hit);
	if(DEBUG_HITCOLL_S)
	  printHIT(j);
	j++;
      }
    
  }
  return;
}

void HITColl_Seg::addHIT(HIT * hit)
{
  if(DEBUG_HITCOLL_S)
    cout << "HITColl_Seg::addHIT" << endl;
  
  _hits_seg.push_back(hit);
  
  return;
}

void HITColl_Seg::eraseHIT(int i)
{
  if(DEBUG_HITCOLL_S)
    cout << "HITColl_Seg::eraseHIT n." << i << endl;

  
  _hits_seg.erase(_hits_seg.begin()+i);
  
  return;
}

HIT * HITColl_Seg::hit(int i)
{	
  HIT * hit_ptr = _hits_seg[i];
  
  return hit_ptr;
}

int HITColl_Seg::Get_NumberHITS()
{
  int _H = _hits_seg.size();
  
  return _H;
}

void HITColl_Seg::printHIT(int i)
{
  HIT * h = hit(i); 
  cout << "HIT " << i << " information" << endl;
  h->print();
  return;
}
