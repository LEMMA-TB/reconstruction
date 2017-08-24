#include "HITColl_Layer.h"


HITColl_Layer::~HITColl_Layer()
{
  for(int h=0; h<Get_NumberHITS(); h++)
    eraseHIT(h);
  if(DEBUG_HITCOLL_L) Printf("I'm deleting hit into HITColl_Layer\n");
  return;
}

void HITColl_Layer::selectHIT(HITColl_Seg * hit_seg,int L) 
{
  if(DEBUG_HITCOLL_L)
    cout << "\n\n *** HITColl_Layer::selectHIT" << endl;
  
  HIT *hit=hit_seg->hit(0);  
  int CH=hit->CH_ID();
  
  int j=0;
  for(int i=0;i<hit_seg->Get_NumberHITS();i++){
    HIT *hit=hit_seg->hit(i);
    int SL=hit->SL_ID();
    if( ((CH==10||CH==11) && (SL==1 || SL==2))
	|| CH==8 || CH==9 )
      if(hit->L_ID()==L){
	addHIT(hit);
	if(!(hit->IsSolved()))
	  addHIT(hit);
	if(DEBUG_HITCOLL_L){
	  printHIT(j);
	  j++;
	  if(!(hit->IsSolved())){
	  printHIT(j);
	  j++;
	  }
	}
      }
    
    if( (CH==10||CH==11) && SL==3 )
      if(hit->L_ID()==(L-4)){
	addHIT(hit);
	if(!(hit->IsSolved()))
	  addHIT(hit);
	if(DEBUG_HITCOLL_L){
	  printHIT(j);
	  j++;
	  if(!(hit->IsSolved())){
	    printHIT(j);
	    j++;
	  }
	}
      }
  }
  return;
}

void HITColl_Layer::addHIT(HIT * hit)
{
  if(DEBUG_HITCOLL_L)
    cout << "HITColl_Layer::addHIT" << endl;
  
  _hits_layer.push_back(hit);
  
  return;
}

void HITColl_Layer::eraseHIT(int i)
{
  
  _hits_layer.erase(_hits_layer.begin()+i);
  
  return;
}

HIT * HITColl_Layer::hit(int i)
{	
  HIT * hit_ptr = _hits_layer[i];
  
  return hit_ptr;
}

int HITColl_Layer::Get_NumberHITS()
{
  int _H = _hits_layer.size();
  
  return _H;
}

void HITColl_Layer::printHIT(int i)
{
  HIT * h = hit(i); 
  cout << "HIT n." << i << " information" << endl;
  h->print();
  
  return;
}
