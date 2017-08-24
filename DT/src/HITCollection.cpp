#include "HITCollection.h"

HITCollection::~HITCollection()
{
  int N=Get_NumberHITS();
  if(DEBUG_HITCOLL){
    printf("\n%d hits into HITCollection:\n\n",N);
    for(int j=0;j<N;j++)
      printHIT(j);
    printf("...Into ~HITCollection() deleting N.%d hits...\n",N);
  }
  for(int h=(N-1);h>=0;h--){
    
    delete _hits[h];
    _hits[h]=NULL;
    
    eraseHIT(h);

  }

  return;
}

void HITCollection::createHIT(int Event, int chamber, int SL, int L, int tube,
			      float x_wire, float y_wire,
			      int TDC, int ROB, int channel,   
			      float raw_time, float t0, float ttrig) 
{
  if(DEBUG_HITCOLL)
    cout << "\n\n *** HITCollection::createHIT" << endl;
  
  pair<bool,bool> code=make_pair(true,true);
  bool isdouble=false;

  // compute drift time in ns
  float drift_time =  T_shift + raw_time - t0 - ttrig;
  
  HIT * hit = new HIT(Event, chamber, SL, L, tube, x_wire, y_wire, TDC, ROB, channel, raw_time, t0, ttrig, drift_time, drift_time, code, isdouble);
  
  addHIT(hit);

  
  return;
}

void HITCollection::addHIT(HIT * hit)
{
  if(DEBUG_HITCOLL)
    cout << "HITCollection::addHIT" << endl;
  
  _hits.push_back(hit);
  
  return;
}


void HITCollection::eraseHIT(int i)
{
  if(DEBUG_HITCOLL)
    cout << "HITCollection::eraseHIT n." << i << endl;

  _hits.erase(_hits.begin()+i);
 
  return;
}

HIT * HITCollection::hit(int i)
{	
  HIT * hit_ptr = _hits[i];
  
  return hit_ptr;
}

void HITCollection::Clean2FE()
{
    for(vector<HIT*>::iterator it = _hits.begin(); it<_hits.end(); it++ )
    {
       HIT *hit1 = *it;
       int n_tube = hit1->wire_ID();
       if(hit1->CH_ID() == 10 && hit1->SL_ID() == 1 && hit1->L_ID() == 4 && n_tube >= 32 && n_tube <= 41) {
           for(vector<HIT*>::iterator it2 = _hits.begin(); it2<_hits.end(); it2++ )
           {
               HIT *hit2 = *it2;
               int n_tube2 = hit2->wire_ID();
               if(hit2->CH_ID() == 10 && hit2->SL_ID() == 1 && hit2->L_ID() == 4 && n_tube2 >= 32 && n_tube2 <= 41 && abs(n_tube-n_tube2) == 1 ) {
                   if(hit1->rtime() < hit2->rtime()) _hits.erase(it2);
               }

           }
       }
    }
}

int HITCollection::Get_NumberHITS()
{
  int _H = _hits.size();
  return _H;
}

void HITCollection::printHIT(int i)
{
  HIT * h = hit(i); 
  cout << "INFO HIT N. " << i << endl;
  h->print();
  
  return;
}

void HITCollection::dumpHITCollection()
{
    cout << "*** HITCollection dump ***" << endl;
    for(int ih=0; ih<_hits.size(); ih++){
        HIT * h = hit(ih);
        h->print();
    }
    return;
}
