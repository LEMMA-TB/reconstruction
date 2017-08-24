#include "HIT.h"


HIT::HIT(){}
HIT::HIT(int Event, int chamber, int SL, int L, int tube,
	 float x_wire, float y_wire,
	 int TDC, int ROB, int channel,   
	 float raw_time, float t0, float ttrig, 
	 float drift_time_in, float drift_time,
	 pair<bool,bool> code, bool isdouble) 
{
  if(DEBUG_HIT)
    cout << "HIT::HIT creating a new HIT " << endl;
  
  _EvNumber = Event;
  _CH_ID = chamber;
  _SL_ID = SL;
  _L_ID  = L;
  _wire_ID = tube;
  _x_wire = x_wire;
  _y_wire = y_wire;
  _TDC_ID = TDC;
  _ROB_ID = ROB;
  _channel = channel;
  _rtime = raw_time;
  _t0 = t0;
  _ttrig = ttrig;
  _dtime_in = drift_time_in;
  _dtime = drift_time;
  _code_ID.first = code.first;
  _code_ID.second = code.second;
  _isdouble = isdouble;
  
  return;
}


HIT::~HIT()
{   
  return;
}

pair<float,float> HIT::Get_XPair(float time_ns){
 
  float x_left= x_wire_ID() - time_ns*0.00547;
  float x_right= x_wire_ID() + time_ns*0.00547;

  return make_pair(x_left,x_right);
}

bool HIT::Change_LCode(bool xleft)
{
  return _code_ID.first=xleft;  
}

bool HIT::Change_RCode(bool xright)
{
  return _code_ID.second=xright;  
}

bool HIT::IsSolved()
{
  bool is_solved=false;
  if( (_code_ID.first==true &&_code_ID.second==false)
      || (_code_ID.first==false &&_code_ID.second==true) )
    is_solved=true;
  
  return is_solved;  
}

bool HIT::Set_IsDouble(bool isdouble)
{
  return _isdouble=isdouble;  
}

pair<bool,bool> HIT::Get_Code_ID()
{
  return make_pair(_code_ID.first,_code_ID.second);  
}

float HIT::Set_CorrTime(float corr_time)
{
  _dtime=corr_time;
  return _dtime;  
}

void HIT::print()
{
  cout <<"Ev: "<<EvNumber()<<", CH: "<<CH_ID()<<", SL: "<<SL_ID()
       <<", L: "<<L_ID()<<", wire: "<<wire_ID()
       <<", x_wire: "<< x_wire_ID()<<", y_wire: "<< y_wire_ID()
       <<", raw time (ns): "<<rtime()<<", t0 (ns): "<<t0()<<", ttrig (ns): "<<ttrig()
       <<", drift time iniziale (ns): "<<dtime_in()
       <<", drift time fin (ns): "<<dtime()
       <<"\nTDC: "<<TDC_ID()<<", ROB: "<<ROB_ID()<<", channel: "<<channel_ID()<<"\n" << endl;
  
  return;
}

