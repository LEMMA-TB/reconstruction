#ifndef HIT_h
#define HIT_h

#include <iostream>
#include "TVectorD.h"
#include "TMath.h"
#include <math.h>
#include <utility>
#include "Constants.h"
#include "Debugs.h"
using namespace std;

// /////////////////////////////////////////////////////////////////////
//                                                                    //
// Silvia Pesente March 2009                                          //
//                                                                    //
// PATTERN RECOGNITION:                                               //
// Class to store single hit and get info on hit                      //
//                                                                    //
// /////////////////////////////////////////////////////////////////////


class HIT {
 public:
  // default constructor
  HIT();
  
  HIT(int Event, int chamber, int SL, int L, int tube,
      float x_wire, float y_wire,
      int TDC, int ROB, int channel,
      float raw_time, float t0, float ttrig,
      float drift_time_in, float drift_time, 
      pair <bool,bool> code, bool isdouble); 
  
  // destructor
  ~HIT(); 
  
  
  // return function
  inline int EvNumber()            {return _EvNumber;}
  inline int CH_ID()               {return _CH_ID;}
  inline int SL_ID()               {return _SL_ID;}
  inline int L_ID()                {return _L_ID;}
  inline int wire_ID()             {return _wire_ID;}
  inline float x_wire_ID()         {return _x_wire;}
  inline float y_wire_ID()         {return _y_wire;}
  inline int TDC_ID()              {return _TDC_ID;}
  inline int ROB_ID()              {return _ROB_ID;}
  inline int channel_ID()          {return _channel;}
  inline float rtime()             {return _rtime;}
  inline float t0()                {return _t0;}
  inline float ttrig()             {return _ttrig;}
  inline float dtime_in()          {return _dtime_in;}
  inline float dtime()             {return _dtime;}
  inline pair<bool,bool> code_ID() {return _code_ID;}
  inline bool IsDouble()           {return _isdouble;}
  
  pair<float,float> Get_XPair(float time_ns);
  bool Change_LCode(bool xleft);
  bool Change_RCode(bool xright);
  bool IsSolved();
  bool Set_IsDouble(bool isdouble);
  pair<bool,bool> Get_Code_ID();
  float Set_CorrTime(float corr_time);
  
  void print();
  
 private:
  int _EvNumber;
  int _CH_ID;   
  int _SL_ID;   
  int _L_ID;   
  int _wire_ID;   
  float _x_wire;
  float _y_wire;
  int _TDC_ID;
  int _ROB_ID;
  int _channel;
  float _rtime;
  float _t0;
  float _ttrig;
  float _dtime_in;
  float _dtime;
  pair<bool,bool> _code_ID;
  bool _isdouble;

};
#endif /*HIT_h*/
