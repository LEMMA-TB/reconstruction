#include "Geom.h"


bool debug_geom  = false;


Geom::Geom()
{   
    // MB3 default for LNL stand
    m_chtype = 3;

    // MB2 LEMMa tb
    m_chtype = 2;

  return;
}

Geom::~Geom()
{
    return;
}

float Geom::get_x_wire(int CH, int SL, int L, int W)
{

  float _x = 0;

//  // LNL MB3 reference system
//  float _x0_phi_lay1 = 148.6;
//  float _x0_phi_lay2 = 150.7;
//  float _x0_theta_lay1 =117.35;
//  float _x0_theta_lay2 = 119.45;

  // LEMMA reference system
  float _x0_phi_lay1 = 120.;
  float _x0_phi_lay2 = _x0_phi_lay1 + 2.1;

  float _x0_theta_lay1 =117.35;
  float _x0_theta_lay2 = _x0_theta_lay1 + 2.1;

  /// PHI SL
  // chamber 11 --> LEMMA tb and LNL upper chamber
  if( (CH==11 && (SL==1||SL==3)) || CH==9 ||CH==8){
    if(L==1 || L==3)
      _x= (W-1)*4.2-_x0_phi_lay1;
    if(L==2 || L==4)
      _x= (W-1)*4.2-_x0_phi_lay2;

    // add SL phi1 staggering
    if(m_chtype==2 && SL==1)
        _x += 4.2;

    // X axis overturn
    if(m_chtype==2)
        _x = - _x;
    }

  /// THETA SL
  if( (CH==11||CH==10) && (SL==2) ){
    if(L==1 || L==3)
      _x= _x0_theta_lay1 - (W-1)*4.2;
    if(L==2 || L==4)
      _x= _x0_theta_lay2 - (W-1)*4.2;

    // X axis overturn
    if(m_chtype==2)
        _x = - _x;
  }

  // chamber 10 --> LNL lower chamber
  if( CH==10 && (SL==1||SL==3) ){
    if(L==1 || L==3)
      _x=-(W-1)*4.2+_x0_phi_lay1;
    if(L==2 || L==4)
      _x=-(W-1)*4.2+_x0_phi_lay2;
  }


  //std::cout << "X wire for ch " << CH << ", SL " << SL << ", LAY " << L << ", wire " << W << " = " << _x << std::endl;
  return _x;
}
  
float Geom::get_y_wire(int CH, int SL, int L, int W)
{
  float _y = 0;
//  float dy_ch1 = 0.1581;	//chamber 1 honeycomb correction for MB3 @ LNL
  float dy_ch1 = 0.;

  float ay_phi_SL1[4] = {10.75, 9.45, 8.15, 6.85};
  float ay_phi_SL3[4] = {12.85, 14.15, 15.45, 16.75};
  float ay_theta[4] = {7.45, 8.75, 10.05, 11.35};
  float ay_phi_34[4] = {-1.95,-0.65,0.65,1.95};
  
  if(CH==11 && SL==1)
    for(int i=0;i<4;i++) _y=-ay_phi_SL1[L-1]-dy_ch1;
  if(CH==11 && SL==3)
    for(int i=0;i<4;i++) _y=ay_phi_SL3[L-1];
  if(CH==11 && SL==2)
    for(int i=0;i<4;i++) _y=ay_theta[L-1];

  if(CH==10 && SL==1)
    for(int i=0;i<4;i++) _y=ay_phi_SL1[L-1];
  if(CH==10 && SL==3)
    for(int i=0;i<4;i++) _y=-ay_phi_SL3[L-1];
  if(CH==10 && SL==2)
    for(int i=0;i<4;i++) _y=-ay_theta[L-1];

  if( CH==8||CH==9 )
    for(int i=0;i<4;i++) _y=ay_phi_34[L-1];
  
  return _y;
}

void Geom::printGeom()
{}

