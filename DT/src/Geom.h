#ifndef Geom_h
#define Geom_h

#include <iostream>
#include <math.h>
//#include <algorithm>

// /////////////////////////////////////////////////////////////////////
//                                                                    //
// Silvia Pesente March 2009                                          //
//                                                                    //
// PATTERN RECOGNITION:                                               //
// Class to store GEOMETRY (wires position and layers position)       //
//                                                                    //
// /////////////////////////////////////////////////////////////////////

class Geom {
 public:
  // default constructor
  Geom();

  // destructor
  ~Geom(); 
  
  
  // * return function
  float get_x_wire(int CH, int SL, int L, int W);
  float get_y_wire(int CH, int SL, int L, int W);
  
  void printGeom();
  
 private:
  int m_chtype;
  
};
#endif /*Geom_h*/
