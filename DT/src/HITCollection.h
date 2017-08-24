#ifndef HITCollection_h
#define HITCollection_h

#include <iostream>
#include <vector>

#include "TObject.h"

#include "HIT.h"

// /////////////////////////////////////////////////////////////////////
//                                                                    //
// Silvia Pesente March 2009                                          //
//                                                                    //
// PATTERN RECOGNITION:                                               //
// Class to store all the HITs of each event                          //
//                                                                    //
// /////////////////////////////////////////////////////////////////////


class HITCollection : public TObject {
 public:
  // default constructor
  HITCollection() { } 
  
  ~HITCollection();
  
  // operations on HITS
  void createHIT(int Event, int chamber, int SL, int L, int tube,
		 float x_wire, float y_wire,
		 int TDC, int ROB, int channel,   
		 float raw_time, float t0, float ttrig); 
  void addHIT(HIT * hit);
  void eraseHIT(int i);
  void printHIT(int i);
  void dumpHITCollection();
  HIT * hit(int i);

  // patch per pulire doppia lettura front end
  void Clean2FE();
  
  // operations on collection
  int Get_NumberHITS();
  
 private:
  // store hits objects in array
  vector< HIT *> _hits;
  
};

#endif /*HITCollection_h*/

