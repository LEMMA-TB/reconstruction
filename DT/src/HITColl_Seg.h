#ifndef HITColl_Seg_h
#define HITColl_Seg_h

#include <iostream>
#include <vector>

#include "TObject.h"

#include "HITCollection.h"


// /////////////////////////////////////////////////////////////////////
//                                                                    //
// Silvia Pesente March 2009                                          //
//                                                                    //
// PATTERN RECOGNITION:                                               //
// Class to store HITs of each segment                                //
//                                                                    //
// /////////////////////////////////////////////////////////////////////


class HITColl_Seg : public TObject {
 public:
  // default constructor
  HITColl_Seg() { } 
  
  ~HITColl_Seg();
  
  // operations on HITS
  void selectHIT(HITCollection * hits, int CH, bool phi); 
  void selectHIT(HITCollection * hits, int CH, int SL); 
  void addHIT(HIT * hit);
  void eraseHIT(int i);
  void printHIT(int i);
  HIT * hit(int i);
  
  // operations on collection
  int Get_NumberHITS();
  
 private:
  // store hit objects in array
  vector< HIT *> _hits_seg;
  
};

#endif /*HITColl_Seg_h*/

