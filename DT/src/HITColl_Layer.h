#ifndef HITColl_Layer_h
#define HITColl_Layer_h

#include <iostream>
#include <vector>

#include "TObject.h"

#include "HITColl_Seg.h"


// /////////////////////////////////////////////////////////////////////
//                                                                    //
// Silvia Pesente March 2009                                          //
//                                                                    //
// PATTERN RECOGNITION:                                               //
// Class to store HITs of each layer                                  //
//                                                                    //
// /////////////////////////////////////////////////////////////////////


class HITColl_Layer : public TObject {
 public:
  // default constructor
  HITColl_Layer() { } 
  
  ~HITColl_Layer();
  
  // operations on HITS
  void selectHIT(HITColl_Seg * hit_seg, int L); 
  void addHIT(HIT * hit);
  void eraseHIT(int i);
  void printHIT(int i);
  HIT * hit(int i);
  
  // operations on collection
  int Get_NumberHITS();
  
 private:
  // store hit objects in array
  vector< HIT *> _hits_layer;
  
};

#endif /*HITColl_Layer_h*/

