////////////////////////////////////////////////////////////////////////////
/// 20170905 Sara Vanini DT LEMMA-TB ntuple
///////////////////////////////////////////////////////////////////////////

#ifndef DTNtuple_h
#define DTNtuple_h

#include <strings.h>

#include "TFile.h"
#include "TTree.h"

#include "HITCollection.h"
#include "Geom.h"

class DTNtuplizer {
public:
    DTNtuplizer(char *fileName=0);
    ~DTNtuplizer();
    void cleanTree();
    void    Init();
    Bool_t   Notify();

    void    fillHits(HITCollection *hits, int numEvent);
    void    write();

    Int_t    Cut(Long64_t entry);

public:
    TFile          *m_file;    //!pointer to the TFile
    TTree        *m_tree;   //!pointer to the TTree

    // Declaration of leaf types
    Int_t         EVENT;
    static const Int_t kMaxHits = 500;
    Int_t        NHits;
    Int_t        lay[kMaxHits];
    Int_t        tube[kMaxHits];
    Float_t    dtime[kMaxHits];
    Float_t    Xwire[kMaxHits];
    Float_t    Zwire[kMaxHits];

    static const Int_t kMaxSeg = 10;
    Int_t         NSeg;
    Float_t     segX[kMaxSeg];
    Float_t     segS[kMaxSeg];
    Float_t     segK[kMaxSeg];
    Float_t     segN[kMaxSeg];
    Float_t     segT0[kMaxSeg];
    Float_t     l1r[kMaxSeg];
    Float_t     l2r[kMaxSeg];
    Float_t     l3r[kMaxSeg];
    Float_t     l4r[kMaxSeg];
    Float_t     l5r[kMaxSeg];
    Float_t     l6r[kMaxSeg];
    Float_t     l7r[kMaxSeg];
    Float_t     l8r[kMaxSeg];

    Float_t     xh1[kMaxSeg];
    Float_t     xh2[kMaxSeg];
    Float_t     xh3[kMaxSeg];
    Float_t     xh4[kMaxSeg];
    Float_t     xh5[kMaxSeg];
    Float_t     xh6[kMaxSeg];
    Float_t     xh7[kMaxSeg];
    Float_t     xh8[kMaxSeg];

    Geom        *m_geom; //! ponter to geometry
};

#endif

