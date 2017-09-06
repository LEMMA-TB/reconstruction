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
    Int_t         EventID;
    static const Int_t kMaxHits = 500;
    Int_t        NHits;
    Int_t        lay[kMaxHits];
    Int_t        tube[kMaxHits];
    Float_t    dtime[kMaxHits];
    Float_t    Xwire[kMaxHits];
    Float_t    Zwire[kMaxHits];

    Geom        *m_geom; //! ponter to geometry
};

#endif

