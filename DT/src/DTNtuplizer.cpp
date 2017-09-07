#include "DTNtuplizer.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

DTNtuplizer::DTNtuplizer(char * fileName) : m_tree(0)
{
    // open file
    cout << "Opening file for output  " << fileName << endl;
    m_file = new TFile(fileName,"RECREATE");
    m_tree = new TTree("LEMMA-DT","dt hits for LEMMA analysis");

    m_geom = new Geom();

   Init();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

DTNtuplizer::~DTNtuplizer()
{
//   if (!m_tree) return;
//   delete m_tree->GetCurrentFile();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DTNtuplizer::Init()
{
   // Set branch addresses and branch pointers
   m_tree->SetMakeClass(1);

   m_tree->Branch("EVENT",        &EVENT,      "EVENT/I");
   m_tree->Branch("DT_NHits",    &NHits,         "NHits/I");
   m_tree->Branch("DT_lay",         lay,                 "lay[NHits]/I");
   m_tree->Branch("DT_tube",      tube,               "tube[NHits]/I");
   m_tree->Branch("DT_dtime",    dtime,             "dtime[NHits]/F");
   m_tree->Branch("DT_Xwire",    Xwire,             "Xwire[NHits]/F");
   m_tree->Branch("DT_Zwire",    Zwire,             "Zwire[NHits]/F");

  m_tree->Branch( "SEG_ns",         &NSeg,  "NSeg/I" );
  m_tree->Branch( "SEG_X",          segX,     "segX[NSeg]/F" );
  m_tree->Branch( "SEG_slope",   segS,     "segS[NSeg]/F" );
  m_tree->Branch( "SEG_chi",        segK,     "segK[NSeg]/F" );
  m_tree->Branch( "SEG_num",     segN,     "segN[NSeg]/I" );
  m_tree->Branch( "SEG_t0",         segT0,    "segT0[NSeg]/F" );
  m_tree->Branch( "SEG_1r",         l1r,   "l1r[NSeg]/F" );
  m_tree->Branch( "SEG_2r",         l2r,   "l2r[NSeg]/F" );
  m_tree->Branch( "SEG_3r",         l3r,   "l3r[NSeg]/F" );
  m_tree->Branch( "SEG_4r",         l4r,   "l4r[NSeg]/F" );
  m_tree->Branch( "SEG_5r",         l5r,   "l5r[NSeg]/F" );
  m_tree->Branch( "SEG_6r",         l6r,   "l6r[NSeg]/F" );
  m_tree->Branch( "SEG_7r",         l7r,   "l7r[NSeg]/F" );
  m_tree->Branch( "SEG_8r",         l8r,   "l8r[NSeg]/F" );

  m_tree->Branch( "SEG_xh1",    xh1,   "xh1[NSeg]/F" );
  m_tree->Branch( "SEG_xh2",    xh2,   "xh2[NSeg]/F" );
  m_tree->Branch( "SEG_xh3",    xh3,   "xh3[NSeg]/F" );
  m_tree->Branch( "SEG_xh4",    xh4,   "xh4[NSeg]/F" );
  m_tree->Branch( "SEG_xh5",    xh5,   "xh5[NSeg]/F" );
  m_tree->Branch( "SEG_xh6",    xh6,   "xh6[NSeg]/F" );
  m_tree->Branch( "SEG_xh7",    xh7,   "xh7[NSeg]/F" );
  m_tree->Branch( "SEG_xh8",    xh8,   "xh8[NSeg]/F" );

  Notify();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t DTNtuplizer::Notify()
{
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Int_t DTNtuplizer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DTNtuplizer::cleanTree() {

  for (int ih = 0; ih < NHits; ih++ ){
      EVENT = -999.;
      lay[ih] = -999;
      tube[ih] = -999;
      dtime[ih] = -999.;
      Xwire[ih] = -999.;
      Zwire[ih] = -999.;
  }

  NHits = -999;

  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DTNtuplizer::fillHits(HITCollection *hits, int numEvent){

    cleanTree();

    EVENT = numEvent;
    NHits = hits->Get_NumberHITS();

    for(int ih=0; ih<NHits; ih++){

        HIT *hit=hits->hit(ih);

        lay[ih]        = hit->L_ID() + (4*(hit->SL_ID() - 1));
        tube[ih]     = hit->wire_ID();
        dtime[ih]      = hit->dtime();
        Xwire[ih]           = m_geom->get_x_wire(11, hit->SL_ID(), hit->L_ID(), hit->wire_ID());
        Zwire[ih]           = m_geom->get_y_wire(11, hit->SL_ID(), hit->L_ID(), hit->wire_ID());
    }

    m_tree->Fill();

    return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DTNtuplizer::write(){

    m_file->cd();
    m_tree->Write();

    return;
}
