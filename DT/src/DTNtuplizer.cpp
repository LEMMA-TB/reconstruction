#include "DTNtuplizer.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

DTNtuplizer::DTNtuplizer(char * fileName) : m_tree(0)
{
    // open file
    cout << "Opening file for output  " << fileName << endl;
    m_file = new TFile(fileName,"RECREATE");
    m_tree = new TTree("DTHITS","dt hits for LEMMA analysis");

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

   m_tree->Branch("EventID",       &EventID,        "EventID/I");
   m_tree->Branch("DT_NHits",    &NHits,         "NHits/I");
   m_tree->Branch("DT_lay",         lay,                 "lay[NHits]/I");
   m_tree->Branch("DT_tube",      tube,               "tube[NHits]/I");
   m_tree->Branch("DT_dtime",    dtime,             "dtime[NHits]/F");
   m_tree->Branch("DT_Xwire",    Xwire,             "Xwire[NHits]/F");
   m_tree->Branch("DT_Zwire",    Zwire,             "Zwire[NHits]/F");

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
      EventID = -999.;
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

    EventID = numEvent;
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
