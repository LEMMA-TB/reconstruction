// NEW VERSION MODIFIED BY A. LORENZON 23.08.2017


#include <algorithm>

#include "EventBuilder.h"
#include "TSystem.h"

using namespace std;

EventBuilder::EventBuilder()
{
    m_debug = false;

    detIDs[0] = 10;
    detIDs[1] = 20;
    detIDs[2] = 30;
    detIDs[3] = 40;
    detIDs[4] = 51;
    detIDs[5] = 50;
    detIDs[6] = 56;
    detIDs[7] = 55;
    detIDs[8] = 70;     

    Int_t iLay=0;
    map_detID[10]=iLay++;
    map_detID[20]=iLay++;
    map_detID[30]=iLay++;
    map_detID[40]=iLay++;
    map_detID[51]=iLay++;
    map_detID[50]=iLay++;
    map_detID[56]=iLay++;
    map_detID[55]=iLay++;
    map_detID[70]=iLay++;
    
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void EventBuilder::setDebug( bool debug ){
    m_debug = debug;
    cout << "Debug flag  " << m_debug << endl;
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void EventBuilder::setAlignment( bool doAli ){
    doAlignment = doAli;
    cout << "Set Local->Global alignment flag  " << doAlignment << endl;
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void EventBuilder::openDataFiles( string inputDTFile, string inputSiFile, string outputFile ){

    // Open DT root file
    m_dtFile = new TFile( ( inputDTFile ).c_str() ); //ALTEA 
    cout << "Opening mu data file " << ( inputDTFile ).c_str() << endl;
    
    m_tree = ( TTree* )m_dtFile->Get( "RADMU" );

    m_tree->SetBranchAddress( "EVENT"  , &dtEvent, &b_EVENT   );
    m_tree->SetBranchAddress( "SEG_ns" , &nseg   , &b_ntes    );
    m_tree->SetBranchAddress( "SEG_sn" ,  segN   , &b_SEG_sn  );
    m_tree->SetBranchAddress( "SEG_xh1",  s1p    , &b_SEG_xh1 );
    m_tree->SetBranchAddress( "SEG_xh2",  s2p    , &b_SEG_xh2 );
    m_tree->SetBranchAddress( "SEG_xh3",  s3p    , &b_SEG_xh3 );
    m_tree->SetBranchAddress( "SEG_xh4",  s4p    , &b_SEG_xh4 );
    m_tree->SetBranchAddress( "SEG_xh5",  s5p    , &b_SEG_xh5 );
    m_tree->SetBranchAddress( "SEG_xh6",  s6p    , &b_SEG_xh6 );
    m_tree->SetBranchAddress( "SEG_xh7",  s7p    , &b_SEG_xh7 );
    m_tree->SetBranchAddress( "SEG_xh8",  s8p    , &b_SEG_xh8 );
    
    m_idt = 0;
    
    // Open Si root file;

    m_siFile = new TFile( ( inputSiFile ).c_str() );
    cout << "Opening Si data file " << ( inputSiFile ).c_str() << endl;
    
    t2 = ( TTree* )m_siFile->Get( "LEMMA" );
    
    t2->SetBranchAddress( "iev",                       &si_iev,                       &b_iev                       );
    t2->SetBranchAddress( "iev_mu",                    &si_ievmu,                     &b_ievmu                     );
    t2->SetBranchAddress( "nhits",                     &si_nhits,                     &b_nhits                     );
    t2->SetBranchAddress( "subdet",                     si_subdet,                    &b_subdet                    );
    t2->SetBranchAddress( "xh",                         si_xh,                        &b_xh                        );
    t2->SetBranchAddress( "xh_PulseHeight",             si_xh_PulseHeight,            &b_xh_PulseHeight            );
    t2->SetBranchAddress( "yh",                         si_yh,                        &b_yh                        );
    t2->SetBranchAddress( "yh_PulseHeight",             si_yh_PulseHeight,            &b_yh_PulseHeight            );
    t2->SetBranchAddress( "zh",                         si_zh,                        &b_zh                        );
    t2->SetBranchAddress( "itrack",                     si_itrack,                    &b_itrack                    );
    t2->SetBranchAddress( "S2_PulseHeight",            &si_S2_PulseHeight,            &b_S2_PulseHeight            );
    t2->SetBranchAddress( "LemmaTrigger_PulseHeight",  &si_LemmaTrigger_PulseHeight,  &b_LemmaTrigger_PulseHeight  );   
    t2->SetBranchAddress( "Deva_PulseHeight",           si_Deva_PulseHeight,          &b_Deva_PulseHeight          );
    t2->SetBranchAddress( "Cherenkov_PulseHeight",      si_Cherenkov_PulseHeight,     &b_Cherenkov_PulseHeight     );
    t2->SetBranchAddress( "S3_PulseHeight",            &si_S3_PulseHeight,            &b_S3_PulseHeight            );
    t2->SetBranchAddress( "PbGlass_PulseHeight",       &si_PbGlass_PulseHeight,       &b_PbGlass_PulseHeight       );
    t2->SetBranchAddress( "S2_TimeSpectrum",           &si_S2_TimeSpectrum,           &b_S2_TimeSpectrum           );
    t2->SetBranchAddress( "LemmaTrigger_TimeSpectrum", &si_LemmaTrigger_TimeSpectrum, &b_LemmaTrigger_TimeSpectrum );
    t2->SetBranchAddress( "Deva_TimeSpectrum",          si_Deva_TimeSpectrum,         &b_Deva_TimeSpectrum         );
    t2->SetBranchAddress( "Cherenkov_TimeSpectrum",     si_Cherenkov_TimeSpectrum,    &b_Cherenkov_TimeSpectrum    );
    t2->SetBranchAddress( "S3_TimeSpectrum",           &si_S3_TimeSpectrum,           &b_S3_TimeSpectrum           );
    t2->SetBranchAddress( "PbGlass_TimeSpectrum",      &si_PbGlass_TimeSpectrum,      &b_PbGlass_TimeSpectrum      );         

    m_isi = 0;

    // Output file
    m_outFile = new TFile(outputFile.c_str(),"recreate");
    m_outTree = new TTree("LEMMA","Event Builder LEMMA data output");
    
    m_outTree->Branch( "siiev",                     &siiev,                           "siiev/I"                                                             );
    m_outTree->Branch( "iev",                       &iev,                             "iev/I"                                                               );
    m_outTree->Branch( "nhits",                     &nhits,                           "nhits/I"                                                             );
    m_outTree->Branch( "subdet",                    &subdet,                          "subdet[100]/I"                                                       );
    m_outTree->Branch( "xh",                        &xh,                              "xh[100]/F"                                                           );
    m_outTree->Branch( "xh_PulseHeight",             xh_PulseHeight,                  "xh_PulseHeight[100]/F"                                               );
    m_outTree->Branch( "yh",                        &yh,                              "yh[100]/F"                                                           );
    m_outTree->Branch( "yh_PulseHeight",             yh_PulseHeight,                  "yh_PulseHeight[100]/F"                                               );
    m_outTree->Branch( "zh",                        &zh,                              "zh[100]/F"                                                           );
    m_outTree->Branch( "itrack",                    &itrack,                          "itrack[100]/I"                                                       );
    m_outTree->Branch( "S2_PulseHeight",            &S2_PulseHeight,                  "S2_PulseHeight/I"                                                    );
    m_outTree->Branch( "LemmaTrigger_PulseHeight",  &LemmaTrigger_PulseHeight,        "LemmaTrigger_PulseHeight/I"                                          );   
    m_outTree->Branch( "Deva_PulseHeight",           Deva_PulseHeight,                "Deva_PulseHeight[6]/I"                                               );
    m_outTree->Branch( "Cherenkov_PulseHeight",      Cherenkov_PulseHeight,           "Cherenkov_PulseHeight[4]/I"                                          );
    m_outTree->Branch( "S3_PulseHeight",            &S3_PulseHeight,                  "S3_PulseHeight/I"                                                    );
    m_outTree->Branch( "PbGlass_PulseHeight",       &PbGlass_PulseHeight,             "PbGlass_PulseHeight/I"                                               );
    m_outTree->Branch( "S2_TimeSpectrum",           &S2_TimeSpectrum,                 "S2_TimeSpectrum/I"                                                   );
    m_outTree->Branch( "LemmaTrigger_TimeSpectrum", &LemmaTrigger_TimeSpectrum,       "LemmaTrigger_TimeSpectrum/I"                                         );
    m_outTree->Branch( "Deva_TimeSpectrum",          Deva_TimeSpectrum,               "Deva_TimeSpectrum[6]/I"                                              );
    m_outTree->Branch( "Cherenkov_TimeSpectrum",     Cherenkov_TimeSpectrum,          "Cherenkov_TimeSpectrum[4]/I"                                         );
    m_outTree->Branch( "S3_TimeSpectrum",           &S3_TimeSpectrum,                 "S3_TimeSpectrum/I"                                                   );
    m_outTree->Branch( "PbGlass_TimeSpectrum",      &PbGlass_TimeSpectrum,            "PbGlass_TimeSpectrum/I"                                              );         
    
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void EventBuilder::matchEvents( int nEvents ){
    
    int nsi=0;
    int ndt=0;
    
    getDtEvent(); //why?
    getDtEvent();
    
    while( m_idt<min( nEvents, int( m_tree->GetEntries() ) ) && int( m_tree->GetEntries() ) &&
           m_isi<min( nEvents, int( t2    ->GetEntries() ) ) && int( t2    ->GetEntries() )    ) {
        
        if ( m_idt%1000 ) cout << "\n --- Matching event " << m_idt << std::endl;
        
        int dtN = getDtEvent(); 
        int siN = getSiEvent();
    
        while( siN < dtN ) {
            siN = getSiEvent();
            nsi++;

            if ( nsi>int( t2->GetEntries() ) || ndt>int( m_tree->GetEntries() ) ) break;
        }
    
        while( siN > dtN ) {
            dtN = getDtEvent();
            ndt++;
            
            if ( nsi>int( t2->GetEntries() ) || ndt>int( m_tree->GetEntries() ) ) break;
        }
    
        if ( siN == dtN ){
            if ( m_idt%1000 ) cout << "\033[1;31mEvent " << dtN << " matched! \033[0m" << std::endl;
            dumpGlobalEvent();
        }
    }
    
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
int EventBuilder::getDtEvent(){
    
    m_tree->GetEntry( m_idt );
    m_idt++;
    return dtEvent;
    
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
int EventBuilder::getSiEvent(){
    
    t2->GetEntry( m_isi );
    m_isi++;
    return si_ievmu;
    
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void EventBuilder::calibrate(){
    
    /// TODO COMPLETE!
    
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void EventBuilder::openAlignments(std::string alignment_file){
    std::cout << "Opening alignment file: " << alignment_file << std::endl;
    
    //   float x0SiAndMu[nSiLayers+1][3];
    //   float phi0Theta0Mu        [2];
    
    std::string line;
    ifstream myfile (alignment_file.c_str());
    std::vector<float> vAli;
    if (myfile.is_open())
    {
        while( getline(myfile, line) )
        {
            std::stringstream ss(line); 
            std::string token;
            while (std::getline(ss, token, ' ')){
                if (!TString(token).IsFloat()) {
                    std::cout << "\033[1;31mERROR. Non float entry in alignment file\033[0m" << std::endl;
                    return;
                }
                vAli.push_back(TString(token).Atof());
            }   
        }
        if ((int)vAli.size() != 3*(nSiLayers+1)+2) {
            std::cout << "\033[1;31mERROR in parsing alignment file\033[0m" << std::endl;
            return;
        }
        for (int iAli=0; iAli<3*(nSiLayers+1); ++iAli) {
            int iLay=iAli/3;
            int xORy=iAli%3;
            x0SiAndMu[iLay][xORy] = vAli.at(iAli);
        }
        phi0Theta0Mu[0] = vAli.at(3*(nSiLayers+1)  );
        phi0Theta0Mu[1] = vAli.at(3*(nSiLayers+1)+1);
        
    }
    
    std::cout << "xh \t yh \t zh" << std::endl;
    for (Int_t iLayer=0; iLayer<nSiLayers+1; iLayer++) {
        std::cout << x0SiAndMu[iLayer][0] << "\t"<< x0SiAndMu[iLayer][1] << "\t"<< x0SiAndMu[iLayer][2] << std::endl;
    }
    std::cout << "phi0mu \t theta0mu" << std::endl;
    std::cout <<  phi0Theta0Mu[0] << "\t" << phi0Theta0Mu[1] << std::endl;
    
    return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventBuilder::dumpGlobalEvent(){
    
    for (int ie=0; ie<100; ie++) {
        
        xh[ie]    = -999;
        yh[ie]    = -999;
        zh[ie]    = -999;
        itrack[ie]=    2;
        subdet[ie]= -999;
        
    }

    /// Si data
    
    int nSiHits = si_nhits;
    
    if ( nSiHits+12>nMaxHits ) { //Does it make sense?
        cout << "\033[1;31m Warning: too many hits. Skipping event. \033[0m" << endl;
        return;
    }
    
    siiev = si_iev;
    
    for (Int_t iSiHit=0; iSiHit<nSiHits; ++iSiHit) {
        subdet[iSiHit] = si_subdet[iSiHit];
        if ( si_xh[iSiHit] == -999 ) { xh[iSiHit] = si_xh[iSiHit]; } else { xh[iSiHit] = -si_xh[iSiHit]; }
        if ( si_yh[iSiHit] == -999 ) { yh[iSiHit] = si_yh[iSiHit]; } else { yh[iSiHit] = -si_yh[iSiHit]; }
        if ( si_zh[iSiHit] == -999 ) { zh[iSiHit] = si_zh[iSiHit]; } else { zh[iSiHit] = -si_zh[iSiHit]; }
        xh_PulseHeight[iSiHit] = si_xh_PulseHeight[iSiHit];
        yh_PulseHeight[iSiHit] = si_yh_PulseHeight[iSiHit];
        itrack [iSiHit] = si_itrack[iSiHit];
    }
    
    S2_PulseHeight            = si_S2_PulseHeight;
    S2_TimeSpectrum           = si_S2_TimeSpectrum;
    LemmaTrigger_PulseHeight  = si_LemmaTrigger_PulseHeight;
    LemmaTrigger_TimeSpectrum = si_LemmaTrigger_TimeSpectrum;
    
    for( Int_t iChannel=0; iChannel<nDevaChannels; ++iChannel ) {
        Deva_PulseHeight [iChannel] = si_Deva_PulseHeight [iChannel];
        Deva_TimeSpectrum[iChannel] = si_Deva_TimeSpectrum[iChannel];
    }
    
    for( Int_t iChannel=0; iChannel<nCherenkovChannels; ++iChannel ) {
        Cherenkov_PulseHeight [iChannel] = si_Cherenkov_PulseHeight [iChannel];
        Cherenkov_TimeSpectrum[iChannel] = si_Cherenkov_TimeSpectrum[iChannel];        
    }
    
    S3_PulseHeight       = si_S3_PulseHeight;
    S3_TimeSpectrum      = si_S3_TimeSpectrum;
    PbGlass_PulseHeight  = si_PbGlass_PulseHeight;
    PbGlass_TimeSpectrum = si_PbGlass_TimeSpectrum;

    /// DT data
    iev   = dtEvent;
    nhits = nSiHits + 12;
    
    // set default for variables x,y, subdet
    for(int ih=nSiHits; ih<12 + nSiHits; ih++){
        xh[ih]     = -999;
        yh[ih]     = -999;
        zh[ih]     = -999;
        subdet[ih] =   70;
    }
    
    // set z of the layers
    float z_layers[12] = {-10.75, -9.45, -8.15, -6.85,  7.45, 8.75, 10.05, 11.35,12.85, 14.15, 15.45, 16.75};
    
    // fill HITS - Z
    for(int ih=0 ; ih<12 ; ih++) {
        zh[ih+nSiHits] = z_layers[ih];
    }
    
    for ( int is = 0; is < nseg; is++ ){
        if( segN[is]<0 ){
            // SL THETA hits
            yh[4 + nSiHits] = s1p[is];
            yh[5 + nSiHits] = s2p[is];
            yh[6 + nSiHits] = s3p[is];
            yh[7 + nSiHits] = s4p[is];
        } else {
            // SL PHI hits
            if ( s1p[is]!=-999 ) { xh[ 0 + nSiHits] = -s1p[is]; } else { xh[ 0 + nSiHits] = s1p[is]; }
            if ( s2p[is]!=-999 ) { xh[ 1 + nSiHits] = -s2p[is]; } else { xh[ 1 + nSiHits] = s2p[is]; }
            if ( s3p[is]!=-999 ) { xh[ 2 + nSiHits] = -s3p[is]; } else { xh[ 2 + nSiHits] = s3p[is]; }
            if ( s4p[is]!=-999 ) { xh[ 3 + nSiHits] = -s4p[is]; } else { xh[ 3 + nSiHits] = s4p[is]; }
            if ( s5p[is]!=-999 ) { xh[ 8 + nSiHits] = -s5p[is]; } else { xh[ 8 + nSiHits] = s5p[is]; }
            if ( s6p[is]!=-999 ) { xh[ 9 + nSiHits] = -s6p[is]; } else { xh[ 9 + nSiHits] = s6p[is]; }
            if ( s7p[is]!=-999 ) { xh[10 + nSiHits] = -s7p[is]; } else { xh[10 + nSiHits] = s7p[is]; }
            if ( s8p[is]!=-999 ) { xh[11 + nSiHits] = -s8p[is]; } else { xh[11 + nSiHits] = s8p[is]; }
        }
    }
    
    // Apply alignment to get global coordinates
    
    if ( doAlignment )
        for (Int_t ih = 0; ih < nhits; ++ih) {
            std::map<Int_t,Int_t>::iterator it_map_detID;
            it_map_detID = map_detID.find( subdet[ih] );
            Int_t iLayer=-1;
            if (it_map_detID != map_detID.end()) iLayer = it_map_detID->second;

            xh[ih] += x0SiAndMu[iLayer][0];
            yh[ih] += x0SiAndMu[iLayer][1];
            zh[ih] += x0SiAndMu[iLayer][2];
        }
        
        
        if( m_debug ){
            cout << "---- EventBuilder::dumpGlobalEvent " << endl;
            cout << "Filling DT data..." << endl;
            cout << " ---> event ID " << dtEvent << endl;
            cout << " ---> muon track nhits " << nhits << endl;
        }
        
        /// fill tree
        m_outTree->Fill();
        
        return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventBuilder::dumpOutput(){
    
    m_outTree->Write();
    m_outFile->Close();
    
    return;
}
