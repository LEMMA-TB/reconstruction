/**************************************************************
 *                       treeSiliconBig.C                     *
 * ************************************************************
 * Adapted by A. Lorenzon from treeSilicon3.C (author R. Rossin)
 * 
 * OBJECT
 * This macro reads data file of Silicon, Scintillators and 
 * Calorimeters Detectors like "run[SiRunNumber]_big.dat"
 * and store in a root file only information from Silicon
 * Detectors.
 * 
 * DATA STRUCTURE
 * A Record is made by 211 words organized in this sequence:
 * 
 * Line 1: 35 words    ->  1x Si Event Number
 *                     ->  1x Mu Event Number (max 65536)
 *                     ->  1x Mu Event NUmber (ri calculated)
 *                     ->  14x Pulse Height:   -->  1x S2 Scintillator
 *                                             -->  1x LEMMA trigger
 *                                             -->  6x Deva Calorimeter
 *                                             -->  4x Cherenkov Calorimeter
 *                                             -->  1x S3 Scintillator
 *                                             -->  1x Pb Glass
 *                     ->  2x not used
 *                     ->  14x Time Spectra    -->  1x S2 Scintillator
 *                                             -->  1x LEMMA trigger
 *                                             -->  6x Deva Calorimeter
 *                                             -->  4x Cherenkov Calorimeter
 *                                             -->  1x S3 Scintillator
 *                                             -->  1x Pb Glass
 *                     ->  2x not used
 *
 * Line 2:        Silicon Detector 10   11 words    ->  1x Number of clusters (max 5)
 *                                                  ->  5x Hit position (cm) x view
 *                                                  ->  5x Pulse Height
 * 
 * Line 3:        Silicon Detector 10   11 words    ->  1x Number of clusters (max 5)
 *                                                  ->  5x Hit position (cm) y view
 *                                                  ->  5x Pulse Height
 * 
 * Lines  4 -  5: Silicon Detector 20
 * Lines  6 -  7: Silicon Detector 30
 * Lines  8 -  9: Silicon Detector 40
 * Lines 10 - 11: Silicon Detector 51
 * Lines 12 - 13: Silicon Detector 50
 * Lines 14 - 15: Silicon Detector 56
 * Lines 16 - 17: Silicon Detector 55
 * 
 * TO COMPILE/RUN
 * $ root -l -q treeSiliconBig.C+(SiRunNumber, MuRunNumber)
 * 
 * The output file will have MuRunNumber in the name, in order to proceed 
 * with Event Building; if there is no muon chamber run associated, then 
 * MuRunNumber=0 and the output file will have SiRunNumber in the name.
 * 
 */

#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

// --- Initialize global variables
const Int_t nSiLayers           = 8;                                           // Si layers
const Int_t nDigis              = 32;                                          // 16 pulse height + 16 time spectrum
const Int_t nMaxClusters        = 5;                                           // maximum number of clusters recorded
const Int_t nTokensInRecord     = 1+1+1+nDigis+2*nSiLayers*(1+2*nMaxClusters); // = 211 = number of words recored in one Si event
const Int_t nMaxHits            = 100;
const float epsilon             = 1E-12;
const Int_t detIDs[nSiLayers+1] = {10,20,30,40,51,50,56,55,70};                // array with detector ID (50 55 mu-) (51 56 mu+)
float       alignments[nSiLayers+1][2];

// --- Define structures
typedef struct {
    
    Int_t    eventID_Si; 
    Int_t    eventID_mu;
    Int_t    eventID_mu_expanded;
    Int_t    digis1[nDigis]; 
    Int_t    nxClusters[nSiLayers]; 
    Float_t  x[nSiLayers][nMaxClusters];
    Float_t  xPulseHeight[nSiLayers][nMaxClusters];
    Int_t    nyClusters[nSiLayers];
    Float_t  y[nSiLayers][nMaxClusters];
    Float_t  yPulseHeight[nSiLayers][nMaxClusters];
    
} Si_t; 

typedef struct {
    
    float x0SiAndMu[nSiLayers+1][3];
    float phi0Theta0Mu          [2];
    
} Alignment_t;

typedef struct {
    
    Int_t    subdet;
    Float_t  xh;
    Float_t  yh;
    Float_t  zh;
    Int_t    itrack;
    
} Hit_t;

// --- Define functions
int  loadAlignments  ( Alignment_t &Ali, bool doAlignment, bool m_debug );
void applyAlignment  ( vector<Hit_t> &vhit, const Si_t &SiEvent, Alignment_t &Ali );
void printAlignments ( Alignment_t &Ali );
int  readSievents    ( Int_t SiRunN, Int_t MuRunN, Alignment_t Ali, bool m_debug );
void fillEvent       ( float *vRecord, Si_t &SiEvent );
void printEvent      ( Si_t &Ev );
void printHit        ( Hit_t &hit );

//\\//\\//\\//\\ MACRO MAIN BODY //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\

void treeSiliconBig( int SiRunN, int MuRunN = 0, bool doAlignment = 0, int m_debug = 0 ) {
    
    Alignment_t Ali;
    int exitcode = 0;
    exitcode = loadAlignments( Ali, doAlignment, m_debug );
    if ( exitcode ) return;
    exitcode = readSievents( SiRunN, MuRunN, Ali, m_debug );
    if ( exitcode ) return;

}

//\\//\\//\\//\\ LOAD ALIGNMENTS //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\

int loadAlignments( Alignment_t &Ali, bool doAlignment, bool m_debug ) { 
    
    string line;
    ifstream myfile ("../utils/alignments.dat");
    vector<float> vAli;
    
    if ( myfile.is_open() ) {
        while( getline( myfile, line ) ){
            stringstream ss( line ); 
            string token;
            while ( getline( ss, token, ' ') ){
                if ( !TString( token ).IsFloat() ) {
                    cout << "ERROR. Non float entry in alignment file" << endl;
                    return 1;
                }
                if ( doAlignment ) {
                    vAli.push_back( TString( token ).Atof() );
                } else { 
                    vAli.push_back( 0 );
                }
            }   
        }
        
        if ( ( int )vAli.size() != 3*( nSiLayers+1 )+2) {
            cout << "ERROR in parsing alignment file" << endl;
            if ( m_debug ) cout << vAli.size() << endl;
            return 2;
        }
        
        for (int iAli=0; iAli<3*(nSiLayers+1); ++iAli) {
            int iLay=iAli/3;
            int xORy=iAli%3;
            Ali.x0SiAndMu[iLay][xORy] = vAli.at( iAli );
        }
        
        Ali.phi0Theta0Mu[0] = vAli.at( 3*( nSiLayers+1 )    );
        Ali.phi0Theta0Mu[1] = vAli.at( 3*( nSiLayers+1 ) +1 );

    }
    
    if ( m_debug ) printAlignments( Ali );
    
    return 0;
}

//\\//\\//\\//\\ APPLY ALIGNMENTS //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\

void applyAlignment( vector<Hit_t> &vhit, const Si_t &SiEvent, Alignment_t &Ali ) {

    vhit.clear();
    
    for ( Int_t iLayer=0; iLayer<nSiLayers; iLayer++ ) { 
        Int_t nxClusters = TMath::Min( SiEvent.nxClusters[iLayer], nMaxClusters);
        Int_t nyClusters = TMath::Min( SiEvent.nyClusters[iLayer], nMaxClusters);
        
        for ( Int_t iCluster=0; iCluster<nxClusters; iCluster++ ) {
            Hit_t hit;
            hit.subdet = detIDs[iLayer];
            hit.xh     = SiEvent.x[iLayer][iCluster] + Ali.x0SiAndMu[iLayer][0];
            hit.yh     = -999;
            hit.zh     =                               Ali.x0SiAndMu[iLayer][2];
            hit.itrack = 0;
            vhit.push_back( hit );
        }
        
        for (Int_t iCluster=0; iCluster<nyClusters; iCluster++) {
            Hit_t hit;
            hit.subdet = detIDs[iLayer];
            hit.xh     = -999;
            hit.yh     = SiEvent.y[iLayer][iCluster] + Ali.x0SiAndMu[iLayer][1];
            hit.zh     =                               Ali.x0SiAndMu[iLayer][2];
            hit.itrack = 0;
            vhit.push_back( hit );
        }
    }
}

//\\//\\//\\//\\ DEBUG: PRINT ALIGNMENTS //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//

void printAlignments( Alignment_t &Ali ) {
    
    cout << "**************** Beginning Print alignments..." << endl;
    cout << "Silicon detectors" << endl;
    for (Int_t iLayer=0; iLayer<nSiLayers+1; iLayer++) {
        cout << "x = " << Ali.x0SiAndMu[iLayer][0] << "\t y = "<< Ali.x0SiAndMu[iLayer][1] << "\t z ="<< Ali.x0SiAndMu[iLayer][2] << endl;
    }
    cout << endl << "Muon Chamber" << endl;
    cout <<  "phi0 = " << Ali.phi0Theta0Mu[0] << "\t theta0 = " << Ali.phi0Theta0Mu[1] << endl;
    
}

//\\//\\//\\//\\ READ EVENTS FROM SILICON ASCII FILE //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\

int readSievents( Int_t SiRunN, Int_t MuRunN, Alignment_t Ali, bool m_debug ) {
    
    Si_t SiEvent;
    stringstream si;
    si << SiRunN;
    string Si_RunN = si.str();
    
    stringstream sm;
    sm << MuRunN;
    string Mu_RunN = sm.str();
    
    TString outFile( TString( "/home/utente/LEMMA/data_Si/Run_") + TString( MuRunN!=0 ? Mu_RunN : Si_RunN ) + TString( "_Si_big.root" ) );
    
    TFile *f = new TFile( outFile,"recreate" );
    TTree *t = new TTree( "t", "a Tree with data from silicon, scintillator and calo detectors" );
    
    t->Branch( "eventID_Si",          &SiEvent.eventID_Si,          "eventID_Si/I"                                               );
    t->Branch( "eventID_mu",          &SiEvent.eventID_mu,          "eventID_mu/I"                                               );
    t->Branch( "eventID_mu_expanded", &SiEvent.eventID_mu_expanded, "eventID_mu_expanded/I"                                      );
    t->Branch( "digis1 ",              SiEvent.digis1,              TString::Format( "digis1[%i]/F", nDigis )                    );
    t->Branch( "x",                    SiEvent.x,                   TString::Format( "x[%i][%i]/F" , nSiLayers,nMaxClusters)     );
    t->Branch( "y",                    SiEvent.y,                   TString::Format( "y[%i][%i]/F" , nSiLayers,nMaxClusters)     );

    Int_t iev    = 0;
    Int_t iev_mu = 0;
    Int_t nhits  = 0;
    Int_t itrack = 0;
    vector <Hit_t> vHits;
    Hit_t   hit;
    Int_t   hitsubdet[nMaxHits];
    Int_t   hititrack[nMaxHits];
    Float_t hitxh    [nMaxHits];
    Float_t hityh    [nMaxHits];
    Float_t hitzh    [nMaxHits];
    
    TTree *t2 = new TTree( "LEMMA", "Event Builder LEMMA data output" );

    t2->Branch( "iev"   , &SiEvent.eventID_Si,          "iev/I"         );
    t2->Branch( "iev_mu", &SiEvent.eventID_mu_expanded, "iev_mu/I"      );
    t2->Branch( "nhits" , &nhits,                       "nhits/I"       );
    t2->Branch( "subdet",  hitsubdet,                   "subdet[100]/I" );
    t2->Branch( "xh",      hitxh,                       "xh[100]/F"     );
    t2->Branch( "yh",      hityh,                       "yh[100]/F"     );
    t2->Branch( "zh",      hitzh,                       "zh[100]/F"     );
    t2->Branch( "itrack",  hititrack,                   "itrack[100]/I" );
    
    string line;
    ifstream myfile( TString( "/home/utente/LEMMA/data_Si/run" ) + TString( Si_RunN ) + TString( "_big.dat" ) );
    //ifstream myfile( TString( "/home/utente/run" ) + TString( Si_RunN ) + TString( "_big_000001.dat" ) );
    
    if ( myfile.is_open() ) {
        bool isNewRecord=true;
        float vRecord[nTokensInRecord];
        Int_t nRecords  = 0;
        Int_t evCounter = 0;
        
        while(  myfile >> vRecord[nRecords] ) {
            
            ++nRecords;
            
            if ( nRecords == nTokensInRecord ) {
                
                for( Int_t i=0; i<nMaxHits; ++i ) {
                    hitsubdet[i] = -999;
                    hititrack[i] = -999;
                    hitxh    [i] = -999;
                    hityh    [i] = -999;
                    hitzh    [i] = -999;
                }
                
                if( evCounter%1000 == 0 ) cout << "Filling Event " << evCounter << endl;
                ++evCounter;
                fillEvent( vRecord, SiEvent );
                if ( m_debug ) printEvent( SiEvent );
                t->Fill();
                applyAlignment( vHits, SiEvent, Ali );
                
                nhits = TMath::Min( ( int )vHits.size(), nMaxHits );
                
                for (Int_t iHits=0; iHits<nhits; ++iHits) {
                    hitsubdet[iHits] = vHits.at( iHits ).subdet; // exactly the same structure in applyAlignment
                    hitxh    [iHits] = vHits.at( iHits ).xh    ;
                    hityh    [iHits] = vHits.at( iHits ).yh    ;
                    hitzh    [iHits] = vHits.at( iHits ).zh    ;
                    hititrack[iHits] = vHits.at( iHits ).itrack;
                    
                    if( m_debug ) printHit( vHits.at( iHits ) );
                }
                
                t2->Fill();
                
                nRecords=0;
            } 
        }
        
        myfile.close();
        
    } //end if (myfile.is_open())
    
    else {
        cout << "ERROR: Unable to open file" << endl;
        return 3;
    }
    
    t->Write();
    t2->Write();
    
    return 0;
    
}

//\\//\\//\\//\\ FILL EVENT //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\

void fillEvent( float *vRecord, Si_t &SiEvent ) {

    Ssiz_t it=0;
    
    SiEvent.eventID_Si           = vRecord[ it++ ];
    SiEvent.eventID_mu           = vRecord[ it++ ];
    SiEvent.eventID_mu_expanded  = vRecord[ it++ ];
        
    for ( Int_t idigi=0; idigi<nDigis; idigi++ ) {
        SiEvent.digis1[idigi]    = vRecord[ it++ ];
    }
    
    for ( Int_t iLayer=0; iLayer<nSiLayers; iLayer++ ) { 
        SiEvent.nxClusters[iLayer] = vRecord[ it++ ];
        for ( Int_t iCluster=0; iCluster<nMaxClusters; iCluster++ ) SiEvent.x[iLayer][iCluster] = vRecord[ it++ ];
        for ( Int_t iCluster=0; iCluster<nMaxClusters; iCluster++ ) SiEvent.xPulseHeight[iLayer][iCluster] = vRecord[ it++ ];
        SiEvent.nyClusters[iLayer] = vRecord[it++];
        for ( Int_t iCluster=0; iCluster<nMaxClusters; iCluster++ ) SiEvent.y[iLayer][iCluster] = vRecord[ it++ ];
        for ( Int_t iCluster=0; iCluster<nMaxClusters; iCluster++ ) SiEvent.yPulseHeight[iLayer][iCluster] = vRecord[ it++ ];
    }
    
    if ( it != nTokensInRecord ) {
        cout << "ERROR in parsing silicon data" << endl;
        return;
    }

}

//\\//\\//\\//\\ DEBUG: PRINT EVENT //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\

void printEvent( Si_t &Ev ) {
    
    cout << "**************** Beginning print event..." << endl;
    cout << "Event ID Si = " << Ev.eventID_Si << endl;
    cout << "Event ID Mu = " << Ev.eventID_mu << endl;
    cout << "Event ID Mu expanded = " << Ev.eventID_mu_expanded << endl;
    
    for ( Int_t idigi=0; idigi<nDigis; idigi++ ) {
        cout << Ev.digis1[idigi] << "\t";
    }
    
    for ( Int_t iLayer=0; iLayer<nSiLayers; iLayer++ ) {
        cout << endl << "Layer = " << iLayer+1 << endl;
        cout << "Number of clusters x view = " << Ev.nxClusters[iLayer] << endl;
        for ( Int_t iCluster=0; iCluster<nMaxClusters; iCluster++ ) {
            cout << "xh = " << Ev.x[iLayer][iCluster] << "\t";
        }
        cout << endl;
        for ( Int_t iCluster=0; iCluster<nMaxClusters; iCluster++ ) {
            cout << "ph = " << Ev.xPulseHeight[iLayer][iCluster] << "\t";
        }
        cout << endl;
        cout << "Number of clusters y view = " << Ev.nyClusters[iLayer] << endl;
        for ( Int_t iCluster=0; iCluster<nMaxClusters; iCluster++ ) {
            cout << "yh = " << Ev.y[iLayer][iCluster] << "\t";
        }
        cout << endl;
        for ( Int_t iCluster=0; iCluster<nMaxClusters; iCluster++ ) {
            cout << "ph = " << Ev.yPulseHeight[iLayer][iCluster] << "\t";
        }
        cout << endl;
    }
}

//\\//\\//\\//\\ DEBUG: PRINT HIT //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\

void printHit( Hit_t &hit ) {
    cout << "**************** Beginning Print Hit..." << endl;
    cout << "subdet = " << hit.subdet << " xh = " <<  hit.xh << " yh = " <<  hit.yh << " zh = " <<  hit.zh << endl;
}
