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
    Int_t    digis1[nDigis]; //TODO per ora questa parte sta qui in Si_t poi bisognerà o cambiare questo nome o costruire un'altra struttura
    Int_t    nxClusters[nSiLayers]; 
    Float_t  x[nSiLayers][nMaxClusters];
    Float_t  xPulseHeight[nSiLayers][nMaxClusters];
    Int_t    nyClusters[nSiLayers];
    Float_t  y[nSiLayers][nMaxClusters];
    Float_t  yPulseHeight[nSiLayers][nMaxClusters];
    
} Si_t; 

typedef struct {
    //TODO aggiungere informazioni anche dei calorimetri e degli scintillatori
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
int  readSievents    ( Int_t runN, Alignment_t Ali, bool m_debug );
void fillEvent       ( vector<float> *vRecord, Si_t &SiEvent );
void printEvent      ( Si_t &Ev );
void printHit        ( Hit_t &hit );

//\\//\\//\\//\\ MACRO MAIN BODY //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\

void treeSiliconBig_bad( int runN, bool doAlignment = 0, int m_debug = 0 ) {
    
    Alignment_t Ali;
    int exitcode = 0;
    exitcode = loadAlignments( Ali, doAlignment, m_debug );
    if ( exitcode ) return;
    exitcode = readSievents( runN, Ali, m_debug );
    if ( exitcode ) return;

}

//\\//\\//\\//\\ LOAD ALIGNMENTS //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\

int loadAlignments( Alignment_t &Ali, bool doAlignment, bool m_debug ) { //TODO modificare con le informazioni dei calorimetri e degli scintillatori
    
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
    
    cout << "vhit.size() = " << vhit.size() << endl;
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
    
    //TODO modificare con le informazioni dei calorimetri e degli scintillatori
    cout << "**************** Beginning Print alignments..." << endl;
    cout << "Silicon detectors" << endl;
    for (Int_t iLayer=0; iLayer<nSiLayers+1; iLayer++) {
        cout << "x = " << Ali.x0SiAndMu[iLayer][0] << "\t y = "<< Ali.x0SiAndMu[iLayer][1] << "\t z ="<< Ali.x0SiAndMu[iLayer][2] << endl;
    }
    cout << endl << "Muon Chamber" << endl;
    cout <<  "phi0 = " << Ali.phi0Theta0Mu[0] << "\t theta0 = " << Ali.phi0Theta0Mu[1] << endl;
    
}

//\\//\\//\\//\\ READ EVENTS FROM SILICON ASCII FILE //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\

int readSievents( Int_t runN, Alignment_t Ali, bool m_debug ) {
    
    Si_t SiEvent;
    stringstream ss;
    ss << runN;
    string srunN = ss.str();
    
    TString outFile( TString( "/home/utente/LEMMA/data/Run_") + TString( srunN ) + TString( "_SiEcal.root" ) ); 
    
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
    Int_t   hitsubdet[100]={-999};
    Int_t   hititrack[100]={-999};
    Float_t hitxh    [100]={-999};
    Float_t hityh    [100]={-999};
    Float_t hitzh    [100]={-999};
    
    TTree *t2 = new TTree( "LEMMA", "Event Builder LEMMA data output" );
    //TODO la struttura va modificata quando si inseriranno le info dei calorimetri
    t2->Branch( "iev"   , &SiEvent.eventID_Si,          "iev/I"         );
    t2->Branch( "iev_mu", &SiEvent.eventID_mu_expanded, "iev_mu/I"      );
    t2->Branch( "nhits" , &nhits,                       "nhits/I"       );
    t2->Branch( "subdet",  hitsubdet,                   "subdet[100]/I" );
    t2->Branch( "xh",      hitxh,                       "xh[100]/F"     );
    t2->Branch( "yh",      hityh,                       "yh[100]/F"     );
    t2->Branch( "zh",      hitzh,                       "zh[100]/F"     );
    t2->Branch( "itrack",  hititrack,                   "itrack[100]/I" );
    
    string line;
//     ifstream myfile( TString( "/home/utente/LEMMA/data/run" ) + TString( srunN ) + TString( "_big.dat" ) );
    ifstream myfile( TString( "/home/utente/run" ) + TString( srunN ) + TString( "_big_000001.dat" ) );
    
    if ( myfile.is_open() ) {
        bool isNewRecord=true;
        vector<float> vRecord;
        Int_t nRecords = 0;
        
        // --- loop on lines 
        while ( getline ( myfile,line ) ) {
            if ( m_debug ) cout << line << endl;
            
            TString sline(line);
            TString token;
            Ssiz_t  from = 0;
            Ssiz_t  iLayer = 0;
            
            TObjArray* slineArray = sline.Tokenize(" ");
            Int_t ntokens = slineArray->GetEntries();
            isNewRecord = (ntokens > 30 ? 1 : 0);
            
            if ( isNewRecord && vRecord.size()>0 ) {
                if ( vRecord.size() != nTokensInRecord ) { 
                    cout << "ERROR in reading Si data. Wrong # of tokens per event" << endl;
                    return 1;
                }
                
                if ( ( nRecords++ )%1000==0 ) cout << "Filling vRecord " << nRecords << endl; 
                
                fillEvent( &vRecord, SiEvent );
                
                if ( m_debug ) printEvent( SiEvent );
                
                cout << "vRecord.size() = " << vRecord.size() << endl;
                vRecord.clear();
                
                isNewRecord=false;
                
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
                
            } // end if (isNewRecord && vRecord.size()>0)
            
            for ( Int_t it=0; it<ntokens; ++it ) { // Store all tokes of each line
                TObjString* objString = ( TObjString* ) ( slineArray->At( it ) );
                TString stoken = objString->String();
                
                if ( !stoken.IsFloat() ) {
                    if ( m_debug ) cout << stoken.Data() << endl;
                    cout << "ERROR in reading. Non float entry in record" << endl;
                    return 2;
                }
                
                vRecord.push_back( stoken.Atof() + epsilon ); 
            } 
            
        } // end loop on lines
        
        myfile.close();
        
        //--- last record added at the end of the parsing
        fillEvent( &vRecord, SiEvent ); 
        
        t->Fill();
        
        applyAlignment( vHits, SiEvent, Ali );
        nhits = TMath::Min( ( int )vHits.size(), nMaxHits );
        for (Int_t iHits=0; iHits < nhits; ++iHits) {
            hitsubdet[iHits] = vHits.at( iHits ).subdet;
            hitxh    [iHits] = vHits.at( iHits ).xh    ;
            hityh    [iHits] = vHits.at( iHits ).yh    ;
            hitzh    [iHits] = vHits.at( iHits ).zh    ;
            hititrack[iHits] = vHits.at( iHits ).itrack;
        }
        
        t2->Fill();
        
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

void fillEvent( vector<float> *vRecord, Si_t &SiEvent ) {

    Ssiz_t it=0;
    
    SiEvent.eventID_Si           = vRecord->at( it++ );
    SiEvent.eventID_mu           = vRecord->at( it++ );
    SiEvent.eventID_mu_expanded  = vRecord->at( it++ );
        
    for ( Int_t idigi=0; idigi<nDigis; idigi++ ) { //TODO in this part you will store calorimeters and scintillator detectors info
        SiEvent.digis1[idigi]    = vRecord->at( it++ );
    }
    
    for ( Int_t iLayer=0; iLayer<nSiLayers; iLayer++ ) { 
        SiEvent.nxClusters[iLayer] = vRecord->at( it++ );
        for ( Int_t iCluster=0; iCluster<nMaxClusters; iCluster++ ) SiEvent.x[iLayer][iCluster] = vRecord->at( it++ );
        for ( Int_t iCluster=0; iCluster<nMaxClusters; iCluster++ ) SiEvent.xPulseHeight[iLayer][iCluster] = vRecord->at( it++ );
        SiEvent.nyClusters[iLayer] = vRecord->at(it++);
        for ( Int_t iCluster=0; iCluster<nMaxClusters; iCluster++ ) SiEvent.y[iLayer][iCluster] = vRecord->at( it++ );
        for ( Int_t iCluster=0; iCluster<nMaxClusters; iCluster++ ) SiEvent.yPulseHeight[iLayer][iCluster] = vRecord->at( it++ );
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
        cout << Ev.digis1[idigi] << "\t"; //TODO modify when introducing ecal e scintillators
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
