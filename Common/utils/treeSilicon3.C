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


// This example illustrates how to make a Tree from variables or arrays
// in a C struct. Use of C structs is strongly discouraged and one should
// use classes instead. However support for C structs is important for 
// legacy applications written in C or Fortran.
//    see tree2a.C for the same example using a class instead of a C-struct.
//
// In this example, we are mapping a C struct to one of the Geant3
// common blocks /gctrak/. In the real life, this common will be filled
// by Geant3 at each step and only the Tree Fill function should be called.
// The example emulates the Geant3 step routines.
//
// to run the example, do:
// .x tree2.C   to execute with the CINT interpreter
// .x tree2.C++ to execute with native compiler
//
//  Author: Mi nonno
// // *******  Data structure ******
// // A new record starts with 55 tokens (in the code the check is implemented with a >30)
// // First token is the eventID from the Si detectors. Following 4 are x1, y1, x2, y2; one cluster per view/detector.
// // The next entry is the eventID from the mu
// // The next is a triggerInfo
// // Then there are 16*3 digis I can store them but I do not know what to do with them
// // Then there are 2*(nLayers-2) lines. Each line is a view/detector. 6 tokens per line. nClusters, 5*localXorY (for the first 5 clusters only)
// // 000000116  0.66560E+00  0.43680E+01 -0.30000E+04 -0.40000E+04 000001 0000 00000 00000 00000 00000 00000 00000 00000 00000 00000 00000 00000 00000 00000 00000 00000 00000 00000 00000 00000 00000 00000 00000 00000 00000 00000 00000 00000 00000 00000 00000 00000 00000 000 000 000 000 000 000 000 000 000 000 000 000 000 000 000 000
// // 00 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02
// // 00 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02
// // 00 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02
// // 00 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02
// // 00 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02
// // 00 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02
// // 00 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02
// // 00 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02
// // 00 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02
// // 00 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02
// // 00 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02
// // 00 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02 -0.99000E+02



const Int_t nSiLayers = 8; // Si layers
const Int_t nDigis  = 48;
const Int_t nMaxClusters = 5;
const Int_t nTokensInRecord = 1+2*2+1+1+nDigis+2*(nSiLayers-2)*(1+nMaxClusters);
const float epsilon = 1E-12;
float alignments[nSiLayers+1][2];
const Int_t detIDs[nSiLayers+1]={10,20,30,40,50,51,55,56,70};
const Int_t nMaxHits = 100;


typedef struct { 
    Int_t    eventID_Si; 
    Float_t  x[nSiLayers][nMaxClusters]; 
    Float_t  y[nSiLayers][nMaxClusters];  
    Int_t    triggerInfo;
    Int_t    eventID_mu; 
    Int_t    digis1[nDigis];
    //   Int_t    digis2[nDigis];
    //   Int_t    digis3[nDigis];
    Int_t    nxClusters[nSiLayers]; 
    Int_t    nyClusters[nSiLayers]; 
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

void printEvent(Si_t &Ev) {
    std::cout << Ev.eventID_Si << std::endl; 
    for (Int_t iLayer=0; iLayer<nSiLayers; iLayer++) {
        std::cout << "Ev.x         ["<< iLayer << "][0]       " << Ev.x[iLayer][0]          << std::endl;
        std::cout << "Ev.y         ["<< iLayer << "][0]       " << Ev.y[iLayer][0]           << std::endl;
        std::cout << "Ev.nxClusters["<< iLayer << "] " << Ev.nxClusters[iLayer]    << std::endl;
        std::cout << "Ev.nyClusters["<< iLayer << "] " << Ev.nyClusters[iLayer]    << std::endl;
    }		       	          		        
    std::cout   << "Ev.triggerInfo	" << Ev.triggerInfo	           << std::endl;
    std::cout   << "Ev.eventID_mu         " << Ev.eventID_mu              << std::endl;
    
}

void printHit(Hit_t &hit) {
    std::cout << "\033[0;2m" << hit.subdet << " " <<  hit.xh << " " <<  hit.yh << " " <<  hit.zh << "\033[0m" << std::endl;
}

void printAlignments(Alignment_t &Ali) {
    cout << "xh \t yh \t zh" << endl;
    for (Int_t iLayer=0; iLayer<nSiLayers+1; iLayer++) {
        cout << Ali.x0SiAndMu[iLayer][0] << "\t"<< Ali.x0SiAndMu[iLayer][1] << "\t"<< Ali.x0SiAndMu[iLayer][2] << endl;
    }
    cout << "phi0mu \t theta0mu" << endl;
    cout <<  Ali.phi0Theta0Mu[0] << "\t" << Ali.phi0Theta0Mu[1] << endl;
}

int loadAlignments(Alignment_t &Ali, bool doIt, bool debug) {
    string line;
    ifstream myfile ("../utils/alignments.dat");
    vector<float> vAli;
    if (myfile.is_open())
    {
        while( getline(myfile, line) )
        {
            std::stringstream ss(line); 
            std::string token;
            while (std::getline(ss, token, ' ')){
                if (!TString(token).IsFloat()) {
                    cout << "\033[1;31mERROR. Non float entry in alignment file\033[0m" << endl;
                    return 1;
                }
                if (doIt) vAli.push_back(TString(token).Atof());
                else vAli.push_back(0);
            }   
        }
        if ((int)vAli.size() != 3*(nSiLayers+1)+2) {
            cout << "\033[1;31mERROR in parsing alignment file\033[0m" << endl;
            if (debug) cout << vAli.size() << endl;
            return 2;
        }
        for (int iAli=0; iAli<3*(nSiLayers+1); ++iAli) {
            int iLay=iAli/3;
            int xORy=iAli%3;
            Ali.x0SiAndMu[iLay][xORy] = vAli.at(iAli);
        }
        Ali.phi0Theta0Mu[0] = vAli.at(3*(nSiLayers+1)  );
        Ali.phi0Theta0Mu[1] = vAli.at(3*(nSiLayers+1)+1);
        
        //     while ( getline (myfile,line) )
        //     {
        //       cout << line << endl;
        //     }
    }
    
    printAlignments(Ali);
    
    return 0;
}

int applyAlignment(vector<Hit_t> & vhit, const Si_t &SiEvent, Alignment_t &Ali) {
    vhit.clear();
    for (Int_t iLayer=0; iLayer<nSiLayers; iLayer++) { // align Silicon detectors
        Int_t nxClusters = TMath::Min(SiEvent.nxClusters[iLayer],nMaxClusters);
        Int_t nyClusters = TMath::Min(SiEvent.nyClusters[iLayer],nMaxClusters);
        for (Int_t iCluster=0; iCluster<nxClusters; iCluster++) {
            Hit_t hit;
            hit.subdet = detIDs[iLayer];
            hit.xh     = SiEvent.x[iLayer][iCluster] + Ali.x0SiAndMu[iLayer][0];
            hit.yh     = -999;
            hit.zh     =                               Ali.x0SiAndMu[iLayer][2];
            hit.itrack = 0; //TO DO. Still not clear how to use it.
            vhit.push_back(hit);
        }
        for (Int_t iCluster=0; iCluster<nyClusters; iCluster++) {
            Hit_t hit;
            hit.subdet = detIDs[iLayer];
            hit.xh     = -999;
            hit.yh     = SiEvent.y[iLayer][iCluster] + Ali.x0SiAndMu[iLayer][1];
            hit.zh     =                               Ali.x0SiAndMu[iLayer][2];
            hit.itrack = 0; //TODO Still not clear how to use it.
            vhit.push_back(hit);
        }
    }
    return 0;
}

int fillEvent( vector<float> *vRecord, Si_t &SiEvent) {
    
    if (0) {
        std::cout << "vRecordSize = " << vRecord->size() << std::endl;
        for (Int_t irec=0; irec<(int)vRecord->size(); ++irec) {
            std::cout << vRecord->at(irec) << " " ;
            if (irec==54||irec==60||irec==66||irec==72||irec==78||irec==84||irec==90||irec==96||irec==102||irec==108||irec==114||irec==120)
                std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    
    Ssiz_t it=0;
    SiEvent.eventID_Si  = vRecord->at(it++);
    SiEvent.x[0][0]     = vRecord->at(it++);
    SiEvent.y[0][0]     = vRecord->at(it++);
    SiEvent.x[1][0]     = vRecord->at(it++);
    SiEvent.y[1][0]     = vRecord->at(it++);
     
        
    for (Int_t iLayer=0; iLayer<2; iLayer++) {
        for (Int_t iCluster=1; iCluster<nMaxClusters; iCluster++) {
            SiEvent.x[iLayer][iCluster] = 0;
            SiEvent.y[iLayer][iCluster] = 0;
        }
    }
    
    SiEvent.nxClusters[0] = (int)(SiEvent.x[0][0]>-99); // subdet 10 & 20 -> single cluster mode
    SiEvent.nyClusters[0] = (int)(SiEvent.y[0][0]>-99);
    SiEvent.nxClusters[1] = (int)(SiEvent.x[1][0]>-99);
    SiEvent.nyClusters[1] = (int)(SiEvent.y[1][0]>-99);
    
    SiEvent.eventID_mu  = vRecord->at(it++);
    SiEvent.triggerInfo = vRecord->at(it++);
        
    for (Int_t idigi=0; idigi<nDigis; idigi++) {
        SiEvent.digis1[idigi] =  vRecord->at(it++);
    }
    
    for (Int_t iLayer=2; iLayer<nSiLayers; iLayer++) { 
        SiEvent.nxClusters[iLayer] = vRecord->at(it++);
        for (Int_t iCluster=0; iCluster<nMaxClusters; iCluster++) SiEvent.x[iLayer][iCluster] = vRecord->at(it++);
        SiEvent.nyClusters[iLayer] = vRecord->at(it++);
        for (Int_t iCluster=0; iCluster<nMaxClusters; iCluster++) SiEvent.y[iLayer][iCluster] = vRecord->at(it++);
    }
    
    if (it != nTokensInRecord) {
        cout << "\033[1;31mERROR in parsing silicon data\033[0m" << endl;
        return 1;
    }
    return 0;
}


int readSievents(Int_t runN, Alignment_t Ali, bool debug) {
    
    Si_t SiEvent;
    stringstream ss;
    ss << runN;
    string srunN = ss.str();
    //  char *intStr = itoa(runN);
    //string srunN = string(intStr);
    //std::string srunN = std::to_string(runN);
    TString outFile(TString("/home/utente/LEMMA/data_Si/Run_")+TString(srunN)+TString("_Si.root"));
    
    TFile f(outFile,"recreate");
    TTree t("t","a Tree with data from TB Silicon");
    t.Branch("eventID_Si"    ,&SiEvent.eventID_Si    ,"eventID_Si/I"                                                  );
    t.Branch("x"             , SiEvent.x             ,TString::Format("x[%i][%i]/F"     , nSiLayers,nMaxClusters)     );
    t.Branch("y"             , SiEvent.y             ,TString::Format("y[%i][%i]/F"     , nSiLayers,nMaxClusters)     );
    t.Branch("triggerInfo"   ,&SiEvent.triggerInfo   ,"triggerInfo/I"                                                 );
    t.Branch("eventID_mu"    ,&SiEvent.eventID_mu    ,"eventID_mu/I"                                                  );
    t.Branch("digis1 "       , SiEvent.digis1        ,TString::Format("digis1[%i]/F", nDigis )                        );
    
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
    
    Int_t ID_mu_expanded;
    
    //   TFile f2("treeSiliconAligned.root","recreate");
    TTree t2("LEMMA","Event Builder LEMMA data output"  );
    t2.Branch("iev"   , &SiEvent.eventID_Si, "iev/I"    );
    t2.Branch("iev_mu" ,&ID_mu_expanded, "iev_mu/I"     );
    t2.Branch("nhits" , &nhits     , "nhits/I"          );
    t2.Branch("subdet",  hitsubdet , "subdet[100]/I"    );
    t2.Branch("xh",      hitxh     , "xh[100]/F"        );
    t2.Branch("yh",      hityh     , "yh[100]/F"        );
    t2.Branch("zh",      hitzh     , "zh[100]/F"        );
    t2.Branch("itrack",  hititrack , "itrack[100]/I"    );
    
    string line;
//     ifstream myfile (TString("/home/utente/LEMMA/data/run")+TString(srunN)+TString("_multi.dat"));
    ifstream myfile (TString("/home/utente/LEMMA/data_Si/run200568_multi.dat"));

    if (myfile.is_open())
    {
        bool isNewRecord=true;
        //    TObjArray objRecord;
        vector<float> vRecord;
        Int_t nRecords = 0;
        
        Int_t previous_ID = -99999;
        Int_t count =0;
        
        // --- loop on lines 
        while ( getline (myfile,line) )
        {
            if (debug) cout << line << endl;
            
            TString sline(line);
            TString token;
            Ssiz_t  from = 0;
            Ssiz_t  iLayer = 0;
            
            TObjArray* slineArray = sline.Tokenize(" ");
            Int_t ntokens = slineArray->GetEntries();
            isNewRecord = (ntokens > 30 ? 1 : 0);
            
            if (isNewRecord && vRecord.size()>0) {
                
                if (vRecord.size() != nTokensInRecord) { 
                    cout << "\033[1;31mERROR in reading Si data. Wrong # of tokens per event \033[0m" << endl;
                    return 1;
                }
                
                if ((nRecords++)%100==0) cout << "Filling vRecord " << nRecords << endl; 
                
                for( Int_t i=0; i<nMaxHits; ++i ) { // ALTEA added 22-08-2017
                    hitsubdet[i] = -999;
                    hititrack[i] = -999;
                    hitxh    [i] = -999;
                    hityh    [i] = -999;
                    hitzh    [i] = -999;
                }
                
                fillEvent(&vRecord, SiEvent);
                
                if (previous_ID > SiEvent.eventID_mu) count++;
                
                ID_mu_expanded = SiEvent.eventID_mu  + count*pow(2,16);
                
                previous_ID = SiEvent.eventID_mu;
                
                // 	    printEvent(SiEvent);
                vRecord.clear();
                isNewRecord=0;
                t.Fill();
                
                applyAlignment(vHits, SiEvent, Ali);
                nhits = TMath::Min((int)vHits.size(),nMaxHits);
                
                for (Int_t iHits=0; iHits < nhits; ++iHits) {
                    hitsubdet[iHits] = vHits.at(iHits).subdet;
                    hitxh    [iHits] = vHits.at(iHits).xh    ;
                    hityh    [iHits] = vHits.at(iHits).yh    ;
                    hitzh    [iHits] = vHits.at(iHits).zh    ;
                    hititrack[iHits] = vHits.at(iHits).itrack;
                    // 	      printHit(vHits.at(iHits));
                    // 	      std::cout << "\033[0;3m" << hitsubdet[iHits] << " " <<  hitxh    [iHits] << " " <<  hityh    [iHits] << " " <<  hitzh    [iHits] << "\033[0m" << std::endl;
                }
                
                t2.Fill();
                
                
                // 	    for (Int_t iHits=0; iHits<nhits; ++iHits) {
                // 	      cout << nhits << endl;
                // 	      hit = vHits.at(iHits);
                // 	      cout << SiEvent.eventID_Si << " " << nhits << " " << hit.subdet << " " << hit.xh << " " << hit.yh << " " << hit.zh <<  " " << hit.itrack << endl;
                // 	      t2.Fill();
                // 	    }
                
            } // end if (isNewRecord && vRecord.size()>0)
            
            for (Int_t it=0; it<ntokens; ++it) { // Store all tokes of each line
                TObjString* objString = (TObjString*) (slineArray->At(it));
                TString stoken = objString->String();
                // 	cout << stoken.Data() << endl;
                
                if (!stoken.IsFloat()) {
                    if (debug) cout << stoken.Data() << endl;
                    cout << "\033[1;31mERROR in reading. Non float entry in record\033[0m" << endl;
                    return 2;
                }
                
                vRecord.push_back(stoken.Atof()+epsilon); 
            } 
            
            //      while (sline.Tokenize(token, from, " ")) {
            //	cout << token.Data() << endl;
            //      }
            
        } // end loop on lines
        
        myfile.close();
        fillEvent(&vRecord, SiEvent); //last record added at the end of the parsing
        t.Fill();
        
        applyAlignment(vHits, SiEvent, Ali);
        nhits = TMath::Min((int)vHits.size(),nMaxHits);
        for (Int_t iHits=0; iHits < nhits; ++iHits) {
            hitsubdet[iHits] = vHits.at(iHits).subdet;
            hitxh    [iHits] = vHits.at(iHits).xh    ;
            hityh    [iHits] = vHits.at(iHits).yh    ;
            hitzh    [iHits] = vHits.at(iHits).zh    ;
            hititrack[iHits] = vHits.at(iHits).itrack;
        }
        t2.Fill();
        
    } //end if (myfile.is_open())
    
    else {
        cout << "Unable to open file" << endl;
        return 3;
    }
    
    t.Write();
    t2.Write();
    
    return 0;
    
}



//void tree2r()
//{
//   //read the Tree generated by tree2w and fill one histogram
//   //we are only interested by the destep branch.
//     
//   //note that we use "new" to create the TFile and TTree objects !
//   //because we want to keep these objects alive when we leave 
//   //this function.
//   TFile *f = new TFile("tree2.root");
//   TTree *t2 = (TTree*)f->Get("t2");
//   static Float_t destep;
//   TBranch *b_destep = t2->GetBranch("destep");
//   b_destep->SetAddress(&destep);
//   
//   //create one histogram
//   TH1F *hdestep   = new TH1F("hdestep","destep in Mev",100,1e-5,3e-5);
//   
//   //read only the destep branch for all entries
//   Long64_t nentries = t2->GetEntries();
//   for (Long64_t i=0;i<nentries;i++) {
//      b_destep->GetEntry(i); 
//      hdestep->Fill(destep);
//   }
//  
//   //we do not close the file. 
//   //We want to keep the generated histograms
//   //We fill a 3-d scatter plot with the particle step coordinates
//   TCanvas *c1 = new TCanvas("c1","c1",600,800);
//   c1->SetFillColor(42);
//   c1->Divide(1,2);
//   c1->cd(1);
//   hdestep->SetFillColor(45);
//   hdestep->Fit("gaus");
//   c1->cd(2);
//   gPad->SetFillColor(37);
//   t2->SetMarkerColor(kRed);
//   t2->Draw("vect[0]:vect[1]:vect[2]");
//   if (gROOT->IsBatch()) return;
//   
//   // invoke the x3d viewer
//   gPad->GetViewer3D("x3d");
//}   
//
void treeSilicon3(int runN, bool doAlignment = 0, int debug = 0) {
    
    Alignment_t Ali;
    int exitcode = 0;
    exitcode = loadAlignments(Ali,doAlignment,debug);
    if (exitcode) return;
    exitcode = readSievents(runN,Ali,debug);
    if (exitcode) return;
    //   tree2r();
}
