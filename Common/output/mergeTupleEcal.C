/********************************************************
 *                  mergeTupleEcal.C                    *
 ********************************************************
 * Modified mergeTuple.C by A. Lorenzon
 * 
 * Merge root file pos e root file neg produced after
 * running the new version of Event Builder.
 * 
 * FIRST  ARGUMENT: neg    file
 * SECOND ARGUMENT: pos    file
 * THIRD  ARGUMENT: merged file
 * 
 */

#include<TTree.h>
#include<TFile.h>
#include<iostream>

using namespace std;

void mergeTupleEcal(TString file_neg, TString file_pos, TString file_new){
    
    // --- leaf types of neg tree
    Int_t           siiev1;
    Int_t           iev1;
    Int_t           nhits1;
    Int_t           subdet1[100];   
    Float_t         xh1[100];
    Float_t         xh1_PulseHeight[100];
    Float_t         yh1[100];
    Float_t         yh1_PulseHeight[100];
    Float_t         zh1[100];   
    Int_t           itrack1[100];
    Int_t           S21_PulseHeight;
    Int_t           S21_TimeSpectrum;
    Int_t           LemmaTrigger1_PulseHeight;
    Int_t           LemmaTrigger1_TimeSpectrum;
    Int_t           Deva1_PulseHeight[6];
    Int_t           Deva1_TimeSpectrum[6];
    Int_t           Cherenkov1_PulseHeight[4];
    Int_t           Cherenkov1_TimeSpectrum[4];
    Int_t           S31_PulseHeight;
    Int_t           S31_TimeSpectrum;
    Int_t           PbGlass1_PulseHeight;
    Int_t           PbGlass1_TimeSpectrum;
    
    // --- leaf typer of pos tree
    Int_t           siiev2;
    Int_t           iev2;
    Int_t           nhits2;
    Int_t           subdet2[100];   
    Float_t         xh2[100];
    Float_t         xh2_PulseHeight[100];
    Float_t         yh2[100];  
    Float_t         yh2_PulseHeight[100];
    Float_t         zh2[100];   
    Int_t           itrack2[100];
    Int_t           S22_PulseHeight;
    Int_t           S22_TimeSpectrum;
    Int_t           LemmaTrigger2_PulseHeight;
    Int_t           LemmaTrigger2_TimeSpectrum;
    Int_t           Deva2_PulseHeight[6];
    Int_t           Deva2_TimeSpectrum[6];
    Int_t           Cherenkov2_PulseHeight[4];
    Int_t           Cherenkov2_TimeSpectrum[4];
    Int_t           S32_PulseHeight;
    Int_t           S32_TimeSpectrum;
    Int_t           PbGlass2_PulseHeight;
    Int_t           PbGlass2_TimeSpectrum;
    
    // --- leaf types of merged tree
    Int_t           siiev_new;
    Int_t           iev_new;
    Int_t           nhits_new;
    Int_t           subdet_new[100];
    Float_t         xh_new[100]; 
    Float_t         xh_new_PulseHeight[100];
    Float_t         yh_new[100]; 
    Float_t         yh_new_PulseHeight[100];
    Float_t         zh_new[100];   
    Int_t           itrack_new[100];
    Int_t           S2_new_PulseHeight;
    Int_t           S2_new_TimeSpectrum;
    Int_t           LemmaTrigger_new_PulseHeight;
    Int_t           LemmaTrigger_new_TimeSpectrum;
    Int_t           Deva_new_PulseHeight[6];
    Int_t           Deva_new_TimeSpectrum[6];
    Int_t           Cherenkov_new_PulseHeight[4];
    Int_t           Cherenkov_new_TimeSpectrum[4];
    Int_t           S3_new_PulseHeight;
    Int_t           S3_new_TimeSpectrum;
    Int_t           PbGlass_new_PulseHeight;
    Int_t           PbGlass_new_TimeSpectrum;
    
    // --- Open neg root file
    TFile *f1 = new TFile( file_neg );
    TTree *t1 = ( TTree* ) f1->Get( "LEMMA" );
    
    t1->SetBranchAddress( "siiev", &siiev1);
    t1->SetBranchAddress( "iev", &iev1);
    t1->SetBranchAddress( "nhits", &nhits1);
    t1->SetBranchAddress( "subdet", subdet1);
    t1->SetBranchAddress( "xh", xh1);
    t1->SetBranchAddress( "xh_PulseHeight",  xh1_PulseHeight );
    t1->SetBranchAddress( "yh", yh1);
    t1->SetBranchAddress( "yh_PulseHeight",  yh1_PulseHeight        );
    t1->SetBranchAddress( "zh", zh1);
    t1->SetBranchAddress( "itrack", itrack1);
    t1->SetBranchAddress( "S2_PulseHeight",            &S21_PulseHeight            );
    t1->SetBranchAddress( "LemmaTrigger_PulseHeight",  &LemmaTrigger1_PulseHeight  );   
    t1->SetBranchAddress( "Deva_PulseHeight",           Deva1_PulseHeight          );
    t1->SetBranchAddress( "Cherenkov_PulseHeight",      Cherenkov1_PulseHeight     );
    t1->SetBranchAddress( "S3_PulseHeight",            &S31_PulseHeight            );
    t1->SetBranchAddress( "PbGlass_PulseHeight",       &PbGlass1_PulseHeight       );
    t1->SetBranchAddress( "S2_TimeSpectrum",           &S21_TimeSpectrum           );
    t1->SetBranchAddress( "LemmaTrigger_TimeSpectrum", &LemmaTrigger1_TimeSpectrum );
    t1->SetBranchAddress( "Deva_TimeSpectrum",          Deva1_TimeSpectrum         );
    t1->SetBranchAddress( "Cherenkov_TimeSpectrum",     Cherenkov1_TimeSpectrum     );
    t1->SetBranchAddress( "S3_TimeSpectrum",           &S31_TimeSpectrum           );
    t1->SetBranchAddress( "PbGlass_TimeSpectrum",      &PbGlass1_TimeSpectrum      );
    
    
    // --- Open pos root file
    TFile *f2 = new TFile( file_pos );
    TTree *t2 = ( TTree* ) f2->Get( "LEMMA" );
    
    t2->SetBranchAddress( "siiev",                     &siiev2                     );
    t2->SetBranchAddress( "iev",                       &iev2                       );
    t2->SetBranchAddress( "nhits",                     &nhits2                     );
    t2->SetBranchAddress( "subdet",                     subdet2                    );
    t2->SetBranchAddress( "xh",                         xh2                        );
    t2->SetBranchAddress( "xh_PulseHeight",             xh2_PulseHeight            );
    t2->SetBranchAddress( "yh",                         yh2                        );
    t2->SetBranchAddress( "yh_PulseHeight",             yh2_PulseHeight            );
    t2->SetBranchAddress( "zh",                         zh2                        );
    t2->SetBranchAddress( "itrack",                     itrack2                    );
    t2->SetBranchAddress( "S2_PulseHeight",            &S22_PulseHeight            );
    t2->SetBranchAddress( "LemmaTrigger_PulseHeight",  &LemmaTrigger2_PulseHeight  );   
    t2->SetBranchAddress( "Deva_PulseHeight",           Deva2_PulseHeight          );
    t2->SetBranchAddress( "Cherenkov_PulseHeight",      Cherenkov2_PulseHeight     );
    t2->SetBranchAddress( "S3_PulseHeight",            &S32_PulseHeight            );
    t2->SetBranchAddress( "PbGlass_PulseHeight",       &PbGlass2_PulseHeight       );
    t2->SetBranchAddress( "S2_TimeSpectrum",           &S22_TimeSpectrum           );
    t2->SetBranchAddress( "LemmaTrigger_TimeSpectrum", &LemmaTrigger2_TimeSpectrum );
    t2->SetBranchAddress( "Deva_TimeSpectrum",          Deva2_TimeSpectrum         );
    t2->SetBranchAddress( "Cherenkov_TimeSpectrum",     Cherenkov2_TimeSpectrum    );
    t2->SetBranchAddress( "S3_TimeSpectrum",           &S32_TimeSpectrum           );
    t2->SetBranchAddress( "PbGlass_TimeSpectrum",      &PbGlass2_TimeSpectrum      );
    
    // --- Open new root file
    TFile *f_new = new TFile(file_new,"RECREATE");
    TTree *t_new = new TTree("LEMMA","LEMMA");
    
    t_new->Branch( "siiev", &siiev_new);
    t_new->Branch( "iev", &iev_new);
    t_new->Branch( "nhits", &nhits_new);
    t_new->Branch( "subdet", subdet_new, "subdet[100]/I");
    t_new->Branch( "xh", xh_new, "xh[100]/F");
    t_new->Branch( "xh_PulseHeight",             xh_new_PulseHeight,                  "xh_PulseHeight[100]/F"                                               );
    t_new->Branch( "yh", yh_new, "yh[100]/F");
    t_new->Branch( "yh_PulseHeight",             yh_new_PulseHeight,                  "yh_PulseHeight[100]/F"                                               );
    t_new->Branch( "zh", zh_new, "zh[100]/F");
    t_new->Branch( "itrack", itrack_new, "itrack[100]/I");
    t_new->Branch( "S2_PulseHeight",            &S2_new_PulseHeight,                  "S2_PulseHeight/I"             );
    t_new->Branch( "LemmaTrigger_PulseHeight",  &LemmaTrigger_new_PulseHeight,        "LemmaTrigger_PulseHeight/I"   );   
    t_new->Branch( "Deva_PulseHeight",           Deva_new_PulseHeight,                "Deva_PulseHeight[6]/I"        );
    t_new->Branch( "Cherenkov_PulseHeight",      Cherenkov_new_PulseHeight,           "Cherenkov_PulseHeight[4]/I"   );
    t_new->Branch( "S3_PulseHeight",            &S3_new_PulseHeight,                  "S3_PulseHeight/I"             );
    t_new->Branch( "PbGlass_PulseHeight",       &PbGlass_new_PulseHeight,             "PbGlass_PulseHeight/I"        );
    t_new->Branch( "S2_TimeSpectrum",           &S2_new_TimeSpectrum,                 "S2_TimeSpectrum/I"            );
    t_new->Branch( "LemmaTrigger_TimeSpectrum", &LemmaTrigger_new_TimeSpectrum,       "LemmaTrigger_TimeSpectrum/I"  );
    t_new->Branch( "Deva_TimeSpectrum",          Deva_new_TimeSpectrum,               "Deva_TimeSpectrum[6]/I"       );
    t_new->Branch( "Cherenkov_TimeSpectrum",     Cherenkov_new_TimeSpectrum,          "Cherenkov_TimeSpectrum[4]/I"  );
    t_new->Branch( "S3_TimeSpectrum",           &S3_new_TimeSpectrum,                 "S3_TimeSpectrum/I"            );
    t_new->Branch( "PbGlass_TimeSpectrum",      &PbGlass_new_TimeSpectrum,            "PbGlass_TimeSpectrum/I"       );        
    
    
    for (int k=0; k < 100; k++) {
        
        xh_new[k] = -999;                  
        yh_new[k] = -999;
        zh_new[k] = -999;
        itrack_new[k] = -999;
        subdet_new[k] = -999;
        xh_new_PulseHeight[k] = -999;
        yh_new_PulseHeight[k] = -999;
        
    }
    
    vector<int> iev_vec ; 
    
    f_new->cd();
    
    for (int i=0; i < t2->GetEntries(); i++) iev_vec.push_back(i);
    
    Int_t evt=0;
    
    for (int i=0; i < t1->GetEntries(); i++){ //loop on neg tree entries
        
        t1->GetEntry( i );
        
        evt++;
        
        if (i%1000==0) cout << "evt = " << evt << endl;
        
        for (int j=0; j < iev_vec.size(); j++){ //loop on pos tree entries
            
            t2->GetEntry( iev_vec.at(j) );
            
            if ( iev1 == iev2 ) {
                                
                for (int k=0; k <= j; k++) {
                    iev_vec.erase( iev_vec.begin() );
                }
                
                siiev_new = siiev1;
                iev_new = iev1;
                
                if (nhits1 > 100 || nhits2 > 100) break;
                
                Int_t imu=0;
                
                Float_t xh2_mu[12];   
                Float_t yh2_mu[12]; 
                Float_t zh2_mu[12]; 
                Int_t itrack2_mu[12]; 
                Int_t subdet2_mu[12];
                
                for (int k=0; k < nhits2; k++) {
                    
                    if (subdet2[k]==70) { 
                        
                        xh2_mu[imu] = xh2[k];   
                        yh2_mu[imu] = yh2[k]; 
                        zh2_mu[imu] = zh2[k]; 
                        itrack2_mu[imu] = itrack2[k]; 
                        subdet2_mu[imu] = subdet2[k];
                        imu++; 
                        
                    }
                }
                
                nhits_new = nhits1+imu;
                
                if (nhits_new > 100) nhits_new = 100;
                
                for (int k=0; k < nhits1; k++) {
                    
                    xh_new[k] = xh1[k];            //if(xh_new[k] > 0) cout << "xh1 neg " << xh1[k] << endl;       
                    yh_new[k] = yh1[k];
                    zh_new[k] = zh1[k];
                    itrack_new[k] = -1;
                    subdet_new[k] = subdet1[k];
                    xh_new_PulseHeight[k] = xh1_PulseHeight[k];
                    yh_new_PulseHeight[k] = yh1_PulseHeight[k];
                
                }
                
                
                for (int k=nhits1; k < nhits_new; k++) {
                    
                    xh_new[k] = xh2_mu[k-nhits1];   //if(xh_new[k] != -999) cout << "xh2 pos " << xh_new[k] << endl;               
                    yh_new[k] = yh2_mu[k-nhits1];
                    zh_new[k] = zh2_mu[k-nhits1];
                    itrack_new[k] = +1;
                    subdet_new[k] = subdet2_mu[k-nhits1];                    
                    
                }
                
                S2_new_PulseHeight  = S21_PulseHeight;
                S2_new_TimeSpectrum = S21_TimeSpectrum;
                LemmaTrigger_new_PulseHeight  = LemmaTrigger1_PulseHeight;
                LemmaTrigger_new_TimeSpectrum = LemmaTrigger1_TimeSpectrum;
                
                for( int k=0; k<6; k++ ) {
                    Deva_new_PulseHeight [k] = Deva1_PulseHeight [k];
                    Deva_new_TimeSpectrum[k] = Deva2_TimeSpectrum[k];
                }
                
                for( int k=0; k<4; k++ ) {
                    Cherenkov_new_PulseHeight [k] = Cherenkov1_PulseHeight [k];
                    Cherenkov_new_TimeSpectrum[k] = Cherenkov1_TimeSpectrum[k];
                }
                
                S3_new_PulseHeight  = S31_PulseHeight;
                S3_new_TimeSpectrum = S31_TimeSpectrum;
                PbGlass_new_PulseHeight  = PbGlass1_PulseHeight;
                PbGlass_new_TimeSpectrum = PbGlass1_TimeSpectrum;
                
                t_new->Fill();
                
                
                for (int k=0; k < 100; k++) {
                    
                    xh_new[k] = -999;                  
                    yh_new[k] = -999;
                    zh_new[k] = -999;
                    itrack_new[k] = -999;
                    subdet_new[k] = -999;
                    xh_new_PulseHeight[k] = -999;
                    yh_new_PulseHeight[k] = -999;
                    
                }
                
                break;
                
            } //close if iev1==iev2
            
        } //end loop on t2 entries
        
    } //end loop on t1 entries
    
    
    t_new->Write();
    f_new->Close();
    
    return;
    
}
