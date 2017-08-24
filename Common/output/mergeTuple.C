
#include<TTree.h>
#include<TFile.h>
#include<iostream>

using namespace std;


// mergeTuple(file_name_neg,file_name_pos, merged)

void mergeTuple(TString file1, TString file2, TString file_new){

  
   Int_t           siiev1;
   Int_t           iev1;
   Int_t           nhits1;
   Int_t           subdet1[100];   //[nhits]
   Float_t         xh1[100];   //[nhits]
   Float_t         yh1[100];   //[nhits]
   Float_t         zh1[100];   //[nhits]
   Int_t           itrack1[100];   //[nhits]

   Int_t           siiev2;
   Int_t           iev2;
   Int_t           nhits2;
   Int_t           subdet2[100];   //[nhits]
   Float_t         xh2[100];   //[nhits]
   Float_t         yh2[100];   //[nhits]
   Float_t         zh2[100];   //[nhits]
   Int_t           itrack2[100];   //[nhits]


   Int_t           siiev_new;
   Int_t           iev_new;
   Int_t           nhits_new;
   Int_t           subdet_new[100];   //[nhits]
   Float_t         xh_new[100];   //[nhits]
   Float_t         yh_new[100];   //[nhits]
   Float_t         zh_new[100];   //[nhits]
   Int_t           itrack_new[100];   //[nhits]


   TFile *f1 = new TFile(file1);

   TTree *t1 = (TTree*) f1->Get("LEMMA");

   t1->SetBranchAddress("siiev", &siiev1);
   t1->SetBranchAddress("iev", &iev1);
   t1->SetBranchAddress("nhits", &nhits1);
   t1->SetBranchAddress("subdet", subdet1);
   t1->SetBranchAddress("xh", xh1);
   t1->SetBranchAddress("yh", yh1);
   t1->SetBranchAddress("zh", zh1);
   t1->SetBranchAddress("itrack", itrack1);

   TFile *f2 = new TFile(file2);

   TTree *t2 = (TTree*) f2->Get("LEMMA");

   t2->SetBranchAddress("siiev", &siiev2);
   t2->SetBranchAddress("iev", &iev2);
   t2->SetBranchAddress("nhits", &nhits2);
   t2->SetBranchAddress("subdet", subdet2);
   t2->SetBranchAddress("xh", xh2);
   t2->SetBranchAddress("yh", yh2);
   t2->SetBranchAddress("zh", zh2);
   t2->SetBranchAddress("itrack", itrack2);

   TFile *f_new = new TFile(file_new,"RECREATE");

   TTree *t_new = new TTree("LEMMA","LEMMA");


   t_new->Branch("siiev", &siiev_new);
   t_new->Branch("iev", &iev_new);
   t_new->Branch("nhits", &nhits_new);
   t_new->Branch("subdet", subdet_new, "subdet[100]/I");
   t_new->Branch("xh", xh_new, "xh[100]/F");
   t_new->Branch("yh", yh_new, "yh[100]/F");
   t_new->Branch("zh", zh_new, "zh[100]/F");
   t_new->Branch("itrack", itrack_new, "itrack[100]/I");

   for (int k=0; k < 100; k++) {

        xh_new[k] = -999;                  
        yh_new[k] = -999;
        zh_new[k] = -999;
        itrack_new[k] = -999;
        subdet_new[k] = -999;

   }


   std::vector<int> iev_vec ; 

   f_new->cd();

   for (int i=0; i <  t2->GetEntries(); i++) iev_vec.push_back(i);

   Int_t evt=0;

   for (int i=0; i < t1->GetEntries(); i++){ //loop on t1 entries

        t1->GetEntry(i);
        
        evt++;
        
        if (i%1000==0) std::cout << "evt = " << evt << std::endl;
        //std::cout << "evt size= " << iev_vec.size() << std::endl;
        //std::cout << "evt1 = " << iev1 << std::endl;
        //getchar();
        
        for (int j=0; j < iev_vec.size(); j++){ //loop on t2 entries
            
            t2->GetEntry(iev_vec.at(j));
            //std::cout << "evt2 = " << iev2 << std::endl;
            
            if (iev1 == iev2) {
                //cout << "found" << endl;
                
                for (int k=0; k <= j; k++) {
                    iev_vec.erase(iev_vec.begin());
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
                //cout << nhits_new << endl;
                
                if (nhits_new > 100) nhits_new = 100;
                
                for (int k=0; k < nhits1; k++) {
 
                      xh_new[k] = xh1[k];            //if(xh_new[k] > 0) cout << "xh1 neg " << xh1[k] << endl;       
                      yh_new[k] = yh1[k];
                      zh_new[k] = zh1[k];
                      itrack_new[k] = -1;
                      subdet_new[k] = subdet1[k];
                      
                }


                for (int k=nhits1; k < nhits_new; k++) {
 
                      xh_new[k] = xh2_mu[k-nhits1];   //if(xh_new[k] != -999) cout << "xh2 pos " << xh_new[k] << endl;               
                      yh_new[k] = yh2_mu[k-nhits1];
                      zh_new[k] = zh2_mu[k-nhits1];
                      itrack_new[k] = +1;
                      subdet_new[k] = subdet2_mu[k-nhits1];                    
                    
                }


                t_new->Fill();


                for (int k=0; k < 100; k++) {
                    
                    xh_new[k] = -999;                  
                    yh_new[k] = -999;
                    zh_new[k] = -999;
                    itrack_new[k] = -999;
                    subdet_new[k] = -999;
                    
                }

                break;
                
            } //close if iev1==iev2
            
        } //end loop on t2 entries

   } //end loop on t1 entries

 
   t_new->Write();
   f_new->Close();
   
   return;

}
