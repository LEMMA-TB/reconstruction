#include "Occupancy.h"

Occupancy::Occupancy(){
 
    // init all
    initVariables();
    initHistos();

    return;

}

Occupancy::~Occupancy(){

    for( int i=0; i<12; ++i ) {

	delete occLay[i];
    }

    return;

}

void Occupancy::initVariables(){

    occLayWidth = 100;
    occLayEdgeL = 0;
    occLayEdgeH = 100;

    return;
}

void Occupancy::initHistos(){

   TString occLayName;
   for( int i=0; i<12; ++i) {
        occLayName = "hOccupancy_SL_";
        if( i<4 ) occLayName += 1;
        if( i>=4 && i<8 ) occLayName += 2;
        if( i>=8 ) occLayName += 3;
        occLayName += "_Lay_";
        occLayName += (i%4 + 1);
        occLay[i] = new TH1F( occLayName, occLayName, occLayWidth, occLayEdgeL, occLayEdgeH );
    
    }
    
    Number_of_hits_per_event = new TH1F( "Number_of_hits_per_event", "Number of hits per event", 150, 0, 150 );
    
    return;

}


void Occupancy::fillOccHistos( HITCollection *hits ) {

    for( int i=0; i<hits->Get_NumberHITS(); ++i ) {

        HIT *hit = hits->hit(i);
	
        int lay=0;
        if( hit->CH_ID()==11 ) {
            if( hit->SL_ID()==1 ) {
            if( hit->L_ID()==1 ) lay=0;
            if( hit->L_ID()==2 ) lay=1;
            if( hit->L_ID()==3 ) lay=2;
            if( hit->L_ID()==4 ) lay=3;
	    }
            if( hit->SL_ID()==2 ) {
            if( hit->L_ID()==1 ) lay=4;
            if( hit->L_ID()==2 ) lay=5;
            if( hit->L_ID()==3 ) lay=6;
            if( hit->L_ID()==4 ) lay=7;
	    }
            if( hit->SL_ID()==3 ) {
            if( hit->L_ID()==1 ) lay=8;
            if( hit->L_ID()==2 ) lay=9;
            if( hit->L_ID()==3 ) lay=10;
            if( hit->L_ID()==4 ) lay=11;
	    }
	}
	//cout << "hit = " << hit->SL_ID() << " " <<  hit->wire_ID() << endl;
	occLay[lay]->Fill( hit->wire_ID() );
	
    }
    
    Number_of_hits_per_event->Fill( hits->Get_NumberHITS() );

    return;
}

void Occupancy::saveOccHistos() {

    for( int i=0; i<12; ++i ){

        occLay[i]->Write();

    }
    
    Number_of_hits_per_event->Write();

}


