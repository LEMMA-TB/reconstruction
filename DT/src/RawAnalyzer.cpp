#include "RawAnalyzer.h"

// output flags
static const bool DUMP_HISTOS = 1;   // fill root histograms
static const bool CREATE_TREE = 1; // fill root file RADMU
static const bool CREATE_HITBANK = 0; // fill HitBank file
static const bool DUMP_STAT = 1;   // save statistics on a txt file

static const double ConvToNs = 25./32.;

//define map array and function to fill it
static const int nros =   1;
static const int nrob =  19;
static const int ntdc =   4;
static const int ncha =  32;

RawAnalyzer::RawAnalyzer() {
  
  // init useful variables
  
  if(DEBUG_RA_FLAG||DEBUG_RA_MAP||DEBUG_RA_MAP_t0||DEBUG_RA_MAP_ttrig||DEBUG_RA_HIT)
    maxwords = 10000;
  else 
    maxwords = 1000000000;  

  nHEADER=0;
  ngheader=0;
  nTRAILER=0;
  ngtrailer=0;
  ngroup=0;
  nERROR=0;
  nDEBUG=0;
  check_header=0;
  
  numEvent=0;
  numEventDAQ=0;
  rawfile=NULL; 
  fo_txt=NULL; 
  fo_Tree=new TFile();
  fo_Tree=NULL; 
  tree=NULL;
  fo_Histo=new TFile();
  fo_Histo=NULL; 

  inout=new Track_IO();
  dump=new Save_HistosAndTree();
  geo=new Geom();        // Geo and reference system
  corr= new TimeCorr();  // Time Correction
  corr->InitSpline();    // Load spline for linear correction
  
  // init TOMTOOL 
  _rawHistos = NULL;
  _ttrigCalib=NULL;
  hits=NULL;
//   _imgAnalyzer =NULL;
    
  return;
}

RawAnalyzer::~RawAnalyzer() {
  
  
  delete geo;  
  delete corr;  
  
  if(DUMP_HISTOS){
    fo_Histo->cd();
    dump->writeHistos();    
    fo_Histo->Close();
    delete fo_Histo;
    fo_Histo = NULL;
    //     dump->resetHistos();	
    //     dump->deleteHistos(); 
  }
  
  if(CREATE_HITBANK){
    // closing raw data file
    char eor[4];
    sprintf(eor,"EOR");
    HBFile->write(eor,sizeof(eor));
    HBFile->close();
    if(DEBUG_HB)
      cout << " Closing HitBank " << eor << endl;
    delete HBFile;
    HBFile=NULL;
  }
  
  if(CREATE_TREE){
    fo_Tree->cd();
    tree->Write();
    fo_Tree->Close();
    delete fo_Tree;
    fo_Tree = NULL;
  }
  
  //   fclose(rawfile);
  
  delete dump;
  dump=NULL;
  
  return;
  
}

void RawAnalyzer::goAnalysis(TString fin, int maxEvent, int runN, int runTrig, bool ttrig) {

  if(CREATE_TREE){
    fo_Tree=inout->openOUTRootFile(runN,maxEvent);
    dump->initTree();
    tree = new TTree("RADMU","radmu analysis");
    dump->bookTree(tree);
  }

  if(CREATE_HITBANK){
    HBFile=new ofstream();
    HBFile=NULL; 
    HBFile=inout->openOUTHBFile(runN,maxEvent);
    
    dump->initHB(runN,numEvent,HBFile);
  }
  
  if(DUMP_HISTOS){
    fo_Histo=inout->openOUTHistoFile(runN,maxEvent);
  }
  
  if(DUMP_STAT){
    dump->init_Statistics();
  }
  
  fillMap();              //fill channel map
  fillMap_t0();           //fill t0 map
  fillMap_ttrig(runTrig); //fill ttrig map
  
  if(DEBUG_RA_MAP)
    checkMap();
  if(DEBUG_RA_MAP_t0)
    checkMap_t0();
  if(DEBUG_RA_MAP_ttrig)
    checkMap_ttrig(runTrig);
  
  if(DUMP_HISTOS){
    dump->initHistos();
  }
  
  if(ttrig==1)
    {  //Set up the Drawing Canvas & Pad Configuration
    if(!fC1)
      fC1 = new TCanvas("ReadbackDisplayProgram", "Tomography Display", 5, 5, 800, 900);
    if(!_ttrigCalib)
      _ttrigCalib = new TTrigCalibration( runTrig );
    _ttrigCalib->setFlagFillHistos(true);
    // pass CalibHistos pointer to CalibAnalyzer
    setTTrigCalibPtr(_ttrigCalib);
    _ttrigCalib->buildCanvasTSL(fC1);
    _ttrigCalib->setCanvas(1);
  }
  
  //NB size_t fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
  //Reads data from the given stream into the array pointed to by ptr.
  //It reads nmemb number of elements of size size (in byte).
  //The total number of bytes read is (size*nmemb).
  //On success the number of elements read is returned.
  //On error or end-of-file the total number of elements successfully
  //read (which may be zero) is returned.
  
  //NB ftell(infile) returns the current file position. For binary stream,
  //then the value is the number of bytes from the beginning of the file.
  
  
  if(DEBUG_RA_FLAG)
    maxwords = 10000;
  int words_read=0;
  
  for(int ID=0; ID<50; ID++){
    
    if(numEvent<=maxEvent) 
      {
    rawfile=inout->openINFile(fin,ID);
	
	if(rawfile==NULL)
	  continue;
	
	for(int words=1; words<maxwords; words++){
	  if(DEBUG_RA_FLAG)
	    cout << "Decoding event " << numEvent << "..." << endl;	
	  
	  if(numEvent<=maxEvent) 
	    {
	      evbu = 0x00000000;
	      type = 0x00000000;
	      TTCcount = 0x00000000;
	      dtime=999.;
	      
	      bool break_flag=0;
	      readWord(break_flag,rawfile);
	      if(break_flag)
		break;
	      
	      if(type == 0x1F) // Event header
		{
		  check_header=1;	  
		  //reset variables
		  Event_Id = 999;
		  Bunch_Id = 999;
		  ROB_Id = 999;

                  numEventDAQ       = TTCcount; // DAQ event numeration

		  numEvent ++;
		  if(numEvent/10000. == int(numEvent/10000.))
		    cout<<"Event number: "<<numEvent<<endl;
		  nHEADER ++;
		  hits = new HITCollection();
		}
	      
	      if(type == 0x3F) // Event trailer
		{
		  //reset variables
		  Event_Id = 999;
		  Bunch_Id = 999;
		  ROB_Id = 999;
		  
		  nTRAILER ++;
		  
          track=new Track();
          if(!hits)
            cout << "Sono in trailer, Hits nullo!!, check header " << check_header << endl;

          if(ttrig==0 && hits)
           {
              //cut sulle hit per doppia lettura front end //
              //hits->Clean2FE();

              track->SelectTrack(hits,corr,1);
              if(track->Track_IsGood()){
                  if(DUMP_HISTOS){

                      dump->dumpHisto(track,hits,numEvent);
                  }
                  if(CREATE_TREE)
                      dump->dumpTree(track,hits,numEventDAQ,tree);
                  if(CREATE_HITBANK)
                      dump->dumpHB(track,hits,numEvent,HBFile);
                  if(DUMP_STAT){
                      dump->compute_Statistics(track,hits,numEvent);
                  }
              }
          }

		  
		  // /////////////////////////////////////////////////////////
		  // SV START
		  
		  
		  // if _rawHistos has been initialized, and fill flag is true, go fill
		  
		  if(_rawHistos && _rawHistos->flagFillHistos() && hits){
		    _rawHistos->fillHistos(hits);
		    if(numEvent%100==0)
		      _rawHistos->updateCanvas(_rawHistos->getCanvas());
		  }
		  
		  // if _ttrigCalib has been initialized, and fill flag is true, go fill

		  if(_ttrigCalib && _ttrigCalib->flagFillHistos() && hits){
		    _ttrigCalib->fillHistos(hits);
		    if(numEvent%10000==0){
		      _ttrigCalib->computeTTrig(_ttrigCalib->getCanvas());
		      
		    }
		  }
		  
		  
		  if(ttrig==0)
		    {
		      delete track;
		      track=NULL;
		    }
              if(hits)
		    delete hits;
		  hits=NULL;

		  check_header=0;
		}
	      
	      if(type == 0xDF)
		nERROR ++;
	      
	      if(type == 0xFF)
		nDEBUG++;
	      
	      //	if(type != 0x1F && type != 0x3F && type != 0xDF && type != 0xFF )
	      if(type != 0x0000001F && type != 0x0000003F && type != 0x000000DF && type != 0x000000FF ){
		if(check_header==0){
		  printf("NEW EVENT WITHOUT EVENT HEADER!!! \n");
		  numEvent ++;
		  if(numEvent/10000. == int(numEvent/10000.))
		    cout<<"Event number: "<<numEvent<<endl;
		  nHEADER ++;
		  hits = new HITCollection();
		}
		
		collectHit();
		
		if(check_header==0)
		  check_header=1;
		
	      }
	      
	      words_read=words;	
	      
	    } // close if(numEvent...)
	  
	  else{
	    cout<<"words read = "<<words_read<<endl; 
	    cout<<"nHEADER = "<<nHEADER<<", nTRAILER = "<<nTRAILER
		<<", ngheader = "<<ngheader
		<<", ngtrailer = "<<ngtrailer
		<<", ngroup = "<<ngroup
		<<", nERROR = " << nERROR 
		<<", nDEBUG = " << nDEBUG << endl;
	    
	    // DO when all events are read:
	    if(_ttrigCalib && _ttrigCalib->flagFillHistos()){
	      // fC1->Clear();
	      // fC1->cd();
	      // _ttrigCalib->buildCanvasTSL(fC1);
	      // _ttrigCalib->setCanvas(1);
	      _ttrigCalib->computeTTrig(_ttrigCalib->getCanvas());
	      fC1->Update();
          char  ttrigName[100];
	      sprintf(ttrigName,"./ttrig/ttrig_%d.txt",runN); 
	      _ttrigCalib->dumpTTrigs(ttrigName);
	    } 
	    
	    // SV 100208 reconstruct image
	    //if(_imgAnalyzer && _imgAnalyzer->flagBuildImg()) ...
	    
	    if(DUMP_STAT){
          fo_txt=inout->openTXTFile(runN,maxEvent);
	      if(fo_txt!=NULL)
		{
		  dump->write_Statistics(fo_txt);
		  fclose (fo_txt);
		}
	    }
	    
	    return;
	  }
	  
	}//end loop on words
      }// close if(numEvent...)
    
    fclose(rawfile);
    rawfile==NULL;
    
  } // close if(ID=0;ID<runID_max;ID++)
  
  if(DUMP_STAT){
    fo_txt=inout->openTXTFile(runN,maxEvent);
    
    if(fo_txt!=NULL){
      dump->write_Statistics(fo_txt);
      fclose (fo_txt);
    }
  }
  
  cout<<"words read = "<<words_read<<endl; 
  cout<<"nHEADER = "<<nHEADER<<", nTRAILER = "<<nTRAILER
      <<", ngheader = "<<ngheader
      <<", ngtrailer = "<<ngtrailer
      <<", ngroup = "<<ngroup
      <<", nERROR = " << nERROR 
      <<", nDEBUG = " << nDEBUG << endl;
  
  // DO when all events are read:
  if(_ttrigCalib && _ttrigCalib->flagFillHistos()){
    // fC1->Clear();
    // fC1->cd();
    // _ttrigCalib->buildCanvasTSL(fC1);
    // _ttrigCalib->setCanvas(1);
    _ttrigCalib->computeTTrig(_ttrigCalib->getCanvas());
    fC1->Update();
    
    char  ttrigName[100];
    sprintf(ttrigName,"./ttrig/ttrig_%d.txt",runN); 
    _ttrigCalib->dumpTTrigs(ttrigName);
  }


  
  return;
}


void RawAnalyzer::fillMap(){
  
  ofstream myfile;
  myfile.open ("out.txt");
  if(DEBUG_RA_FLAG) 
    myfile << "reading map" << endl;
  
  ifstream file("utils/legnaro2ROS25.txt");//"legnaro2ROS25-prototipoFRANCO.txt");
  if(!file.is_open()){
    cout << "ERROR : no channel map file found - exiting ! " << endl;
    exit(1);
  }  
  
  int ros;
  int rob;
  int tdc;
  int cha;
  
  int se;
  int sl;
  int lay;
  int tube;
  
  int ddu;
  int wh;
  int st;
  
  if(DEBUG_RA_FLAG) 
    cout << "map loop ... " << endl;
  while( file >> ddu >> ros >> rob >> tdc >> cha >> wh >> st >> se >> sl >> lay >> tube ) 
    {
      if(DEBUG_RA_FLAG) myfile << "fillMap "
			    << ddu << " " << ros << " " << rob << " " << tdc << " " << cha  << " => "
			    << wh << " " << st << " " << se << " " << " " <<  sl << " " << lay << " " << tube << endl;
      
      int idTDC = getTDCid(ros,rob,tdc,cha);
      int idTube = getTubeId(se,sl,lay,tube);
      
      chmap.insert( std::make_pair( idTDC, idTube ) ); 
      
    }
  if(DEBUG_RA_FLAG) 
    myfile << "TDC map read" << endl;
  return;
}

void RawAnalyzer::checkMap(){

  ofstream myfile;
  myfile.open ("check.txt");
  if(DEBUG_RA_FLAG) 
    myfile << "checking map" << endl;
  
  ifstream file("legnaro2ROS25.txt");//"legnaro2ROS25-prototipoFRANCO.txt");
  if(!file.is_open()){
    cout << "ERROR : no channel map file found - exiting ! " << endl;
    exit(1);
  }  
  
  int ros;
  int rob;
  int tdc;
  int cha;
  
  int se;
  int sl;
  int lay;
  int tube;
  
  int ddu;
  int wh;
  int st;
  
  map<int,int>::iterator iter;
  
  if(DEBUG_RA_FLAG) 
    cout << "map loop ... " << endl;
  while( file >> ddu >> ros >> rob >> tdc >> cha >> wh >> st >> se >> sl >> lay >> tube ) 
    {
      int rawdata = getTDCid(1,rob,tdc,cha);
      //cout << "TDCid " << rawdata << endl;
      
      iter = chmap.find(rawdata);
      if( iter != chmap.end() )
	{
	  if(DEBUG_RA_FLAG)
	    {
	      cout << "SE " << getSe(iter->second) <<  
		" sl " << getSL(iter->second) << 
		" lay " << getLay(iter->second) << 
		" tube " << getTube(iter->second) << endl;
	    }
	  
	  if(DEBUG_RA_FLAG) 
	    myfile << "fillMap " << ddu << " " << ros << " " << rob << " " << 
	      tdc << " " << cha  << " => " << wh << " " << st << " " << 
	      getSe(iter->second) << " " << " " <<  getSL(iter->second) << " " << 
	      getLay(iter->second) << " " << getTube(iter->second) << endl;
	  
	}
    }
  if(DEBUG_RA_FLAG) 
    myfile << "end map check" << endl;
  return;
}

void RawAnalyzer::fillMap_t0()
{
  
  ofstream myfile;
  myfile.open ("out2.txt");
  if(DEBUG_RA_FLAG) 
    myfile << "reading t0" << endl;
  
  ifstream file("utils/t0.txt");
  if(!file.is_open()){
    cout << "WARNING : no t0 file found... " << endl;
  }  
  
  int se;
  int sl;
  int lay;
  int tube;
  float t0;
  float s_t0;
  
  
  if(DEBUG_RA_FLAG) 
    cout << "map_t0 loop ... " << endl;
  while( file >> se >> sl >> lay >> tube >> t0 >> s_t0) 
    {
      if(DEBUG_RA_FLAG) myfile << "fillMap t0 "
			    << se << " " << " " <<  sl << " " << lay << " " << tube << " " << t0 << endl;
      
      int idTube = getTubeId(se,sl,lay,tube);
      
      t0map.insert( std::make_pair( idTube , t0) ); 
      
    }
  if(DEBUG_RA_FLAG) 
    myfile << "t0 map read" << endl;
  return;
}

void RawAnalyzer::checkMap_t0()
{
  
  ofstream myfile;
  myfile.open ("check2.txt");
  if(DEBUG_RA_FLAG || DEBUG_RA_MAP_t0) 
    myfile << "checking map t0" << endl;
  
  ifstream file("utils/t0.txt");
  if(!file.is_open()){
    cout << "WARNING : no t0 file found... " << endl;
  }  
  
  
  int se;
  int sl;
  int lay;
  int tube;
  float t0;
  float s_t0;
  
  
  map<int,float>::iterator iter_t0;
  
  if(DEBUG_RA_FLAG || DEBUG_RA_MAP_t0) 
    cout << "map_t0 loop ... " << endl;
  while( file >> se >> sl >> lay >> tube >> t0 >> s_t0) 
    {
      int t0data = getTubeId(se,sl,lay,tube);
      //cout << "TDCid " << rawdata << endl;
      
      iter_t0 = t0map.find(t0data);
      if( iter_t0 != t0map.end() )
	{
	  if(DEBUG_RA_FLAG || DEBUG_RA_MAP_t0)
	    {
	      cout << "SE " << getSe(iter_t0->first) <<  
		" sl " << getSL(iter_t0->first) << 
		" lay " << getLay(iter_t0->first) << 
		" tube " << getTube(iter_t0->first) << 
		" time "<< iter_t0->second << endl;
	      //" time "<< t0 << endl;
	    }	  
	}
    }
  if(DEBUG_RA_FLAG || DEBUG_RA_MAP_t0) 
    myfile << "end map_t0 check" << endl;
  return;
}


void RawAnalyzer::fillMap_ttrig(int runTrig)
{
  
  ofstream myfile;
  myfile.open ("out3.txt");
  if(DEBUG_RA_FLAG) 
    myfile << "reading ttrig" << endl;
  
  char fn[100];
  sprintf(fn,"./ttrig/ttrig_%d.txt",runTrig);
  ifstream file(fn);
  if(!file.is_open()){
    cout << "WARNING : no ttrig file found... " << endl;
  }  
  
  int se;
  int sl;
  float ttrig;
  float s_ttrig;
  
  
  if(DEBUG_RA_FLAG) 
    cout << "map_ttrig loop ... " << endl;
  while( file >> se >> sl >> ttrig >> s_ttrig) 
    {
      if(DEBUG_RA_FLAG) myfile << "fillMap ttrig "
			    << se << " " << " " <<  sl << " " << ttrig << endl;
      
      int idSL = getSLId(se,sl);
      
      ttrigmap.insert( std::make_pair( idSL , ttrig) ); 
      
    }
  if(DEBUG_RA_FLAG) 
    myfile << "ttrig map read" << endl;
  return;
}

void RawAnalyzer::checkMap_ttrig(int runTrig)
{
  
  ofstream myfile;
  myfile.open ("check2.txt");
  if(DEBUG_RA_FLAG || DEBUG_RA_MAP_ttrig) 
    myfile << "checking map ttrig" << endl;
  
  char fn[100];
  sprintf(fn,"./ttrig/ttrig_%d.txt",runTrig);
  ifstream file(fn);
  if(!file.is_open()){
    cout << "WARNING : no ttrig file found... " << endl;
  }  
  
  
  int se;
  int sl;
  float ttrig;
  float s_ttrig;
  
  
  map<int,float>::iterator iter_ttrig;
  
  if(DEBUG_RA_FLAG || DEBUG_RA_MAP_ttrig) 
    cout << "map_ttrig loop ... " << endl;
  while( file >> se >> sl >> ttrig >> s_ttrig) 
    {
      int ttrigdata = getSLId(se,sl);
      //cout << "TDCid " << rawdata << endl;
      
      iter_ttrig = ttrigmap.find(ttrigdata);
      if( iter_ttrig != ttrigmap.end() )
	{
	  if(DEBUG_RA_FLAG || DEBUG_RA_MAP_ttrig)
	    {
	      cout << "SE " << getSe(iter_ttrig->first) <<  
		" sl " << getSL(iter_ttrig->first) << 
		" time "<< iter_ttrig->second << endl;
	    }
	  	  
	}
    }
  if(DEBUG_RA_FLAG || DEBUG_RA_MAP_ttrig) 
    myfile << "end map_ttrig check" << endl;
  return;
}


void RawAnalyzer::readWord(bool &break_flag, FILE *infile){
  
  // get the 32-bit word
  if(fread(&evbu,4,1,infile)==0)
    break_flag=1;
  
  // DMA: swap 16 lsb with 16 msb
  long word = 0x00000000;
  word = ((evbu>>16) & 0x0000FFFF) | ((evbu<<16) & 0xFFFF0000);
  
  if(DEBUG_RA_FLAG){
    cout << "Reading evbu ---> " << hex<<evbu<<dec << endl;
    cout << "Reading word ---> " <<hex<<word<<dec << endl;
    cout << "Reading type ---> " <<hex<<type<<dec << endl;
  }
  // type of packet and TTCcounts
  type = word & 0xFF000000;
  type = type >> 24;			// type of data packet: bit 24-31
  TTCcount = word & 0x00FFFFFF;	// TTC count: bit 0-23

  if(DEBUG_RA_FLAG){
    cout << "Reading word ---> " <<hex<<word<<dec << endl;
    cout << "Reading type ---> " <<hex<<type<<dec << endl;
  }
  
  if(DEBUG_RA_FLAG)
    {
      if(type == 0x1F){
	cout << " ---> Event header ! " << endl;
      }
      if(type == 0x3F){
	cout << " ---> Event trailer ! " << endl;
      }
      if(type == 0xDF)
	cout << " ---> Error Flag ! " << endl;
      
      if(type == 0xFF)
	cout << " ---> Debugging data ! " << endl;
    }

  return;
  
}

void RawAnalyzer::collectHit(){
  
  long group = type & 0xE0;
  group = group >> 5;
  
  if(DEBUG_RA_FLAG)
    cout << "group " << group << endl;
  switch(group)
    {
    case 5:
      if(DEBUG_RA_FLAG)
	cout << " Trailing measurement " << endl;
      break;	
      
    case 6:
      if(DEBUG_RA_FLAG)
	cout << " Errors " << endl;
      break;
      
    case 7:	
      if(DEBUG_RA_FLAG)
	cout << " Debugging data " << endl;	
      break;
      
    case 0:
      Event_Id = TTCcount & 0xFFF000; 
      Event_Id  = Event_Id >> 12;
      Bunch_Id = TTCcount & 0xFFF;
      ROB_Id = type & 0x1F;
      
      if(DEBUG_RA_FLAG)
	{
	  cout << "Group header: ";
	  cout << " Event_Id " << setbase(10) << Event_Id << 
	    " bunch " << Bunch_Id << 
	    " rob " << ROB_Id << endl;	
	}
      
      if(numEvent==0) ngheader=-1;
      ngheader ++;	
      
      break;
      
    case 1:
      Event_Id = TTCcount & 0xFFF000; 
      Event_Id  = Event_Id >> 12;
      Bunch_Id = TTCcount & 0xFFF;
      if(DEBUG_RA_FLAG)
	{
	  cout << "Group trailer: ";
	  cout << " Event_Id " << setbase(10) << Event_Id << 
	    " bunch " << Bunch_Id << 
	    " rob " << ROB_Id << endl;	
	}
      //reset variables
      Event_Id = 999;
      Bunch_Id = 999;
      ROB_Id = 999;					
      
      ngtrailer ++;
      
      break;
      
    case 2:
      Event_Id = TTCcount & 0xFFF000; 
      Event_Id  = Event_Id >> 12;
      Bunch_Id = TTCcount & 0xFFF;
      TDC_Id = ROB_Id & 0x3;	
      if(DEBUG_RA_FLAG)
	{
	  cout << "TDC Id " << setbase(10) << TDC_Id << " header: ";
	  cout << " Event_Id " << setbase(10) << Event_Id << 
	    " bunch " << Bunch_Id << 
	    " rob " << ROB_Id << endl;	
	}
      break;
      
    case 3:
      Event_Id = TTCcount & 0xFFF000; 
      Event_Id  = Event_Id >> 12;
      Bunch_Id = TTCcount & 0xFFF;
      TDC_Id = ROB_Id & 0x3;		
      if(DEBUG_RA_FLAG)
	{
	  cout << "TDC Id " << setbase(10) << TDC_Id << " trailer: ";
	  cout << " Event_Id " << setbase(10) << Event_Id << 
	    " bunch " << Bunch_Id << 
	    " rob " << ROB_Id << endl;	
	}
      break;
      
    case 4:	
      if(ROB_Id==999){
	if(DEBUG_RA_FLAG) 
	  cout<<"ROB_Id==999"<<endl; 
	break;
      }
      
      TDC_Id = type & 0x3;	
      channel = TTCcount & 0xF80000;
      channel = channel >> 19;
      rawtime = TTCcount & 0x7FFFF;
      
      // ROB: 0..24;  TDC:0..3;  ch:0..31
      int tdcflag = TDC_Id + 4*ROB_Id;
      int chflag = channel + 32*tdcflag;
      if(DEBUG_RA_FLAG)
	cout <<	"Event_Id " << Event_Id << " tdc " << TDC_Id << 
	  " rob " << ROB_Id << " ch " << channel << 
	  " tdcflag " << tdcflag << " chflag " << chflag << 
	  " numEvent " << numEvent << endl;
      
      int rawdata = getTDCid(1,ROB_Id,TDC_Id,channel);
      iter = chmap.find(rawdata);
      // SV 100203 add severe failure in case non channel is found in map 
      if(iter==chmap.end()) {
	cout << "ERROR : no channel is found in map - exit ! " << endl;
	cout <<	"Event_Id " << Event_Id << " tdc " << TDC_Id << 
	  " rob " << ROB_Id << " ch " << channel << endl;
        exit;
      }
      
      // SV 100203 get t0: set t0 from map, otherwise value is 0
      int t0data  = getTubeId(getSe(iter->second),getSL(iter->second),getLay(iter->second),getTube(iter->second));
      iter_t0 = t0map.find(t0data);
      float t0 = 0; 
      if( iter_t0 != t0map.end())
	if(iter_t0->second > -100 && (iter_t0->second)< 100) 
	  t0 = convToNs*(iter_t0->second);		      
      
      // SV 100203 get ttrig: set ttrig from map, otherwise value is 0
      int ttrigdata  = getSLId(getSe(iter->second),getSL(iter->second));
      iter_ttrig = ttrigmap.find(ttrigdata);
      float ttrig = 0; 
      if( iter_ttrig != ttrigmap.end())
        ttrig = iter_ttrig->second;	      
      
      if( iter != chmap.end() )
	{
	  
	  float x_wire = geo->get_x_wire(getSe(iter->second),getSL(iter->second),getLay(iter->second),getTube(iter->second));
	  float y_wire = geo->get_y_wire(getSe(iter->second),getSL(iter->second),getLay(iter->second),getTube(iter->second));
	  
	  // create HIT and add HIT to HITCollection	  
	  if(numEvent>=0 && (getSe(iter->second) !=1 )) {
	    hits->createHIT(numEvent,getSe(iter->second),getSL(iter->second),getLay(iter->second),getTube(iter->second),x_wire,y_wire,TDC_Id,ROB_Id,channel,convToNs*(float(rawtime/4.)),t0,ttrig);
	    
	    if(DEBUG_RA_HIT){
	      cout<<"\n Event: "<<numEvent<<", 1 HIT created"<<endl;
	    }		    
	  }
	  
	  if(DEBUG_RA_FLAG || DEBUG_RA_HIT)
	    {
	      cout << "SE " << getSe(iter->second) <<  
		" sl " << getSL(iter->second) << 
		" lay " << getLay(iter->second) << 
		" tube " << getTube(iter->second) << 
		" x_tube " << x_wire << 
		" y_tube " << y_wire << 
		" TDC_Id " << rawdata <<
		" ROB_Id " << ROB_Id <<
		" channel " << channel <<
		" time (ns) " << convToNs*(float(rawtime)/4.) << 
		" t0 (ns) " << convToNs*(iter_t0->second) << 
		" ttrig (ns) " << iter_ttrig->second << 
		" drift_time (ns) " << T_shift + convToNs*(float(rawtime)/4.) - t0 - ttrig << endl;
	    }		        
	} // close if( iter != chmap.end() )
      
      if(DEBUG_RA_FLAG)
	cout 	<< " Leading measurment " << " TDC " << setbase(10) << TDC_Id 
		<< " ch " << channel << " time " << rawtime << endl;
      //reset
      TDC_Id = 999;	
      channel = 999;
      rawtime = 999;
      
      ngroup ++;
      
      break;
      
    }//end switch
  
  
  return;
}


int RawAnalyzer::getTDCid( int ros, int rob, int tdc, int cha ) {
  return ( ( ( ( ( ros * nrob ) + rob ) * ntdc ) + tdc ) * ncha ) + cha;
}

int RawAnalyzer::getTubeId( int se, int sl, int lay, int tube ) {
  return (se * 10000 + sl * 1000 + lay * 100 + tube);
}

int RawAnalyzer::getTube(int idTube) {
  return (idTube - int(idTube/100.)*100.);
}

int RawAnalyzer::getLay(int idTube) {
  return (int(idTube/100.) - int(idTube/1000.)*10);
}

int RawAnalyzer::getSL(int idTube) {
  return (int(idTube/1000.) - int(idTube/10000.)*10);
}

int RawAnalyzer::getSe(int idTube) {
  return int(idTube/10000.);
}

int RawAnalyzer::getSLId( int se, int sl ) {
  return (se * 10000 + sl * 1000);
}



