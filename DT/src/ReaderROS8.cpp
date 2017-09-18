#include "ReaderROS8.h"

// output flags
static const bool DUMP_HISTOS = 0;   // fill root histograms
static const bool CREATE_TREE = 1;   // fill root file RADMU
static const bool DUMP_STAT   = 0;   // save statistics on a txt file

static const double ConvToNs = 25./32.;

//define map array and function to fill it
static const int nros =   1;
static const int nrob =  19;
static const int ntdc =   4;
static const int ncha =  32;

ReaderROS8::ReaderROS8() {

  fo_txt=NULL;
  fo_Tree=new TFile();
  fo_Tree=NULL;
  tree=NULL;
  fo_Histo=new TFile();
  fo_Histo=NULL;

  inout=new Track_IO();
  geo=new Geom();        // Geo and reference system
  corr= new TimeCorr();  // Time Correction
  corr->InitSpline();    // Load spline for linear correction

  // init TOMTOOL
  _rawHistos = NULL;
  _ttrigCalib = NULL;

  //ALTEA
  _occupancy = NULL;

  // SV new ntuplizer
  m_ntuplizer = NULL;

  return;
}

ReaderROS8::~ReaderROS8() {


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


  if(CREATE_TREE){
    fo_Tree->cd();
    tree->Write();
    fo_Tree->Close();
    delete fo_Tree;
    fo_Tree = NULL;
  }

  delete dump;
  dump=NULL;

  return;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ReaderROS8::goUnpack(TString fin, int maxEvent, int runN, int runTrig, bool ttrig, bool n2chambers, int chside) {
    
    /// NB chside=-1 reconstruct only negative side of phi orizontal coordinate, chside=1 only positive, chside=0 all chamber
    m_chside = chside;

    /// --------- dumps initialization
    dump=new Save_HistosAndTree();
    if(!n2chambers)
        dump->set1chamber();

    if(CREATE_TREE){
        fo_Tree=inout->openOUTRootFile(runN,maxEvent,m_chside);
        dump->initTree();
        tree = new TTree("RADMU","radmu analysis");
        dump->bookTree(tree,n2chambers);
    }

    if(DUMP_HISTOS){
        fo_Histo=inout->openOUTHistoFile(runN,maxEvent);
    }

    if(DUMP_STAT){
        dump->init_Statistics();
    }

    /// --------- maps initialization
    fillMap();                        //fill channel map
    fillMap_t0();                     //fill t0 map
    fillMap_ttrig(runTrig);           //fill ttrig map

    if(DEBUG_RA_MAP)
        checkMap();
    if(DEBUG_RA_MAP_t0)
        checkMap_t0();
    if(DEBUG_RA_MAP_ttrig)
        checkMap_ttrig(runTrig);
    if(DUMP_HISTOS)
        dump->initHistos();

    /// --------- ttrig calibration 
    if(ttrig==1)
    {
        if(!_ttrigCalib)
        _ttrigCalib = new TTrigCalibration( runN );  //ALTEA
        _ttrigCalib->setFlagFillHistos(true);	
        setTTrigCalibPtr(_ttrigCalib);

        //ALTEA
        _occupancy = new Occupancy();

        // SV 20170613 FIX not sure this is needed....
        TCanvas * fC1 = new TCanvas("ReadbackDisplayProgram", "Tomography Display", 5, 5, 800, 900);
        _ttrigCalib->buildCanvasTSL(fC1);
        _ttrigCalib->setCanvas(1); //ALTEA IMPORTANTE!
        //in TTrigCalibration::setCanvas(int flag) viene inizializzata a 1 la variabile _nCanvas
        //che poi viene restituita dalla funzione TTrigCalibration::getCanvas()
        //con cui viene inizializzata la flag di TTrigCalibration::computeTTrig(int flag)
        //per cui alla fine vengono calcolati i ttrig dei SL e non quelli delle rob.
				   
    }

    /// --------- read raw file
    const char *fileName = fin.Data();
    FILE * rawfile = fopen(fileName,"rb");
    printf("Opening file: %s\n",fileName);
    m_integrity = 1;
    int numEvent = 0;

    if(rawfile==NULL) {
        printf("ERROR file %s doesn't exist ! \n",fileName);
    }
    else {
        while ( ! feof (rawfile) && m_integrity ){
            if(m_debug || numEvent/10000. == int(numEvent/10000.))
                cout << "Unpacking event " << dec << numEvent << "..." << endl;

            // create hit collection
            HITCollection * hits = new HITCollection();

            // unpack event
            int daqEvNum = readEvent(rawfile,hits);

            if(m_ntuplizer && hits->Get_NumberHITS()>0)
                m_ntuplizer->fillHits(hits,daqEvNum);
            else{
                 // track reconstruction
                Track *track = new Track();
                bool n2chambers=0;
                track->SelectTrack(hits,corr,n2chambers);

                if(DUMP_HISTOS)
                    dump->dumpHisto(track,hits,daqEvNum);
                if(CREATE_TREE)
                    dump->dumpTree(track,hits,daqEvNum,tree);
                if(DUMP_STAT)
                    dump->compute_Statistics(track,hits,daqEvNum);

                // fill ttrig calibration histos + occupancy(ALTEA)
                if(_ttrigCalib && _ttrigCalib->flagFillHistos() && hits) {
                    _ttrigCalib->fillHistos(hits);
                    _occupancy->fillOccHistos(hits);
                }
                delete track;
            }

            // delete collections
            delete hits;

            // increment event number and check integrity
            numEvent ++;
            if(numEvent >= maxEvent) m_integrity = false;
        }
        fclose (rawfile);
    } // end file reading

    /// --------- dump statistics
    if(DUMP_STAT){
      fo_txt=inout->openTXTFile(runN,maxEvent);
      if(fo_txt!=NULL){
        dump->write_Statistics(fo_txt);
        fclose (fo_txt);
      }
    }

    /// --------- ttrig dump
    if(_ttrigCalib && _ttrigCalib->flagFillHistos()){
        _ttrigCalib->computeTTrig(_ttrigCalib->getCanvas());
        char  ttrigName[100];
        sprintf(ttrigName,"./ttrig/ttrig_%d.txt",runN);
        _ttrigCalib->dumpTTrigs(ttrigName);
        //ALTEA
    	_occupancy->saveOccHistos();

    }
    cout << "END file unpacking, " << dec << numEvent << " events read" << endl;
    return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ReaderROS8::fillMap(){

  ifstream file("utils/mappaMB2.txt");

  if(!file.is_open()){
    cout << "ERROR : no channel map file found - exiting ! " << endl;
    exit(1);
  }

  //ALTEA per la nuova mappa devo cambiare queste variabili
  int ros=1;
  int rob;
  int tdc;
  int cha;

  int se=11;
  int sl;
  int lay;
  int tube;

  while( file >> rob >> tdc >> cha >> sl >> lay >> tube )
  {
      int idTDC = getTDCid(ros,rob,tdc,cha);
      int idTube = getTubeId(se,sl,lay,tube);

      chmap.insert( std::make_pair( idTDC, idTube ) );
      
  } 
  
  return;
}



void ReaderROS8::checkMap(){

  ofstream myfile;
  myfile.open ("check.txt");
  if(m_debug)
    myfile << "checking map" << endl;

  ifstream file("utils/mappaMB2.txt");
  
  if(!file.is_open()){
    cout << "ERROR : no channel map file found - exiting ! " << endl;
    exit(1);
  }

  int ros=1;
  int rob;
  int tdc;
  int cha;

  int se=11;
  int sl;
  int lay;
  int tube;
  
  map<int,int>::iterator iter;

  if(m_debug)
    cout << "map loop ... " << endl;
  while( file >> rob >> tdc >> cha >> sl >> lay >> tube  )
  {
      int rawdata = getTDCid(1,rob,tdc,cha);

      iter = chmap.find(rawdata);
      
      if( iter != chmap.end() )
      {
          if(m_debug)
          {
          	cout << "SE " << dec << getSe(iter->second) <<
                " sl " << getSL(iter->second) <<
                " lay " << getLay(iter->second) <<
                " tube " << getTube(iter->second) << endl;
          }

          if(m_debug)
          {
        	myfile << "fillMap " << ros << " " << rob << " " << tdc << " " << cha  << " => " <<	getSe(iter->second) << " " 
            << " " <<  getSL(iter->second) << " " << getLay(iter->second) << " " << getTube(iter->second) << endl;
          }
       }
    }

  if(m_debug)
    myfile << "end map check" << endl;
  
  return;
}



void ReaderROS8::fillMap_t0()
{
  ifstream file("utils/t0.txt");
  
  if(!file.is_open())
    cout << "WARNING : no t0 file found... " << endl;

  int se;
  int sl;
  int lay;
  int tube;
  float t0;
  float s_t0;


  if(m_debug)
    cout << "map_t0 loop ... " << endl;
  
  while( file >> se >> sl >> lay >> tube >> t0 >> s_t0 )
  {
      int idTube = getTubeId(se,sl,lay,tube);
      t0map.insert( std::make_pair( idTube , t0) );
  }
  
  return;
}



void ReaderROS8::checkMap_t0()
{

  ofstream myfile;
  myfile.open ("check2.txt");
  if(m_debug || DEBUG_RA_MAP_t0)
    myfile << "checking map t0" << endl;

  ifstream file("utils/t0.txt");
  if(!file.is_open())
    cout << "WARNING : no t0 file found... " << endl;

  int se;
  int sl;
  int lay;
  int tube;
  float t0;
  float s_t0;

  map<int,float>::iterator iter_t0;

  if(m_debug || DEBUG_RA_MAP_t0)
    cout << "map_t0 loop ... " << endl;
  
  while( file >> se >> sl >> lay >> tube >> t0 >> s_t0)
  {
      int t0data = getTubeId(se,sl,lay,tube);
      
      iter_t0 = t0map.find(t0data);
      
      if( iter_t0 != t0map.end() )
      {
          if(m_debug || DEBUG_RA_MAP_t0)
          {
              cout << "SE " << dec << getSe(iter_t0->first) <<
              " sl " << getSL(iter_t0->first) <<
              " lay " << getLay(iter_t0->first) <<
              " tube " << getTube(iter_t0->first) <<
              " time "<< iter_t0->second << endl; 
          }
      }
  }
  
  if(m_debug || DEBUG_RA_MAP_t0)
    myfile << "end map_t0 check" << endl;
  
  return;
}



void ReaderROS8::fillMap_ttrig(int runTrig)
{
  char fn[100];
  sprintf(fn,"./ttrig/ttrig_%d.txt",runTrig);
  ifstream file(fn);
  if(!file.is_open())
    cout << "WARNING : no ttrig file found... " << endl;

  int se;
  int sl;
  float ttrig;
  float s_ttrig;

  if(m_debug)
    cout << "map_ttrig loop ... " << endl;
  
  while( file >> se >> sl >> ttrig >> s_ttrig)
  {
      int idSL = getSLId(se,sl);

      ttrigmap.insert( std::make_pair( idSL , ttrig) );

  }

  return;
}



void ReaderROS8::checkMap_ttrig(int runTrig)
{
  ofstream myfile;
  myfile.open ("check2.txt");
  if(m_debug || DEBUG_RA_MAP_ttrig)
    myfile << "checking map ttrig" << endl;

  char fn[100];
  sprintf(fn,"./ttrig/ttrig_%d.txt",runTrig);
  ifstream file(fn);
  if(!file.is_open())
    cout << "WARNING : no ttrig file found... " << endl;

  int se;
  int sl;
  float ttrig;
  float s_ttrig;

  map<int,float>::iterator iter_ttrig;

  if(m_debug || DEBUG_RA_MAP_ttrig)
    cout << "map_ttrig loop ... " << endl;
  
  while( file >> se >> sl >> ttrig >> s_ttrig)
  {
      int ttrigdata = getSLId(se,sl);
      //cout << "TDCid " << rawdata << endl;

      iter_ttrig = ttrigmap.find(ttrigdata);
      if( iter_ttrig != ttrigmap.end() )
      {
          if(m_debug || DEBUG_RA_MAP_ttrig)
          {
              cout << "SE " << getSe(iter_ttrig->first) <<
              " sl " << getSL(iter_ttrig->first) <<
              " time "<< iter_ttrig->second << endl;
          }
      } 
  }
  
  if(m_debug || DEBUG_RA_MAP_ttrig)
    myfile << "end map_ttrig check" << endl;
  
  return;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int ReaderROS8::readEvent(FILE *infile, HITCollection * hits){

    // total word counter
    unsigned int wordCount = 0;

    /// first word is event word count
    unsigned int eventWordCountHeader = readWord(infile,wordCount);
    
    if(m_debug)
        cout << "\nEVENT HEADER, word count header " << dec << eventWordCountHeader << endl;
    
    /// event header: word from 1 to 9
    unsigned int runNum     =  readWord(infile,wordCount);
    unsigned int evNum      =  readWord(infile,wordCount);
    unsigned int spillNum   =  readWord(infile,wordCount);
    unsigned int res1       =  readWord(infile,wordCount);
    unsigned int rosOffset  =  readWord(infile,wordCount);
    unsigned int puOffset   =  readWord(infile,wordCount);
    unsigned int res2       =  readWord(infile,wordCount);
    unsigned int res3       =  readWord(infile,wordCount);

    if(m_debug){
        cout << "Reading event header " << endl;
        cout << "Run " << dec << runNum << endl;
        cout << "event number " << dec << evNum << endl;
        cout << "rosOffset " << dec << rosOffset << endl;
        cout << "puOffset " << dec << puOffset << endl;
        cout << "wordCount " << dec << wordCount << endl;
    }

    /// ROS board data
    readROS(infile,wordCount,rosOffset,hits);
    
    /// PU data 
    readPU(infile,wordCount,puOffset);

    /// last word is event word count
    unsigned int eventWordCountTrailer = readWord(infile,wordCount);
    
    if(m_debug)
        cout << "EVENT TRAILER, word count " << dec << eventWordCountTrailer << endl;

    // end event cross-check
    if(eventWordCountHeader != eventWordCountTrailer) {
        cout << "ERROR in data unpacking : event header word counter " << dec << eventWordCountHeader << " != event trailer word counter "
             << eventWordCountTrailer << endl;
        if (!m_integrity) {
            cout << "-> Don't panic, file ended before" << endl;
        }
        else {
            cout << "-> Serious problem in this file" << endl;
        }
    }
    return evNum;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ReaderROS8::readROS(FILE * infile, unsigned int & wordCount, unsigned int & rosOffset, HITCollection * hits){
       
    if(m_debug)
        cout << "START ROS board data unpacking!" << endl;

    unsigned int rosWordCount =  readWord(infile,wordCount);
    
    if(rosOffset+1 != wordCount)
        cout << "ERROR in data unpacking : rosOffset " << dec << rosOffset << " != word counter " << wordCount << endl;

    unsigned int rosWordLock =  readWord(infile,wordCount);
    
    if(m_debug)
        cout << "ROS board data, ROS word count " << dec << rosWordCount << ", lock status " << hex <<  rosWordLock << endl;
    
    int finalROSWordCount = wordCount+rosWordCount-2;

    if(m_debug)
        cout << "finalROSWordCount = " << dec << finalROSWordCount << endl;
    
    do {
      unsigned int rosData, rosID, rosChID;
      rosData = readWord(infile,wordCount);
      unsigned int typeRos = rosData & 0xF0000000;
      typeRos = typeRos >> 24;
      
      if (typeRos==0xF0) {
        rosID = rosData & 0x000007F8;
        rosID = rosID >> 3;
        rosChID = rosData & 0x00000007;
        if(m_debug) cout << "   ---> typeRos " << hex << typeRos << ", rosID " << dec << rosID <<", rosChID " << rosChID << ", wordCount " << wordCount << endl;
      }
      else {
        /// read TDC
        readTDCGroup(finalROSWordCount, rosData, infile, wordCount, rosChID, hits);
      }
    } while(wordCount < finalROSWordCount && m_integrity);

    if(finalROSWordCount != wordCount) {
        cout << "ERROR in data unpacking : ros word count " << dec << finalROSWordCount << " != word counter " << wordCount << endl;
        if (!m_integrity) {
            cout << "-> Don't panic, file ended before" << endl;
        }
        else {
            cout << "-> Serious problem in this file" << endl;
        }
    }

    if(m_debug)
        cout << "END ROS board data unpacking!" << endl;
    return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ReaderROS8::readTDCGroup(int finalROSWordCount, unsigned int tdcData, FILE * infile, unsigned int & wordCount, unsigned int rosChID, HITCollection * hits) {

    /// read TDC group type
    unsigned int typeTDC = tdcData & 0xF0000000;
	typeTDC = typeTDC >> 28;
	
	int evID = -1;
        /// group header
        if(typeTDC==0x0) {
            unsigned int tdcID = tdcData & 0x0F000000;
            tdcID = tdcID >> 24;
            evID = tdcData & 0x00FFF000;
            evID = evID >> 12;
            if(m_debug)
                cout << "      ---> TDC group header" << ", TDC ID " << dec << tdcID << ", event ID " << dec << evID << ", wordCount " << wordCount << endl;

            tdcData =  readWord(infile,wordCount);
            typeTDC = tdcData & 0xF0000000;
            typeTDC = typeTDC >> 28;

            /// leading measurement
            while(typeTDC == 0x4 || typeTDC==0x6 || typeTDC==0x7 || typeTDC==0xF){
                unsigned int tdcID, chID, rawtime;
                tdcID = tdcData & 0x0F000000;
                tdcID = tdcID >> 24;
                chID = tdcData & 0x00F80000;
                chID = chID >> 19;
                rawtime = tdcData & 0x0007FFFF;
                if(m_debug)
                    cout << "         ---> TDC data, TDC ID " << dec << tdcID << ", channel ID " << chID << ", time " << rawtime 
                    << ", wordCount " << wordCount << endl;

                // 20170725 LEMMA suppress no physical channel
//                 if(rosChID==4 && tdcID==3 && chID==0 && m_debug) {
//                     cout << "ERROR : no physical channel ! " << endl;
//                    return;
//                 }
                
                /// collect hit
                int rawdata = getTDCid(1,rosChID,tdcID,chID);
                map<int,int>::iterator iter = chmap.find(rawdata);

                // SV 100203 add severe failure in case no channel is found in map
                if(iter==chmap.end() && m_debug) {
                    cout << "ERROR : no channel is found in map -> rosChID " << dec << rosChID << ", tdcID " << tdcID << ", channel ID " << chID << " -> exit ! " << endl;
                    return;
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

                if( iter != chmap.end() && typeTDC!=0x6 && typeTDC!=0x7 && typeTDC!=0xF ) //<<<<<< QUA E' GIUSTO??
                {
                    float x_wire = geo->get_x_wire(getSe(iter->second),getSL(iter->second),getLay(iter->second),getTube(iter->second));
                    float y_wire = geo->get_y_wire(getSe(iter->second),getSL(iter->second),getLay(iter->second),getTube(iter->second));

                    // create HIT and add HIT to HITCollection
                    if(getSe(iter->second) !=1 ) {
                        int sl = getSL(iter->second);
                        if(m_chside==0) {
                            hits->createHIT(evID,getSe(iter->second),getSL(iter->second),getLay(iter->second),getTube(iter->second),x_wire,y_wire,tdcID,rosChID,chID,convToNs*(float(rawtime/4.)),t0,ttrig);
                        }
                        else if (sl==2) {
                            hits->createHIT(evID,getSe(iter->second),getSL(iter->second),getLay(iter->second),getTube(iter->second),x_wire,y_wire,tdcID,rosChID,chID,convToNs*(float(rawtime/4.)),t0,ttrig);
                        }
                        else if (m_chside*x_wire >0 && (sl==1 || sl==3)){
                            hits->createHIT(evID,getSe(iter->second),getSL(iter->second),getLay(iter->second),getTube(iter->second),x_wire,y_wire,tdcID,rosChID,chID,convToNs*(float(rawtime/4.)),t0,ttrig);
                        }
			
//                         if(m_debug) {
//                           cout<<"HIT created !"<< endl;
//                           cout << "SE " << dec << getSe(iter->second) <<
//                           " sl " << getSL(iter->second) <<
//                           " lay " << getLay(iter->second) <<
//                           " tube " << getTube(iter->second) <<
//                           " x_tube " << x_wire <<
//                           " y_tube " << y_wire <<
//                           " ROB_Id " << rosChID <<
//                           " TDC_Id " << tdcID <<
//                           " channel " << chID <<
//                           " time (ns) " << convToNs*(float(rawtime)/4.) <<
//                           " t0 (ns) " << convToNs*(iter_t0->second) <<
//                           " ttrig (ns) " << iter_ttrig->second <<
//                           " drift_time (ns) " << T_shift + convToNs*(float(rawtime)/4.) - t0 - ttrig << endl << endl;
//                         }
                    }
                } // close if( iter != chmap.end() )
                
                if(wordCount >= finalROSWordCount) return;
                
                tdcData = readWord(infile,wordCount);
                typeTDC = tdcData & 0xF0000000;
                typeTDC = typeTDC >> 28;
                
            } // close while

            /// now trailer comes...
            if(typeTDC==0x1){
                tdcID = tdcData & 0x0F000000;
                tdcID = tdcID >> 24;
                evID = tdcData & 0x00FFF000;
                evID = evID >> 12;
                if(m_debug)
                        cout << "      ---> TDC group trailer" << ", TDC ID " << dec <<tdcID << ", event ID " << dec << evID 
                        << " wordCount " << wordCount << endl;
            } 
            else if(m_debug) {
                cout << "I read this word " << hex << tdcData << " and I don't know what it means" << endl;
            }
        } // end if TDC group
        else if (typeTDC==0x6 || typeTDC==0x7 || typeTDC==0xF ){
            if(m_debug)
                cout << "  ---> (MARCO) Error, debugging data or empty ROS!" << endl;
        }

        return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ReaderROS8::readPU(FILE * infile, unsigned int & wordCount, unsigned int & puOffset){
    
    if(m_debug)
        cout << "START PU data unpacking!" << endl;

    unsigned int puWordCount =  readWord(infile,wordCount);
    
    if(puOffset+1 != wordCount){
        cout << "ERROR in data unpacking : puOffset " << dec << puOffset << " != word counter " << wordCount << endl;
        if (!m_integrity) {
            cout << "-> Don't panic, file ended before" << endl;
        }
        else {
            cout << "-> Serious problem in this file" << endl;
        }
    }
    
    if(m_debug)
        cout << "PU data, PU word count " << dec << puWordCount << endl;

    for(int ipw=1; ipw<puWordCount; ipw++)
        unsigned int puData =  readWord(infile,wordCount);

    if(m_debug)
        cout << "END reading PU data " << endl;

    return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

unsigned int ReaderROS8::readWord(FILE * infile, unsigned int & wordCount){

    //NB size_t fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
    //Reads data from the given stream into the array pointed to by ptr.
    //It reads nmemb number of elements of size size (in byte).
    //The total number of bytes read is (size*nmemb).
    //On success the number of elements read is returned.
    //On error or end-of-file the total number of elements successfully read (which may be zero) is returned.

    // get the 32-bit event buffer
    unsigned int word = 0x00000000;
    if(fread(&word,4,1,infile)==0)
        m_integrity = 0;
    else{
        wordCount++;
        m_integrity = 1;

//         if(m_debug) {
//             cout << "Reading word.... ---> (hex) " << hex << word << dec << "   " << endl;
//             cout << ", (dec) " << dec << word << endl;
//         }
    }

    return word;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int ReaderROS8::getTDCid( int ros, int rob, int tdc, int cha ) {
  return ( ( ( ( ( ros * nrob ) + rob ) * ntdc ) + tdc ) * ncha ) + cha;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int ReaderROS8::getTubeId( int se, int sl, int lay, int tube ) {
  return (se * 10000 + sl * 1000 + lay * 100 + tube);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int ReaderROS8::getTube(int idTube) {
  return (idTube - int(idTube/100.)*100.);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int ReaderROS8::getLay(int idTube) {
  return (int(idTube/100.) - int(idTube/1000.)*10);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int ReaderROS8::getSL(int idTube) {
  return (int(idTube/1000.) - int(idTube/10000.)*10);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int ReaderROS8::getSe(int idTube) {
  return int(idTube/10000.);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int ReaderROS8::getSLId( int se, int sl ) {
  return (se * 10000 + sl * 1000);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



