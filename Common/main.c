////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 20170720 Sara Vanini Common software for LEMMA testbeam
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//---- STL
#include <iostream>
#include <fstream>
#include <map>

//---- Core
#include "src/Options.h"

//---- ROOT
#include <TError.h>

//---- LEMMA classes
#include "src/EventBuilder.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// RECIPE PARAMETERS //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace {
    static struct Parameters : Options {
        std::string rootDTFile;
        std::string inputSiFile;
        std::string outputFile;
        int runNum;
        int nEvents;
        std::string execute;
        bool alignment;
        bool debug;

        Parameters(const char *hello = "Program options") : Options(hello) {
        add_options()
        ("help", "printout help")

        // GENERAL //
        ( "rootDTFile",  &rootDTFile,  std::string("test"),         "DT  chamber reconstructed data file"                )
        ( "inputSiFile", &inputSiFile, std::string("test"),         "Si detectors raw data file"                         )
        ( "outputFile",  &outputFile,  std::string("test"),         "Output file name"                                   )
        ( "runNum",      &runNum,      (int)0.,                     "Run number ID"                                      )
        ( "nEvents",     &nEvents,     (int)10.,                    "number of events to read"                           )
        ( "execute",     &execute,     std::string("eventbuilder"), "Execution options: eventbuilder, PUT MORE IF NEEDED")
        ( "alignment",   &alignment,   (bool)0,                     "Alignment flag"                                     )
        ( "debug",       &debug,       (bool)0,                     "Debugging flag"                                     )
        ;
        }
    } p;   // <-- INSTANCE //
} // namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv) {

    //--- Avoid annoying ROOT reader warnings
    gErrorIgnoreLevel=kError;

    //--- Parameters
    std::string config_file("utils/parameters.cfg");
    p.add_options() ("config",&config_file,"set config file");
    p.parse_command_line(argc,argv);
    p.parse_config_file(config_file);
    p.parse_command_line(argc,argv);

    /// event builder
    if (p.execute == "eventbuilder"){
        std::cout << "\n\n *** RUNNING event builder *** \n\n";

        EventBuilder *builder = new EventBuilder();
        std::cout << "debug " << p.debug << std::endl;

        builder->setDebug(p.debug);
        builder->setAlignment(p.alignment);
        builder->openDataFiles(p.rootDTFile,p.inputSiFile,p.outputFile);
        builder->matchEvents(p.nEvents);
        builder->dumpOutput();

        delete builder;
    }

    return 0;
}





