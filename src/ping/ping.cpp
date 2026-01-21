#include "CTI.hpp"
using namespace CTI;

#include "Journal.hpp"
#include "ImageTools.hpp"

class Tycho : public ImageTools {

public:

  void run_journal(istringstream &instream) {
    if (mpi_rank==0) logger->setKillFilename("killping");

    /*
     *  Read Input File to input journal vector<string>, match input file exactly
     *  Use Params to parse and execute input journal
     *  Copy executed lines of input journal to output journal vector<string>
     */
    using namespace Journal;
    readInputJournal(instream);

    Param journalParam;
    bool areMoreJournalLines = parseParamFromJournal(journalParam);
    if (!areMoreJournalLines){
      Param helpParam("HELP");
      parse(&helpParam);
    }
    while (areMoreJournalLines){
      try {
        logCommandLine(0);
        parse(&journalParam);
        journalParam = Param();
        areMoreJournalLines = parseParamFromJournal(journalParam);
      }
      catch (...) {
        // if a catatrophic (i.e., surface not constructed) exception is thrown, exit surfer
        areMoreJournalLines = false;
      }
    }
    dumpOutputJournal();

  }

};

int main(int argc, char * argv[]) {

  try {

    // this call to CTI_Init will update the input file based on any command line input
    CTI_Init(argc,argv,"ping.in");

    /*
    FOR_ALL_PARAM {
      cout << "got param: " << param->str() << endl;
    }
    getchar();
    */

    {

      //Build a new string list of params including any command line
      //params for the Journal class.
      ostringstream ofullParamStream;
      dumpParams(ofullParamStream);
      //Next, delete Param list, params will be processed sequentially from the
      //ping.in file in run_journal().  Ultimately, CTI_Init should
      //be changed not to build the list in the first place; however,
      //to avoid changes to cti_core for now just delete and start over locally.
      
      // skip param delete as the first step in transition to param-driven
      // surfer-like architecture...
      /*
      FOR_ALL_PARAM {
        param->setFlag();
      }
      deleteFlaggedParams();
      */
      //END Delete Params list


      Tycho tb;
      istringstream ifullParamStream(ofullParamStream.str());
      tb.run_journal(ifullParamStream);

    }

    CTI_Finalize();

  }
  catch (int e) {
    if (e >= 0)
      CTI_Finalize();
    else
      CTI_Abort();
  }
  catch(...) {
    CTI_Abort();
  }

  return(0);

}
