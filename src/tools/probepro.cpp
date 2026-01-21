
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <assert.h>
#include <unistd.h>
#include <sys/stat.h> 
#include <string.h>
#include <ctime>
#include <vector>
//#include "Defs.hpp"

using namespace std;

int main( const int argc, const char* argv[] ) { 

  try {

    // each arg should contain a probe filename...

    if (argc == 1) {
      cout << "Usage:\n\nprobepro.exe [args] <probe-filename> [<more-probe-filenames> ...]\n" << endl;
      return 0;
    }

    int first = -1;

    for (int iarg = 1; iarg < argc; ++iarg) {
      
      cout << "Working on probe file \"" << argv[iarg] << "\" step 1..." << endl;
      
      int first_step = -1;
      int last_step = -1;
      int inc = 0;
      bool got_hash = false;

      const int buffer_max = 256;
      char buffer[buffer_max];
      vector<int> last_step_vec;
      
      ifstream infile(argv[iarg]);
      while (infile) {
	infile.getline(buffer,buffer_max);
	if ((buffer[0] == '#')||(strlen(buffer) == 0)) {
	  got_hash = true;
	  continue;
	}
	// this is a valid line: now get the step number...
	if (first == -1) {
	  // flux probes start with "probe"...
	  if (buffer[0] == 'p')
	    first = 6;
	  else
	    first = 0;
	}
	int step = 0;
	char * str = buffer + first;
	while(*str != ' ') {
	  step = step*10 + (*str++ - '0');
	}

	if (got_hash) {
	  // coming out of a hash, things may be misaligned...
	  cout << " > block starts at step " << step << endl;
	  if (first_step != -1) {
	    assert(last_step != -1);
	    cout << " > data restarts at step: " << step << " step-last_step: " << step-last_step << " buffer \"" << buffer << "\"" << endl;
	    if ((step > first_step)&&(step <= last_step+inc)) {
	      cout << " > THIS IS A RESTART" << endl;
	      last_step_vec.push_back(step);
	    }
	    else {
	      cout << " > THIS IS NEW DATA" << endl;
	      last_step_vec.push_back(-1);
	    }
	  }
	  first_step = step;
	  got_hash = false;
	}
	else if (inc == 0) {
	  // when inc is zero, the step increment can be calculated...
	  inc = step - first_step;
	  cout << " > inc " << inc << endl; 
	}
	else {
	  if (step != last_step + inc) {
	    cout << "Error: inc changes: " << step - last_step << endl;
	    assert(0);
	  }
	}
	last_step = step;
      }

      cout << " > step2..." << endl;

      first_step = -1;
      last_step = -1;
      inc = 0;
      got_hash = false;
      
      infile.clear();
      infile.seekg(0);
      int iblock = 0;
      int ifile = 0;

      char ofname[256];
      sprintf(ofname,"%s-probepro",argv[iarg]);
      ofstream outfile(ofname);
      int outfile_step = -1; // for checking
      
      while (infile) {
	infile.getline(buffer,buffer_max);
	if ((buffer[0] == '#')||(strlen(buffer) == 0)) {
	  got_hash = true;
	  continue;
	}
	// this is a valid line: now get the step number...
	int step = 0;
	char * str = buffer + first;
	while(*str != ' ') {
	  step = step*10 + (*str++ - '0');
	}
	if (got_hash) {
	  // coming out of a hash, things may be misaligned...
	  cout << " > block starts at step " << step << endl;
	  if (first_step != -1) {
	    assert(iblock < last_step_vec.size());
	    if (last_step_vec[iblock] == -1) {
	      outfile.close();
	      ++ifile;
	      sprintf(ofname,"%s-probepro-%d",argv[iarg],ifile);
	      outfile.open(ofname);
	      outfile_step = -1; // reset check
	    }
	    ++iblock;
	  }
	  first_step = step;
	  got_hash = false;
	}
	else if (inc == 0) {
	  // when inc is zero, the step increment can be calculated...
	  inc = step - first_step;
	  cout << " > inc " << inc << endl; 
	}
	else {
	  if (step != last_step + inc) {
	    cout << "Error: inc changes: " << step - last_step << endl;
	    assert(0);
	  }
	}
	last_step = step;
	// now decide if we are going to write this piece of data...
	if ((iblock == last_step_vec.size())||(last_step_vec[iblock] == -1)||(step < last_step_vec[iblock])) {
	  // ensure contiguousness...
	  if (outfile_step == -1) {
	    outfile_step = step;
	  }
	  else {
	    assert(inc > 0);
	    assert(step = outfile_step + inc);
	    outfile_step = step;
	  }
	  // and write...
	  outfile << buffer << "\n";
	}
      }

      outfile.close(); 
      
    }

  } 
  catch(...) { 

    return 1;

  }

  return 0;

}



