/* April/August 2015.
   This file is part of MACMA
   Class MacmaNonGui: Command line version of MACMA (author: leyaouanq@cervval.com)
  
   author: cecile.grigne@univ-brest.fr
*/

#include "MACMA/graphics/macmaNonGui.h"
#include "MACMA/physics/earth.h"

namespace fs = boost::filesystem;

using namespace std;


// ==================================================

MacmaNonGui::MacmaNonGui() : Macma()
{

  if(!checkFolder()) // prepare workspace
    exit(EXIT_FAILURE);

  _earth = new Earth(200.0);

  // not drawing cells in command line mode
  _earth->setDrawCells(false);

  // set to false for numerous runs done automatically
  _earth->setShowEvents(true);
}

MacmaNonGui::~MacmaNonGui()
{
}

// ==================================================
bool
MacmaNonGui::checkFolder()
{
  bool ready = false;

  fs::path current(fs::current_path());
  if(fs::exists(current))
    {
      bool empty = true;
      vector<fs::path> toRemove;
      for(fs::directory_iterator it(current); it!=fs::directory_iterator(); it++)
	{
	  if( is_directory(it->status()) 
	      && (it->path().filename() == "ages" || 
		  it->path().filename() == "logs" || 
		  it->path().filename() == "elements") )
	    {
	      empty = false;
	      toRemove.push_back(it->path());
	    }
	}

      if(empty)
	ready = prepareWorkspace();
      else
	{
	  cout << endl;
	  cout << "************************ WARNING **************************** " << endl;
	  cout << "Current workspace already contains folders ages, logs and elements." << endl;
	  cout << "These folders will be erased. Confirm [y or n]:\t ";
	  string input;
	  while(true)
	    {
	      getline(cin, input);
	      if(input.compare("y") == 0)
		{
		  for(unsigned int i=0; i<toRemove.size(); i++)
		    fs::remove_all(toRemove[i]); 
		  break;
		}
	      else if(input.compare("n") == 0)
		{
		  cout << "Abort... " << endl;
		  exit(EXIT_SUCCESS);
		}
	      else
		cout << "Please answer y or n." << endl;
	    }
	  ready = prepareWorkspace();
	}
    }
  else
    {
      cerr << "ERROR: cannot open current workspace " << current << endl;
      exit(EXIT_FAILURE);
    }

  return ready;
}
// ------------------------------------------------- 
bool
MacmaNonGui::prepareWorkspace()
{
  bool ready = true;

  fs::path dirAges("ages");
  if(fs::create_directory(dirAges))
    cout << " - - - folder ages/ was created" << endl; 
  else
    ready = false;

  fs::path dirLogs("logs");
  if(fs::create_directory(dirLogs))
    cout << " - - - folder logs/ was created" << endl; 
  else
    ready = false;

  fs::path dirElements("elements");
  if(fs::create_directory(dirElements))
    cout << " - - - folder elements/ was created" << endl;
  else
    ready = false;


  return ready;
}
// --------------------------------------------------
void
MacmaNonGui::run(string filename)
{
  openXmlFile(filename); /* read the configuration 
			    -> Makes _state = PLATES */

  if(_state < READY) 
    {
      _earth->makeReady();
      setMACMAState(READY);
    }

  if(_state == READY && _earth->isReady())
    {
      setMACMAState(RUNNING);
      _earth->setRunning(true);
      _earth->setRunStarted(true);

      _thread = boost::thread(&MacmaNonGui::_earthThread, this);
      _thread.join(); 

      _earth->setRunning(false);
      _earth->clear();
      setMACMAState(EMPTY);
    }
  else
    {
      cerr << "ERROR while initializing the structure" << endl;
      exit(EXIT_FAILURE);
    }
}
// --------------------------------------------------
void
MacmaNonGui::_earthThread()
{
  double writeDiagnosis = 50; // write a diagnosis every XX Myr

  while(true)
    {
      if(_state == RUNNING && _earth->isRunning())
	{
	  _earth->compute();

	  // diagnosis
	  if(_earth->getTimeMyr() >= _nextTimeWriteDiagnosis)
	    {
	      _earth->diagnosis(); 
	      _nextTimeWriteDiagnosis += writeDiagnosis;
	    }

	  if(_earth->getTimeMyr() >= _nextTimeWriteXml)
	    {
	      ostringstream oss;
	      oss << _workspace_ << "logs/param_"
		  << _earth->getTimeString();
	      writeXmlFile(oss.str());

	      oss.str("");

	      _nextTimeWriteXml += _earth->getWriteXml(); 
	    }

	  if( _earth->final_condition_is_reached() )
	    break;
	}
    }
}
// --------------------------------------------------
