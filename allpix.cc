// ********************************************************************
//                                                     AllPix Geant4  *
//                 Generic Geant4 implementation for pixel detectors  *
//    Laboratoire de l'Acc�l�rateur Lin�aire Universit� Paris-Sud 11  *
//                                                                    *
//                           John Idarraga <idarraga@lal.in2p3.fr>    *
//                            Mathieu Benoit <benoit@lal.in2p3.fr>    *
// ********************************************************************

// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: exampleN06.cc,v 1.14 2006/06/29 17:53:52 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "AllPixDebug.h"

#include <time.h>
#include <string>

using namespace std;

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4ios.hh"

#include "AllPixPhysicsList.hh"
#include "AllPixDetectorConstruction.hh"
#include "AllPixPostDetConstruction.hh"
#include "AllPixSteppingVerbose.hh"
#include "AllPixRunAction.hh"
#include "AllPixPrimaryGeneratorAction.hh"
#include "AllPixEventAction.hh"
#include "AllPixRun.hh"
#include "AllPixRunAction.hh"

// digits, frames
#include "AllPix_Frames_WriteToEntuple.h"
// hits
#include "AllPix_Hits_WriteToEntuple.h"
#include "AllPixWriteROOTFile.hh" //nalipour
#include "Randomize.hh"


#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

/*
#if defined(G4UI_USE_TCSH)
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#elif defined(G4UI_USE_XM)
#include "G4UIXm.hh"
#elif defined(G4UI_USE_WIN32)
#include "G4UIWin32.hh"
#elif defined(G4UI_USE_QT)
#include "G4UIQt.hh"
#include "G4Qt.hh"
#else
#include "G4UIterminal.hh"
#endif
 */

#include "TH1F.h"

//#include "G4Qt.hh"
//#include "G4UIQt.hh"

// Input parameters
// indexes for argv[]
typedef enum { 
	_MACRO = 1,
	_RUN_BATCH // run batch mode in this case
} inputPars;

DebugLevel debug;

void checkflags(int,char**);
void SplashWindow();

// Check the options given on the command line and retrieve the values
bool cmdOptionExists(char** begin, char** end, const std::string& option);
char* getCmdOption(char** begin, char** end, const std::string & option);


#include "TMatrix.h"
#include "TRotMatrix.h"

// geometry
#include "ReadGeoDescription.hh"

int main(int argc, char** argv)
{

	// Flags
	checkflags(argc, argv);
	G4String fileName = argv[_MACRO];

	// Choose debug level, default error (0)
	debug = ERROR;
	if (cmdOptionExists(argv, argv + argc, "--debug")) {
		switch (std::stoi(getCmdOption(argv, argv+argc, "--debug")))
		{
					case 5:
						debug = INFO;
						break;
					case 10:
						debug = DEBUG;
						break;
					default:
						debug = ERROR;
		}

	}
	G4cout << " DEBUG level " << debug << endl;

	// Seed the random number generator manually
	time_t rawtime;
	time(&rawtime);
	G4long myseed = G4long(rawtime);
	SplashWindow();
	G4cout << "The random seed (localtime): " << myseed << G4endl;


	// User Verbose output class
	G4VSteppingVerbose* verbosity = new AllPixSteppingVerbose;
	G4VSteppingVerbose::SetInstance(verbosity);

	G4RunManager* runManager = new G4RunManager;

	// UserInitialization classes - mandatory;
	AllPixDetectorConstruction * detector =
			new AllPixDetectorConstruction();

	runManager->SetUserInitialization(detector);

	G4VUserPhysicsList * physics = new AllPixPhysicsList;
	runManager->SetUserInitialization(physics);

	// Hits ! --> Ntuple to store hits
	// creates AllPixRun to analyze hits at the end of event
	TString dataset = "allPix";
	TString tempdir = "";
	AllPixRunAction * run_action = new AllPixRunAction(detector, dataset, tempdir,
			"lciobridge_allpix.txt",
			"lciobridge_allpix_dut.txt"); // dataset
	runManager->SetUserAction(run_action);

	// Particle gun
	SourceType st = _GeneralParticleSource;
	//SourceType st = _HEPEvtInterface;
	AllPixPrimaryGeneratorAction * gen_action = new AllPixPrimaryGeneratorAction(st);
	runManager->SetUserAction(gen_action);

	// Get the PrimaryGeneratorMessenger to get a hold of some of its member variables
	run_action->GetPrimaryGeneratorMessenger(gen_action);

	// Digits ! --> calls Digitize, and makes ntuple to store digits
	AllPixEventAction * event_action = new AllPixEventAction(run_action);
	runManager->SetUserAction(event_action);

	// Initialize G4 kernel
	//
	// ATTENTION !  runManager->Initialize() calls AllPixDetectorConstruction::Construct()
	//  but the macro has not been read and the position of the Medipixes is decided in the macro.
	//  This automatic call won't have any effect.  It only build devices when the use calls
	//  /allpix/det/update from the macro.  AllPixDetectorConstruction::Construct() is called from
	//  AllPixDetectorConstruction::UpdateGeometry()
	//runManager->Initialize();

	// now done through the DetectorConstructor Messenger via the macro
	// Set digitizers.  I need to do it after
	//  DetectorConstruction::Construct() is called
	//  event_action->SetupDigitizers();
	//event_action->SetDetectorDigitInput(8.*keV); // thl !!!

	CLHEP::HepRandom::setTheSeed(myseed);
	//G4cout << *(G4Material::GetMaterialTable()) << G4endl;
#ifdef G4VIS_USE
	// Initialize visualization
	//
	G4VisManager* visManager = new G4VisExecutive;
	visManager->Initialize();
#endif

	// Get the pointer to the User Interface manager
	//
	G4UImanager* UI = G4UImanager::GetUIpointer();
	//UI->ApplyCommand("/tracking/verbose 1"); //nalipour: get information from secondary vertices???


	G4String command = "/control/execute ";

	//if (argc-1 == _RUN_BATCH)   // batch mode
	if (cmdOptionExists(argv, argv+argc, "--batch"))   // batch mode
	{
		//G4String command = "/control/execute ";
		//G4String fileName = argv[_MACRO];
		UI->ApplyCommand(command+fileName);
	}
	else
	{
		G4UIsession * session = 0;

		/*
#if defined(G4UI_USE_TCSH)
      session = new G4UIterminal(new G4UItcsh);      
#elif defined(G4UI_USE_XM)
      session = new G4UIXm(argc,argv);
      UI->ApplyCommand("/control/execute visTutor/gui.mac");      
#elif defined(G4UI_USE_WIN32)
      session = new G4UIWin32();
      UI->ApplyCommand("/control/execute visTutor/gui.mac");      
#elif defined(G4UI_USE_QT)
      session = new G4UIQt(argc,argv);
      UI->ApplyCommand("/control/execute visTutor/gui.mac");      
#else
      session = new G4UIterminal();
#endif


#ifdef G4VIS_USE
      UI->ApplyCommand(command+fileName);
      //AllPixPostDetConstruction::GetInstance()->WriteTracks("tracks_G4.root");
#endif
		 */


		//session = new G4UIQt(argc,argv);
		//session = new G4UIXm(argc, argv);
		session = new G4UIterminal();

		UI->ApplyCommand(command+fileName);
		//UI->ApplyCommand("/control/execute visTutor/gui.mac");

		session->SessionStart();
		delete session;
	}
	
	time_t simulationTime;
	time(&simulationTime);
	G4cout << " [INFO] Time to simulate " << simulationTime - rawtime << " s." << G4endl;

#ifdef G4VIS_USE
	delete visManager;
#endif

	// Geo description
	extern ReadGeoDescription * g_GeoDsc; // already loaded ! :)
	map<int, AllPixGeoDsc *> * geoMap = g_GeoDsc->GetDetectorsMap();
	map<int, AllPixGeoDsc *>::iterator detItr;

	// Frames ntuple closing
	// G4int nDigitizers = event_action->GetNumberOfDigitizers();
	for( detItr = geoMap->begin() ; detItr != geoMap->end() ; detItr++) {
	  WriteToNtuple::GetInstance("", "", "", (int)geoMap->size(), (*detItr).first)->closeNtuple();

	  //nalipour: Close and write the ROOT files
	  if(run_action->writeROOTFile !=NULL)
	    {
	      run_action->writeROOTFile[run_action->ReturnAllPixRun()->return_detIdToIndex((*detItr).first)]->AllPixCloseROOTFile();
	    }
	}

	// hits ntuple closing
	G4int nHC = event_action->GetNumberOfHC();
	for(G4int i = 0 ; i < nHC ; i++)
		Hits_WriteToNtuple::GetInstance("", "", "", nHC, i)->closeNtuple();

	G4cout << "[DONE]" << G4endl;

	delete runManager;
	delete verbosity;

	return 0;
}

void checkflags(int argc, char** argv){

	if(argc < 2){
		G4cout << "use: " << G4endl;
		//G4cout << "     " << argv[0] << " macro[filename]" << G4endl;
		G4cout << "     " << argv[0] << " macro/<filename> [options]" << G4endl;
		G4cout << "     additional options: --debug <level> where level is 0 (ERROR), 10 (DEBUG), or smth between (INFO)" << G4endl;
		G4cout << "                         --batch Runs simulation without graphical " << G4endl <<
				"                           window/driver." << G4endl;
		//G4cout << "                         --nEvents <n> Number of events to create." << G4endl; // not yet implemented
		exit(1);
	}

}

void SplashWindow(){

	G4cout << "*************************************************************" << G4endl;
	G4cout << "                                    AllPix Geant4" << G4endl;
	G4cout << "            Generic Geant4 implementation for pixel detectors" << G4endl;
	G4cout << "                                                             " << G4endl;
	G4cout << "                      John Idarraga, NIKHEF  <idarraga@CERN.CH>" << G4endl;
	G4cout << "                       Mathieu Benoit, UNIGE  <mbenoit@CERN.CH>" << G4endl;
	G4cout << "*************************************************************" << G4endl;
	G4cout << G4endl;

}

bool
cmdOptionExists(char** begin, char** end, const std::string& option)
{
  return std::find(begin, end, option) != end;
}

char*
getCmdOption(char** begin, char** end, const std::string & option)
{
  char** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)
    {
      return *itr;
    }
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
