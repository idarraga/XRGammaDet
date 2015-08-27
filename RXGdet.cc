//
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

#include "TFile.h"
#include "TTree.h"

#include <unistd.h>
#include <sstream>
#include <string.h>

#include "RXGDetectorConstruction.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UIterminal.hh"

#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "QGSP_BERT_HP.hh"

//#include "G4StepLimiterPhysics.hh"
#include "RXGSteppingAction.hh"
#include "RXGActionInitialization.hh"

#include "RXGTrackerSD.hh"

#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "G4PhysListFactory.hh"

// Prototypes
void SaveRandomSeed(long int seed, TString name);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{

	// Seed the random number generator manually
	time_t rawtime;
	time(&rawtime);
	G4long myseed = G4long(rawtime);
	G4cout << "[INFO] The local time : " << myseed << G4endl;
	// Get the process id and use it to offset the time
	G4long pid = (G4long) getpid();
	G4cout << "[INFO] The pid { getpid() } : " << pid << G4endl;
	// Offset the local time
	myseed += pid;
	// Finally this will be the seed
	G4cout << "[INFO] The random seed (local time + pid): " << myseed << G4endl;
	// Save the seed to double check later
	SaveRandomSeed(myseed, "run");
	//choose the Random engine
	CLHEP::RanecuEngine * reng = new CLHEP::RanecuEngine();
	reng->setSeed(myseed);
	G4Random::setTheEngine(reng);

	Outdata * oDataNaI = new Outdata;
	Outdata * oDataChoroid = new Outdata;
	Outdata * oDataRetina = new Outdata;
	Outdata * oDataVitreous = new Outdata;

	TFile * oF = new TFile("output.root","RECREATE");
	TTree * oT_NaI = new TTree("TNaI","TNaI");
	TTree * oTChoroid = new TTree("TChoroid","TChoroid");
	TTree * oTRetina = new TTree("TRetina","TRetina");
	TTree * oTVitreous = new TTree("TVitreous","TVitreous");

	oT_NaI->Branch("edep", &(oDataNaI->edep));
	oT_NaI->Branch("edepSec", &(oDataNaI->edepSec));
	oT_NaI->Branch("pdgId", &(oDataNaI->pdgId));
	oT_NaI->Branch("nReach", &(oDataNaI->nReach));
	oT_NaI->Branch("x", &(oDataNaI->x));
	oT_NaI->Branch("y", &(oDataNaI->y));
	oT_NaI->Branch("z", &(oDataNaI->z));

	oTChoroid->Branch("edep", &(oDataChoroid->edep));
	oTChoroid->Branch("edepSec", &(oDataChoroid->edepSec));
	oTChoroid->Branch("pdgId", &(oDataChoroid->pdgId));

	oTRetina->Branch("edep", &(oDataRetina->edep));
	oTRetina->Branch("edepSec", &(oDataRetina->edepSec));
	oTRetina->Branch("pdgId", &(oDataRetina->pdgId));

	oTVitreous->Branch("edep", &(oDataVitreous->edep));
	oTVitreous->Branch("edepSec", &(oDataVitreous->edepSec));
	oTVitreous->Branch("pdgId", &(oDataVitreous->pdgId));
	oTVitreous->Branch("x", &(oDataVitreous->x));
	oTVitreous->Branch("y", &(oDataVitreous->y));
	oTVitreous->Branch("z", &(oDataVitreous->z));

	// Put them together to pass them through the DetectorConstruction
	TTree * oT[4] = { oT_NaI, oTChoroid, oTRetina, oTVitreous };
	Outdata * oData[4] = { oDataNaI, oDataChoroid, oDataRetina, oDataVitreous };

	// Construct the default run manager
	G4RunManager * runManager = new G4RunManager;

	// Set mandatory initialization classes
	RXGDetectorConstruction * det = new RXGDetectorConstruction();
	det->SetOutputTree(oT, oData);
	runManager->SetUserInitialization( det );
	// Certain properties for detector building
	det->SetOuterVolumesTransparent( true );
	det->SetOverlapCheck( false );

	G4VModularPhysicsList* physicsList = new QGSP_BERT_HP;
	//physicsList->RegisterPhysics(new G4StepLimiterPhysics());
	runManager->SetUserInitialization(physicsList);

	// Physics List
	//G4PhysListFactory factory;
	//G4VModularPhysicsList * phys = 0;
	//G4String physName = "FTFP_BERT";
	//G4String physName = "QGSP_BERT_EMY";
	// reference PhysicsList via its name
	//phys = factory.GetReferencePhysList(physName);
	//runManager->SetUserInitialization(phys);

	G4UserSteppingAction * stepping_action = new RXGSteppingAction();
	runManager->SetUserAction(stepping_action);

	runManager->SetUserInitialization( new RXGActionInitialization() );


	// Initialize G4 kernel
	runManager->Initialize();


#ifdef G4VIS_USE
	// Initialize visualization
	G4VisManager* visManager = new G4VisExecutive;
	// G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
	// G4VisManager* visManager = new G4VisExecutive("Quiet");
	visManager->Initialize();
#endif

	// Get the pointer to the User Interface manager
	G4UImanager* UImanager = G4UImanager::GetUIpointer();

	if (argc!=1)   // batch mode
	{
		G4String command = "/control/execute ";
		G4String fileName = argv[1];
		UImanager->ApplyCommand(command+fileName);
	}
	else
	{  // interactive mode : define UI session
#ifdef G4UI_USE
		//G4UIExecutive* ui = new G4UIExecutive(argc, argv);
		G4UIterminal * ui = new G4UIterminal();
#ifdef G4VIS_USE
		UImanager->ApplyCommand("/control/execute vis_gps.mac");
#else
		UImanager->ApplyCommand("/control/execute init.mac");
#endif
		//if (ui->IsGUI())
		//UImanager->ApplyCommand("/control/execute vis.mac");
		//ui->SessionStart();
		ui->SessionStart();
		delete ui;
#endif
	}

	// Job termination
	// Free the store: user actions, physics_list and detector_description are
	// owned and deleted by the run manager, so they should not be deleted
	// in the main() program !


	oF->cd();

	// Write Tree
	oT_NaI->Write();
	oTChoroid->Write();
	oTRetina->Write();
	oTVitreous->Write();

	// Close the File
	oF->Close();


#ifdef G4VIS_USE
	delete visManager;
#endif
	delete runManager;

	return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void SaveRandomSeed(long int seed, TString name) {

	TString ofn = "randseed_";
	ofn += name;
	ofn += ".txt";
	std::ofstream ofs (ofn.Data(), std::ofstream::out);
	ofs << seed << endl;
	ofs.close();

}
