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
// $Id: B2TrackerSD.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file B2TrackerSD.cc
/// \brief Implementation of the B2TrackerSD class

#include "RXGTrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "TTree.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RXGTrackerSD::RXGTrackerSD(const G4String& name,
		const G4String& hitsCollectionName)
: G4VSensitiveDetector(name),
  fHitsCollection(NULL)
{

	_primariesCount = 0;

	collectionName.insert(hitsCollectionName);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RXGTrackerSD::~RXGTrackerSD() 
{



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RXGTrackerSD::Initialize(G4HCofThisEvent* hce)
{

	// Create hits collection
	fHitsCollection
	= new RXGTrackerHitsCollection(SensitiveDetectorName, collectionName[0]);

	// Add this collection in hce

	G4int hcID
	= G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
	hce->AddHitsCollection( hcID, fHitsCollection );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool RXGTrackerSD::ProcessHits(G4Step* aStep, 
		G4TouchableHistory*)
{

	// energy deposit
	G4double edep = aStep->GetTotalEnergyDeposit();

	//if (edep==0.) return false;

	// Add edep
	G4Track * tr = aStep->GetTrack();
	G4int trackId = tr->GetTrackID();

	//_od->x.push_back( stepPoint->GetPosition().x() / CLHEP::mm );

	if( trackId == 1 ) {

		_od->edep.push_back( edep / CLHEP::keV );

		if ( _od->x.empty() ) _od->x.push_back( aStep->GetPreStepPoint()->GetPosition().x() / CLHEP::mm );
		if ( _od->y.empty() ) _od->y.push_back( aStep->GetPreStepPoint()->GetPosition().y() / CLHEP::mm );
		if ( _od->z.empty() ) _od->z.push_back( aStep->GetPreStepPoint()->GetPosition().z() / CLHEP::mm );

		//int si = _od->edep.size();
		//cout << si << endl;
		_primariesCount++;
	} else {
		_od->edepSec.push_back( edep / CLHEP::keV );
	}
	// store the pdgId info
	_od->pdgId.push_back( tr->GetParticleDefinition()->GetPDGEncoding( ) );

	RXGTrackerHit* newHit = new RXGTrackerHit();
	newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
	newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchableHandle()
			->GetCopyNumber());
	newHit->SetEdep(edep);
	newHit->SetPos (aStep->GetPostStepPoint()->GetPosition());
	fHitsCollection->insert( newHit );

	return true;
}
void RXGTrackerSD::EndOfEvent(G4HCofThisEvent*)
{
	if ( verboseLevel>1 ) {
		G4int nofHits = fHitsCollection->entries();
		G4cout << "\n-------->Hits Collection: in this event they are " << nofHits
				<< " hits in the tracker chambers: " << G4endl;
		for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
	}
	// Fill the Tree

	if ( ! _od->edep.empty() ) _od->nReach = 1;
	else _od->nReach = 0;

	_T->Fill();

	// Clean vars for next event
	_od->edep.clear();
	_od->edepSec.clear();
	_od->pdgId.clear();
	_od->x.clear();
	_od->y.clear();
	_od->z.clear();
	_od->nReach = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
