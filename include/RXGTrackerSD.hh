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
// $Id: B2TrackerSD.hh 69706 2013-05-13 09:12:40Z gcosmo $
//
/// \file B2TrackerSD.hh
/// \brief Definition of the B2TrackerSD class

#ifndef RXGTrackerSD_h
#define RXGTrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "RXGTrackerHit.hh"

#include <vector>
using namespace std;

class G4Step;
class G4HCofThisEvent;
class TTree;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// B2Tracker sensitive detector class
///
/// The hits are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step. A hit is created with each step with non zero 
/// energy deposit.

typedef struct {
	vector<double> edep;
	vector<double> edepSec;
	vector<int> pdgId;
	vector<double> x;
	vector<double> y;
	vector<double> z;
	int nReach;
} Outdata;

class RXGTrackerSD : public G4VSensitiveDetector
{
  public:
    RXGTrackerSD(const G4String& name, 
                const G4String& hitsCollectionName);
    virtual ~RXGTrackerSD();
  
    // methods from base class
    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);
    void SetOutputTree(TTree * T, Outdata * od) { _T = T ; _od = od; };

  private:
    RXGTrackerHitsCollection* fHitsCollection;
    TTree * _T;
    Outdata * _od;

    // Primaries count
    int _primariesCount;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
