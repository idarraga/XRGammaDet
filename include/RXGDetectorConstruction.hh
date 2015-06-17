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
// $Id: B2aDetectorConstruction.hh 73722 2013-09-09 10:23:05Z gcosmo $
//
/// \file B2aDetectorConstruction.hh
/// \brief Definition of the B2aDetectorConstruction class

#ifndef RXGDetectorConstruction_h
#define RXGDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "tls.hh"

#include "RXGTrackerSD.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4GlobalMagFieldMessenger;
class G4Polyhedra;

class RXGDetectorMessenger;
class TTree;

/// Detector construction class to define materials, geometry
/// and global uniform magnetic field.

class RXGDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    RXGDetectorConstruction();
    virtual ~RXGDetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();
    void SetOutputTree(TTree ** T, Outdata ** od) { _T = T ; _od = od; };
    G4LogicalVolume * BuildCollimator(G4LogicalVolume * motherL, G4double centerZ);
    G4LogicalVolume * PlaceAndHcc(G4Polyhedra *, G4ThreeVector);

  private:
    // methods
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();

    // Materials used
    G4Material * _air;
    G4Material * _NaI;
    G4Material * _Pb;
    G4Material * _RXReadOutChipMaterial;
    G4Material * _RXSensorMaterial;
    G4Material * _beamColimatorMaterial;

    // Logical volumes
    G4LogicalVolume * _NaILogic;
    G4LogicalVolume * _CollimatorLogic;
    G4LogicalVolume * _RXBackscatterShieldLogic;
    G4LogicalVolume * _RXReadOutLogic;
    G4LogicalVolume * _RXSensorLogic;

    // Output to pass to the SD
    Outdata ** _od;
    TTree ** _T;

};


//  Collimator geometry
//  Separation: 0.15 mm
//  Lenght: 30 mm
//  Diameter: 1.5 mm



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
