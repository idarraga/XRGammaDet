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
// $Id: B2aDetectorConstruction.cc 77603 2013-11-26 17:11:49Z gcosmo $
//
/// \file B2aDetectorConstruction.cc
/// \brief Implementation of the B2aDetectorConstruction class

#include "TTree.h"

#include "RXGDetectorConstruction.hh"
#include "RXGDetectorMessenger.hh"
#include "RXGTrackerSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Sphere.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"

#include "G4SubtractionSolid.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


RXGDetectorConstruction::RXGDetectorConstruction()
: G4VUserDetectorConstruction()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RXGDetectorConstruction::~RXGDetectorConstruction()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* RXGDetectorConstruction::Construct()
{
	// Define materials
	DefineMaterials();

	// Define volumes
	G4VPhysicalVolume * py = DefineVolumes();

	return py;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RXGDetectorConstruction::DefineMaterials()
{
	// Material definition

	G4NistManager* nistManager = G4NistManager::Instance();
	nistManager->ListMaterials("all");

	// Air defined using NIST Manager
	_air = nistManager->FindOrBuildMaterial("G4_AIR");
	_NaI = nistManager->FindOrBuildMaterial("G4_SODIUM_IODIDE"); // Eventually try LSO too
	//_NaI = nistManager->FindOrBuildMaterial("G4_AIR");
	_RXReadOutChipMaterial = nistManager->FindOrBuildMaterial("G4_Si");
	_RXSensorMaterial = nistManager->FindOrBuildMaterial("G4_CESIUM_IODIDE");
	_beamColimatorMaterial = nistManager->FindOrBuildMaterial("G4_Pb"); // ("G4_Pb");
	//_beamColimatorMaterial = nistManager->FindOrBuildMaterial("G4_AIR");
	_Pb = nistManager->FindOrBuildMaterial("G4_Pb");

	// Print materials
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* RXGDetectorConstruction::DefineVolumes()
{

	G4double outSideHalfX = 10*cm;
	G4double outSideHalfY = 10*cm;
	G4double inSideHalfX = 8*cm;
	G4double inSideHalfY = 8*cm;

	G4double NaIHalfZ = 9*mm;
	G4double CollimatorHalfZ = 2*cm;
	G4double RXBackscatterShieldfHalfZ = (127/2.)*um;
	G4double ReadOutHalfZ = 1000*um;
	G4double RXSensorHalfZ = 500*um;  // FIXME try 300 too

	// The center of the NaI crystal is at the origin.  detectroStrecht is the distance from the
	//  origin to the first surface, i.e. the X-Ray detector.  The X-Ray source will be placed
	//  at 100cm from this point
	G4double detectroStrecht = ( NaIHalfZ + 2*CollimatorHalfZ + 2*RXBackscatterShieldfHalfZ + 2*ReadOutHalfZ + 2*RXSensorHalfZ );
	//G4double detectroStrecht = ( NaIHalfZ + 2*CollimatorHalfZ + 2*ReadOutHalfZ + 2*RXSensorHalfZ );
	G4double worldLength = 200*cm;  		// A hundred centimeters on each side
	//G4double worldLength = 400*cm;
	//		+
	//		(2 * detectroStrecht);			// twice the thickness of the detector.


	G4cout.precision(8);
	G4cout << "[INFO] For the source to be located at a 100cm from the surface of the X-Ray detector" << G4endl;
	G4cout << " it must be placed at " << (detectroStrecht + 100*cm)/mm << " mm" <<  endl;

	///////////////////////////////////////////////////////////
	// Definitions of Solids, Logical Volumes, Physical Volumes
	// World
	G4GeometryManager::GetInstance()->SetWorldMaximumExtent( worldLength );
	G4cout << "Computed tolerance = "
			<< G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
			<< " mm" << G4endl;
	G4Box* worldS
	= new G4Box("world",                                    // its name
			worldLength/2,worldLength/2,worldLength/2);     // its size
	G4LogicalVolume* worldLV
	= new G4LogicalVolume(
			worldS,    // its solid
			_air,      // its material
			"World");  // its name

	//  Must place the World Physical volume unrotated at (0,0,0).
	//
	G4VPhysicalVolume * worldPV
	= new G4PVPlacement(
			0,               // no rotation
			G4ThreeVector(), // at (0,0,0)
			worldLV,         // its logical volume
			"World",         // its name
			0,               // its mother  volume
			false,           // no boolean operations
			0,               // copy number
			true); // checking overlaps



	///////////////////////////////////////////////////////////
	// NaI Crystal
	G4Box  * NaISolid = new G4Box("NaISolid", outSideHalfX, outSideHalfY, NaIHalfZ);

	_NaILogic = new G4LogicalVolume(
			NaISolid,    // its solid
			_NaI,      // its material
			"NaILogic");  // its name

	new G4PVPlacement(
			0,               // no rotation
			G4ThreeVector(), // at (0,0,0)
			_NaILogic,         // its logical volume
			"Nai",         // its name
			worldLV,               // its mother  volume
			false,           // no boolean operations
			0,               // copy number
			true); // checking overlaps

	// Visualization attributes
	G4VisAttributes* NaIVisAtt = new G4VisAttributes( G4Colour::Yellow() );
	NaIVisAtt->SetForceSolid( true );
	//vitVisAtt->SetForceWireframe( true );
	NaIVisAtt->SetForceAuxEdgeVisible( true );
	_NaILogic->SetVisAttributes( NaIVisAtt );


	///////////////////////////////////////////////////////////
	// Collimator
	_CollimatorLogic = BuildCollimator( worldLV, 50*cm);//NaIHalfZ + CollimatorHalfZ ); //placement of collimator

	///////////////////////////////////////////////////////////
	// RX Backscatter shield
	// FIXME ... test and remove
	G4Box  * RXBackscatterShieldSolid = new G4Box("RXBackscatterShieldSolid", inSideHalfX, inSideHalfY, RXBackscatterShieldfHalfZ);
	_RXBackscatterShieldLogic = new G4LogicalVolume(
			RXBackscatterShieldSolid,    				// its solid
			_Pb,      									// its material
			"RXBackscatterShieldLogic");  			// its name
	G4ThreeVector bshieldPos(0,0, NaIHalfZ + 2*CollimatorHalfZ + RXBackscatterShieldfHalfZ);

	new G4PVPlacement(
			0,					// no rotation
			bshieldPos,
			_RXBackscatterShieldLogic,	// its logical volume
			"RXBackscatterShield",        // its name
			worldLV,            // its mother  volume
			false,           	// no boolean operations
			0,               	// copy number
			true); 				// checking overlaps

	// Visualization attributes
	G4VisAttributes* RXBShieldVisAtt = new G4VisAttributes( G4Colour::Blue() );
	RXBShieldVisAtt->SetForceSolid( true );
	//vitVisAtt->SetForceWireframe( true );
	RXBShieldVisAtt->SetForceAuxEdgeVisible( true );
	_RXBackscatterShieldLogic->SetVisAttributes( RXBShieldVisAtt );

	///////////////////////////////////////////////////////////
	// Si Read Out
	G4Box  * RXReadOutSolid = new G4Box("RXReadOutSolid", inSideHalfX, inSideHalfY, ReadOutHalfZ);
	_RXReadOutLogic = new G4LogicalVolume(
			RXReadOutSolid,    				// its solid
			_RXReadOutChipMaterial,      	// its material
			"RXReadOutLogic");  			// its name
	G4ThreeVector roPos(0,0, NaIHalfZ + 2*CollimatorHalfZ + 2*RXBackscatterShieldfHalfZ + ReadOutHalfZ);
	new G4PVPlacement(
			0,					// no rotation
			roPos,
			_RXReadOutLogic,	// its logical volume
			"RXReadOut",        // its name
			worldLV,            // its mother  volume
			false,           	// no boolean operations
			0,               	// copy number
			true); 				// checking overlaps

	// Visualization attributes
	G4VisAttributes* RXReadOutVisAtt = new G4VisAttributes( G4Colour::Black() );
	RXReadOutVisAtt->SetForceSolid( true );
	//vitVisAtt->SetForceWireframe( true );
	RXReadOutVisAtt->SetForceAuxEdgeVisible( true );
	_RXReadOutLogic->SetVisAttributes( RXReadOutVisAtt );


	///////////////////////////////////////////////////////////
	// RX Sensor
	G4Box  * RXSensorSolid = new G4Box("RXSensorSolid", inSideHalfX, inSideHalfY, RXSensorHalfZ);
	_RXSensorLogic = new G4LogicalVolume(
			RXSensorSolid,    				// its solid
			_RXSensorMaterial,      	// its material
			"RXSensorLogic");  			// its name
	G4ThreeVector rxSensPos(0,0, NaIHalfZ + 2*CollimatorHalfZ + 2*RXBackscatterShieldfHalfZ + 2*ReadOutHalfZ + RXSensorHalfZ);
	new G4PVPlacement(
			0,					// no rotation
			rxSensPos,
			_RXSensorLogic,	// its logical volume
			"RXSensor",        // its name
			worldLV,            // its mother  volume
			false,           	// no boolean operations
			0,               	// copy number
			true); 				// checking overlaps

	// Visualization attributes
	G4VisAttributes* RXSensorVisAtt = new G4VisAttributes( G4Colour::Red() );
	RXSensorVisAtt->SetForceSolid( true );
	//vitVisAtt->SetForceWireframe( true );
	RXSensorVisAtt->SetForceAuxEdgeVisible( true );
	_RXSensorLogic->SetVisAttributes( RXSensorVisAtt );


	// Always return the physical world

	return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RXGDetectorConstruction::ConstructSDandField()
{

	// Sensitive NaI
	// _NaILogic
	G4SDManager * SDman = G4SDManager::GetSDMpointer();
	RXGTrackerSD * NaI_SD = new RXGTrackerSD("NaI_SD", "HitsCollection_NaI");
	NaI_SD->SetOutputTree(_T[0], _od[0]);
	SDman->AddNewDetector( NaI_SD );
	_NaILogic->SetSensitiveDetector( NaI_SD );

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

