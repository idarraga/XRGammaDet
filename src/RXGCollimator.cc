/**
 *  John Idarraga <idarraga@cern.ch>, 2015
 *
 */

#include "RXGDetectorConstruction.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Trd.hh"
#include "G4Sphere.hh"
#include "G4Polyhedra.hh"
#include "G4SubtractionSolid.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "TString.h"

using namespace CLHEP;

G4LogicalVolume * RXGDetectorConstruction::BuildCollimator(G4LogicalVolume * motherL, G4double centerZ) {

	double hc_rOuter_0 = 2*mm;
	double hc_rOuter_1 = hc_rOuter_0 * 0.5; // TODO ... calculate according to focus
	double hc_sep = 150*um;

	///////////////////////////////////////////////////////////////////////////////////////
	// Collimator Pb body
	G4Cons * collimatorBody = new G4Cons("CollimatorBody", 0, 10*cm, 0, 8*cm, 2*cm, pi, 2*pi);
	_CollimatorLogic = new G4LogicalVolume(
			collimatorBody,
			_beamColimatorMaterial,
			"CollimatorBodyLogic");

	G4ThreeVector collPos(25*cm,25*cm, centerZ);

	new G4PVPlacement(
			0,                // no rotation
			collPos,
			_CollimatorLogic, // its logical volume
			"CollimatorBody", // its name
			motherL,          // its mother  volume
			true,            // no boolean operations
			0,                // copy number
			true); 			  // checking overlaps

	// Visualization attributes
	G4VisAttributes* CollVisAtt = new G4VisAttributes( G4Colour(1,1,1, 1) );
	//CollVisAtt->SetForceSolid( true );
	CollVisAtt->SetForceWireframe( true );
	CollVisAtt->SetForceAuxEdgeVisible( true );
	_CollimatorLogic->SetVisAttributes( CollVisAtt );
	///////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////
	// Honeycomb core
	const G4double z[] = { 0, 4*cm };
	const G4double rInner[] = { 0, 0 };
	const G4double rOuter[] = { hc_rOuter_0, hc_rOuter_1 };
	G4Polyhedra * hcc = new G4Polyhedra(
			"oneTube",
			0.,		// initial phi starting angle
			2*pi,	// total phi angle
			6,		// number sides
			2,		// number of z planes
			z,
			rInner,
			rOuter
	);
	G4LogicalVolume * hccLogic = new G4LogicalVolume(
			hcc,
			_air,
			"HoneycombElement");

	// Visualization attributes
	G4VisAttributes* hccVisAtt = new G4VisAttributes( G4Colour::Blue() );
	hccVisAtt->SetForceSolid( true );
	//CollVisAtt->SetForceWireframe( true );
	hccVisAtt->SetForceAuxEdgeVisible( true );
	hccLogic->SetVisAttributes( hccVisAtt );

	// Substract the honey comb structure
	// Loop in y
	G4VSolid * substractColl = 0x0;
	G4ThreeVector orgPos(0,0, -2*cm);
	G4ThreeVector runPos = orgPos;
//	double inStep = 2*hc_rOuter_0 + hc_sep;
	double xStep = cos(30*pi/180) * (2*hc_rOuter_0 + hc_sep);
	double yStep = 2*hc_rOuter_0 + hc_sep;
	double posx = 0., posy = 0.;

	for ( int i = -20 ; i <= 20 ; i++ ) { // loop in x


		for ( int j = 20 ; j <= 20 ; j++ ) { // loop in y

//			posy = (double)j * inStep + ( (double)i*inStep/2. );
//			posx = i * inStep;
			if (i % 2 == 0) {
				posy = (double)j * yStep;
			}
			else {
				posy = (double)j * yStep + ( hc_rOuter_0 + hc_sep/2);
			}
			posx = i * xStep;
			runPos.setY( posy );
			runPos.setX( posx );

			TString hcc_name = "HoneycombElement_";
			hcc_name += i; hcc_name += "_"; hcc_name += j;

			// Condition to place or not a HoneyComp element
			if (sqrt(posy*posy + posx*posx) > 7.5*cm) {
				continue;
			}

			// Define rotation
//			G4RotationMatrix * pRot = new G4RotationMatrix;
//			pRot->rotateX( 20*deg );

			new G4PVPlacement(
					0,                // no rotation
					runPos,
					hccLogic, // its logical volume
					hcc_name.Data(), // its name
					_CollimatorLogic,          // its mother  volume
					true,            // no boolean operations
					0,                // copy number
					true); 			  // checking overlaps

			G4cout << "[COLL] Building collimator element " << i << "," << j << " | " << hcc_name << G4endl;


			/*
			// The first one
			if ( ! substractColl ) {
				substractColl = new G4SubtractionSolid("Coll", collimatorBody, hcc, 0x0, runPos );
			} else {
				substractColl = new G4SubtractionSolid("Coll", substractColl, hcc, 0x0, runPos );
			}
			 */

		}

	}

	/*
	_CollimatorLogic = new G4LogicalVolume(
			substractColl,
			_beamColimatorMaterial,
			"CollimatorBodyLogic");

	G4ThreeVector collPos(0,0, centerZ);
	new G4PVPlacement(
			0,                // no rotation
			collPos,
			_CollimatorLogic, // its logical volume
			"CollimatorBody", // its name
			motherL,          // its mother  volume
			true,            // no boolean operations
			0,                // copy number
			true); 			  // checking overlaps



	vector<G4LogicalVolume *> hcc_logic_v;
	 */


	/*

	// Build first
	G4ThreeVector hccCenterPos(0, 0, -2*cm);
	G4ThreeVector hccPos(0, 0, -2*cm);
	hcc_logic_v.push_back( PlaceAndHcc(hcc, hccPos) );

	// Build sides

	hccPos.setY( 2*hc_rOuter_0 + hc_sep );
	hcc_logic_v.push_back( PlaceAndHcc(hcc, hccPos) );
	hccPos.setY( -1 * hccPos.y() );
	hcc_logic_v.push_back( PlaceAndHcc(hcc, hccPos) );

	// up
	hccPos = hccCenterPos;
	hccPos.setY( (2*hc_rOuter_0 + hc_sep) / 2. );
	hccPos.setX( 2*hc_rOuter_0 + hc_sep );
	hcc_logic_v.push_back( PlaceAndHcc(hcc, hccPos) );
	hccPos.setY( -1 * hccPos.y() );
	hcc_logic_v.push_back( PlaceAndHcc(hcc, hccPos) );

	// down
	hccPos = hccCenterPos;
	hccPos.setY( (2*hc_rOuter_0 + hc_sep) / 2. );
	hccPos.setX( -1 * (2*hc_rOuter_0 + hc_sep) );
	hcc_logic_v.push_back( PlaceAndHcc(hcc, hccPos) );
	hccPos.setY( -1 * hccPos.y() );
	hcc_logic_v.push_back( PlaceAndHcc(hcc, hccPos) );




	G4VisAttributes* hccVisAtt = new G4VisAttributes( G4Colour(0.1, 0.1, 0.9, 0.5) );
	hccVisAtt->SetForceSolid( true );
	//hccVisAtt->SetForceWireframe( true );
	hccVisAtt->SetForceAuxEdgeVisible( true );


	vector<G4LogicalVolume *>::iterator hccItr = hcc_logic_v.begin();
	vector<G4LogicalVolume *>::iterator hccItrE = hcc_logic_v.end();

	for ( ; hccItr != hccItrE ; hccItr++ ) {
		(*hccItr)->SetVisAttributes( hccVisAtt );
	}

	 */


	return _CollimatorLogic;


}

G4LogicalVolume * RXGDetectorConstruction::PlaceAndHcc(G4Polyhedra * hcc, G4ThreeVector hccPos) {

	G4LogicalVolume * hccLogic = new G4LogicalVolume(
			hcc,
			_air,
			"hccLogic");

	new G4PVPlacement(
			0,                // no rotation
			hccPos,
			hccLogic, 	      // its logical volume
			"hcc",            // its name
			_CollimatorLogic, // its mother  volume
			false,            // no boolean operations
			0,                // copy number
			true); 			  // checking overlaps

	return hccLogic;
}
