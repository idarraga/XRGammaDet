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
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4QuadrangularFacet.hh"
#include "G4SubtractionSolid.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "TString.h"

using namespace CLHEP;

G4LogicalVolume * RXGDetectorConstruction::BuildCollimator(G4LogicalVolume * motherL, G4double centerZ) {

	///////////////////////////////////////////////////////////////////////////////////////
	// Collimator Pb body
	//G4Cons * collimatorBody = new G4Cons("CollimatorBody", 0, 10*cm, 0, 8*cm, 2*cm, pi, 2*pi);
	G4Trd * collimatorBody = new G4Trd("CollimatorBody", 10*cm, 8*cm, 10*cm, 8*cm, 2*cm);
	_CollimatorLogic = new G4LogicalVolume(
			collimatorBody,
			_beamColimatorMaterial,
			"CollimatorBodyLogic");

	G4ThreeVector collPos(0, 0, centerZ);

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
//	CollVisAtt->SetForceSolid( true );
	CollVisAtt->SetForceWireframe( true );
	CollVisAtt->SetForceAuxEdgeVisible( true );
	_CollimatorLogic->SetVisAttributes( CollVisAtt );
	///////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////
	// Honeycomb core
	double hc_rOuter_0 = 1*mm;
	double hc_sep = 150*um;
	G4ThreeVector orgPos(0,0, 2*cm);
	double xStep = cos(30*pi/180) * (2*hc_rOuter_0 + hc_sep);
	double yStep = 2*hc_rOuter_0 + hc_sep;
	G4ThreeVector runPos = orgPos;

	for ( int i = -75; i <= 75 ; i++ ) { // loop in x
			for ( int j = -75 ; j <= 75 ; j++ ) { // loop in y

				double dy = hc_rOuter_0;
				double dx = tan(30*deg) * dy;
				double dz = 4*cm;
				double dmpx = 0*mm;
				double dmpy = 0*mm;

				if (i % 2 == 0) {
						dmpy = (double)j * yStep;
					}
					else {
						dmpy = (double)j * yStep + ( hc_rOuter_0 + hc_sep/2);
					}
				dmpx = i * xStep;

				runPos.setY( dmpy );
				runPos.setX( dmpx );

				// Check if honeycomb element is still within borders of collimator
				if ((abs(dmpx) > 7.75*cm) || (abs(dmpy) > 7.75*cm)){
					continue;
				}

				double mpxUp = 0;
				double mpyUp = 0;

				// Calculate offset of lower hexagon due to focussing
				double dmpxLow = dmpx * (100*cm + 3127*um + 4*cm) / (100*cm + 3127*um);
				double dmpyLow = dmpy * (100*cm + 3127*um + 4*cm) / (100*cm + 3127*um);
				double mpxLow = dmpxLow - dmpx;
				double mpyLow = dmpyLow - dmpy;

				// Make hcc as tessellated solid
				G4TessellatedSolid *hcc = new G4TessellatedSolid("hcc");

				// lower hexagon
				G4QuadrangularFacet *facet1 = new
				G4QuadrangularFacet (G4ThreeVector(mpxLow+dx, mpyLow+dy, -dz),
								   G4ThreeVector(mpxLow+dx, mpyLow-dy, -dz),
								   G4ThreeVector(mpxLow-dx, mpyLow-dy,-dz),
								   G4ThreeVector(mpxLow-dx, mpyLow+dy, -dz),

								   ABSOLUTE);
				G4TriangularFacet *facet2 = new
				G4TriangularFacet (G4ThreeVector(mpxLow+dx, mpyLow+dy, -dz),
						   	   	   G4ThreeVector(mpxLow+2*dx, mpyLow, -dz),
								   G4ThreeVector(mpxLow+dx, mpyLow-dy, -dz),
								   ABSOLUTE);
				G4TriangularFacet *facet3 = new
				G4TriangularFacet (G4ThreeVector(mpxLow-dx, mpyLow+dy, -dz),
								   G4ThreeVector(mpxLow-dx, mpyLow-dy, -dz),
								   G4ThreeVector(mpxLow-2*dx, mpyLow, -dz),
								   ABSOLUTE);


				// upper hexagon
				G4QuadrangularFacet *facet4 = new
				G4QuadrangularFacet (G4ThreeVector(mpxUp+dx, mpyUp+dy, 0),
								   G4ThreeVector(mpxUp-dx, mpyUp+dy,0),
								   G4ThreeVector(mpxUp-dx, mpyUp-dy, 0),
								   G4ThreeVector(mpxUp+dx, mpyUp-dy, 0),
								   ABSOLUTE);
				G4TriangularFacet *facet5 = new
				G4TriangularFacet (G4ThreeVector(mpxUp+dx, mpyUp+dy, 0),
								   G4ThreeVector(mpxUp+dx, mpyUp-dy, 0),
								   G4ThreeVector(mpxUp+2*dx, mpyUp, 0),
								   ABSOLUTE);
				G4TriangularFacet *facet6 = new
				G4TriangularFacet (G4ThreeVector(mpxUp-dx, mpyUp+dy, 0),
								   G4ThreeVector(mpxUp-2*dx, mpyUp, 0),
								   G4ThreeVector(mpxUp-dx, mpyUp-dy, 0),
								   ABSOLUTE);

				//side panels
				G4QuadrangularFacet *facet7 = new
				G4QuadrangularFacet (G4ThreeVector(mpxLow+dx, mpyLow+dy, -dz),
									 G4ThreeVector(mpxLow-dx, mpyLow+dy, -dz),
									 G4ThreeVector(mpxUp-dx, mpyUp+dy, 0),
									 G4ThreeVector(mpxUp+dx, mpyUp+dy, 0),
									 ABSOLUTE);
				G4QuadrangularFacet *facet8 = new
				G4QuadrangularFacet (G4ThreeVector(mpxLow-dx, mpyLow+dy, -dz),
									 G4ThreeVector(mpxLow-2*dx, mpyLow, -dz),
									 G4ThreeVector(mpxUp-2*dx, mpyUp, 0),
									 G4ThreeVector(mpxUp-dx, mpyUp+dy, 0),
									 ABSOLUTE);
				G4QuadrangularFacet *facet9 = new
				G4QuadrangularFacet (G4ThreeVector(mpxLow-2*dx, mpyLow, -dz),
									 G4ThreeVector(mpxLow-dx, mpyLow-dy, -dz),
									 G4ThreeVector(mpxUp-dx, mpyUp-dy, 0),
									 G4ThreeVector(mpxUp-2*dx, mpyUp, 0),
									 ABSOLUTE);
				G4QuadrangularFacet *facet10 = new
				G4QuadrangularFacet (G4ThreeVector(mpxLow-dx, mpyLow-dy, -dz),
									 G4ThreeVector(mpxLow+dx, mpyLow-dy, -dz),
									 G4ThreeVector(mpxUp+dx, mpyUp-dy, 0),
									 G4ThreeVector(mpxUp-dx, mpyUp-dy, 0),
									 ABSOLUTE);
				G4QuadrangularFacet *facet11 = new
				G4QuadrangularFacet (G4ThreeVector(mpxLow+dx, mpyLow-dy, -dz),
									 G4ThreeVector(mpxLow+2*dx, mpyLow, -dz),
									 G4ThreeVector(mpxUp+2*dx, mpyUp, 0),
									 G4ThreeVector(mpxUp+dx, mpyUp-dy, 0),
									 ABSOLUTE);
				G4QuadrangularFacet *facet12 = new
				G4QuadrangularFacet (G4ThreeVector(mpxLow+2*dx, mpyLow, -dz),
									 G4ThreeVector(mpxLow+dx, mpyLow+dy, -dz),
									 G4ThreeVector(mpxUp+dx, mpyUp+dy, 0),
									 G4ThreeVector(mpxUp+2*dx, mpyUp, 0),
									 ABSOLUTE);

				// Now add the facets to the solid
				//
				hcc->AddFacet((G4VFacet*) facet1);
				hcc->AddFacet((G4VFacet*) facet2);
				hcc->AddFacet((G4VFacet*) facet3);
				hcc->AddFacet((G4VFacet*) facet4);
				hcc->AddFacet((G4VFacet*) facet5);
				hcc->AddFacet((G4VFacet*) facet6);
				hcc->AddFacet((G4VFacet*) facet7);
				hcc->AddFacet((G4VFacet*) facet8);
				hcc->AddFacet((G4VFacet*) facet9);
				hcc->AddFacet((G4VFacet*) facet10);
				hcc->AddFacet((G4VFacet*) facet11);
				hcc->AddFacet((G4VFacet*) facet12);

				//Finally declare the solid is complete
				//
				hcc->SetSolidClosed(true);

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

				TString hcc_name = "HoneycombElement_";
				hcc_name += i; hcc_name += "_"; hcc_name += j;

				new G4PVPlacement(
						0,                // no rotation
						runPos,
						hccLogic, // its logical volume
						hcc_name.Data(), // its name
						_CollimatorLogic,          // its mother  volume
						true,            // no boolean operations
						0,                // copy number
						false);//true); 			  // checking overlaps

				G4cout << "[COLL] Building collimator element " << i << "," << j << " | " << hcc_name << G4endl;
			}
		}




	// Substract the honey comb structure
	// Loop in y
//	G4VSolid * substractColl = 0x0;

//	double inStep = 2*hc_rOuter_0 + hc_sep;
//	double xStep = cos(30*pi/180) * (2*hc_rOuter_0 + hc_sep);
//	double yStep = 2*hc_rOuter_0 + hc_sep;
//	double posx = 0., posy = 0.;
//
//	for ( int i = 0; i <= 0 ; i++ ) { // loop in x
//
//
//		for ( int j = 0 ; j <= 0 ; j++ ) { // loop in y
//
////			posy = (double)j * inStep + ( (double)i*inStep/2. );
////			posx = i * inStep;
//			if (i % 2 == 0) {
//				posy = (double)j * yStep;
//			}
//			else {
//				posy = (double)j * yStep + ( hc_rOuter_0 + hc_sep/2);
//			}
//			posx = i * xStep;
//			runPos.setY( posy );
//			runPos.setX( posx );
//
//			TString hcc_name = "HoneycombElement_";
//			hcc_name += i; hcc_name += "_"; hcc_name += j;
//
//			// Condition to place or not a HoneyComp element
//			//if (sqrt(posy*posy + posx*posx) > 7.5*cm) {
//			if ((abs(posy) > 6.5*cm) || (abs(posx) > 6.5*cm)){
//				continue;
//			}
//
//			// Define rotation
//			G4RotationMatrix * pRot = new G4RotationMatrix;
//			double thetay = atan2(posy,100);
//			pRot->rotateX( -thetay );
//			double thetax = atan2(posx,100);
//			pRot->rotateY( thetax );
//
//
//			new G4PVPlacement(
//					pRot,                // no rotation
//					runPos,
//					hccLogic, // its logical volume
//					hcc_name.Data(), // its name
//					_CollimatorLogic,          // its mother  volume
//					true,            // no boolean operations
//					0,                // copy number
//					true); 			  // checking overlaps
//
//			G4cout << "[COLL] Building collimator element " << i << "," << j << " | " << hcc_name << G4endl;
//
//
//
//			// The first one
////			if ( ! substractColl ) {
////				substractColl = new G4SubtractionSolid("Coll", collimatorBody, hcc, 0x0, runPos );
////			} else {
////				substractColl = new G4SubtractionSolid("Coll", substractColl, hcc, 0x0, runPos );
////			}
//
//
//		}
//
//	}

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
