/*
 * RXGCollimatorParameterization.cc
 *
 *  Created on: Jan 14, 2015
 *      Author: idarraga
 */

#include "RXGCollimatorParameterization.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Polyhedra.hh"
#include "G4SystemOfUnits.hh"

RXGCollimatorParameterization::RXGCollimatorParameterization(
		G4int nSpiral,
		G4int nCoresX,
		G4int nCoresY,
		G4double spacing
) : G4VPVParameterisation()
{

	_nSpiral = nSpiral;
	_nCoresX = nCoresX;
	_nCoresY = nCoresY;
	_spacing = spacing;

}

RXGCollimatorParameterization::~RXGCollimatorParameterization()
{ }

void RXGCollimatorParameterization::ComputeTransformation
(const G4int copyNo, G4VPhysicalVolume* physVol) const
{

	// X,Y positioning
	// copyNo is the ID of the hcc (honeycomb core) in the mesh.
	// I need to map copyNo --> x,y (center of a given hcc).
	// First find the number of hcc in x and y for this spiral
	//G4int Nx = GetNElementsSpiral(0, __spiral_x);
	//G4int Ny = GetNElementsSpiral(0, __spiral_y);

	//G4double xposition = fStartZ + copyNo * fSpacing;


}

void RXGCollimatorParameterization::ComputeDimensions
(G4Polyhedra & hcc, const G4int copyNo, const G4VPhysicalVolume*) const
{



}

G4int RXGCollimatorParameterization::GetNElementsSpiral(G4int idx, spiral_axis axis) {

	// The formula for x and y is --> N_{i} = N_{i-1} + 1
	// Where for spiral0 	Nx = __spiral0_nelements_x = 3
	//  and  				Ny = __spiral0_nelements_y = 2

	G4int N = __spiral0_nelements_x;

	if( axis == __spiral_y ) N = __spiral0_nelements_y;

	if ( idx == 0 ) return N;

	// If idx is > 0
	for ( int i = 1 ; i <= idx ; i++ ) {
		N = 2*N + 1;
	}

	return N;
}

