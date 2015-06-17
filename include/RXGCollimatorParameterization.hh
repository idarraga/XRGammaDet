/*
 * RXGCollimatorParameterization.hh
 *
 *  Created on: Jan 14, 2015
 *      Author: idarraga
 */

#ifndef RXGCOLLIMATORPARAMETERIZATION_HH_
#define RXGCOLLIMATORPARAMETERIZATION_HH_

#include "globals.hh"
#include "G4VPVParameterisation.hh"

#define __spiral0_nelements_x 3
#define __spiral0_nelements_y 2

typedef enum {
	__spiral_x = 0,
	__spiral_y
} spiral_axis;

class RXGCollimatorParameterization : public G4VPVParameterisation {

public:
	RXGCollimatorParameterization(
			G4int nSpiral,
			G4int nCoresX,
			G4int nCoresY,
			G4double spacing
	);
	~RXGCollimatorParameterization();

	void ComputeTransformation (const G4int copyNo,
			G4VPhysicalVolume* physVol) const;

	void ComputeDimensions (G4Polyhedra & hcc,const G4int copyNo,
			const G4VPhysicalVolume * physVol) const;

	G4int GetNElementsSpiral(G4int idx, spiral_axis);

private:  // Dummy declarations to get rid of warnings ...

	void ComputeDimensions (G4Tubs&, const G4int,
			const G4VPhysicalVolume*) const {}
	void ComputeDimensions (G4Box&,const G4int,
			const G4VPhysicalVolume*) const {}
	void ComputeDimensions (G4Trd&,const G4int,
			const G4VPhysicalVolume*) const {}
	void ComputeDimensions (G4Trap&,const G4int,
			const G4VPhysicalVolume*) const {}
	void ComputeDimensions (G4Cons&,const G4int,
			const G4VPhysicalVolume*) const {}
	void ComputeDimensions (G4Sphere&,const G4int,
			const G4VPhysicalVolume*) const {}
	void ComputeDimensions (G4Orb&,const G4int,
			const G4VPhysicalVolume*) const {}
	void ComputeDimensions (G4Ellipsoid&,const G4int,
			const G4VPhysicalVolume*) const {}
	void ComputeDimensions (G4Torus&,const G4int,
			const G4VPhysicalVolume*) const {}
	void ComputeDimensions (G4Para&,const G4int,
			const G4VPhysicalVolume*) const {}
	void ComputeDimensions (G4Hype&,const G4int,
			const G4VPhysicalVolume*) const {}
	void ComputeDimensions (G4Polycone&,const G4int,
			const G4VPhysicalVolume*) const {}

private:

	G4int		_nSpiral;
	G4int    	_nCoresX;
	G4int    	_nCoresY;
	G4double 	_spacing;

};

#endif /* RXGCOLLIMATORPARAMETERIZATION_HH_ */
