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
// Code developed by:
//  S.Larsson and J. Generowicz.

#include "G4SystemOfUnits.hh"
#include "G4AutoLock.hh"

#include "TVector2.h"
#include "TVector3.h"

#include "PrtManager.h"
#include "PrtRun.h"
#include "PrtField.h"

namespace {
G4Mutex myPrtFieldLock = G4MUTEX_INITIALIZER;
}

using namespace std;

PrtField::PrtField(const char *filename, double zOffset)
  : fZoffset(zOffset), invertX(false), invertY(false), invertZ(false) {

  double lenUnit;
  double fieldUnit;

  //G4cout << "\n-----------------------------------------------------------"
  //       << "\n      Magnetic field"
  //       << "\n-----------------------------------------------------------";
	G4cout<<"PrtField::PrtField -- Reading field map: "<<filename<<G4endl;
	auto fRun	= PrtManager::Instance()->getRun();
	int kField	= fRun->getField();
	if (kField<1 || kField>3){
		G4cout<<"PrtField::PrtField -- Do not understand field parameter! \t field = "<<kField<<G4endl;
		exit(0);
	} else {
		if (kField==1){		// default field.tab is in units of Gauss
			lenUnit		= cm;
			fieldUnit	= gauss;
		} else if (kField==2||kField==3){
			lenUnit		= cm;
			fieldUnit	= tesla;
		}
	}

  // This is a thread-local class and we have to avoid that all workers open the
  // file at the same time
  G4AutoLock lock(&myPrtFieldLock);

  ifstream file(filename); // Open the file for reading.

  if (!file.is_open()) {
    G4ExceptionDescription ed;
    ed << "Could not open input file " << filename << std::endl;
    G4Exception("PrtField::PrtField", "pugmag001", FatalException, ed);
  }
  char buffer[256];

  // Read table dimensions
  file >> nz >> ny >> nx;		 // Note dodgy order
  G4cout << "PrtField::PrtField -- Number of values x,y,z: "<<nx<<" "<<ny<<" "<<nz<<G4endl;

  // Set up storage space for table
  xField.resize(nx);
  yField.resize(nx);
  zField.resize(nx);
  int ix, iy, iz;
  for (ix = 0; ix < nx; ix++) {
    xField[ix].resize(ny);
    yField[ix].resize(ny);
    zField[ix].resize(ny);
    for (iy = 0; iy < ny; iy++) {
      xField[ix][iy].resize(nz);
      yField[ix][iy].resize(nz);
      zField[ix][iy].resize(nz);
    }
  }

	//---- read in the dummy lines at top of field.tab
	if (kField==1){
	  // Ignore other header information
	  // The first line whose second character is '0' is considered to
	  // be the last line of the header.
	  do {
		file.getline(buffer, 256);
	  } while (buffer[1] != '0');
	}	// end kField==1
	
  // Read in the data
  double xval, yval, zval, bx, by, bz;
  for (iy = 0; iy < ny; iy++) {
    for (iz = 0; iz < nz; iz++) {
      file >> xval >> yval >> zval >> bx >> by >> bz;		// xval always zero! bx BETTER always be zero!!
      if (kField==1){
	      if (fabs(bx) < 0.01) bx = 0;		// MISSES SIGN...
	      if (fabs(by) < 0.01) by = 0;		// MISSES SIGN...
	      if (fabs(bz) < 0.01) bz = 0;		// MISSES SIGN...
      }
      if (iy == 0 && iz == 0) {
        minx = xval * lenUnit;
        miny = yval * lenUnit;
        minz = zval * lenUnit;
      }
      xField[0][iy][iz] = bx * fieldUnit;
      yField[0][iy][iz] = by * fieldUnit;
      zField[0][iy][iz] = bz * fieldUnit;
    }
  }
  file.close();
  lock.unlock();
  G4cout<<"PrtField::PrtField -- done reading"<<G4endl;

  maxx = xval * lenUnit;
  maxy = yval * lenUnit;
  maxz = zval * lenUnit;

  // G4cout << " Read values of field from file " << filename << G4endl;
  G4cout << " ---> assumed the order:  x, y, z, Bx, By, Bz "
         << "\n ---> Min values x,y,z: " << minx / cm << " " << miny / cm << " " << minz / cm
         << " cm "
         << "\n ---> Max values x,y,z: " << maxx / cm << " " << maxy / cm << " " << maxz / cm
         << " cm "
         << "\n ---> The field will be offset by " << zOffset / cm << " cm " << G4endl;

  // Should really check that the limits are not the wrong way around.
  if (maxx < minx) {
    swap(maxx, minx);
    invertX = true;
  }
  if (maxy < miny) {
    swap(maxy, miny);
    invertY = true;
  }
  if (maxz < minz) {
    swap(maxz, minz);
    invertZ = true;
  }
  G4cout << "\nAfter reordering if neccesary"
         << "\n ---> Min values x,y,z: " << minx / cm << " " << miny / cm << " " << minz / cm
         << " cm "
         << " \n ---> Max values x,y,z: " << maxx / cm << " " << maxy / cm << " " << maxz / cm
         << " cm ";

  dx = maxx - minx;
  dx = 1;
  dy = maxy - miny;
  dz = maxz - minz;
  G4cout << "\n ---> Dif values x,y,z (range): " << dx / cm << " " << dy / cm << " " << dz / cm
         << " cm in z "
         << "\n-----------------------------------------------------------" << G4endl;
}

//
//	WJL	NOTE default map is only defined in the y-z plane, then 
//	WJL		a rotation about Z is used to get all other values
//
void PrtField::GetFieldValue(const double point[4], double *Bfield) const {

  TVector2 t(point[0],point[1]);	// (x,y)
  double z = point[2];				//  z
  double r = t.Mod();				//  rho

  // Check that the point is within the defined region
  if (r >= miny && r <= maxy && z >= minz && z <= maxz) {

    // Position of given point within region, normalized to the range  [0,1]
    double rfraction = (r - miny) / dy;
    double zfraction = (z - minz) / dz;

    // Need addresses of these to pass to modf below.
    // modf uses its second argument as an OUTPUT argument.
    double rdindex, zdindex;

    // Position of the point within the cuboid defined by the
    // nearest surrounding tabulated points
    std::modf(rfraction * (ny - 1), &rdindex);
    std::modf(zfraction * (nz - 1), &zdindex);

    // The indices of the nearest tabulated point whose coordinates
    // are all less than those of the given point
    int rindex = static_cast<int>(rdindex);
    int zindex = static_cast<int>(zdindex);

    TVector3 vfield(xField[0][rindex][zindex], 
                    yField[0][rindex][zindex],
                    zField[0][rindex][zindex]);

    vfield.RotateZ(t.Phi()-TMath::PiOver2());

    Bfield[0] = vfield.X();
    Bfield[1] = vfield.Y();
    Bfield[2] = vfield.Z();

  } else {
    Bfield[0] = 0.0;
    Bfield[1] = 0.0;
    Bfield[2] = 0.0;
  }
}
