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
// $Id: PrtActionInitialization.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file PrtActionInitialization.cc
/// \brief Implementation of the PrtActionInitialization class

#include "PrtActionInitialization.h"
#include "PrtPrimaryGeneratorAction.h"
#include "PrtRunAction.h"
#include "PrtSteppingAction.h"
#include "PrtStackingAction.h"
#include "PrtSteppingVerbose.h"
#include "PrtTrackingAction.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrtActionInitialization::PrtActionInitialization()
  : G4VUserActionInitialization(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrtActionInitialization::~PrtActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrtActionInitialization::BuildForMaster() const
{
  SetUserAction(new PrtRunAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrtActionInitialization::Build() const
{
  SetUserAction(new PrtPrimaryGeneratorAction());
  SetUserAction(new PrtRunAction());
  SetUserAction(new PrtSteppingAction());
  SetUserAction(new PrtStackingAction());
  SetUserAction(new PrtTrackingAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VSteppingVerbose*
PrtActionInitialization::InitializeSteppingVerbose() const
{
  return new PrtSteppingVerbose();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
