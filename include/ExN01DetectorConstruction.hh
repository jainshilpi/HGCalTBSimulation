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
//
// $Id: ExN01DetectorConstruction.hh,v 1.6 2006-06-29 17:47:13 gunter Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//

#include "G4SystemOfUnits.hh"
#ifndef ExN01DetectorConstruction_H
#define ExN01DetectorConstruction_H 1

#include "TString.h"

class G4LogicalVolume;
class G4VPhysicalVolume;

const double x0_pb = 0.0056*CLHEP::m;
//const double x0_pb = 0.003504*CLHEP::m;

#include "G4VUserDetectorConstruction.hh"

class ExN01DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  ExN01DetectorConstruction();
  ExN01DetectorConstruction(TString config);
  ~ExN01DetectorConstruction();
  
  G4VPhysicalVolume* Construct();
  G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
  int    nsensorLayer;
  G4double abs_hz;
  G4double abs_hx;
  G4double abs_hy;
  G4double sensorWidth;
  G4double sensorThickness[18];
  double thickness_cover;
  
  static ExN01DetectorConstruction* Instance () { return fInstance ; } ;
  static ExN01DetectorConstruction* fInstance ;
  
private:
  
  // Logical volumes
  //
  G4LogicalVolume* worldLog;
  G4LogicalVolume* absLog;
  G4LogicalVolume* detLog;
  G4LogicalVolume* outLog;        
  G4LogicalVolume* mcpLog;        
  
  // Physical volumes
  //
  G4VPhysicalVolume* worldPhys;
  G4VPhysicalVolume* absPhys;
  G4VPhysicalVolume* detPhys;    
  G4VPhysicalVolume* outPhys;        
  G4VPhysicalVolume* mcpPhys;   
  
  G4LogicalVolume*  fScoringVolume;         
  


  
  
};

#endif

