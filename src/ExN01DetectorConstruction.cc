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
// $Id: ExN01DetectorConstruction.cc,v 1.9 2006-06-29 17:47:19 gunter Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//

#include "ExN01DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "ExN01CreateTree.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4NistManager.hh"
//pre-shower = |abs = X0| |mcp| |	abs = 2*X0	| |mcp|
   		
/*
float n = 1;    //absorber thickness
float b = 20;    //silica window thickness
float mcp_z = 0.8; //mcp thickness

//Pb
//double X_0 = 0.0143;       //Cu
double X_0 = 0.0056 ;       //Pb (X_0 = 0.0056, but absorber thickness is 0.004)
G4double abs_hx = 0.3*m;
G4double abs_hy = 0.3*m;
*/


//Cu
/*
double X_0 = 0.0143;       //Cu
G4double abs_hx = 0.1*m;
G4double abs_hy = 0.1*m;
*/

ExN01DetectorConstruction* ExN01DetectorConstruction::fInstance = NULL ;

int verbosity = 1;
///SJ
///world
G4double world_hx = 0.3*m;
G4double world_hy = 0.3*m;
G4double world_hz = 2.*m;

///Absorber

//double x0_pb = 0.0056*CLHEP::m;
//polyethylene sheet




/*G4double sensorThickness[6] = {15*CLHEP::mm, 15*CLHEP::mm,
			       10*CLHEP::mm, 10*CLHEP::mm,
			       5*CLHEP::mm, 5*CLHEP::mm
}; 
*/

//G4double sensorThickness = 300*CLHEP::mm/1000.;





ExN01DetectorConstruction::ExN01DetectorConstruction()
 :  worldLog(0), absLog(0), outLog(0),
    worldPhys(0), absPhys(0), outPhys(0)
{

  if ( fInstance )
    {
      return;
    }

  this->fInstance = this ;


  thickness_cover = 1.5*CLHEP::cm;   
  //G4double abs_hz = 5*CLHEP::cm;    //absorber thickness
  abs_hz = 3*x0_pb;    //absorber thickness
  //abs_hz = 6*x0_pb;    //absorber thickness
  //abs_hz = 30*x0_pb;    //absorber thickness
  abs_hx = 20*CLHEP::cm;
  abs_hy = 20*CLHEP::cm;

  
  ///sensor
  nsensorLayer = 6;
  sensorWidth     = 5*CLHEP::mm;
  //G4double sensorWidth     = 0.3*CLHEP::m;
  sensorThickness[0] = 300*CLHEP::mm/1000.;
  sensorThickness[1] = 300*CLHEP::mm/1000;
  sensorThickness[2] = 300*CLHEP::mm/1000; 
  sensorThickness[3] = 300*CLHEP::mm/1000;
  sensorThickness[4] = 300*CLHEP::mm/1000;
  sensorThickness[5] = 300*CLHEP::mm/1000;
  
  
  
}

ExN01DetectorConstruction::~ExN01DetectorConstruction()
{
}

G4VPhysicalVolume* ExN01DetectorConstruction::Construct()
{

  //------------------------------------------------------ materials

G4String name, symbol;
G4double z, a, density,fractionmass;
G4int ncomponents, natoms;


//defining Pb
density = 11.34*g/cm3;
a = 207.2*g/mole;
G4Element* elPb  = new G4Element(name="Lead"  ,symbol="Pb" , z= 82., a);
G4Material* Pb = new G4Material(name="Lead", z=82., a, density);

//defining Cu
density = 8.96*g/cm3;
a = 63.546*g/mole;
G4Element* elCu  = new G4Element(name="Copper"  ,symbol="Cu" , z= 29., a);
G4Material* Cu = new G4Material(name="Cu", z=29., a, density);

//defining O and O2
density = 0.001429*g/cm3;
a = 15.9994*g/mole;
G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

//defining Si
density = 2.33*g/cm3;
a = 28.0855*g/mole;
G4Element* elSi  = new G4Element(name="Silicon"  ,symbol="Si" , z= 14., a);
G4Material* matSi = new G4Material(name="Silicon",density,ncomponents=1);
matSi->AddElement(elSi, natoms=1);

//defining N and N2
density = 0.001251*g/cm3;
a = 14.0067*g/mole;
G4Element* elN  = new G4Element(name="Nitrogen"  ,symbol="N" , z= 7., a);

//defining air
density = 0.001275*g/cm3;
G4Material* Air = new G4Material(name="Air",density,ncomponents=2);

Air->AddElement(elN, fractionmass=70.*perCent);
Air->AddElement(elO, fractionmass=30.*perCent);

//defining Si02
density = 2.196*g/cm3;
G4Material* SiO2 = new G4Material(name="Silica",density,ncomponents=2);
SiO2->AddElement(elSi, natoms=1);
SiO2->AddElement(elO, natoms=2);
G4MaterialPropertiesTable* SiO2_MPT = new G4MaterialPropertiesTable();
G4double ppckov[9] = {0.62*eV,0.69*eV,0.775*eV,0.8856*eV,1.0332*eV,1.2398*eV,1.5498*eV,2.0664*eV,3.0996*eV};
G4double rindex[9] = {1.52   ,1.524  ,1.529   ,1.531    ,1.532    ,1.536    ,1.539    ,1.543    ,1.558};
SiO2_MPT->AddProperty("RINDEX",ppckov,rindex,9);
SiO2 -> SetMaterialPropertiesTable(SiO2_MPT);

//defining PbO
density = 9.5*g/cm3;
G4Material* PbO = new G4Material(name="Lead Oxide",density,ncomponents=2);
PbO->AddElement(elPb, natoms=1);
PbO->AddElement(elO, natoms=1);
  
//defining lead glass
density = 5.12*g/cm3;
G4Material* LeadGlass = new G4Material(name="LeadGlass",density,ncomponents=2);

//defining silicon material


LeadGlass->AddMaterial(PbO, fractionmass=40.*perCent);
LeadGlass->AddMaterial(SiO2, fractionmass=60.*perCent);  


/////Define polyethylene material
//G4_POLYETHYLENE
  
 G4NistManager* man_mat = G4NistManager::Instance();
 G4Material* c2h4_n  = man_mat->FindOrBuildMaterial("G4_POLYETHYLENE");
  //------------------------------------------------------ experimental apparatus

//defining world environment

G4Box* worldBox = new G4Box("World", 0.5*world_hx, 0.5*world_hy, 0.5*world_hz);
G4LogicalVolume* worldLog = new G4LogicalVolume(worldBox, Air, "World");
G4PVPlacement* worldPhys = new G4PVPlacement(0,G4ThreeVector(),worldLog,"World",0,false,0);
// worldLog->SetVisAttributes (G4VisAttributes::Invisible);

//defining first absorber: X0

//G4double abs_hz = 2*0.5*X_0*m;

G4Box* absBox = new G4Box("Abs", 0.5*abs_hx, 0.5*abs_hy, 0.5*abs_hz);

//G4LogicalVolume* absLog = new G4LogicalVolume(absBox, Cu, "Abs");
G4LogicalVolume* absLog = new G4LogicalVolume(absBox, Pb, "Abs");

 G4double pos_x =  0.0*CLHEP::m;
 G4double pos_y =  0.0*CLHEP::m;
//G4double pos_z =  4.*0.5*n*X_0;//+0.1*meter;
 G4double pos_z =  0.1 * CLHEP::m;
//G4double pos_z =  0.20*meter;

G4PVPlacement* absPhys = new G4PVPlacement(0,G4ThreeVector(pos_x, pos_y, pos_z),absLog,"Abs",worldLog,false,0);


///polyethylene attached to the absorber
//c2h4_n


 G4Box* coverBox = new G4Box("coverBox", 0.5*abs_hx, 0.5*abs_hy, 0.5*thickness_cover);
 G4LogicalVolume* coverLog = new G4LogicalVolume(coverBox, Pb, "cover");
 
 pos_x =  0.0*CLHEP::m;
 pos_y =  0.0*CLHEP::m;
 pos_z =  pos_z + 0.5*abs_hz + 0.5*thickness_cover;
 //G4double pos_z =  0.20*meter;
 
 //for now - remove it
 ///G4PVPlacement* coverPhys = new G4PVPlacement(0,G4ThreeVector(pos_x, pos_y, pos_z),coverLog,"cover",worldLog,false,0);
 
////sensors now - square shaped - place 4 layers of sensor
 //G4Box* solid = new G4Box("ECalSensitive", 0.5*sensorWidth, 0.5*sensorWidth, 0.5*sensorThickness);
 /*G4Box* solid = new G4Box("ECalSensitive", 0.5*sensorWidth, 0.5*sensorWidth, 0.5*sensorThickness[0]);
 G4LogicalVolume *logSens = new G4LogicalVolume(solid, matSi, "ECalSensitive");
 */

 /*G4VisAttributes* atb = new G4VisAttributes(G4Colour(1.0,1.0,0.5));
 atb->SetForceSolid(true);  
 logSens->SetVisAttributes(atb);
 */
 //if (verbosity > 0) G4cout << "Create " << solid << G4endl;

 if (verbosity > 0) std::cout<<"Earlier position of Z "<<pos_z<<std::endl;
 pos_x =  0.0*CLHEP::m;
 pos_y =  0.0*CLHEP::m;
 //pos_z =  pos_z + 0.5*thickness_cover + 0.5*sensorThickness + 7.2*CLHEP::cm; //position of first layer
 pos_z =  pos_z + 0.5*thickness_cover + 0.5*sensorThickness[0] + 7.2*CLHEP::cm; //position of first layer
 if (verbosity > 0) std::cout<<"Now position of Z "<<pos_z<<std::endl;

 G4ThreeVector trans;
 
 for (int i=0; i<nsensorLayer; i++) {
   
   ////sensors now - square shaped - place 4 layers of sensor
   G4Box* solid = new G4Box("ECalSensitive", 0.5*sensorWidth, 0.5*sensorWidth, 0.5*sensorThickness[i]);
   G4LogicalVolume *logSens = new G4LogicalVolume(solid, matSi, "ECalSensitive");
   
   
   if (verbosity > 0) G4cout << "Create " << solid << G4endl;
   

   
   trans = G4ThreeVector(0,0,pos_z);
   
   
   if (verbosity > 0) std::cout<<"For "<<i<<"th layer, Z pos is "<<pos_z<<std::endl;
   //new G4PVPlacement(0,trans,logSens,Form("sensitiveLayer_%d",i),worldLog,false,i+1); ///last argument is for the copy number - useful later
   new G4PVPlacement(0,trans,logSens,"sensitiveLayer",worldLog,false,i+1); ///last argument is for the copy number - useful later
   //pos_z += 4*CLHEP::cm + 2*0.5*sensorThickness;
   
   //pos_z += 4*CLHEP::cm + 2*0.5*sensorThickness[0];
   
   if(i != (nsensorLayer-1)) pos_z += 2*CLHEP::cm + 0.5*sensorThickness[i] + 0.5*sensorThickness[i+1];
 }


 if(verbosity > 0) std::cout<<"Done defining the geometry "<<std::endl;


  //------------------------------------------------------------------

  return worldPhys;
}
