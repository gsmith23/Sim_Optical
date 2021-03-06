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
/// \file ACYL/src/ACYLDetectorConstruction.cc
/// \brief Implementation of the ACYLDetectorConstruction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ACYLDetectorConstruction.hh"

#include "G4NistManager.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ACYLDetectorConstruction::ACYLDetectorConstruction()
 : G4VUserDetectorConstruction()
{
  // half lengths for boxes
  fExpHall_x = fExpHall_y = fExpHall_z = 0.55*m;
  fDarkBox_x = fDarkBox_y = fDarkBox_z = 0.5*m;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ACYLDetectorConstruction::~ACYLDetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ACYLDetectorConstruction::Construct()
{

  // --------------------------------------
  // ------------- G4Material -------------
  
  G4NistManager* nist = G4NistManager::Instance();
  
   G4double a, z, density;
   //G4int nelements;
  
  //-------------
  // Air
  //
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
  //-------------
  
  //-------------
  // glass 
  G4Element * H = new G4Element("H", "H", z=1., a=1.01*g/mole);
  G4Element * C = new G4Element("C", "C", z=6., a=12.01*g/mole);
  
  G4Material* glass;
  glass = new G4Material("Glass", density=1.032*g/cm3,2);
  glass->AddElement(C,91.533*perCent);
  glass->AddElement(H,8.467*perCent);
  
  //-------------
  
  // array elements for optics
  static const G4int nElements = 2;
  G4double photon_Energies[] = {1.7712027*eV,6.1992093*eV}; //[700,200] nm
  G4double air_R_Indices[]   = {1.0,1.0};
  G4double glass_R_Indices[] = {1.0,2.0}; // exagerated for maximum effect
  
  G4MaterialPropertiesTable* air_Table = new G4MaterialPropertiesTable();
  air_Table->AddProperty("RINDEX", photon_Energies, air_R_Indices, nElements);
  
  G4MaterialPropertiesTable* glass_Table = new G4MaterialPropertiesTable();
  glass_Table->AddProperty("RINDEX", photon_Energies, glass_R_Indices, nElements);
  
  G4cout << "\n Air G4MaterialPropertiesTable" << G4endl;
  air_Table->DumpTable();

  G4cout << "\n Glass G4MaterialPropertiesTable" << G4endl;
  glass_Table->DumpTable();
  
  air->SetMaterialPropertiesTable(air_Table);
  glass->SetMaterialPropertiesTable(glass_Table);

  // --------------------------------------
  // ------------- Volumes ----------------
  
  // The Experimental Hall
  //
  G4Box* expHall_box = new G4Box("World",fExpHall_x,fExpHall_y,fExpHall_z);
  
  G4LogicalVolume* expHall_log
    = new G4LogicalVolume(expHall_box,air,"World",0,0,0);
  
  G4VPhysicalVolume* expHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0);
  
  // The Dark Box
  //
  G4Box* darkBox_box = new G4Box("DarkBox",fDarkBox_x,fDarkBox_y,fDarkBox_z);
  
  G4LogicalVolume* darkBox_log
    = new G4LogicalVolume(darkBox_box,air,"DarkBox",0,0,0);
  
  G4VPhysicalVolume* darkBox_phys
    = new G4PVPlacement(0,G4ThreeVector(),darkBox_log,"DarkBox",
                        expHall_log,false,0);
  

  // Prism
  //  http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/ForApplicationDeveloper/BackupVersions/V10.4/html/Detector/Geometry/geomSolids.html

  G4double prism_x1 = 0.25*m;
  G4double prism_x2 = 0.00*m;
  G4double prism_y1 = 0.25*m;
  G4double prism_y2 = 0.25*m;
  G4double prism_z  = 0.25*m;
  
  G4Trd * prism = new G4Trd("Prism",
			    prism_x1,
			    prism_x2,
			    prism_y1,
			    prism_y2,
			    prism_z);
  
  G4LogicalVolume* prism_log
    = new G4LogicalVolume(prism,
			  glass,
			  "Prism",
			  0,0,0);

  G4RotationMatrix* yRot = new G4RotationMatrix;
  yRot->rotateY(-M_PI/2.*rad);
  
  G4VPhysicalVolume* prism_phys =
    new G4PVPlacement(yRot,G4ThreeVector(0,0,0),
		      prism_log,"Prism",
		      darkBox_log,false,0);
  
  
  // ------------- Surfaces --------------
  //
  // DarkBox
  //
  G4OpticalSurface* opDarkBoxSurface = new G4OpticalSurface("DarkBoxSurface");
  opDarkBoxSurface->SetType(dielectric_metal);
  opDarkBoxSurface->SetFinish(polished);
  opDarkBoxSurface->SetModel(glisur);
  opDarkBoxSurface->SetPolish(0.0); // 0 is max roughness
  
  G4LogicalBorderSurface* darkBoxSurface =
    new G4LogicalBorderSurface("DarkBoxSurface",
			       darkBox_phys,
			       expHall_phys,
			       opDarkBoxSurface); // opTyvekSurface
 
  // print info
  G4OpticalSurface* opticalSurface = dynamic_cast <G4OpticalSurface*>
    (darkBoxSurface->GetSurface(darkBox_phys,expHall_phys)->
     GetSurfaceProperty());
  
  if (opticalSurface) {
    G4cout << "\n Dark Box surface properties " << G4endl;
    opticalSurface->DumpInfo();
  }


    

  //------------------------------------------------------
  //---------- Surfaces G4MaterialPropertiesTables -------
  //
  // Generate & Add Material Properties Table attached to the optical surfaces
  //
  
  //--------------------------
  //  Dark Box Surface (border)
  G4double refractiveIndex[nElements] = {1.5, 1.5};
  G4double darkBox_reflectivity[nElements]  = {0.0, 0.0};
  G4double darkBox_efficiency[nElements]    = {0.0, 0.0}; 
  
  G4MaterialPropertiesTable* darkBoxSurf_Table = new G4MaterialPropertiesTable();

  darkBoxSurf_Table->AddProperty("REFLECTIVITY",photon_Energies,
				 darkBox_reflectivity ,nElements);
  darkBoxSurf_Table->AddProperty("EFFICIENCY",photon_Energies,
				 darkBox_efficiency   ,nElements);
  darkBoxSurf_Table->AddProperty("RINDEX",photon_Energies,
				 refractiveIndex, nElements);
  
  G4cout << "\n Dark Box Surface G4MaterialPropertiesTable" << G4endl;
  darkBoxSurf_Table->DumpTable();
  G4cout << G4endl;
  
  opDarkBoxSurface->SetMaterialPropertiesTable(darkBoxSurf_Table);
  
  // always return the physical World
  return expHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
