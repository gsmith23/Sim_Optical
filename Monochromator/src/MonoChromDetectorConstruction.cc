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
/// \file MonoChrom/src/MonoChromDetectorConstruction.cc
/// \brief Implementation of the MonoChromDetectorConstruction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "MonoChromDetectorConstruction.hh"

#include "G4NistManager.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MonoChromDetectorConstruction::MonoChromDetectorConstruction()
 : G4VUserDetectorConstruction()
{
  // half lengths for boxes
  fExpHall_x = fExpHall_y = fExpHall_z = 0.55*m;
  fDarkBox_x = fDarkBox_y = fDarkBox_z = 0.5*m;

  fTyvek_x   = fTyvek_y   = 0.25*m;
  // http://protectiontechnologies.dupont.com/tyvek-graphics-na-product-selector
  // 1082D 257 um thick, 105 gsm  
  fTyvek_z   = 0.01375*cm; 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MonoChromDetectorConstruction::~MonoChromDetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* MonoChromDetectorConstruction::Construct()
{

  // --------------------------------------
  // ------------- G4Material -------------
  
  G4NistManager* nist = G4NistManager::Instance();
  
  //-------------
  // Air
  //
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
  
  //-------------
  // Polyethelene
  //  /material/nist/listMaterials

  G4Material* polyEth = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  
  // 1082D 257 um thick, 105 gsm
  // 275 / 105 = 0.382
  
  // http://protectiontechnologies.dupont.com/tyvek-graphics-na-product-selector
  //  G4double tyvek_thickness = 0.275 * m;
  //  G4double tyvek_gsm       = 105.0 * g / m2;
  G4double tyvek_density   = 0.382 * g / cm3;   
  
  G4Material* tyvek = 
    nist->BuildMaterialWithNewDensity("TYVEK",
				      "G4_POLYETHYLENE",
				      tyvek_density);
  //-------------
  
  // array elements for optics
  static const G4int nElements = 2;
  G4double photon_Energies[] = {1.7712027*eV,6.1992093*eV}; //[700,200] nm
  
  //----------------------------------------------------
  
  G4double air_R_Indices[] = {1.0,1.0};  
  G4double pe_R_Indices[]  = {1.5,1.5};
  G4double absorption[]    = {1.*mm,1.*mm};
  
  G4MaterialPropertiesTable* air_Table = new G4MaterialPropertiesTable();
  air_Table->AddProperty("RINDEX", photon_Energies, air_R_Indices, nElements);
  
  G4MaterialPropertiesTable* pe_Table = new G4MaterialPropertiesTable();
  pe_Table->AddProperty("RINDEX",photon_Energies, pe_R_Indices,nElements);
  pe_Table->AddProperty("ABSLENGTH",photon_Energies,absorption,nElements);

  G4cout << "\n Air G4MaterialPropertiesTable" << G4endl;
  air_Table->DumpTable();
  
  G4cout << "\n PE G4MaterialPropertiesTable" << G4endl;
  pe_Table->DumpTable();
  
  // Assign materials with refractive indices 
  air->SetMaterialPropertiesTable(air_Table);
  tyvek->SetMaterialPropertiesTable(pe_Table);
  polyEth->SetMaterialPropertiesTable(pe_Table);
  
  // --------------------------------------
  // ------------- Volumes ----------------
  

  
  // --------------------------------------
  // The Experimental Hall
  //
  
  // Solid
  G4Box* expHall_box = new G4Box("World",fExpHall_x,fExpHall_y,fExpHall_z);
  
  // Logical
  G4LogicalVolume* expHall_log
    = new G4LogicalVolume(expHall_box,air,"World",0,0,0);
  
  // Physical
  G4VPhysicalVolume* expHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0);
  
  // --------------------------------------
  // The Dark Box
  //
  G4Box* darkBox_box = new G4Box("DarkBox",fDarkBox_x,fDarkBox_y,fDarkBox_z);
  
  G4LogicalVolume* darkBox_log
    = new G4LogicalVolume(darkBox_box,air,"DarkBox",0,0,0);
  
  G4VPhysicalVolume* darkBox_phys
    = new G4PVPlacement(0,G4ThreeVector(),darkBox_log,"DarkBox",
                        expHall_log,false,0);
  
  // --------------------------------------
  // The Tyvek Panel
  //
  G4Box* tyvekSheet_box = new G4Box("Tyvek",fTyvek_x,fTyvek_y,fTyvek_z);
  
  G4LogicalVolume* tyvekSheet_log
    = new G4LogicalVolume(tyvekSheet_box,
			  tyvek,"Tyvek",0,
			  0,
			  0);
  
  G4RotationMatrix* yRot = new G4RotationMatrix;
  yRot->rotateY(-M_PI/4.*rad);
  
  G4VPhysicalVolume* tyvekSheet_phys = nullptr;
  tyvekSheet_phys = new G4PVPlacement(yRot,G4ThreeVector(0,0,0),
				      tyvekSheet_log,"Tyvek",
				      darkBox_log,false,0); 
  
  // --------------------------------------
  // Photodiode
  G4double  pRMin = 0.;
  G4double  pRMax = 5.*cm;
  G4double  pDz   = 0.1*cm;
  G4double  pSPhi = 0.0;
  G4double  pDPhi = 2*M_PI;
  
  G4Tubs * detector = new G4Tubs("Detector",
				 pRMin,pRMax,pDz,
				 pSPhi,pDPhi);
  
  G4LogicalVolume* detector_log
    = new G4LogicalVolume(detector,
			  polyEth,
			  "Detector",
			  0,0,0);

  G4RotationMatrix* yRot2 = new G4RotationMatrix;
  yRot2->rotateY(-M_PI*3/8.*rad);
  
  G4VPhysicalVolume* detector_phys = nullptr;
  detector_phys = new G4PVPlacement(yRot2,G4ThreeVector(-5*cm,0,-3*cm),
				    detector_log,"Detector",
				    darkBox_log,false,0); 

  // --------------------------------------
  // ------------- Surfaces --------------
  //
  // DarkBox
  //
  G4OpticalSurface* opDarkBoxSurface = new G4OpticalSurface("DarkBoxSurface");
  //  opDarkBoxSurface->SetType(dielectric_dielectric);
  //  opDarkBoxSurface->SetFinish(ground);
  //  opDarkBoxSurface->SetModel(unified);
  
  // G4cout << G4endl;
//   opDarkBoxSurface->SetType(dielectric_LUTDAVIS);
//   opDarkBoxSurface->SetFinish(Rough_LUT);
//   opDarkBoxSurface->SetModel(DAVIS);

  opDarkBoxSurface->SetType(dielectric_metal);
  opDarkBoxSurface->SetFinish(polished);
  opDarkBoxSurface->SetModel(glisur);
  opDarkBoxSurface->SetPolish(0.0); // 0 is max roughness
  
  //
  // Tyvek
  //
  G4OpticalSurface* opTyvekSurface = new G4OpticalSurface("TyvekSurface");
  
  // // below are the defaults anyway
  //   opTyvekSurface->SetType(dielectric_dielectric);
  //   //opTyvekSurface->SetType(dielectric_metal);
  //   opTyvekSurface->SetFinish(polished);
  //   opTyvekSurface->SetModel(glisur);
  //   // 0.0 max roughness - Lambertian
  //   // 1.0     smooth - Snell's 
  //   opTyvekSurface->SetPolish(0.0); 
  
  // UNIFIED MODEL
  //http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/ForApplicationDeveloper/BackupVersions/V10.4/html/_images/UNIFIED_model_diagram.png
     
  opTyvekSurface->SetModel(unified);
  
  
  /*
    dielectric_metal
     polished             - no refraction, if not absorbed then specular spike reflection 
     ground               - no refraction, if reflected then coefficients
  */
  
  // Apply coefficients but no transmission
  // dielectric_metal gives no transmission ie no refraction
  opTyvekSurface->SetType(dielectric_metal);
  opTyvekSurface->SetFinish(ground);
  
  /* 
     dielectric_dielectric
     polished             - if not absorbed then apply Snell's law - ie refract or reflect 
     ground               - as above but if reflected use coefficients
     polishedfrontpainted - no refraction, if not absorbed then specular spike reflection (yes)
     groundfrontpainted   - lambertian reflection or absorption only  
     polishedbackpainted  - out-in: absorption or specular spike, in-out: absorption or coefficients
     groundbackpainted    - out-in: absorption or lambertian, in-out: absorption or coefficients    (ID option)
  */

  // Apply coefficients
  //   opTyvekSurface->SetType(dielectric_dielectric);
  //   opTyvekSurface->SetFinish(polishedbackpainted);
  
  opTyvekSurface->SetSigmaAlpha(0.23); // 0.0 rad means no lobe
  
  G4LogicalBorderSurface* darkBoxSurface =
    new G4LogicalBorderSurface("DarkBoxSurface",
			       darkBox_phys,
			       expHall_phys,
			       opDarkBoxSurface); // opTyvekSurface
 
//   G4LogicalSkinSurface* darkBoxSurface =
//     new G4LogicalSkinSurface("DarkBoxSurface",
// 			     darkBox_phys,
// 			     expHall_phys,
// 			     opDarkBoxSurface); 
  // print info
  G4OpticalSurface* opticalSurface = dynamic_cast <G4OpticalSurface*>
    (darkBoxSurface->GetSurface(darkBox_phys,expHall_phys)->
     GetSurfaceProperty());
  
  if (opticalSurface) {
    G4cout << "\n Dark Box surface properties " << G4endl;
    opticalSurface->DumpInfo();
  }


  //   G4LogicalSkinSurface* tyvekSurface =
  //     new G4LogicalSkinSurface("TyvekSurface",
  // 			     tyvekSheet_log,
  // 			     opTyvekSurface);
  

  if(tyvekSheet_phys){
    G4LogicalBorderSurface* tyvekSurface =
      new G4LogicalBorderSurface("TyvekSurface",
				 darkBox_phys,
				 tyvekSheet_phys,
				 opTyvekSurface);
  
  
    // print info
    opticalSurface = dynamic_cast <G4OpticalSurface*>
      (tyvekSurface->GetSurface(darkBox_phys,tyvekSheet_phys)->
       GetSurfaceProperty());
    //(tyvekSurface->GetSurface(tyvekSheet_log)->GetSurfaceProperty());
    
    if (opticalSurface) {
      G4cout << "\n Tyvek Sheet surface properties " << G4endl;
      opticalSurface->DumpInfo();
    }
    
  }
  //------------------------------------------------------
  //---------- Surfaces G4MaterialPropertiesTables -------
  //
  // Generate & Add Material Properties Table attached to the optical surfaces
  //
  
  //--------------------------
  //  Dark Box Surface (border)
  G4double refractiveIndex[nElements] = {1.5, 1.5};
  // G4double specularLobe[nElements]    = {0.0, 0.0};
  // G4double specularSpike[nElements]   = {0.0, 0.0};
  // G4double backScatter[nElements]     = {0.0, 0.0};  

  G4double darkBox_reflectivity[nElements]  = {0.0, 0.0};
  G4double darkBox_efficiency[nElements]    = {0.0, 0.0}; 
  
  G4MaterialPropertiesTable* darkBoxSurf_Table = new G4MaterialPropertiesTable();

  darkBoxSurf_Table->AddProperty("REFLECTIVITY",photon_Energies,
				 darkBox_reflectivity ,nElements);
  darkBoxSurf_Table->AddProperty("EFFICIENCY",photon_Energies,
				 darkBox_efficiency   ,nElements);
  darkBoxSurf_Table->AddProperty("RINDEX",photon_Energies,
				 refractiveIndex, nElements);
  
  // darkBoxSurf_Table->AddProperty("SPECULARLOBECONSTANT",  photon_Energies, specularLobe,    nElements);
  // darkBoxSurf_Table->AddProperty("SPECULARSPIKECONSTANT", photon_Energies, specularSpike,   nElements);
  // darkBoxSurf_Table->AddProperty("BACKSCATTERCONSTANT",   photon_Energies, backScatter,     nElements);
  
  G4cout << "\n Dark Box Surface G4MaterialPropertiesTable" << G4endl;
  darkBoxSurf_Table->DumpTable();
  G4cout << G4endl;
  
  opDarkBoxSurface->SetMaterialPropertiesTable(darkBoxSurf_Table);
  
  //--------------------------
  //  Tyvek Surface 
  
  G4double tyvek_R[nElements]   = {0.9, 0.9};
  // no transmission for dielectric_metal
  //  G4double tyvek_Eff[nElements] = {1.0, 1.0}; 
  
  // diffuse Lobe = 1 - (Csl + Css + cBS)
  G4double tyvek_Csl[nElements] = {0.2, 0.2}; // spec lobe
  G4double tyvek_Css[nElements] = {0.0, 0.0}; // spec spike
  G4double tyvek_Cbs[nElements] = {0.0, 0.0}; // back spike
  
  G4MaterialPropertiesTable *tyvekSurf_Table = new G4MaterialPropertiesTable();
  
  tyvekSurf_Table->AddProperty("REFLECTIVITY"         ,photon_Energies,
			       tyvek_R ,nElements);
//   tyvekSurf_Table->AddProperty("EFFICIENCY"           ,photon_Energies,
// 			       tyvek_Eff   ,nElements);
  tyvekSurf_Table->AddProperty("RINDEX"               ,photon_Energies,
			       pe_R_Indices ,nElements);
  tyvekSurf_Table->AddProperty("SPECULARLOBECONSTANT" ,photon_Energies,
			       tyvek_Csl ,nElements);
  tyvekSurf_Table->AddProperty("SPECULARSPIKECONSTANT",photon_Energies,
			       tyvek_Css,nElements);
  tyvekSurf_Table->AddProperty("BACKSCATTERCONSTANT"  ,photon_Energies,
			       tyvek_Cbs  ,nElements);
  
  G4cout << "\n Tyvek Surface G4MaterialPropertiesTable" << G4endl;
  tyvekSurf_Table->DumpTable();
  G4cout << G4endl;
  
  opTyvekSurface->SetMaterialPropertiesTable(tyvekSurf_Table);

  // always return the physical World
  return expHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MonoChromDetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  
  // declare photodiode as a MultiFunctionalDetector scorer
  //  
  G4MultiFunctionalDetector* photodiode;
  photodiode = new G4MultiFunctionalDetector("photodiode");
  
  G4SDManager::GetSDMpointer()->AddNewDetector(photodiode);

  G4VPrimitiveScorer* PS = new G4PSEnergyDeposit("edep");
  photodiode->RegisterPrimitive(PS);
  SetSensitiveDetector("Detector",photodiode);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
