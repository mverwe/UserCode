// -*- C++ -*-
//
// Package:    SiPixelAnalyzer
// Class:      SiPixelAnalyzer
//  
/**\class SiPixelAnalyzer SiPixelAnalyzer.cc Validation/SiPixelAnalyzer/src/SiPixelAnalyzer.cc

 Description: <one line class summary>
  
 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Ivan Amos Cali
//         Created:  Wed Jun  4 18:45:43 CEST 2008
// $Id$
//
//


// system include files
#include <memory>
#include <vector>
#include <iostream>
#include <ctime>        // std::clock()
#include <algorithm>    // std::find()

// user include files 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/SiPixelRawData/interface/SiPixelRawDataError.h"   // for error vs. centrality plots

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

#include "CondFormats/SiPixelObjects/interface/SiPixelFedCablingMap.h"
#include "CondFormats/SiPixelObjects/interface/SiPixelFedCablingTree.h"
#include "CondFormats/DataRecord/interface/SiPixelFedCablingMapRcd.h"
#include "CondFormats/SiPixelObjects/interface/PixelFEDCabling.h"
#include "CondFormats/SiPixelObjects/interface/SiPixelFrameConverter.h"
#include "CondFormats/SiPixelObjects/interface/PixelFEDLink.h"
#include "CondFormats/SiPixelObjects/interface/PixelROC.h"
#include "CondFormats/SiPixelObjects/interface/ElectronicIndex.h"
#include "CondFormats/SiPixelObjects/interface/DetectorIndex.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

//ROOT inclusion
#include "TROOT.h"
#include "TFile.h"
//#include "TNtuple.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace std;
using namespace edm;
//
// class decleration
//

class SiPixelAnalyzer : public edm::EDAnalyzer {
   public:
      explicit SiPixelAnalyzer(const edm::ParameterSet&);
      ~SiPixelAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
  


  // ----------member data ---------------------------
  std::string outputFile_;
  std::string outTreeName_;
  edm::InputTag src_;

  edm::EDGetTokenT<HFRecHitCollection> srcHFhits_;

  unsigned int BarrelColumnsOffset_;
  unsigned int gnHModuleBarrel_;
  unsigned int EndcapColumnsOffset_;
  unsigned int gnHModuleEndcap_;
  const TrackerGeometry* trGeo_;
  TFile* oFile_;
  TFile* oTree_;

  TTree* treeEvent_;
  int run_;
  int event_;
  int lumi_;
  int orbitNumber_;
  int bunchCrossing_;
  double HFRecHitSum_;

  TTree* treePixelDigi_;
  int nPixelDigi_;
  std::vector<int> pixelDigi_id_;
  std::vector<int> pixelDigi_detLayer_;
  std::vector<int> pixelDigi_fedID_;
  std::vector<int> pixelDigi_nHits_B_L1_;
  std::vector<int> pixelDigi_nHits_B_L2_;
  std::vector<int> pixelDigi_nHits_B_L3_;
  std::vector<int> pixelDigi_nHits_EC_D1_;
  std::vector<int> pixelDigi_nHits_EC_D2_;
  std::vector<int> pixelDigi_nHits_B_L123_;
  std::vector<int> pixelDigi_nHits_Total_;
  std::vector<float> moduleGeom_phi_;
  std::vector<float> moduleGeom_eta_;
  std::vector<float> moduleGeom_z_;

  TTree* treePixelRawDataError_;
  int nError_;
  std::vector<int> error_id_;
  std::vector<int> error_detLayer_;
  std::vector<int> error_fedID_;
  std::vector<int> error_type_;

  TTree* treeFED_;
  int nFED_;
  std::vector<int>   FED_id_;
  std::vector<double> FED_occu_B_L1_;
  std::vector<double> FED_occu_B_L2_;
  std::vector<double> FED_occu_B_L3_;
  std::vector<double> FED_occu_EC_D1_;
  std::vector<double> FED_occu_EC_D2_;
  std::vector<double> FED_occu_B_L123_;   // occupancy over Barrel Layer 1,2,3

  // TNtuple* BarrelModuleNt_;
  // TNtuple* BarrelDigisNt_;
  // TNtuple* BarrelColumnsNt_;
  // TNtuple* BarrelDColumnsNt_;
  

  // TNtuple* EndcapModuleNt_;
  // TNtuple* EndcapDigisNt_;
  // TNtuple* EndcapColumnsNt_;
  // TNtuple* EndcapDColumnsNt_;  
    
  // TNtuple* FEDNt_;
  // TNtuple* LinksNt_;
  // TNtuple* ROCOccNt_;
  
  // TH1F* DcolumnHits_; 
  // TH1D* Occupancy_;
  // TH2D* OccupancyZ_;
//  TH1F * BarrelSlinkHitsDistrib_;

  /////
  TH2D* FEDOcc_[6];
  TH2D* FEDErrors_[6];
  TH2D* FEDErrors_Type40_[6];

  TH2D* FEDOccPerFED_[41];

  TH2D* FEDOcc_10xBin_[6];
  TH2D* FEDErrors_10xBin_[6];
  TH2D* FEDErrors_Type40_10xBin_[6];

  TH1D* numErrorsPerEvent_[7];
  TH1D* numErrorsPerEvent_10xRange_[7];
  TH2D* numErrorsPerEvent_vs_HF_[7];
  TH2D* numErrorsPerEvent_vs_HF_10xBin_[7];

  TH1D* numErrorsPerEvent_no3740_[7];
  TH1D* numErrorsPerEvent_no3740_10xRange_[7];
  TH2D* numErrorsPerEvent_no3740_vs_HF_[7];
  TH2D* numErrorsPerEvent_no3740_vs_HF_10xBin_[7];

  //FED errors vs HFE for each FED seperate
  TH2D *hFEDerrorHFE[41];
  
  TH1D* HF_;
  TH1D* HF_10xBin_;
  /////

  //  TH1D* LinkByLinkOcc_;  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SiPixelAnalyzer::SiPixelAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed 
   outputFile_ = iConfig.getUntrackedParameter<string>("outputFile", "pixeldigihisto.root");
   outTreeName_ = iConfig.getUntrackedParameter<string>("outTreeName", "pixeldigitree.root");
   src_ =  iConfig.getParameter<edm::InputTag>( "src" );
   
   srcHFhits_ = consumes<HFRecHitCollection>(iConfig.getParameter<edm::InputTag>("srcHFhits"));

    
     //FEDs studies---------------------------------------------------------------------
//     FEDNt_ = new TNtuple("FED", "FED","fedid:links:roc",100000);
//     LinksNt_= new TNtuple("Links", "Links","fedid:linkn:nHits",100000);
//     ROCOccNt_= new TNtuple("ROCOcc", "ROCOcc","idH:idL:ROCHits:SubDet:gHModule:ROCn:ModuleHits:nROCs",100000);
     
     //Barrel Ntuples----------------------------------------------------------------------------
//     BarrelModuleNt_ = new TNtuple("BarrelPixelDistrib", "BarrelPixelDistri","idH:idL:layer:ladder:ring:z:R:phi:eta:Mhits:HMod0Hits:HMod1Hits:ncolumns:nrows",100000);
//     BarrelDigisNt_ = new TNtuple("BarrelDigis", "BarrelPixelDigisDist","idH:idL:layer:ladder:ring:adc:column:row",100000);
//     BarrelColumnsNt_ = new TNtuple("BarrelColumnOcc", "BarrelColumnOcc","idH:idL:layer:ladder:ring:gHmodule:column:colHits:gcolumn",100000);
//     BarrelDColumnsNt_ = new TNtuple("BarrelDColumnOcc", "BarrelDColumnOcc","idH:idL:layer:ladder:ring:gHmodule:Dcolumn:DcolHits:gDcolumn",100000);
   
     //histos---------------------------------------------------------
//     BarrelSlinkHitsDistrib_ = new TH1F("BarrelSlinkHitsDistrib", "Slink Hits -Barrel-", 801, 0., 800.);
//     BarrelSlinkHitsDistrib_->GetYaxis()->SetTitle("Entries");
//     BarrelSlinkHitsDistrib_->GetXaxis()->SetTitle("Hits per RO Link");
//     BarrelSlinkHitsDistrib_->SetDirectory(oFile_->GetDirectory(0));
 
     /////
     oTree_ = new TFile((const char*)outTreeName_.c_str(), "RECREATE");
   
     treeEvent_ = new TTree("Event", "event information");
     treeEvent_->Branch("run",&run_);
     treeEvent_->Branch("event",&event_);
     treeEvent_->Branch("lumi",&lumi_);
     treeEvent_->Branch("orbitNumber",&orbitNumber_);
     treeEvent_->Branch("bunchCrossing",&bunchCrossing_);
     treeEvent_->Branch("HFRecHitSum",&HFRecHitSum_);

     treePixelDigi_ = new TTree("PixelDigi", "PixelDigi information");
     treePixelDigi_->Branch("nPixelDigi",&nPixelDigi_);
     treePixelDigi_->Branch("pixelDigi_id",&pixelDigi_id_);
     treePixelDigi_->Branch("pixelDigi_detLayer",&pixelDigi_detLayer_);
     treePixelDigi_->Branch("pixelDigi_fedID",&pixelDigi_fedID_);
     treePixelDigi_->Branch("pixelDigi_nHits_L1",&pixelDigi_nHits_B_L1_);
     treePixelDigi_->Branch("pixelDigi_nHits_L2",&pixelDigi_nHits_B_L2_);
     treePixelDigi_->Branch("pixelDigi_nHits_L3",&pixelDigi_nHits_B_L3_);
     treePixelDigi_->Branch("pixelDigi_nHits_EC_D1",&pixelDigi_nHits_EC_D1_);
     treePixelDigi_->Branch("pixelDigi_nHits_EC_D2",&pixelDigi_nHits_EC_D2_);
     treePixelDigi_->Branch("pixelDigi_nHits_B_L123",&pixelDigi_nHits_B_L123_);
     treePixelDigi_->Branch("pixelDigi_nHits_Total",&pixelDigi_nHits_Total_);
     treePixelDigi_->Branch("moduleGeom_phi",&moduleGeom_phi_);
     treePixelDigi_->Branch("moduleGeom_eta",&moduleGeom_eta_);
     treePixelDigi_->Branch("moduleGeom_z",&moduleGeom_z_);

     treePixelRawDataError_ = new TTree("PixelRawDataError", "PixelRawDataError information");
     treePixelRawDataError_->Branch("nError",&nError_);
     treePixelRawDataError_->Branch("error_id",&error_id_);
     treePixelRawDataError_->Branch("error_detLayer",&error_detLayer_);
     treePixelRawDataError_->Branch("error_fedID",&error_fedID_);
     treePixelRawDataError_->Branch("error_type",&error_type_);

     treeFED_ = new TTree("FED", "FED information");
     treeFED_->Branch("nFED",&nFED_);
     treeFED_->Branch("FED_id",&FED_id_);
     treeFED_->Branch("FED_occu_B_L1",&FED_occu_B_L1_);
     treeFED_->Branch("FED_occu_B_L2",&FED_occu_B_L2_);
     treeFED_->Branch("FED_occu_B_L3",&FED_occu_B_L3_);
     treeFED_->Branch("FED_occu_EC_D1",&FED_occu_EC_D1_);
     treeFED_->Branch("FED_occu_EC_D2",&FED_occu_EC_D2_);
     treeFED_->Branch("FED_occu_B_L123",&FED_occu_B_L123_);

     //histo inits
     oFile_ = new TFile((const char*)outputFile_.c_str(), "RECREATE");
     
     for(int i = 0; i<6; i++)
     {    
       if(i<3)      FEDOcc_[i] = new TH2D(Form("FEDOccupancy_B_L%d",i+1),Form("FED Occupancy (Barrel Layer%d);HF Energy Sum;Average FED Occupancy",i+1),40,0,300000,100,0,3);
       else if(i<5) FEDOcc_[i] = new TH2D(Form("FEDOccupancy_EC_D%d",i-2),Form("FED Occupancy (Endcap Disk%d);HF Energy Sum;Average FED Occupancy",i-2),40,0,300000,100,0,3);
       else         FEDOcc_[i] = new TH2D(Form("FEDOccupancy_B_L123"),Form("FED Occupancy (Barrel Layers : 1, 2, 3);HF Energy Sum;Average FED Occupancy"),40,0,300000,100,0,3);
       FEDOcc_[i]->SetDirectory(oFile_->GetDirectory(0));
     }

     //     https://github.com/cms-sw/cmssw/blob/CMSSW_7_5_X/DQM/SiPixelMonitorRawData/doc/HistogramMeanings.doc#L3-L17
     //     errorType - a number (25-38) indicating the type of error recorded.
     //         25 indicates an invalid ROC of 25
     //         26 indicates a gap word
     //         27 indicates a dummy word
     //         28 indicates a FIFO full error
     //         29 indicates a timeout error
     //         30 indicates a TBM error trailer
     //         31 indicates an event number error (TBM and FED event number mismatch)
     //         32 indicates an incorrectly formatted Slink Header
     //         33 indicates an incorrectly formatted Slink Trailer
     //         34 indicates the event size encoded in the Slink Trailer is different than the size found at raw to digi conversion
     //         35 indicates an invalid FED channel number
     //         36 indicates an invalid ROC value
     //         37 indicates an invalid dcol or pixel value
     //         38 indicates the pixels on a ROC weren't read out from lowest to highest row and dcol value
     int errorBin1 = 25;
     int errorBin2 = 38+1;
     int errorBin2_Type40 = 40+1;
     int numErrorBins = errorBin2-errorBin1;
     int numErrorBins_Type40 = errorBin2_Type40-errorBin1;
//     di->errorMessage_ = Error: Unknown error type
//     di->getFedId()    = 37
//     di->getType()     = 40
     for(int i = 0; i<6; i++)
     {
       if(i<3)      FEDErrors_[i] = new TH2D(Form("FEDErrors_B_L%d",i+1) ,Form("SiPixel Errors (Barrel Layer%d);HF Energy Sum;errorType",i+1),40,0,300000,    numErrorBins, errorBin1, errorBin2);
       else if(i<5) FEDErrors_[i] = new TH2D(Form("FEDErrors_EC_D%d",i-2),Form("SiPixel Errors (Endcap Disk%d);HF Energy Sum;errorType",i-2),40,0,300000,     numErrorBins, errorBin1, errorBin2);
       else         FEDErrors_[i] = new TH2D(Form("FEDErrors_B_L123"),    Form("SiPixel Errors (Barrel Layers : 1, 2, 3);HF Energy Sum;errorType"),40,0,300000, numErrorBins, errorBin1, errorBin2);
       FEDErrors_[i]->SetDirectory(oFile_->GetDirectory(0));
     }
     for(int i = 0; i<6; i++)
     {
       if(i<3)      FEDErrors_Type40_[i] = new TH2D(Form("FEDErrors_Type40_B_L%d",i+1) ,Form("SiPixel Errors (Barrel Layer%d);HF Energy Sum;errorType",i+1),40,0,300000,    numErrorBins_Type40, errorBin1, errorBin2_Type40);
       else if(i<5) FEDErrors_Type40_[i] = new TH2D(Form("FEDErrors_Type40_EC_D%d",i-2),Form("SiPixel Errors (Endcap Disk%d);HF Energy Sum;errorType",i-2),40,0,300000,     numErrorBins_Type40, errorBin1, errorBin2_Type40);
       else         FEDErrors_Type40_[i] = new TH2D(Form("FEDErrors_Type40_B_L123"),    Form("SiPixel Errors (Barrel Layers : 1, 2, 3);HF Energy Sum;errorType"),40,0,300000, numErrorBins_Type40, errorBin1, errorBin2_Type40);
       FEDErrors_Type40_[i]->SetDirectory(oFile_->GetDirectory(0));
     }

     // same histograms with 10x bins
     for(int i = 0; i<6; i++)
     {
       if(i<3)      FEDOcc_10xBin_[i] = new TH2D(Form("FEDOccupancy_10xBin_B_L%d",i+1),Form("FED Occupancy (Barrel Layer%d);HF Energy Sum;Average FED Occupancy",i+1),40*10,0,300000,100*10,0,3);
       else if(i<5) FEDOcc_10xBin_[i] = new TH2D(Form("FEDOccupancy_10xBin_EC_D%d",i-2),Form("FED Occupancy (Endcap Disk%d);HF Energy Sum;Average FED Occupancy",i-2),40*10,0,300000,100*10,0,3);
       else         FEDOcc_10xBin_[i] = new TH2D(Form("FEDOccupancy_10xBin_B_L123"),Form("FED Occupancy (Barrel Layers : 1, 2, 3);HF Energy Sum;Average FED Occupancy"),40*10,0,300000,100*10,0,3);
       FEDOcc_10xBin_[i]->SetDirectory(oFile_->GetDirectory(0));
     }
     for(int i = 0; i<6; i++)
     {
       if(i<3)      FEDErrors_10xBin_[i] = new TH2D(Form("FEDErrors_10xBin_B_L%d",i+1) ,Form("SiPixel Errors (Barrel Layer%d);HF Energy Sum;errorType",i+1),40*10,0,300000,    numErrorBins, errorBin1, errorBin2);
       else if(i<5) FEDErrors_10xBin_[i] = new TH2D(Form("FEDErrors_10xBin_EC_D%d",i-2),Form("SiPixel Errors (Endcap Disk%d);HF Energy Sum;errorType",i-2),40*10,0,300000,     numErrorBins, errorBin1, errorBin2);
       else         FEDErrors_10xBin_[i] = new TH2D(Form("FEDErrors_10xBin_B_L123"),    Form("SiPixel Errors (Barrel Layers : 1, 2, 3);HF Energy Sum;errorType"),40*10,0,300000, numErrorBins, errorBin1, errorBin2);
       FEDErrors_10xBin_[i]->SetDirectory(oFile_->GetDirectory(0));
     }
     for(int i = 0; i<6; i++)
     {
       if(i<3)      FEDErrors_Type40_10xBin_[i] = new TH2D(Form("FEDErrors_Type40_10xBin_B_L%d",i+1) ,Form("SiPixel Errors (Barrel Layer%d);HF Energy Sum;errorType",i+1),40*10,0,300000,    numErrorBins_Type40, errorBin1, errorBin2_Type40);
       else if(i<5) FEDErrors_Type40_10xBin_[i] = new TH2D(Form("FEDErrors_Type40_10xBin_EC_D%d",i-2),Form("SiPixel Errors (Endcap Disk%d);HF Energy Sum;errorType",i-2),40*10,0,300000,     numErrorBins_Type40, errorBin1, errorBin2_Type40);
       else         FEDErrors_Type40_10xBin_[i] = new TH2D(Form("FEDErrors_Type40_10xBin_B_L123"),    Form("SiPixel Errors (Barrel Layers : 1, 2, 3);HF Energy Sum;errorType"),40*10,0,300000, numErrorBins_Type40, errorBin1, errorBin2_Type40);
       FEDErrors_Type40_10xBin_[i]->SetDirectory(oFile_->GetDirectory(0));
     }

     for(int i = 0; i<7; i++)
     {
       if(i<3)      numErrorsPerEvent_[i] = new TH1D(Form("numErrorsPerEvent_B_L%d",i+1) ,Form("number of SiPixel Errors (Barrel Layer%d);nSiPixelError",i+1), 50, 0, 50);
       else if(i<5) numErrorsPerEvent_[i] = new TH1D(Form("numErrorsPerEvent_EC_D%d",i-2),Form("number of SiPixel Errors (Endcap Disk%d);nSiPixelError",i-2),  50, 0, 50);
       else if(i<6) numErrorsPerEvent_[i] = new TH1D(Form("numErrorsPerEvent_B_L123"),    Form("number of SiPixel Errors (Barrel Layers : 1, 2, 3);nSiPixelError"), 50, 0, 50);
       else         numErrorsPerEvent_[i] = new TH1D(Form("numErrorsPerEvent_"),          Form("number of SiPixel Errors (Barrel Layers : 1, 2, 3);nSiPixelError"), 50, 0, 50);
       numErrorsPerEvent_[i]->SetDirectory(oFile_->GetDirectory(0));
     }
     for(int i = 0; i<7; i++)
     {
       if(i<3)      numErrorsPerEvent_10xRange_[i] = new TH1D(Form("numErrorsPerEvent_10xRange_B_L%d",i+1) ,Form("number of SiPixel Errors (Barrel Layer%d);nSiPixelError",i+1), 100, 0, 100);
       else if(i<5) numErrorsPerEvent_10xRange_[i] = new TH1D(Form("numErrorsPerEvent_10xRange_EC_D%d",i-2),Form("number of SiPixel Errors (Endcap Disk%d);nSiPixelError",i-2),  100, 0, 100);
       else if(i<6) numErrorsPerEvent_10xRange_[i] = new TH1D(Form("numErrorsPerEvent_10xRange_B_L123"),    Form("number of SiPixel Errors (Barrel Layers : 1, 2, 3);nSiPixelError"), 100, 0, 100);
       else         numErrorsPerEvent_10xRange_[i] = new TH1D(Form("numErrorsPerEvent_10xRange_"),          Form("number of SiPixel Errors (Barrel Layers : 1, 2, 3);nSiPixelError"), 100, 0, 100);
       numErrorsPerEvent_10xRange_[i]->SetDirectory(oFile_->GetDirectory(0));
     }
     for(int i = 0; i<7; i++)
     {
       if(i<3)      numErrorsPerEvent_vs_HF_[i] = new TH2D(Form("numErrorsPerEvent_vs_HF_B_L%d",i+1) ,Form("number of SiPixel Errors (Barrel Layer%d);HF Energy Sum;nSiPixelError",i+1),40,0,300000, 50, 0, 50);
       else if(i<5) numErrorsPerEvent_vs_HF_[i] = new TH2D(Form("numErrorsPerEvent_vs_HF_EC_D%d",i-2),Form("number of SiPixel Errors (Endcap Disk%d);HF Energy Sum;nSiPixelError",i-2),40,0,300000,  50, 0, 50);
       else if(i<6) numErrorsPerEvent_vs_HF_[i] = new TH2D(Form("numErrorsPerEvent_vs_HF_B_L123"),    Form("number of SiPixel Errors (Barrel Layers : 1, 2, 3);HF Energy Sum;nSiPixelError"),40,0,300000, 50, 0, 50);
       else         numErrorsPerEvent_vs_HF_[i] = new TH2D(Form("numErrorsPerEvent_vs_HF_"),          Form("number of SiPixel Errors (Barrel Layers : 1, 2, 3);HF Energy Sum;nSiPixelError"),40,0,300000, 50, 0, 50);
       numErrorsPerEvent_vs_HF_[i]->SetDirectory(oFile_->GetDirectory(0));
     }
     for(int i = 0; i<7; i++)
     {
       if(i<3)      numErrorsPerEvent_vs_HF_10xBin_[i] = new TH2D(Form("numErrorsPerEvent_vs_HF_10xBin_B_L%d",i+1) ,Form("number of SiPixel Errors (Barrel Layer%d);HF Energy Sum;nSiPixelError",i+1),40*10, 0, 300000, 50, 0, 50);
       else if(i<5) numErrorsPerEvent_vs_HF_10xBin_[i] = new TH2D(Form("numErrorsPerEvent_vs_HF_10xBin_EC_D%d",i-2),Form("number of SiPixel Errors (Endcap Disk%d);HF Energy Sum;nSiPixelError", i-2), 40*10, 0, 300000, 50, 0, 50);
       else if(i<6) numErrorsPerEvent_vs_HF_10xBin_[i] = new TH2D(Form("numErrorsPerEvent_vs_HF_10xBin_B_L123"),    Form("number of SiPixel Errors (Barrel Layers : 1, 2, 3);HF Energy Sum;nSiPixelError"),40*10, 0, 300000, 50, 0, 50);
       else         numErrorsPerEvent_vs_HF_10xBin_[i] = new TH2D(Form("numErrorsPerEvent_vs_HF_10xBin_"),          Form("number of SiPixel Errors (Barrel Layers : 1, 2, 3);HF Energy Sum;nSiPixelError"),40*10, 0, 300000, 50, 0, 50);
       numErrorsPerEvent_vs_HF_10xBin_[i]->SetDirectory(oFile_->GetDirectory(0));
     }
     // "_no3740" block
     for(int i = 0; i<7; i++)
     {
       if(i<3)      numErrorsPerEvent_no3740_[i] = new TH1D(Form("numErrorsPerEvent_no3740_B_L%d",i+1) ,Form("number of SiPixel Errors (without error 37 and 40) (Barrel Layer%d);nSiPixelError",i+1), 50, 0, 50);
       else if(i<5) numErrorsPerEvent_no3740_[i] = new TH1D(Form("numErrorsPerEvent_no3740_EC_D%d",i-2),Form("number of SiPixel Errors (without error 37 and 40) (Endcap Disk%d);nSiPixelError",i-2),  50, 0, 50);
       else if(i<6) numErrorsPerEvent_no3740_[i] = new TH1D(Form("numErrorsPerEvent_no3740_B_L123"),    Form("number of SiPixel Errors (without error 37 and 40) (Barrel Layers : 1, 2, 3);nSiPixelError"), 50, 0, 50);
       else         numErrorsPerEvent_no3740_[i] = new TH1D(Form("numErrorsPerEvent_no3740_"),          Form("number of SiPixel Errors (without error 37 and 40) (Barrel Layers : 1, 2, 3);nSiPixelError"), 50, 0, 50);
       numErrorsPerEvent_no3740_[i]->SetDirectory(oFile_->GetDirectory(0));
     }
     for(int i = 0; i<7; i++)
     {
       if(i<3)      numErrorsPerEvent_no3740_10xRange_[i] = new TH1D(Form("numErrorsPerEvent_no3740_10xRange_B_L%d",i+1) ,Form("number of SiPixel Errors (without error 37 and 40) (Barrel Layer%d);nSiPixelError",i+1), 100, 0, 100);
       else if(i<5) numErrorsPerEvent_no3740_10xRange_[i] = new TH1D(Form("numErrorsPerEvent_no3740_10xRange_EC_D%d",i-2),Form("number of SiPixel Errors (without error 37 and 40) (Endcap Disk%d);nSiPixelError",i-2),  100, 0, 100);
       else if(i<6) numErrorsPerEvent_no3740_10xRange_[i] = new TH1D(Form("numErrorsPerEvent_no3740_10xRange_B_L123"),    Form("number of SiPixel Errors (without error 37 and 40) (Barrel Layers : 1, 2, 3);nSiPixelError"), 100, 0, 100);
       else         numErrorsPerEvent_no3740_10xRange_[i] = new TH1D(Form("numErrorsPerEvent_no3740_10xRange_"),          Form("number of SiPixel Errors (without error 37 and 40) (Barrel Layers : 1, 2, 3);nSiPixelError"), 100, 0, 100);
       numErrorsPerEvent_no3740_10xRange_[i]->SetDirectory(oFile_->GetDirectory(0));
     }
     for(int i = 0; i<7; i++)
     {
       if(i<3)      numErrorsPerEvent_no3740_vs_HF_[i] = new TH2D(Form("numErrorsPerEvent_no3740_vs_HF_B_L%d",i+1) ,Form("number of SiPixel Errors (without error 37 and 40)(Barrel Layer%d);HF Energy Sum;nSiPixelError",i+1),40,0,300000, 50, 0, 50);
       else if(i<5) numErrorsPerEvent_no3740_vs_HF_[i] = new TH2D(Form("numErrorsPerEvent_no3740_vs_HF_EC_D%d",i-2),Form("number of SiPixel Errors (without error 37 and 40)(Endcap Disk%d);HF Energy Sum;nSiPixelError",i-2),40,0,300000,  50, 0, 50);
       else if(i<6) numErrorsPerEvent_no3740_vs_HF_[i] = new TH2D(Form("numErrorsPerEvent_no3740_vs_HF_B_L123"),    Form("number of SiPixel Errors (without error 37 and 40)(Barrel Layers : 1, 2, 3);HF Energy Sum;nSiPixelError"),40,0,300000, 50, 0, 50);
       else         numErrorsPerEvent_no3740_vs_HF_[i] = new TH2D(Form("numErrorsPerEvent_no3740_vs_HF_"),          Form("number of SiPixel Errors (without error 37 and 40)(Barrel Layers : 1, 2, 3);HF Energy Sum;nSiPixelError"),40,0,300000, 50, 0, 50);
       numErrorsPerEvent_no3740_vs_HF_[i]->SetDirectory(oFile_->GetDirectory(0));
     }
     for(int i = 0; i<7; i++)
     {
       if(i<3)      numErrorsPerEvent_no3740_vs_HF_10xBin_[i] = new TH2D(Form("numErrorsPerEvent_no3740_vs_HF_10xBin_B_L%d",i+1) ,Form("number of SiPixel Errors (without error 37 and 40)(Barrel Layer%d);HF Energy Sum;nSiPixelError",i+1),40*10, 0, 300000, 50, 0, 50);
       else if(i<5) numErrorsPerEvent_no3740_vs_HF_10xBin_[i] = new TH2D(Form("numErrorsPerEvent_no3740_vs_HF_10xBin_EC_D%d",i-2),Form("number of SiPixel Errors (without error 37 and 40)(Endcap Disk%d);HF Energy Sum;nSiPixelError", i-2), 40*10, 0, 300000, 50, 0, 50);
       else if(i<6) numErrorsPerEvent_no3740_vs_HF_10xBin_[i] = new TH2D(Form("numErrorsPerEvent_no3740_vs_HF_10xBin_B_L123"),    Form("number of SiPixel Errors (without error 37 and 40)(Barrel Layers : 1, 2, 3);HF Energy Sum;nSiPixelError"),40*10, 0, 300000, 50, 0, 50);
       else         numErrorsPerEvent_no3740_vs_HF_10xBin_[i] = new TH2D(Form("numErrorsPerEvent_no3740_vs_HF_10xBin_"),          Form("number of SiPixel Errors (without error 37 and 40)(Barrel Layers : 1, 2, 3);HF Energy Sum;nSiPixelError"),40*10, 0, 300000, 50, 0, 50);
       numErrorsPerEvent_no3740_vs_HF_10xBin_[i]->SetDirectory(oFile_->GetDirectory(0));
     }
     HF_        = new TH1D("HF_", ";HF Energy Sum;", 30, 0, 300000);
     HF_->SetDirectory(oFile_->GetDirectory(0));
     HF_10xBin_ = new TH1D("HF_10xBin_", ";HF Energy Sum;", 30*10, 0, 300000);
     HF_10xBin_->SetDirectory(oFile_->GetDirectory(0));
     /////

  //FEDerror vs HFE for each FED
  for(int i = 0; i<41; ++i) {
    hFEDerrorHFE[i] = new TH2D(Form("hFEDerrorHFE_FED%d",i),Form("hFEDerrorHFE_FED%d;E_{HF};FED Error",i),40,0,300000,numErrorBins, errorBin1, errorBin2);
    hFEDerrorHFE[i]->SetDirectory(oFile_->GetDirectory(0));

    FEDOccPerFED_[i] = new TH2D(Form("FEDOccupancyPerFED_%d",i),Form("FED Occupancy (FED %d);HF Energy Sum;Average FED Occupancy",i),40,0,300000,100,0,3);
    FEDOccPerFED_[i]->SetDirectory(oFile_->GetDirectory(0));
    for(int j = 0; j<6; ++j) {
      // FEDOccPerFED_[i][j] = new TH2D(Form("FEDOccupancyPerFED_%d_L%d",i,j+1),Form("FED Occupancy (FED %d L %d);HF Energy Sum;Average FED Occupancy",i,j+1),40,0,300000,100,0,3);
      // FEDOccPerFED_[i][j]->SetDirectory(oFile_->GetDirectory(0));
    }
  }

  //  LinkByLinkOcc_ = new TH1D("LinkByLinkOcc",";36*detid+linkid;hits",1500,0,1500);
  //  LinkByLinkOcc_->SetDirectory(oFile_->GetDirectory(0));
 
  //Endcap Ntuples----------------------------------------------------------------------------
  // EndcapModuleNt_ = new TNtuple("EndcapPixelDistrib", "EndcapPixelDistri","idH:idL:side:disk:blade:panel:module:z:R:phi:eta:HMod0Hits:HMod1Hits:ncolumns:nrows",100000);
  // EndcapDigisNt_ = new TNtuple("EndcapDigis", "EndcapPixelDigisDist","idH:idL:side:disk:blade:panel:module:adc:column:row",100000);
  // EndcapColumnsNt_ = new TNtuple("EndcapColumnOcc", "EndcapColumnOcc","idH:idL:side:disk:blade:panel:gHmodule:column:colHits:gcolumn",100000);
  // EndcapDColumnsNt_ = new TNtuple("EndcapDColumnOcc", "EndcapDColumnOcc","idH:idL:side:disk:blade:panel:gHmodule:Dcolumn:DcolHits:gDcolumn",100000);

  // //General histos---------------------------------------------------------
  // DcolumnHits_ = new TH1F("DcolumnHits", "Hits per Dcolumn", 161, 0, 160);
  // DcolumnHits_ ->GetYaxis()->SetTitle("Entries");
  // DcolumnHits_->GetXaxis()->SetTitle("Hits per Dcolumn");
  // DcolumnHits_->SetDirectory(oFile_->GetDirectory(0));
 
  // Occupancy_ = new TH1D("Occupancy", "Occupancy", 3001, 0., 3.);
  // Occupancy_->GetYaxis()->SetTitle("Entries");
  // Occupancy_->GetXaxis()->SetTitle("Occupancy [%]");
  // Occupancy_->SetDirectory(oFile_->GetDirectory(0));

  // OccupancyZ_ = new TH2D("OccupancyZ", "Occupancy", 1020, -51., 51., 3001, 0., 3.);
  // OccupancyZ_->GetYaxis()->SetTitle("Occupancy [%]");
  // OccupancyZ_->GetXaxis()->SetTitle("Z");
  // OccupancyZ_->SetDirectory(oFile_->GetDirectory(0));
}


SiPixelAnalyzer::~SiPixelAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
SiPixelAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    pixelDigi_id_.clear();
    pixelDigi_detLayer_.clear();
    pixelDigi_fedID_.clear();
    pixelDigi_nHits_B_L1_.clear();
    pixelDigi_nHits_B_L2_.clear();
    pixelDigi_nHits_B_L3_.clear();
    pixelDigi_nHits_EC_D1_.clear();
    pixelDigi_nHits_EC_D2_.clear();
    pixelDigi_nHits_B_L123_.clear();
    pixelDigi_nHits_Total_.clear();
    moduleGeom_phi_.clear();
    moduleGeom_eta_.clear();
    moduleGeom_z_.clear();

    error_id_.clear();
    error_detLayer_.clear();
    error_fedID_.clear();
    error_type_.clear();

    FED_id_.clear();
    FED_occu_B_L1_.clear();
    FED_occu_B_L2_.clear();
    FED_occu_B_L3_.clear();
    FED_occu_EC_D1_.clear();
    FED_occu_EC_D2_.clear();
    FED_occu_B_L123_.clear();

//    std::clock_t    start_analyze, end_analyze;
//    std::clock_t    start_analyze1, start_analyze2;
//    start_analyze = std::clock();
//    std::cout << "started, CLOCK_REALTIME         : " << (double)CLOCK_REALTIME << std::endl;

    run_   = iEvent.id().run();
    event_ = iEvent.id().event();
    lumi_  = iEvent.id().luminosityBlock();
    orbitNumber_   = iEvent.orbitNumber();
    bunchCrossing_ = iEvent.bunchCrossing();
    std::cout<<"orbitNumber   = "<<orbitNumber_<<std::endl;
    std::cout<<"bunchCrossing = "<<bunchCrossing_<<std::endl;

   //centrality first
   double HFRecHitSum = 0;
   Handle<HFRecHitCollection> hits;
   iEvent.getByToken(srcHFhits_,hits);
   for( size_t ihit = 0; ihit<hits->size(); ++ ihit){
     const HFRecHit & rechit = (*hits)[ ihit ];
     HFRecHitSum += rechit.energy();
   }
   HFRecHitSum_ = HFRecHitSum;
   
   if(HFRecHitSum>150000) {
       int event = iEvent.id().event();
       int run = iEvent.id().run();
       int lumi = iEvent.id().luminosityBlock();
       std::cout << "------ run: " << run << "  lumi: " << lumi << " event: " << event << " HFESum: " << HFRecHitSum << std::endl;
   }
//   std::cout << "hits->size() = " <<hits->size() << std::endl;
//   std::cout << "HFRecHitSum  = " <<HFRecHitSum << std::endl;

   using namespace sipixelobjects;

   gnHModuleBarrel_=0;
   BarrelColumnsOffset_ =0;
   gnHModuleEndcap_=0;
   EndcapColumnsOffset_ =0;

//   start_analyze1 = std::clock();

   edm::Handle<edm::DetSetVector<PixelDigi> > pixelDigis;
   iEvent.getByLabel(src_, pixelDigis);

   // get pixel digi errors from the event
   // https://github.com/cms-sw/cmssw/blob/CMSSW_7_5_X/EventFilter/SiPixelRawToDigi/plugins/SiPixelRawToDigi.cc#L187
   // https://github.com/cms-sw/cmssw/blob/CMSSW_7_5_X/EventFilter/SiPixelRawToDigi/plugins/SiPixelRawToDigi.cc#L285
   edm::Handle<edm::DetSetVector<SiPixelRawDataError>> pixelRawDataErrors;
   iEvent.getByLabel(src_, pixelRawDataErrors);

   edm::ESHandle<SiPixelFedCablingMap> map;
   iSetup.get<SiPixelFedCablingMapRcd>().get( map );
//   const SiPixelFedCablingMap * theCablingMap = map.product();
//   std::cout <<"Map number " << map->version() << std::endl;
   
   edm::ESHandle<TrackerGeometry> tracker;
   iSetup.get<TrackerDigiGeometryRecord>().get( tracker );    
   trGeo_ = tracker.product();

   //looping over FED and Links
   //-------------------------------------------------------------------------
   //EDIT HERE FOR CONVERSION ERROR UNIQUE POINTER TO SIPIXELFEDCBLINGTREE
   std::unique_ptr<SiPixelFedCablingTree> CablingTree = map->cablingTree();
   const std::vector<const PixelFEDCabling *>  cabling = CablingTree->fedList();
   typedef std::vector<const PixelFEDCabling *>::const_iterator FEDiter;
   edm::DetSetVector<PixelDigi>::const_iterator DSViter;   

//   start_analyze2 = std::clock();

   int numFEDs = cabling.size();
//   std::cout << "numFEDs = " << numFEDs << std::endl;
//   int numPixelDigis_global = pixelDigis->size();

   std::vector<sipixelobjects::PixelROC> pixelROCs;
   std::vector<unsigned int> rawIDsROC;
   std::vector<unsigned int> fedIDs;
   for (std::vector<const PixelFEDCabling *>::const_iterator ifed=cabling.begin();
    ifed != cabling.end(); ifed++) {
     unsigned int fed = (**ifed).id();
     unsigned int numLink = (**ifed).numberOfLinks();
     for (unsigned int link=1; link <= numLink; link++) {
       const PixelFEDLink * pLink = (**ifed).link(link);
       if (pLink==0) continue;
       unsigned int linkId = pLink->id();
       if (linkId != 0 && linkId!= link)
           std::cout << "PROBLEM WITH LINK NUMBER!!!!" << std::endl;
       unsigned int numberROC = pLink->numberOfROCs();
       for (unsigned int roc=1; roc <= numberROC; roc++) {
         const PixelROC * pROC = pLink->roc(roc);
         if (pROC==0) continue;
         if (pROC->idInLink() != roc)
             std::cout << "PROBLEM WITH ROC NUMBER!!!!" << std::endl;
//         Key key = {fed, link, roc};
//         theMap[key] = (*pROC);
         pixelROCs.push_back(*pROC);
         rawIDsROC.push_back(pROC->rawId());
         fedIDs.push_back(fed);
       }
     }
   }

//   int numRawIDsROC = rawIDsROC.size();
//   std::cout<< "numRawIDsROC = " << numRawIDsROC <<std::endl;
   uint32_t nhits[numFEDs][6] = { {0} };
   uint32_t totalPix[numFEDs][6] = { {0} };

   // check if FED occupancy exceeds the given threshold
   bool monsterEventCandidate = false;
   bool belowDiagonal         = false;
   double FEDoccupancyThresholds[5]       = {1.5, 0.8, 0.3, 0.3, 0.3};
   double FEDoccupancyThresholds_lower[1] = {0.1};
   int HFRecHitSumThresh = 70000;

   double FEDoccupancies[numFEDs][6] = { {0} };
   int numErrors[7] = { 0 };
   int numErrors_no3740[7] = { 0 };

   nPixelDigi_ = pixelDigis->size();
   for( DSViter = pixelDigis->begin() ; DSViter != pixelDigis->end(); DSViter++) {
       unsigned int id = DSViter->id;
       unsigned int fedID = 999999;

       std::vector<unsigned int>::iterator it = std::find(rawIDsROC.begin(),rawIDsROC.end(), id);
       if(it == rawIDsROC.end()) continue;
       auto pos = it - rawIDsROC.begin();
       fedID = fedIDs.at(pos);

       edm::DetSet<PixelDigi>::const_iterator  begin = DSViter->data.begin();
       edm::DetSet<PixelDigi>::const_iterator  end   = DSViter->data.end();
       edm::DetSet<PixelDigi>::const_iterator iter;

       DetId  detId(id);
       int detLayer = -1;
       if(detId.subdetId() ==PixelSubdetector::PixelBarrel ) {                //selcting barrel modules
           PXBDetId  bdetid(id);
       detLayer  = bdetid.layer()-1;   // Layer:1,2,3. -1
       }
       else{
           PXFDetId  fdetid(id);
           detLayer  = fdetid.disk()+2; //1, 2, 3 +2 for endcap
       }

       const PixelGeomDetUnit* PixelModuleGeom = dynamic_cast<const PixelGeomDetUnit*> (trGeo_->idToDet(id));   //detector geometry -> it returns the center of the module
       uint32_t ncolumns = PixelModuleGeom->specificTopology().ncolumns(); //n of columns
       uint32_t nrows = PixelModuleGeom->specificTopology().nrows();       //n of rows

//       // temporary verbose
//       // position of pixels
//       std::cout<<"PixelModuleGeom->position().barePhi() = "<<PixelModuleGeom->position().barePhi()<<std::endl;
//       std::cout<<"PixelModuleGeom->position().eta() = "<<PixelModuleGeom->position().eta()<<std::endl;
//       std::cout<<"fedID = "<<fedID<<std::endl;

       int nhits_PixelDigi[6] = {0};
       int nhits_PixelDigi_Total = 0;

       totalPix[fedID][detLayer] += ncolumns*nrows;
       for ( iter = begin ; iter != end; iter++ ){  //llop over digi
           GlobalPixel global = {(*iter).row(), (*iter).column()};

           for(std::vector<sipixelobjects::PixelROC>::const_iterator it_pixelROCs = pixelROCs.begin(); it_pixelROCs!=pixelROCs.end(); ++it_pixelROCs){
               PixelROC roc = *it_pixelROCs;
               LocalPixel local = roc.toLocal(global);

               if(local.valid())
               {
                   nhits[fedID][detLayer] ++;
                   nhits_PixelDigi[detLayer]++;
                   nhits_PixelDigi_Total++;
                   break;
               }
           }
       }

       pixelDigi_id_.push_back(id);
       pixelDigi_detLayer_.push_back(detLayer);
       pixelDigi_fedID_.push_back(fedID);
       pixelDigi_nHits_B_L1_.push_back(nhits_PixelDigi[0]);
       pixelDigi_nHits_B_L2_.push_back(nhits_PixelDigi[1]);
       pixelDigi_nHits_B_L3_.push_back(nhits_PixelDigi[2]);
       pixelDigi_nHits_EC_D1_.push_back(nhits_PixelDigi[3]);
       pixelDigi_nHits_EC_D2_.push_back(nhits_PixelDigi[4]);
       pixelDigi_nHits_B_L123_.push_back(nhits_PixelDigi[0]+nhits_PixelDigi[1]+nhits_PixelDigi[2]);
       pixelDigi_nHits_Total_.push_back(nhits_PixelDigi_Total);
       moduleGeom_phi_.push_back(PixelModuleGeom->position().phi());
       moduleGeom_eta_.push_back(PixelModuleGeom->position().eta());
       moduleGeom_z_.push_back(PixelModuleGeom->position().z());
   }

 //Loop over FEDs and calculate FED occupancy per layer
   nFED_ = numFEDs;
   for(int j=0; j<numFEDs; ++j)
   {
//       std::cout << "Fed id: " << j << std::endl;
       for(int i = 0; i<5; i++){
           if(totalPix[j][i]!=0) {

               FEDoccupancies[j][i] = 100.0*nhits[j][i]/((double)totalPix[j][i]);

               if (FEDoccupancies[j][i] > FEDoccupancyThresholds[i])
               {
                   monsterEventCandidate = true;

                   std::cout<<"monsterEventCandidate = "<< monsterEventCandidate <<std::endl;
                   std::cout<<"fedID = j              = "<< j <<std::endl;
                   std::cout<<"detLayer = i           = "<< i <<std::endl;
                   std::cout<<"FEDoccupancies[j][i]     = "<< FEDoccupancies[j][i] <<std::endl;
                   std::cout<<"nhits[j][i]              = "<< nhits[j][i] <<std::endl;
                   std::cout<<"totalPix[j][i]           = "<< totalPix[j][i] <<std::endl;
                   std::cout<<"HFRecHitSum              = "<< HFRecHitSum <<std::endl;
               }
               if (i<1){
                   if (FEDoccupancies[j][i] < FEDoccupancyThresholds_lower[i] && HFRecHitSum > HFRecHitSumThresh)
                   {
                       belowDiagonal = true;

                       std::cout<<"belowDiagonal          = "<< belowDiagonal <<std::endl;
                       std::cout<<"fedID = j              = "<< j <<std::endl;
                       std::cout<<"detLayer = i           = "<< i <<std::endl;
                       std::cout<<"FEDoccupancies[j][i]     = "<< FEDoccupancies[j][i] <<std::endl;
                       std::cout<<"nhits[j][i]              = "<< nhits[j][i] <<std::endl;
                       std::cout<<"totalPix[j][i]           = "<< totalPix[j][i] <<std::endl;
                       std::cout<<"HFRecHitSum              = "<< HFRecHitSum <<std::endl;
                   }
               }

               FEDOcc_[i]->Fill(HFRecHitSum, FEDoccupancies[j][i]);
               FEDOcc_10xBin_[i]->Fill(HFRecHitSum, FEDoccupancies[j][i]);
//               std::cout << "Layer: " << i << " Occupancy: " << 100.0*nhits[j][i]/((double)totalPix[j][i]) << std::endl;
           }
       }
       if ((totalPix[j][0]+totalPix[j][1]+totalPix[j][2])!=0) {
           FEDoccupancies[j][5] = 100.0*(nhits[j][0]+nhits[j][1]+nhits[j][2])/((double)(totalPix[j][0]+totalPix[j][1]+totalPix[j][2]));

           FEDOcc_[5]->Fill(HFRecHitSum, FEDoccupancies[j][5]);
           FEDOcc_10xBin_[5]->Fill(HFRecHitSum, FEDoccupancies[j][5]);
       }

       FED_id_.push_back(j);
       FED_occu_B_L1_.push_back(FEDoccupancies[j][0]);
       FED_occu_B_L2_.push_back(FEDoccupancies[j][1]);
       FED_occu_B_L3_.push_back(FEDoccupancies[j][2]);
       FED_occu_EC_D1_.push_back(FEDoccupancies[j][3]);
       FED_occu_EC_D2_.push_back(FEDoccupancies[j][4]);
       FED_occu_B_L123_.push_back(FEDoccupancies[j][5]);
   }

   // fill Errors vs. HFRecHitSum histogram
//   std::cout<<"pixelDigis->size()         = "         <<pixelDigis->size()<<std::endl;
//   std::cout<<"pixelRawDataErrors->size() = " <<pixelRawDataErrors->size()<<std::endl;
   edm::DetSetVector<SiPixelRawDataError>::const_iterator iter_Errors;
   nError_ = 0;
   for(iter_Errors = pixelRawDataErrors->begin(); iter_Errors!=pixelRawDataErrors->end(); ++iter_Errors)
   {
       DetId  detId(iter_Errors->id);
       int detLayer = -1;
       if(detId.subdetId() ==PixelSubdetector::PixelBarrel ) {                //selcting barrel modules
           PXBDetId  bdetid(iter_Errors->id);
           detLayer  = bdetid.layer()-1;   // Layer:1,2,3. -1
       }
       else if (PixelSubdetector::PixelEndcap){
           PXFDetId  fdetid(iter_Errors->id);
           detLayer  = fdetid.disk()+2; //1, 2, 3 +2 for endcap
       }

//       std::cout<<"iter_Errors->id          = "<<iter_Errors->id<<std::endl;
//       std::cout<<"iter_Errors->data.size() = "<<iter_Errors->data.size()<<std::endl;
       edm::DetSet<SiPixelRawDataError>::const_iterator di;
       for(di = iter_Errors->data.begin(); di != iter_Errors->data.end(); di++) {
           nError_++;
           int FedId = di->getFedId();                  // FED the error came from
           int errorType = di->getType();               // type of error
           numErrors[6]++;
           if (errorType != 37 && errorType != 40)
           {
               numErrors_no3740[6]++;
           }

           std::cout<<"di->errorMessage_ = "<< di->getMessage() <<std::endl;
           std::cout<<"di->getFedId()    = "<< FedId <<std::endl;
           std::cout<<"di->getType()     = "<< errorType <<std::endl;

           if (detLayer >= 0 && detLayer<=5){
               numErrors[detLayer]++;
               if (errorType != 37 && errorType != 40)
               {
                   numErrors_no3740[detLayer]++;
               }

               FEDErrors_[detLayer]->Fill(HFRecHitSum,errorType);
               FEDErrors_10xBin_[detLayer]->Fill(HFRecHitSum,errorType);

                         if(FedId<0 || FedId>39)
            std::cout << "!!!!!!!!!!! FedId: " << FedId << std::endl;
               if(FedId>-1 && FedId<41)
            hFEDerrorHFE[FedId]->Fill(HFRecHitSum,errorType);
               
               FEDErrors_Type40_[detLayer]->Fill(HFRecHitSum,errorType);
               FEDErrors_Type40_10xBin_[detLayer]->Fill(HFRecHitSum,errorType);
               if(detLayer<3)      {
                   numErrors[5]++;
                   if (errorType != 37 && errorType != 40)
                   {
                       numErrors_no3740[5]++;
                   }
                   FEDErrors_[5]->Fill(HFRecHitSum,errorType);
                   FEDErrors_10xBin_[5]->Fill(HFRecHitSum,errorType);

                   FEDErrors_Type40_[5]->Fill(HFRecHitSum,errorType);
                   FEDErrors_Type40_10xBin_[5]->Fill(HFRecHitSum,errorType);
               }
           }
           else
           {
               std::cout<<"error-prone detLayer = "<<detLayer<<std::endl;
           }

           error_id_.push_back(iter_Errors->id);
           error_detLayer_.push_back(detLayer);
           error_fedID_.push_back(FedId);
           error_type_.push_back(errorType);
       }
   }
   HF_->Fill(HFRecHitSum);
   HF_10xBin_->Fill(HFRecHitSum);

   for(int i=0; i<7; ++i) {
       numErrorsPerEvent_[i]->Fill(numErrors[i]);
       numErrorsPerEvent_10xRange_[i]->Fill(numErrors[i]);

       numErrorsPerEvent_vs_HF_[i]->Fill(HFRecHitSum,numErrors[i]);
       numErrorsPerEvent_vs_HF_10xBin_[i]->Fill(HFRecHitSum,numErrors[i]);

       numErrorsPerEvent_no3740_[i]->Fill(numErrors_no3740[i]);
       numErrorsPerEvent_no3740_10xRange_[i]->Fill(numErrors_no3740[i]);

       numErrorsPerEvent_no3740_vs_HF_[i]->Fill(HFRecHitSum,numErrors_no3740[i]);
       numErrorsPerEvent_no3740_vs_HF_10xBin_[i]->Fill(HFRecHitSum,numErrors_no3740[i]);

       if(numErrors[i]>50){
           std::cout << "numErrors[i]>50" << std::endl;
           std::cout<<"detLayer = i = "<< i << "  numErrors[i] = "<< numErrors[i] <<std::endl;
       }
   }
   if(monsterEventCandidate) {
       std::cout<<"monsterEventCandidate = "<< monsterEventCandidate <<std::endl;
       for(int i=0; i<7; ++i) {

           std::cout<<"detLayer = i = "<< i << "  numErrors[i] = "<< numErrors[i] <<std::endl;
       }
   }
   if(belowDiagonal){
       std::cout<<"belowDiagonal = "<< belowDiagonal <<std::endl;
       for(int i=0; i<7; ++i) {

           std::cout<<"detLayer = i = "<< i << "  numErrors[i] = "<< numErrors[i] <<std::endl;
       }
   }

   // fill trees
   treeEvent_->Fill();
   treePixelDigi_->Fill();
   treePixelRawDataError_->Fill();
   treeFED_->Fill();

//   std::cout << "numPixelDigis_global       = " << numPixelDigis_global << std::endl;

//   end_analyze = std::clock();
   std::cout.precision(6);      // get back to default precision
//   std::cout << "analyzer finished in             : " << (end_analyze - start_analyze) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
//   std::cout << "analyzer1 finished in            : " << (end_analyze - start_analyze1) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
//   std::cout << "analyzer2 finished in            : " << (end_analyze - start_analyze2) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
//   std::cout << "finished, CLOCK_REALTIME         : " << (double)CLOCK_REALTIME << std::endl;
     
}

// ------------ method called once each job just before starting event loop  ------------
void 
SiPixelAnalyzer::beginJob(const edm::EventSetup&)
{
    
  //  oFile_ = new TFile((const char*)outputFile_.c_str(), "RECREATE");
  //FEDs studies---------------------------------------------------------------------
  // FEDNt_ = new TNtuple("FED", "FED","fedid:links:roc",100000);
  // LinksNt_= new TNtuple("Links", "Links","fedid:linkn:nHits",100000);
  // ROCOccNt_= new TNtuple("ROCOcc", "ROCOcc","idH:idL:ROCHits:SubDet:gHModule:ROCn:ModuleHits:nROCs",100000);
     
     
  //Barrel Ntuples----------------------------------------------------------------------------
  // BarrelModuleNt_ = new TNtuple("BarrelPixelDistrib", "BarrelPixelDistri","idH:idL:layer:ladder:ring:z:R:phi:eta:Mhits:HMod0Hits:HMod1Hits:ncolumns:nrows",100000);
  // BarrelDigisNt_ = new TNtuple("BarrelDigis", "BarrelPixelDigisDist","idH:idL:layer:ladder:ring:adc:column:row",100000);
  // BarrelColumnsNt_ = new TNtuple("BarrelColumnOcc", "BarrelColumnOcc","idH:idL:layer:ladder:ring:gHmodule:column:colHits:gcolumn",100000);
  // BarrelDColumnsNt_ = new TNtuple("BarrelDColumnOcc", "BarrelDColumnOcc","idH:idL:layer:ladder:ring:gHmodule:Dcolumn:DcolHits:gDcolumn",100000);
   
  // //histos---------------------------------------------------------
  // BarrelSlinkHitsDistrib_ = new TH1F("BarrelSlinkHitsDistrib", "Slink Hits -Barrel-", 801, 0., 800.);
  // BarrelSlinkHitsDistrib_->GetYaxis()->SetTitle("Entries");
  // BarrelSlinkHitsDistrib_->GetXaxis()->SetTitle("Hits per RO Link");
  // BarrelSlinkHitsDistrib_->SetDirectory(oFile_->GetDirectory(0));
 
  // //Endcap Ntuples----------------------------------------------------------------------------
  // EndcapModuleNt_ = new TNtuple("EndcapPixelDistrib", "EndcapPixelDistri","idH:idL:side:disk:blade:panel:module:z:R:phi:eta:HMod0Hits:HMod1Hits:ncolumns:nrows",100000);
  // EndcapDigisNt_ = new TNtuple("EndcapDigis", "EndcapPixelDigisDist","idH:idL:side:disk:blade:panel:module:adc:column:row",100000);
  // EndcapColumnsNt_ = new TNtuple("EndcapColumnOcc", "EndcapColumnOcc","idH:idL:side:disk:blade:panel:gHmodule:column:colHits:gcolumn",100000);
  // EndcapDColumnsNt_ = new TNtuple("EndcapDColumnOcc", "EndcapDColumnOcc","idH:idL:side:disk:blade:panel:gHmodule:Dcolumn:DcolHits:gDcolumn",100000);

  // //General histos---------------------------------------------------------
  // DcolumnHits_ = new TH1F("DcolumnHits", "Hits per Dcolumn", 161, 0, 160);
  // DcolumnHits_ ->GetYaxis()->SetTitle("Entries");
  // DcolumnHits_->GetXaxis()->SetTitle("Hits per Dcolumn");
  // DcolumnHits_->SetDirectory(oFile_->GetDirectory(0));
 
  // Occupancy_ = new TH1D("Occupancy", "Occupancy", 3001, 0., 3.);
  // Occupancy_->GetYaxis()->SetTitle("Entries");
  // Occupancy_->GetXaxis()->SetTitle("Occupancy [%]");
  // Occupancy_->SetDirectory(oFile_->GetDirectory(0));

  // OccupancyZ_ = new TH2D("OccupancyZ", "Occupancy", 1020, -51., 51., 3001, 0., 3.);
  // OccupancyZ_->GetYaxis()->SetTitle("Occupancy [%]");
  // OccupancyZ_->GetXaxis()->SetTitle("Z");
  // OccupancyZ_->SetDirectory(oFile_->GetDirectory(0));
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SiPixelAnalyzer::endJob() {
   oFile_->Write();
   oFile_->Close();
}


//define this as a plug-in
DEFINE_FWK_MODULE(SiPixelAnalyzer);
