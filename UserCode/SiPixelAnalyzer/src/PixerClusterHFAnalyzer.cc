// -*- C++ -*-
//
// Package:    UserCode/PixerClusterHFAnalyzer
// Class:      PixerClusterHFAnalyzer
// 
/**\class PixerClusterHFAnalyzer PixerClusterHFAnalyzer.cc UserCode/PixerClusterHFAnalyzer/plugins/PixerClusterHFAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marta Verweij
//         Created:  Fri, 13 Nov 2015 13:17:14 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "CommonTools/Utils/interface/TFileDirectory.h"

//ROOT inclusion
#include "TROOT.h"
//#include "TFile.h"
//#include "TNtuple.h"
//#include "TH1F.h"
#include "TH2F.h"

class PixerClusterHFAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit PixerClusterHFAnalyzer(const edm::ParameterSet&);
      ~PixerClusterHFAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
  edm::InputTag inputTag_;          // input tag identifying product containing pixel digis
  edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster> > inputToken_;
  //edm::EDGetTokenT<edm::DetSetVector<SiPixelDigi> > siPixelToken_;
  edm::EDGetTokenT<HFRecHitCollection> HFHitsToken_;
  edm::InputTag HFHits_;
  double eCut_HF_;

//  TFile* oFile_;
  //edm::Service<TFileService> fs_;
  TH2F *fh2NClusSumHF;
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
PixerClusterHFAnalyzer::PixerClusterHFAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   //usesResource("TFileService");

   inputTag_     = iConfig.getParameter<edm::InputTag>("inputTag");
   HFHits_       = iConfig.getParameter<edm::InputTag>("HFHitCollection");
   inputToken_   = consumes<edmNew::DetSetVector<SiPixelCluster> >(inputTag_);
   //siPixelToken_ = consumes<edmNew::DetSetVector<SiPixelDigi> >(siPixelTag_);
   HFHitsToken_  = consumes<HFRecHitCollection>(HFHits_);
   eCut_HF_      = iConfig.getParameter<double>("eCut_HF");

   //oFile_ = new TFile("histos.root", "RECREATE");//
   //fh2NClusSumHF = new TH2F("fh2NClusSumHF","fh2NClusSumHF;sum E HF;N pix clusters",50,0.,50000.,50,0.,50000.);
   //fh2NClusSumHF->SetDirectory(oFile_->GetDirectory(0));
//  fh2NClusSumHF = fs_->make<TH2F>( "fh2NClusSumHF"  , "sum E HF", "N pix clusters", 50, 0. ,50000., 50, 0. ,50000.);// 100,  0., 100. );
}


PixerClusterHFAnalyzer::~PixerClusterHFAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PixerClusterHFAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


 // get hold of products from Event
  edm::Handle<edmNew::DetSetVector<SiPixelCluster> > clusterColl;
  iEvent.getByToken(inputToken_, clusterColl);

  //edm::Handle<edm::DetSetVector<PixelDigi> > siPixelColl;
  //iEvent.getByToken(siPixelToken_, siPixelColl);

  unsigned int clusterSize = clusterColl->dataSize();//size();//
  LogDebug("") << "Number of clusters accepted: " << clusterSize;

  edm::Handle<HFRecHitCollection> HFRecHitsH;
  iEvent.getByToken(HFHitsToken_,HFRecHitsH);

  double sumE = 0.;

  for (HFRecHitCollection::const_iterator it=HFRecHitsH->begin(); it!=HFRecHitsH->end(); it++) {
    if (it->energy()>eCut_HF_) {
      sumE += it->energy();
    }
  }
  //  Printf("clusterSize: %d  sumE: %f  ",clusterSize,sumE);
  fh2NClusSumHF->Fill(sumE,clusterSize);

}


// ------------ method called once each job just before starting event loop  ------------
void 
PixerClusterHFAnalyzer::beginJob()
{
  edm::Service<TFileService> fs_;
  //TFileDirectory subDir = fs_->mkdir( "mySubDirectory" );
  fh2NClusSumHF = fs_->make<TH2F>("sumE_HF", "N pix clusters", 100, 0. ,100000., 50, 0. ,50000.);
  //fh2NClusSumHF = subDir.make<TH2F>("sumE_HF", "N pix clusters", 50, 0. ,50000., 50, 0. ,50000.);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PixerClusterHFAnalyzer::endJob() 
{
   //oFile_->Write();
  // oFile_->Close();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PixerClusterHFAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PixerClusterHFAnalyzer);
