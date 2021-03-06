// -*- C++ -*-
//
// Package:    PixelHitAnalyzer
// Class:      PixelHitAnalyzer
// 
/**\class PixelHitAnalyzer PixelHitAnalyzer.cc MitHig/PixelHitAnalyzer/src/PixelHitAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     Prepare the Hit Tree for analysis
*/
//
// Original Author:  Yilmaz Yetkin, Yen-Jie 
//         Created:  Tue Sep 30 15:14:28 CEST 2008
// $Id: PixelHitAnalyzer.cc,v 1.28 2010/11/09 00:14:40 yjlee Exp $
//
//


// system include files
#include <memory>
#include <iostream>
#include <vector>
#include <string>
#include <map>

// CMSSW user include files
#include "DataFormats/Common/interface/DetSetAlgorithm.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
//#include "Geometry/TrackerGeometryBuilder/interface/TrackerLayerIdAccessor.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

//#include "L1Trigger/GlobalTrigger/interface/L1GlobalTrigger.h"

#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"

// Heavyion
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
//#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"


// Root include files
#include "TTree.h"

using namespace std;
using namespace edm;
using namespace reco;

//
// class decleration
//

#define PI 3.14159265358979

#define MAXPARTICLES 500000
#define MAXHITS 500000
#define MAXVTX 100
#define MAXHLTBITS 256

struct PixelEvent{

   int nRun;
   int nEv;
   int nLumi;
   int nBX;
   
   int nhits1;
   int nhits2;
   int nhits3;
   int nhits4;
   int nhits5;
   int ntrks;
   int ntrksCut;
   int nHFp;
   int nHFn;
   
   int mult;
   
   // Vertex information
   int nv;
   float vx[MAXVTX];
   float vy[MAXVTX];
   float vz[MAXVTX];
   float vxError[MAXVTX];
   float vyError[MAXVTX];
   float vzError[MAXVTX];

   float beamSpotX;
   float beamSpotY;
   float beamSpotZ;

   // First layer hit
   float eta1[MAXHITS];
   float phi1[MAXHITS];
   float r1[MAXHITS];
   float cs1[MAXHITS];
   float ch1[MAXHITS];
   //int id1[MAXHITS];
   //int gp1[MAXHITS];
   //int type1[MAXHITS];

   // Second layer hit
   float eta2[MAXHITS];
   float phi2[MAXHITS];
   float r2[MAXHITS];
   float cs2[MAXHITS];
   float ch2[MAXHITS];
   //int id2[MAXHITS];
   //int gp2[MAXHITS];
   //int type2[MAXHITS];

   // Third layer hit
   float eta3[MAXHITS];
   float phi3[MAXHITS];
   float r3[MAXHITS];
   float cs3[MAXHITS];
   float ch3[MAXHITS];
   //int id3[MAXHITS];
   //int gp3[MAXHITS];
   //int type3[MAXHITS];

   // First disk hit
   float eta4[MAXHITS];
   float phi4[MAXHITS];
   float r4[MAXHITS];
   float cs4[MAXHITS];
   float ch4[MAXHITS];
   //int id4[MAXHITS];
   //int gp4[MAXHITS];
   //int type4[MAXHITS];

   // Second disk hit
   float eta5[MAXHITS];
   float phi5[MAXHITS];
   float r5[MAXHITS];
   float cs5[MAXHITS];
   float ch5[MAXHITS];
   //int id5[MAXHITS];
   //int gp5[MAXHITS];
   //int type5[MAXHITS];

   // genparticle
   int nparticle;
   float pt[MAXPARTICLES];
   float eta[MAXPARTICLES];
   float phi[MAXPARTICLES];
   int pdg[MAXPARTICLES];
   int chg[MAXPARTICLES];
   int evtType;
   float xi;
   bool passDS;
   bool passSingleTrack;
   //float x[MAXPARTICLES];
   //float y[MAXPARTICLES];
   //float z[MAXPARTICLES];
   
   // hlt
   int nHLTBit;
   bool hltBit[MAXHLTBITS];

   // l1
   int nL1TBit;
   bool l1TBit[MAXHLTBITS];
   int nL1ABit;
   bool l1ABit[MAXHLTBITS];

   // HI
   int cBin;
   int nbins;
   int binsize;
   float hf;
   float hftp;
   float hftm;
   float hfp;
   float hfm;
   float eb;
   float eep;
   float eem;
   float etMid;
   float nparti;
   float npartiSigma;
   float ncoll;
   float ncollSigma;
   float nhard;
   float nhardSigma;
   float b;
   float bSigma;
   float pixel;
   float zdcp;
   float zdcm;
};

class PixelHitAnalyzer : public edm::EDAnalyzer {
   public:
      explicit PixelHitAnalyzer(const edm::ParameterSet&);
      ~PixelHitAnalyzer();

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

   void fillVertices(const edm::Event& iEvent);
   void fillBeamSpot(const edm::Event& iEvent);
   void fillHits(const edm::Event& iEvent);
   void fillParticles(const edm::Event& iEvent);
   void fillPixelTracks(const edm::Event& iEvent);
   void fillL1Bits(const edm::Event& iEvent);
   void fillHLTBits(const edm::Event& iEvent);
   void fillHF(const edm::Event& iEvent);
   void fillTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup);
   void fillCentrality(const edm::Event& iEvent, const edm::EventSetup& iSetup);

   static bool pairCompare(pair<double,HepMC::FourVector> i, pair<double,HepMC::FourVector> j);
   
  //   int associateSimhitToTrackingparticle(unsigned int trid );
   bool checkprimaryparticle(const TrackingParticle* tp);

      // ----------member data ---------------------------

   bool doMC_;
   bool doHF_;
   bool doCentrality_;
   bool doTrackingParticle_;
   bool doPixel_;
   bool doTrack_;
   bool doBeamSpot_;

  //vector<string> vertexSrc_;
  //edm::InputTag trackSrc_;
  //edm::InputTag TowerSrc_;  // > 2.87 for HF
  //edm::InputTag beamSpotProducer_;

   vector<edm::EDGetTokenT<reco::VertexCollection> > vertexSrc_;
   edm::EDGetTokenT<vector<Track> > trackSrc_;
   edm::EDGetTokenT<CaloTowerCollection> TowerSrc_;
   edm::EDGetTokenT<reco::BeamSpot> beamSpotProducer_;
   edm::EDGetTokenT<SiPixelRecHitCollection> siPixelSrc_;

   edm::EDGetTokenT<TrackingVertexCollection> simVertexSrc_;
  
  //   edm::InputTag L1gtReadout_; 
   edm::InputTag CentralityTag_; 
   double etaMult_;

   const TrackerGeometry* geo_;
   edm::Service<TFileService> fs;           
   edm::ESHandle < ParticleDataTable > pdt;
   edm::Handle<TrackingParticleCollection> trackingParticles;

   map<int,int> tpmap_;

   edm::InputTag                 hltResName_;         //HLT trigger results name
   std::vector<std::string>      hltProcNames_;       //HLT process name(s)
   std::vector<std::string>      hltTrgNames_;        //HLT trigger name(s)

   std::vector<int>              hltTrgBits_;         //HLT trigger bit(s)
   std::vector<bool>             hltTrgDeci_;         //HLT trigger descision(s)
   std::vector<std::string>      hltTrgUsedNames_;    //HLT used trigger name(s)
   std::string                   hltUsedResName_;     //used HLT trigger results name

//   CentralityProvider * centrality_;

   // Root object
   TTree* pixelTree_;

   PixelEvent pev_;

};

//--------------------------------------------------------------------------------------------------
// PixelHitAnalyzer::PixelHitAnalyzer(const edm::ParameterSet& iConfig)

// {
//    doMC_             = iConfig.getUntrackedParameter<bool>  ("doMC",true);
//    doCentrality_             = iConfig.getUntrackedParameter<bool>  ("doCentrality",false);
//    doHF_             = iConfig.getUntrackedParameter<bool>  ("doHF",true);
//    doTrackingParticle_             = iConfig.getUntrackedParameter<bool>  ("doTrackingParticle",false);
//    doPixel_             = iConfig.getUntrackedParameter<bool>  ("doPixel",true);
//    doTrack_             = iConfig.getUntrackedParameter<bool>  ("doTrack",true);
//    doBeamSpot_             = iConfig.getUntrackedParameter<bool>  ("doBeamSpot",true);
//    //vertexSrc_ = iConfig.getParameter<vector<string> >("vertexSrc");
//    std::vector<std::string> vertexSrcString_ = iConfig.getParameter<vector<string> >("vertexSrc");
//    for(unsigned i = 0; i < vertexSrcString_.size(); i++) {
//      vertexSrc_.push_back(consumes<reco::VertexCollection>(edm::InputTag(vertexSrcString_[i])));
//    }
   
//    etaMult_ = iConfig.getUntrackedParameter<double>  ("nHitsRegion",1.);
//    trackSrc_ = consumes<vector<Track> > (iConfig.getUntrackedParameter<edm::InputTag>("trackSrc",edm::InputTag("generalTracks")));
//    siPixelSrc_ = consumes<SiPixelRecHitCollection>(edm::InputTag("siPixelRecHits"));
   
//    //trackSrc_ =    iConfig.getUntrackedParameter<edm::InputTag>("trackSrc",edm::InputTag("generalTracks"));
//    //L1gtReadout_ = iConfig.getUntrackedParameter<edm::InputTag>("L1gtReadout",edm::InputTag("gtDigis"));

//    //edm::InputTag TowerSrc_;  // > 2.87 for HF
//    //edm::InputTag beamSpotProducer_;

//    //edm::EDGetTokenT<vector<Track> > trackSrc_;
//    // edm::EDGetTokenT<reco::BeamSpot> beamSpotProducer_;
//    // edm::EDGetTokenT<CaloTowerCollection> TowerSrc_;
  
//    //edm::InputTag L1gtReadout_; 
//    edm::InputTag CentralityTag_; 
//    //double etaMult_;

//    //const TrackerGeometry* geo_;
//    edm::Service<TFileService> fs;           
//    edm::ESHandle < ParticleDataTable > pdt;
//    edm::Handle<TrackingParticleCollection> trackingParticles;

//    map<int,int> tpmap_;

//    edm::InputTag                 hltResName_;         //HLT trigger results name
//    std::vector<std::string>      hltProcNames_;       //HLT process name(s)
//    std::vector<std::string>      hltTrgNames_;        //HLT trigger name(s)

//    std::vector<int>              hltTrgBits_;         //HLT trigger bit(s)
//    std::vector<bool>             hltTrgDeci_;         //HLT trigger descision(s)
//    std::vector<std::string>      hltTrgUsedNames_;    //HLT used trigger name(s)
//    std::string                   hltUsedResName_;     //used HLT trigger results name

// //   CentralityProvider * centrality_;

//    // Root object
//    TTree* pixelTree_;

//    PixelEvent pev_;

// };

//--------------------------------------------------------------------------------------------------
PixelHitAnalyzer::PixelHitAnalyzer(const edm::ParameterSet& iConfig)

{
   doMC_             = iConfig.getUntrackedParameter<bool>  ("doMC",true);
   doCentrality_             = iConfig.getUntrackedParameter<bool>  ("doCentrality",false);
   doHF_             = iConfig.getUntrackedParameter<bool>  ("doHF",true);
   doTrackingParticle_             = iConfig.getUntrackedParameter<bool>  ("doTrackingParticle",false);
   doPixel_             = iConfig.getUntrackedParameter<bool>  ("doPixel",true);
   doTrack_             = iConfig.getUntrackedParameter<bool>  ("doTrack",true);
   doBeamSpot_             = iConfig.getUntrackedParameter<bool>  ("doBeamSpot",true);
   std::vector<std::string> vertexSrcString_ = iConfig.getParameter<vector<string> >("vertexSrc");
   for(unsigned i = 0; i < vertexSrcString_.size(); i++) {
     vertexSrc_.push_back(consumes<reco::VertexCollection>(edm::InputTag(vertexSrcString_[i])));
   }
   //   vertexSrc_ = iConfig.getParameter<vector<string> >("vertexSrc");
   etaMult_ = iConfig.getUntrackedParameter<double>  ("nHitsRegion",1.);
   trackSrc_ = consumes<vector<Track> > (iConfig.getUntrackedParameter<edm::InputTag>("trackSrc",edm::InputTag("generalTracks")));
   //trackSrc_ =    iConfig.getUntrackedParameter<edm::InputTag>("trackSrc",edm::InputTag("generalTracks"));
   //L1gtReadout_ = iConfig.getUntrackedParameter<edm::InputTag>("L1gtReadout",edm::InputTag("gtDigis"));
   TowerSrc_ = consumes<CaloTowerCollection> (iConfig.getUntrackedParameter<edm::InputTag>("towersSrc",edm::InputTag("towerMaker")));
   //TowerSrc_ =   iConfig.getUntrackedParameter<edm::InputTag>("towersSrc",edm::InputTag("towerMaker"));
   hltResName_ = iConfig.getUntrackedParameter<edm::InputTag>("hltTrgResults",edm::InputTag("TriggerResults"));
   beamSpotProducer_  = consumes<reco::BeamSpot>(iConfig.getUntrackedParameter<edm::InputTag>("beamSpotSrc",edm::InputTag("offlineBeamSpot")));
   //beamSpotProducer_  = iConfig.getUntrackedParameter<edm::InputTag>("towersSrc",edm::InputTag("offlineBeamSpot"));   
   CentralityTag_ = iConfig.getUntrackedParameter<edm::InputTag>("CentralityTag",edm::InputTag("pACentrality"));

   siPixelSrc_ = consumes<SiPixelRecHitCollection>(edm::InputTag("siPixelRecHits"));

   // if it's not MC, don't do TrackingParticle
   if (doMC_ == false) doTrackingParticle_ = false;
   if(doMC_==true)   simVertexSrc_ =  consumes<TrackingVertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tpVtxSrc",edm::InputTag("mix","MergedTrackTruth")));
   if (iConfig.exists("hltTrgNames"))
    hltTrgNames_ = iConfig.getUntrackedParameter<vector<string> >("hltTrgNames");

   if (iConfig.exists("hltProcNames"))
      hltProcNames_ = iConfig.getUntrackedParameter<vector<string> >("hltProcNames");
   else {
      hltProcNames_.push_back("FU");
      hltProcNames_.push_back("HLT");
   }

   
}

//--------------------------------------------------------------------------------------------------
PixelHitAnalyzer::~PixelHitAnalyzer()
{
}

//--------------------------------------------------------------------------------------------------
void
PixelHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   edm::ESHandle<TrackerGeometry> tGeo;
   iSetup.get<TrackerDigiGeometryRecord>().get(tGeo);
   geo_ = tGeo.product();
   iSetup.getData(pdt);

   tpmap_.clear();
   pev_.nhits1 = 0;
   pev_.nhits2 = 0;
   pev_.nhits3 = 0;
   pev_.nhits4 = 0;
   pev_.nhits5 = 0;
   pev_.ntrks = 0;
   pev_.ntrksCut = 0;
   pev_.mult = 0;
   pev_.nparticle = 0;
   pev_.nEv = (int)iEvent.id().event();
   pev_.nRun = (int)iEvent.id().run();
   pev_.nLumi = (int)iEvent.luminosityBlock();
   pev_.nBX = (int)iEvent.bunchCrossing();
   pev_.xi = -1;
   pev_.nv = 0;

   if (doTrack_) fillTracks(iEvent, iSetup);

   //if (doMC_) fillParticles(iEvent);
//   cout <<"Fill Vtx"<<endl;
   fillVertices(iEvent);
//   cout <<"Fill BeamSpot"<<endl;
   fillBeamSpot(iEvent);
//   cout <<"Fill Hits"<<endl;
   if (doPixel_) fillHits(iEvent);
   fillTracks(iEvent, iSetup);
//   cout <<"Fill L1"<<endl;
   fillL1Bits(iEvent);
//   cout <<"Fill HLT"<<endl;
   fillHLTBits(iEvent);
//   cout <<"Fill HF"<<endl;
   fillHF(iEvent);
//   cout <<"Fill Centrality"<<endl;
   if (doCentrality_) fillCentrality(iEvent, iSetup);
   map<int,int>::iterator begin = tpmap_.begin();
   map<int,int>::iterator end = tpmap_.end();

   pixelTree_->Fill();
}

//--------------------------------------------------------------------------------------------------
void
PixelHitAnalyzer::fillTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Handle<vector<Track> > etracks;
  iEvent.getByToken(trackSrc_, etracks);
  //iEvent.getByLabel(trackSrc_, etracks);
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(beamSpotProducer_,recoBeamSpotHandle);
  //iEvent.getByLabel(beamSpotProducer_,recoBeamSpotHandle);
  beamSpot = *recoBeamSpotHandle;

  pev_.ntrks=0;
  pev_.ntrksCut=0;
  for(unsigned it=0; it<etracks->size(); ++it){
    const reco::Track & etrk = (*etracks)[it];
    reco::TrackRef trackRef=reco::TrackRef(etracks,it);
    pev_.ntrks++;

    if (etrk.pt()<0.4||fabs(etrk.eta())>0.8) continue;
//    if(fiducialCut_ && hitDeadPXF(etrk)) continue; // if track hits the dead region, igonore it;

    if(etrk.quality(reco::TrackBase::qualityByName("highPurity")) != 1) continue;


    math::XYZPoint v1(pev_.vx[1],pev_.vy[1], pev_.vz[1]);
    double dz=etrk.dz(v1);
    double dzError=sqrt(etrk.dzError()*etrk.dzError()+pev_.vzError[1]*pev_.vzError[1]);
    double dxy=etrk.dxy(v1);
    double dxyError=sqrt(etrk.dxyError()*etrk.dxyError()+pev_.vxError[1]*pev_.vxError[1]+pev_.vyError[1]*pev_.vyError[1]);

   if (fabs(dz/dzError)>3||fabs(dxy/dxyError)>3) continue; 
   
    pev_.ntrksCut++;
  }

}

//--------------------------------------------------------------------------------------------------
void
PixelHitAnalyzer::fillVertices(const edm::Event& iEvent){

   // Vertex 0 : pev_vz[0] MC information from TrackingVertexCollection
   // Vertex 1 - n : Reconstructed Vertex from various of algorithms
   if(doMC_){
     // unsigned int daughter = 0;
//      int nVertex;
//      nVertex = 0;
     //int greatestvtx = 0;
      //  simVertexSrc_ =  consumes<TrackingVertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tpVtxSrc",edm::InputTag("mix","MergedTrackTruth")));
      Handle<TrackingVertexCollection> simvertices;
      iEvent.getByToken(simVertexSrc_, simvertices);
      //edm::Handle<SimVertexContainer> SimVtx;
      //iEvent.getByLabel("g4SimHits",SimVtx);
 
      if (doTrackingParticle_) {
         // Handle<TrackingVertexCollection> vertices;
         // iEvent.getByLabel("mergedtruth","MergedTrackTruth", vertices);
         // //nVertex = vertices->size();
         // for (unsigned int i = 0 ; i< vertices->size(); ++i){
   	 //    daughter = (*vertices)[i].nDaughterTracks();
   	 //    if(  daughter >(*vertices)[greatestvtx].nDaughterTracks()
   	 //       &&fabs((*vertices)[i].position().z())<30) greatestvtx = i;
         // }
      
         // if(vertices->size()>0&&fabs((*vertices)[greatestvtx].position().z())<30){
   	 //    pev_.vx[pev_.nv] = (*vertices)[greatestvtx].position().x();
   	 //    pev_.vy[pev_.nv] = (*vertices)[greatestvtx].position().y();
	 //    pev_.vz[pev_.nv] = (*vertices)[greatestvtx].position().z();
	    
         // }else{
	    pev_.vx[pev_.nv] =  -99; 
	    pev_.vy[pev_.nv] =  -99; 
	    pev_.vz[pev_.nv] =  -99; 
            // }
      } else {
         pev_.vx[pev_.nv] = -99;
         pev_.vy[pev_.nv] = -99;
         pev_.vz[pev_.nv] = -99;

         if(simvertices->size()>0) {
           if (fabs((*simvertices)[0].position().z())<30) {
             pev_.vx[0]=(*simvertices)[0].position().x();//simVtxItr->position().X();
             pev_.vy[0]=(*simvertices)[0].position().y();
             pev_.vz[0]=(*simvertices)[0].position().z();
           }
         }
         
        //  edm::SimVertexContainer::const_iterator simVtxItr= SimVtx->begin();
        // // for (unsigned int iv=0; iv<SimVtx->size(); iv++) {
        //  if (SimVtx->size()>0) {
   	//     for (unsigned int iv=0; iv<=0; iv++) {
	//        //itr->momentum();
        //        if (fabs(simVtxItr->position().Z())<30) {
   	//           pev_.vx[0]=simVtxItr->position().X();
        //           pev_.vy[0]=simVtxItr->position().Y();
        //           pev_.vz[0]=simVtxItr->position().Z();
        //        }
	//        simVtxItr++;
	//        //if (fabs(pev_.vz[0]>100) cout <<"!!!"<<endl;
	//        //cout <<pev_.vx[0]<<" "<<pev_.vy[0]<<" "<<pev_.vz[0]<<" "<<simVtxItr->position().T()<<endl;
	//     }
	//  }  
	//  //cout <<"============"<<endl;
      }

      pev_.nv++;
   } else {
      // Fill a dummy MC information
      pev_.vz[pev_.nv] = -99;
      pev_.nv++;
   }
   
   // Fill reconstructed vertices.   
   for(unsigned int iv = 0; iv < vertexSrc_.size(); ++iv){
      const reco::VertexCollection * recoVertices;
      edm::Handle<reco::VertexCollection> vertexCollection;
      iEvent.getByToken(vertexSrc_[iv],vertexCollection);
      //iEvent.getByLabel(vertexSrc_[iv],vertexCollection);
      recoVertices = vertexCollection.product();
      unsigned int daughter = 0;
      int greatestvtx = 0;
      
      for (unsigned int i = 0 ; i< recoVertices->size(); ++i){
	 daughter = (*recoVertices)[i].tracksSize();
	 if( daughter > (*recoVertices)[greatestvtx].tracksSize()) greatestvtx = i;
	 //         cout <<"Vertex: "<< (*recoVertices)[i].position().z()<<" "<<daughter<<endl;
      }
      
      if(recoVertices->size()>0){
	 pev_.vx[pev_.nv] = (*recoVertices)[greatestvtx].position().x();
	 pev_.vy[pev_.nv] = (*recoVertices)[greatestvtx].position().y();
	 pev_.vz[pev_.nv] = (*recoVertices)[greatestvtx].position().z();
	 pev_.vxError[pev_.nv] = (*recoVertices)[greatestvtx].xError();
	 pev_.vyError[pev_.nv] = (*recoVertices)[greatestvtx].yError();
	 pev_.vzError[pev_.nv] = (*recoVertices)[greatestvtx].zError();
      }else{
	 pev_.vx[pev_.nv] =  -99;
	 pev_.vy[pev_.nv] =  -99;
	 pev_.vz[pev_.nv] =  -99;
	 pev_.vxError[pev_.nv] = 0;
	 pev_.vyError[pev_.nv] = 0;
	 pev_.vzError[pev_.nv] = 0;
       }
      pev_.nv++;
   }

}

//--------------------------------------------------------------------------------------------------
void
PixelHitAnalyzer::fillBeamSpot(const edm::Event& iEvent){

  // Get the Beam Spot
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(beamSpotProducer_,recoBeamSpotHandle);
  //iEvent.getByLabel(beamSpotProducer_,recoBeamSpotHandle);
  beamSpot = *recoBeamSpotHandle;

  pev_.beamSpotX = beamSpot.x0();
  pev_.beamSpotY = beamSpot.y0();
  pev_.beamSpotZ = beamSpot.z0();

}


//--------------------------------------------------------------------------------------------------
void
PixelHitAnalyzer::fillHits(const edm::Event& iEvent){

//   double matchEtaMax = 0.005;
//   double matchPhiMax = 0.01;

  //if(doMC_&&doTrackingParticle_) iEvent.getByLabel("mergedtruth","MergedTrackTruth",trackingParticles);
   
   const SiPixelRecHitCollection* rechits;
   Handle<SiPixelRecHitCollection> rchts;
   //iEvent.getByLabel("siPixelRecHits",rchts);
   iEvent.getByToken(siPixelSrc_,rchts);
   rechits = rchts.product();

   for (SiPixelRecHitCollection::const_iterator it = rechits->begin(); it!=rechits->end();it++)
   {
      SiPixelRecHitCollection::DetSet hits = *it;
      DetId detId = DetId(hits.detId());
      SiPixelRecHitCollection::const_iterator recHitMatch = rechits->find(detId);
      const SiPixelRecHitCollection::DetSet recHitRange = *recHitMatch;
      unsigned int detType=detId.det();    // det type, tracker=1
      unsigned int subid=detId.subdetId(); //subdetector type, barrel=1, fpix=2
      
      PXBDetId pdetIdB = PXBDetId(detId);
      PXFDetId pdetIdF = PXFDetId(detId);
      unsigned int layer=0;
      layer=pdetIdB.layer();
      unsigned int disk=0;
      disk=pdetIdF.disk();
      //cout <<subid<<" "<<layer<<" "<<disk<<endl;
      if (detType!=1) continue;

      for ( SiPixelRecHitCollection::DetSet::const_iterator recHitIterator = recHitRange.begin(); 
	 recHitIterator != recHitRange.end(); ++recHitIterator) {
         const SiPixelRecHit * recHit = &(*recHitIterator);

         // SIM INFO
	 /*
         bool isprimary    = false;
         bool issecondary  = false;
         bool isbackground = false;
         int ptype = -99;
//         int gpid = -9999;
 //        int trid = -9999;
*/
         const PixelGeomDetUnit* pixelLayer = dynamic_cast<const PixelGeomDetUnit*> (geo_->idToDet(recHit->geographicalId()));
         GlobalPoint gpos = pixelLayer->toGlobal(recHit->localPosition());
         math::XYZVector rechitPos(gpos.x(),gpos.y(),gpos.z());

         // position
         double eta = rechitPos.eta();
         double phi = rechitPos.phi();
         double r   = rechitPos.rho();

         /*
         if (doMC_&&doTrackingParticle_) {
            
            TrackerHitAssociator theHitAssociator(iEvent);
            vector<PSimHit> simHits1 = theHitAssociator.associateHit(*recHitIterator);
            const PSimHit * bestSimHit1 = 0;
            int simIdx =0;

            //gets the primary simhit and its specifications for the rechit	       
            for(vector<PSimHit>::const_iterator simHit1 = simHits1.begin(); simHit1!= simHits1.end(); simHit1++){  
  	       simIdx++;
 
	       GlobalPoint simpos = pixelLayer->toGlobal((*simHit1).localPosition());
	       math::XYZVector simhitPos(simpos.x(),simpos.y(),simpos.z());

	       double simrecdeltaphi = fabs(simhitPos.phi()-phi)>matchPhiMax;
	       if(simrecdeltaphi > 2*PI) simrecdeltaphi -= 2*PI;
	       if(simrecdeltaphi > PI) simrecdeltaphi = PI -simrecdeltaphi;
               if(fabs(simhitPos.eta()-eta)>matchEtaMax) continue;
               if(simrecdeltaphi>matchPhiMax) continue;

	       int associatedTPID = associateSimhitToTrackingparticle((*simHit1).trackId());
	       ptype = (&(*simHit1))->processType();

	       if (associatedTPID == -1){
	          isbackground = true;
	          continue;    // doesn't match to any Trackingparticle
	       }
	    
	       const TrackingParticle* associatedTP = &(*trackingParticles)[associatedTPID];
	    
	       TrackingParticle::genp_iterator itb = associatedTP->genParticle_begin();
	       TrackingParticle::genp_iterator itend = associatedTP->genParticle_end();

	       isprimary = checkprimaryparticle(associatedTP);
               issecondary = itb == itend; 

	       if(itb == itend){
	          //	       cout<<"This is a secondary particle"<<endl;
	       }

	       //  	    cout<<"TP eta : "<<associatedTP->eta()<<" phi : "<<associatedTP->phi()<<endl;

	       for(TrackingParticle::genp_iterator itp = itb; itp != itend; ++itp){
	          gpid = tpmap_[(*itp)->barcode()];
	          //	       cout<<" Particle : "<<gpid<<endl;
	       }
	       
	       if (isprimary && bestSimHit1==0){ 
	          bestSimHit1 = &(*simHit1);
	          break;
	       }
            } 
         
            if(bestSimHit1!=0){
	       trid = bestSimHit1->trackId();  

               GlobalPoint simpos = pixelLayer->toGlobal((*bestSimHit1).localPosition());
            }
         }
         */
/*
         int type=-99;
	 
	 type = -99;
	 
	 if(isbackground) type = 0;
	 if(isprimary) type = 1;
	 if(ptype != 2) type = 2;
         if(issecondary) type = 3;
*/
	
	 if(layer == 1 && subid==1){ 
	    pev_.eta1[pev_.nhits1] = eta;
	    pev_.phi1[pev_.nhits1] = phi;
	    pev_.r1[pev_.nhits1] = r;
	    //pev_.id1[pev_.nhits1] = trid;
	    pev_.cs1[pev_.nhits1] = recHit->cluster()->size(); //Cluster Size
            pev_.ch1[pev_.nhits1] = recHit->cluster()->charge(); //Cluster Charge
	    //pev_.gp1[pev_.nhits1] = gpid;
	    //pev_.type1[pev_.nhits1] = type;
	    pev_.nhits1++;
	    if(fabs(gpos.eta()) < etaMult_ ) pev_.mult++;
	 }
	 
	 if(layer == 2 && subid==1){
	    pev_.eta2[pev_.nhits2] = eta;
	    pev_.phi2[pev_.nhits2] = phi;
	    pev_.r2[pev_.nhits2] = r;
            //pev_.id2[pev_.nhits2] = trid;
	    pev_.cs2[pev_.nhits2] = recHit->cluster()->size(); //Cluster Size
            pev_.ch2[pev_.nhits2] = recHit->cluster()->charge(); //Cluster Charge
	    //pev_.gp2[pev_.nhits2] = gpid;
	    //pev_.type2[pev_.nhits2] = type;
	    pev_.nhits2++;
	 } 

	 if(layer == 3 && subid==1){
	    pev_.eta3[pev_.nhits3] = eta;
	    pev_.phi3[pev_.nhits3] = phi;
	    pev_.r3[pev_.nhits3] = r;
            //pev_.id3[pev_.nhits3] = trid;
	    pev_.cs3[pev_.nhits3] = recHit->cluster()->size(); //Cluster Size
            pev_.ch3[pev_.nhits3] = recHit->cluster()->charge(); //Cluster Charge
	    //pev_.gp3[pev_.nhits3] = gpid;
	    //pev_.type3[pev_.nhits3] = type;
	    pev_.nhits3++;
	 } 
	 if(disk == 1 && subid==2){
	    pev_.eta4[pev_.nhits4] = eta;
	    pev_.phi4[pev_.nhits4] = phi;
	    pev_.r4[pev_.nhits4] = r;
            //pev_.id4[pev_.nhits4] = trid;
	    pev_.cs4[pev_.nhits4] = recHit->cluster()->size(); //Cluster Size
            pev_.ch4[pev_.nhits4] = recHit->cluster()->charge(); //Cluster Charge
	    //pev_.gp4[pev_.nhits4] = gpid;
	    //pev_.type4[pev_.nhits4] = type;
	    pev_.nhits4++;
	 } 
	 if(disk == 2 && subid==2){
	    pev_.eta5[pev_.nhits5] = eta;
	    pev_.phi5[pev_.nhits5] = phi;
	    pev_.r5[pev_.nhits5] = r;
            //pev_.id5[pev_.nhits5] = trid;
	    pev_.cs5[pev_.nhits5] = recHit->cluster()->size(); //Cluster Size
            pev_.ch5[pev_.nhits5] = recHit->cluster()->charge(); //Cluster Charge
	    //pev_.gp5[pev_.nhits5] = gpid;
	    //pev_.type5[pev_.nhits5] = type;
	    pev_.nhits5++;
	 } 
      }
   }
}

//--------------------------------------------------------------------------------------------------
void
PixelHitAnalyzer::fillParticles(const edm::Event& iEvent)
{
   // Handle<HepMCProduct> mc;
   // iEvent.getByLabel("generator",mc);
   // const HepMC::GenEvent* evt = mc->GetEvent();

   // int evtType = evt->signal_process_id();
   
   // pev_.evtType = evtType;
   // pev_.passDS = 0;
   // pev_.passSingleTrack = 0;
   
   // HepMC::GenEvent::particle_const_iterator begin = evt->particles_begin();
   // HepMC::GenEvent::particle_const_iterator end = evt->particles_end();
   
   // int nPos=0;
   // int nNeg=0;
   // for(HepMC::GenEvent::particle_const_iterator it = begin; it != end; ++it){
   //    if((*it)->status() != 1) continue;
   //       if ((*it)->momentum().eta()>3&&(*it)->momentum().eta()<5&&(*it)->momentum().e()>3) nPos++;
   //       if ((*it)->momentum().eta()>-5&&(*it)->momentum().eta()<-3&&(*it)->momentum().e()>3) nNeg++;
	 
   //       tpmap_[(*it)->barcode()] = pev_.nparticle;
   //       pev_.pdg[pev_.nparticle] = (*it)->pdg_id();
   //       pev_.eta[pev_.nparticle] = (*it)->momentum().eta();
   //       pev_.phi[pev_.nparticle] = (*it)->momentum().phi();
   //       pev_.pt[pev_.nparticle] = (*it)->momentum().perp();
   //       const ParticleData * part = pdt->particle(pev_.pdg[pev_.nparticle]);
   //       pev_.chg[pev_.nparticle] = (int)part->charge();
   //       if (fabs((*it)->momentum().eta())<0.8&&(*it)->momentum().perp()>0.4&&pev_.chg[pev_.nparticle]!=0) pev_.passSingleTrack=1;
   //       //pev_.x[pev_.nparticle] = (*it)->production_vertex()->position().x();
   //       //pev_.y[pev_.nparticle] = (*it)->production_vertex()->position().y();
   //       //pev_.z[pev_.nparticle] = (*it)->production_vertex()->position().z();
   //       if (pev_.chg[pev_.nparticle]==0) continue;
   //       if (fabs(pev_.eta[pev_.nparticle])>3) continue;
   //       pev_.nparticle++;
   // }
   
   // if (nPos>0&&nNeg>0) pev_.passDS=1;
   
   // // xi calculation
   // pev_.xi = -1;

   // vector<pair<double,HepMC::FourVector> > particles;
   // for(HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p!= evt->particles_end(); p++) {
   //   if((*p)->status() == 1) {
   //     double e  = (*p)->momentum().e();
   //     double pz = (*p)->momentum().pz();
   //     double y = 1./2*log((e+pz)/(e-pz));
   //     particles.push_back(pair<double,HepMC::FourVector> (y, (*p)->momentum()) );
   //   }
   // }

   // sort(particles.begin(), particles.end(), pairCompare);

   // double detamax = 0.;
   // vector<pair<double,HepMC::FourVector> >::const_iterator plast;
   // for(vector<pair<double,HepMC::FourVector> >::const_iterator p = particles.begin(); p!= particles.end() - 1; p++) {
   //   double deta = (p+1)->first - p->first;
   //   if(deta > detamax) { detamax = deta; plast = p; }
   // }

   // HepMC::FourVector sum_low, sum_hig;
   // for(vector<pair<double,HepMC::FourVector> >::iterator p = particles.begin(); p!= particles.end(); p++) {
   //   if(p <= plast) {
   //      sum_low.setPx(sum_low.px() + p->second.px());
   //      sum_low.setPy(sum_low.py() + p->second.py());
   //      sum_low.setPz(sum_low.pz() + p->second.pz());
   //      sum_low.setE (sum_low.e()  + p->second.e() );
   //   }
   //   else {
   //      sum_hig.setPx(sum_hig.px() + p->second.px());
   //      sum_hig.setPy(sum_hig.py() + p->second.py());
   //      sum_hig.setPz(sum_hig.pz() + p->second.pz());
   //      sum_hig.setE (sum_hig.e()  + p->second.e() );
   //   }
   // }

   // double m_low = sum_low.m();
   // double m_hig = sum_hig.m();

   // double Mx = (m_low > m_hig ? m_low : m_hig);

   // pev_.xi = Mx*Mx / (13e+3*13e+3); // 13 TeV

}
/*
//--------------------------------------------------------------------------------------------------
int PixelHitAnalyzer::associateSimhitToTrackingparticle(unsigned int trid )
{
   int ref=-1;
   const TrackingParticleCollection* TPCProd = trackingParticles.product();
   for (TrackingParticleCollection::size_type i=0; i<TPCProd->size(); i++){
      const TrackingParticle* tp = &(*TPCProd)[i];

      vector <PSimHit> particlesimhits = tp->trackPSimHit();
      for(vector<PSimHit>::const_iterator simhit = particlesimhits.begin(); simhit != particlesimhits.end(); ++simhit)
	 {
	    //cout <<"       matching TP: "<<i<<" TPsimhitid: "<<simhit->trackId()<<" simhitId: "<<trid<<endl;
	    if(simhit->trackId()==trid)//  checkprimaryparticle(tp))
	       {
		  ref=i;
		  break;
	       }
	 }
      if (ref!=-1) break;
   }

   return ref;
}
*/

//--------------------------------------------------------------------------------------------------   
bool PixelHitAnalyzer::checkprimaryparticle(const TrackingParticle* tp)
{
   int primarycheck=2;

   if(((tp->charge()==1)||(tp->charge()==-1))&&(tp->vertex().Rho()<0.2))
      {
	 primarycheck=1;
      } else {
	 primarycheck=0;
      }
   return primarycheck;
}       



// ------------ method called once each job just before starting event loop  ------------
void 
PixelHitAnalyzer::beginJob()
{

  pixelTree_ = fs->make<TTree>("PixelTree","Tree of Pixel Hits");
  pixelTree_->Branch("nEv",&pev_.nEv,"nEv/I");
  pixelTree_->Branch("nLumi",&pev_.nLumi,"nLumi/I");
  pixelTree_->Branch("nBX",&pev_.nBX,"nBX/I");
  pixelTree_->Branch("nRun",&pev_.nRun,"nRun/I");

  // Beam spot  
  pixelTree_->Branch("beamSpotX",&pev_.beamSpotX,"beamSpotX/F");
  pixelTree_->Branch("beamSpotY",&pev_.beamSpotY,"beamSpotY/F");
  pixelTree_->Branch("beamSpotZ",&pev_.beamSpotZ,"beamSpotZ/F");

  pixelTree_->Branch("nhits1",&pev_.nhits1,"nhits1/I");
  pixelTree_->Branch("nhits2",&pev_.nhits2,"nhits2/I");
  pixelTree_->Branch("nhits3",&pev_.nhits3,"nhits3/I");
  pixelTree_->Branch("nhits4",&pev_.nhits4,"nhits4/I");
  pixelTree_->Branch("nhits5",&pev_.nhits5,"nhits5/I");
  pixelTree_->Branch("ntrks",&pev_.ntrks,"ntrks/I");
  pixelTree_->Branch("ntrksCut",&pev_.ntrksCut,"ntrksCut/I");
  pixelTree_->Branch("mult",&pev_.mult,"mult/I");
  pixelTree_->Branch("nv",&pev_.nv,"nv/I");
  pixelTree_->Branch("vx",pev_.vx,"vx[nv]/F");
  pixelTree_->Branch("vy",pev_.vy,"vy[nv]/F");
  pixelTree_->Branch("vz",pev_.vz,"vz[nv]/F");
  pixelTree_->Branch("eta1",pev_.eta1,"eta1[nhits1]/F");
  pixelTree_->Branch("phi1",pev_.phi1,"phi1[nhits1]/F");
  pixelTree_->Branch("r1",pev_.r1,"r1[nhits1]/F");
  //pixelTree_->Branch("id1",pev_.id1,"id1[nhits1]/I");
  pixelTree_->Branch("cs1",pev_.cs1,"cs1[nhits1]/F");
  pixelTree_->Branch("ch1",pev_.ch1,"ch1[nhits1]/F");
  //pixelTree_->Branch("gp1",pev_.gp1,"gp1[nhits1]/I");
  //pixelTree_->Branch("type1",pev_.type1,"type1[nhits1]/I");

  pixelTree_->Branch("eta2",pev_.eta2,"eta2[nhits2]/F");
  pixelTree_->Branch("phi2",pev_.phi2,"phi2[nhits2]/F");
  pixelTree_->Branch("r2",pev_.r2,"r2[nhits2]/F");
  //pixelTree_->Branch("id2",pev_.id2,"id2[nhits2]/I");
  pixelTree_->Branch("cs2",pev_.cs2,"cs2[nhits2]/F");
  pixelTree_->Branch("ch2",pev_.ch2,"ch2[nhits2]/F");
  //pixelTree_->Branch("gp2",pev_.gp2,"gp2[nhits2]/I");
  //pixelTree_->Branch("type2",pev_.type2,"type2[nhits2]/I");

  pixelTree_->Branch("eta3",pev_.eta3,"eta3[nhits3]/F");
  pixelTree_->Branch("phi3",pev_.phi3,"phi3[nhits3]/F");
  pixelTree_->Branch("r3",pev_.r3,"r3[nhits3]/F");
  //pixelTree_->Branch("id3",pev_.id3,"id3[nhits3]/I");
  pixelTree_->Branch("cs3",pev_.cs3,"cs3[nhits3]/F");
  pixelTree_->Branch("ch3",pev_.ch3,"ch3[nhits3]/F");
  //pixelTree_->Branch("gp3",pev_.gp3,"gp3[nhits3]/I");
  //pixelTree_->Branch("type3",pev_.type3,"type3[nhits3]/I");

  pixelTree_->Branch("eta4",pev_.eta4,"eta4[nhits4]/F");
  pixelTree_->Branch("phi4",pev_.phi4,"phi4[nhits4]/F");
  pixelTree_->Branch("r4",pev_.r4,"r4[nhits4]/F");
  //pixelTree_->Branch("id4",pev_.id4,"id4[nhits4]/I");
  pixelTree_->Branch("cs4",pev_.cs4,"cs4[nhits4]/F");
  pixelTree_->Branch("ch4",pev_.ch4,"ch4[nhits4]/F");
  //pixelTree_->Branch("gp4",pev_.gp4,"gp4[nhits4]/I");
  //pixelTree_->Branch("type4",pev_.type4,"type4[nhits4]/I");

  pixelTree_->Branch("eta5",pev_.eta5,"eta5[nhits5]/F");
  pixelTree_->Branch("phi5",pev_.phi5,"phi5[nhits5]/F");
  pixelTree_->Branch("r5",pev_.r5,"r5[nhits5]/F");
  //pixelTree_->Branch("id5",pev_.id5,"id5[nhits5]/I");
  pixelTree_->Branch("cs5",pev_.cs5,"cs5[nhits5]/F");
  pixelTree_->Branch("ch5",pev_.ch5,"ch5[nhits5]/F");
  //pixelTree_->Branch("gp5",pev_.gp5,"gp5[nhits5]/I");
  //pixelTree_->Branch("type5",pev_.type5,"type5[nhits5]/I");


  pixelTree_->Branch("evtType",&pev_.evtType,"evtType/I");
  pixelTree_->Branch("npart",&pev_.nparticle,"npart/I");
  pixelTree_->Branch("pt",pev_.pt,"pt[npart]/F");
  pixelTree_->Branch("eta",pev_.eta,"eta[npart]/F");
  pixelTree_->Branch("phi",pev_.phi,"phi[npart]/F");
  pixelTree_->Branch("pdg",pev_.pdg,"pdg[npart]/I");
  pixelTree_->Branch("chg",pev_.chg,"chg[npart]/I");
  pixelTree_->Branch("xi",&pev_.xi,"xi/F");
  pixelTree_->Branch("passDS",&pev_.passDS,"passDS/O");
  pixelTree_->Branch("passSingleTrack",&pev_.passSingleTrack,"passSingleTrack/O");

  /* Not needed anymore 
  pixelTree_->Branch("x",pev_.x,"x[npart]/F");
  pixelTree_->Branch("y",pev_.y,"y[npart]/F");
  pixelTree_->Branch("z",pev_.z,"z[npart]/F");
  */

  pixelTree_->Branch("nHLTBit",&pev_.nHLTBit,"nHLTBit/I");
  pixelTree_->Branch("hltBit",pev_.hltBit,"hltBit[nHLTBit]/O");

  pixelTree_->Branch("nL1T",&pev_.nL1TBit,"nL1T/I");
  pixelTree_->Branch("L1T",pev_.l1TBit,"L1T[nL1T]/O");

  pixelTree_->Branch("nL1A",&pev_.nL1ABit,"nL1A/I");
  pixelTree_->Branch("L1A",pev_.l1ABit,"L1A[nL1A]/O");

  // HI related
  pixelTree_->Branch("zdcp",&pev_.zdcp,"zdcp/F");
  pixelTree_->Branch("zdcm",&pev_.zdcm,"zdcm/F");
  pixelTree_->Branch("hf",&pev_.hf,"hf/F");
  pixelTree_->Branch("hftp",&pev_.hftp,"hftp/F");
  pixelTree_->Branch("hftm",&pev_.hftm,"hftm/F");
  pixelTree_->Branch("hfp",&pev_.hfp,"hfp/F");
  pixelTree_->Branch("hfm",&pev_.hfm,"hfm/F");
  pixelTree_->Branch("etMid",&pev_.etMid,"etMid/F");
  pixelTree_->Branch("eb",&pev_.eb,"eb/F");
  pixelTree_->Branch("eep",&pev_.eep,"eep/F");
  pixelTree_->Branch("eem",&pev_.eem,"eem/F");
  pixelTree_->Branch("cBin",&pev_.cBin,"cBin/I");
  pixelTree_->Branch("nbins",&pev_.nbins,"nbins/I"); 
  pixelTree_->Branch("binsize",&pev_.binsize,"binsize/I");
  pixelTree_->Branch("nparti",&pev_.nparti,"nparti/F");
  pixelTree_->Branch("npartiSigma",&pev_.npartiSigma,"npartiSigma/F");
  pixelTree_->Branch("ncoll",&pev_.ncoll,"ncoll/F");
  pixelTree_->Branch("ncollSigma",&pev_.ncollSigma,"ncollSigma/F");
  pixelTree_->Branch("nhard",&pev_.nhard,"nhard/F");
  pixelTree_->Branch("nhardSigma",&pev_.nhardSigma,"nhardSigma/F");
  pixelTree_->Branch("b",&pev_.b,"b/F");
  pixelTree_->Branch("bSigma",&pev_.bSigma,"bSigma/F");
  pixelTree_->Branch("pixel",&pev_.pixel,"pixel/F");
  pixelTree_->Branch("nHFp",&pev_.nHFp,"nHFp/I");
  pixelTree_->Branch("nHFn",&pev_.nHFn,"nHFn/I");

  HLTConfigProvider hltConfig;
  /*
  cout <<"Configure hlt"<<endl;
  string teststr;

  // setup "Any" bit
  hltTrgBits_.clear();
  hltTrgBits_.push_back(-1);
  hltTrgDeci_.clear();
  hltTrgDeci_.push_back(true);
  hltTrgUsedNames_.clear();
  hltTrgUsedNames_.push_back("Any");

  // figure out relation of trigger name to trigger bit and store used trigger names/bits
  for(size_t i=0;i<hltTrgNames_.size();++i) {
    const string &n1(hltTrgNames_.at(i));
    bool found = 0;
    for(size_t j=0;j<hltConfig.size();++j) {
      const string &n2(hltConfig.triggerName(j));
      if (n2==n1) {
        hltTrgBits_.push_back(j);
        hltTrgUsedNames_.push_back(n1);
        hltTrgDeci_.push_back(false);
        cout <<Form("Added trigger %d with name %s for bit %d", 
                     hltTrgBits_.size()-1, n1.c_str(), j)<<endl;
        found = 1;
        break;
      }
    }      
    if (!found) {
      cout <<Form("Could not find trigger bit for %s", n1.c_str())<<endl;
    }
  }

  // ensure that trigger collections are of same size
  if (hltTrgBits_.size()!=hltTrgUsedNames_.size())
    cout <<Form("Size of trigger bits not equal used names: %d %d",
                 hltTrgBits_.size(), hltTrgUsedNames_.size())<<endl;
  if (hltTrgDeci_.size()!=hltTrgUsedNames_.size())
    cout <<Form("Size of decision bits not equal names: %d %d",
                 hltTrgDeci_.size(), hltTrgUsedNames_.size())<<endl;

  */
}

//--------------------------------------------------------------------------------------------------
void
PixelHitAnalyzer::fillPixelTracks(const edm::Event& iEvent){
/*
  // First fish the pixel tracks out of the event
  edm::Handle<reco::TrackCollection> trackCollection;
  edm::InputTag trackCollName = trackSrc_; 
  iEvent.getByLabel(trackCollName,trackCollection);
  const reco::TrackCollection tracks = *(trackCollection.product());
  reco::TrackRefVector trks;
  for (unsigned int i=0; i<tracks.size(); i++) {
    if (tracks[i].pt() > 0.15)     
      trks.push_back( reco::TrackRef(trackCollection, i) );
  }
  pev_.ntrksCut = trks.size();  
  pev_.ntrks = tracks.size();  
*/
}


//--------------------------------------------------------------------------------------------------
void PixelHitAnalyzer::fillL1Bits(const edm::Event &iEvent)
{
  // edm::Handle< L1GlobalTriggerReadoutRecord >  L1GlobalTrigger;

  // iEvent.getByLabel(L1gtReadout_, L1GlobalTrigger);
  // const TechnicalTriggerWord&  technicalTriggerWordBeforeMask = L1GlobalTrigger->technicalTriggerWord();

  // for (int i=0; i<64;i++)
  // {
  //   pev_.l1TBit[i] = technicalTriggerWordBeforeMask.at(i);
  // }
  // pev_.nL1TBit = 64;

  // int ntrigs = L1GlobalTrigger->decisionWord().size();
  // pev_.nL1ABit = ntrigs;

  // for (int i=0; i != ntrigs; i++) {
  //   bool accept = L1GlobalTrigger->decisionWord()[i];
  //   pev_.l1ABit[i] = (accept == true)? 1:0;
  // }
}

//--------------------------------------------------------------------------------------------------
void PixelHitAnalyzer::fillHLTBits(const edm::Event &iEvent)
{
  // // Fill HLT trigger bits.
  // Handle<TriggerResults> triggerResultsHLT;
  // iEvent.getByLabel(hltResName_,triggerResultsHLT);

  // const TriggerResults *hltResults = triggerResultsHLT.product();
  // const TriggerNames & triggerNames = iEvent.triggerNames(*hltResults);

  // pev_.nHLTBit = triggerNames.size();

  // for(size_t i=0;i<triggerNames.size();++i) {
  //     cout <<triggerNames.triggerName(i)<<endl;
  //     pev_.hltBit[i] = triggerResultsHLT->accept(i);
  // }
}

//--------------------------------------------------------------------------------------------------

void PixelHitAnalyzer::fillCentrality(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<reco::Centrality> centrality;
  iEvent.getByLabel(CentralityTag_, centrality);

/*
  hiNpix = centrality->multiplicityPixel();
  hiNpixelTracks = centrality->NpixelTracks();
  hiNtracks = centrality->Ntracks();
  hiNtracksPtCut = centrality->NtracksPtCut();
  hiNtracksEtaCut = centrality->NtracksEtaCut();
  hiNtracksEtaPtCut = centrality->NtracksEtaPtCut();

  hiHF = centrality->EtHFtowerSum();
  hiHFplus = centrality->EtHFtowerSumPlus();
  hiHFminus = centrality->EtHFtowerSumMinus();
  hiHFplusEta4 = centrality->EtHFtruncatedPlus();
  hiHFminusEta4 = centrality->EtHFtruncatedMinus();
  hiHFhit = centrality->EtHFhitSum();
  hiHFhitPlus = centrality->EtHFhitSumPlus();
  hiHFhitMinus = centrality->EtHFhitSumMinus();

  hiZDC = centrality->zdcSum();
  hiZDCplus = centrality->zdcSumPlus();
  hiZDCminus = centrality->zdcSumMinus();

  hiEEplus = centrality->EtEESumPlus();
  hiEEminus = centrality->EtEESumMinus();
  hiEE = centrality->EtEESum();
  hiEB = centrality->EtEBSum();
  hiET = centrality->EtMidRapiditySum();

  centrality_->newEvent(iEvent,iSetup);
*/
  pev_.hf = centrality->EtHFhitSum();
  pev_.hftp = (double)centrality->EtHFtowerSumPlus();
  pev_.hftm = (double)centrality->EtHFtowerSumMinus();
  pev_.zdcp = (double)centrality->zdcSumPlus();
  pev_.zdcm = (double)centrality->zdcSumMinus();
 
  pev_.hfp = (double)centrality->EtHFhitSumPlus();
  pev_.hfm = (double)centrality->EtHFhitSumMinus();
  pev_.etMid = (double)centrality->EtMidRapiditySum();
  pev_.eb = (double)centrality->EtEBSum();
  pev_.eep = (double)centrality->EtEESumPlus();
  pev_.eem = (double)centrality->EtEESumMinus();
/*
  pev_.cBin = (int)centrality->getBin();
  pev_.nbins = (int)centrality->getNbins(); 
  pev_.binsize = (int)(100/centrality->getNbins() );
  pev_.nparti = (double)centrality->NpartMean();
  pev_.npartiSigma = (double)centrality->NpartSigma();
  pev_.ncoll = (double)centrality->NcollMean();
  pev_.ncollSigma = (double)centrality->NcollSigma();
  pev_.nhard = (double)centrality->NhardMean();
  pev_.nhardSigma = (double)centrality->NhardSigma();
  pev_.b = (double)centrality->bMean();
  pev_.bSigma = (double)centrality->bSigma();
*/
  pev_.pixel = (double)centrality->multiplicityPixel();

}


void PixelHitAnalyzer::fillHF(const edm::Event& iEvent)
{

 edm::Handle<CaloTowerCollection> towers;
 iEvent.getByToken(TowerSrc_,towers);
 //iEvent.getByLabel(TowerSrc_,towers);

 int negTowers = 0;
 int posTowers = 0;
 for(CaloTowerCollection::const_iterator cal = towers->begin(); cal != towers->end(); ++cal) {
    for(unsigned int i = 0; i < cal->constituentsSize(); i++) {
       const DetId id = cal->constituent(i);
       if(id.det() == DetId::Hcal) {
         HcalSubdetector subdet=(HcalSubdetector(id.subdetId()));
         if(subdet == HcalForward) {
           if(cal->energy()>3. && cal->eta()<-3.)
             negTowers++;
           if(cal->energy()>3. && cal->eta()>3.)
             posTowers++;
         }
       }
    }
 }
 
 pev_.nHFp = posTowers;
 pev_.nHFn = negTowers;
 
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PixelHitAnalyzer::endJob() {
}

// --------------------------------------------------------------------------------------------------
bool
PixelHitAnalyzer::pairCompare(pair<double,HepMC::FourVector> i, pair<double,HepMC::FourVector> j)
{
  return (i.first < j.first);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PixelHitAnalyzer);
