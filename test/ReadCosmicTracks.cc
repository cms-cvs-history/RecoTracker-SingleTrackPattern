// File: ReadSeeds.cc
// Description:  see ReadSeeds.h
// Author:  O. Gutsche
// Creation Date:  OGU Aug. 1 2005 Initial version.
//
//--------------------------------------------
#include <memory>
#include <string>
#include <iostream>

#include "RecoTracker/SingleTrackPattern/test/ReadCosmicTracks.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "Geometry/Vector/interface/GlobalPoint.h"
#include "Geometry/Vector/interface/GlobalVector.h"
#include "Geometry/Vector/interface/LocalVector.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
using namespace std;
using namespace edm;
ReadCosmicTracks::ReadCosmicTracks(edm::ParameterSet const& conf) : 
  conf_(conf)
{
}
void ReadCosmicTracks::beginJob(const edm::EventSetup& c){
  hFile = new TFile ( "ptRes.root", "RECREATE" );
  hptres = new TH1F("hptres","Pt resolution",100,-0.1,0.1);
  hptres1 = new TH1F("hptres1","Pt resolution  0-10 GeV",100,-0.005,0.005);
  hptres2 = new TH1F("hptres2","Pt resolution 10-15 GeV",100,-0.005,0.005);
  hptres3 = new TH1F("hptres3","Pt resolution 15-20 GeV",100,-0.005,0.005);
  hptres4 = new TH1F("hptres4","Pt resolution 20-25 GeV",100,-0.005,0.005);
  hptres5 = new TH1F("hptres5","Pt resolution 25-30 GeV",100,-0.005,0.005);
  hptres6 = new TH1F("hptres6","Pt resolution 30-35 GeV",100,-0.005,0.005);
  hptres7 = new TH1F("hptres7","Pt resolution 35-40 GeV",100,-0.005,0.005);
  hptres8 = new TH1F("hptres8","Pt resolution 40-45 GeV",100,-0.005,0.005);
  hptres9 = new TH1F("hptres9","Pt resolution 45-50 GeV",100,-0.005,0.005);
  hptres10 = new TH1F("hptres10","Pt resolution >50 GeV",100,-0.005,0.005);
  hrespt= new TH1F("hrespt","resolution as a function of Pt", 10, 5, 55);
  heffpt= new TH1F("heffpt","efficiency as a function of Pt", 10, 5, 55);
  heffhit= new TH1F("heffhit","efficiency as a function of number of hits", 25, 4.5, 29.5);
  for (uint ik=0;ik<10;ik++){
    inum[ik]=0;
    iden[ik]=0;
  }
  for (uint ik=0;ik<30;ik++){
    inum2[ik]=0;
    iden2[ik]=0;
  }
}
// Virtual destructor needed.
ReadCosmicTracks::~ReadCosmicTracks() {  }  

// Functions that gets called by framework every event
void ReadCosmicTracks::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  using namespace edm;
  //  std::cout<<"EV "<<e.id()<<std::endl; 
  // Step A: Get Inputs 
  
  bool trackable_cosmic=false;
  //  const reco::TrackCollection tC = *(trackCollection.product());
  //
  theStripHits.clear();
  edm::Handle<edm::PSimHitContainer> TOBHitsLowTof;
  edm::Handle<edm::PSimHitContainer> TOBHitsHighTof;
  edm::Handle<edm::PSimHitContainer> TIBHitsLowTof;
  edm::Handle<edm::PSimHitContainer> TIBHitsHighTof;
  //
  e.getByLabel("SimG4Object","TrackerHitsTOBLowTof", TOBHitsLowTof);
  e.getByLabel("SimG4Object","TrackerHitsTOBHighTof", TOBHitsHighTof);
  e.getByLabel("SimG4Object","TrackerHitsTIBLowTof", TIBHitsLowTof);
  e.getByLabel("SimG4Object","TrackerHitsTIBHighTof", TIBHitsHighTof);

  theStripHits.insert(theStripHits.end(), TOBHitsLowTof->begin(), TOBHitsLowTof->end());
  theStripHits.insert(theStripHits.end(), TOBHitsHighTof->begin(), TOBHitsHighTof->end());
  theStripHits.insert(theStripHits.end(), TIBHitsLowTof->begin(), TIBHitsLowTof->end());
  theStripHits.insert(theStripHits.end(), TIBHitsHighTof->begin(), TIBHitsHighTof->end());
  edm::ESHandle<TrackerGeometry> tracker;
  es.get<TrackerDigiGeometryRecord>().get(tracker);
  uint nshit=theStripHits.size();
  edm::Handle<TrajectorySeedCollection> seedcoll;
  e.getByType(seedcoll);
  bool seed_plus= false;
  if( (seedcoll.product()->size()>0) &&(nshit>4)) {
    trackable_cosmic=true;
    unsigned int iraw=(*(*(*seedcoll).begin()).recHits().first).geographicalId().rawId();
    LocalPoint lp=(*(*(*seedcoll).begin()).recHits().first).localPosition();
    seed_plus=(tracker->idToDet(DetId(iraw))->surface().toGlobal(lp).y()>0.);
  }
  if (nshit>30) nshit=30;

  
  stable_sort(theStripHits.begin(),theStripHits.end(),CompareTOF());
  PSimHit isimfirst=(*theStripHits.begin());
  DetId tmp1=DetId(isimfirst.detUnitId());
  GlobalVector gvs= tracker->idToDet(tmp1)->surface().toGlobal(isimfirst.localDirection());
  float ptsims=gvs.perp()*isimfirst.pabs();
  unsigned int iptsims= uint(ptsims/5 -1);
  if (iptsims>10) iptsims=10;  

  if (trackable_cosmic){
    iden[iptsims]++;
    iden2[nshit]++;
  }


  edm::Handle<reco::TrackCollection> trackCollection;
  //    event.getByLabel("trackp", trackCollection);
  e.getByLabel("cosmictrackfinder",trackCollection);
  

  edm::Handle<TrackingRecHitCollection> trackrechitCollection;
  //    event.getByLabel("trackp", trackCollection);
  e.getByType(trackrechitCollection);
  

  const   reco::TrackCollection *tracks=trackCollection.product();
  if (tracks->size()>0){
    reco::TrackCollection::const_iterator ibeg=tracks->begin();
    if (trackable_cosmic) {
      inum[iptsims]++;
      inum2[nshit]++;
    }
    
    GlobalPoint gp((*ibeg).outerPosition().x(),
		   (*ibeg).outerPosition().y(),
		   (*ibeg).outerPosition().z());
  
    float ptrec = (*ibeg).outerPt();
  


    //    cout<<"PTREC "<<ptrec<<endl;


    std::vector<PSimHit>::iterator ihit;
    
    float magmag=10000;
    PSimHit isim;
    DetId idet;
    for (ihit=theStripHits.begin();ihit!=theStripHits.end();ihit++){
      DetId tmp=DetId((*ihit).detUnitId());
      GlobalPoint gp1 =tracker->idToDet(tmp)->surface().toGlobal((*ihit).localPosition());
      
      if ((gp1-gp).mag()<magmag){
	magmag=(gp1-gp).mag();
	isim=(*ihit);
	idet=DetId(tmp);
      }
    }
  

    //    GlobalPoint max =tracker->idToDet(idet)->surface().toGlobal(isim.localPosition());
    GlobalVector gv= tracker->idToDet(idet)->surface().toGlobal(isim.localDirection());
    float PP=isim.pabs();
    float ptsim=gv.perp()*PP;
    //    cout<<"PTSIM "<<ptsim<<endl;
    float ptresrel=((1./ptrec)-(1./ptsim));
    if (seed_plus){  
      hptres->Fill(ptresrel);
      unsigned int iptsim= uint(ptsim/5 -1);
      if (iptsim>10) iptsim=10;
      
      if (iptsim==1) hptres1->Fill(ptresrel);
      if (iptsim==2) hptres2->Fill(ptresrel);
      if (iptsim==3) hptres3->Fill(ptresrel);
      if (iptsim==4) hptres4->Fill(ptresrel);
      if (iptsim==5) hptres5->Fill(ptresrel);
      if (iptsim==6) hptres6->Fill(ptresrel);
      if (iptsim==7) hptres7->Fill(ptresrel);
      if (iptsim==8) hptres8->Fill(ptresrel);
      if (iptsim==9) hptres9->Fill(ptresrel);
      if (iptsim==10) hptres10->Fill(ptresrel);
    }
  
  }
}
void ReadCosmicTracks::endJob(){
  for  (uint ik=0;ik<10;ik++)  heffpt->Fill(7.5+(ik*5),float(inum[ik])/float(iden[ik]));
  for  (uint ik=0;ik<30;ik++)  heffhit->Fill(ik,float(inum2[ik])/float(iden2[ik]));
  hFile->Write();
  hFile->Close();
}
