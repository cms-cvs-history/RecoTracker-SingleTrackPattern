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
using namespace std;
ReadCosmicTracks::ReadCosmicTracks(edm::ParameterSet const& conf) : 
  conf_(conf)
{
}
void ReadCosmicTracks::beginJob(const edm::EventSetup& c){
  hFile = new TFile ( "ptRes.root", "RECREATE" );
  hptres = new TH1F("hptres","Pt resolution",100,-2.,2.);
}
// Virtual destructor needed.
ReadCosmicTracks::~ReadCosmicTracks() {  }  

// Functions that gets called by framework every event
void ReadCosmicTracks::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  using namespace edm;
  std::cout<<"EV "<<e.id()<<std::endl; 
  // Step A: Get Inputs 
  
  edm::Handle<reco::TrackCollection> trackCollection;
  //    event.getByLabel("trackp", trackCollection);
  e.getByType(trackCollection);
  

  edm::Handle<TrackingRecHitCollection> trackrechitCollection;
  //    event.getByLabel("trackp", trackCollection);
  e.getByType(trackrechitCollection);
  
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
 
  const   reco::TrackCollection *tracks=trackCollection.product();
  reco::TrackCollection::const_iterator ibeg=tracks->begin();
 

  if (tracks->size()>0){
    TrackingRecHitCollection::const_iterator irec;
    GlobalPoint gp;
    cout<<endl;
    for(irec=trackrechitCollection.product()->begin();
	irec!=trackrechitCollection.product()->end();
	irec++)   gp =tracker->idToDet((*irec).geographicalId())
		    ->surface().toGlobal((*irec).localPosition());
    //   cout<<gp<<endl;
    

    float ptrec = (*ibeg).outerPt();
  




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
    hptres->Fill(((1./ptrec)-(1./ptsim))*ptsim);
    //   std::cout<<max<<" max "<<gv*PP<<" "<<isim.particleType()<<std::endl;	
  }
}
void ReadCosmicTracks::endJob(){
  hFile->Write();
  hFile->Close();
}
