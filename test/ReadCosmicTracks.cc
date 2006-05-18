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

ReadCosmicTracks::ReadCosmicTracks(edm::ParameterSet const& conf) : 
  conf_(conf)
{
}

// Virtual destructor needed.
ReadCosmicTracks::~ReadCosmicTracks() { }  

// Functions that gets called by framework every event
void ReadCosmicTracks::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  using namespace edm;
  std::cout<<"EV "<<e.id()<<std::endl; 
  // Step A: Get Inputs 
  
  edm::Handle<reco::TrackCollection> trackCollection;
  //    event.getByLabel("trackp", trackCollection);
  e.getByType(trackCollection);
  
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
   reco::TrackCollection::const_iterator itr;
   const   reco::TrackCollection *tracks=trackCollection.product();
   reco::TrackCollection::const_iterator ibeg=tracks->begin();
   reco::TrackCollection::const_iterator iend=tracks->end();
//     //  reco::track_iterator ibeg=tracks->begin();
//   //  reco::track_iterator ibeg2=coll.begin();
//   //  std::cout <<" FOUND "<<(coll.product())->size()<<" Seeds."<<std::endl;

  for (itr=ibeg;itr!=iend;itr++){
    //  std::cout<<std::endl<<std::endl<<"EVENT ANAL "<<e.id()<<std::endl;
    //  std::cout<<"Mom of the track "<<(*itr).outerMomentum()<<std::endl;
  }

  std::vector<PSimHit>::iterator ihit;

  float ymin=-1000;
  PSimHit isim;
  DetId idet;
  for (ihit=theStripHits.begin();ihit!=theStripHits.end();ihit++){
    
    DetId tmp=DetId((*ihit).detUnitId());
    //    GlobalPoint gp = tracker->idToDet((*ihit).detUnitId().rawId())
    GlobalPoint gp =tracker->idToDet(tmp)->surface().toGlobal((*ihit).localPosition());
    std::cout<<"SIMHIT POSITION "<<gp<<std::endl;  
  //  GlobalVector gv= tracker->idToDet(pippo)->surface().toGlobal((*ihit).localDirection());
    if (gp.y()>ymin){
      ymin=gp.y();
      isim=(*ihit);
      idet=DetId(tmp);
    }
    // std::cout<<gp<<" "<<std::endl;
  }
   GlobalPoint max =tracker->idToDet(idet)->surface().toGlobal(isim.localPosition());
   //   LocalVector loc=isim.localDirection();
   GlobalVector gv= tracker->idToDet(idet)->surface().toGlobal(isim.localDirection());
   float PP=isim.pabs();
   //   if (tracks->size()>0)
     //   std::cout<<max<<" max "<<gv*PP<<std::endl;	
}
