
#include <memory>
#include <string>
#include <iostream>

#include "RecoTracker/SingleTrackPattern/test/AnalyzeMTCCTracks.h"

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
AnalyzeMTCCTracks::AnalyzeMTCCTracks(edm::ParameterSet const& conf) : 
  conf_(conf)
{
}
void AnalyzeMTCCTracks::beginJob(const edm::EventSetup& c){
  hFile = new TFile ( "trackhisto.root", "RECREATE" );
  hphi = new TH1F("hphi","Phi distribution",20,-3.14,3.14);
  hnhit = new TH1F("hnhit","Number of Hits per Track ",5,2.5,7.5);

}
// Virtual destructor needed.
AnalyzeMTCCTracks::~AnalyzeMTCCTracks() {  }  

// Functions that gets called by framework every event
void AnalyzeMTCCTracks::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  using namespace edm;
  // Step A: Get Inputs 
  
  edm::Handle<reco::TrackCollection> trackCollection;
  //    event.getByLabel("trackp", trackCollection);
  e.getByType(trackCollection);
  

 
  edm::ESHandle<TrackerGeometry> tracker;
  es.get<TrackerDigiGeometryRecord>().get(tracker);
 
  const   reco::TrackCollection *tracks=trackCollection.product();
 
  if (tracks->size()>0){
    reco::TrackCollection::const_iterator ibeg=trackCollection.product()->begin();
    hphi->Fill((*ibeg).outerPhi());
    hnhit->Fill((*ibeg).recHitsSize() );
  }
}
void AnalyzeMTCCTracks::endJob(){
  hFile->Write();
  hFile->Close();
}
