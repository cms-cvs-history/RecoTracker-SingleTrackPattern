#include <memory>
#include <string>
#include <iostream>

#include "RecoTracker/SingleTrackPattern/test/AnalyzeMTCCTracks.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
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
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h" 
#include "TrackingTools/TrajectoryState/interface/BasicSingleTrajectoryState.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateWithArbitraryError.h"
#include "RecoTracker/SingleTrackPattern/test/TrajectoryMeasurementResidual.h"

using namespace std;
AnalyzeMTCCTracks::AnalyzeMTCCTracks(edm::ParameterSet const& conf) : 
  conf_(conf)
{
  trinevents=conf_.getParameter<bool>("TrajInEvents");
}
void AnalyzeMTCCTracks::beginJob(const edm::EventSetup& c){
  hFile = new TFile ( "trackhisto.root", "RECREATE" );
  hphi = new     TH1F("hphi",    "Phi distribution",         100,-3.14,  3.14);
  heta = new     TH1F("heta",    "Eta distribution",         100,  -5.,  5.   ); 
  hnhit = new    TH1F("hnhit",   "Number of Hits per Track ", 15,  2.5, 17.5  );
  hchi = new     TH1F("hchi",    "Chi squared of the track", 300,    0, 60    );
  hresTIB = new  TH1F("hresTIB", "TIB residual",             300,-1.5,1.5 );
  hresTOB = new  TH1F("hresTOB", "TOB residual",             300,-1.5,1.5 );
  hresTIB1 = new  TH1F("hresTIB1", "TIB residual1",             300,-1.5,1.5 );
  hresTOB1 = new  TH1F("hresTOB1", "TOB residual1",             300,-1.5,1.5 );
  hresTIB2 = new  TH1F("hresTIB2", "TIB residual2",             300,-1.5,1.5 );
  hresTOB2 = new  TH1F("hresTOB2", "TOB residual2",             300,-1.5,1.5 );
  hpt = new      TH1F("hpt"    , "Transv. moment",           100,  0.0,100    );
  hpx = new      TH1F("hpx"    , "Px",                       100, -50.0,50    );
  hpy = new      TH1F("hpy"    , "Py",                       100, -50.0,50    );
  hpz = new      TH1F("hpz"    , "Pz",                       100, -50.0,50    );
  hq = new       TH1F("hq"     , "charge",                     2, -1.5,  1.5  );
}
// Virtual destructor needed.
AnalyzeMTCCTracks::~AnalyzeMTCCTracks() {  }  

// Functions that gets called by framework every event
void AnalyzeMTCCTracks::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  //  std::cout<<"P"<<std::endl;

  using namespace edm;
  // Step A: Get Inputs 
  edm::InputTag seedTag = conf_.getParameter<edm::InputTag>("cosmicSeeds");
  edm::InputTag TkTag = conf_.getParameter<edm::InputTag>("cosmicTracks");
  edm::Handle<TrajectorySeedCollection> seedcoll;
  e.getByLabel(seedTag,seedcoll);

  edm::Handle<reco::TrackCollection> trackCollection;
  e.getByLabel(TkTag,trackCollection);

  edm::Handle<TrackingRecHitCollection> trackerchitCollection;
  e.getByLabel(TkTag,trackerchitCollection);
  edm::Handle<std::vector<Trajectory> > TrajectoryCollection;
  if (trinevents){
    e.getByLabel(TkTag,TrajectoryCollection);
  }
  edm::ESHandle<TrackerGeometry> tracker;
  es.get<TrackerDigiGeometryRecord>().get(tracker);
 
  const   reco::TrackCollection *tracks=trackCollection.product();

  if (tracks->size()>0){
    AnalHits(*trackerchitCollection);
  
   
    if (nlay>2){
      reco::TrackCollection::const_iterator ibeg=trackCollection.product()->begin();
      hphi->Fill((*ibeg).outerPhi());
      heta->Fill((*ibeg).outerEta());
      hnhit->Fill((*ibeg).recHitsSize() );
      hchi->Fill((*ibeg).chi2());
      hpt->Fill((*ibeg).pt());
      hpx->Fill((*ibeg).px());
      hpy->Fill((*ibeg).py());
      hpz->Fill((*ibeg).pz());
      hq->Fill((*ibeg).charge());

      if((*ibeg).chi2()<200){

	if (trinevents)
	  makeResiduals(*(TrajectoryCollection.product()->begin()));
	else
	  makeResiduals((*(*seedcoll).begin()),
			*trackerchitCollection,
			e,
			es);
      }
    }
  }
}
void AnalyzeMTCCTracks::endJob(){
  hFile->Write();
  hFile->Close();
}
void AnalyzeMTCCTracks::makeResiduals(const Trajectory traj){
  std::vector<TrajectoryMeasurement> TMeas=traj.measurements();
  vector<TrajectoryMeasurement>::iterator itm;
  for (itm=TMeas.begin();itm!=TMeas.end();itm++){
    TrajectoryMeasurementResidual*  TMR=new TrajectoryMeasurementResidual(*itm);
    uint iidd =(*itm).recHit()->detUnit()->geographicalId().rawId();
    StripSubdetector iid=StripSubdetector(iidd);
    unsigned int subid=iid.subdetId();

    if    (subid==  StripSubdetector::TIB) hresTIB->Fill(TMR->measurementXResidual());
    if    (subid==  StripSubdetector::TOB) hresTOB->Fill(TMR->measurementXResidual());
    int lay=(iidd>>16) & 0xF;  int sub=(iidd>>25)&0x7 ;
    if ((lay==1)&&(sub==3)) hresTIB1->Fill(TMR->measurementXResidual());
    if ((lay==2)&&(sub==3)) hresTIB2->Fill(TMR->measurementXResidual());
    if ((lay==1)&&(sub==5)) hresTOB1->Fill(TMR->measurementXResidual());
    if ((lay==2)&&(sub==5)) hresTOB2->Fill(TMR->measurementXResidual());
    delete TMR;
  }
}
void AnalyzeMTCCTracks::makeResiduals(const TrajectorySeed& seed,
				     const TrackingRecHitCollection &hits,
				     const edm::Event& e, 
				     const edm::EventSetup& es){



  //services
  es.get<IdealMagneticFieldRecord>().get(magfield);
  es.get<TrackerDigiGeometryRecord>().get(tracker);
  
  bool  seed_plus=(seed.direction()==alongMomentum);
  
  if (seed_plus) { 	 
    thePropagator=      new PropagatorWithMaterial(alongMomentum,0.1057,&(*magfield) ); 	 
    thePropagatorOp=    new PropagatorWithMaterial(oppositeToMomentum,0.1057,&(*magfield) );} 	 
  else {
    thePropagator=      new PropagatorWithMaterial(oppositeToMomentum,0.1057,&(*magfield) ); 	
    thePropagatorOp=    new PropagatorWithMaterial(alongMomentum,0.1057,&(*magfield) );
  }

  theUpdator=       new KFUpdator();
  theEstimator=     new Chi2MeasurementEstimator(30);
  
  
  edm::ESHandle<TransientTrackingRecHitBuilder> theBuilder;
  std::string builderName = conf_.getParameter<std::string>("TTRHBuilder");   
  es.get<TransientRecHitRecord>().get(builderName,theBuilder);  
  RHBuilder=   theBuilder.product();


  theFitter=        new KFTrajectoryFitter(*thePropagator,
					   *theUpdator,	
					   *theEstimator) ;
  theSmoother=      new KFTrajectorySmoother(*thePropagatorOp,
 					     *theUpdator,	
 					     *theEstimator);
 


  Trajectory traj=createStartingTrajectory(*(&seed));
  TransientTrackingRecHit::RecHitContainer trans_hits;
  for (unsigned int icosmhit=hits.size()-1;icosmhit+1>0;icosmhit--){
    TransientTrackingRecHit::RecHitPointer tmphit=RHBuilder->build(&(hits[icosmhit]));

    trans_hits.push_back(&(*tmphit));
    if (icosmhit<hits.size()-1){
      TSOS prSt= 
	thePropagator->propagate(traj.lastMeasurement().updatedState(),
				 tracker->idToDet(hits[icosmhit].geographicalId())->surface());

      if (prSt.isValid()){

	if(theUpdator->update( prSt, *tmphit).isValid()){
	  traj.push(TrajectoryMeasurement(prSt,
					  theUpdator->update( prSt, *tmphit),
					  tmphit,
					  theEstimator->estimate(prSt, *tmphit).second));
	}
      }
    }
  }

  if (thePropagatorOp->propagate(traj.lastMeasurement().updatedState(),
				 tracker->idToDet((*trans_hits.begin())->geographicalId())->surface()).isValid()){

    TSOS startingState=  TrajectoryStateWithArbitraryError()
      (thePropagatorOp->propagate(traj.lastMeasurement().updatedState(),
				  tracker->idToDet((*trans_hits.begin())->geographicalId())->surface()));

    std::vector<Trajectory> fittraj=theFitter->fit(seed,trans_hits,startingState);
    if (fittraj.size()>0){

      const Trajectory ifitted= *(fittraj.begin());

      std::vector<Trajectory> smoothtraj=theSmoother->trajectories(ifitted);
      if (smoothtraj.size()>0){
	const Trajectory smooth=*(smoothtraj.begin());


	std::vector<TrajectoryMeasurement> TMeas=smooth.measurements();
	vector<TrajectoryMeasurement>::iterator itm;
	for (itm=TMeas.begin();itm!=TMeas.end();itm++){
	  TrajectoryMeasurementResidual*  TMR=new TrajectoryMeasurementResidual(*itm);
	  uint iidd =(*itm).recHit()->detUnit()->geographicalId().rawId();
	  StripSubdetector iid=StripSubdetector(iidd);
	  unsigned int subid=iid.subdetId();
	  if    (subid==  StripSubdetector::TIB) hresTIB->Fill(TMR->measurementXResidual());
	  if    (subid==  StripSubdetector::TOB) hresTOB->Fill(TMR->measurementXResidual());
	  int lay=(iidd>>16) & 0xF;  int sub=(iidd>>25)&0x7 ;
	  if ((lay==1)&&(sub==3)) hresTIB1->Fill(TMR->measurementXResidual());
	  if ((lay==2)&&(sub==3)) hresTIB2->Fill(TMR->measurementXResidual());
	  if ((lay==1)&&(sub==5)) hresTOB1->Fill(TMR->measurementXResidual());
	  if ((lay==2)&&(sub==5)) hresTOB2->Fill(TMR->measurementXResidual());
	  delete TMR;
	}
      }
    }
  }
  delete thePropagator;
  delete thePropagatorOp;
  delete theUpdator;
  delete theEstimator;
  delete theFitter;
  delete theSmoother;


}
TrajectoryStateOnSurface
AnalyzeMTCCTracks::startingTSOS(const TrajectorySeed& seed)const
{
  PTrajectoryStateOnDet pState( seed.startingState());
  const GeomDet* gdet  = (&(*tracker))->idToDet(DetId(pState.detId()));
  TSOS  State= tsTransform.transientState( pState, &(gdet->surface()), 
					   &(*magfield));
  return State;

}
Trajectory AnalyzeMTCCTracks::createStartingTrajectory( const TrajectorySeed& seed) const
{
  Trajectory result( seed, seed.direction());
  std::vector<TrajectoryMeasurement> seedMeas = seedMeasurements(seed);
  if ( !seedMeas.empty()) {
    for (std::vector<TrajectoryMeasurement>::const_iterator i=seedMeas.begin(); i!=seedMeas.end(); i++){
      result.push(*i);
    }
  }
 
  return result;
}

void AnalyzeMTCCTracks::AnalHits(const TrackingRecHitCollection &hits){
  using namespace std;
  ltob1=false; ltob2=false; ltib1=false; ltib2=false;

  TrackingRecHitCollection::const_iterator hit;
  for(hit=hits.begin();hit!=hits.end();hit++){
  
    uint iid=(*hit).geographicalId().rawId();
    int sub=(iid>>25)&0x7 ;
    int lay=(iid>>16) & 0xF;
    if ((lay==1)&&(sub==3)) ltib1=true;
    if ((lay==2)&&(sub==3)) ltib2=true;
    if ((lay==1)&&(sub==5)) ltob1=true;
    if ((lay==2)&&(sub==5)) ltob2=true;
  
  }
  nlay=ltib1+ltib2+ltob1+ltob2;
 
}
std::vector<TrajectoryMeasurement> 
AnalyzeMTCCTracks::seedMeasurements(const TrajectorySeed& seed) const
{
  std::vector<TrajectoryMeasurement> result;
  TrajectorySeed::range hitRange = seed.recHits();
  for (TrajectorySeed::const_iterator ihit = hitRange.first; 
       ihit != hitRange.second; ihit++) {
    //RC TransientTrackingRecHit* recHit = RHBuilder->build(&(*ihit));
    TransientTrackingRecHit::RecHitPointer recHit = RHBuilder->build(&(*ihit));
    const GeomDet* hitGeomDet = (&(*tracker))->idToDet( ihit->geographicalId());
    TSOS invalidState( new BasicSingleTrajectoryState( hitGeomDet->surface()));

    if (ihit == hitRange.second - 1) {
      TSOS  updatedState=startingTSOS(seed);
      result.push_back(TrajectoryMeasurement( invalidState, updatedState, recHit));

    } 
    else {
      result.push_back(TrajectoryMeasurement( invalidState, recHit));
    }
    
  }

  return result;
};
