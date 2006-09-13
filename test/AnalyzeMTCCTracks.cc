
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
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h" 
#include "TrackingTools/TrajectoryState/interface/BasicSingleTrajectoryState.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateWithArbitraryError.h"
#include "RecoTracker/SingleTrackPattern/test/TrajectoryMeasurementResidual.h"

using namespace std;
AnalyzeMTCCTracks::AnalyzeMTCCTracks(edm::ParameterSet const& conf) : 
  conf_(conf)
{
}
void AnalyzeMTCCTracks::beginJob(const edm::EventSetup& c){
  hFile = new TFile ( "trackhisto.root", "RECREATE" );
  hphi = new TH1F("hphi","Phi distribution",20,-3.14,3.14);
  hnhit = new TH1F("hnhit","Number of Hits per Track ",5,2.5,7.5);
  hchi = new TH1F("hchi","Chi squared of the track",100,0,100);
}
// Virtual destructor needed.
AnalyzeMTCCTracks::~AnalyzeMTCCTracks() {  }  

// Functions that gets called by framework every event
void AnalyzeMTCCTracks::analyze(const edm::Event& e, const edm::EventSetup& es)
{

  cout<<"b1"<<endl;
  using namespace edm;
  // Step A: Get Inputs 
    edm::Handle<TrajectorySeedCollection> seedcoll;
    e.getByType(seedcoll);

  edm::Handle<reco::TrackCollection> trackCollection;
  e.getByLabel("cosmictrackfinder",trackCollection);

  edm::Handle<TrackingRecHitCollection> trackrechitCollection;
  e.getByType(trackrechitCollection);

 
  edm::ESHandle<TrackerGeometry> tracker;
  es.get<TrackerDigiGeometryRecord>().get(tracker);
 
  const   reco::TrackCollection *tracks=trackCollection.product();
   cout<<"b2"<<endl;
  if (tracks->size()>0){
    reco::TrackCollection::const_iterator ibeg=trackCollection.product()->begin();
    hphi->Fill((*ibeg).outerPhi());
    hnhit->Fill((*ibeg).recHitsSize() );
    hchi->Fill((*ibeg).chi2());
    cout<<"b3"<<endl;
    makeResiduals((*(*seedcoll).begin()),
		  *trackrechitCollection,
		  e,
		  es);
    cout<<"b4"<<endl;
  }
}
void AnalyzeMTCCTracks::endJob(){
  hFile->Write();
  hFile->Close();
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
 
  TransientTrackingRecHit::RecHitPointer recHit = RHBuilder->build(&(*(seed.recHits().first)));
  const GeomDet* hitGeomDet = (&(*tracker))->idToDet( recHit->geographicalId());
  TSOS invalidState( new BasicSingleTrajectoryState( hitGeomDet->surface()));
  PTrajectoryStateOnDet pState( seed.startingState());
  TSOS  State= tsTransform.transientState( pState, &(hitGeomDet->surface()), 
					   &(*magfield));
  cout<<"a000"<<endl;
  //  const  CosmicTrajectoryBuilder *pippo= new  CosmicTrajectoryBuilder(conf_);
  cout<<"a001"<<endl;
  Trajectory traj=pippo->createStartingTrajectory(*(&seed));
  cout<<"a002"<<endl;
  //  Trajectory traj=CosmicTrajectoryBuilder::createStartingTrajectory(&seed);
  // Trajectory traj2=CosmicTrajectoryBuilder::createStartingTrajectory(*seed);
  // Trajectory traj=const CosmicTrajectoryBuilder::createStartingTrajectory(*(&seed));
  //  const Trajectory traj4=CosmicTrajectoryBuilder::createStartingTrajectory(*(&seed));
  // Trajectory traj( seed, seed.direction());
  // cout<<"a1"<<endl;
  //traj.push(TrajectoryMeasurement(State,recHit));
  cout<<"a2"<<endl;
  // traj.push(TrajectoryMeasurement(invalidState,recHit));
  cout<<"a3"<<endl;
 
  TransientTrackingRecHit::RecHitContainer trans_hits;
  for (unsigned int icosmhit=0;icosmhit<hits.size();icosmhit++){
 
    TransientTrackingRecHit::RecHitPointer tmphit=RHBuilder->build(&(hits[icosmhit]));
    // traj.push(TrajectoryMeasurement(invalidState,tmphit));
    cout<<"a4"<<endl;
    trans_hits.push_back(&(*tmphit));
    if (icosmhit>0){
      // traj.push(TrajectoryMeasurement(invalidState,tmphit));
      //  TSOS puppolo =traj.lastMeasurement().updatedState();
      // cout<<"a5"<<endl;
      TSOS prSt;
      if (icosmhit==1)
	prSt      = thePropagator->propagate(TrajectoryStateWithArbitraryError()(traj.lastMeasurement().updatedState()),
					     tracker->idToDet(hits[icosmhit].geographicalId())->surface());
      else  prSt=thePropagator->propagate(traj.lastMeasurement().updatedState(),
					  tracker->idToDet(hits[icosmhit].geographicalId())->surface());

      cout<<"a6"<<prSt.isValid()<<endl;
      //      TSOS UpdatedState= theUpdator->update( prSt, *tmphit);
     //   float chidue=theEstimator->estimate(prSt, *tmphit).second ;
     cout<<"a7"<<endl;
 //       TrajectoryMeasurement pippo(prSt,
// 				   theUpdator->update( prSt, *tmphit)
// 				   ,tmphit,chidue);
//       cout<<"aa "<<traj.isValid()<<" "<<prSt.isValid()<<" "<<UpdatedState.isValid()<<endl;
      //      if (traj.isValid())     traj.push(pippo);
     if (prSt.isValid()){
       cout<<"DF"<<endl;
     traj.push(TrajectoryMeasurement(prSt,
				     theUpdator->update( prSt, *tmphit),
				     tmphit,
				     theEstimator->estimate(prSt, *tmphit).second));
     }
     // traj.push(TrajectoryMeasurement(invalidState,tmphit));
    }
  }
//   TSOS startingState=  TrajectoryStateWithArbitraryError()
//     (thePropagatorOp->propagate(traj.lastMeasurement().updatedState(),
// 				tracker->idToDet((*trans_hits.begin())->geographicalId())->surface()));
//   const Trajectory ifitted= *(theFitter->fit(seed,trans_hits,startingState).begin());
//   const Trajectory smooth=*(theSmoother->trajectories(ifitted).begin());

//    std::vector<TrajectoryMeasurement> TMeas=smooth.measurements();
//    vector<TrajectoryMeasurement>::iterator itm;
//    for (itm=TMeas.begin();itm!=TMeas.end();itm++){
//      TrajectoryMeasurementResidual*  TMR=new TrajectoryMeasurementResidual(*itm);
// //     if    (subid==  StripSubdetector::TIB) hresTIB->Fill(r[0]);
// //     if    (subid==  StripSubdetector::TOB) hresTOB->Fill(r[0]);
// //     if    (subid==  StripSubdetector::TID) hresTID->Fill(r[0]);
// //     if    (subid==  StripSubdetector::TEC) hresTEC->Fill(r[0]);
// //     //cout<<r<<" "<<r[0]<<endl;
// //     //  cout<<"tm "<<(*itm).updatedState().globalPosition()<<" "<<(*itm).recHit()->globalPosition()<<endl;
//      delete TMR;
//  }


}
TrajectoryStateOnSurface
CosmicTrajectoryBuilder::startingTSOS(const TrajectorySeed& seed)const
{
  PTrajectoryStateOnDet pState( seed.startingState());
  const GeomDet* gdet  = (&(*tracker))->idToDet(DetId(pState.detId()));
  TSOS  State= tsTransform.transientState( pState, &(gdet->surface()), 
					   &(*magfield));
  return State;

}
Trajectory CosmicTrajectoryBuilder::createStartingTrajectory( const TrajectorySeed& seed) const
{
  Trajectory result( seed, seed.direction());
  std::vector<TM> seedMeas = seedMeasurements(seed);
  if ( !seedMeas.empty()) {
    for (std::vector<TM>::const_iterator i=seedMeas.begin(); i!=seedMeas.end(); i++){
      result.push(*i);
    }
  }
 
  return result;
}


std::vector<TrajectoryMeasurement> 
CosmicTrajectoryBuilder::seedMeasurements(const TrajectorySeed& seed) const
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
      result.push_back(TM( invalidState, updatedState, recHit));

    } 
    else {
      result.push_back(TM( invalidState, recHit));
    }
    
  }

  return result;
};
