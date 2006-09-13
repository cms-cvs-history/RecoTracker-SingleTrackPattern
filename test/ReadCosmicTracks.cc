// File: ReadCosmicTracks.cc
// Description:  see ReadCosmicTracks.h
// Author:  M.Pioppi Univ. & INFN Perugia
//--------------------------------------------
#include <memory>
#include <string>
#include <iostream>

#include "RecoTracker/SingleTrackPattern/test/ReadCosmicTracks.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "Geometry/Vector/interface/GlobalPoint.h"
#include "Geometry/Vector/interface/GlobalVector.h"
#include "Geometry/Vector/interface/LocalVector.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h" 
#include "RecoTracker/SingleTrackPattern/interface/CosmicTrajectoryBuilder.h"
#include "TrackingTools/PatternTools/interface/MeasurementExtractor.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "TrackingTools/TrajectoryState/interface/BasicSingleTrajectoryState.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateWithArbitraryError.h"
#include "RecoTracker/SingleTrackPattern/test/TrajectoryMeasurementResidual.h"
using namespace std;
using namespace edm;
ReadCosmicTracks::ReadCosmicTracks(edm::ParameterSet const& conf) : 
  conf_(conf)
{
}
void ReadCosmicTracks::beginJob(const edm::EventSetup& c){

  //HISTOGRAM BOOKING
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
  hchiSq = new TH1F("hchiSq","Chi2",50,0 ,300);
  hchiSq1 = new TH1F("hchiSq1","Chi2  0-10 GeV",50, 0 ,300);
  hchiSq2 = new TH1F("hchiSq2","Chi2 10-15 GeV",50, 0 ,300);
  hchiSq3 = new TH1F("hchiSq3","Chi2 15-20 GeV",50, 0 ,300);
  hchiSq4 = new TH1F("hchiSq4","Chi2 20-25 GeV",50, 0 ,300);
  hchiSq5 = new TH1F("hchiSq5","Chi2 25-30 GeV",50, 0 ,300);
  hchiSq6 = new TH1F("hchiSq6","Chi2 30-35 GeV",50, 0 ,300);
  hchiSq7 = new TH1F("hchiSq7","Chi2 35-40 GeV",50, 0 ,300);
  hchiSq8 = new TH1F("hchiSq8","Chi2 40-45 GeV",50, 0 ,300);
  hchiSq9 = new TH1F("hchiSq9","Chi2 45-50 GeV",50, 0 ,300);
  hchiSq10 = new TH1F("hchiSq10","Chi2 >50 GeV",50, 0 ,300);
  hchipt= new TH1F("hchipt","Chi2 as a function of Pt", 10, 5, 55);
  hchiR= new TH1F("hchiR","Chi2 as a function of R", 10, 5, 105);
  hrespt= new TH1F("hrespt","resolution as a function of Pt", 10, 5, 55);
  heffpt= new TH1F("heffpt","efficiency as a function of Pt", 10, 5, 55);
  heffhit= new TH1F("heffhit","efficiency as a function of number of hits", 9, 5, 50);
  hcharge=new TH1F("hcharge","1=SIM+REC+,2=SIM+REC-,3=SIM-REC-,4=SIM-REC+",4,0.5,4.5);
  hresTIB=new TH1F("hresTIB","residuals for TIB modules",100,-0.05,0.05);
  hresTOB=new TH1F("hresTOB","residuals for TOB modules",100,-0.05,0.05);
  hresTID=new TH1F("hresTID","residuals for TID modules",100,-0.05,0.05);
  hresTEC=new TH1F("hresTEC","residuals for TEC modules",100,-0.05,0.05);

  // COUNTERS INITIALIZATION
  for (uint ik=0;ik<10;ik++){
    inum[ik]=0;
    iden[ik]=0;
    ichiR[ik]=0;
    iden3[ik]=0;
  }
  for (uint ik=0;ik<9;ik++){
    inum2[ik]=0;
    iden2[ik]=0;
  }
}
// Virtual destructor needed.
ReadCosmicTracks::~ReadCosmicTracks() {  }  

// Functions that gets called by framework every event
void ReadCosmicTracks::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  
 
  
  bool trackable_cosmic=false;
  //SIM HITS
  theStripHits.clear();
  edm::Handle<edm::PSimHitContainer> TOBHitsLowTof;
  edm::Handle<edm::PSimHitContainer> TOBHitsHighTof;
  edm::Handle<edm::PSimHitContainer> TIBHitsLowTof;
  edm::Handle<edm::PSimHitContainer> TIBHitsHighTof;
  edm::Handle<edm::PSimHitContainer> PixelBarrelLowTof;
  edm::Handle<edm::PSimHitContainer> PixelBarrelHighTof;
  
  e.getByLabel("g4SimHits","TrackerHitsTOBLowTof", TOBHitsLowTof);
  e.getByLabel("g4SimHits","TrackerHitsTOBHighTof", TOBHitsHighTof);
  e.getByLabel("g4SimHits","TrackerHitsTIBLowTof", TIBHitsLowTof);
  e.getByLabel("g4SimHits","TrackerHitsTIBHighTof", TIBHitsHighTof);
  e.getByLabel("g4SimHits","TrackerHitsPixelBarrelLowTof", PixelBarrelLowTof);
  e.getByLabel("g4SimHits","TrackerHitsPixelBarrelHighTof", PixelBarrelHighTof);
  theStripHits.insert(theStripHits.end(), TOBHitsLowTof->begin(), TOBHitsLowTof->end());
  theStripHits.insert(theStripHits.end(), TOBHitsHighTof->begin(), TOBHitsHighTof->end());
  theStripHits.insert(theStripHits.end(), TIBHitsLowTof->begin(), TIBHitsLowTof->end());
  theStripHits.insert(theStripHits.end(), TIBHitsHighTof->begin(), TIBHitsHighTof->end());
  theStripHits.insert(theStripHits.end(), PixelBarrelLowTof->begin(), PixelBarrelLowTof->end());
  theStripHits.insert(theStripHits.end(), PixelBarrelHighTof->begin(), PixelBarrelHighTof->end());

  //RECONSTRUCTED OBJECTS
  //SEEDS AND TRACKS
  edm::Handle<TrajectorySeedCollection> seedcoll;
  e.getByType(seedcoll);

  edm::Handle<reco::TrackCollection> trackCollection;
  e.getByLabel("cosmictrackfinder",trackCollection);

  edm::Handle<TrackingRecHitCollection> trackrechitCollection;
  e.getByType(trackrechitCollection);

  //TRACKER GEOMETRY

  es.get<TrackerDigiGeometryRecord>().get(tracker);

  //SIMULATED INFORMATION:
  //NUMBER OF SIMHIT ,CHARGE  AND PT
  uint nshit=theStripHits.size();
  stable_sort(theStripHits.begin(),theStripHits.end(),CompareTOF());
  PSimHit isimfirst=(*theStripHits.begin());
  int chsim=-1*isimfirst.particleType()/abs(isimfirst.particleType());
  DetId tmp1=DetId(isimfirst.detUnitId());
  GlobalVector gvs= tracker->idToDet(tmp1)->surface().toGlobal(isimfirst.localDirection());
  float ptsims=gvs.perp()*isimfirst.pabs();


 
  //CRITERION TO SAY IF A MUON IS TRACKABLE OR NOT
  if( (seedcoll.product()->size()>0) &&(nshit>4)) {
    trackable_cosmic=true;
    seed_plus=((*(*seedcoll).begin()).direction()==alongMomentum);
   }

  if (nshit>50) nshit=50;
  uint inshit=(nshit/5)-1;
  unsigned int iptsims= uint(ptsims/5 -1);
  if (iptsims>10) iptsims=10;  

  //INCREASE THE NUMBER OF TRACKABLE MUONS
  if (trackable_cosmic){
    iden[iptsims]++;
    iden2[inshit]++;
  }
  
  //RECONSTRUCTED INFORMATION
  const   reco::TrackCollection *tracks=trackCollection.product();
  //INCREASE THE NUMBER OF TRACKED MUONS
  if (tracks->size()>0){
    reco::TrackCollection::const_iterator ibeg=tracks->begin();
    if (trackable_cosmic) {
      inum[iptsims]++;
      inum2[inshit]++;
    }
    
    GlobalPoint gp((*ibeg).outerPosition().x(),
		   (*ibeg).outerPosition().y(),
		   (*ibeg).outerPosition().z());

    //PT,Chi2,CHARGE RECONSTRUCTED
    float ptrec = (*ibeg).outerPt();
    float chiSquared =(*ibeg).chi2();
    int chrec=(*ibeg).charge();
    //cout<<"CHI1 "<< chiSquared<<endl;
    //FILL CHARGE HISTOGRAM 
    if(chrec==1 && chsim==1) hcharge->Fill(1);
    if(chrec==-1 && chsim==1) hcharge->Fill(2);
    if(chrec==-1 && chsim==-1) hcharge->Fill(3);
    if(chrec==1 && chsim==-1) hcharge->Fill(4);
  

    //FIND THE SIMULATED HIT CLOSEST TO THE 
    //FIRST REC HIT OF THE TRACK  
    //AND THE RADIAL DISTANCE OF THE TRACK
    std::vector<PSimHit>::iterator ihit;
    float MAGR=10000;
    float magmag=10000;
    DetId idet;
    for (ihit=theStripHits.begin();ihit!=theStripHits.end();ihit++){
      DetId tmp=DetId((*ihit).detUnitId());
      GlobalPoint gp1 =tracker->idToDet(tmp)->surface().toGlobal((*ihit).localPosition());

     float RR=sqrt((gp1.x()*gp1.x())+(gp1.y()*gp1.y()));
      if (RR<MAGR) MAGR=RR;
      if ((gp1-gp).mag()<magmag){
	magmag=(gp1-gp).mag();
	isim=&(*ihit);
	idet=DetId(tmp);
      }
    }
//     //CHI2 vs RADIAL DISTANCE
    uint iRR=uint(MAGR/10);
    iden3[iRR]++;
    ichiR[iRR]+=chiSquared;


    //REEVALUATION OF THE SIMULATED PT
    GlobalVector gv= tracker->idToDet(idet)->surface().toGlobal(isim->localDirection());
    float PP=isim->pabs();
    float ptsim=gv.perp()*PP;
    //PT^-1 RESOLUTION
    float ptresrel=((1./ptrec)-(1./ptsim));
    //   if (seed_plus){  

    //FILL PT RESOLUTION AND CHI2 HISTO 
    //IN PT BIN
      hptres->Fill(ptresrel);
      hchiSq->Fill(chiSquared);
      unsigned int iptsim= uint(ptsim/5 -1);
      if (iptsim>10) iptsim=10;
      
      if (iptsim==1) {
	hptres1->Fill(ptresrel);
	hchiSq1->Fill(chiSquared);
      }
      if (iptsim==2) {
	hchiSq2->Fill(chiSquared);
	hptres2->Fill(ptresrel);
      }
      if (iptsim==3) {
	hchiSq3->Fill(chiSquared);
	hptres3->Fill(ptresrel);
      }
      if (iptsim==4) {
	hchiSq4->Fill(chiSquared);
	hptres4->Fill(ptresrel);
      }
      if (iptsim==5) {
	hchiSq5->Fill(chiSquared);
	hptres5->Fill(ptresrel);
      }
      if (iptsim==6) {
	hchiSq6->Fill(chiSquared);
	hptres6->Fill(ptresrel);
      }
      if (iptsim==7) {
	hptres7->Fill(ptresrel);
	hchiSq7->Fill(chiSquared);
      }
      if (iptsim==8) {
	hchiSq8->Fill(chiSquared);
	hptres8->Fill(ptresrel);
      }
      if (iptsim==9) {
	hchiSq9->Fill(chiSquared);
	hptres9->Fill(ptresrel);
      }
      if (iptsim==10) {
	hchiSq10->Fill(chiSquared);
	hptres10->Fill(ptresrel);
      }
      // }
  
      makeResiduals((*(*seedcoll).begin()),
		    *trackrechitCollection,
		    e,
		    es);
  }
}
void ReadCosmicTracks::makeResiduals(const TrajectorySeed& seed,
				     const TrackingRecHitCollection &hits,
				     const edm::Event& e, 
				     const edm::EventSetup& es){



  //services
  es.get<IdealMagneticFieldRecord>().get(magfield);
  es.get<TrackerDigiGeometryRecord>().get(tracker);
  
 
  
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
  Trajectory traj( seed, seed.direction());
  traj.push(TrajectoryMeasurement(invalidState,State,recHit));


 
  TransientTrackingRecHit::RecHitContainer trans_hits;
  for (unsigned int icosmhit=0;icosmhit<hits.size();icosmhit++){

    TransientTrackingRecHit::RecHitPointer tmphit=RHBuilder->build(&(hits[icosmhit]));
    trans_hits.push_back(&(*tmphit));
    if (icosmhit>0){
      TSOS prSt= thePropagator->propagate(traj.lastMeasurement().updatedState(),
					  tracker->idToDet(hits[icosmhit].geographicalId())->surface());
      TSOS UpdatedState= theUpdator->update( prSt, *tmphit);
      traj.push(TrajectoryMeasurement(prSt,UpdatedState,tmphit
				      ,theEstimator->estimate(prSt, *tmphit).second ));
    }
  }
  TSOS startingState=  TrajectoryStateWithArbitraryError()
    (thePropagatorOp->propagate(traj.lastMeasurement().updatedState(),
				tracker->idToDet((*trans_hits.begin())->geographicalId())->surface()));
  const Trajectory ifitted= *(theFitter->fit(seed,trans_hits,startingState).begin());
  const Trajectory smooth=*(theSmoother->trajectories(ifitted).begin());

   std::vector<TrajectoryMeasurement> TMeas=smooth.measurements();
   vector<TrajectoryMeasurement>::iterator itm;
   for (itm=TMeas.begin();itm!=TMeas.end();itm++){
     TrajectoryMeasurementResidual*  TMR=new TrajectoryMeasurementResidual(*itm);
//     if    (subid==  StripSubdetector::TIB) hresTIB->Fill(r[0]);
//     if    (subid==  StripSubdetector::TOB) hresTOB->Fill(r[0]);
//     if    (subid==  StripSubdetector::TID) hresTID->Fill(r[0]);
//     if    (subid==  StripSubdetector::TEC) hresTEC->Fill(r[0]);
//     //cout<<r<<" "<<r[0]<<endl;
//     //  cout<<"tm "<<(*itm).updatedState().globalPosition()<<" "<<(*itm).recHit()->globalPosition()<<endl;
     delete TMR;
   }


}

void ReadCosmicTracks::endJob(){
  //FILL EFFICIENCY vs PT
  //AND CHI2 vs R
  for  (uint ik=0;ik<10;ik++) {
    heffpt->Fill(7.5+(ik*5),float(inum[ik])/float(iden[ik]));
    hchiR->Fill(10+(ik*10),ichiR[ik]/float(iden3[ik]));
  }
  //FILL EFFICIENCY vs NHIT
  for  (uint ik=0;ik<9;ik++) {
    heffhit->Fill(7.5+(ik*5),float(inum2[ik])/float(iden2[ik]));
  }
  //FILL CHI2 vs PT 
  hchipt->Fill(7.5,hchiSq1->GetMean());
  hchipt->Fill(12.5,hchiSq2->GetMean());
  hchipt->Fill(17.5,hchiSq3->GetMean());
  hchipt->Fill(22.5,hchiSq4->GetMean());
  hchipt->Fill(27.5,hchiSq5->GetMean());
  hchipt->Fill(32.5,hchiSq6->GetMean());
  hchipt->Fill(37.5,hchiSq7->GetMean());
  hchipt->Fill(42.5,hchiSq8->GetMean());
  hchipt->Fill(47.5,hchiSq9->GetMean());
  hchipt->Fill(52.5,hchiSq10->GetMean());

  //WRITE ROOT FILE
  hFile->Write();
  hFile->Close();
}
TrajectoryStateOnSurface
ReadCosmicTracks::startingTSOS(const TrajectorySeed& seed)const
{
  PTrajectoryStateOnDet pState( seed.startingState());
  const GeomDet* gdet  = (&(*tracker))->idToDet(DetId(pState.detId()));
  TSOS  State= tsTransform.transientState( pState, &(gdet->surface()), 
					   &(*magfield));
  return State;

}
