#include "PluginManager/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoTracker/SingleTrackPattern/test/ReadCosmicTracks.h"
#include "RecoTracker/SingleTrackPattern/test/AnalyzeMTCCTracks.h"
#include "RecoTracker/SingleTrackPattern/test/CombTrack.h"
DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(ReadCosmicTracks);
DEFINE_ANOTHER_FWK_MODULE(AnalyzeMTCCTracks);
DEFINE_ANOTHER_FWK_MODULE(CombTrack);
