/// \class AlCaRecoSimTruthSelector
///
/// Selects the sim truth needed for the gen and sim-hit association of the
/// track refit (ResidualGlobalCorrectionMaker*) in full-event MC ALCARECO:
/// the SimTracks that come from a generator particle (genpartIndex >= 0) of
/// one of the configured |pdgId| (all if empty), and the PSimHits of those
/// SimTracks from the configured hit collections (same instance names).
/// SimVertices are not written: the vertIndex of the selected SimTracks
/// refers to the original SimVertex collection.
///
/// In a full Pythia event the tracker PSimHits and SimTracks are O(200 kB)
/// and O(70 kB) per event; for the muons of a J/psi this is O(1 kB).

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include <cstdlib>
#include <memory>
#include <unordered_set>
#include <vector>

class AlCaRecoSimTruthSelector : public edm::global::EDProducer<> {
public:
  explicit AlCaRecoSimTruthSelector(const edm::ParameterSet& cfg)
      : simTracksToken_(consumes<edm::SimTrackContainer>(cfg.getParameter<edm::InputTag>("simTracks"))) {
    for (const auto& pdgId : cfg.getParameter<std::vector<int>>("pdgIds")) {
      pdgIds_.insert(std::abs(pdgId));
    }
    for (const auto& tag : cfg.getParameter<std::vector<edm::InputTag>>("simHits")) {
      simHitsTokens_.push_back(consumes<edm::PSimHitContainer>(tag));
      simHitsPutTokens_.push_back(produces<edm::PSimHitContainer>(tag.instance()));
    }
    simTracksPutToken_ = produces<edm::SimTrackContainer>();
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("simTracks", edm::InputTag("g4SimHits", "", "SIM"));
    desc.add<std::vector<edm::InputTag>>("simHits", {});
    desc.add<std::vector<int>>("pdgIds", {13});
    descriptions.add("alcaRecoSimTruthSelector", desc);
  }

  void produce(edm::StreamID, edm::Event& evt, const edm::EventSetup&) const override {
    const auto& simTracks = evt.get(simTracksToken_);

    auto outTracks = std::make_unique<edm::SimTrackContainer>();
    std::unordered_set<unsigned int> trackIds;
    for (const auto& st : simTracks) {
      if (st.noGenpart())
        continue;
      if (!pdgIds_.empty() && !pdgIds_.count(std::abs(st.type())))
        continue;
      outTracks->push_back(st);
      trackIds.insert(st.trackId());
    }

    for (unsigned int i = 0; i < simHitsTokens_.size(); ++i) {
      auto outHits = std::make_unique<edm::PSimHitContainer>();
      for (const auto& hit : evt.get(simHitsTokens_[i])) {
        if (trackIds.count(hit.trackId()))
          outHits->push_back(hit);
      }
      evt.put(simHitsPutTokens_[i], std::move(outHits));
    }
    evt.put(simTracksPutToken_, std::move(outTracks));
  }

private:
  const edm::EDGetTokenT<edm::SimTrackContainer> simTracksToken_;
  std::vector<edm::EDGetTokenT<edm::PSimHitContainer>> simHitsTokens_;
  std::vector<edm::EDPutTokenT<edm::PSimHitContainer>> simHitsPutTokens_;
  edm::EDPutTokenT<edm::SimTrackContainer> simTracksPutToken_;
  std::unordered_set<int> pdgIds_;
};

DEFINE_FWK_MODULE(AlCaRecoSimTruthSelector);
