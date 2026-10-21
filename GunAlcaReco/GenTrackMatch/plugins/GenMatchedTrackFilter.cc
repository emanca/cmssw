// GenMatchedTrackFilter
//
// Loose gen-matched track selection for particle-gun ALCARECO production:
// keeps only tracks loosely matched (deltaR + relative pT) to a stable,
// charged generator-level particle. Filter semantics (EDFilter):
//   * emits a reco::TrackCollection with the matched clones (preserving
//     TrackExtra references, so downstream refits keep hits),
//   * the event passes iff it contains at least one matched track,
// so pure-pileup events are dropped and files carry only matched content.
// No ValueMap products (maps removed by design decision of 2026-10-21).
//
// Motivation: drop PU tracks in gun+premix samples while keeping the signal
// track and its decay/interaction daughters, per the campaign decision
// (see MuKPiGun_UL16_Premix_e2e_v1/CAMPAIGN_NOTES.md).

#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <cmath>
#include <string>
#include <vector>

class GenMatchedTrackFilter : public edm::stream::EDFilter<> {
public:
  explicit GenMatchedTrackFilter(const edm::ParameterSet &iConfig);
  bool filter(edm::Event &iEvent, const edm::EventSetup &iSetup) override;

private:
  edm::EDGetTokenT<reco::TrackCollection> srcToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  const double dRMatchMax_;
  const double dPtRelMatchMax_;
  const double genMinPt_;
  const double genMaxAbsEta_;
  const bool keepUnmatchedAll_;
  const unsigned keepEveryNthUnmatched_;
};

GenMatchedTrackFilter::GenMatchedTrackFilter(const edm::ParameterSet &iConfig)
    : srcToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      genToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
      dRMatchMax_(iConfig.getParameter<double>("dRMatchMax")),
      dPtRelMatchMax_(iConfig.getParameter<double>("dPtRelMatchMax")),
      genMinPt_(iConfig.getParameter<double>("genMinPt")),
      genMaxAbsEta_(iConfig.getParameter<double>("genMaxAbsEta")),
      keepUnmatchedAll_(iConfig.getParameter<bool>("keepUnmatchedAll")),
      keepEveryNthUnmatched_(iConfig.getParameter<unsigned>("keepEveryNthUnmatched")) {
  produces<reco::TrackCollection>();
}

bool GenMatchedTrackFilter::filter(edm::Event &iEvent, const edm::EventSetup &) {
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(srcToken_, tracks);
  edm::Handle<reco::GenParticleCollection> gens;
  iEvent.getByToken(genToken_, gens);

  std::vector<const reco::GenParticle *> stable;
  if (gens.isValid()) {
    for (const auto &g : *gens) {
      if (g.status() != 1)
        continue;
      if (g.charge() == 0)
        continue;
      if (g.pt() < genMinPt_)
        continue;
      if (std::abs(g.eta()) > genMaxAbsEta_)
        continue;
      stable.push_back(&g);
    }
  }

  auto outTracks = std::make_unique<reco::TrackCollection>();
  bool eventMatched = false;

  unsigned trackIndex = 0;
  for (const auto &trk : *tracks) {
    float bestR = dRMatchMax_;
    bool matched = false;
    for (const auto *g : stable) {
      const double dr = reco::deltaR(trk.eta(), trk.phi(), g->eta(), g->phi());
      if (dr >= bestR)
        continue;
      const float dptrel = std::abs(trk.pt() - g->pt()) / std::max(g->pt(), 1e-9);
      if (dptrel > dPtRelMatchMax_)
        continue;
      bestR = static_cast<float>(dr);
      matched = true;
    }
    const bool keep = matched || keepUnmatchedAll_ ||
                      (keepEveryNthUnmatched_ > 0 && (trackIndex % keepEveryNthUnmatched_) == 0);
    if (keep) {
      outTracks->push_back(trk);
      eventMatched = true;
    }
    ++trackIndex;
  }

  iEvent.put(std::move(outTracks));
  return eventMatched;
}

DEFINE_FWK_MODULE(GenMatchedTrackFilter);
