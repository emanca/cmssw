#ifndef GunAlcaReco_GenTrackMatch_GenMatchedTrackSelector_h
#define GunAlcaReco_GenTrackMatch_GenMatchedTrackSelector_h

// Selects tracks loosely matched (deltaR + relative pT) to a stable, charged
// generator-level particle. Packaged as a CommonTools/UtilAlgos selector
// functor so that plugging it into ObjectSelectorStream gives the full stock
// ALCARECO clone family (tracks + TrackExtras + owned RecHits + SiPixel and
// SiStripClusters, references recursively re-pointed), exactly like every
// TkAl stream. Replaces the plain EDFilter version of 2026-10-21, which kept
// matched-track clones but left references dangling.

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <cmath>
#include <vector>

class GenMatchedTrackSelector {
public:
  typedef std::vector<const reco::Track *> container;
  typedef container input;
  GenMatchedTrackSelector(const edm::ParameterSet &cfg, edm::ConsumesCollector &iC)
      : genToken_(iC.consumes<reco::GenParticleCollection>(cfg.getParameter<edm::InputTag>("genParticles"))),
        dRMatchMax_(cfg.getParameter<double>("dRMatchMax")),
        dPtRelMatchMax_(cfg.getParameter<double>("dPtRelMatchMax")),
        genMinPt_(cfg.getParameter<double>("genMinPt")),
        genMaxAbsEta_(cfg.getParameter<double>("genMaxAbsEta")) {}

  container select(const container &tracks, const edm::Event &iEvent, const edm::EventSetup &) const {
    edm::Handle<reco::GenParticleCollection> gens;
    iEvent.getByToken(genToken_, gens);

    std::vector<const reco::GenParticle *> stable;
    if (gens.isValid()) {
      for (const auto &g : *gens) {
        if (g.status() != 1 || g.charge() == 0)
          continue;
        if (g.pt() < genMinPt_ || std::abs(g.eta()) > genMaxAbsEta_)
          continue;
        stable.push_back(&g);
      }
    }

    container out;
    for (const auto *trk : tracks) {
      float bestR = dRMatchMax_;
      bool matched = false;
      for (const auto *g : stable) {
        const double dr = reco::deltaR(trk->eta(), trk->phi(), g->eta(), g->phi());
        if (dr >= bestR)
          continue;
        const float dptrel = std::abs(trk->pt() - g->pt()) / std::max(g->pt(), 1e-9);
        if (dptrel > dPtRelMatchMax_)
          continue;
        bestR = static_cast<float>(dr);
        matched = true;
      }
      if (matched)
        out.push_back(trk);
    }
    return out;
  }

private:
  const edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  const double dRMatchMax_;
  const double dPtRelMatchMax_;
  const double genMinPt_;
  const double genMaxAbsEta_;
};

#endif
