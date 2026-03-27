#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

namespace {
struct FitQuality {
  bool isValid = false;
  double chi2 = -1.;
  double dof = -1.;
  double pval = -1.;
  bool isGood = false;
};

struct CandidateSummary {
  float dstPt;
  float d0Mass;
  float dstMass;
  float deltaMass;
  float d0Pval;
  float dstPval;
  float kPt;
  float kEta;
  float kPhi;
  int kCharge;
  float piPt;
  float piEta;
  float piPhi;
  int piCharge;
  float pisPt;
  float pisEta;
  float pisPhi;
  int pisCharge;
};
}

class DstToD0PiValidationAnalyzer : public edm::EDAnalyzer {
public:
  explicit DstToD0PiValidationAnalyzer(const edm::ParameterSet& iConfig);
  ~DstToD0PiValidationAnalyzer() override = default;

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  bool checkCharge(const reco::Track* kaon, const reco::Track* pion, const reco::Track* softPion) const;
  TLorentzVector makeP4(const reco::Track& track, double mass) const;
  void resetBranches();
  FitQuality fitQuality(RefCountedKinematicTree tree) const;
  RefCountedKinematicTree fitD0(const edm::EventSetup&, const reco::Track&, const reco::Track&) const;
  RefCountedKinematicTree fitDst(const edm::EventSetup&, const reco::Track&, const RefCountedKinematicParticle&) const;

  edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vertexToken_;

  double kaonMass_;
  double pionMass_;
  double softPionMass_;
  double kaonMassErr_;
  double pionMassErr_;
  double minD0Mass_;
  double maxD0Mass_;
  double minDstMass_;
  double maxDstMass_;
  double minDeltaMass_;
  double maxDeltaMass_;
  double pvalMin_;
  bool applyChargeFilter_;
  bool useUnsignedCharge_;
  int targetCharge_;

  TTree* tree_;
  int run_;
  int lumi_;
  long long event_;
  int nTracks_;
  int nCandidates_;
  float pvX_;
  float pvY_;
  float pvZ_;
  std::vector<float> massD0_;
  std::vector<float> massDst_;
  std::vector<float> deltaM_;
  std::vector<float> d0Pval_;
  std::vector<float> dstPval_;
  std::vector<float> kPt_;
  std::vector<float> kEta_;
  std::vector<float> kPhi_;
  std::vector<int> kCharge_;
  std::vector<float> piPt_;
  std::vector<float> piEta_;
  std::vector<float> piPhi_;
  std::vector<int> piCharge_;
  std::vector<float> pisPt_;
  std::vector<float> pisEta_;
  std::vector<float> pisPhi_;
  std::vector<int> pisCharge_;
};

DstToD0PiValidationAnalyzer::DstToD0PiValidationAnalyzer(const edm::ParameterSet& iConfig)
    : trackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      vertexToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
      kaonMass_(iConfig.getParameter<double>("kaonMass")),
      pionMass_(iConfig.getParameter<double>("pionMass")),
      softPionMass_(iConfig.getParameter<double>("softPionMass")),
      kaonMassErr_(1.e-6),
      pionMassErr_(1.e-6),
      minD0Mass_(iConfig.getParameter<double>("minD0Mass")),
      maxD0Mass_(iConfig.getParameter<double>("maxD0Mass")),
      minDstMass_(iConfig.getParameter<double>("minDstMass")),
      maxDstMass_(iConfig.getParameter<double>("maxDstMass")),
      minDeltaMass_(iConfig.getParameter<double>("minDeltaMass")),
      maxDeltaMass_(iConfig.getParameter<double>("maxDeltaMass")),
      pvalMin_(iConfig.getParameter<double>("pvalMin")),
      applyChargeFilter_(iConfig.getParameter<bool>("applyChargeFilter")),
      useUnsignedCharge_(iConfig.getParameter<bool>("useUnsignedCharge")),
      targetCharge_(iConfig.getParameter<int>("charge")),
      tree_(nullptr) {
  edm::Service<TFileService> fs;
  TFile& outputFile = fs->file();
  outputFile.cd();
  tree_ = new TTree("DstToD0PiTree", "D* -> D0 pi validation tree");
  tree_->SetDirectory(&outputFile);
  tree_->Branch("run", &run_, "run/I");
  tree_->Branch("lumi", &lumi_, "lumi/I");
  tree_->Branch("event", &event_, "event/L");
  tree_->Branch("nTracks", &nTracks_, "nTracks/I");
  tree_->Branch("nCandidates", &nCandidates_, "nCandidates/I");
  tree_->Branch("pv_x", &pvX_, "pv_x/F");
  tree_->Branch("pv_y", &pvY_, "pv_y/F");
  tree_->Branch("pv_z", &pvZ_, "pv_z/F");
  tree_->Branch("mass_D0", &massD0_);
  tree_->Branch("mass_Dst", &massDst_);
  tree_->Branch("deltaM", &deltaM_);
  tree_->Branch("d0_pval", &d0Pval_);
  tree_->Branch("dst_pval", &dstPval_);
  tree_->Branch("K_pt", &kPt_);
  tree_->Branch("K_eta", &kEta_);
  tree_->Branch("K_phi", &kPhi_);
  tree_->Branch("K_charge", &kCharge_);
  tree_->Branch("pi_pt", &piPt_);
  tree_->Branch("pi_eta", &piEta_);
  tree_->Branch("pi_phi", &piPhi_);
  tree_->Branch("pi_charge", &piCharge_);
  tree_->Branch("pis_pt", &pisPt_);
  tree_->Branch("pis_eta", &pisEta_);
  tree_->Branch("pis_phi", &pisPhi_);
  tree_->Branch("pis_charge", &pisCharge_);
}

void DstToD0PiValidationAnalyzer::resetBranches() {
  nTracks_ = 0;
  nCandidates_ = 0;
  pvX_ = 0.f;
  pvY_ = 0.f;
  pvZ_ = 0.f;
  massD0_.clear();
  massDst_.clear();
  deltaM_.clear();
  d0Pval_.clear();
  dstPval_.clear();
  kPt_.clear();
  kEta_.clear();
  kPhi_.clear();
  kCharge_.clear();
  piPt_.clear();
  piEta_.clear();
  piPhi_.clear();
  piCharge_.clear();
  pisPt_.clear();
  pisEta_.clear();
  pisPhi_.clear();
  pisCharge_.clear();
}

TLorentzVector DstToD0PiValidationAnalyzer::makeP4(const reco::Track& track, double mass) const {
  TLorentzVector p4;
  p4.SetXYZT(track.px(), track.py(), track.pz(), std::sqrt(track.p() * track.p() + mass * mass));
  return p4;
}

bool DstToD0PiValidationAnalyzer::checkCharge(const reco::Track* kaon,
                                              const reco::Track* pion,
                                              const reco::Track* softPion) const {
  if (!applyChargeFilter_)
    return true;

  const bool validD0Charges = (kaon->charge() + pion->charge() == 0);
  const bool validSoftCharge = (pion->charge() == softPion->charge());
  if (!validD0Charges || !validSoftCharge)
    return false;

  int totalCharge = kaon->charge() + pion->charge() + softPion->charge();
  if (useUnsignedCharge_)
    totalCharge = std::abs(totalCharge);

  return totalCharge == targetCharge_;
}

FitQuality DstToD0PiValidationAnalyzer::fitQuality(RefCountedKinematicTree tree) const {
  FitQuality quality;
  if (!tree || !tree->isValid())
    return quality;

  tree->movePointerToTheTop();
  RefCountedKinematicVertex vertex = tree->currentDecayVertex();
  if (!vertex || !vertex->vertexIsValid())
    return quality;

  quality.isValid = true;
  quality.chi2 = vertex->chiSquared();
  quality.dof = vertex->degreesOfFreedom();
  quality.pval = TMath::Prob(quality.chi2, quality.dof);
  quality.isGood = quality.pval > pvalMin_;
  return quality;
}

RefCountedKinematicTree DstToD0PiValidationAnalyzer::fitD0(const edm::EventSetup& iSetup,
                                                           const reco::Track& kaon,
                                                           const reco::Track& pion) const {
  edm::ESHandle<TransientTrackBuilder> ttBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttBuilder);

  reco::TransientTrack kaonTT = ttBuilder->build(kaon);
  reco::TransientTrack pionTT = ttBuilder->build(pion);

  KinematicParticleFactoryFromTransientTrack particleFactory;
  std::vector<RefCountedKinematicParticle> parts;
  float chi = 0.f;
  float ndf = 0.f;
  float kaonMassErr = kaonMassErr_;
  float pionMassErr = pionMassErr_;
  parts.push_back(particleFactory.particle(kaonTT, kaonMass_, chi, ndf, kaonMassErr));
  parts.push_back(particleFactory.particle(pionTT, pionMass_, chi, ndf, pionMassErr));

  KinematicParticleVertexFitter fitter;
  return fitter.fit(parts);
}

RefCountedKinematicTree DstToD0PiValidationAnalyzer::fitDst(const edm::EventSetup& iSetup,
                                                            const reco::Track& softPion,
                                                            const RefCountedKinematicParticle& d0) const {
  edm::ESHandle<TransientTrackBuilder> ttBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttBuilder);

  reco::TransientTrack softPionTT = ttBuilder->build(softPion);
  reco::TransientTrack d0TT = d0->refittedTransientTrack();

  KinematicParticleFactoryFromTransientTrack particleFactory;
  std::vector<RefCountedKinematicParticle> parts;
  float chi = 0.f;
  float ndf = 0.f;
  float d0Mass = d0->currentState().mass();
  float d0MassErr = std::sqrt(std::abs(d0->currentState().kinematicParametersError().matrix()(6, 6)));
  float pionMassErr = pionMassErr_;
  parts.push_back(particleFactory.particle(softPionTT, softPionMass_, chi, ndf, pionMassErr));
  parts.push_back(particleFactory.particle(d0TT, d0Mass, chi, ndf, d0MassErr));

  KinematicParticleVertexFitter fitter;
  return fitter.fit(parts);
}

void DstToD0PiValidationAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  resetBranches();

  run_ = iEvent.id().run();
  lumi_ = iEvent.luminosityBlock();
  event_ = iEvent.id().event();

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(trackToken_, tracks);

  edm::Handle<std::vector<reco::Vertex>> vertices;
  iEvent.getByToken(vertexToken_, vertices);

  if (tracks.isValid())
    nTracks_ = tracks->size();

  if (vertices.isValid() && !vertices->empty()) {
    const reco::Vertex& pv = (*vertices)[0];
    pvX_ = pv.x();
    pvY_ = pv.y();
    pvZ_ = pv.z();
  }

  std::vector<CandidateSummary> candidates;

  if (tracks.isValid() && tracks->size() >= 3) {
    for (unsigned int iK = 0; iK < tracks->size(); ++iK) {
      const reco::Track* kaon = &(*tracks)[iK];
      const int kaonCharge = kaon->charge();

      for (unsigned int iPi = 0; iPi < tracks->size(); ++iPi) {
        if (iPi == iK)
          continue;

        const reco::Track* pion = &(*tracks)[iPi];
        if (pion->charge() != -kaonCharge)
          continue;

        RefCountedKinematicTree d0Tree = fitD0(iSetup, *kaon, *pion);
        FitQuality d0Quality = fitQuality(d0Tree);
        if (!d0Quality.isGood)
          continue;

        d0Tree->movePointerToTheTop();
        RefCountedKinematicParticle d0 = d0Tree->currentParticle();
        const double d0Mass = d0->currentState().mass();
        if (d0Mass < minD0Mass_ || d0Mass > maxD0Mass_)
          continue;

        for (unsigned int iPis = 0; iPis < tracks->size(); ++iPis) {
          if (iPis == iK || iPis == iPi)
            continue;

          const reco::Track* softPion = &(*tracks)[iPis];
          if (softPion->charge() != pion->charge())
            continue;
          if (!checkCharge(kaon, pion, softPion))
            continue;

          RefCountedKinematicTree dstTree = fitDst(iSetup, *softPion, d0);
          FitQuality dstQuality = fitQuality(dstTree);
          if (!dstQuality.isGood)
            continue;

          dstTree->movePointerToTheTop();
          RefCountedKinematicParticle dst = dstTree->currentParticle();
          const double dstMass = dst->currentState().mass();
          const double deltaMass = dstMass - d0Mass;

          if (dstMass < minDstMass_ || dstMass > maxDstMass_)
            continue;
          if (deltaMass < minDeltaMass_ || deltaMass > maxDeltaMass_)
            continue;

          candidates.push_back(CandidateSummary{static_cast<float>(dst->currentState().globalMomentum().perp()),
                                                static_cast<float>(d0Mass),
                                                static_cast<float>(dstMass),
                                                static_cast<float>(deltaMass),
                                                static_cast<float>(d0Quality.pval),
                                                static_cast<float>(dstQuality.pval),
                                                static_cast<float>(kaon->pt()),
                                                static_cast<float>(kaon->eta()),
                                                static_cast<float>(kaon->phi()),
                                                kaon->charge(),
                                                static_cast<float>(pion->pt()),
                                                static_cast<float>(pion->eta()),
                                                static_cast<float>(pion->phi()),
                                                pion->charge(),
                                                static_cast<float>(softPion->pt()),
                                                static_cast<float>(softPion->eta()),
                                                static_cast<float>(softPion->phi()),
                                                softPion->charge()});
        }
      }
    }
  }

  std::sort(candidates.begin(), candidates.end(), [](const CandidateSummary& a, const CandidateSummary& b) {
    return a.dstPt > b.dstPt;
  });

  if (!candidates.empty())
    candidates.resize(1);

  nCandidates_ = candidates.size();
  for (const auto& candidate : candidates) {
    massD0_.push_back(candidate.d0Mass);
    massDst_.push_back(candidate.dstMass);
    deltaM_.push_back(candidate.deltaMass);
    d0Pval_.push_back(candidate.d0Pval);
    dstPval_.push_back(candidate.dstPval);
    kPt_.push_back(candidate.kPt);
    kEta_.push_back(candidate.kEta);
    kPhi_.push_back(candidate.kPhi);
    kCharge_.push_back(candidate.kCharge);
    piPt_.push_back(candidate.piPt);
    piEta_.push_back(candidate.piEta);
    piPhi_.push_back(candidate.piPhi);
    piCharge_.push_back(candidate.piCharge);
    pisPt_.push_back(candidate.pisPt);
    pisEta_.push_back(candidate.pisEta);
    pisPhi_.push_back(candidate.pisPhi);
    pisCharge_.push_back(candidate.pisCharge);
  }

  tree_->Fill();
}

DEFINE_FWK_MODULE(DstToD0PiValidationAnalyzer);
