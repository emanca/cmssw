#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "TLorentzVector.h"
#include "TMath.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace {
struct FitQuality {
  bool isValid = false;
  bool isGood = false;
  double chi2 = -1.;
  double dof = -1.;
  double pval = -1.;
};

struct BestCandidate {
  bool found = false;
  unsigned int iK = 0;
  unsigned int iPi = 0;
  unsigned int iPis = 0;
  double pt = -1.;
};

constexpr unsigned int kD0FitExceptionCounter = 11;
constexpr unsigned int kDstFitExceptionCounter = 12;

TLorentzVector p4FromTrack(const reco::Track& track, const double mass) {
  TLorentzVector p4;
  p4.SetXYZM(track.px(), track.py(), track.pz(), mass);
  return p4;
}

TLorentzVector p4FromKinParticle(const RefCountedKinematicParticle& particle) {
  const auto momentum = particle->currentState().globalMomentum();
  TLorentzVector p4;
  p4.SetXYZM(momentum.x(), momentum.y(), momentum.z(), particle->currentState().mass());
  return p4;
}

double massWithHypotheses(const TLorentzVector& p1, const double m1, const TLorentzVector& p2, const double m2) {
  const double e1 = std::hypot(m1, p1.P());
  const double e2 = std::hypot(m2, p2.P());
  double p1p2 = p1.Pt() * p2.Pt();
  p1p2 *= std::cos(p1.Phi() - p2.Phi()) + std::sinh(p1.Eta()) * std::sinh(p2.Eta());
  const double mass2 = m1 * m1 + m2 * m2 + 2. * (e1 * e2 - p1p2);
  return mass2 > 0. ? std::sqrt(mass2) : 0.;
}

double safeRatio(const double numerator, const double denominator) {
  return denominator > 0. ? numerator / denominator : -99.;
}

double pointingCos(const reco::Vertex& pv,
                   const RefCountedKinematicVertex& vertex,
                   const RefCountedKinematicParticle& particle) {
  const TVector3 displacement(vertex->position().x() - pv.x(),
                              vertex->position().y() - pv.y(),
                              vertex->position().z() - pv.z());
  const auto momentum = particle->currentState().globalMomentum();
  const TVector3 p3(momentum.x(), momentum.y(), momentum.z());
  if (displacement.Mag() == 0. || p3.Mag() == 0.)
    return -99.;
  return std::cos(displacement.Angle(p3));
}

double pointingCosTransverse(const reco::Vertex& pv,
                             const RefCountedKinematicVertex& vertex,
                             const RefCountedKinematicParticle& particle) {
  const TVector3 displacement(vertex->position().x() - pv.x(), vertex->position().y() - pv.y(), 0.);
  const auto momentum = particle->currentState().globalMomentum();
  const TVector3 p3(momentum.x(), momentum.y(), 0.);
  if (displacement.Mag() == 0. || p3.Mag() == 0.)
    return -99.;
  return std::cos(displacement.Angle(p3));
}

std::pair<double, double> vertexDistance(const reco::Vertex& pv, const RefCountedKinematicVertex& vertex) {
  const double dx = vertex->position().x() - pv.x();
  const double dy = vertex->position().y() - pv.y();
  const double dz = vertex->position().z() - pv.z();
  const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
  if (distance <= 0.)
    return std::make_pair(distance, 0.);

  const double direction[3] = {dx / distance, dy / distance, dz / distance};
  const auto pvCov = pv.covariance();
  const auto vertexCov = vertex->error().matrix();
  double error2 = 0.;
  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j) {
      error2 += direction[i] * (pvCov.At(i, j) + vertexCov.At(i, j)) * direction[j];
    }
  }
  return std::make_pair(distance, error2 > 0. ? std::sqrt(error2) : 0.);
}

std::pair<double, double> transverseDistanceFromBeamSpot(const reco::BeamSpot& beamSpot,
                                                         const RefCountedKinematicVertex& vertex) {
  const double dx = vertex->position().x() - beamSpot.x0();
  const double dy = vertex->position().y() - beamSpot.y0();
  const double distance = std::sqrt(dx * dx + dy * dy);
  if (distance <= 0.)
    return std::make_pair(distance, 0.);

  const double direction[2] = {dx / distance, dy / distance};
  const auto bsCov = beamSpot.covariance();
  const auto vertexCov = vertex->error().matrix();
  double error2 = 0.;
  for (unsigned int i = 0; i < 2; ++i) {
    for (unsigned int j = 0; j < 2; ++j) {
      error2 += direction[i] * (bsCov.At(i, j) + vertexCov.At(i, j)) * direction[j];
    }
  }
  return std::make_pair(distance, error2 > 0. ? std::sqrt(error2) : 0.);
}

std::pair<double, double> transverseDistance(const reco::Vertex& pv, const RefCountedKinematicVertex& vertex) {
  const double dx = vertex->position().x() - pv.x();
  const double dy = vertex->position().y() - pv.y();
  const double distance = std::sqrt(dx * dx + dy * dy);
  if (distance <= 0.)
    return std::make_pair(distance, 0.);

  const double direction[2] = {dx / distance, dy / distance};
  const auto pvCov = pv.covariance();
  const auto vertexCov = vertex->error().matrix();
  double error2 = 0.;
  for (unsigned int i = 0; i < 2; ++i) {
    for (unsigned int j = 0; j < 2; ++j) {
      error2 += direction[i] * (pvCov.At(i, j) + vertexCov.At(i, j)) * direction[j];
    }
  }
  return std::make_pair(distance, error2 > 0. ? std::sqrt(error2) : 0.);
}
}  // namespace

class DstToD0PiProducer : public edm::EDProducer {
public:
  explicit DstToD0PiProducer(const edm::ParameterSet&);
  ~DstToD0PiProducer() override;

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  RefCountedKinematicTree fitD0(const edm::EventSetup&, const reco::Track&, const reco::Track&) const;
  RefCountedKinematicTree fitDst(const edm::EventSetup&,
                                 const reco::Track&,
                                 const RefCountedKinematicParticle&) const;
  FitQuality fitQuality(RefCountedKinematicTree) const;
  bool passTrackPreselection(const reco::Track&) const;
  bool passSoftPionPreselection(const reco::Track&) const;
  void clearTreeVectors();
  void addP4(const TLorentzVector&,
             std::vector<float>&,
             std::vector<float>&,
             std::vector<float>&,
             std::vector<float>&,
             std::vector<float>&);
  bool fillBestCandidate(const edm::EventSetup&,
                         const reco::Track&,
                         const reco::Track&,
                         const reco::Track&,
                         const reco::Vertex&,
                         const reco::BeamSpot&);

  edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vertexToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken_;

  std::vector<std::string> hltPaths_;
  bool doTrigger_ = false;
  bool deferTreeFill_ = false;
  int verbose_ = 0;

  double massPi_ = 0.139570;
  double massK_ = 0.493677;
  double massD0_ = 1.86483;
  double massDst_ = 2.01026;
  double massErrPi_ = 1.e-6;
  double massErrK_ = 1.e-6;

  double hadPtMin_ = 0.6;
  // A negative value preserves the historical common threshold for the soft pion.
  double softPionPtMin_ = -1.0;
  double pvalChi2VtxMin_ = 0.005;
  double dmD0Max_ = 0.05;
  double dmDstMax_ = 0.15;
  double dmDstMinusD0Max_ = 0.003;

  TTree* tree_ = nullptr;
  unsigned int run_ = 0;
  unsigned int lumi_ = 0;
  unsigned long long event_ = 0;
  int n_K_ = 0;
  int n_pi_ = 0;
  int n_D0_ = 0;
  int n_pis_ = 0;
  int n_Dst_ = 0;
  int n_BestDst_ = 0;

  std::vector<float> K_dxy_BS_, K_dxyErr_BS_, K_sigdxy_BS_, K_dxy_PV_, K_dxyErr_PV_, K_sigdxy_PV_;
  std::vector<float> K_norm_chi2_, K_N_valid_hits_, K_lostInnerHits_, K_pt_, K_pz_, K_E_, K_eta_, K_phi_;
  std::vector<float> pi_dxy_BS_, pi_dxyErr_BS_, pi_sigdxy_BS_, pi_dxy_PV_, pi_dxyErr_PV_, pi_sigdxy_PV_;
  std::vector<float> pi_norm_chi2_, pi_N_valid_hits_, pi_lostInnerHits_, pi_pt_, pi_pz_, pi_E_, pi_eta_, pi_phi_;
  std::vector<float> pis_dxy_BS_, pis_dxyErr_BS_, pis_sigdxy_BS_, pis_dxy_PV_, pis_dxyErr_PV_, pis_sigdxy_PV_;
  std::vector<float> pis_norm_chi2_, pis_N_valid_hits_, pis_lostInnerHits_, pis_pt_, pis_pz_, pis_E_, pis_eta_, pis_phi_;
  std::vector<float> massb_piK_, massb_KK_, massb_Kpi_conj_;
  std::vector<float> chi2_piK_, dof_piK_, pval_piK_, mass_piK_;
  std::vector<float> D0_pt_, D0_pz_, D0_E_, D0_eta_, D0_phi_;
  std::vector<float> K_refitpiK_pt_, K_refitpiK_pz_, K_refitpiK_E_, K_refitpiK_eta_, K_refitpiK_phi_;
  std::vector<float> pi_refitpiK_pt_, pi_refitpiK_pz_, pi_refitpiK_E_, pi_refitpiK_eta_, pi_refitpiK_phi_;
  std::vector<float> mass_piK_hKK_, mass_piK_hpipi_, mass_piK_hKpi_;
  std::vector<float> cos_D0_PV_, cosT_D0_PV_, d_vtxD0_PV_, sigd_vtxD0_PV_, dxy_vtxD0_BS_, sigdxy_vtxD0_BS_;
  std::vector<float> chi2_D0pis_, dof_D0pis_, pval_D0pis_, mass_D0pis_;
  std::vector<float> cos_Dst_PV_, cosT_Dst_PV_, d_vtxDst_PV_, sigd_vtxDst_PV_, dxy_vtxDst_PV_, sigdxy_vtxDst_PV_;
  std::vector<float> pis_refitD0pis_pt_, pis_refitD0pis_pz_, pis_refitD0pis_E_, pis_refitD0pis_eta_, pis_refitD0pis_phi_;
  std::vector<float> D0_refitD0pis_pt_, D0_refitD0pis_pz_, D0_refitD0pis_E_, D0_refitD0pis_eta_, D0_refitD0pis_phi_;
  std::vector<float> deltaMass_D0pis_;
  std::vector<int> K_index_, pi_index_, pis_index_;

  // Fit helpers are const because fitting does not alter candidate state; the
  // exception counters are diagnostic state updated only when a fit throws.
  mutable std::vector<unsigned long long> counters_;
};

DstToD0PiProducer::DstToD0PiProducer(const edm::ParameterSet& iConfig)
    : trackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      vertexToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
      beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))) {
  verbose_ = iConfig.existsAs<int>("verbose") ? iConfig.getParameter<int>("verbose") : 0;
  deferTreeFill_ = iConfig.existsAs<bool>("deferTreeFill") ? iConfig.getParameter<bool>("deferTreeFill") : false;
  doTrigger_ = iConfig.existsAs<bool>("doTrigger") ? iConfig.getParameter<bool>("doTrigger") : false;
  if (doTrigger_) {
    triggerBitsToken_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"));
    hltPaths_ = iConfig.getParameter<std::vector<std::string>>("hltPaths");
  }

  // Keep the transiently selected best triplet available to downstream refits.
  // The existing diagnostic tree and its candidate selection are unchanged.
  produces<reco::TrackCollection>("dstarTracks");

  if (iConfig.existsAs<double>("hadPtMin"))
    hadPtMin_ = iConfig.getParameter<double>("hadPtMin");
  if (iConfig.existsAs<double>("softPionPtMin"))
    softPionPtMin_ = iConfig.getParameter<double>("softPionPtMin");
  if (iConfig.existsAs<double>("pvalChi2VtxMin"))
    pvalChi2VtxMin_ = iConfig.getParameter<double>("pvalChi2VtxMin");
  if (iConfig.existsAs<double>("dmD0Max"))
    dmD0Max_ = iConfig.getParameter<double>("dmD0Max");
  if (iConfig.existsAs<double>("dmDstMax"))
    dmDstMax_ = iConfig.getParameter<double>("dmDstMax");
  if (iConfig.existsAs<double>("dmDstMinusD0Max"))
    dmDstMinusD0Max_ = iConfig.getParameter<double>("dmDstMinusD0Max");
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("Events", "BPH-like D* -> D0 pi validation tree");
  tree_->Branch("Run", &run_);
  tree_->Branch("LumiBlock", &lumi_);
  tree_->Branch("eventNumber", &event_);
  tree_->Branch("n_K", &n_K_);
  tree_->Branch("n_pi", &n_pi_);
  tree_->Branch("n_D0", &n_D0_);
  tree_->Branch("n_pis", &n_pis_);
  tree_->Branch("n_Dst", &n_Dst_);
  tree_->Branch("n_BestDst", &n_BestDst_);

  tree_->Branch("K_index", &K_index_);
  tree_->Branch("pi_index", &pi_index_);
  tree_->Branch("pis_index", &pis_index_);
  tree_->Branch("K_dxy_BS", &K_dxy_BS_);
  tree_->Branch("K_dxyErr_BS", &K_dxyErr_BS_);
  tree_->Branch("K_sigdxy_BS", &K_sigdxy_BS_);
  tree_->Branch("K_dxy_PV", &K_dxy_PV_);
  tree_->Branch("K_dxyErr_PV", &K_dxyErr_PV_);
  tree_->Branch("K_sigdxy_PV", &K_sigdxy_PV_);
  tree_->Branch("K_norm_chi2", &K_norm_chi2_);
  tree_->Branch("K_N_valid_hits", &K_N_valid_hits_);
  tree_->Branch("K_lostInnerHits", &K_lostInnerHits_);
  tree_->Branch("K_pt", &K_pt_);
  tree_->Branch("K_pz", &K_pz_);
  tree_->Branch("K_E", &K_E_);
  tree_->Branch("K_eta", &K_eta_);
  tree_->Branch("K_phi", &K_phi_);
  tree_->Branch("pi_dxy_BS", &pi_dxy_BS_);
  tree_->Branch("pi_dxyErr_BS", &pi_dxyErr_BS_);
  tree_->Branch("pi_sigdxy_BS", &pi_sigdxy_BS_);
  tree_->Branch("pi_dxy_PV", &pi_dxy_PV_);
  tree_->Branch("pi_dxyErr_PV", &pi_dxyErr_PV_);
  tree_->Branch("pi_sigdxy_PV", &pi_sigdxy_PV_);
  tree_->Branch("pi_norm_chi2", &pi_norm_chi2_);
  tree_->Branch("pi_N_valid_hits", &pi_N_valid_hits_);
  tree_->Branch("pi_lostInnerHits", &pi_lostInnerHits_);
  tree_->Branch("pi_pt", &pi_pt_);
  tree_->Branch("pi_pz", &pi_pz_);
  tree_->Branch("pi_E", &pi_E_);
  tree_->Branch("pi_eta", &pi_eta_);
  tree_->Branch("pi_phi", &pi_phi_);
  tree_->Branch("massb_piK", &massb_piK_);
  tree_->Branch("massb_KK", &massb_KK_);
  tree_->Branch("massb_Kpi_conj", &massb_Kpi_conj_);
  tree_->Branch("chi2_piK", &chi2_piK_);
  tree_->Branch("dof_piK", &dof_piK_);
  tree_->Branch("pval_piK", &pval_piK_);
  tree_->Branch("mass_piK", &mass_piK_);
  tree_->Branch("D0_pt", &D0_pt_);
  tree_->Branch("D0_pz", &D0_pz_);
  tree_->Branch("D0_E", &D0_E_);
  tree_->Branch("D0_eta", &D0_eta_);
  tree_->Branch("D0_phi", &D0_phi_);
  tree_->Branch("K_refitpiK_pt", &K_refitpiK_pt_);
  tree_->Branch("K_refitpiK_pz", &K_refitpiK_pz_);
  tree_->Branch("K_refitpiK_E", &K_refitpiK_E_);
  tree_->Branch("K_refitpiK_eta", &K_refitpiK_eta_);
  tree_->Branch("K_refitpiK_phi", &K_refitpiK_phi_);
  tree_->Branch("pi_refitpiK_pt", &pi_refitpiK_pt_);
  tree_->Branch("pi_refitpiK_pz", &pi_refitpiK_pz_);
  tree_->Branch("pi_refitpiK_E", &pi_refitpiK_E_);
  tree_->Branch("pi_refitpiK_eta", &pi_refitpiK_eta_);
  tree_->Branch("pi_refitpiK_phi", &pi_refitpiK_phi_);
  tree_->Branch("mass_piK_hKK", &mass_piK_hKK_);
  tree_->Branch("mass_piK_hpipi", &mass_piK_hpipi_);
  tree_->Branch("mass_piK_hKpi", &mass_piK_hKpi_);
  tree_->Branch("cos_D0_PV", &cos_D0_PV_);
  tree_->Branch("cosT_D0_PV", &cosT_D0_PV_);
  tree_->Branch("d_vtxD0_PV", &d_vtxD0_PV_);
  tree_->Branch("sigd_vtxD0_PV", &sigd_vtxD0_PV_);
  tree_->Branch("dxy_vtxD0_BS", &dxy_vtxD0_BS_);
  tree_->Branch("sigdxy_vtxD0_BS", &sigdxy_vtxD0_BS_);
  tree_->Branch("pis_dxy_BS", &pis_dxy_BS_);
  tree_->Branch("pis_dxyErr_BS", &pis_dxyErr_BS_);
  tree_->Branch("pis_sigdxy_BS", &pis_sigdxy_BS_);
  tree_->Branch("pis_dxy_PV", &pis_dxy_PV_);
  tree_->Branch("pis_dxyErr_PV", &pis_dxyErr_PV_);
  tree_->Branch("pis_sigdxy_PV", &pis_sigdxy_PV_);
  tree_->Branch("pis_norm_chi2", &pis_norm_chi2_);
  tree_->Branch("pis_N_valid_hits", &pis_N_valid_hits_);
  tree_->Branch("pis_lostInnerHits", &pis_lostInnerHits_);
  tree_->Branch("pis_pt", &pis_pt_);
  tree_->Branch("pis_pz", &pis_pz_);
  tree_->Branch("pis_E", &pis_E_);
  tree_->Branch("pis_eta", &pis_eta_);
  tree_->Branch("pis_phi", &pis_phi_);
  tree_->Branch("chi2_D0pis", &chi2_D0pis_);
  tree_->Branch("dof_D0pis", &dof_D0pis_);
  tree_->Branch("pval_D0pis", &pval_D0pis_);
  tree_->Branch("mass_D0pis", &mass_D0pis_);
  tree_->Branch("deltaMass_D0pis", &deltaMass_D0pis_);
  tree_->Branch("cos_Dst_PV", &cos_Dst_PV_);
  tree_->Branch("cosT_Dst_PV", &cosT_Dst_PV_);
  tree_->Branch("d_vtxDst_PV", &d_vtxDst_PV_);
  tree_->Branch("sigd_vtxDst_PV", &sigd_vtxDst_PV_);
  tree_->Branch("dxy_vtxDst_PV", &dxy_vtxDst_PV_);
  tree_->Branch("sigdxy_vtxDst_PV", &sigdxy_vtxDst_PV_);
  tree_->Branch("pis_refitD0pis_pt", &pis_refitD0pis_pt_);
  tree_->Branch("pis_refitD0pis_pz", &pis_refitD0pis_pz_);
  tree_->Branch("pis_refitD0pis_E", &pis_refitD0pis_E_);
  tree_->Branch("pis_refitD0pis_eta", &pis_refitD0pis_eta_);
  tree_->Branch("pis_refitD0pis_phi", &pis_refitD0pis_phi_);
  tree_->Branch("D0_refitD0pis_pt", &D0_refitD0pis_pt_);
  tree_->Branch("D0_refitD0pis_pz", &D0_refitD0pis_pz_);
  tree_->Branch("D0_refitD0pis_E", &D0_refitD0pis_E_);
  tree_->Branch("D0_refitD0pis_eta", &D0_refitD0pis_eta_);
  tree_->Branch("D0_refitD0pis_phi", &D0_refitD0pis_phi_);

  counters_.assign(13, 0);
}

DstToD0PiProducer::~DstToD0PiProducer() {
  std::cout << "DstToD0PiProducer counters:";
  for (const auto counter : counters_)
    std::cout << " " << counter;
  std::cout << std::endl;
}

bool DstToD0PiProducer::passTrackPreselection(const reco::Track& track) const {
  return track.pt() >= hadPtMin_;
}

bool DstToD0PiProducer::passSoftPionPreselection(const reco::Track& track) const {
  const double threshold = softPionPtMin_ >= 0. ? softPionPtMin_ : hadPtMin_;
  return track.pt() >= threshold;
}

RefCountedKinematicTree DstToD0PiProducer::fitD0(const edm::EventSetup& iSetup,
                                                const reco::Track& pi,
                                                const reco::Track& K) const {
  edm::ESHandle<TransientTrackBuilder> ttBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttBuilder);

  reco::TransientTrack piTT = ttBuilder->build(pi);
  reco::TransientTrack KTT = ttBuilder->build(K);

  KinematicParticleFactoryFromTransientTrack particleFactory;
  std::vector<RefCountedKinematicParticle> parts;
  float chi = 0.f;
  float ndf = 0.f;
  float dmK = massErrK_;
  float dmPi = massErrPi_;
  parts.push_back(particleFactory.particle(KTT, massK_, chi, ndf, dmK));
  parts.push_back(particleFactory.particle(piTT, massPi_, chi, ndf, dmPi));

  KinematicParticleVertexFitter fitter;
  try {
    return fitter.fit(parts);
  } catch (const std::exception& exception) {
    ++counters_[kD0FitExceptionCounter];
    edm::LogWarning("DstToD0PiProducer")
        << "D0 vertex fit threw std::exception in event " << run_ << ":" << lumi_ << ":" << event_
        << " (K pt=" << K.pt() << ", pi pt=" << pi.pt() << "): " << exception.what()
        << ". Candidate rejected.";
    return RefCountedKinematicTree();
  }
}

RefCountedKinematicTree DstToD0PiProducer::fitDst(const edm::EventSetup& iSetup,
                                                 const reco::Track& softPi,
                                                 const RefCountedKinematicParticle& D0) const {
  edm::ESHandle<TransientTrackBuilder> ttBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttBuilder);

  reco::TransientTrack softPiTT = ttBuilder->build(softPi);
  reco::TransientTrack D0TT = D0->refittedTransientTrack();

  KinematicParticleFactoryFromTransientTrack particleFactory;
  std::vector<RefCountedKinematicParticle> parts;
  float chi = 0.f;
  float ndf = 0.f;
  float dmPi = massErrPi_;
  float mD0 = D0->currentState().mass();
  float dmD0 = std::sqrt(std::abs(D0->currentState().kinematicParametersError().matrix()(6, 6)));
  parts.push_back(particleFactory.particle(softPiTT, massPi_, chi, ndf, dmPi));
  parts.push_back(particleFactory.particle(D0TT, mD0, chi, ndf, dmD0));

  KinematicParticleVertexFitter fitter;
  try {
    return fitter.fit(parts);
  } catch (const std::exception& exception) {
    ++counters_[kDstFitExceptionCounter];
    edm::LogWarning("DstToD0PiProducer")
        << "D* vertex fit threw std::exception in event " << run_ << ":" << lumi_ << ":" << event_
        << " (soft-pion pt=" << softPi.pt() << "): " << exception.what() << ". Candidate rejected.";
    return RefCountedKinematicTree();
  }
}

FitQuality DstToD0PiProducer::fitQuality(RefCountedKinematicTree tree) const {
  FitQuality quality;
  if (!tree || !tree->isValid())
    return quality;

  tree->movePointerToTheTop();
  const RefCountedKinematicVertex vertex = tree->currentDecayVertex();
  if (!vertex || !vertex->vertexIsValid())
    return quality;

  quality.isValid = true;
  quality.chi2 = vertex->chiSquared();
  quality.dof = vertex->degreesOfFreedom();
  if (quality.chi2 < 0.) {
    quality.pval = -1.;
    quality.isGood = false;
    return quality;
  }
  quality.pval = TMath::Prob(quality.chi2, quality.dof);
  quality.isGood = pvalChi2VtxMin_ <= 0. ? true : quality.pval > pvalChi2VtxMin_;
  return quality;
}

void DstToD0PiProducer::clearTreeVectors() {
  n_K_ = 0;
  n_pi_ = 0;
  n_D0_ = 0;
  n_pis_ = 0;
  n_Dst_ = 0;
  n_BestDst_ = 0;
  K_index_.clear();
  pi_index_.clear();
  pis_index_.clear();
  K_dxy_BS_.clear();
  K_dxyErr_BS_.clear();
  K_sigdxy_BS_.clear();
  K_dxy_PV_.clear();
  K_dxyErr_PV_.clear();
  K_sigdxy_PV_.clear();
  K_norm_chi2_.clear();
  K_N_valid_hits_.clear();
  K_lostInnerHits_.clear();
  K_pt_.clear();
  K_pz_.clear();
  K_E_.clear();
  K_eta_.clear();
  K_phi_.clear();
  pi_dxy_BS_.clear();
  pi_dxyErr_BS_.clear();
  pi_sigdxy_BS_.clear();
  pi_dxy_PV_.clear();
  pi_dxyErr_PV_.clear();
  pi_sigdxy_PV_.clear();
  pi_norm_chi2_.clear();
  pi_N_valid_hits_.clear();
  pi_lostInnerHits_.clear();
  pi_pt_.clear();
  pi_pz_.clear();
  pi_E_.clear();
  pi_eta_.clear();
  pi_phi_.clear();
  pis_dxy_BS_.clear();
  pis_dxyErr_BS_.clear();
  pis_sigdxy_BS_.clear();
  pis_dxy_PV_.clear();
  pis_dxyErr_PV_.clear();
  pis_sigdxy_PV_.clear();
  pis_norm_chi2_.clear();
  pis_N_valid_hits_.clear();
  pis_lostInnerHits_.clear();
  pis_pt_.clear();
  pis_pz_.clear();
  pis_E_.clear();
  pis_eta_.clear();
  pis_phi_.clear();
  massb_piK_.clear();
  massb_KK_.clear();
  massb_Kpi_conj_.clear();
  chi2_piK_.clear();
  dof_piK_.clear();
  pval_piK_.clear();
  mass_piK_.clear();
  D0_pt_.clear();
  D0_pz_.clear();
  D0_E_.clear();
  D0_eta_.clear();
  D0_phi_.clear();
  K_refitpiK_pt_.clear();
  K_refitpiK_pz_.clear();
  K_refitpiK_E_.clear();
  K_refitpiK_eta_.clear();
  K_refitpiK_phi_.clear();
  pi_refitpiK_pt_.clear();
  pi_refitpiK_pz_.clear();
  pi_refitpiK_E_.clear();
  pi_refitpiK_eta_.clear();
  pi_refitpiK_phi_.clear();
  mass_piK_hKK_.clear();
  mass_piK_hpipi_.clear();
  mass_piK_hKpi_.clear();
  cos_D0_PV_.clear();
  cosT_D0_PV_.clear();
  d_vtxD0_PV_.clear();
  sigd_vtxD0_PV_.clear();
  dxy_vtxD0_BS_.clear();
  sigdxy_vtxD0_BS_.clear();
  chi2_D0pis_.clear();
  dof_D0pis_.clear();
  pval_D0pis_.clear();
  mass_D0pis_.clear();
  deltaMass_D0pis_.clear();
  cos_Dst_PV_.clear();
  cosT_Dst_PV_.clear();
  d_vtxDst_PV_.clear();
  sigd_vtxDst_PV_.clear();
  dxy_vtxDst_PV_.clear();
  sigdxy_vtxDst_PV_.clear();
  pis_refitD0pis_pt_.clear();
  pis_refitD0pis_pz_.clear();
  pis_refitD0pis_E_.clear();
  pis_refitD0pis_eta_.clear();
  pis_refitD0pis_phi_.clear();
  D0_refitD0pis_pt_.clear();
  D0_refitD0pis_pz_.clear();
  D0_refitD0pis_E_.clear();
  D0_refitD0pis_eta_.clear();
  D0_refitD0pis_phi_.clear();
}

void DstToD0PiProducer::addP4(const TLorentzVector& p4,
                             std::vector<float>& pt,
                             std::vector<float>& pz,
                             std::vector<float>& energy,
                             std::vector<float>& eta,
                             std::vector<float>& phi) {
  pt.push_back(p4.Pt());
  pz.push_back(p4.Pz());
  energy.push_back(p4.E());
  eta.push_back(p4.Eta());
  phi.push_back(p4.Phi());
}

bool DstToD0PiProducer::fillBestCandidate(const edm::EventSetup& iSetup,
                                          const reco::Track& K,
                                          const reco::Track& pi,
                                          const reco::Track& pis,
                                          const reco::Vertex& primaryVtx,
                                          const reco::BeamSpot& beamSpot) {
  RefCountedKinematicTree D0KinTree = fitD0(iSetup, pi, K);
  const FitQuality resPiK = fitQuality(D0KinTree);
  if (!resPiK.isGood)
    return false;

  D0KinTree->movePointerToTheTop();
  const double massPiK = D0KinTree->currentParticle()->currentState().mass();
  const RefCountedKinematicParticle D0 = D0KinTree->currentParticle();
  const RefCountedKinematicVertex vtxD0 = D0KinTree->currentDecayVertex();

  RefCountedKinematicTree DstKinTree = fitDst(iSetup, pis, D0);
  const FitQuality resD0Pis = fitQuality(DstKinTree);
  if (!resD0Pis.isGood)
    return false;

  DstKinTree->movePointerToTheTop();
  const RefCountedKinematicParticle Dst = DstKinTree->currentParticle();
  const RefCountedKinematicVertex vtxDst = DstKinTree->currentDecayVertex();
  const double massD0Pis = Dst->currentState().mass();

  const auto KdxyBS = K.dxy(beamSpot.position());
  const auto KdzBS = K.dz(beamSpot.position());
  const auto KdxyErrBS = K.dxyError();
  const auto KdxyPV = K.dxy(primaryVtx.position());
  const auto KdxyErrPV = K.dxyError();
  const auto pidxyBS = pi.dxy(beamSpot.position());
  const auto pidzBS = pi.dz(beamSpot.position());
  const auto pidxyErrBS = pi.dxyError();
  const auto pidxyPV = pi.dxy(primaryVtx.position());
  const auto pidxyErrPV = pi.dxyError();
  const auto pisdxyBS = pis.dxy(beamSpot.position());
  const auto pisdxyErrBS = pis.dxyError();
  const auto pisdxyPV = pis.dxy(primaryVtx.position());
  const auto pisdxyErrPV = pis.dxyError();

  K_dxy_BS_.push_back(KdxyBS);
  K_dxyErr_BS_.push_back(KdxyErrBS);
  K_sigdxy_BS_.push_back(safeRatio(std::abs(KdxyBS), KdxyErrBS));
  K_dxy_PV_.push_back(KdxyPV);
  K_dxyErr_PV_.push_back(KdxyErrPV);
  K_sigdxy_PV_.push_back(safeRatio(std::abs(KdxyPV), KdxyErrPV));
  K_norm_chi2_.push_back(K.normalizedChi2());
  K_N_valid_hits_.push_back(K.numberOfValidHits());
  K_lostInnerHits_.push_back(K.numberOfLostHits());
  addP4(p4FromTrack(K, massK_), K_pt_, K_pz_, K_E_, K_eta_, K_phi_);

  pi_dxy_BS_.push_back(pidxyBS);
  pi_dxyErr_BS_.push_back(pidxyErrBS);
  pi_sigdxy_BS_.push_back(safeRatio(std::abs(pidxyBS), pidxyErrBS));
  pi_dxy_PV_.push_back(pidxyPV);
  pi_dxyErr_PV_.push_back(pidxyErrPV);
  pi_sigdxy_PV_.push_back(safeRatio(std::abs(pidxyPV), pidxyErrPV));
  pi_norm_chi2_.push_back(pi.normalizedChi2());
  pi_N_valid_hits_.push_back(pi.numberOfValidHits());
  pi_lostInnerHits_.push_back(pi.numberOfLostHits());
  addP4(p4FromTrack(pi, massPi_), pi_pt_, pi_pz_, pi_E_, pi_eta_, pi_phi_);

  massb_piK_.push_back((p4FromTrack(pi, massPi_) + p4FromTrack(K, massK_)).M());
  massb_KK_.push_back((p4FromTrack(pi, massK_) + p4FromTrack(K, massK_)).M());
  massb_Kpi_conj_.push_back((p4FromTrack(pi, massK_) + p4FromTrack(K, massPi_)).M());

  chi2_piK_.push_back(resPiK.chi2);
  dof_piK_.push_back(resPiK.dof);
  pval_piK_.push_back(resPiK.pval);
  mass_piK_.push_back(massPiK);
  addP4(p4FromKinParticle(D0), D0_pt_, D0_pz_, D0_E_, D0_eta_, D0_phi_);

  D0KinTree->movePointerToTheFirstChild();
  const TLorentzVector refitK = p4FromKinParticle(D0KinTree->currentParticle());
  addP4(refitK, K_refitpiK_pt_, K_refitpiK_pz_, K_refitpiK_E_, K_refitpiK_eta_, K_refitpiK_phi_);
  D0KinTree->movePointerToTheNextChild();
  const TLorentzVector refitPi = p4FromKinParticle(D0KinTree->currentParticle());
  addP4(refitPi, pi_refitpiK_pt_, pi_refitpiK_pz_, pi_refitpiK_E_, pi_refitpiK_eta_, pi_refitpiK_phi_);
  mass_piK_hKK_.push_back(massWithHypotheses(refitPi, massK_, refitK, massK_));
  mass_piK_hpipi_.push_back(massWithHypotheses(refitPi, massPi_, refitK, massPi_));
  mass_piK_hKpi_.push_back(massWithHypotheses(refitPi, massK_, refitK, massPi_));

  const auto dD0PV = vertexDistance(primaryVtx, vtxD0);
  const auto dxyD0BS = transverseDistanceFromBeamSpot(beamSpot, vtxD0);
  cos_D0_PV_.push_back(pointingCos(primaryVtx, vtxD0, D0));
  cosT_D0_PV_.push_back(pointingCosTransverse(primaryVtx, vtxD0, D0));
  d_vtxD0_PV_.push_back(dD0PV.first);
  sigd_vtxD0_PV_.push_back(safeRatio(dD0PV.first, dD0PV.second));
  dxy_vtxD0_BS_.push_back(dxyD0BS.first);
  sigdxy_vtxD0_BS_.push_back(safeRatio(dxyD0BS.first, dxyD0BS.second));

  pis_dxy_BS_.push_back(pisdxyBS);
  pis_dxyErr_BS_.push_back(pisdxyErrBS);
  pis_sigdxy_BS_.push_back(safeRatio(std::abs(pisdxyBS), pisdxyErrBS));
  pis_dxy_PV_.push_back(pisdxyPV);
  pis_dxyErr_PV_.push_back(pisdxyErrPV);
  pis_sigdxy_PV_.push_back(safeRatio(std::abs(pisdxyPV), pisdxyErrPV));
  pis_norm_chi2_.push_back(pis.normalizedChi2());
  pis_N_valid_hits_.push_back(pis.numberOfValidHits());
  pis_lostInnerHits_.push_back(pis.numberOfLostHits());
  addP4(p4FromTrack(pis, massPi_), pis_pt_, pis_pz_, pis_E_, pis_eta_, pis_phi_);

  chi2_D0pis_.push_back(resD0Pis.chi2);
  dof_D0pis_.push_back(resD0Pis.dof);
  pval_D0pis_.push_back(resD0Pis.pval);
  mass_D0pis_.push_back(massD0Pis);
  deltaMass_D0pis_.push_back(massD0Pis - massPiK);
  const auto dDstPV = vertexDistance(primaryVtx, vtxDst);
  const auto dxyDstPV = transverseDistance(primaryVtx, vtxDst);
  cos_Dst_PV_.push_back(pointingCos(primaryVtx, vtxDst, Dst));
  cosT_Dst_PV_.push_back(pointingCosTransverse(primaryVtx, vtxDst, Dst));
  d_vtxDst_PV_.push_back(dDstPV.first);
  sigd_vtxDst_PV_.push_back(safeRatio(dDstPV.first, dDstPV.second));
  dxy_vtxDst_PV_.push_back(dxyDstPV.first);
  sigdxy_vtxDst_PV_.push_back(safeRatio(dxyDstPV.first, dxyDstPV.second));

  DstKinTree->movePointerToTheFirstChild();
  addP4(p4FromKinParticle(DstKinTree->currentParticle()),
        pis_refitD0pis_pt_,
        pis_refitD0pis_pz_,
        pis_refitD0pis_E_,
        pis_refitD0pis_eta_,
        pis_refitD0pis_phi_);
  DstKinTree->movePointerToTheNextChild();
  addP4(p4FromKinParticle(DstKinTree->currentParticle()),
        D0_refitD0pis_pt_,
        D0_refitD0pis_pz_,
        D0_refitD0pis_E_,
        D0_refitD0pis_eta_,
        D0_refitD0pis_phi_);

  (void)KdzBS;
  (void)pidzBS;
  n_BestDst_ = 1;
  return true;
}

void DstToD0PiProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  run_ = iEvent.run();
  lumi_ = iEvent.luminosityBlock();
  event_ = iEvent.id().event();
  clearTreeVectors();
  counters_[0]++;
  auto selectedDstarTracks = std::make_unique<reco::TrackCollection>();

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(trackToken_, tracks);
  edm::Handle<std::vector<reco::Vertex>> vertices;
  iEvent.getByToken(vertexToken_, vertices);
  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByToken(beamSpotToken_, beamSpot);

  if (!tracks.isValid() || !vertices.isValid() || vertices->empty() || !beamSpot.isValid()) {
    if (!deferTreeFill_)
      tree_->Fill();
    iEvent.put(std::move(selectedDstarTracks), "dstarTracks");
    return;
  }

  const reco::Vertex* primaryVtx = nullptr;
  for (const auto& vertex : *vertices) {
    if (!vertex.isFake()) {
      primaryVtx = &vertex;
      break;
    }
  }
  if (!primaryVtx)
    primaryVtx = &vertices->front();

  std::vector<char> valid(tracks->size(), 0);
  for (unsigned int i = 0; i < tracks->size(); ++i) {
    if (passTrackPreselection((*tracks)[i]))
      valid[i] = 1;
  }

  BestCandidate bestCandidate;
  for (unsigned int iK = 0; iK < tracks->size(); ++iK) {
    if (!valid[iK])
      continue;
    const reco::Track& K = (*tracks)[iK];
    const int chargeK = K.charge();
    ++n_K_;
    counters_[3]++;

    for (unsigned int iPi = 0; iPi < tracks->size(); ++iPi) {
      if (!valid[iPi] || iPi == iK)
        continue;
      const reco::Track& pi = (*tracks)[iPi];
      if (pi.charge() != -chargeK)
        continue;

      ++n_pi_;
      counters_[4]++;

      RefCountedKinematicTree D0KinTree = fitD0(iSetup, pi, K);
      const FitQuality resPiK = fitQuality(D0KinTree);
      if (!resPiK.isGood)
        continue;
      counters_[5]++;

      D0KinTree->movePointerToTheTop();
      const double massPiK = D0KinTree->currentParticle()->currentState().mass();
      if (std::abs(massPiK - massD0_) > dmD0Max_)
        continue;
      counters_[6]++;

      const RefCountedKinematicParticle D0 = D0KinTree->currentParticle();
      counters_[7]++;
      ++n_D0_;

      for (unsigned int iPis = 0; iPis < tracks->size(); ++iPis) {
        if (iPis == iK || iPis == iPi)
          continue;
        const reco::Track& pis = (*tracks)[iPis];
        if (!passSoftPionPreselection(pis))
          continue;
        if (pis.charge() != pi.charge())
          continue;

        counters_[8]++;
        ++n_pis_;

        RefCountedKinematicTree DstKinTree = fitDst(iSetup, pis, D0);
        const FitQuality resD0Pis = fitQuality(DstKinTree);
        if (!resD0Pis.isGood)
          continue;
        counters_[9]++;

        DstKinTree->movePointerToTheTop();
        const RefCountedKinematicParticle Dst = DstKinTree->currentParticle();
        const double massD0Pis = Dst->currentState().mass();
        if (std::abs(massD0Pis - massDst_) > dmDstMax_)
          continue;
        if (std::abs(massD0Pis - massPiK - (massDst_ - massD0_)) > dmDstMinusD0Max_)
          continue;

        counters_[10]++;
        ++n_Dst_;
        const double dstPt = Dst->currentState().globalMomentum().perp();
        if (!bestCandidate.found || dstPt > bestCandidate.pt) {
          bestCandidate = {true, iK, iPi, iPis, dstPt};
        }
      }
    }
  }

  if (bestCandidate.found) {
    const bool bestCandidateFilled = fillBestCandidate(iSetup,
                                                       (*tracks)[bestCandidate.iK],
                                                       (*tracks)[bestCandidate.iPi],
                                                       (*tracks)[bestCandidate.iPis],
                                                       *primaryVtx,
                                                       *beamSpot);
    if (bestCandidateFilled) {
      K_index_.push_back(bestCandidate.iK);
      pi_index_.push_back(bestCandidate.iPi);
      pis_index_.push_back(bestCandidate.iPis);
      selectedDstarTracks->push_back((*tracks)[bestCandidate.iK]);
      selectedDstarTracks->push_back((*tracks)[bestCandidate.iPi]);
      selectedDstarTracks->push_back((*tracks)[bestCandidate.iPis]);
    }
  }

  if (verbose_ > 1) {
    std::cout << "DstToD0PiProducer event " << run_ << ":" << lumi_ << ":" << event_ << " nDst=" << n_Dst_
              << " best=" << n_BestDst_ << std::endl;
  }

  if (!deferTreeFill_)
    tree_->Fill();
  iEvent.put(std::move(selectedDstarTracks), "dstarTracks");
}

DEFINE_FWK_MODULE(DstToD0PiProducer);
