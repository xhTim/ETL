#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "TH1D.h"
#include "TH2D.h"

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <vector>

class MTDSimHitPrimaryAnalyzer : public edm::one::EDAnalyzer<> {
public:
  explicit MTDSimHitPrimaryAnalyzer(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  struct TrackAccum {
    unsigned nETL = 0;
    double tMin = 1e99;
    double tMax = -1e99;
    bool hasHit = false;

    // local-position proxy (midpoint of entry/exit in local coords)
    double x0 = 0., y0 = 0., z0 = 0.;
    bool hasRefPos = false;
    double maxDrLocalFromFirst = 0.;
    double zLocalMin = 1e99;
    double zLocalMax = -1e99;
  };

  edm::EDGetTokenT<std::vector<PSimHit>> tokETL_;
  edm::EDGetTokenT<std::vector<PSimHit>> tokBTL_;  // keep interface, do not count
  edm::EDGetTokenT<edm::SimTrackContainer> tokSimTrk_;
  edm::EDGetTokenT<edm::SimVertexContainer> tokSimVtx_;

  TH1D* h_nETLPrimary_;
  TH1D* h_timeSpreadETLPrimary_;
  TH1D* h_spaceSpreadDrLocalETLPrimary_;
  TH1D* h_spaceSpreadDzLocalETLPrimary_;
  TH2D* h_timeVsDrLocalETLPrimary_;
};

MTDSimHitPrimaryAnalyzer::MTDSimHitPrimaryAnalyzer(const edm::ParameterSet& iConfig) {
  tokETL_ = consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("etlSimHits"));
  tokBTL_ = consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("btlSimHits"));
  tokSimTrk_ = consumes<edm::SimTrackContainer>(iConfig.getParameter<edm::InputTag>("simTracks"));
  tokSimVtx_ = consumes<edm::SimVertexContainer>(iConfig.getParameter<edm::InputTag>("simVertices"));

  edm::Service<TFileService> fs;
  h_nETLPrimary_ = fs->make<TH1D>(
      "h_nETLPrimary",
      "Primary-track ETL PSimHit multiplicity;N_{ETL}^{primary};Tracks",
      20, -0.5, 19.5);

  h_timeSpreadETLPrimary_ = fs->make<TH1D>(
      "h_timeSpreadETLPrimary",
      "Primary-track ETL time spread;#Delta t=t_{max}-t_{min} [ns];Tracks with N_{ETL}#geq1",
      200, 0.0, 10.0);

  h_spaceSpreadDrLocalETLPrimary_ = fs->make<TH1D>(
      "h_spaceSpreadDrLocalETLPrimary",
      "Primary-track ETL local-space spread proxy;max #Delta r_{local} from first hit;Tracks with N_{ETL}#geq1",
      200, 0.0, 50.0);

  h_spaceSpreadDzLocalETLPrimary_ = fs->make<TH1D>(
      "h_spaceSpreadDzLocalETLPrimary",
      "Primary-track ETL local-z spread proxy;#Delta z_{local}=z_{max}-z_{min};Tracks with N_{ETL}#geq1",
      200, 0.0, 20.0);

  h_timeVsDrLocalETLPrimary_ = fs->make<TH2D>(
      "h_timeVsDrLocalETLPrimary",
      "Primary-track ETL time vs local-space spread proxy;#Delta t [ns];max #Delta r_{local}",
      200, 0.0, 10.0, 200, 0.0, 50.0);
}

void MTDSimHitPrimaryAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  edm::Handle<edm::SimTrackContainer> hTrk;
  edm::Handle<std::vector<PSimHit>> hETL;
  edm::Handle<std::vector<PSimHit>> hBTL;

  iEvent.getByToken(tokSimTrk_, hTrk);
  iEvent.getByToken(tokETL_, hETL);
  iEvent.getByToken(tokBTL_, hBTL);  // keep interface intentionally (not counted)

  // trackId -> (pdgId, isPrimaryApprox)
  std::unordered_map<unsigned, std::pair<int, bool>> trkInfo;
  trkInfo.reserve(hTrk->size());

  // primary trackId -> accumulators (ETL only)
  // IMPORTANT: pre-create all primary tracks so nETL=0 tracks are included.
  std::unordered_map<unsigned, TrackAccum> acc;
  acc.reserve(hTrk->size());

  for (auto const& trk : *hTrk) {
    bool isPrimary = (trk.genpartIndex() >= 0);  // particle gun下通常够用
    const unsigned tid = trk.trackId();
    trkInfo.emplace(tid, std::make_pair(trk.type(), isPrimary));
    if (isPrimary) {
      acc.emplace(tid, TrackAccum{});  // create even if it gets zero ETL hits
    }
  }

  // Fill ETL hit information for primary-assigned tracks
  for (auto const& hit : *hETL) {
    unsigned tid = hit.trackId();
    auto it = trkInfo.find(tid);
    if (it == trkInfo.end()) continue;
    if (!it->second.second) continue;  // only primary (current definition)

    auto ait = acc.find(tid);
    if (ait == acc.end()) continue;  // defensive (should not happen)
    auto& a = ait->second;

    a.nETL++;
    a.hasHit = true;

    const double t = hit.timeOfFlight();
    a.tMin = std::min(a.tMin, t);
    a.tMax = std::max(a.tMax, t);

    // local midpoint proxy (entry/exit)
    const auto& pIn  = hit.entryPoint();
    const auto& pOut = hit.exitPoint();
    const double x = 0.5 * (pIn.x() + pOut.x());
    const double y = 0.5 * (pIn.y() + pOut.y());
    const double z = 0.5 * (pIn.z() + pOut.z());

    if (!a.hasRefPos) {
      a.x0 = x; a.y0 = y; a.z0 = z;
      a.hasRefPos = true;
    }

    const double dr = std::sqrt((x - a.x0) * (x - a.x0) + (y - a.y0) * (y - a.y0));
    a.maxDrLocalFromFirst = std::max(a.maxDrLocalFromFirst, dr);
    a.zLocalMin = std::min(a.zLocalMin, z);
    a.zLocalMax = std::max(a.zLocalMax, z);
  }

  edm::LogPrint("MTDSimHitPrimaryAnalyzer")
      << "Run " << iEvent.id().run() << " Event " << iEvent.id().event()
      << " : primary-track ETL simhit multiplicity + spreads (trackId: pdgId, nETL, dT[ns], dRlocal, dZlocal)";

  if (acc.empty()) {
    edm::LogPrint("MTDSimHitPrimaryAnalyzer") << "  (no primary tracks found by current definition)";
    return;
  }

  for (auto const& kv : acc) {
    const unsigned tid = kv.first;
    const auto& a = kv.second;
    const int pdg = trkInfo[tid].first;

    // Multiplicity includes nETL=0 primaries (important!)
    h_nETLPrimary_->Fill(a.nETL);

    if (a.hasHit) {
      const double dT = a.tMax - a.tMin;
      const double dRlocal = a.maxDrLocalFromFirst;
      const double dZlocal = a.hasRefPos ? (a.zLocalMax - a.zLocalMin) : 0.0;

      h_timeSpreadETLPrimary_->Fill(dT);
      h_spaceSpreadDrLocalETLPrimary_->Fill(dRlocal);
      h_spaceSpreadDzLocalETLPrimary_->Fill(dZlocal);
      h_timeVsDrLocalETLPrimary_->Fill(dT, dRlocal);

      edm::LogPrint("MTDSimHitPrimaryAnalyzer")
          << "  trackId=" << tid << " pdgId=" << pdg
          << " nETL=" << a.nETL
          << " dT=" << dT
          << " dRlocal=" << dRlocal
          << " dZlocal=" << dZlocal;
    } else {
      edm::LogPrint("MTDSimHitPrimaryAnalyzer")
          << "  trackId=" << tid << " pdgId=" << pdg
          << " nETL=0 (no ETL simhits)";
    }
  }
}

DEFINE_FWK_MODULE(MTDSimHitPrimaryAnalyzer);