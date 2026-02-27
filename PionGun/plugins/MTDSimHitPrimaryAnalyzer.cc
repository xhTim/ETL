#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include <vector>
#include <unordered_map>

class MTDSimHitPrimaryAnalyzer : public edm::one::EDAnalyzer<> {
public:
  explicit MTDSimHitPrimaryAnalyzer(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  edm::EDGetTokenT<std::vector<PSimHit>> tokETL_;
  edm::EDGetTokenT<std::vector<PSimHit>> tokBTL_;
  edm::EDGetTokenT<edm::SimTrackContainer> tokSimTrk_;
  edm::EDGetTokenT<edm::SimVertexContainer> tokSimVtx_;
};

MTDSimHitPrimaryAnalyzer::MTDSimHitPrimaryAnalyzer(const edm::ParameterSet& iConfig) {
  tokETL_ = consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("etlSimHits"));
  tokBTL_ = consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("btlSimHits"));
  tokSimTrk_ = consumes<edm::SimTrackContainer>(iConfig.getParameter<edm::InputTag>("simTracks"));
  tokSimVtx_ = consumes<edm::SimVertexContainer>(iConfig.getParameter<edm::InputTag>("simVertices"));
}

void MTDSimHitPrimaryAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  edm::Handle<edm::SimTrackContainer> hTrk;
  edm::Handle<std::vector<PSimHit>> hETL;
  edm::Handle<std::vector<PSimHit>> hBTL;

  iEvent.getByToken(tokSimTrk_, hTrk);
  iEvent.getByToken(tokETL_, hETL);
  iEvent.getByToken(tokBTL_, hBTL);

  // trackId -> (pdgId, isPrimaryApprox)
  std::unordered_map<unsigned, std::pair<int,bool>> trkInfo;
  trkInfo.reserve(hTrk->size());
  for (auto const& trk : *hTrk) {
    bool isPrimary = (trk.genpartIndex() >= 0); // particle gun 主粒子通常满足
    trkInfo.emplace(trk.trackId(), std::make_pair(trk.type(), isPrimary));
  }

  // primary trackId -> (nETL, nBTL)
  std::unordered_map<unsigned, std::pair<unsigned,unsigned>> cnt;

  auto countHits = [&](std::vector<PSimHit> const& hits, bool isETL) {
    for (auto const& hit : hits) {
      unsigned tid = hit.trackId();
      auto it = trkInfo.find(tid);
      if (it == trkInfo.end()) continue;
      if (!it->second.second) continue; // only primary
      auto& c = cnt[tid];
      if (isETL) c.first++;
      else c.second++;
    }
  };

  countHits(*hETL, true);
  countHits(*hBTL, false);

  edm::LogPrint("MTDSimHitPrimaryAnalyzer")
    << "Run " << iEvent.id().run() << " Event " << iEvent.id().event()
    << " : primary-track simhit counts (trackId: pdgId, nETL, nBTL)";

  if (cnt.empty()) {
    edm::LogPrint("MTDSimHitPrimaryAnalyzer") << "  (no primary-track simhits in ETL/BTL)";
    return;
  }

  for (auto const& kv : cnt) {
    unsigned tid = kv.first;
    int pdg = trkInfo[tid].first;
    edm::LogPrint("MTDSimHitPrimaryAnalyzer")
      << "  trackId=" << tid << " pdgId=" << pdg
      << " nETL=" << kv.second.first << " nBTL=" << kv.second.second;
  }
}

DEFINE_FWK_MODULE(MTDSimHitPrimaryAnalyzer);
