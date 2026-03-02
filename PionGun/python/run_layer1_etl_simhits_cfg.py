import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")

process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring("file:step3.root")
)

process.MessageLogger = cms.Service("MessageLogger",
    destinations=cms.untracked.vstring("cout"),
    cout=cms.untracked.PSet(threshold=cms.untracked.string("INFO")),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("etl_primary_simhit_spreads.root")
)

process.mtdSimHitPrimary = cms.EDAnalyzer("MTDSimHitPrimaryAnalyzer",
    etlSimHits   = cms.InputTag("g4SimHits", "FastTimerHitsEndcap", "SIM"),
    btlSimHits   = cms.InputTag("g4SimHits", "FastTimerHitsBarrel", "SIM"),
    simTracks    = cms.InputTag("g4SimHits", "", "SIM"),
    simVertices  = cms.InputTag("g4SimHits", "", "SIM"),
)

process.p = cms.Path(process.mtdSimHitPrimary)