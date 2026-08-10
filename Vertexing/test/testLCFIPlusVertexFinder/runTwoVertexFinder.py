from Gaudi.Configuration import INFO
from Configurables import EventDataSvc, UniqueIDGenSvc
from k4FWCore import ApplicationMgr, IOSvc


svc = IOSvc("IOSvc")
svc.Input = "two_vertex_input.root"
svc.Output = "two_vertex_output.root"

from Configurables import LCFIPlusVertexFinder

vertex_finder = LCFIPlusVertexFinder(
    "LCFIPlusVertexFinder",
    InputFittedTracks=["InputTracks"],
    OutputVerticesCandidates=["VertexCandidates"],
    OutputV0Vertices=["V0Vertices"],
    PrimaryTrackChi2Cut=25.0,
    Chi2Cut=25.0,
    AddedTrackChi2Cut=25.0,
    RejectV0s=True,
    ReconstructV0Vertices=True,
    IncludeV0InVertexCandidates=True,
    OutputLevel=INFO,
)

ApplicationMgr(
    TopAlg=[vertex_finder],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc"), UniqueIDGenSvc("uidSvc")],
    StopOnSignal=True,
)
