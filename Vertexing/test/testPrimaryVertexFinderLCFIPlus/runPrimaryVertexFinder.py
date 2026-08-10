from Configurables import EventDataSvc, PrimaryVertexFinderLCFIPlus, UniqueIDGenSvc
from Gaudi.Configuration import INFO
from k4FWCore import ApplicationMgr, IOSvc


io_svc = IOSvc("IOSvc")
io_svc.Input = "primary_vertex_input.root"
io_svc.Output = "primary_vertex_output.root"

vertex_finder = PrimaryVertexFinderLCFIPlus(
    "PrimaryVertexFinderLCFIPlus",
    InputTracks=["InputTracks"],
    OutputPrimaryVertex=["PrimaryVertex"],
    TrackChi2Cut=25.0,
    OutputLevel=INFO,
)

ApplicationMgr(
    TopAlg=[vertex_finder],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc"), UniqueIDGenSvc("uidSvc")],
    StopOnSignal=True,
)
