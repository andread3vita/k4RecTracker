import os

from Gaudi.Configuration import INFO, WARNING

from Configurables import AuditorSvc, ChronoAuditor
from Configurables import EventDataSvc
from Configurables import HiveSlimEventLoopMgr, HiveWhiteBoard, AvalancheSchedulerSvc

# evtslots = 10
# threads = 12

# whiteboard = HiveWhiteBoard(
#     "EventDataSvc",
#     EventSlots=evtslots,
#     ForceLeaves=True,
# )

# slimeventloopmgr = HiveSlimEventLoopMgr(
#     "HiveSlimEventLoopMgr", SchedulerName="AvalancheSchedulerSvc", OutputLevel=WARNING
# )

# scheduler = AvalancheSchedulerSvc(ThreadPoolSize=threads, ShowDataFlow=True, OutputLevel=WARNING)

################ parser
from k4FWCore.parseArgs import parser

parser.add_argument("--inputFile", default="/afs/cern.ch/user/a/adevita/public/workDir/test/dataset/muons_eta0017_mom10/eval/output_IDEA_tracking_gun_10.root", help="InputFile")
parser.add_argument("--outputFile", default="output_eff.root", help="OutputFile")
args = parser.parse_args()

################ input & output
from k4FWCore import IOSvc

io_svc = IOSvc("IOSvc")
io_svc.input = args.inputFile
io_svc.output = args.outputFile


# pattern recognition over digitized hits
from Configurables import GGTF_tracking_dbscan

GGTF = GGTF_tracking_dbscan(
    "GGTF_tracking_dbscan",
    inputHits_CDC=["CDCHDigis"],
    inputHits_VTXD=["VTXDDigis"],
    inputHits_VTXIB=["VTXIBDigis"],
    inputHits_VTXOB=["VTXOBDigis"],
    inputHits_CDC_sim=["CDCHHits"],
    inputHits_VTXD_sim=["VTXDCollection"],
    inputHits_VTXIB_sim=["VTXIBCollection"],
    inputHits_VTXOB_sim=["VTXOBCollection"],
    outputTracks=["CDCHTracks"],
    outputHits=["outputHits"],
    OutputLevel=INFO,
)
GGTF.modelPath = "/afs/cern.ch/user/a/adevita/public/workDir/k4RecTracker/Tracking/model_multivector_1_input.onnx"

################ Application
from k4FWCore import ApplicationMgr

chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]

# ApplicationMgr(
#     TopAlg=[GGTF],
#     EvtSel="NONE",
#     StopOnSignal=True,
#     EvtMax=-1,
#     ExtSvc=[whiteboard],
#     EventLoop=slimeventloopmgr,
#     MessageSvcType="InertMessageSvc",
#     OutputLevel=INFO,
# )

ApplicationMgr(
    TopAlg=[GGTF],
    EvtSel="NONE",
    ExtSvc=[EventDataSvc("EventDataSvc"), audsvc],
    StopOnSignal=True,
    EvtMax=-1,
    OutputLevel=INFO,
)