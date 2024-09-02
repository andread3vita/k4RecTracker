import os

from Gaudi.Configuration import INFO, WARNING

from Configurables import AuditorSvc, ChronoAuditor
from Configurables import EventDataSvc
from Configurables import HiveSlimEventLoopMgr, HiveWhiteBoard, AvalancheSchedulerSvc

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

GGTF_tracking = GGTF_tracking_dbscan(
    "GGTF_tracking_dbscan_eval",
    inputHits_CDC=["CDCHDigis"],
    inputHits_VTXD=["VTXDDigis"],
    inputHits_VTXIB=["VTXIBDigis"],
    inputHits_VTXOB=["VTXOBDigis"],
    outputTracks=["CDCHTracks"],
    OutputLevel=INFO,
)


GGTF_tracking.modelPath = "/afs/cern.ch/user/a/adevita/public/workDir/k4RecTracker/Tracking/model_multivector_1_input.onnx"

################ Application
from k4FWCore import ApplicationMgr

chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]

ApplicationMgr(
    TopAlg=[GGTF_tracking],
    EvtSel="NONE",
    ExtSvc=[EventDataSvc("EventDataSvc"), audsvc],
    StopOnSignal=True,
    EvtMax=2,
    OutputLevel=INFO,
)