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

# # pattern recognition over digitized hits
# from Configurables import GGTF_tracking_dbscan_check

# GGTF_tracking = GGTF_tracking_dbscan_check(
#     "GGTF_tracking_dbscan_eval",
#     inputHits_CDC=["CDCHDigis"],
#     inputHits_VTXD=["VTXDDigis"],
#     inputHits_VTXIB=["VTXIBDigis"],
#     inputHits_VTXOB=["VTXOBDigis"],
#     inputHits_CDC_sim=["CDCHHits"],
#     inputHits_VTXD_sim=["VTXDCollection"],
#     inputHits_VTXIB_sim=["VTXIBCollection"],
#     inputHits_VTXOB_sim=["VTXOBCollection"],
#     inputModel=["inputModel"],
#     outputModel=["outputModel"],
#     outputClustering=["outputClustering"],
#     OutputLevel=INFO,
# )

# pattern recognition over digitized hits
from Configurables import GGTF_check

GGTF_c = GGTF_check(
    "GGTF_check",
    outputTracks=["CDCHTracks"],
    test=["test"],
    OutputLevel=INFO,
)

################ Application
from k4FWCore import ApplicationMgr

chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]

ApplicationMgr(
    TopAlg=[GGTF_c],
    EvtSel="NONE",
    ExtSvc=[EventDataSvc("EventDataSvc"), audsvc],
    StopOnSignal=True,
    EvtMax=2,
    OutputLevel=INFO,
)