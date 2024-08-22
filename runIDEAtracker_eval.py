import os
from Gaudi.Configuration import *

# Loading the input SIM file
from Configurables import k4DataSvc, PodioInput
from k4FWCore.parseArgs import parser

parser.add_argument("--inputFile", default="../dataset/muons_eta0017_mom100/output_IDEA_DIGI_gun_1.root", help="InputFile")
parser.add_argument("--outputFile", default="output_tracking.root", help="OutputFile")
args = parser.parse_args()

evtsvc = k4DataSvc("EventDataSvc")
evtsvc.input = args.inputFile
inp = PodioInput("InputReader")


# pattern recognition over digitized hits
from Configurables import GGTF_tracking_dbscan_eval

GGTF_tracking = GGTF_tracking_dbscan_eval(
    "GGTF_tracking_dbscan_eval",
    modelPath="/afs/cern.ch/user/a/adevita/public/workDir/k4RecTracker/Tracking/model_multivector_1_input.onnx",
    inputHits_CDC="CDCHDigis",
    inputHits_VTXD="VTXDDigis",
    inputHits_VTXIB="VTXIBDigis",
    inputHits_VTXOB="VTXOBDigis",
    inputHits_CDC_sim="CDCHHits",
    inputHits_VTXD_sim="VTXDCollection",
    inputHits_VTXIB_sim="VTXIBCollection",
    inputHits_VTXOB_sim="VTXOBCollection",
    inputAssociation_CDC_sim="CDCHDigisAssociation",
    outputTracks="CDCHTracks",
    clustering_space_tracks="clustering_space_tracks",
    OutputLevel=INFO,
)

################ Output
from Configurables import PodioOutput

out = PodioOutput("out", OutputLevel=INFO)
out.outputCommands = ["keep *"]

out.filename = args.outputFile

# CPU information
from Configurables import AuditorSvc, ChronoAuditor

chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
out.AuditExecute = True

from Configurables import ApplicationMgr

ApplicationMgr(
    TopAlg=[
        inp,
        GGTF_tracking,
        out,
    ],
    EvtSel="NONE",
    ExtSvc=[evtsvc, audsvc],
    StopOnSignal=True,
)