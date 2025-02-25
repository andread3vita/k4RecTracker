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
from Configurables import GGTF_tracking_dbscan_IDEAv3

GGTF = GGTF_tracking_dbscan_IDEAv3(
    "GGTF_tracking_dbscan",
    inputHits_CDC=["DCH_DigiCollection"],
    inputHits_VTXD=["VertexEndcapCollection_digi"],
    inputHits_VTXB=["VertexBarrelCollection_digi"],
    outputTracks=["CDCHTracks"],
    OutputLevel=INFO,
)
GGTF.modelPath = "/eos/user/a/adevita/workDir/k4RecTracker/Tracking/model_multivector_input_011124_v2.onnx"

################ Detector geometry
from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")
path_to_detector = os.environ.get("K4GEO", "")
print(path_to_detector)
detectors_to_use=['FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml']
geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
geoservice.OutputLevel = INFO
    
################ Application
from k4FWCore import ApplicationMgr

chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]

ApplicationMgr(
    TopAlg=[GGTF],
    EvtSel="NONE",
    ExtSvc=[geoservice,EventDataSvc("EventDataSvc"), audsvc],
    StopOnSignal=True,
    EvtMax=-1,
    OutputLevel=INFO,
)