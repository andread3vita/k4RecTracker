from Gaudi.Configuration import INFO, WARNING
from k4FWCore import ApplicationMgr
from Configurables import k4DataSvc, PodioOutput, PodioInput
from Configurables import AuditorSvc, ChronoAuditor
from k4FWCore.parseArgs import parser
from Configurables import HiveSlimEventLoopMgr, HiveWhiteBoard, AvalancheSchedulerSvc

import os


################ parser

parser.add_argument("--inputFile", default="/afs/cern.ch/user/a/adevita/public/workDir/test/dataset/muons_eta0017_mom10/eval/output_IDEA_tracking_gun_10.root", help="InputFile")
parser.add_argument("--outputFile", default="output_eff.root", help="OutputFile")
args = parser.parse_args()

################ input & output

from k4FWCore import IOSvc
io_svc = IOSvc("IOSvc")
io_svc.input = args.inputFile
io_svc.output = args.outputFile

################## Simulation setup
# Detector geometry
from Configurables import GeoSvc

geoservice = GeoSvc("GeoSvc")
path_to_detector = os.environ.get("K4GEO", "")
print(path_to_detector)
detectors_to_use = ["/afs/cern.ch/user/a/adevita/public/workDir/k4geo/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml"]
# prefix all xmls with path_to_detector
geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
geoservice.OutputLevel = INFO

################ eff
from Configurables import DCHdigi
DCHdigi = DCHdigi("DCHdigi")
DCHdigi.DCH_simhits=["CDCHHits"]
DCHdigi.DCH_name="DCH_v2"
# DCHdigi.fileDataAlg="DataAlgFORGEANT.root"
DCHdigi.zResolution_mm=1
DCHdigi.xyResolution_mm=0.1
DCHdigi.OutputLevel=INFO

################ Application
from Configurables import EventDataSvc, UniqueIDGenSvc

chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]

ApplicationMgr(
    TopAlg=[DCHdigi],
    EvtSel="NONE",
    ExtSvc=[geoservice,EventDataSvc("EventDataSvc"),UniqueIDGenSvc("uidSvc")],
    StopOnSignal=True,
    EvtMax=-1,
    OutputLevel=INFO,
)