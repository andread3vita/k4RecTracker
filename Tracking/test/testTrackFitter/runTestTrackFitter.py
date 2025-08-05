import os

from Gaudi.Configuration import INFO, WARNING, DEBUG

from Configurables import AuditorSvc, ChronoAuditor
from Configurables import EventDataSvc


################ parser
from k4FWCore.parseArgs import parser

parser.add_argument("--inputFile", default="output_tracks.root", help="InputFile")
parser.add_argument("--outputFile", default="testFitter.root", help="OutputFile")
parser.add_argument("--Beta_init", default=100, help="OutputFile")
parser.add_argument("--Beta_final", default=0.1, help="OutputFile")
parser.add_argument("--Beta_steps", default=10, help="OutputFile")
args = parser.parse_args()

################ input & output
from k4FWCore import IOSvc

io_svc = IOSvc("IOSvc")
io_svc.Input = args.inputFile
io_svc.Output = args.outputFile

# pattern recognition over digitized hits
from Configurables import GenfitTrackFitter

trackFitter = GenfitTrackFitter(
    "GenfitTrackFitter",
    
    tracks_input=["CDCHTracks"],
    Fitted_tracks_electron=["Fitted_tracks_electron"],
    Fitted_tracks_muon=["Fitted_tracks_muon"],
    Fitted_tracks_pion=["Fitted_tracks_pion"],
    Fitted_tracks_kaon=["Fitted_tracks_kaon"],
    Fitted_tracks_proton=["Fitted_tracks_proton"],
    OutputLevel=DEBUG,
)

trackFitter.Beta_init = args.Beta_init
trackFitter.Beta_final = args.Beta_final
trackFitter.Beta_steps = args.Beta_steps
trackFitter.debug_lvl = 0

# ################ Detector geometry
from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [os.environ["K4GEO"]+"/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml"]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False
    
################ Application
from k4FWCore import ApplicationMgr
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
ApplicationMgr(
    TopAlg=[trackFitter],
    EvtSel="NONE",
    ExtSvc=[geoservice,EventDataSvc("EventDataSvc"), audsvc],
    StopOnSignal=True,
    EvtMax=-1,
    OutputLevel=INFO,
)