import os
import math 

from Gaudi.Configuration import INFO,DEBUG

from Configurables import EventDataSvc,UniqueIDGenSvc
from Configurables import RndmGenSvc
from Configurables import SimG4SaveTrackerHits
from Configurables import GeoSvc

from k4FWCore import IOSvc,ApplicationMgr
from k4FWCore.parseArgs import parser


################## Parser
parser.add_argument("--inputFile", default="out_sim_edm4hep.root", help="InputFile")
parser.add_argument("--outputFile", default="output_digi_v3.root", help="OutputFile")
args = parser.parse_args()

# ################## InputOutput
svc = IOSvc("IOSvc")
svc.input = args.inputFile
svc.output = args.outputFile

from Configurables import VTXdigitizer
from Configurables import DCHdigi_v01
    
################## Vertex sensor resolutions
innerVertexResolution_x = 0.003 # [mm], assume 3 µm resolution for ARCADIA sensor
innerVertexResolution_y = 0.003 # [mm], assume 3 µm resolution for ARCADIA sensor
innerVertexResolution_t = 1000 # [ns]
outerVertexResolution_x = 0.050/math.sqrt(12) # [mm], assume ATLASPix3 sensor with 50 µm pitch
outerVertexResolution_y = 0.150/math.sqrt(12) # [mm], assume ATLASPix3 sensor with 150 µm pitch
outerVertexResolution_t = 1000 # [ns]
    
################ Detector geometry
geoservice = GeoSvc("GeoSvc")
path_to_detector = os.environ.get("K4GEO", "")
print(path_to_detector)
detectors_to_use=['FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml']
geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
geoservice.OutputLevel = INFO
    
    
############### Vertex Digitizer
idea_vtxb_digitizer = VTXdigitizer("VTXBdigitizer",
    inputSimHits = "VertexBarrelCollection",
    outputDigiHits = "VertexBarrelCollection_digi",
    outputSimDigiAssociation = "VertexBarrel_Association",
    detectorName = "Vertex",
    readoutName = "VertexBarrelCollection",
    xResolution = [innerVertexResolution_x, innerVertexResolution_x, innerVertexResolution_x, outerVertexResolution_x, outerVertexResolution_x], # mm, r-phi direction
    yResolution = [innerVertexResolution_y, innerVertexResolution_y, innerVertexResolution_y, outerVertexResolution_y, outerVertexResolution_y], # mm, z direction
    tResolution = [innerVertexResolution_t, innerVertexResolution_t, innerVertexResolution_t, outerVertexResolution_t, outerVertexResolution_t], # ns
    forceHitsOntoSurface = False,
    OutputLevel = INFO
)

idea_vtxd_digitizer  = VTXdigitizer("VTXDdigitizer",
    inputSimHits = "VertexEndcapCollection",
    outputDigiHits = "VertexEndcapCollection_digi",
    outputSimDigiAssociation = "VertexEndcap_Association",
    detectorName = "Vertex",
    readoutName = "VertexEndcapCollection",
    xResolution = [outerVertexResolution_x, outerVertexResolution_x, outerVertexResolution_x], # mm, r direction
    yResolution = [outerVertexResolution_y, outerVertexResolution_y, outerVertexResolution_y], # mm, phi direction
    tResolution = [outerVertexResolution_t, outerVertexResolution_t, outerVertexResolution_t], # ns
    forceHitsOntoSurface = False,
    OutputLevel = INFO
)
        
################ DC digitizer
DCHdigi = DCHdigi_v01("DCHdigi")
DCHdigi.DCH_simhits=["DCHCollection"]
DCHdigi.DCH_name="DCH_v2"
# DCHdigi.fileDataAlg="DataAlgFORGEANT.root"
DCHdigi.calculate_dndx=True
DCHdigi.create_debug_histograms=False
DCHdigi.zResolution_mm=1
DCHdigi.xyResolution_mm=0.1
DCHdigi.OutputLevel=INFO
    
rndm_gen_svc = RndmGenSvc()
    
mgr = ApplicationMgr(TopAlg=[DCHdigi,idea_vtxb_digitizer,idea_vtxd_digitizer],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[geoservice,EventDataSvc("EventDataSvc"),UniqueIDGenSvc("uidSvc"),RndmGenSvc()],
    OutputLevel=INFO,
    )
