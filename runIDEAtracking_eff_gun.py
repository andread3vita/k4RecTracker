from Gaudi.Configuration import INFO, WARNING
from k4FWCore import ApplicationMgr
from Configurables import k4DataSvc, PodioOutput, PodioInput
from Configurables import Efficiency_calc_gun
from Configurables import AuditorSvc, ChronoAuditor
from k4FWCore.parseArgs import parser
from Configurables import HiveSlimEventLoopMgr, HiveWhiteBoard, AvalancheSchedulerSvc

import os

evtslots = 10
threads = 12

whiteboard = HiveWhiteBoard(
    "EventDataSvc",
    EventSlots=evtslots,
    ForceLeaves=True,
)

slimeventloopmgr = HiveSlimEventLoopMgr(
    "HiveSlimEventLoopMgr", SchedulerName="AvalancheSchedulerSvc", OutputLevel=WARNING
)

scheduler = AvalancheSchedulerSvc(ThreadPoolSize=threads, ShowDataFlow=True, OutputLevel=WARNING)

################ parser

parser.add_argument("--inputFile", default="/afs/cern.ch/user/a/adevita/public/workDir/test/dataset/muons_eta0017_mom10/eval/output_IDEA_tracking_gun_10.root", help="InputFile")
parser.add_argument("--outputFile", default="output_eff.root", help="OutputFile")
args = parser.parse_args()

################ input & output

from k4FWCore import IOSvc
io_svc = IOSvc("IOSvc")
io_svc.input = args.inputFile
io_svc.output = args.outputFile

################ eff

trackin_eff = Efficiency_calc_gun(
    "Evaluation_efficiency",
    InputCollectionTracks=["CDCHTracks"],
    InputCollectionParticles=["MCParticles"],
    inputHits_CDC_sim=["CDCHHits"],
    inputHits_VTXD_sim=["VTXDCollection"],
    inputHits_VTXIB_sim=["VTXIBCollection"],
    inputHits_VTXOB_sim=["VTXOBCollection"],
    pdg_MCParticles=["pdg_MCParticles"],
    pt_MCParticles=["pt_MCParticles"],
    index_MCParticles=["index_MCParticles"],
    isReco=["isReco"],
    isPrimary=["isPrimary"],
    outputTrackingEff_primary=["outputTrackingEff_primary"],
    outputTrackingEff_secondary=["outputTrackingEff_secondary"],
    OutputLevel=INFO,
)

################ Application
from Configurables import EventDataSvc

chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]

ApplicationMgr(
    TopAlg=[trackin_eff],
    EvtSel="NONE",
    # ExtSvc=[EventDataSvc("EventDataSvc"), audsvc],
    StopOnSignal=True,
    EvtMax=-1,
    ExtSvc=[whiteboard],
    EventLoop=slimeventloopmgr,
    MessageSvcType="InertMessageSvc",
    OutputLevel=INFO,
)

