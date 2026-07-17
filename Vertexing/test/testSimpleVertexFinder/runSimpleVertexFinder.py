from Gaudi.Configuration import INFO
from Configurables import EventDataSvc, SimpleVertexFinder
from k4FWCore import ApplicationMgr, IOSvc
from k4FWCore.parseArgs import parser


parser.add_argument("--inputFile", help="Input file containing fitted tracks")
parser.add_argument("--outputFile", help="Output file containing vertex candidates")
args = parser.parse_args()

io_svc = IOSvc("IOSvc")
io_svc.Input = args.inputFile
io_svc.Output = args.outputFile

vertex_finder = SimpleVertexFinder(
    "SimpleVertexFinder",
    InputFittedTracks=["FittedTracks"],
    OutputVerticesCandidates=["VertexCandidates"],
    OutputLevel=INFO,
)

ApplicationMgr(
    TopAlg=[vertex_finder],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc")],
    StopOnSignal=True,
)
