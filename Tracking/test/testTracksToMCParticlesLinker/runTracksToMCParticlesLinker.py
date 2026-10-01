"""Run TracksToMCParticlesLinker on the IDEA_o1_v03 track-finder output."""

from Configurables import TracksToMCParticlesLinker
from k4FWCore import ApplicationMgr, IOSvc
from k4FWCore.parseArgs import parser


parser.add_argument("--input", required=True, help="Track-finder EDM4hep input file")
parser.add_argument("--output", required=True, help="Linked EDM4hep output file")
args = parser.parse_args()

io_svc = IOSvc("IOSvc")
io_svc.Input = args.input
io_svc.Output = args.output

linker = TracksToMCParticlesLinker(
    "TracksToMCParticlesLinker",
    TrackCollection=["GGTFTracks"],
    DigiToSimHitsLinks=[
        "VTXBSimDigiLinks",
        "VTXDSimDigiLinks",
        "SiWrBSimDigiLinks",
        "SiWrDSimDigiLinks",
        "DCH_DigiSimAssociationCollection",
    ],
    LinksTracksMCParticles=["TracksMCParticlesLinks"],
)

ApplicationMgr(
    TopAlg=[linker],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[io_svc],
)
