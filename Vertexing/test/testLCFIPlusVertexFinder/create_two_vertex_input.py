"""Create one idealised event with a primary, secondary, and fitted K-short."""

import math
import edm4hep
import podio
from podio import root_io


def add_track(collection, z0, phi, tan_lambda, omega):
    track = collection.create()
    state = edm4hep.TrackState()
    state.location = edm4hep.TrackState.AtIP
    state.D0 = 0.0
    state.Z0 = z0
    state.phi = phi
    state.tanLambda = tan_lambda
    state.omega = omega

    # EDM4hep TrackState covariance ordering is lower triangular.  These
    # non-zero uncertainties make the toy states usable by the vertex fitter.
    covariance = [0.0] * 21
    covariance[0] = 1.0e-4
    covariance[2] = 1.0e-6
    covariance[5] = 1.0e-10
    covariance[9] = 1.0e-4
    covariance[14] = 1.0e-6
    state.covMatrix = covariance
    track.addToTrackStates(state)


tracks = edm4hep.TrackCollection()

# Primary vertex at (0, 0, 0) mm.
add_track(tracks, 0.0, -2.0, 0.2, 0.0010)
add_track(tracks, 0.0, 0.2, -0.1, 0.0011)
add_track(tracks, 0.0, 2.2, 0.3, 0.0009)

# K-short candidate at (0, 0, 20) mm.  Equal and opposite transverse
# momenta give a pion-pair mass of 0.497611 GeV; the positive longitudinal
# components make its momentum point away from the primary vertex.
ks_omega = 0.0029110082929612977
add_track(tracks, 20.0, 0.0, 1.0, ks_omega)
add_track(tracks, 20.0, math.pi, 1.0, -ks_omega)

# Secondary vertex at (0, 0, 10) mm.  The positive longitudinal momenta make
# the summed momentum point from the primary vertex towards this vertex.
add_track(tracks, 10.0, -1.0, 0.8, 0.0012)
add_track(tracks, 10.0, 0.7, 0.9, 0.0010)
add_track(tracks, 10.0, 2.5, 0.7, 0.0008)

frame = podio.Frame()
frame.put(tracks, "InputTracks")
writer = root_io.Writer("two_vertex_input.root")
writer.write_frame(frame, "events")
writer.finish()
