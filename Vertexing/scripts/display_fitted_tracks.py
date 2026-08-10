#!/usr/bin/env python3
"""Create XY and YZ event displays for fitted EDM4hep tracks."""

import argparse
import math
import os
from pathlib import Path
import sys
import tempfile

default_matplotlib_config = Path.home() / ".config" / "matplotlib"
if "MPLCONFIGDIR" not in os.environ and not os.access(default_matplotlib_config, os.W_OK):
    os.environ["MPLCONFIGDIR"] = str(
        Path(tempfile.gettempdir()) / f"matplotlib-{os.getuid()}"
    )

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from podio import root_io  # noqa: E402


DEFAULT_INPUT = Path(
    "/afs/cern.ch/work/a/adevita/public/workDir/vertexingMeeting/"
    "out_vertex_lcfiplus.root"
)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Draw the XY and YZ projections of FittedTracks and reconstructed "
            "vertices for one EDM4hep event."
        )
    )
    parser.add_argument(
        "event_id",
        type=int,
        help="EventHeader event number (or zero-based frame entry with --entry)",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Input ROOT file (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--entry",
        action="store_true",
        help="Interpret event_id as a zero-based frame entry instead of an event number",
    )
    parser.add_argument(
        "--tracks",
        default="FittedTracks",
        help="Track collection name (default: FittedTracks)",
    )
    parser.add_argument(
        "--vertices",
        default="VertexCandidates",
        help="Vertex collection name (default: VertexCandidates)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output PNG path (default: event_display_<event_id>.png)",
    )
    return parser.parse_args()


def event_number(frame, entry):
    collections = frame.getAvailableCollections()
    if "EventHeader" not in collections:
        return entry

    headers = frame.get("EventHeader")
    if len(headers) == 0:
        return entry
    return int(headers[0].getEventNumber())


def select_event(input_path, requested_id, select_by_entry):
    reader = root_io.Reader(str(input_path))
    available_ids = []

    for entry, frame in enumerate(reader.get("events")):
        current_event_number = event_number(frame, entry)
        available_ids.append(entry if select_by_entry else current_event_number)
        if (select_by_entry and entry == requested_id) or (
            not select_by_entry and current_event_number == requested_id
        ):
            return frame, entry, current_event_number

    selector = "entry" if select_by_entry else "event number"
    shown_ids = ", ".join(str(value) for value in available_ids[:20])
    if len(available_ids) > 20:
        shown_ids += ", ..."
    raise ValueError(
        f"No event with {selector} {requested_id}. Available values: {shown_ids or 'none'}"
    )


def coordinate(value, name, index):
    """Read an EDM vector component across EDM4hep/PyROOT binding variants."""
    if hasattr(value, name):
        return float(getattr(value, name))
    return float(value[index])


def track_points(track):
    """Return ordered trajectory points, prepending the fitted IP perigee."""
    points = []

    for state in track.getTrackStates():
        if state.location != 1:
            continue
        reference = state.referencePoint
        reference_x = coordinate(reference, "x", 0)
        reference_y = coordinate(reference, "y", 1)
        reference_z = coordinate(reference, "z", 2)
        points.append(
            (
                reference_x - float(state.D0) * math.sin(float(state.phi)),
                reference_y + float(state.D0) * math.cos(float(state.phi)),
                reference_z + float(state.Z0),
            )
        )
        break

    for hit in track.getTrackerHits():
        position = hit.getPosition()
        points.append(
            (
                coordinate(position, "x", 0),
                coordinate(position, "y", 1),
                coordinate(position, "z", 2),
            )
        )

    return points


def vertex_label(vertex, index):
    if hasattr(vertex, "isPrimary") and vertex.isPrimary():
        return f"PV {index}", "*", "red"
    if hasattr(vertex, "isSecondary") and vertex.isSecondary():
        return f"SV {index}", "X", "darkorange"
    return f"V{index}", "D", "black"


def draw_event(
    frame, entry, selected_event_number, track_collection, vertex_collection, output_path
):
    available_collections = frame.getAvailableCollections()
    if track_collection not in available_collections:
        raise KeyError(
            f"Track collection '{track_collection}' is absent. "
            f"Available collections: {', '.join(available_collections)}"
        )

    tracks = frame.get(track_collection)
    vertices = None
    if vertex_collection in available_collections:
        vertices = frame.get(vertex_collection)

    figure, (axis_xy, axis_yz) = plt.subplots(1, 2, figsize=(14, 6.5))
    color_map = plt.get_cmap("tab20")
    drawn_tracks = 0
    skipped_tracks = 0
    show_track_legend = len(tracks) <= 20

    for index, track in enumerate(tracks):
        points = track_points(track)
        if len(points) < 2:
            skipped_tracks += 1
            continue

        x = [point[0] for point in points]
        y = [point[1] for point in points]
        z = [point[2] for point in points]
        color = color_map(index % color_map.N)
        label = f"track {index}" if show_track_legend else None

        axis_xy.plot(x, y, color=color, linewidth=1.0, alpha=0.9, label=label)
        axis_yz.plot(y, z, color=color, linewidth=1.0, alpha=0.9, label=label)
        axis_xy.scatter(x[1:], y[1:], color=[color], s=3, alpha=0.45)
        axis_yz.scatter(y[1:], z[1:], color=[color], s=3, alpha=0.45)
        drawn_tracks += 1

    drawn_vertices = 0
    if vertices is not None:
        for index, vertex in enumerate(vertices):
            position = vertex.getPosition()
            x = coordinate(position, "x", 0)
            y = coordinate(position, "y", 1)
            z = coordinate(position, "z", 2)
            label, marker, color = vertex_label(vertex, index)

            axis_xy.scatter(
                [x], [y], marker=marker, color=color, edgecolors="white", s=180, zorder=10,
                linewidths=0.8, label=label,
            )
            axis_yz.scatter(
                [y], [z], marker=marker, color=color, edgecolors="white", s=180, zorder=10,
                linewidths=0.8, label=label,
            )
            axis_xy.annotate(label, (x, y), xytext=(6, 6), textcoords="offset points", color=color)
            axis_yz.annotate(label, (y, z), xytext=(6, 6), textcoords="offset points", color=color)
            drawn_vertices += 1

    for axis in (axis_xy, axis_yz):
        axis.axhline(0.0, color="0.75", linewidth=0.7, linestyle="--", zorder=0)
        axis.axvline(0.0, color="0.75", linewidth=0.7, linestyle="--", zorder=0)
        axis.grid(True, color="0.9", linewidth=0.5)
        axis.set_aspect("equal", adjustable="datalim")

    axis_xy.set_xlabel("x [mm]")
    axis_xy.set_ylabel("y [mm]")
    axis_xy.set_title("XY projection")
    axis_yz.set_xlabel("y [mm]")
    axis_yz.set_ylabel("z [mm]")
    axis_yz.set_title("YZ projection")

    if show_track_legend or drawn_vertices:
        axis_xy.legend(loc="best", fontsize="small", ncols=2)
        axis_yz.legend(loc="best", fontsize="small", ncols=2)

    figure.suptitle(
        f"Event {selected_event_number} (entry {entry}) — {track_collection}", fontsize=14
    )
    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(figure)

    print(
        f"Wrote {output_path} with {drawn_tracks} tracks and {drawn_vertices} vertices"
        + (
            f" ({skipped_tracks} tracks had fewer than two drawable points)"
            if skipped_tracks
            else ""
        )
    )


def main():
    arguments = parse_arguments()
    if not arguments.input.is_file():
        raise FileNotFoundError(f"Input file does not exist: {arguments.input}")

    output_path = arguments.output
    if output_path is None:
        output_path = Path(f"event_display_{arguments.event_id}.png")

    frame, entry, selected_event_number = select_event(
        arguments.input, arguments.event_id, arguments.entry
    )
    draw_event(
        frame,
        entry,
        selected_event_number,
        arguments.tracks,
        arguments.vertices,
        output_path,
    )


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, KeyError, RuntimeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1) from error
