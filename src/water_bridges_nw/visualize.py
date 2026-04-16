import json
import logging
import os

logger = logging.getLogger(__name__)

def export_vmd_script(data, output_file="draw_pathways.tcl", mode="frame", frame_idx=None):
    """
    Generates a Tcl script to visualize the network in VMD.
    """
    frames = data.get("frames", {})

    if mode == "frame":
        if frame_idx is None or str(frame_idx) not in frames:
            # If not specified, pick the first available frame
            frame_idx = list(frames.keys())[0] if frames else None

        if frame_idx is None:
            logger.error("No frames available for visualization.")
            return

        paths = frames[str(frame_idx)]

        with open(output_file, 'w') as f:
            f.write("graphics top delete all\n")
            f.write("materials change opacity Ghost 0.7\n")
            f.write("graphics top material Ghost\n")
            f.write("graphics top color cyan\n\n")

            for path in paths:
                coords = path["coords"]
                if len(coords) < 2:
                    continue

                for i in range(len(coords) - 1):
                    p1 = coords[i]
                    p2 = coords[i+1]
                    f.write(f"graphics top cylinder {{{p1[0]:.3f} {p1[1]:.3f} {p1[2]:.3f}}} {{{p2[0]:.3f} {p2[1]:.3f} {p2[2]:.3f}}} radius 0.2\n")
                    f.write(f"graphics top sphere {{{p1[0]:.3f} {p1[1]:.3f} {p1[2]:.3f}}} radius 0.3\n")

                # Last node sphere
                last_p = coords[-1]
                f.write(f"graphics top sphere {{{last_p[0]:.3f} {last_p[1]:.3f} {last_p[2]:.3f}}} radius 0.3\n")

        logger.info(f"VMD script written to {output_file} for frame {frame_idx}")

    elif mode == "density":
        # Draw all paths across all frames with low opacity to show density
        with open(output_file, 'w') as f:
            f.write("graphics top delete all\n")
            f.write("materials change opacity Ghost 0.1\n")
            f.write("graphics top material Ghost\n")
            f.write("graphics top color cyan\n\n")

            for f_idx, paths in frames.items():
                for path in paths:
                    coords = path["coords"]
                    for i in range(len(coords) - 1):
                        p1 = coords[i]
                        p2 = coords[i+1]
                        f.write(f"graphics top cylinder {{{p1[0]:.3f} {p1[1]:.3f} {p1[2]:.3f}}} {{{p2[0]:.3f} {p2[1]:.3f} {p2[2]:.3f}}} radius 0.1\n")

        logger.info(f"VMD density script written to {output_file}")


def export_pymol_script(data, output_file="draw_pathways.py", mode="frame", frame_idx=None):
    """
    Generates a Python script to visualize the network in PyMOL as CGO (Compiled Graphics Objects).
    """
    frames = data.get("frames", {})

    with open(output_file, 'w') as f:
        f.write("from pymol import cmd\n")
        f.write("from pymol.cgo import *\n\n")
        f.write("obj = []\n")

        def write_cylinder(p1, p2, radius=0.2, color=[0.0, 1.0, 1.0], alpha=1.0):
            f.write(f"obj.extend([CYLINDER, {p1[0]:.3f}, {p1[1]:.3f}, {p1[2]:.3f}, {p2[0]:.3f}, {p2[1]:.3f}, {p2[2]:.3f}, {radius:.3f}, ")
            f.write(f"{color[0]:.1f}, {color[1]:.1f}, {color[2]:.1f}, {color[0]:.1f}, {color[1]:.1f}, {color[2]:.1f}])\n")

        def write_sphere(p, radius=0.3, color=[0.0, 1.0, 1.0]):
            f.write(f"obj.extend([COLOR, {color[0]:.1f}, {color[1]:.1f}, {color[2]:.1f}])\n")
            f.write(f"obj.extend([SPHERE, {p[0]:.3f}, {p[1]:.3f}, {p[2]:.3f}, {radius:.3f}])\n")

        if mode == "frame":
            if frame_idx is None or str(frame_idx) not in frames:
                frame_idx = list(frames.keys())[0] if frames else None

            if frame_idx is None:
                logger.error("No frames available for visualization.")
                return

            paths = frames[str(frame_idx)]
            for path in paths:
                coords = path["coords"]
                if len(coords) < 2:
                    continue
                for i in range(len(coords) - 1):
                    write_cylinder(coords[i], coords[i+1], radius=0.2)
                    write_sphere(coords[i], radius=0.3)
                write_sphere(coords[-1], radius=0.3)

            logger.info(f"PyMOL script written to {output_file} for frame {frame_idx}")

        elif mode == "density":
            # For PyMOL, transparency in CGO can be tricky, we'll draw thin lines/cylinders
            for f_idx, paths in frames.items():
                for path in paths:
                    coords = path["coords"]
                    for i in range(len(coords) - 1):
                        write_cylinder(coords[i], coords[i+1], radius=0.05, color=[0.0, 0.8, 0.8])
            logger.info(f"PyMOL density script written to {output_file}")

        f.write("\ncmd.load_cgo(obj, 'water_network')\n")

def run_visualization(data_file, format="vmd", mode="density", frame_idx=None, output_file=None):
    if not os.path.exists(data_file):
        logger.error(f"Data file {data_file} not found.")
        return

    with open(data_file, 'r') as f:
        data = json.load(f)

    if format == "vmd":
        out = output_file or ("vmd_density.tcl" if mode == "density" else "vmd_frame.tcl")
        export_vmd_script(data, output_file=out, mode=mode, frame_idx=frame_idx)
    elif format == "pymol":
        out = output_file or ("pymol_density.py" if mode == "density" else "pymol_frame.py")
        export_pymol_script(data, output_file=out, mode=mode, frame_idx=frame_idx)
    else:
        logger.error(f"Unknown format: {format}")
