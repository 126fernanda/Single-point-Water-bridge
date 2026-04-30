import json
import logging
import os

logger = logging.getLogger(__name__)

def read_cluster_json(data_file):
    with open(data_file, 'r') as f:
        return json.load(f)

def read_jsonl(data_file):
    frames = {}
    with open(data_file, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            obj = json.loads(line)
            if obj.get('type') == 'frame':
                f_idx = str(obj['frame_idx'])
                if f_idx not in frames:
                    frames[f_idx] = []
                frames[f_idx].extend(obj['paths'])
    return frames

def export_vmd_script(data_file, output_file="draw_pathways.tcl", mode="frame", frame_idx=None):
    """
    Generates a Tcl script to visualize the network in VMD using dynamic selections.
    """
    if mode in ["frame", "density"]:
        frames = read_jsonl(data_file)
        
    if mode == "frame":
        if frame_idx is None or str(frame_idx) not in frames:
            frame_idx = list(frames.keys())[0] if frames else None

        if frame_idx is None:
            logger.error("No frames available for visualization.")
            return

        paths = frames[str(frame_idx)]

        # Collect all atom indices
        all_indices = set()
        for path in paths:
            all_indices.update(path["nodes"])

        indices_str = " ".join(str(i) for i in all_indices)

        if not indices_str:
            logger.warning(f"No paths found for frame {frame_idx}. Skipping VMD script generation to prevent Tcl crash.")
            return

        with open(output_file, 'w') as f:
            f.write("mol color Name\n")
            f.write("mol representation Licorice 0.3 12 12\n")
            f.write(f"mol selection \"index {indices_str}\"\n")
            f.write("mol addrep top\n")
            f.write("mol delrep 0 top\n")

            # Optional: Move to the frame
            f.write(f"animate goto {frame_idx}\n")

        logger.info(f"VMD script written to {output_file} for frame {frame_idx}")

    elif mode == "density":
        with open(output_file, 'w') as f:
            f.write("graphics top delete all\n")
            f.write("material change opacity Ghost 0.1\n")
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
        
    elif mode == "cluster":
        clusters = read_cluster_json(data_file)
        with open(output_file, 'w') as f:
            f.write("graphics top delete all\n")
            f.write("material change opacity Opaque 1.0\n")
            f.write("graphics top material Opaque\n")
            f.write("graphics top color orange\n\n")

            for cluster in clusters:
                coords = cluster.get("medoid_coords", [])
                if len(coords) < 2:
                    continue
                for i in range(len(coords) - 1):
                    p1 = coords[i]
                    p2 = coords[i+1]
                    f.write(f"graphics top cylinder {{{p1[0]:.3f} {p1[1]:.3f} {p1[2]:.3f}}} {{{p2[0]:.3f} {p2[1]:.3f} {p2[2]:.3f}}} radius 0.2\n")
                    f.write(f"graphics top sphere {{{p1[0]:.3f} {p1[1]:.3f} {p1[2]:.3f}}} radius 0.3\n")
                
                # write the last sphere
                p_last = coords[-1]
                f.write(f"graphics top sphere {{{p_last[0]:.3f} {p_last[1]:.3f} {p_last[2]:.3f}}} radius 0.3\n")

        logger.info(f"VMD cluster script written to {output_file}")


def export_pymol_script(data_file, output_file="draw_pathways.py", mode="frame", frame_idx=None):
    """
    Generates a Python script to visualize the network in PyMOL as CGO.
    """
    if mode in ["frame", "density"]:
        frames = read_jsonl(data_file)

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
                f.write("\ncmd.load_cgo(obj, 'water_network')\n")
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
            for f_idx, paths in frames.items():
                for path in paths:
                    coords = path["coords"]
                    for i in range(len(coords) - 1):
                        write_cylinder(coords[i], coords[i+1], radius=0.05, color=[0.0, 0.8, 0.8])
            logger.info(f"PyMOL density script written to {output_file}")
            
        elif mode == "cluster":
            clusters = read_cluster_json(data_file)
            for cluster in clusters:
                coords = cluster.get("medoid_coords", [])
                if len(coords) < 2:
                    continue
                for i in range(len(coords) - 1):
                    write_cylinder(coords[i], coords[i+1], radius=0.2, color=[1.0, 0.5, 0.0], alpha=1.0)
                    write_sphere(coords[i], radius=0.3, color=[1.0, 0.5, 0.0])
                write_sphere(coords[-1], radius=0.3, color=[1.0, 0.5, 0.0])
            logger.info(f"PyMOL cluster script written to {output_file}")

        f.write("\ncmd.load_cgo(obj, 'water_network')\n")
        if mode == "density":
            f.write("cmd.set('cgo_transparency', 0.8, 'water_network')\n")


def export_chimera_script(data_file, output_file="draw_pathways.py", mode="frame", frame_idx=None):
    """
    Generates a Python script to visualize the network in UCSF Chimera.
    """
    if mode in ["frame", "density"]:
        frames = read_jsonl(data_file)

    if mode == "frame":
        if frame_idx is None or str(frame_idx) not in frames:
            frame_idx = list(frames.keys())[0] if frames else None

        if frame_idx is None:
            logger.error("No frames available for visualization.")
            return

        paths = frames[str(frame_idx)]

        all_ids = set()
        for path in paths:
            # Use atom_ids if available (from analysis.py upgrade), otherwise fallback to indices+1
            if "atom_ids" in path:
                all_ids.update(path["atom_ids"])
            else:
                all_ids.update(i+1 for i in path["nodes"])

        with open(output_file, 'w') as f:
            f.write("import chimera\n")
            f.write("from chimera import runCommand\n")
            sel_str = ", ".join(str(i) for i in sorted(all_ids))
            f.write(f"serial_set = set([{sel_str}])\n")
            f.write("atoms = [a for m in chimera.openModels.list() for a in m.atoms if a.serialNumber in serial_set]\n")
            f.write("chimera.selection.setCurrent(atoms)\n")
            f.write("runCommand('show sel')\n")
            f.write("runCommand('repr stick sel')\n")

        logger.info(f"Chimera script written to {output_file} for frame {frame_idx}")

    elif mode == "density":
        bild_file = os.path.abspath(output_file.replace('.py', '.bild'))
        with open(bild_file, 'w') as bild_f:
            bild_f.write('.color cyan\n')
            
            for f_idx, paths in frames.items():
                for path in paths:
                    coords = path["coords"]
                    for i in range(len(coords) - 1):
                        p1 = coords[i]
                        p2 = coords[i+1]
                        bild_f.write(f".cylinder {p1[0]:.3f} {p1[1]:.3f} {p1[2]:.3f} {p2[0]:.3f} {p2[1]:.3f} {p2[2]:.3f} 0.1\n")
        with open(output_file, 'w') as f:
            f.write("from chimera import runCommand, openModels\n")
            f.write("before_models = set(openModels.list())\n")
            f.write(f"runCommand(\"open '{bild_file}'\")\n")
            f.write("new_models = set(openModels.list()) - before_models\n")
            f.write("for m in new_models:\n")
            f.write("    runCommand('transparency 50 models #{}'.format(m.id))\n")
        logger.info(f"Chimera density script written to {output_file} (BILD geometry: {bild_file})")
        
    elif mode == "cluster":
        clusters = read_cluster_json(data_file)
        bild_file = os.path.abspath(output_file.replace('.py', '.bild'))
        with open(bild_file, 'w') as bild_f:
            bild_f.write('.color orange\n')

            for cluster in clusters:
                coords = cluster.get("medoid_coords", [])
                if len(coords) < 2:
                    continue
                for i in range(len(coords) - 1):
                    p1 = coords[i]
                    p2 = coords[i+1]
                    bild_f.write(f".cylinder {p1[0]:.3f} {p1[1]:.3f} {p1[2]:.3f} {p2[0]:.3f} {p2[1]:.3f} {p2[2]:.3f} 0.2\n")
                    bild_f.write(f".sphere {p1[0]:.3f} {p1[1]:.3f} {p1[2]:.3f} 0.3\n")

                # write the last sphere
                p_last = coords[-1]
                bild_f.write(f".sphere {p_last[0]:.3f} {p_last[1]:.3f} {p_last[2]:.3f} 0.3\n")

        with open(output_file, 'w') as f:
            f.write("from chimera import runCommand, openModels\n")
            f.write(f"runCommand(\"open '{bild_file}'\")\n")
        logger.info(f"Chimera cluster script written to {output_file} (BILD geometry: {bild_file})")


def run_visualization(data_file, format="vmd", mode="density", frame_idx=None, output_file="pathways_viz"):
    if not os.path.exists(data_file):
        logger.error(f"Data file {data_file} not found.")
        return

    # Append appropriate extension based on format if not provided by user
    if format == "vmd" and not output_file.endswith(".tcl"):
        out = f"{output_file}.tcl"
    elif format in ["pymol", "chimera"] and not output_file.endswith(".py"):
        out = f"{output_file}.py"
    else:
        out = output_file

    if format == "vmd":
        export_vmd_script(data_file, output_file=out, mode=mode, frame_idx=frame_idx)
    elif format == "pymol":
        export_pymol_script(data_file, output_file=out, mode=mode, frame_idx=frame_idx)
    elif format == "chimera":
        export_chimera_script(data_file, output_file=out, mode=mode, frame_idx=frame_idx)
    else:
        logger.error(f"Unknown format: {format}")
