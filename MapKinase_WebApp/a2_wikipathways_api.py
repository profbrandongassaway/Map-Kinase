import os
import re
import requests
import xml.etree.ElementTree as ET
from PIL import Image
import io
from pathlib import Path
from MapKinase_WebApp.a1_base_api import BasePathwayAPI
from pywikipathways import get_pathway, get_pathway_info
from MapKinase_WebApp.d3_entrez_to_uniprot import entrez_to_uniprot, ensembl_to_uniprot

def _normalize_species_folder(species: str) -> str:
    cleaned = re.sub(r"[^a-z0-9]+", "_", (species or "").strip().lower())
    return cleaned.strip("_") or "unknown"


def _extract_species_from_gpml(gpml_content: str) -> str:
    try:
        root = ET.fromstring(gpml_content)
        for key in ("Organism", "organism", "Species", "species"):
            value = root.attrib.get(key)
            if value:
                return value.strip()
    except Exception:
        return ""
    return ""


class WikiPathwaysAPI(BasePathwayAPI):
    def download_pathway_data(self, pathway_id, species_hint=None):
        base_dir = Path(__file__).resolve().parent.parent / "stored_pathways" / "wikipathways"
        if species_hint:
            species_folder = _normalize_species_folder(str(species_hint))
            cached = base_dir / species_folder / f"{pathway_id}.gpml"
            if cached.exists():
                print(f"Using cached file: {cached}")
                return str(cached)
        if base_dir.exists():
            for cached in base_dir.rglob(f"{pathway_id}.gpml"):
                print(f"Using cached file: {cached}")
                return str(cached)
        candidates = [
            os.path.join(os.getcwd(), f"{pathway_id}.gpml"),
            os.path.join(os.path.dirname(os.getcwd()), f"{pathway_id}.gpml"),
        ]
        fallback_local = None
        for local_file in candidates:
            if os.path.exists(local_file):
                print(f"Using local file: {local_file}")
                fallback_local = local_file
                return local_file
        try:
            print(f"Fetching pathway {pathway_id} using pywikipathways")
            gpml_content = get_pathway(pathway_id)
            species_label = _extract_species_from_gpml(gpml_content)
            if not species_label:
                species_label = str(species_hint or "").strip()
            if not species_label:
                try:
                    info = get_pathway_info(pathway_id)
                    if isinstance(info, dict):
                        species_label = str(info.get("species") or info.get("organism") or "").strip()
                except Exception:
                    species_label = ""
            species_folder = _normalize_species_folder(species_label)
            target_dir = base_dir / species_folder
            target_dir.mkdir(parents=True, exist_ok=True)
            file_path = target_dir / f"{pathway_id}.gpml"
            with file_path.open("w", encoding="utf-8") as f:
                f.write(gpml_content)
            print(f"Successfully downloaded {file_path}")
            return str(file_path)
        except Exception as e:
            print(f"Error downloading {pathway_id} with pywikipathways: {e}")
            try:
                print(f"Fetching pathway info for debugging: {pathway_id}")
                info = get_pathway_info(pathway_id)
                print(f"Pathway info: {info}")
            except Exception as info_e:
                print(f"Failed to get pathway info: {info_e}")
            if fallback_local and os.path.exists(fallback_local):
                print(f"Falling back to local file: {fallback_local}")
                return fallback_local
            raise Exception(f"Failed to download pathway {pathway_id} from WikiPathways")

    def download_pathway_image(self, pathway_id):
        try:
            url = f"https://www.wikipathways.org/wpi/wpi.php?action=downloadFile&type=png&pwTitle=Pathway:{pathway_id}"
            response = requests.get(url)
            response.raise_for_status()
            return Image.open(io.BytesIO(response.content))
        except requests.exceptions.RequestException as e:
            print(f"Error downloading image for {pathway_id}: {e}")
            return None

    def parse_pathway(self, file_path):
        try:
            tree = ET.parse(file_path)
            root = tree.getroot()
            ns = {'gpml': 'http://pathvisio.org/GPML/2013a'}
            print("Top-level tags in GPML:", [elem.tag for elem in root])

            entries = []
            groups = []
            arrows = []
            entrez_ids: set[str] = set()
            ensembl_ids: set[str] = set()

            def _normalize_color(val: str, default: str = "#000000") -> str:
                if not val:
                    return default
                v = val.strip()
                if not v.startswith("#") and len(v) in {6, 8}:
                    v = f"#{v}"
                return v

            # Parse DataNode elements
            datanodes = root.findall(f".//gpml:DataNode", namespaces=ns)
            print(f"Found {len(datanodes)} DataNode elements")
            for datanode in datanodes:
                graph_id = datanode.get("GraphId", "Unknown")
                print(f"Processing DataNode with GraphId: {graph_id}")

                graphics = datanode.find(f"gpml:Graphics", namespaces=ns)
                if graphics is None:
                    print(f"Warning: DataNode {graph_id} is missing Graphics element, skipping")
                    continue

                center_x = graphics.get("CenterX", 0)
                center_y = graphics.get("CenterY", 0)
                width = graphics.get("Width", 0)
                height = graphics.get("Height", 0)

                if not all([center_x, center_y, width, height]):
                    print(f"Warning: DataNode {graph_id} is missing required Graphics attributes, skipping")
                    continue

                node_type = datanode.get("Type", "GeneProduct").lower()
                if node_type in ["geneproduct", "protein", "molecule"]:
                    entry_type = 'prot_box'
                elif node_type == "pathway":
                    entry_type = 'map'  # treat pathway nodes as text boxes with outline
                elif node_type == "metabolite":
                    entry_type = 'compound'
                else:
                    entry_type = node_type

                text_label = datanode.get("TextLabel", "").strip()
                graphics_type = graphics.get("ShapeType") or graphics.get("Shape") or ""
                fg_color_norm = _normalize_color(graphics.get("Color", "#000000"), "#000000")
                entry = {
                    "id": graph_id,
                    "name": text_label,
                    "type": entry_type,
                    "backup_label": text_label,
                    "x": float(center_x) - float(width) / 2,
                    "y": float(center_y) - float(height) / 2,
                    "width": float(width),
                    "height": float(height),
                    "first_name": text_label.split(",")[0].strip() if text_label else "",
                    "fgcolor": fg_color_norm,
                    "bgcolor": _normalize_color(graphics.get("FillColor", "#FFFFFF"), "#FFFFFF"),
                    "group_ref": datanode.get("GroupRef"),
                    "graphics_type": graphics_type,
                    "border_color": fg_color_norm if entry_type == 'map' else None,
                    "border_width": 1 if entry_type == 'map' else None,
                    "xref": {}
                }
                # Apply manual offset for metabolites to better align downstream
                if entry_type == 'compound':
                    entry["x"] -= 15
                    entry["y"] -= 10

                # Add Xref information
                xref = datanode.find(f"gpml:Xref", namespaces=ns)
                if xref is not None:
                    db = xref.get("Database", "")
                    xid = xref.get("ID", "")
                    if db.lower() == "entrez gene" and xid:
                        entrez_ids.add(xid)
                    if db.lower() == "ensembl" and xid:
                        ensembl_ids.add(xid)
                    entry["xref"] = {
                        "Database": db,
                        "ID": xid
                    }

                entries.append(entry)
                print(f"Parsed DataNode {graph_id}: {entry}")

            # Parse Group elements
            group_elements = root.findall(f".//gpml:Group", namespaces=ns)
            print(f"Found {len(group_elements)} Group elements")
            for group in group_elements:
                group_id = group.get("GraphId", "")
                if not group_id:
                    print("Warning: Group element missing GraphId, skipping")
                    continue

                graphics = group.find(f"gpml:Graphics", namespaces=ns)
                if graphics is None:
                    print(f"Warning: Group {group_id} is missing Graphics element, skipping")
                    continue

                group_data = {
                    "id": group_id,
                    "type": "group",
                    "x": float(graphics.get("CenterX", 0)) - float(graphics.get("Width", 0)) / 2,
                    "y": float(graphics.get("CenterY", 0)) - float(graphics.get("Height", 0)) / 2,
                    "width": float(graphics.get("Width", 0)),
                    "height": float(graphics.get("Height", 0)),
                    "fgcolor": graphics.get("Color", "#000000"),
                    "bgcolor": graphics.get("FillColor", "#FFFFFF"),
                    "components": []
                }
                groups.append(group_data)
                print(f"Added group {group_id}")

            # Link DataNodes to Groups via GroupRef
            for entry in entries:
                group_ref = entry.get("group_ref")
                if group_ref:
                    for group in groups:
                        if group["id"] == group_ref:
                            group["components"].append(entry["id"])
                            print(f"Linked DataNode {entry['id']} to group {group_ref}")

            # Parse Interaction elements for arrows
            interactions = root.findall(f".//gpml:Interaction", namespaces=ns)
            print(f"Found {len(interactions)} Interaction elements")
            # Precompute bounding boxes of existing entries to help with elbow routing
            rects = []
            for e in entries:
                try:
                    rects.append({
                        "x1": float(e.get("x", 0)),
                        "y1": float(e.get("y", 0)),
                        "x2": float(e.get("x", 0)) + float(e.get("width", 0)),
                        "y2": float(e.get("y", 0)) + float(e.get("height", 0)),
                    })
                except Exception:
                    continue

            def _segment_intersects_rect(p1, p2, rect):
                """Simple AABB vs segment intersection check."""
                x1, y1 = p1
                x2, y2 = p2
                rx1, ry1, rx2, ry2 = rect["x1"], rect["y1"], rect["x2"], rect["y2"]
                # quick reject if both points on one side
                if max(x1, x2) < rx1 or min(x1, x2) > rx2 or max(y1, y2) < ry1 or min(y1, y2) > ry2:
                    return False
                # Check if either endpoint inside
                if rx1 <= x1 <= rx2 and ry1 <= y1 <= ry2:
                    return True
                if rx1 <= x2 <= rx2 and ry1 <= y2 <= ry2:
                    return True
                # Check intersection with rect edges
                def _intersects(p1, p2, p3, p4):
                    def ccw(a, b, c):
                        return (c[1]-a[1])*(b[0]-a[0]) > (b[1]-a[1])*(c[0]-a[0])
                    return (ccw(p1, p3, p4) != ccw(p2, p3, p4)) and (ccw(p1, p2, p3) != ccw(p1, p2, p4))
                corners = [(rx1, ry1), (rx2, ry1), (rx2, ry2), (rx1, ry2)]
                edges = [ (corners[i], corners[(i+1)%4]) for i in range(4) ]
                for e1, e2 in edges:
                    if _intersects(p1, p2, e1, e2):
                        return True
                return False

            for interaction in interactions:
                graphics = interaction.find(f"gpml:Graphics", namespaces=ns)
                if graphics is not None:
                    points = graphics.findall(f"gpml:Point", namespaces=ns)
                    if len(points) >= 2:
                        connector_type = (graphics.get("ConnectorType", "") if graphics is not None else "").lower()
                        line_style_attr = (graphics.get("LineStyle", "") if graphics is not None else "").lower()

                        def _build_arrow(p_start, p_end, force_line=None):
                            try:
                                x1 = float(p_start.get("X", "0"))
                                y1 = float(p_start.get("Y", "0"))
                                x2 = float(p_end.get("X", "0"))
                                y2 = float(p_end.get("Y", "0"))
                            except Exception:
                                x1 = y1 = x2 = y2 = None
                            arrow_head = p_end.get("ArrowHead", "") or p_start.get("ArrowHead", "")
                            line_style = force_line or "arrow"
                            if not force_line and arrow_head.lower() == "tbar":
                                line_style = "inhibition"
                            arrow = {
                                "entry1": p_start.get("GraphRef"),
                                "entry2": p_end.get("GraphRef"),
                                "type": interaction.get("Type", "interaction"),
                                "line": line_style,
                                "x1": x1,
                                "y1": y1,
                                "x2": x2,
                                "y2": y2,
                                "connector_type": connector_type,
                                "control_points": [],
                            }
                            if line_style_attr == "broken":
                                if arrow["line"] == "arrow":
                                    arrow["line"] = "dashed_arrow"
                                elif arrow["line"] == "inhibition":
                                    arrow["line"] = "dashed_inhibition"
                                else:
                                    arrow["line"] = "dashed_line"
                            return arrow

                        if connector_type == "elbow":
                            if len(points) >= 3:
                                # Use provided three points: first->middle, middle->end
                                arrow_on_start = (points[0].get("ArrowHead", "") or "").lower() in {"tbar", "arrow"}
                                if arrow_on_start:
                                    arrow_seg = _build_arrow(points[0], points[1])
                                    line_seg = _build_arrow(points[1], points[2], force_line="line")
                                else:
                                    line_seg = _build_arrow(points[0], points[1], force_line="line")
                                    arrow_seg = _build_arrow(points[1], points[2])
                                for seg in (line_seg, arrow_seg):
                                    seg["connector_type"] = "elbow"
                                    arrows.append(seg)
                                    print(f"Added arrow (elbow): {seg}")
                            else:
                                # Only two points: fabricate an elbow bend that minimally intersects existing boxes
                                p_start, p_end = points[0], points[1]
                                try:
                                    x1 = float(p_start.get("X", "0")); y1 = float(p_start.get("Y", "0"))
                                    x2 = float(p_end.get("X", "0")); y2 = float(p_end.get("Y", "0"))
                                except Exception:
                                    continue
                                bend_candidates = [ (x2, y1), (x1, y2) ]  # |_ or _|
                                def overlap_score(cp):
                                    score = 0
                                    for rect in rects:
                                        if _segment_intersects_rect((x1,y1), cp, rect):
                                            score += 1
                                        if _segment_intersects_rect(cp, (x2,y2), rect):
                                            score += 1
                                    return score
                                scores = [overlap_score(cp) for cp in bend_candidates]
                                best_cp = bend_candidates[scores.index(min(scores))]
                                # First segment: line; second: arrow
                                mid_point = ET.Element("Point", {"X": str(best_cp[0]), "Y": str(best_cp[1])})
                                line_seg = _build_arrow(p_start, mid_point, force_line="line")
                                arrow_seg = _build_arrow(mid_point, p_end)
                                for seg in (line_seg, arrow_seg):
                                    seg["connector_type"] = "elbow"
                                    arrows.append(seg)
                                    print(f"Added arrow (elbow fabricated): {seg}")
                        else:
                            start_point = points[0]
                            end_point = points[-1]
                            control_points = []
                            # Use middle point as control for curved connectors
                            if connector_type == "curved" and len(points) >= 3:
                                mid = points[1]
                                try:
                                    cx = float(mid.get("X", "0"))
                                    cy = float(mid.get("Y", "0"))
                                    control_points.append({"x": cx, "y": cy})
                                except Exception:
                                    control_points = []
                            arrow = _build_arrow(start_point, end_point)
                            arrow["control_points"] = control_points
                            arrows.append(arrow)
                            print(f"Added arrow: {arrow}")

            def _label_to_entry(label_elem):
                text_label = (label_elem.get("TextLabel") or "").strip()
                graph_id = label_elem.get("GraphId", "") or text_label or ""
                g = label_elem.find(f"gpml:Graphics", namespaces=ns) or label_elem.find(".//{*}Graphics")
                if g is None:
                    return None
                cx = g.get("CenterX", 0)
                cy = g.get("CenterY", 0)
                width = g.get("Width", 0)
                height = g.get("Height", 0)
                if not all([cx, cy, width, height]):
                    return None
                try:
                    font_size_raw = g.get("FontSize")
                    font_size = float(font_size_raw) if font_size_raw is not None else None
                    font_weight = g.get("FontWeight", "").lower()
                    font_style = g.get("FontStyle", "").lower()
                    h_align = (g.get("Align") or "").lower()  # left/center/right
                    v_align = (g.get("Valign") or g.get("VAlign") or "").lower()  # top/middle/bottom
                    text_style = {}
                    if font_size:
                        text_style["fontSize"] = font_size
                    if font_weight == "bold":
                        text_style["bold"] = True
                    if font_style == "italic":
                        text_style["italic"] = True
                    if h_align in {"left", "center", "right"}:
                        text_style["align"] = h_align
                    if v_align in {"top", "middle", "bottom", "center"}:
                        text_style["vertical"] = "center" if v_align == "middle" else v_align

                    return {
                        "id": graph_id,
                        "name": text_label,
                        "first_name": text_label.split(",")[0].strip() if text_label else "",
                        "type": "map",  # handled as text box downstream
                        "graphics_type": g.get("Shape", "label"),
                        "x": float(cx) - float(width) / 2,
                        "y": float(cy) - float(height) / 2,
                        "width": float(width),
                        "height": float(height),
                        "fgcolor": _normalize_color(g.get("Color", "#000000"), "#000000"),
                        "bgcolor": _normalize_color(g.get("FillColor", "#FFFFFF"), "#FFFFFF"),
                        "text_style": text_style,
                    }
                except Exception:
                    return None

            # Parse Label elements into text entries (namespace-aware + namespace-agnostic fallback)
            label_elements = list(root.findall(f".//gpml:Label", namespaces=ns))
            if not label_elements:
                label_elements = list(root.findall(".//{*}Label"))
            print(f"Found {len(label_elements)} Label elements")
            for label in label_elements:
                entry = _label_to_entry(label)
                if entry:
                    entries.append(entry)

            # Parse Shape elements into shape/text-box entries
            shape_elements = list(root.findall(f".//gpml:Shape", namespaces=ns))
            if not shape_elements:
                shape_elements = list(root.findall(".//{*}Shape"))
            print(f"Found {len(shape_elements)} Shape elements")
            for shape in shape_elements:
                graph_id = shape.get("GraphId", "") or ""
                text_label = (shape.get("TextLabel") or "").strip()
                g = shape.find(f"gpml:Graphics", namespaces=ns) or shape.find(".//{*}Graphics")
                if g is None:
                    continue
                cx = g.get("CenterX", 0)
                cy = g.get("CenterY", 0)
                width = g.get("Width", 0)
                height = g.get("Height", 0)
                if not all([cx, cy, width, height]):
                    continue
                shape_type_raw = g.get("ShapeType") or g.get("Shape") or ""
                shape_type = shape_type_raw.lower()
                if shape_type == "roundedrectangle":
                    shape_type = "rounded"
                elif shape_type == "circle":
                    shape_type = "circle"
                elif shape_type == "rectangle":
                    shape_type = "rect"
                fgcolor = _normalize_color(g.get("Color", "#000000"), "#000000")
                bgcolor = _normalize_color(g.get("FillColor", "#FFFFFF"), "#FFFFFF")
                font_size_raw = g.get("FontSize")
                font_size = float(font_size_raw) if font_size_raw is not None else None
                font_weight = (g.get("FontWeight") or "").lower()
                font_style = (g.get("FontStyle") or "").lower()
                h_align = (g.get("Align") or "").lower()  # left/center/right
                v_align = (g.get("Valign") or g.get("VAlign") or "").lower()  # top/middle/bottom
                text_style = {}
                if font_size:
                    text_style["fontSize"] = font_size
                if font_weight == "bold":
                    text_style["bold"] = True
                if font_style == "italic":
                    text_style["italic"] = True
                if h_align in {"left", "center", "right"}:
                    text_style["align"] = h_align
                if v_align in {"top", "middle", "bottom", "center"}:
                    text_style["vertical"] = "center" if v_align == "middle" else v_align
                border_width = None
                try:
                    lw = g.get("LineThickness")
                    if lw is not None:
                        border_width = float(lw)
                except Exception:
                    border_width = None
                entry = {
                    "id": graph_id or text_label or "",
                    "name": text_label,
                    "first_name": text_label.split(",")[0].strip() if text_label else "",
                    "type": "map",  # render as text-box with shape
                    "graphics_type": shape_type_raw,
                    "shape_type": shape_type,
                    "x": float(cx) - float(width) / 2,
                    "y": float(cy) - float(height) / 2,
                    "width": float(width),
                    "height": float(height),
                    "fgcolor": fgcolor,
                    "bgcolor": bgcolor,
                    "border_color": fgcolor,
                    "border_width": border_width if border_width is not None else 1,
                    "text_style": text_style,
                    "link": "",
                }
                entries.append(entry)
                # If shape has ANY Attribute with Value="Double", add a duplicate shape slightly larger (width/height +10)
                attrs = list(shape.findall(f"gpml:Attribute", namespaces=ns)) or list(shape.findall(".//{*}Attribute"))
                is_cell_bg = any((attr.get("Value") or "").lower() == "cell" and "cellularcomponentproperty" in (attr.get("Key") or "").lower() for attr in attrs)
                has_double = any((attr.get("Value") or "").lower() == "double" for attr in attrs)
                if has_double:
                    dup = dict(entry)
                    dup["id"] = f"{entry['id']}_double" if entry.get("id") else f"{graph_id}_double"
                    dup["width"] = entry["width"] + 10
                    dup["height"] = entry["height"] + 10
                    dup["x"] = entry["x"] - 5
                    dup["y"] = entry["y"] - 5
                    dup["is_background"] = is_cell_bg
                    entries.append(dup)
                entry["is_background"] = is_cell_bg

            # Fallback scan if nothing was added (handles unexpected namespaces or casing)
            if not any(e.get("type") == "map" or e.get("graphics_type") in {"label", "roundrectangle"} for e in entries):
                for elem in root.iter():
                    if not (isinstance(elem.tag, str) and elem.tag.lower().endswith("label")):
                        continue
                    entry = _label_to_entry(elem)
                    if entry:
                        entries.append(entry)

            # Map Entrez/Ensembl IDs to UniProt IDs when possible (batch using helper)
            if entrez_ids or ensembl_ids:
                id_to_uniprot: dict[str, str] = {}
                id_db_map: dict[str, str] = {}  # track which DB each id belongs to
                try:
                    for eid in sorted(entrez_ids):
                        try:
                            uni = entrez_to_uniprot(eid)
                            if uni:
                                id_to_uniprot[eid] = uni
                                id_db_map[eid] = "entrez gene"
                        except Exception as one_exc:
                            print(f"Entrez->UniProt lookup failed for {eid}: {one_exc}")
                    for ens in sorted(ensembl_ids):
                        try:
                            uni = ensembl_to_uniprot(ens)
                            if uni:
                                id_to_uniprot[ens] = uni
                                id_db_map[ens] = "ensembl"
                        except Exception as one_exc:
                            print(f"Ensembl->UniProt lookup failed for {ens}: {one_exc}")
                    if id_to_uniprot:
                        for entry in entries:
                            xref = entry.get("xref") or {}
                            db = str(xref.get("Database", "")).lower()
                            xid = str(xref.get("ID", "")).strip()
                            if xid and xid in id_to_uniprot and ((db == id_db_map.get(xid)) or (db in {"entrez gene", "ensembl"})):
                                uni = id_to_uniprot[xid]
                                entry.setdefault("xref", {})["UniProt"] = uni
                                entry["uniprot_id"] = uni
                                # Normalize the primary identifier to UniProt to align with KEGG-style processing
                                entry.setdefault("label", entry.get("name", ""))
                                entry.setdefault("backup_label", entry.get("name", ""))
                                entry["name"] = uni
                                # Keep display-friendly first_name for fallback labels
                                entry["first_name"] = entry.get("backup_label") or entry.get("label") or entry.get("first_name") or uni
                except Exception as map_exc:
                    print(f"Warning: failed to map IDs to UniProt: {map_exc}")

            print(f"Parsed {len(entries)} entries, {len(groups)} groups, and {len(arrows)} arrows")
            return entries, groups, arrows

        except ET.ParseError as e:
            print(f"XML Parse Error for pathway file {file_path}: {e}")
            return [], [], []
        except Exception as e:
            print(f"Error parsing pathway file {file_path}: {e}")
            return [], [], []
