import json
import os
from pathlib import Path
from typing import Optional
from shiny import App, ui, render, reactive

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR / "output" / "testing_file_001"
DEFAULT_JSON_PATH = Path(
    r"C:\Users\clayt\OneDrive - Brigham Young University\Desktop\MapKinase\MapKinase_WebApp\hsa04010_pathway_data_20250907_210715.json"
)


def _find_latest_json(directory: Path) -> Optional[Path]:
    try:
        if not directory.exists():
            return None
        candidates = sorted(
            directory.glob("*_pathway_data*.json"),
            key=lambda p: p.stat().st_mtime,
            reverse=True,
        )
        return candidates[0] if candidates else None
    except OSError:
        return None


def _resolve_json_path() -> Optional[Path]:
    env_path = os.environ.get("MAPKINASE_JSON")
    if env_path:
        env_candidate = Path(env_path).expanduser()
        if env_candidate.exists():
            return env_candidate
    latest = _find_latest_json(OUTPUT_DIR)
    if latest:
        return latest
    if DEFAULT_JSON_PATH.exists():
        return DEFAULT_JSON_PATH
    return None


_json_path_obj = _resolve_json_path()
if _json_path_obj:
    json_file_path = str(_json_path_obj)
    print(f"Loading pathway JSON from: {json_file_path}")
else:
    json_file_path = ""
    print("Warning: No MAPKINASE_JSON provided and no generated pathway JSON found; viewer will start without data.")

def load_json_data():
    try:
        if not os.path.exists(json_file_path):
            return None
        with open(json_file_path, 'r') as f:
            json_data = json.load(f)
        return json_data
    except:
        return None

def save_json_data(json_data):
    try:
        with open(json_file_path, 'w') as f:
            json.dump(json_data, f, indent=4)
    except:
        pass

def create_pathway_svg(json_data, show_kegg_bg=False):
    if not json_data:
        return ui.div("Error: Could not load JSON data.")
    settings = json_data.get('general_data', {}).get('settings', {})
    protein_data = json_data.get('protein_data', {})
    protbox_data = json_data.get('protbox_data', [])
    groups = json_data.get('groups', [])
    arrows = json_data.get('arrows', [])
    compound_data = json_data.get('compound_data', [])
    text_data = json_data.get('text_data', [])
    show_arrows = bool(settings.get('show_arrows', True))
    show_text_boxes = bool(settings.get('show_text_boxes', True))
    catalog_data = {}
    catalog_info = json_data.get('_global_protein_catalog') or {}
    catalog_path = catalog_info.get('path') or os.environ.get('GLOBAL_PROTEIN_CATALOG_PATH')
    if catalog_path and os.path.exists(catalog_path):
        try:
            with open(catalog_path, 'r', encoding='utf-8') as catalog_fh:
                catalog_payload = json.load(catalog_fh)
            if isinstance(catalog_payload, dict):
                catalog_data = catalog_payload.get('protein_catalog', {}) or {}
        except Exception as catalog_err:
            print(f"Warning: failed to load global protein catalog from {catalog_path}: {catalog_err}")
    if not protbox_data:
        max_x, max_y = 800, 600
    else:
        max_x = max(pb.get('x', 0) + pb.get('width', 0) for pb in protbox_data if pb.get('x') is not None and pb.get('width') is not None) + 50
        max_y = max(pb.get('y', 0) + pb.get('height', 0) for pb in protbox_data if pb.get('y') is not None and pb.get('height') is not None) + 50
        max_x = max(max_x, 800)
        max_y = max(max_y, 600)
    # Embed the full JSON (including any 'kegg_bg_image' and preview settings) for the client
    data_script = f"""<script type="application/json" id="pathway-data">{json.dumps(json_data)}</script>"""
    catalog_script = f"""<script type="application/json" id="global-protein-catalog">{json.dumps(catalog_data)}</script>"""
    svg_js = f"""
        <script src="https://cdnjs.cloudflare.com/ajax/libs/svg.js/3.2.0/svg.min.js"></script>
    <script>
    // Use `var` for globals so re-rendering the HTML (which runs this script again)
    // won't throw a 'Identifier ... has already been declared' error.
    var selectedElement = null, selectedType = null, selectedId = null, selectedProtboxId = null, selectedVisual = null;
    var zoomLevel = 1, minZoom = 0.5, maxZoom = 2;
    var viewBox = {{ x: 0, y: 0, width: {max_x}, height: {max_y} }};
    var isPanning = false, startPanX = null, startPanY = null, startViewBoxX = null, startViewBoxY = null;
    var arrowHandleGroups = arrowHandleGroups || {{}};
    var snapRadius = 3;
    var AXIS_SNAP_TOLERANCE_DEG = typeof AXIS_SNAP_TOLERANCE_DEG === 'number' ? AXIS_SNAP_TOLERANCE_DEG : 2;
    var AXIS_SNAP_ANGLES = Array.isArray(AXIS_SNAP_ANGLES) ? AXIS_SNAP_ANGLES : [0, 90, 180, 270];
    var protboxMap = protboxMap || {{}};
    var attachments = attachments || {{}};
    var attachedBy = attachedBy || {{}};
    var protboxSnapPoints = protboxSnapPoints || {{}};
    var protboxHandleDists = protboxHandleDists || {{}};
    var labelDefaults = labelDefaults || {{'N1': [-17, -10, 'right'],'N2': [-6, -12, 'center'],'N3': [5, -10, 'left'],'S1': [-18, 4, 'right'],'S2': [0, 10, 'center'],'S3': [6, 2, 'left'],'W1': [-18, -6, 'right'],'W2': [-18, 2, 'right'],'E1': [5, -6, 'left'],'E2': [5, 2, 'left']}};
    var anchorMap = anchorMap || {{'left': 'start','center': 'middle','right': 'end'}};
    var mkHistory = typeof mkHistory === 'object' && mkHistory ? mkHistory : null;
    var activeSnapKey = typeof activeSnapKey === 'string' ? activeSnapKey : null;
    var proteinCatalog = (function() {{
        try {{
            var catalogNode = document.getElementById('global-protein-catalog');
            if (catalogNode) {{
                var parsedPayload = JSON.parse(catalogNode.textContent || '{{}}') || {{}};
                if (typeof parsedPayload === 'object' && !Array.isArray(parsedPayload)) {{
                    return parsedPayload;
                }}
            }}
        }} catch (err) {{
            console.warn('m3: failed to parse global protein catalog', err);
        }}
        return {{}};
    }})();
    // Pick the first unoccupied PTM snap point for a protbox
    var ptmSnapRadius = 12;
    var currentSelected = currentSelected || {{}};
    function initializeSvg() {{
            const container = document.querySelector('.svg-container');
            const canvas = document.getElementById('svgCanvas');
            if (!container || !canvas) return;
            while (canvas.firstChild) {{canvas.removeChild(canvas.firstChild);}}
            // Reset client-side caches so switching pathways/modes doesn't keep stale geometry data
            selectedElement = null;
            selectedType = null;
            selectedId = null;
            selectedProtboxId = null;
            selectedVisual = null;
            protboxMap = {{}};
            attachments = {{}};
            arrowHandleGroups = {{}};
            protboxSnapPoints = {{}};
            protboxHandleDists = {{}};
            currentSelected = {{}};
            // If a KEGG background image exists in the JSON and preview is enabled, embed it into the SVG
            try {{
                const showBg = {str(show_kegg_bg).lower()};
                const dataPreview = null;
                if (showBg) {{
                    // data will be parsed below; we'll operate after parsing
                }}
            }} catch (e) {{ console.log('bg pre-check failed', e); }}
            const draw = SVG().addTo('#svgCanvas').size({max_x}, {max_y}).viewbox(viewBox.x, viewBox.y, viewBox.width, viewBox.height);
            const pickFreeSnap = (protboxId, radius = ptmSnapRadius) => {{
              const snaps = protboxSnapPoints[protboxId];
              if (!snaps) return null;
            
              // All PTM shapes already attached to this protbox
              const ptms = draw.find(`[data-type="ptm-shape"][data-protbox-id="${{protboxId}}"]`);
            
              const occupied = new Set();
              for (const key in snaps) {{
                const sp = snaps[key];
                for (let i = 0; i < ptms.length; i++) {{
                  const el = ptms[i];
                  const isCircle = el.type === 'circle';
                  const x = parseFloat(isCircle ? el.cx() : (el.x() || 0));
                  const y = parseFloat(isCircle ? el.cy() : (el.y() || 0));
                  const d = Math.hypot(x - sp.x, y - sp.y);
                  if (d <= radius) {{ occupied.add(key); break; }}
                }}
              }}
            
              // Preferred order (tweak if you like)
              const order = ['N2','N1','N3','S2','S1','S3','E1','E2','W1','W2'];
              for (const key of order) {{
                if (snaps[key] && !occupied.has(key)) return {{ key, x: snaps[key].x, y: snaps[key].y }};
              }}
              return null;
            }};
            const data = JSON.parse(document.getElementById('pathway-data').textContent);
            const preview = data._kegg_preview || {{}};
            const offx = 22.3;
            const offy = 8;
            const fgOffsetX = -offx;
            const fgOffsetY = -offy;
            // Embed KEGG image into the SVG using SVG.js so it shares the same transforms (pan/zoom)
            try {{
                const existingBg = draw.findOne('image[data-kegg-bg="1"]');
                if (data && data.kegg_bg_image && data._kegg_preview && data._kegg_preview.show) {{
                    const opacity = typeof preview.opacity === 'number' ? preview.opacity : 0.9;
                    const size = data.kegg_bg_size || {{}};
                    const fw = Number(size.width);
                    const fh = Number(size.height);
                    const hasFallbackSize = Number.isFinite(fw) && fw > 0 && Number.isFinite(fh) && fh > 0;
                    const fallbackSize = hasFallbackSize ? {{ width: fw, height: fh }} : null;
                    try {{
                        const scale = typeof preview.scale === 'number' ? preview.scale : 1.0;
                        console.log('m3: embedding/updating kegg bg image (canvas)', {max_x}, {max_y}, 'scale', scale, 'offsets', offx, offy, 'opacity', opacity);
                        // Size the image using the provided scale so users can zoom the background separately
                        const targetScale = Number.isFinite(scale) ? scale : 1;
                        let imgEl = existingBg;
                        const applySizing = (element, naturalWidth, naturalHeight) => {{
                            const baseWidth = Number(naturalWidth) || (fallbackSize ? fallbackSize.width : null) || Number(element.width()) || {max_x};
                            const baseHeight = Number(naturalHeight) || (fallbackSize ? fallbackSize.height : null) || Number(element.height()) || {max_y};
                            const width = baseWidth * targetScale;
                            const height = baseHeight * targetScale;
                            const offsetX = 0;
                            const offsetY = 0;
                            console.log('m3: applySizing offsets', offsetX, offsetY, 'width/height', width, height);
                            try {{ element.size(width, height); element.move(offsetX, offsetY).opacity(opacity); }} catch (uerr) {{ console.log('update bg sizing failed', uerr); }}
                            try {{ element.node.style.pointerEvents = 'none'; }} catch (pe) {{ /* ignore */ }}
                            try {{ element.attr({{ preserveAspectRatio: 'xMidYMid meet', 'data-natural-width': baseWidth, 'data-natural-height': baseHeight }}); }} catch (pa) {{ /* ignore */ }}
                            element.back();
                        }};
                        if (imgEl) {{
                            const existingWidth = parseFloat(imgEl.attr('data-natural-width'));
                            const existingHeight = parseFloat(imgEl.attr('data-natural-height'));
                            const naturalWidth = Number.isFinite(existingWidth) ? existingWidth : (fallbackSize ? fallbackSize.width : Number(imgEl.width()));
                            const naturalHeight = Number.isFinite(existingHeight) ? existingHeight : (fallbackSize ? fallbackSize.height : Number(imgEl.height()));
                            applySizing(imgEl, naturalWidth, naturalHeight);
                            console.log('m3: updated existing kegg bg image');
                        }} else {{
                            imgEl = draw.image(data.kegg_bg_image);
                            try {{ imgEl.attr({{ 'data-kegg-bg': '1' }}); }} catch (aerr) {{ /* ignore */ }}
                            imgEl.loaded(function(loader) {{
                                const candidateWidth = Number(loader && loader.width);
                                const candidateHeight = Number(loader && loader.height);
                                const naturalWidth = Number.isFinite(candidateWidth) && candidateWidth > 0 ? candidateWidth : (fallbackSize ? fallbackSize.width : Number(imgEl.width()) || {max_x});
                                const naturalHeight = Number.isFinite(candidateHeight) && candidateHeight > 0 ? candidateHeight : (fallbackSize ? fallbackSize.height : Number(imgEl.height()) || {max_y});
                                applySizing(imgEl, naturalWidth, naturalHeight);
                                console.log('m3: kegg bg image embedded');
                            }});
                        }}
                    }} catch (eimg) {{ console.log('embedding/updating bg image failed', eimg); }}
                }} else {{
                    // preview disabled or no image in data: remove existing background if present
                    if (existingBg) {{
                        try {{ existingBg.remove(); console.log('m3: removed existing kegg bg image because preview is off'); }} catch (rerr) {{ console.log('failed to remove existing bg', rerr); }}
                    }}
                }}
            }} catch (e) {{ console.log('bg embed check failed', e); }}
            // Client-side diagnostic: write protbox count into the debug area and console
            try {{
                const dbgEl = document.getElementById('debug_json');
                const protCount = (data && data.protbox_data) ? data.protbox_data.length : 0;
                if (dbgEl) dbgEl.textContent = 'Client protbox count: ' + protCount;
                console.log('m3 client: protbox count =', protCount);
            }} catch (dbgErr) {{
                console.log('m3 client debug write failed', dbgErr);
            }}
            const {{ 
              protbox_data: protBoxes, 
              protein_data: proteinData, 
              groups, 
              arrows, 
              compound_data: compounds = [], 
              text_data: textBlocks = [],   
              general_data: {{ settings = {{}} }} 
            }} = data;
            protBoxes.forEach(pb => {{
                if (!pb || typeof pb !== 'object') {{
                    return;
                }}
                if (!pb.ptm_overrides || typeof pb.ptm_overrides !== 'object') {{
                    pb.ptm_overrides = {{}};
                }}
            }});
            const showArrows = {str(show_arrows).lower()};
            const showTextBoxes = {str(show_text_boxes).lower()};
            const colorPreview = data._color_preview_override || {{}};
            const searchSource = Object.keys(proteinCatalog).length ? proteinCatalog : (proteinData || {{}});
            const proteinSearchIndex = Object.entries(searchSource).map(([uniprot, protein]) => {{
                const geneSymbol = protein?.gene_symbol || protein?.label || protein?.backup_label || uniprot || 'Unknown';
                return {{
                    uniprot,
                    geneSymbol,
                    searchText: `${{uniprot}}|${{geneSymbol}}`
                }};
            }});
            let protboxCounter = protBoxes.reduce((max, pb) => {{
                const val = Number(pb?.protbox_id);
                if (Number.isFinite(val)) {{
                    return Math.max(max, val);
                }}
                return max;
            }}, 0);
            let compoundAutoCounter = 0;
            let textAutoCounter = 0;
            const registerCompoundAutoId = (value) => {{
                const match = /^compound_auto_(\\d+)$/.exec(value || '');
                if (match) {{
                    const num = Number(match[1]);
                    if (Number.isFinite(num)) {{
                        compoundAutoCounter = Math.max(compoundAutoCounter, num);
                    }}
                }}
            }};
            const registerTextAutoId = (value) => {{
                const match = /^text_auto_(\\d+)$/.exec(value || '');
                if (match) {{
                    const num = Number(match[1]);
                    if (Number.isFinite(num)) {{
                        textAutoCounter = Math.max(textAutoCounter, num);
                    }}
                }}
            }};
            const protboxDefaults = {{
                width: settings?.prot_box_width || 46,
                height: settings?.prot_box_height || 17
            }};
            const prioritizedPtmPositions = ['N1', 'N3', 'S1', 'S3'];
            const normalizeProtboxId = (value) => (value === undefined || value === null ? null : String(value));
            const findProtboxById = (protboxId) => {{
                const targetId = normalizeProtboxId(protboxId);
                if (targetId === null) return null;
                return protBoxes.find(pb => normalizeProtboxId(pb?.protbox_id) === targetId) || null;
            }};
            const ensureProtboxOverride = (protboxId, uniprot) => {{
                const pb = findProtboxById(protboxId);
                if (!pb || !uniprot) return null;
                if (!pb.ptm_overrides || typeof pb.ptm_overrides !== 'object') {{
                    pb.ptm_overrides = {{}};
                }}
                if (!pb.ptm_overrides[uniprot] || typeof pb.ptm_overrides[uniprot] !== 'object') {{
                    pb.ptm_overrides[uniprot] = {{}};
                }}
                return pb.ptm_overrides[uniprot];
            }};
            const getPtmOverride = (protboxOrId, uniprot, ptmKey) => {{
                if (!uniprot || !ptmKey) return null;
                const pb = typeof protboxOrId === 'object' && protboxOrId !== null ? protboxOrId : findProtboxById(protboxOrId);
                if (!pb) return null;
                const map = pb.ptm_overrides && pb.ptm_overrides[uniprot];
                return map && typeof map === 'object' ? map[ptmKey] || null : null;
            }};
            const recordPtmOverride = (protboxId, uniprot, ptmKey, payload = {{}}) => {{
                const overrides = ensureProtboxOverride(protboxId, uniprot);
                if (!overrides || !ptmKey || !payload || typeof payload !== 'object') return;
                overrides[ptmKey] = {{ ...(overrides[ptmKey] || {{}}), ...payload }};
            }};
            const getProtboxesForUniprot = (uniprot) => {{
                if (!uniprot) return [];
                return protBoxes.filter(pb => Array.isArray(pb?.proteins) && pb.proteins.includes(uniprot));
            }};
            const isPrimaryProtbox = (protboxId, uniprot) => {{
                if (!uniprot) return true;
                const matches = getProtboxesForUniprot(uniprot);
                if (!matches.length) return true;
                const firstId = normalizeProtboxId(matches[0]?.protbox_id);
                return firstId === normalizeProtboxId(protboxId);
            }};
            const arrowIdPattern = /^arrow_(\d+)$/;
            const parseArrowIndex = (arrowId) => {{
                if (typeof arrowId !== 'string') return null;
                const match = arrowIdPattern.exec(arrowId);
                if (!match) return null;
                const idx = Number(match[1]);
                return Number.isFinite(idx) ? idx : null;
            }};
            const getArrowById = (arrowId) => {{
                const idx = parseArrowIndex(arrowId);
                if (idx === null) return null;
                return arrows[idx] || null;
            }};
            const arrowIdFromIndex = (index) => (Number.isInteger(index) && index >= 0 ? `arrow_${{index}}` : null);
            const ensureNumber = (value, fallback) => {{
                const num = Number(value);
                return Number.isFinite(num) ? num : fallback;
            }};
            const toFiniteNumber = (value) => {{
                if (value === null || value === undefined) return null;
                const num = Number(value);
                return Number.isFinite(num) ? num : null;
            }};
            const resolveCoordinate = (overrideVal, baseVal) => {{
                const overrideNum = toFiniteNumber(overrideVal);
                if (overrideNum !== null) return overrideNum;
                return toFiniteNumber(baseVal);
            }};
            const parsePtmElementId = (elementId) => {{
                if (!elementId) return null;
                const parts = String(elementId).split('_');
                if (parts.length < 3) return null;
                const suffix = parts.pop();
                const uniprot = parts.shift();
                const ptmKey = parts.join('_');
                if (!uniprot || !ptmKey) return null;
                return {{ uniprot, ptmKey, suffix }};
            }};
            const pickSettingValue = (key, fallback = null) => {{
                if (colorPreview[key] !== undefined && colorPreview[key] !== null) {{
                    return colorPreview[key];
                }}
                if (settings && settings[key] !== undefined && settings[key] !== null) {{
                    return settings[key];
                }}
                return fallback;
            }};
            const defaultGray = 'rgb(128,128,128)';
            const clampChannel = (value) => Math.min(255, Math.max(0, Math.round(value)));
            const toRgbString = (arr, fallback = defaultGray) => {{
                if (!Array.isArray(arr) || arr.length !== 3) return fallback;
                const parsed = arr.map((v) => Number(v));
                if (parsed.some((v) => !Number.isFinite(v))) return fallback;
                const clamped = parsed.map(clampChannel);
                return `rgb(${{clamped.join(',')}})`;
            }};
            const parseFoldChange = (entity, idx = 1) => {{
                if (!entity) return null;
                const key = `fold_change_${{idx}}`;
                const raw = entity[key];
                if (raw === null || raw === undefined) return null;
                const num = Number(raw);
                return Number.isFinite(num) ? num : null;
            }};
            const gradientConfig = (() => {{
                const normalizeColor = (arr, fallback) => {{
                    if (!Array.isArray(arr) || arr.length !== 3) return fallback;
                    const parsed = arr.map((v) => Number(v));
                    return parsed.some((v) => !Number.isFinite(v)) ? fallback : parsed;
                }};
                const toNumberOr = (value, fallback) => {{
                    const num = Number(value);
                    return Number.isFinite(num) ? num : fallback;
                }};
                const negLimit = toNumberOr(pickSettingValue('max_negative', -2), -2) || -2;
                const posLimit = toNumberOr(pickSettingValue('max_positive', 2), 2) || 2;
                return {{
                    negColor: normalizeColor(pickSettingValue('negative_color', [255, 0, 0]), [255, 0, 0]),
                    posColor: normalizeColor(pickSettingValue('positive_color', [0, 0, 255]), [0, 0, 255]),
                    maxNeg: negLimit === 0 ? -2 : negLimit,
                    maxPos: posLimit === 0 ? 2 : posLimit,
                }};
            }})();
            const gradientColorFromFold = (foldValue) => {{
                if (!Number.isFinite(foldValue)) return defaultGray;
                const white = [255, 255, 255];
                if (foldValue < 0) {{
                    const denom = Math.abs(gradientConfig.maxNeg) > 0 ? Math.abs(gradientConfig.maxNeg) : 1;
                    const t = Math.min(Math.abs(foldValue) / denom, 1);
                    const blended = white.map((w, idx) => clampChannel((1 - t) * w + t * gradientConfig.negColor[idx]));
                    return `rgb(${{blended.join(',')}})`;
                }}
                const denom = gradientConfig.maxPos > 0 ? gradientConfig.maxPos : (Math.abs(gradientConfig.maxPos) || 1);
                const t = Math.min(foldValue / denom, 1);
                const blended = white.map((w, idx) => clampChannel((1 - t) * w + t * gradientConfig.posColor[idx]));
                return `rgb(${{blended.join(',')}})`;
            }};
            const entityColor = (entity, idx = 1) => {{
                const foldVal = parseFoldChange(entity, idx);
                if (foldVal !== null) {{
                    return gradientColorFromFold(foldVal);
                }}
                return toRgbString(entity && entity[`fc_color_${{idx}}`]);
            }};
            const trackableMoveTypes = new Set(['prot-box','ptm-shape','ptm-label','ptm-symbol','compound','text-box','arrow','arrow-start','arrow-end']);
            const toCoordinateNumber = (value, fallback = 0) => {{
                const num = Number(value);
                return Number.isFinite(num) ? num : fallback;
            }};
            const ensureElementForId = (elementId) => {{
                if (!elementId) return null;
                return draw.find(`[data-id="${{elementId}}"]`)[0] || null;
            }};
            const captureElementPosition = (element, type) => {{
                if (!element) return null;
                if (type === 'prot-box') {{
                    return {{
                        x: toCoordinateNumber(element.x && element.x()),
                        y: toCoordinateNumber(element.y && element.y()),
                        width: toCoordinateNumber(element.width && element.width()),
                        height: toCoordinateNumber(element.height && element.height())
                    }};
                }}
                if (type === 'arrow') {{
                    const readAttr = (attr) => toCoordinateNumber(element.attr ? element.attr(attr) : null, 0);
                    return {{
                        x: readAttr('x1'),
                        y: readAttr('y1'),
                        x2: readAttr('x2'),
                        y2: readAttr('y2')
                    }};
                }}
                if (type === 'arrow-start' || type === 'arrow-end') {{
                    if (typeof element.cx === 'function' && typeof element.cy === 'function') {{
                        return {{
                            x: toCoordinateNumber(element.cx()),
                            y: toCoordinateNumber(element.cy())
                        }};
                    }}
                }}
                if (type === 'ptm-shape') {{
                    const isCircle = element.type === 'circle';
                    return {{
                        x: toCoordinateNumber(isCircle ? element.cx() : element.x()),
                        y: toCoordinateNumber(isCircle ? element.cy() : element.y()),
                        posKey: element.attr ? element.attr('data-pos-key') || null : null
                    }};
                }}
                if (typeof element.x === 'function') {{
                    return {{
                        x: toCoordinateNumber(element.x()),
                        y: toCoordinateNumber(element.y())
                    }};
                }}
                if (typeof element.cx === 'function') {{
                    return {{
                        x: toCoordinateNumber(element.cx()),
                        y: toCoordinateNumber(element.cy())
                    }};
                }}
                return null;
            }};
            const positionsDiffer = (a, b) => {{
                if (!a || !b) return false;
                if (Math.abs(a.x - b.x) > 0.4 || Math.abs(a.y - b.y) > 0.4) return true;
                if ((a.x2 !== undefined || b.x2 !== undefined) && (Math.abs((a.x2 || 0) - (b.x2 || 0)) > 0.4)) return true;
                if ((a.y2 !== undefined || b.y2 !== undefined) && (Math.abs((a.y2 || 0) - (b.y2 || 0)) > 0.4)) return true;
                if ((a.posKey || null) !== (b.posKey || null)) return true;
                return false;
            }};
            const runWithTemporarySelection = (element, type, id, protboxId, fn) => {{
                if (!element || typeof fn !== 'function') return;
                const prevSelected = selectedElement;
                const prevType = selectedType;
                const prevId = selectedId;
                const prevProtbox = selectedProtboxId;
                const prevVisual = selectedVisual;
                try {{
                    let resolvedElement = element;
                    let resolvedVisual = null;
                    if (type === 'arrow' && id) {{
                        resolvedVisual = draw.find(`[data-id="${{id}}"]`)[0] || null;
                        const hit = draw.find(`[data-id="${{id}}_hit"]`)[0];
                        if (hit) {{
                            resolvedElement = hit;
                        }} else if (resolvedVisual) {{
                            resolvedElement = resolvedVisual;
                        }}
                    }}
                    selectedElement = resolvedElement;
                    selectedType = type;
                    selectedId = id;
                    selectedProtboxId = protboxId || null;
                    selectedVisual = resolvedVisual;
                    fn();
                }} finally {{
                    selectedElement = prevSelected;
                    selectedType = prevType;
                    selectedId = prevId;
                    selectedProtboxId = prevProtbox;
                    selectedVisual = prevVisual;
                }}
            }};
            const cloneData = (payload) => {{
                if (!payload) return null;
                try {{
                    return JSON.parse(JSON.stringify(payload));
                }} catch (err) {{
                    console.warn('mkHistory: clone failed', err);
                    return null;
                }}
            }};
            const resolveLabelCenteringFromElement = (element) => {{
                if (!element || !element.attr) return 'center';
                const anchorAttr = element.attr('text-anchor');
                if (anchorAttr === 'start') return 'left';
                if (anchorAttr === 'end') return 'right';
                return 'center';
            }};
            const captureArrowSnapshot = (arrowId) => {{
                if (!arrowId) return null;
                const arrowElement = draw.find(`[data-id="${{arrowId}}"]`)[0];
                if (!arrowElement) return null;
                const readCoord = (attr) => parseFloat(arrowElement.attr(attr) || 0);
                return {{
                    arrowId,
                    x1: readCoord('x1'),
                    y1: readCoord('y1'),
                    x2: readCoord('x2'),
                    y2: readCoord('y2'),
                    lineType: arrowElement.attr('data-line-type') || 'arrow',
                }};
            }};
            const captureAttachedArrowSnapshots = (protboxId) => {{
                if (!protboxId || !attachments[protboxId]) return [];
                const seen = new Set();
                const snapshots = [];
                Object.values(attachments[protboxId]).forEach(sideAttachments => {{
                    sideAttachments.forEach(att => {{
                        const arrowId = att?.arrowId;
                        const endType = att?.type;
                        if (!arrowId || !endType) return;
                        const key = `${{arrowId}}_${{endType}}`;
                        if (seen.has(key)) return;
                        seen.add(key);
                        const arrowElement = draw.find(`[data-id="${{arrowId}}"]`)[0];
                        if (!arrowElement) return;
                        const coordXAttr = endType === 'start' ? 'x1' : 'x2';
                        const coordYAttr = endType === 'start' ? 'y1' : 'y2';
                        snapshots.push({{
                            arrowId,
                            endType,
                            x: toCoordinateNumber(arrowElement.attr(coordXAttr), 0),
                            y: toCoordinateNumber(arrowElement.attr(coordYAttr), 0),
                            lineType: arrowElement.attr('data-line-type') || 'arrow'
                        }});
                    }});
                }});
                return snapshots;
            }};
            const restoreArrowSnapshots = (snapshots) => {{
                if (!Array.isArray(snapshots) || !snapshots.length) return;
                snapshots.forEach(snap => {{
                    if (!snap?.arrowId) return;
                    const arrowElement = draw.find(`[data-id="${{snap.arrowId}}"]`)[0];
                    if (!arrowElement) return;
                    let x1 = toCoordinateNumber(arrowElement.attr('x1'), 0);
                    let y1 = toCoordinateNumber(arrowElement.attr('y1'), 0);
                    let x2 = toCoordinateNumber(arrowElement.attr('x2'), 0);
                    let y2 = toCoordinateNumber(arrowElement.attr('y2'), 0);
                    let appliedPartial = false;
                    if (typeof snap.endType === 'string') {{
                        if (snap.endType === 'start') {{
                            x1 = toCoordinateNumber(snap.x, x1);
                            y1 = toCoordinateNumber(snap.y, y1);
                            appliedPartial = true;
                        }} else if (snap.endType === 'end') {{
                            x2 = toCoordinateNumber(snap.x, x2);
                            y2 = toCoordinateNumber(snap.y, y2);
                            appliedPartial = true;
                        }}
                    }}
                    if (!appliedPartial) {{
                        x1 = toCoordinateNumber(snap.x1, x1);
                        y1 = toCoordinateNumber(snap.y1, y1);
                        x2 = toCoordinateNumber(snap.x2, x2);
                        y2 = toCoordinateNumber(snap.y2, y2);
                    }}
                    arrowElement.plot(x1, y1, x2, y2);
                    const arrowHitbox = draw.find(`[data-id="${{snap.arrowId}}_hit"]`)[0];
                    if (arrowHitbox) {{
                        arrowHitbox.plot(x1, y1, x2, y2);
                    }}
                    const startHandle = arrowHandleGroups[snap.arrowId]?.startHandle;
                    if (startHandle) {{
                        startHandle.cx(x1).cy(y1);
                    }}
                    const endHandle = arrowHandleGroups[snap.arrowId]?.endHandle;
                    if (endHandle) {{
                        endHandle.cx(x2).cy(y2);
                    }}
                    const arrowHead = draw.find(`[data-id="${{snap.arrowId}}_head"]`)[0];
                    if (!arrowHead) return;
                    const dx = x2 - x1;
                    const dy = y2 - y1;
                    const angle = Math.atan2(dy, dx);
                    const arrowSize = 5;
                    const lineType = snap.lineType || arrowElement.attr('data-line-type') || 'arrow';
                    if (lineType === 'inhibition') {{
                        const arrowData = getArrowById(snap.arrowId);
                        const endSide = arrowData?.protbox_id_2_side;
                        let barX1, barY1, barX2, barY2;
                        if (endSide === 'North' || endSide === 'South') {{
                            barX1 = x2 - arrowSize;
                            barX2 = x2 + arrowSize;
                            barY1 = barY2 = y2;
                        }} else if (endSide === 'West' || endSide === 'East') {{
                            barY1 = y2 - arrowSize;
                            barY2 = y2 + arrowSize;
                            barX1 = barX2 = x2;
                        }} else {{
                            const perp = angle + Math.PI / 2;
                            barX1 = x2 + Math.cos(perp) * arrowSize;
                            barY1 = y2 + Math.sin(perp) * arrowSize;
                            barX2 = x2 - Math.cos(perp) * arrowSize;
                            barY2 = y2 - Math.sin(perp) * arrowSize;
                        }}
                        arrowHead.plot(barX1, barY1, barX2, barY2);
                    }} else {{
                        arrowHead.plot([
                            [x2, y2],
                            [x2 - arrowSize * Math.cos(angle + Math.PI / 6), y2 - arrowSize * Math.sin(angle + Math.PI / 6)],
                            [x2 - arrowSize * Math.cos(angle - Math.PI / 6), y2 - arrowSize * Math.sin(angle - Math.PI / 6)],
                        ]);
                    }}
                }});
            }};
            const createHistoryManager = () => {{
                const undoStack = [];
                const redoStack = [];
                let suppress = false;
                let activeMove = null;
                let undoButton = null;
                let redoButton = null;
                const detachButtons = () => {{
                    if (undoButton) {{
                        undoButton.removeEventListener('click', undoAction);
                        undoButton = null;
                    }}
                    if (redoButton) {{
                        redoButton.removeEventListener('click', redoAction);
                        redoButton = null;
                    }}
                }};
                const updateButtons = () => {{
                    if (undoButton) undoButton.disabled = !undoStack.length;
                    if (redoButton) redoButton.disabled = !redoStack.length;
                }};
                const applyMoveEntry = (entry, direction) => {{
                    const target = direction === 'undo' ? entry.before : entry.after;
                    const element = ensureElementForId(entry.id);
                    if (!element || !target) return;
                    runWithTemporarySelection(element, entry.targetType, entry.id, entry.protboxId, () => {{
                        suppress = true;
                        try {{
                            const current = captureElementPosition(element, entry.targetType);
                            if (!current) return;
                            const deltaX = target.x - current.x;
                            const deltaY = target.y - current.y;
                            if (Math.abs(deltaX) > 0.05 || Math.abs(deltaY) > 0.05) {{
                                moveSelectedElement(deltaX, deltaY);
                            }}
                            if (entry.targetType === 'ptm-shape') {{
                                if (target.posKey) {{
                                    element.attr && element.attr('data-pos-key', target.posKey);
                                }} else if (element.attr) {{
                                    element.attr('data-pos-key', null);
                                    try {{ element.node.removeAttribute('data-pos-key'); }} catch (ignore) {{}}
                                }}
                            }}
                        }} finally {{
                            suppress = false;
                        }}
                    }});
                    if (entry.attachedArrows) {{
                        const snaps = direction === 'undo' ? entry.attachedArrows.before : entry.attachedArrows.after;
                        restoreArrowSnapshots(snaps);
                    }}
                }};
                const applyEntry = (entry, direction) => {{
                    if (!entry) return;
                    if (entry.kind === 'move') {{
                        applyMoveEntry(entry, direction);
                        return;
                    }}
                    if (entry.handlers) {{
                        const handler = direction === 'undo' ? entry.handlers.undo : entry.handlers.redo;
                        if (typeof handler === 'function') {{
                            handler();
                        }}
                    }}
                }};
                const recordAction = (entry) => {{
                    if (!entry) return;
                    undoStack.push(entry);
                    redoStack.length = 0;
                    updateButtons();
                }};
                const undoAction = () => {{
                    if (!undoStack.length) return;
                    const entry = undoStack.pop();
                    suppress = true;
                    try {{
                        applyEntry(entry, 'undo');
                    }} finally {{
                        suppress = false;
                    }}
                    redoStack.push(entry);
                    updateButtons();
                }};
                const redoAction = () => {{
                    if (!redoStack.length) return;
                    const entry = redoStack.pop();
                    suppress = true;
                    try {{
                        applyEntry(entry, 'redo');
                    }} finally {{
                        suppress = false;
                    }}
                    undoStack.push(entry);
                    updateButtons();
                }};
                const beginMoveSession = (meta) => {{
                    if (suppress || !meta || !trackableMoveTypes.has(meta.type)) return;
                    const element = meta.element || ensureElementForId(meta.id);
                    if (!element) return;
                    const before = captureElementPosition(element, meta.type);
                    const attachedBefore = (meta.type === 'prot-box' && meta.protboxId) ? captureAttachedArrowSnapshots(meta.protboxId) : null;
                    activeMove = before ? {{
                        kind: 'move',
                        targetType: meta.type,
                        id: meta.id,
                        protboxId: meta.protboxId || null,
                        before,
                        attachedArrows: attachedBefore ? {{ before: attachedBefore }} : null,
                    }} : null;
                }};
                const finalizeMoveSession = () => {{
                    if (!activeMove) return;
                    const element = ensureElementForId(activeMove.id);
                    const after = element ? captureElementPosition(element, activeMove.targetType) : null;
                    if (!after || !positionsDiffer(activeMove.before, after)) {{
                        activeMove = null;
                        return;
                    }}
                    const entry = {{
                        kind: 'move',
                        targetType: activeMove.targetType,
                        id: activeMove.id,
                        protboxId: activeMove.protboxId,
                        before: activeMove.before,
                        after
                    }};
                    if (activeMove.attachedArrows && activeMove.protboxId) {{
                        const afterSnapshots = captureAttachedArrowSnapshots(activeMove.protboxId);
                        const beforeSnapshots = activeMove.attachedArrows.before || [];
                        if ((beforeSnapshots && beforeSnapshots.length) || (afterSnapshots && afterSnapshots.length)) {{
                            entry.attachedArrows = {{
                                before: beforeSnapshots,
                                after: afterSnapshots
                            }};
                        }}
                    }}
                    recordAction(entry);
                    activeMove = null;
                }};
                const captureInstantMove = (meta, mover) => {{
                    if (!meta || typeof mover !== 'function' || !trackableMoveTypes.has(meta.type)) {{
                        mover && mover();
                        return;
                    }}
                    const element = meta.element || ensureElementForId(meta.id);
                    if (!element) {{
                        mover();
                        return;
                    }}
                    const before = captureElementPosition(element, meta.type);
                    let attachedSnapshots = null;
                    if (meta.type === 'prot-box' && meta.protboxId) {{
                        attachedSnapshots = {{ before: captureAttachedArrowSnapshots(meta.protboxId) }};
                    }}
                    suppress = true;
                    try {{
                        mover();
                    }} finally {{
                        suppress = false;
                    }}
                    const after = captureElementPosition(element, meta.type);
                    if (!before || !after || !positionsDiffer(before, after)) {{
                        return;
                    }}
                    if (attachedSnapshots && meta.protboxId) {{
                        attachedSnapshots.after = captureAttachedArrowSnapshots(meta.protboxId);
                        const beforeCount = attachedSnapshots.before ? attachedSnapshots.before.length : 0;
                        const afterCount = attachedSnapshots.after ? attachedSnapshots.after.length : 0;
                        if (!beforeCount && !afterCount) {{
                            attachedSnapshots = null;
                        }}
                    }}
                    const entry = {{
                        kind: 'move',
                        targetType: meta.type,
                        id: meta.id,
                        protboxId: meta.protboxId || null,
                        before,
                        after
                    }};
                    if (attachedSnapshots) {{
                        entry.attachedArrows = attachedSnapshots;
                    }}
                    recordAction(entry);
                }};
                const runWithoutRecording = (fn) => {{
                    if (suppress || typeof fn !== 'function') {{
                        return fn?.();
                    }}
                    suppress = true;
                    try {{
                        return fn();
                    }} finally {{
                        suppress = false;
                    }}
                }};
                const attachButtons = (buttons = {{}}) => {{
                    if (buttons.undo) {{
                        if (undoButton) undoButton.removeEventListener('click', undoAction);
                        undoButton = buttons.undo;
                        undoButton.addEventListener('click', undoAction);
                    }}
                    if (buttons.redo) {{
                        if (redoButton) redoButton.removeEventListener('click', redoAction);
                        redoButton = buttons.redo;
                        redoButton.addEventListener('click', redoAction);
                    }}
                    updateButtons();
                }};
                const dispose = () => {{
                    detachButtons();
                    undoStack.length = 0;
                    redoStack.length = 0;
                    activeMove = null;
                }};
                return {{
                    beginMoveSession,
                    finalizeMoveSession,
                    captureInstantMove,
                    runWithoutRecording,
                    attachButtons,
                    recordAction,
                    undo: undoAction,
                    redo: redoAction,
                    isSuppressed: () => suppress,
                    dispose
                }};
            }};

            // Create or reuse a single tooltip element so repeated renders don't append duplicates
            let tooltip = document.getElementById('mk-tooltip');
            if (!tooltip) {{
                tooltip = document.createElement('div');
                tooltip.id = 'mk-tooltip';
            }}
            // offsets (in px) for tooltip and context menus to position them near the cursor
            // negative Y values move the element upward relative to the cursor
            // increase negative Y offsets to move elements further up
            // adjusted: make tooltips appear higher (more negative Y)
            const tooltipOffsetX = 6, tooltipOffsetY = -36, menuOffsetX = 4, menuOffsetY = -24;
            const toContainerPosition = (evt) => {{
                const point = evt && evt.touches && evt.touches[0] ? evt.touches[0] : evt;
                if (!point) {{
                    return {{ x: 0, y: 0 }};
                }}
                const rect = container.getBoundingClientRect();
                const scrollLeft = container.scrollLeft || 0;
                const scrollTop = container.scrollTop || 0;
                return {{
                    x: (point.clientX - rect.left) + scrollLeft,
                    y: (point.clientY - rect.top) + scrollTop,
                }};
            }};
            const removeExistingContextMenu = () => {{
                const activeMenu = document.querySelector('.context-menu');
                if (activeMenu) activeMenu.remove();
            }};
            const showPtmContextMenu = (evt, protboxId, meta) => {{
                if (!meta || !meta.uniprot || !meta.ptmKey || !protboxId) {{
                    return;
                }}
                evt.preventDefault();
                removeExistingContextMenu();
                const menu = document.createElement('div');
                menu.className = 'context-menu';
                menu.style.position = 'absolute';
                const menuPos = toContainerPosition(evt);
                menu.style.left = `${{menuPos.x + menuOffsetX}}px`;
                menu.style.top = `${{menuPos.y + menuOffsetY}}px`;
                menu.style.backgroundColor = 'white';
                menu.style.border = '1px solid #ccc';
                menu.style.padding = '5px';
                menu.style.zIndex = '1000';
                const removeItem = document.createElement('div');
                removeItem.textContent = 'Remove';
                removeItem.style.padding = '5px 10px';
                removeItem.style.cursor = 'pointer';
                const cleanupMenu = () => {{
                    if (menu.parentNode) menu.remove();
                    document.removeEventListener('click', cleanupMenu);
                }};
                removeItem.addEventListener('click', () => {{
                    deletePtmFromProtbox(protboxId, meta);
                    cleanupMenu();
                }});
                menu.appendChild(removeItem);
                container.appendChild(menu);
                document.addEventListener('click', cleanupMenu);
            }};
            const bindPtmContextMenu = (svgElement, protboxId, meta) => {{
                if (!svgElement || !svgElement.node) return;
                svgElement.node.addEventListener('contextmenu', evt => showPtmContextMenu(evt, protboxId, meta));
            }};
            const isBackgroundTarget = (target) => {{
                if (!target) return false;
                if (target === draw.node) return true;
                if (target.getAttribute && target.getAttribute('data-kegg-bg') === '1') return true;
                return false;
            }};
            Object.assign(tooltip.style, {{position: 'absolute', backgroundColor: 'rgba(0,0,0,0.8)', color: 'white',padding: '5px 10px', borderRadius: '4px', fontSize: '12px', pointerEvents: 'none',display: 'none', maxWidth: '300px', whiteSpace: 'normal', wordWrap: 'break-word'}});
            if (!document.getElementById('mk-tooltip')) container.appendChild(tooltip);
            const setTooltipContent = (content, useHtml = false) => {{
                if (!tooltip) return;
                if (useHtml) {{
                    tooltip.innerHTML = content || '';
                }} else {{
                    tooltip.textContent = content || '';
                }}
                if (!content) {{
                    tooltip.innerHTML = '';
                }}
            }};
            const clearTooltipContent = () => {{
                if (!tooltip) return;
                tooltip.innerHTML = '';
            }};
            const elementGroups = {{}};
            const foregroundGroup = draw.group();
            const handleGroup = foregroundGroup.group();
            const arrowGroup = foregroundGroup.group();
            const protboxGroup = foregroundGroup.group();
            const compoundGroup = draw.group();  
            const textGroup = foregroundGroup.group();     
            // Keep the background fixed; shift all interactive layers by the inverse offset
            foregroundGroup.translate(fgOffsetX, fgOffsetY);
            if (!showArrows) {{
                arrowGroup.clear();
                handleGroup.clear();
                arrowHandleGroups = {{}};
            }}

            const updateHandleDists = (protboxId) => {{
                let northHas = false, southHas = false, westHas = false, eastHas = false;
                const snapPoints = protboxSnapPoints[protboxId];
                if (!snapPoints) return;
                elementGroups[protboxId].forEach(el => {{
                    if (el.element.attr('data-type') === 'ptm-shape') {{
                        const isCircle = el.element.type === 'circle';
                        const ptmX = isCircle ? el.element.cx() : el.element.x() + (el.element.width() / 2);
                        const ptmY = isCircle ? el.element.cy() : el.element.y() + (el.element.height() / 2);
                        for (const key in snapPoints) {{
                            const sp = snapPoints[key];
                            const dist = Math.sqrt((ptmX - sp.x)**2 + (ptmY - sp.y)**2);
                            if (dist < 1) {{
                                if (key[0] === 'N') northHas = true;
                                else if (key[0] === 'S') southHas = true;
                                else if (key[0] === 'W') westHas = true;
                                else if (key[0] === 'E') eastHas = true;
                            }}
                        }}
                    }}
                }});
                protboxHandleDists[protboxId] = {{North: northHas ? 10 : 5,South: southHas ? 10 : 5,West: westHas ? 10 : 5,East: eastHas ? 10 : 5}};
                elementGroups[protboxId].forEach(el => {{
                    const type = el.element.attr('data-type');
                    const id_ = el.element.attr('data-id');
                    if (type === 'handle-line') {{
                        let dist = 5;
                        const pb = protboxMap[protboxId];
                        const x = pb.x, y = pb.y, width = pb.width, height = pb.height;
                        if (id_.endsWith('_north_line')) {{dist = protboxHandleDists[protboxId].North;el.element.plot(x, y - dist, x + width, y - dist);}} 
                        else if (id_.endsWith('_south_line')) {{dist = protboxHandleDists[protboxId].South;el.element.plot(x, y + height + dist, x + width, y + height + dist);}} 
                        else if (id_.endsWith('_west_line')) {{dist = protboxHandleDists[protboxId].West;el.element.plot(x - dist, y, x - dist, y + height);}} 
                        else if (id_.endsWith('_east_line')) {{dist = protboxHandleDists[protboxId].East;el.element.plot(x + width + dist, y, x + width + dist, y + height);}}
                    }}
                }});
            }};
            const showRelatedHandles = (arrowId, visited = new Set()) => {{
                if (visited.has(arrowId)) return;
                visited.add(arrowId);
                const group = arrowHandleGroups[arrowId];
                if (group) {{
                    group.startHandle?.show();
                    group.endHandle?.show();
                }}
                const arrow = getArrowById(arrowId);
                if (arrow) {{
                    if (arrow.attached_arrow_1) {{
                        const mh = arrowHandleGroups[arrow.attached_arrow_1]?.[`${{arrow.attached_end_1}}Handle`];
                        mh?.show();
                        showRelatedHandles(arrow.attached_arrow_1, visited);
                    }}
                    if (arrow.attached_arrow_2) {{
                        const mh = arrowHandleGroups[arrow.attached_arrow_2]?.[`${{arrow.attached_end_2}}Handle`];
                        mh?.show();
                        showRelatedHandles(arrow.attached_arrow_2, visited);
                    }}
                }}
                if (attachedBy[arrowId]) {{
                    ['start', 'end'].forEach(endType => {{
                        if (attachedBy[arrowId][endType]) {{
                            attachedBy[arrowId][endType].forEach(slaveInfo => {{
                                const sId = slaveInfo.slave;
                                const freeEnd = slaveInfo.slaveEnd === 'start' ? 'endHandle' : 'startHandle';
                                arrowHandleGroups[sId]?.[freeEnd]?.show();
                                showRelatedHandles(sId, visited);
                            }});
                        }}
                    }});
                }}
            }};
            const hideRelatedHandles = (arrowId, visited = new Set()) => {{
                if (visited.has(arrowId)) return;
                visited.add(arrowId);
                const group = arrowHandleGroups[arrowId];
                if (group) {{
                    group.startHandle?.hide();
                    group.endHandle?.hide();
                }}
                const arrow = getArrowById(arrowId);
                if (arrow) {{
                    if (arrow.attached_arrow_1) {{
                        const mh = arrowHandleGroups[arrow.attached_arrow_1]?.[`${{arrow.attached_end_1}}Handle`];
                        mh?.hide();
                        hideRelatedHandles(arrow.attached_arrow_1, visited);
                    }}
                    if (arrow.attached_arrow_2) {{
                        const mh = arrowHandleGroups[arrow.attached_arrow_2]?.[`${{arrow.attached_end_2}}Handle`];
                        mh?.hide();
                        hideRelatedHandles(arrow.attached_arrow_2, visited);
                    }}
                }}
                if (attachedBy[arrowId]) {{
                    ['start', 'end'].forEach(endType => {{
                        if (attachedBy[arrowId][endType]) {{
                            attachedBy[arrowId][endType].forEach(slaveInfo => {{
                                const sId = slaveInfo.slave;
                                const freeEnd = slaveInfo.slaveEnd === 'start' ? 'endHandle' : 'startHandle';
                                arrowHandleGroups[sId]?.[freeEnd]?.hide();
                                hideRelatedHandles(sId, visited);
                            }});
                        }}
                    }});
                }}
            }};
            const setGroupStrokeColor = (group, role, color = null) => {{
                if (!group || typeof group.findOne !== 'function') return;
                const target = group.findOne(`[data-role="${{role}}"]`);
                if (!target) return;
                const originalColor = target.attr('data-original-stroke') || '#000';
                const originalWidth = parseFloat(target.attr('data-original-stroke-width')) || 1;
                if (color) {{
                    target.stroke({{ color, width: originalWidth }});
                }} else {{
                    target.stroke({{ color: originalColor, width: originalWidth }});
                }}
            }};
            const setGroupLabelColor = (group, labelType, color = null) => {{
                if (!group || typeof group.findOne !== 'function') return;
                const label = group.findOne(`[data-type="${{labelType}}"]`);
                if (!label) return;
                const original = label.attr('data-original-fill') || '#000';
                label.fill(color || original);
            }};
            const updateCompoundSelectionStyle = (group, active) => {{
                setGroupStrokeColor(group, 'compound-shape', active ? 'red' : null);
                setGroupLabelColor(group, 'compound-label', active ? 'red' : null);
            }};
            const updateTextBoxSelectionStyle = (group, active) => {{
                setGroupStrokeColor(group, 'text-rect', active ? 'red' : null);
                setGroupLabelColor(group, 'text-label', active ? 'red' : null);
            }};
            const selectElement = (element, type, id, protboxId = null) => {{
                if (selectedElement) {{
                    if (selectedType === 'compound') {{
                        updateCompoundSelectionStyle(selectedElement, false);
                    }} else if (selectedType === 'text-box') {{
                        updateTextBoxSelectionStyle(selectedElement, false);
                    }} else if (selectedType === 'arrow') {{
                        selectedVisual.stroke({{ color: 'black', width: 1 }});
                    }} else if (selectedType === 'arrow-start' || selectedType === 'arrow-end') {{
                        selectedElement.stroke({{ color: 'black', width: 1 }});
                    }} else if (selectedType === 'ptm-label') {{
                        selectedElement.stroke({{ color: 'none' }});
                    }} else {{
                        selectedElement.stroke({{ color: 'black', width: 1 }});
                    }}
                    if ((selectedType === 'arrow' || selectedType === 'arrow-start' || selectedType === 'arrow-end') && selectedId) {{
                        const arrowId = selectedType === 'arrow' ? selectedId : selectedId.replace(/_(start|end)$/, '');
                        hideRelatedHandles(arrowId);
                    }}
                    if (selectedType === 'prot-box' && selectedId) {{
                        elementGroups[selectedId].forEach(el => (el.element.attr('data-type').startsWith('handle-') || el.element.attr('data-type') === 'ptm-snap-circle') && el.element.hide());
                    }}
                }}
                selectedElement = element;
                selectedType = type;
                selectedId = id;
                if (protboxId) {{
                    selectedProtboxId = protboxId;
                    if (elementGroups[protboxId]) {{
                        elementGroups[protboxId].forEach(el => el.element.attr('data-type') === 'ptm-snap-circle' && el.element.show());
                    }}
                }}
                if (type === 'arrow') {{
                    selectedVisual = draw.findOne(`[data-id="${{id}}"]`);
                    selectedVisual.stroke({{ color: 'red', width: 1 }});
                }} else if (type === 'arrow-start' || type === 'arrow-end') {{
                    element.fill('red').stroke({{ color: 'red', width: 2 }});
                    selectedVisual = null;
                }} else if (type === 'compound') {{
                    updateCompoundSelectionStyle(element, true);
                    selectedVisual = null;
                }} else if (type === 'text-box') {{
                    updateTextBoxSelectionStyle(element, true);
                    selectedVisual = null;
                }} else if (type !== 'ptm-label') {{
                    element.stroke({{ color: 'red', width: 1 }});
                    selectedVisual = null;
                }}
                if (type === 'prot-box' && id) {{
                    elementGroups[id].forEach(el => (el.element.attr('data-type').startsWith('handle-') || el.element.attr('data-type') === 'ptm-snap-circle') && el.element.show());
                }}
                if (type === 'arrow' || type === 'arrow-start' || type === 'arrow-end') {{
                    Object.values(elementGroups).forEach(group => group.forEach(el => el.element.attr('data-type').startsWith('handle-') && el.element.show()));
                    const arrowId = type === 'arrow' ? id : id.replace(/_(start|end)$/, '');
                    showRelatedHandles(arrowId);
                }}
            }};
            const deselectElement = () => {{
                if (selectedElement) {{
                    if (selectedType === 'compound') {{
                        updateCompoundSelectionStyle(selectedElement, false);
                    }} else if (selectedType === 'text-box') {{
                        updateTextBoxSelectionStyle(selectedElement, false);
                    }} else if (selectedType === 'arrow') {{
                        selectedVisual?.stroke({{ color: 'black', width: 1 }});
                    }} else if (selectedType === 'arrow-start' || selectedType === 'arrow-end') {{
                        selectedElement.stroke({{ color: 'black', width: 1 }});
                    }} else if (selectedType === 'ptm-label') {{
                        selectedElement.stroke({{ color: 'none' }});
                    }} else {{
                        selectedElement.stroke({{ color: 'black', width: 1 }});
                    }}
                    if ((selectedType === 'arrow' || selectedType === 'arrow-start' || selectedType === 'arrow-end') && selectedId) {{
                        const arrowId = selectedType === 'arrow' ? selectedId : selectedId.replace(/_(start|end)$/, '');
                        hideRelatedHandles(arrowId);
                        const arrow = getArrowById(arrowId);
                        if (arrow) {{
                            if (arrow.protbox_id_1 && arrow.protbox_id_1_side) updateArrowPositions(arrow.protbox_id_1, arrow.protbox_id_1_side);
                            if (arrow.protbox_id_2 && arrow.protbox_id_2_side) updateArrowPositions(arrow.protbox_id_2, arrow.protbox_id_2_side);
                        }}
                    }}
                    if (selectedType === 'prot-box' && selectedId) {{
                        const selectedGroup = elementGroups[selectedId];
                        if (selectedGroup) {{
                            selectedGroup.forEach(el => (el.element.attr('data-type').startsWith('handle-') || el.element.attr('data-type') === 'ptm-snap-circle') && el.element.hide());
                        }}
                    }}
                }}
                if (selectedProtboxId) {{
                    const protGroup = elementGroups[selectedProtboxId];
                    if (protGroup) {{
                        protGroup.forEach(el => el.element.attr('data-type') === 'ptm-snap-circle' && el.element.hide());
                    }}
                }}
                Object.values(elementGroups).forEach(group => group.forEach(el => el.element.attr('data-type').startsWith('handle-') && el.element.hide()));
                selectedElement = selectedType = selectedId = selectedProtboxId = selectedVisual = null;
            }};
            const handleWheel = e => {{
                e.preventDefault();
                if (e.ctrlKey) {{
                    const delta = e.deltaY < 0 ? 1.1 : 0.9;
                    const newZoom = zoomLevel * delta;
                    if (newZoom < minZoom || newZoom > maxZoom) return;
                    zoomLevel = newZoom;
                    const rect = document.getElementById('svgCanvas').getBoundingClientRect();
                    const mouseX = (e.clientX - rect.left) / zoomLevel;
                    const mouseY = (e.clientY - rect.top) / zoomLevel;
                    const newWidth = {max_x} / zoomLevel;
                    const newHeight = {max_y} / zoomLevel;
                    viewBox = {{ 
                        x: mouseX - (mouseX - viewBox.x) * (newWidth / viewBox.width),
                        y: mouseY - (mouseY - viewBox.y) * (newHeight / viewBox.height),
                        width: newWidth, 
                        height: newHeight 
                    }};
                    draw.viewbox(viewBox.x, viewBox.y, viewBox.width, viewBox.height);
                }} else {{
                    const panX = e.deltaX;
                    const panY = e.deltaY;
                    viewBox.x += panX / zoomLevel;
                    viewBox.y += panY / zoomLevel;
                    draw.viewbox(viewBox.x, viewBox.y, viewBox.width, viewBox.height);
                }}
            }};
            const startPanning = e => {{
                if (isBackgroundTarget(e.target)) {{
                    e.preventDefault();
                    isPanning = true;
                    startPanX = e.clientX;
                    startPanY = e.clientY;
                    startViewBoxX = viewBox.x;
                    startViewBoxY = viewBox.y;
                    document.getElementById('svgCanvas').style.cursor = 'grab';
                }}
            }};
            const pan = e => {{
                if (!isPanning) return;
                e.preventDefault();
                viewBox.x = startViewBoxX - (e.clientX - startPanX) / zoomLevel;
                viewBox.y = startViewBoxY - (e.clientY - startPanY) / zoomLevel;
                draw.viewbox(viewBox.x, viewBox.y, viewBox.width, viewBox.height);
            }};
            const endPanning = () => {{
                if (isPanning) {{
                    isPanning = false;
                    document.getElementById('svgCanvas').style.cursor = 'default';
                }}
            }};
            const resetView = () => {{
                zoomLevel = 1;
                viewBox = {{ x: 0, y: 0, width: {max_x}, height: {max_y} }};
                draw.viewbox(viewBox.x, viewBox.y, viewBox.width, viewBox.height);
            }};
            if (mkHistory && typeof mkHistory.dispose === 'function') {{
                mkHistory.dispose();
            }}
            mkHistory = null;
            mkHistory = createHistoryManager();
            const snapToHandle = (x, y, protboxId, side, isDashed = false) => {{
                const pb = protboxMap[protboxId];
                if (!pb) return {{ x, y }};
                const handleDist = isDashed ? 0 : protboxHandleDists[protboxId]?.[side] || 5;
                const px = pb.x || 0, py = pb.y || 0, width = pb.width || 46, height = pb.height || 17;
                let snapX = x, snapY = y, isSnapped = false;
                let handleX1, handleY1, handleX2, handleY2;
                if (side === 'North') {{handleX1 = px; handleY1 = py - handleDist; handleX2 = px + width; handleY2 = py - handleDist;}} 
                else if (side === 'South') {{handleX1 = px; handleY1 = py + height + handleDist; handleX2 = px + width; handleY2 = py + height + handleDist;}} 
                else if (side === 'West') {{handleX1 = px - handleDist; handleY1 = py; handleX2 = px - handleDist; handleY2 = py + height;}} 
                else if (side === 'East') {{handleX1 = px + width + handleDist; handleY1 = py; handleX2 = px + width + handleDist; handleY2 = py + height;}}
                if (side === 'North' || side === 'South') {{
                    if (Math.abs(y - handleY1) <= snapRadius && x >= handleX1 && x <= handleX2) {{
                        snapY = handleY1; snapX = Math.max(handleX1, Math.min(handleX2, x)); isSnapped = true;
                    }}
                }} else if (side === 'West' || side === 'East') {{
                    if (Math.abs(x - handleX1) <= snapRadius && y >= handleY1 && y <= handleY2) {{
                        snapX = handleX1; snapY = Math.max(handleY1, Math.min(handleY2, y)); isSnapped = true;
                    }}
                }}
                return {{ x: snapX, y: snapY, isSnapped, protbox_id: protboxId, side }};
            }};
            var snapPointToPreferredAxis = typeof snapPointToPreferredAxis === 'function' ? snapPointToPreferredAxis : function(point, anchor) {{
                if (!anchor) return point;
                const dx = point.x - anchor.x;
                const dy = point.y - anchor.y;
                if (!Number.isFinite(dx) || !Number.isFinite(dy)) {{
                    return point;
                }}
                if (Math.abs(dx) < 1e-6 && Math.abs(dy) < 1e-6) {{
                    return point;
                }}
                const angleDeg = (Math.atan2(dy, dx) * 180 / Math.PI + 360) % 360;
                let closestAxis = null;
                let smallestDiff = Infinity;
                AXIS_SNAP_ANGLES.forEach(axisAngle => {{
                    let diff = Math.abs(angleDeg - axisAngle);
                    if (diff > 180) diff = 360 - diff;
                    if (diff < smallestDiff) {{
                        smallestDiff = diff;
                        closestAxis = axisAngle;
                    }}
                }});
                if (closestAxis !== null && smallestDiff <= AXIS_SNAP_TOLERANCE_DEG) {{
                    if (closestAxis % 180 === 0) {{
                        return {{ x: point.x, y: anchor.y }};
                    }}
                    return {{ x: anchor.x, y: point.y }};
                }}
                return point;
            }};
            var isPointNearProtboxSnapZone = typeof isPointNearProtboxSnapZone === 'function' ? isPointNearProtboxSnapZone : function(x, y, arrow) {{
                if (!arrow) return false;
                for (const protboxId in protboxMap) {{
                    for (const side of ['North', 'South', 'West', 'East']) {{
                        const candidate = snapToHandle(x, y, protboxId, side, arrow.line === 'dashed_arrow');
                        if (candidate.isSnapped) {{
                            return true;
                        }}
                    }}
                }}
                return false;
            }};
            var isPointNearArrowSnapZone = typeof isPointNearArrowSnapZone === 'function' ? isPointNearArrowSnapZone : function(x, y, currentArrowId = null) {{
                if (!Array.isArray(arrows) || !arrows.length || !draw) {{
                    return false;
                }}
                for (let idx = 0; idx < arrows.length; idx++) {{
                    const arrowId = `arrow_${{idx}}`;
                    if (arrowId === currentArrowId) continue;
                    const otherLine = draw.find(`[data-id="${{arrowId}}"]`)[0];
                    if (!otherLine) continue;
                    const ox1 = parseFloat(otherLine.attr('x1'));
                    const oy1 = parseFloat(otherLine.attr('y1'));
                    const ox2 = parseFloat(otherLine.attr('x2'));
                    const oy2 = parseFloat(otherLine.attr('y2'));
                    if (Number.isFinite(ox1) && Number.isFinite(oy1) && Math.hypot(x - ox1, y - oy1) <= snapRadius) {{
                        return true;
                    }}
                    if (Number.isFinite(ox2) && Number.isFinite(oy2) && Math.hypot(x - ox2, y - oy2) <= snapRadius) {{
                        return true;
                    }}
                }}
                return false;
            }};
            const updateAttachedArrows = (arrowId, endType) => {{
                if (!attachedBy[arrowId] || !attachedBy[arrowId][endType]) return;
                const masterLine = draw.find(`[data-id="${{arrowId}}"]`)[0];
                const masterX = endType === 'start' ? parseFloat(masterLine.attr('x1')) : parseFloat(masterLine.attr('x2'));
                const masterY = endType === 'start' ? parseFloat(masterLine.attr('y1')) : parseFloat(masterLine.attr('y2'));
                attachedBy[arrowId][endType].forEach(slaveInfo => {{
                    const slaveArrowId = slaveInfo.slave;
                    const slaveEnd = slaveInfo.slaveEnd;
                    const slaveLine = draw.find(`[data-id="${{slaveArrowId}}"]`)[0];
                    const slaveHitbox = draw.find(`[data-id="${{slaveArrowId}}_hit"]`)[0];
                    if (slaveEnd === 'start') {{
                        slaveLine.attr({{ x1: masterX, y1: masterY }});
                    }} else {{
                        slaveLine.attr({{ x2: masterX, y2: masterY }});
                    }}
                    if (slaveHitbox) slaveHitbox.plot(slaveLine.attr('x1'), slaveLine.attr('y1'), slaveLine.attr('x2'), slaveLine.attr('y2'));
                    const slaveHead = draw.find(`[data-id="${{slaveArrowId}}_head"]`)[0];
                    if (slaveHead) {{
                        const sX1 = parseFloat(slaveLine.attr('x1')), sY1 = parseFloat(slaveLine.attr('y1'));
                        const sX2 = parseFloat(slaveLine.attr('x2')), sY2 = parseFloat(slaveLine.attr('y2'));
                        const dx = sX2 - sX1, dy = sY2 - sY1, angle = Math.atan2(dy, dx), arrowSize = 5;
                        const lineType = slaveLine.attr('data-line-type');
                        if (lineType === 'arrow') {{
                            slaveHead.plot([[sX2, sY2],[sX2 - arrowSize * Math.cos(angle + Math.PI / 6), sY2 - arrowSize * Math.sin(angle + Math.PI / 6)],[sX2 - arrowSize * Math.cos(angle - Math.PI / 6), sY2 - arrowSize * Math.sin(angle - Math.PI / 6)]]);
                        }} else if (lineType === 'inhibition') {{
                            const perp = angle + Math.PI / 2;
                            const barX1 = sX2 + Math.cos(perp) * arrowSize;
                            const barY1 = sY2 + Math.sin(perp) * arrowSize;
                            const barX2 = sX2 - Math.cos(perp) * arrowSize;
                            const barY2 = sY2 - Math.sin(perp) * arrowSize;
                            slaveHead.plot(barX1, barY1, barX2, barY2);
                        }}
                    }}
                    const slaveHandle = arrowHandleGroups[slaveArrowId]?.[`${{slaveEnd}}Handle`];
                    if (slaveHandle) {{
                        slaveHandle.cx(masterX).cy(masterY);
                    }}
                }});
            }};
            const updateArrowAttachments = (arrowId, endType, x, y) => {{
                const arrow = getArrowById(arrowId);
                if (!arrow) return {{ x, y, isSnapped: false }};
                const num = endType === 'start' ? 1 : 2;
                const oldMaster = arrow[`attached_arrow_${{num}}`];
                const oldMasterEnd = arrow[`attached_end_${{num}}`];
                if (oldMaster && oldMasterEnd) {{
                    if (attachedBy[oldMaster] && attachedBy[oldMaster][oldMasterEnd]) {{
                        attachedBy[oldMaster][oldMasterEnd] = attachedBy[oldMaster][oldMasterEnd].filter(s => s.slave !== arrowId);
                        if (attachedBy[oldMaster][oldMasterEnd].length === 0) delete attachedBy[oldMaster][oldMasterEnd];
                        if (Object.keys(attachedBy[oldMaster]).length === 0) delete attachedBy[oldMaster];
                    }}
                }}
                let snapped = {{ x, y, isSnapped: false }};
                for (let protboxId in protboxMap) {{
                    for (let side of ['North', 'South', 'West', 'East']) {{
                        const result = snapToHandle(x, y, protboxId, side, arrow.line === 'dashed_arrow');
                        if (result.isSnapped) {{
                            snapped = result;
                            break;
                        }}
                    }}
                    if (snapped.isSnapped) break;
                }}
                if (!snapped.isSnapped && arrow.line === 'line') {{
                    for (let otherIndex = 0; otherIndex < arrows.length; otherIndex++) {{
                        if (`arrow_${{otherIndex}}` === arrowId) continue;
                        const otherArrowId = `arrow_${{otherIndex}}`;
                        const otherLine = draw.find(`[data-id="${{otherArrowId}}"]`)[0];
                        if (!otherLine) continue;
                        const otherX1 = parseFloat(otherLine.attr('x1')), otherY1 = parseFloat(otherLine.attr('y1'));
                        const otherX2 = parseFloat(otherLine.attr('x2')), otherY2 = parseFloat(otherLine.attr('y2'));
                        const dist1 = Math.sqrt((x - otherX1)**2 + (y - otherY1)**2);
                        if (dist1 <= snapRadius) {{
                            snapped = {{x: otherX1, y: otherY1, isSnapped: true, arrow_id: otherArrowId, arrow_end: 'start'}};
                            break;
                        }}
                        const dist2 = Math.sqrt((x - otherX2)**2 + (y - otherY2)**2);
                        if (dist2 <= snapRadius) {{
                            snapped = {{x: otherX2, y: otherY2, isSnapped: true, arrow_id: otherArrowId, arrow_end: 'end'}};
                            break;
                        }}
                    }}
                }}
                const arrowElement = draw.find(`[data-id="${{arrowId}}"]`)[0];
                const arrowHitbox = draw.find(`[data-id="${{arrowId}}_hit"]`)[0];
                const arrowHead = draw.find(`[data-id="${{arrowId}}_head"]`)[0];
                const startHandle = arrowHandleGroups[arrowId]?.startHandle;
                const endHandle = arrowHandleGroups[arrowId]?.endHandle;
                let x1 = parseFloat(arrowElement.attr('x1') || 0);
                let y1 = parseFloat(arrowElement.attr('y1') || 0);
                let x2 = parseFloat(arrowElement.attr('x2') || 0);
                let y2 = parseFloat(arrowElement.attr('y2') || 0);
                if (endType === 'start') {{
                    x1 = snapped.x;
                    y1 = snapped.y;
                    if (startHandle) startHandle.cx(x1).cy(y1);
                }} else {{
                    x2 = snapped.x;
                    y2 = snapped.y;
                    if (endHandle) endHandle.cx(x2).cy(y2);
                }}
                arrowElement.plot(x1, y1, x2, y2);
                if (arrowHitbox) arrowHitbox.plot(x1, y1, x2, y2);
                if (snapped.isSnapped) {{
                    if (snapped.protbox_id) {{
                        const {{ protbox_id, side }} = snapped;
                        attachments[protbox_id] = attachments[protbox_id] || {{}};
                        attachments[protbox_id][side] = attachments[protbox_id][side] || [];
                        if (!attachments[protbox_id][side].some(att => att.arrow === arrow && att.type === endType)) {{
                            attachments[protbox_id][side].push({{ type: endType, arrow, arrowId }});
                        }}
                        arrow[`protbox_id_${{num}}`] = protbox_id;
                        arrow[`protbox_id_${{num}}_side`] = side;
                        delete arrow[`attached_arrow_${{num}}`];
                        delete arrow[`attached_end_${{num}}`];
                    }} else if (snapped.arrow_id) {{
                        arrow[`attached_arrow_${{num}}`] = snapped.arrow_id;
                        arrow[`attached_end_${{num}}`] = snapped.arrow_end;
                        delete arrow[`protbox_id_${{num}}`];
                        delete arrow[`protbox_id_${{num}}_side`];
                        const masterId = snapped.arrow_id;
                        const masterEnd = snapped.arrow_end;
                        attachedBy[masterId] = attachedBy[masterId] || {{start: [], end: []}};
                        attachedBy[masterId][masterEnd].push({{slave: arrowId, slaveEnd: endType}});
                    }}
                }} else {{
                    delete arrow[`attached_arrow_${{num}}`];
                    delete arrow[`attached_end_${{num}}`];
                    delete arrow[`protbox_id_${{num}}`];
                    delete arrow[`protbox_id_${{num}}_side`];
                    Object.keys(attachments).forEach(protboxId => {{
                        Object.keys(attachments[protboxId]).forEach(side => {{
                            attachments[protboxId][side] = attachments[protboxId][side].filter(att => !(att.arrow === arrow && att.type === endType));
                            if (attachments[protboxId][side].length === 0) delete attachments[protboxId][side];
                        }});
                        if (Object.keys(attachments[protboxId]).length === 0) delete attachments[protboxId];
                    }});
                }}
                if (arrowHead) {{
                    const lineType = arrowElement.attr('data-line-type') || 'arrow';
                    const dx = x2 - x1, dy = y2 - y1, angle = Math.atan2(dy, dx), arrowSize = 5;
                    if (lineType === 'arrow') {{
                        arrowHead.plot([[x2, y2],[x2 - arrowSize * Math.cos(angle + Math.PI / 6), y2 - arrowSize * Math.sin(angle + Math.PI / 6)],[x2 - arrowSize * Math.cos(angle - Math.PI / 6), y2 - arrowSize * Math.sin(angle - Math.PI / 6)]]);
                    }} else if (lineType === 'inhibition') {{
                        const endSide = snapped.isSnapped ? snapped.side : arrow.protbox_id_2_side;
                        let barX1, barY1, barX2, barY2;
                        if (endSide && (endSide === 'North' || endSide === 'South')) {{barX1 = x2 - 5;barY1 = y2;barX2 = x2 + 5;barY2 = y2;}} 
                        else if (endSide && (endSide === 'West' || endSide === 'East')) {{barX1 = x2;barY1 = y2 - 5;barX2 = x2;barY2 = y2 + 5;}} 
                        else {{const perp = angle + Math.PI / 2;barX1 = x2 + Math.cos(perp) * arrowSize;barY1 = y2 + Math.sin(perp) * arrowSize;barX2 = x2 - Math.cos(perp) * arrowSize;barY2 = y2 - Math.sin(perp) * arrowSize;}}
                        arrowHead.plot(barX1, barY1, barX2, barY2);
                    }}
                }}
                const startProtboxId = arrow.protbox_id_1;
                const startSide = arrow.protbox_id_1_side;
                const endProtboxId = arrow.protbox_id_2;
                const endSide = arrow.protbox_id_2_side;
                if (startProtboxId && startSide) updateArrowPositions(startProtboxId, startSide);
                if (endProtboxId && endSide) updateArrowPositions(endProtboxId, endSide);
                Shiny?.setInputValue('arrow_moved', {{ id: arrowId, end: endType, x: endType === 'start' ? x1 : x2, y: endType === 'start' ? y1 : y2 }}, {{ priority: 'event' }});
                return snapped;
            }};
            const updateArrowPositions = (protboxId, side) => {{
                if (!attachments[protboxId]?.[side]) return;
                const attachList = attachments[protboxId][side];
                const pb = protboxMap[protboxId];
                if (!pb) return;
                attachList.forEach(attachment => {{
                    const {{ arrow, type: endType, arrowId }} = attachment;
                    const arrowElement = draw.find(`[data-id="${{arrowId}}"]`)[0];
                    const arrowHitbox = draw.find(`[data-id="${{arrowId}}_hit"]`)[0];
                    const arrowHead = draw.find(`[data-id="${{arrowId}}_head"]`)[0];
                    const id1 = arrow.protbox_id_1, side1 = arrow.protbox_id_1_side, box1 = protboxMap[id1];
                    const id2 = arrow.protbox_id_2, side2 = arrow.protbox_id_2_side, box2 = protboxMap[id2];
                    let x1, y1;
                    if (arrow.attached_arrow_1) {{
                        const attId = arrow.attached_arrow_1;
                        const attEnd = arrow.attached_end_1;
                        const attLine = draw.find(`[data-id="${{attId}}"]`)[0];
                        x1 = parseFloat(attLine.attr(attEnd === 'start' ? 'x1' : 'x2'));
                        y1 = parseFloat(attLine.attr(attEnd === 'start' ? 'y1' : 'y2'));
                    }} else if (id1 && side1 && box1) {{
                        let frac1 = arrow.line === 'dashed_arrow' ? 0.5 : (attachments[id1][side1].findIndex(att => att.arrow === arrow) + 1) / (attachments[id1][side1].length + 1);
                        const width1 = box1.width || 46, height1 = box1.height || 17;
                        const handleDist1 = arrow.line === 'dashed_arrow' ? 0 : protboxHandleDists[id1]?.[side1] || 5;
                        if (side1 === 'North') {{ x1 = box1.x + frac1 * width1; y1 = box1.y - handleDist1; }}
                        else if (side1 === 'South') {{ x1 = box1.x + frac1 * width1; y1 = box1.y + height1 + handleDist1; }}
                        else if (side1 === 'West') {{ x1 = box1.x - handleDist1; y1 = box1.y + frac1 * height1; }}
                        else if (side1 === 'East') {{ x1 = box1.x + width1 + handleDist1; y1 = box1.y + frac1 * height1; }}
                    }} else {{
                        x1 = arrow.x1 || parseFloat(arrowElement.attr('x1') || 0);
                        y1 = arrow.y1 || parseFloat(arrowElement.attr('y1') || 0);
                    }}
                    let x2, y2;
                    if (arrow.attached_arrow_2) {{
                        const attId = arrow.attached_arrow_2;
                        const attEnd = arrow.attached_end_2;
                        const attLine = draw.find(`[data-id="${{attId}}"]`)[0];
                        x2 = parseFloat(attLine.attr(attEnd === 'start' ? 'x1' : 'x2'));
                        y2 = parseFloat(attLine.attr(attEnd === 'start' ? 'y1' : 'y2'));
                    }} else if (id2 && side2 && box2) {{
                        let frac2 = arrow.line === 'dashed_arrow' ? 0.5 : (attachments[id2][side2].findIndex(att => att.arrow === arrow) + 1) / (attachments[id2][side2].length + 1);
                        const width2 = box2.width || 46, height2 = box2.height || 17;
                        const handleDist2 = arrow.line === 'dashed_arrow' ? 0 : protboxHandleDists[id2]?.[side2] || 5;
                        if (side2 === 'North') {{ x2 = box2.x + frac2 * width2; y2 = box2.y - handleDist2; }}
                        else if (side2 === 'South') {{ x2 = box2.x + frac2 * width2; y2 = box2.y + height2 + handleDist2; }}
                        else if (side2 === 'West') {{ x2 = box2.x - handleDist2; y2 = box2.y + frac2 * height2; }}
                        else if (side2 === 'East') {{ x2 = box2.x + width2 + handleDist2; y2 = box2.y + frac2 * height2; }}
                    }} else {{
                        x2 = arrow.x2 || parseFloat(arrowElement.attr('x2') || 0);
                        y2 = arrow.y2 || parseFloat(arrowElement.attr('y2') || 0);
                    }}
                    arrowElement.plot(x1, y1, x2, y2);
                    if (arrowHitbox) arrowHitbox.plot(x1, y1, x2, y2);
                    const startHandle = arrowHandleGroups[arrowId]?.startHandle;
                    const endHandle = arrowHandleGroups[arrowId]?.endHandle;
                    if (startHandle) startHandle.cx(x1).cy(y1);
                    if (endHandle) endHandle.cx(x2).cy(y2);
                    if (arrowHead) {{
                        const dx = x2 - x1, dy = y2 - y1, angle = Math.atan2(dy, dx), arrowSize = 5;
                        const lineType = arrowElement.attr('data-line-type') || 'arrow';
                        if (lineType === 'arrow') {{
                            arrowHead.plot([[x2, y2],[x2 - arrowSize * Math.cos(angle + Math.PI / 6), y2 - arrowSize * Math.sin(angle + Math.PI / 6)],[x2 - arrowSize * Math.cos(angle - Math.PI / 6), y2 - arrowSize * Math.sin(angle - Math.PI / 6)]]);
                        }} else if (lineType === 'inhibition') {{
                            const endSide = arrow.protbox_id_2_side;
                            let barX1, barY1, barX2, barY2;
                            if (endSide === 'North' || endSide === 'South') {{barX1 = x2 - arrowSize;barY1 = y2;barX2 = x2 + arrowSize;barY2 = y2;}} 
                            else if (endSide === 'West' || endSide === 'East') {{barX1 = x2;barY1 = y2 - arrowSize;barX2 = x2;barY2 = y2 + arrowSize;}} 
                            else {{const perp = angle + Math.PI / 2;barX1 = x2 + Math.cos(perp) * arrowSize;barY1 = y2 + Math.sin(perp) * arrowSize;barX2 = x2 - Math.cos(perp) * arrowSize;barY2 = y2 - Math.sin(perp) * arrowSize;}}
                            arrowHead.plot(barX1, barY1, barX2, barY2);
                        }}
                    }}
                }});
            }};
            const moveSelectedElement = (deltaX, deltaY, pointerX = null, pointerY = null) => {{
                if (!selectedElement) return;
                const isCircle = selectedType === 'ptm-shape';
                let currentX = parseFloat(isCircle ? selectedElement.cx() : selectedElement.x() || 0);
                let currentY = parseFloat(isCircle ? selectedElement.cy() : selectedElement.y() || 0);
                let newX = currentX + deltaX, newY = currentY + deltaY;
                let arrowHandleContext = null;
                if (selectedType === 'arrow-start' || selectedType === 'arrow-end') {{
                    const arrowId = selectedId.replace(/_(start|end)$/, '');
                    const arrowElement = draw.find(`[data-id="${{arrowId}}"]`)[0];
                    const arrowHitbox = draw.find(`[data-id="${{arrowId}}_hit"]`)[0];
                    const arrowHead = draw.find(`[data-id="${{arrowId}}_head"]`)[0];
                    const arrowData = getArrowById(arrowId);
                    arrowHandleContext = {{
                        arrowId,
                        arrowElement,
                        arrowHitbox,
                        arrowHead,
                        arrowData
                    }};
                    if (arrowElement) {{
                        const anchorAttrX = selectedType === 'arrow-start' ? 'x2' : 'x1';
                        const anchorAttrY = selectedType === 'arrow-start' ? 'y2' : 'y1';
                        arrowHandleContext.anchor = {{
                            x: parseFloat(arrowElement.attr(anchorAttrX) || 0),
                            y: parseFloat(arrowElement.attr(anchorAttrY) || 0)
                        }};
                        if (arrowData && arrowHandleContext.anchor) {{
                            const nearProtbox = isPointNearProtboxSnapZone(newX, newY, arrowData);
                            const nearArrow = arrowData.line === 'line' ? isPointNearArrowSnapZone(newX, newY, arrowId) : false;
                            if (!nearProtbox && !nearArrow) {{
                                const snappedPoint = snapPointToPreferredAxis({{ x: newX, y: newY }}, arrowHandleContext.anchor);
                                newX = snappedPoint.x;
                                newY = snappedPoint.y;
                            }}
                        }}
                    }}
                }}
                if (selectedType === 'arrow') {{
                    const arrowId = selectedId;
                    const hitbox = selectedElement;
                    const visual = selectedVisual;
                    const arrowHead = draw.find(`[data-id="${{arrowId}}_head"]`)[0];
                    const startHandle = arrowHandleGroups[arrowId]?.startHandle;
                    const endHandle = arrowHandleGroups[arrowId]?.endHandle;
                    let x1 = parseFloat(hitbox.attr('x1') || 0) + deltaX;
                    let y1 = parseFloat(hitbox.attr('y1') || 0) + deltaY;
                    let x2 = parseFloat(hitbox.attr('x2') || 0) + deltaX;
                    let y2 = parseFloat(hitbox.attr('y2') || 0) + deltaY;
                    hitbox.plot(x1, y1, x2, y2);
                    visual.plot(x1, y1, x2, y2);
                    if (startHandle) startHandle.cx(x1).cy(y1);
                    if (endHandle) endHandle.cx(x2).cy(y2);
                    if (arrowHead) {{
                        const lineType = visual.attr('data-line-type') || 'arrow';
                        const dx = x2 - x1, dy = y2 - y1, angle = Math.atan2(dy, dx), arrowSize = 5;
                        if (lineType === 'arrow') {{
                            arrowHead.plot([[x2, y2],[x2 - arrowSize * Math.cos(angle + Math.PI / 6), y2 - arrowSize * Math.sin(angle + Math.PI / 6)],[x2 - arrowSize * Math.cos(angle - Math.PI / 6), y2 - arrowSize * Math.sin(angle - Math.PI / 6)]]);
                        }} else if (lineType === 'inhibition') {{
                            const arrow = getArrowById(arrowId);
                            const endSide = arrow?.protbox_id_2_side;
                            let barX1, barY1, barX2, barY2;
                            if (endSide === 'North' || endSide === 'South') {{barX1 = x2 - arrowSize;barY1 = y2;barX2 = x2 + arrowSize;barY2 = y2;}} 
                            else if (endSide === 'West' || endSide === 'East') {{barX1 = x2;barY1 = y2 - arrowSize;barX2 = x2;barY2 = y2 + arrowSize;}} 
                            else {{const perp = angle + Math.PI / 2;barX1 = x2 + Math.cos(perp) * arrowSize;barY1 = y2 + Math.sin(perp) * arrowSize;barX2 = x2 - Math.cos(perp) * arrowSize;barY2 = y2 - Math.sin(perp) * arrowSize;}}
                            arrowHead.plot(barX1, barY1, barX2, barY2);
                        }}
                    }}
                    const arrow = getArrowById(arrowId);
                    ['1', '2'].forEach(num => {{
                        const attachedKey = `attached_arrow_${{num}}`;
                        const endKey = `attached_end_${{num}}`;
                        if (arrow[attachedKey]) {{
                            const masterId = arrow[attachedKey];
                            const masterEnd = arrow[endKey];
                            const masterElement = draw.find(`[data-id="${{masterId}}"]`)[0];
                            const masterHitbox = draw.find(`[data-id="${{masterId}}_hit"]`)[0];
                            const masterHead = draw.find(`[data-id="${{masterId}}_head"]`)[0];
                            const masterHandle = arrowHandleGroups[masterId]?.[`${{masterEnd}}Handle`];
                            let mx1 = parseFloat(masterElement.attr('x1')) + (masterEnd === 'start' ? deltaX : 0);
                            let my1 = parseFloat(masterElement.attr('y1')) + (masterEnd === 'start' ? deltaY : 0);
                            let mx2 = parseFloat(masterElement.attr('x2')) + (masterEnd === 'end' ? deltaX : 0);
                            let my2 = parseFloat(masterElement.attr('y2')) + (masterEnd === 'end' ? deltaY : 0);
                            masterElement.plot(mx1, my1, mx2, my2);
                            if (masterHitbox) masterHitbox.plot(mx1, my1, mx2, my2);
                            if (masterHandle) masterHandle.cx(masterEnd === 'start' ? mx1 : mx2).cy(masterEnd === 'start' ? my1 : my2);
                            if (masterHead) {{
                                const dx = mx2 - mx1, dy = my2 - my1, angle = Math.atan2(dy, dx), arrowSize = 5;
                                const lineType = masterElement.attr('data-line-type') || 'arrow';
                                if (lineType === 'arrow') {{
                                    masterHead.plot([[mx2, my2],[mx2 - arrowSize * Math.cos(angle + Math.PI / 6), my2 - arrowSize * Math.sin(angle + Math.PI / 6)],[mx2 - arrowSize * Math.cos(angle - Math.PI / 6), my2 - arrowSize * Math.sin(angle - Math.PI / 6)]]);
                                }} else if (lineType === 'inhibition') {{
                                    const masterArrow = getArrowById(masterId);
                                    const endSide = masterArrow?.protbox_id_2_side;
                                    let barX1, barY1, barX2, barY2;
                                    if (endSide === 'North' || endSide === 'South') {{barX1 = mx2 - arrowSize;barY1 = my2;barX2 = mx2 + arrowSize;barY2 = my2;}} 
                                    else if (endSide === 'West' || endSide === 'East') {{barX1 = mx2;barY1 = my2 - arrowSize;barX2 = mx2;barY2 = my2 + arrowSize;}} 
                                    else {{const perp = angle + Math.PI / 2;barX1 = mx2 + Math.cos(perp) * arrowSize;barY1 = my2 + Math.sin(perp) * arrowSize;barX2 = mx2 - Math.cos(perp) * arrowSize;barY2 = my2 - Math.sin(perp) * arrowSize;}}
                                    masterHead.plot(barX1, barY1, barX2, barY2);
                                }}
                            }}
                            updateAttachedArrows(masterId, masterEnd);
                            Shiny?.setInputValue('arrow_moved', {{ id: masterId, end: masterEnd, x: masterEnd === 'start' ? mx1 : mx2, y: masterEnd === 'start' ? my1 : my2 }}, {{ priority: 'event' }});
                        }}
                    }});
                    updateAttachedArrows(arrowId, 'start');
                    updateAttachedArrows(arrowId, 'end');
                }} else if (isCircle) {{
                    // compute the normal target from the delta
                    let targetX = newX;
                    let targetY = newY;
                    
                    // Live snap using pointer coords (if available)
                    if (pointerX != null && pointerY != null) {{
                      const snaps = protboxSnapPoints[selectedProtboxId];
                      if (snaps) {{
                        // find nearest snap key
                        let nearestKey = null;
                        let nearestDist = Infinity;
                        for (const key in snaps) {{
                          const sp = snaps[key];
                          const d = Math.hypot(pointerX - sp.x, pointerY - sp.y);
                          if (d < nearestDist) {{ nearestDist = d; nearestKey = key; }}
                        }}
                    
                        // release snap if pointer leaves the active point
                        if (activeSnapKey) {{
                          const sp = snaps[activeSnapKey];
                          if (!sp || Math.hypot(pointerX - sp.x, pointerY - sp.y) > ptmSnapRadius) {{
                            activeSnapKey = null;
                          }}
                        }}
                    
                        // acquire snap if close enough
                        if (!activeSnapKey && nearestKey && nearestDist <= ptmSnapRadius) {{
                          activeSnapKey = nearestKey;
                        }}
                    
                        // apply snap position
                        if (activeSnapKey) {{
                          const sp = snaps[activeSnapKey];
                          targetX = sp.x;
                          targetY = sp.y;
                        }}
                      }}
                    }}  
                    // Move the PTM to the final (possibly snapped) spot
                    const appliedDX = targetX - currentX;
                    const appliedDY = targetY - currentY;
                    selectedElement.cx(targetX).cy(targetY);
                    if (activeSnapKey) {{
                        selectedElement.attr('data-pos-key', activeSnapKey);
                    }} else {{
                        selectedElement.attr('data-pos-key', null);
                        try {{ selectedElement.node.removeAttribute('data-pos-key'); }} catch (ignore) {{}}
                    }}
                    
                    // Make the labels/symbols follow the *actual* applied movement (including snap)
                    deltaX = appliedDX;
                    deltaY = appliedDY;
                    }} else {{
                        selectedElement.move(newX, newY);
                    }}
                if (selectedType === 'prot-box' && selectedProtboxId && elementGroups[selectedProtboxId]) {{
                    elementGroups[selectedProtboxId].forEach(assocElement => {{
                        const assocType = assocElement.element.attr('data-type');
                        if (assocType === 'prot-box') return;  // ← don't re-move the rectangle
                        if (assocType === 'handle-line') {{
                            assocElement.element.plot(parseFloat(assocElement.element.attr('x1') || 0) + deltaX,parseFloat(assocElement.element.attr('y1') || 0) + deltaY,parseFloat(assocElement.element.attr('x2') || 0) + deltaX,parseFloat(assocElement.element.attr('y2') || 0) + deltaY);
                        }} else {{
                        assocElement.element.dmove(deltaX, deltaY);
                        }}
                    }});
                    if (selectedProtboxId in attachments) {{
                        Object.keys(attachments[selectedProtboxId]).forEach(side => updateArrowPositions(selectedProtboxId, side));
                    }}
                    const spacing = settings.ptm_circle_spacing || 4;
                    const width = selectedElement.width();
                    const height = selectedElement.height();
                    const pbX = selectedElement.x();
                    const pbY = selectedElement.y();
                    protboxSnapPoints[selectedProtboxId] = {{'N1': {{x: pbX + width * 0.2, y: pbY - spacing}},'N2': {{x: pbX + width * 0.5, y: pbY - spacing}},'N3': {{x: pbX + width * 0.8, y: pbY - spacing}},'S1': {{x: pbX + width * 0.2, y: pbY + height + spacing}},'S2': {{x: pbX + width * 0.5, y: pbY + height + spacing}},'S3': {{x: pbX + width * 0.8, y: pbY + height + spacing}},'W1': {{x: pbX - spacing, y: pbY + height * 0.33 - 2}},'W2': {{x: pbX - spacing, y: pbY + height * 0.66 + 2}},'E1': {{x: pbX + width + spacing, y: pbY + height * 0.33 - 2}},'E2': {{x: pbX + width + spacing, y: pbY + height * 0.66 + 2}}}};
                    elementGroups[selectedProtboxId].filter(el => el.element.attr('data-type') === 'ptm-snap-circle').forEach(snapC => {{
                        const key = snapC.element.attr('data-pos-key');
                        snapC.element.cx(protboxSnapPoints[selectedProtboxId][key].x).cy(protboxSnapPoints[selectedProtboxId][key].y);
                    }});
                    protboxMap[selectedProtboxId].x = pbX;
                    protboxMap[selectedProtboxId].y = pbY;
                }}
                if (selectedType === 'ptm-shape') {{
                    const baseId = selectedId.replace('_shape', '');
                    const labelElement = draw.find(`[data-id="${{baseId}}_label"]`)[0];
                    const symbolElement = draw.find(`[data-id="${{baseId}}_symbol"]`)[0];
                    if (labelElement) labelElement.move(parseFloat(labelElement.x() || 0) + deltaX, parseFloat(labelElement.y() || 0) + deltaY);
                    if (symbolElement) symbolElement.move(parseFloat(symbolElement.x() || 0) + deltaX, parseFloat(symbolElement.y() || 0) + deltaY);
                }}
                if (selectedType === 'arrow-start' || selectedType === 'arrow-end') {{
                    const context = arrowHandleContext;
                    const arrowId = context?.arrowId || selectedId.replace(/_(start|end)$/, '');
                    const arrowElement = context?.arrowElement || draw.find(`[data-id="${{arrowId}}"]`)[0];
                    const arrowHitbox = context?.arrowHitbox || draw.find(`[data-id="${{arrowId}}_hit"]`)[0];
                    const arrowHead = context?.arrowHead || draw.find(`[data-id="${{arrowId}}_head"]`)[0];
                    if (arrowElement) {{
                        let x1 = parseFloat(arrowElement.attr('x1') || 0);
                        let y1 = parseFloat(arrowElement.attr('y1') || 0);
                        let x2 = parseFloat(arrowElement.attr('x2') || 0);
                        let y2 = parseFloat(arrowElement.attr('y2') || 0);
                        if (selectedType === 'arrow-start') {{ x1 = newX; y1 = newY; }} else {{ x2 = newX; y2 = newY; }}
                        arrowElement.plot(x1, y1, x2, y2);
                        if (arrowHitbox) arrowHitbox.plot(x1, y1, x2, y2);
                        if (arrowHead) {{
                            const dx = x2 - x1, dy = y2 - y1, angle = Math.atan2(dy, dx), arrowSize = 5;
                            const lineType = arrowElement.attr('data-line-type') || 'arrow';
                            if (lineType === 'arrow') {{
                                arrowHead.plot([[x2, y2],[x2 - arrowSize * Math.cos(angle + Math.PI / 6), y2 - arrowSize * Math.sin(angle + Math.PI / 6)],[x2 - arrowSize * Math.cos(angle - Math.PI / 6), y2 - arrowSize * Math.sin(angle - Math.PI / 6)]]);
                            }} else if (lineType === 'inhibition') {{
                                const arrow = context?.arrowData || getArrowById(arrowId);
                                const endSide = arrow?.protbox_id_2_side;
                                let barX1, barY1, barX2, barY2;
                                if (endSide === 'North' || endSide === 'South') {{barX1 = x2 - arrowSize;barY1 = y2;barX2 = x2 + arrowSize;barY2 = y2;}} 
                                else if (endSide === 'West' || endSide === 'East') {{barX1 = x2;barY1 = y2 - arrowSize;barX2 = x2;barY2 = y2 + arrowSize;}} 
                                else {{const perp = angle + Math.PI / 2;barX1 = x2 + Math.cos(perp) * arrowSize;barY1 = y2 + Math.sin(perp) * arrowSize;barX2 = x2 - Math.cos(perp) * arrowSize;barY2 = y2 - Math.sin(perp) * arrowSize;}}
                                arrowHead.plot(barX1, barY1, barX2, barY2);
                            }}
                        }}
                    }}
                }}
            }};
            const makeDraggable = (element, type, id, protboxId = null) => {{
                let isDragging = false;
                let prevPointX, prevPointY;
                element.node.style.cursor = 'pointer';
                element.node.addEventListener('mousedown', e => {{
                    e.preventDefault();
                    e.stopPropagation();
                    isDragging = true;
                    selectElement(element, type, id, protboxId);
                    const point = draw.point(e.clientX, e.clientY);
                    prevPointX = point.x;
                    prevPointY = point.y;
                    if (mkHistory && trackableMoveTypes.has(type)) {{
                        mkHistory.beginMoveSession({{
                            element,
                            type,
                            id,
                            protboxId: protboxId || (type === 'prot-box' ? id : (type.startsWith('ptm-') ? selectedProtboxId : null))
                        }});
                    }}
                    if (type === 'arrow') {{
                        const arrowId = id;
                        const arrow = getArrowById(arrowId);
                        if (arrow) {{
                            if (attachedBy[arrowId]) {{
                                ['start', 'end'].forEach(endType => {{
                                    if (attachedBy[arrowId][endType]) {{
                                        attachedBy[arrowId][endType].forEach(slaveInfo => {{
                                            const sArrow = getArrowById(slaveInfo.slave);
                                            const sNum = slaveInfo.slaveEnd === 'start' ? 1 : 2;
                                            delete sArrow[`attached_arrow_${{sNum}}`];
                                            delete sArrow[`attached_end_${{sNum}}`];
                                        }});
                                    }}
                                }});
                                delete attachedBy[arrowId];
                            }}
                            ['1', '2'].forEach(num => {{
                                const protboxId = arrow[`protbox_id_${{num}}`];
                                const side = arrow[`protbox_id_${{num}}_side`];
                                if (protboxId && side) {{
                                    attachments[protboxId][side] = attachments[protboxId][side].filter(att => att.arrow !== arrow);
                                    if (attachments[protboxId][side].length === 0) delete attachments[protboxId][side];
                                    if (Object.keys(attachments[protboxId]).length === 0) delete attachments[protboxId];
                                    delete arrow[`protbox_id_${{num}}`];
                                    delete arrow[`protbox_id_${{num}}_side`];
                                }}
                            }});
                        }}
                    }}
                }});
                const dragging = e => {{
                    if (!isDragging) return;
                    e.preventDefault();
                    e.stopPropagation();
                    const point = draw.point(e.clientX, e.clientY);
                    const deltaX = point.x - prevPointX;
                    const deltaY = point.y - prevPointY;
                    prevPointX = point.x;
                    prevPointY = point.y;
                    moveSelectedElement(deltaX, deltaY, point.x, point.y);
                }};
                const endDragging = e => {{
                    if (!isDragging) return;
                    isDragging = false;
                    if (mkHistory && trackableMoveTypes.has(type)) {{
                        mkHistory.finalizeMoveSession();
                    }}
                    let currentX, currentY;
                    if (selectedType === 'arrow-start' || selectedType === 'arrow-end') {{
                        const arrowId = selectedId.replace(/_(start|end)$/, '');
                        const end = selectedType.replace('arrow-', '');
                        const handle = selectedType === 'arrow-start' ? arrowHandleGroups[arrowId]?.startHandle : arrowHandleGroups[arrowId]?.endHandle;
                        const snapped = updateArrowAttachments(arrowId, end, handle.cx(), handle.cy());
                        if (snapped.isSnapped) {{
                            if (snapped.protbox_id) {{
                                arrow.protbox_id = snapped.protbox_id;
                                arrow.side = snapped.side;
                            }} else if (snapped.arrow_id) {{
                                arrow.attached_arrow = snapped.arrow_id;
                                arrow.attached_end = snapped.arrow_end;
                            }}
                        }}
                        updateAttachedArrows(arrowId, end);
                        if (handle) {{
                            Shiny?.setInputValue('arrow_moved', {{ id: arrowId, end, x: handle.cx(), y: handle.cy() }}, {{ priority: 'event' }});
                        }}
                        const arrow = getArrowById(arrowId);
                        if (arrow && arrow.line === 'inhibition') {{
                            const endSide = end === 'end' ? arrow.protbox_id_2_side : arrow.protbox_id_1_side;
                            const protboxId = end === 'end' ? arrow.protbox_id_2 : arrow.protbox_id_1;
                            if (endSide && protboxId) {{
                                updateArrowPositions(protboxId, endSide);
                            }}
                        }}
                    }} else if (selectedType === 'arrow') {{
                        const arrowId = selectedId;
                        const startHandle = arrowHandleGroups[arrowId]?.startHandle;
                        const endHandle = arrowHandleGroups[arrowId]?.endHandle;
                        if (startHandle && endHandle) {{
                            Shiny?.setInputValue('arrow_moved', {{ id: arrowId, end: 'both', x1: startHandle.cx(), y1: startHandle.cy(), x2: endHandle.cx(), y2: endHandle.cy() }}, {{ priority: 'event' }});
                        }}
                }} else {{
                    const isCircle = (selectedType === 'ptm-shape' && selectedElement && selectedElement.type === 'circle');
                    const rawX = selectedElement ? (isCircle ? selectedElement.cx() : selectedElement.x()) : 0;
                    const rawY = selectedElement ? (isCircle ? selectedElement.cy() : selectedElement.y()) : 0;
                    currentX = parseFloat(rawX || 0);
                    currentY = parseFloat(rawY || 0);
                    let reportedY = currentY;
                    if (selectedType === 'ptm-label' || selectedType === 'ptm-symbol') {{
                        reportedY = currentY - textOffsetY;
                    }}
                    const eventPayload = {{ type: selectedType, id: selectedId, x: currentX, y: reportedY }};
                    if (selectedProtboxId) {{
                        eventPayload.protbox_id = selectedProtboxId;
                    }}
                    if (selectedType === 'ptm-shape') {{
                        const posKeyAttr = selectedElement?.attr ? selectedElement.attr('data-pos-key') : null;
                        if (posKeyAttr) {{
                            eventPayload.ptm_position = posKeyAttr;
                        }}
                    }}
                    if (selectedType === 'ptm-label') {{
                        eventPayload.label_centering = resolveLabelCenteringFromElement(selectedElement);
                    }}
                    Shiny?.setInputValue('element_moved', eventPayload, {{ priority: 'event' }});
                    if (selectedType === 'ptm-shape') {{
                        updateHandleDists(selectedProtboxId);
                        if (selectedProtboxId in attachments) {{
                            Object.keys(attachments[selectedProtboxId]).forEach(side => updateArrowPositions(selectedProtboxId, side));
                        }}
                    }}
                    if (selectedType === 'ptm-shape' || selectedType === 'ptm-label' || selectedType === 'ptm-symbol') {{
                        const meta = parsePtmElementId(selectedId);
                        if (meta && selectedProtboxId) {{
                            const overridePayload = {{}};
                            if (selectedType === 'ptm-shape') {{
                                overridePayload.shape_x = currentX;
                                overridePayload.shape_y = currentY;
                                overridePayload.ptm_position = eventPayload.ptm_position || null;
                            }} else if (selectedType === 'ptm-label') {{
                                overridePayload.label_x = currentX;
                                overridePayload.label_y = reportedY;
                                overridePayload.label_centering = eventPayload.label_centering || resolveLabelCenteringFromElement(selectedElement);
                            }} else if (selectedType === 'ptm-symbol') {{
                                overridePayload.symbol_x = currentX;
                                overridePayload.symbol_y = reportedY;
                            }}
                            recordPtmOverride(selectedProtboxId, meta.uniprot, meta.ptmKey, overridePayload);
                        }}
                    }}
                    activeSnapKey = null;
                }}
                }};
                document.addEventListener('mousemove', dragging);
                document.addEventListener('mouseup', endDragging);
            }};
            const removeElementByDataId = (elementId) => {{
                if (!elementId) return;
                try {{
                    const node = draw.find(`[data-id="${{elementId}}"]`)[0];
                    if (node) {{
                        node.remove();
                    }}
                }} catch (err) {{}}
            }};
            const cleanupProtboxAttachmentsForArrow = (arrowId) => {{
                Object.keys(attachments).forEach(pbId => {{
                    const perSide = attachments[pbId];
                    if (!perSide) return;
                    Object.keys(perSide).forEach(side => {{
                        perSide[side] = perSide[side].filter(entry => entry && entry.arrowId !== arrowId);
                        if (!perSide[side].length) delete perSide[side];
                    }});
                    if (!Object.keys(perSide).length) delete attachments[pbId];
                }});
            }};
            const ensureArrowHandle = (arrowId, whichEnd) => {{
                const handleKey = whichEnd === 'start' ? 'startHandle' : 'endHandle';
                arrowHandleGroups[arrowId] = arrowHandleGroups[arrowId] || {{}};
                if (arrowHandleGroups[arrowId][handleKey]) return;
                const line = draw.find(`[data-id="${{arrowId}}"]`)[0];
                if (!line) return;
                const cx = parseFloat(line.attr(whichEnd === 'start' ? 'x1' : 'x2') || 0);
                const cy = parseFloat(line.attr(whichEnd === 'start' ? 'y1' : 'y2') || 0);
                const handleId = `${{arrowId}}_${{whichEnd}}`;
                const handle = arrowGroup.circle(6).cx(cx).cy(cy).fill('red').stroke({{ color: 'black', width: 1 }}).attr({{ 'data-id': handleId, 'data-type': whichEnd === 'start' ? 'arrow-start' : 'arrow-end' }}).hide();
                arrowHandleGroups[arrowId][handleKey] = handle;
                makeDraggable(handle, whichEnd === 'start' ? 'arrow-start' : 'arrow-end', handleId);
            }};
            const detachDependentArrows = (arrowId) => {{
                const dependents = attachedBy[arrowId];
                if (!dependents) return;
                ['start', 'end'].forEach(endType => {{
                    (dependents[endType] || []).forEach(slaveInfo => {{
                        const slaveId = slaveInfo.slave;
                        const slaveArrow = getArrowById(slaveId);
                        if (!slaveArrow) return;
                        const slot = slaveInfo.slaveEnd === 'start' ? '1' : '2';
                        delete slaveArrow[`attached_arrow_${{slot}}`];
                        delete slaveArrow[`attached_end_${{slot}}`];
                        ensureArrowHandle(slaveId, slaveInfo.slaveEnd);
                    }});
                }});
            }};
            const pruneAttachedByReferences = (arrowId) => {{
                Object.keys(attachedBy).forEach(masterId => {{
                    ['start', 'end'].forEach(endType => {{
                        if (attachedBy[masterId]?.[endType]) {{
                            attachedBy[masterId][endType] = attachedBy[masterId][endType].filter(s => s.slave !== arrowId);
                            if (!attachedBy[masterId][endType].length) delete attachedBy[masterId][endType];
                        }}
                    }});
                    if (attachedBy[masterId] && !Object.keys(attachedBy[masterId]).length) {{
                        delete attachedBy[masterId];
                    }}
                }});
            }};
            const detachArrowsFromProtbox = (targetId) => {{
                if (!targetId) return [];
                const detachedArrowPayloads = [];
                const processed = new Set();
                const detachArrowEnd = (arrowId, arrow, endType, arrowIndex, sideHint = null) => {{
                    if (!arrowId || !arrow) return;
                    const key = `${{arrowId}}_${{endType}}`;
                    if (processed.has(key)) return;
                    processed.add(key);
                    const coordKeyX = endType === 'start' ? 'x1' : 'x2';
                    const coordKeyY = endType === 'start' ? 'y1' : 'y2';
                    const pbIndex = endType === 'start' ? '1' : '2';
                    const pbField = `protbox_id_${{pbIndex}}`;
                    const sideField = `${{pbField}}_side`;
                    const prevSide = arrow[sideField] || sideHint || null;
                    if (normalizeProtboxId(arrow[pbField]) !== targetId) {{
                        delete arrow[pbField];
                        delete arrow[sideField];
                    }}
                    const arrowElement = draw.find(`[data-id="${{arrowId}}"]`)[0];
                    let coordX = Number(arrow[coordKeyX]) || 0;
                    let coordY = Number(arrow[coordKeyY]) || 0;
                    if (arrowElement) {{
                        const attrX = arrowElement.attr(endType === 'start' ? 'x1' : 'x2');
                        const attrY = arrowElement.attr(endType === 'start' ? 'y1' : 'y2');
                        const attrXNum = Number(attrX);
                        const attrYNum = Number(attrY);
                        if (Number.isFinite(attrXNum)) coordX = attrXNum;
                        if (Number.isFinite(attrYNum)) coordY = attrYNum;
                    }}
                    arrow[coordKeyX] = coordX;
                    arrow[coordKeyY] = coordY;
                    delete arrow[pbField];
                    delete arrow[sideField];
                    ensureArrowHandle(arrowId, endType);
                    const handleKey = endType === 'start' ? 'startHandle' : 'endHandle';
                    const handle = arrowHandleGroups[arrowId]?.[handleKey];
                    if (handle) handle.cx(coordX).cy(coordY);
                    detachedArrowPayloads.push({{
                        arrow_id: arrowId,
                        arrow_index: arrowIndex,
                        end: endType,
                        x: coordX,
                        y: coordY,
                        protbox_id: targetId,
                        side: prevSide
                    }});
                }};
                const targetAttachments = attachments[targetId];
                if (targetAttachments) {{
                    Object.keys(targetAttachments).forEach(side => {{
                        const entries = targetAttachments[side] || [];
                        entries.forEach(entry => {{
                            const {{ arrow, arrowId, type }} = entry || {{}};
                            if (!arrow || !arrowId) return;
                            const arrowIndex = parseArrowIndex(arrowId);
                            const endType = type === 'end' ? 'end' : 'start';
                            detachArrowEnd(arrowId, arrow, endType, arrowIndex, side);
                        }});
                    }});
                    delete attachments[targetId];
                }}
                try {{
                    arrows.forEach((arrow, index) => {{
                        if (!arrow) return;
                        const arrowId = arrowIdFromIndex(index);
                        if (!arrowId) return;
                        ['start', 'end'].forEach(endType => {{
                            const pbIndex = endType === 'start' ? '1' : '2';
                            const pbField = `protbox_id_${{pbIndex}}`;
                            const sideField = `${{pbField}}_side`;
                            if (normalizeProtboxId(arrow[pbField]) === targetId) {{
                                detachArrowEnd(arrowId, arrow, endType, index, arrow[sideField]);
                            }}
                        }});
                    }});
                }} catch (err) {{
                    console.warn('Unable to detach arrows for protbox', targetId, err);
                }}
                return detachedArrowPayloads;
            }};
            const reattachArrowPayloads = (payloads) => {{
                if (!Array.isArray(payloads) || !payloads.length) return;
                payloads.forEach(payload => {{
                    if (!payload || !payload.arrow_id || !payload.protbox_id) return;
                    const arrow = getArrowById(payload.arrow_id);
                    if (!arrow) return;
                    const num = payload.end === 'start' ? '1' : '2';
                    const coordKeyX = payload.end === 'start' ? 'x1' : 'x2';
                    const coordKeyY = payload.end === 'start' ? 'y1' : 'y2';
                    arrow[`protbox_id_${{num}}`] = payload.protbox_id;
                    if (payload.side) {{
                        arrow[`protbox_id_${{num}}_side`] = payload.side;
                    }}
                    arrow[coordKeyX] = payload.x;
                    arrow[coordKeyY] = payload.y;
                    const arrowElement = draw.find(`[data-id="${{payload.arrow_id}}"]`)[0];
                    const arrowHitbox = draw.find(`[data-id="${{payload.arrow_id}}_hit"]`)[0];
                    if (arrowElement) {{
                        if (payload.end === 'start') {{
                            arrowElement.attr({{ x1: payload.x, y1: payload.y }});
                        }} else {{
                            arrowElement.attr({{ x2: payload.x, y2: payload.y }});
                        }}
                    }}
                    if (arrowHitbox && arrowElement) {{
                        arrowHitbox.plot(arrowElement.attr('x1'), arrowElement.attr('y1'), arrowElement.attr('x2'), arrowElement.attr('y2'));
                    }}
                    ensureArrowHandle(payload.arrow_id, payload.end);
                    const handleKey = payload.end === 'start' ? 'startHandle' : 'endHandle';
                    const handle = arrowHandleGroups[payload.arrow_id]?.[handleKey];
                    if (handle) handle.cx(payload.x).cy(payload.y);
                    attachments[payload.protbox_id] = attachments[payload.protbox_id] || {{}};
                    if (payload.side) {{
                        attachments[payload.protbox_id][payload.side] = attachments[payload.protbox_id][payload.side] || [];
                        attachments[payload.protbox_id][payload.side].push({{
                            type: payload.end === 'end' ? 'end' : 'start',
                            arrow,
                            arrowId: payload.arrow_id,
                            side: payload.side
                        }});
                        updateArrowPositions(payload.protbox_id, payload.side);
                    }}
                }});
            }};
            const deleteArrowById = (arrowId, options = {{}}) => {{
                if (!arrowId) return false;
                const idx = parseArrowIndex(arrowId);
                if (idx === null || idx < 0 || idx >= arrows.length) return false;
                const arrow = arrows[idx];
                if (!arrow) return false;
                detachDependentArrows(arrowId);
                cleanupProtboxAttachmentsForArrow(arrowId);
                pruneAttachedByReferences(arrowId);
                removeElementByDataId(`${{arrowId}}_head`);
                removeElementByDataId(`${{arrowId}}_hit`);
                removeElementByDataId(`${{arrowId}}`);
                removeElementByDataId(`${{arrowId}}_start`);
                removeElementByDataId(`${{arrowId}}_end`);
                if (arrowHandleGroups[arrowId]) {{
                    Object.values(arrowHandleGroups[arrowId]).forEach(handle => {{
                        try {{ handle?.remove(); }} catch (err) {{}}
                    }});
                    delete arrowHandleGroups[arrowId];
                }}
                delete attachedBy[arrowId];
                arrows[idx] = null;
                if (!options.silent) {{
                    Shiny?.setInputValue('delete_element', {{ type: 'arrow', id: arrowId, arrow_index: idx, timestamp: Date.now() }}, {{ priority: 'event' }});
                }}
                return true;
            }};
            const deleteCompoundById = (compoundId, options = {{}}) => {{
                if (!compoundId) return false;
                removeElementByDataId(compoundId);
                removeElementByDataId(`${{compoundId}}_label`);
                let payloadId = null;
                const idx = compounds.findIndex(c => c && (c._client_id === compoundId || `compound_${{c.compound_id}}` === compoundId));
                const snapshot = idx >= 0 ? cloneData(compounds[idx]) : null;
                if (snapshot && !snapshot._client_id) {{
                    snapshot._client_id = compoundId;
                }}
                if (idx >= 0) {{
                    payloadId = compounds[idx]?.compound_id ?? null;
                    compounds.splice(idx, 1);
                }}
                Shiny?.setInputValue('delete_element', {{ type: 'compound', id: compoundId, compound_id: payloadId, timestamp: Date.now() }}, {{ priority: 'event' }});
                if (mkHistory && options.suppressHistory !== true && snapshot) {{
                    const entry = {{
                        kind: 'compound-delete',
                        snapshot,
                        index: idx
                    }};
                    entry.handlers = {{
                        undo: () => mkHistory.runWithoutRecording(() => {{
                            restoreCompoundSnapshot(snapshot, entry.index);
                        }}),
                        redo: () => mkHistory.runWithoutRecording(() => {{
                            deleteCompoundById(snapshot._client_id || compoundId, {{ suppressHistory: true }});
                        }})
                    }};
                    mkHistory.recordAction(entry);
                }}
                return true;
            }};
            const deleteTextBoxById = (textId, options = {{}}) => {{
                if (!textId) return false;
                removeElementByDataId(textId);
                removeElementByDataId(`${{textId}}_label`);
                let payloadId = null;
                const idx = textBlocks.findIndex(tb => tb && (tb._client_id === textId || `text_${{tb.text_id}}` === textId));
                const snapshot = idx >= 0 ? cloneData(textBlocks[idx]) : null;
                if (snapshot && !snapshot._client_id) {{
                    snapshot._client_id = textId;
                }}
                if (idx >= 0) {{
                    payloadId = textBlocks[idx]?.text_id ?? null;
                    textBlocks.splice(idx, 1);
                }}
                Shiny?.setInputValue('delete_element', {{ type: 'text-box', id: textId, text_id: payloadId, timestamp: Date.now() }}, {{ priority: 'event' }});
                if (mkHistory && options.suppressHistory !== true && snapshot) {{
                    const entry = {{
                        kind: 'text-delete',
                        snapshot,
                        index: idx
                    }};
                    entry.handlers = {{
                        undo: () => mkHistory.runWithoutRecording(() => {{
                            restoreTextSnapshot(snapshot, entry.index);
                        }}),
                        redo: () => mkHistory.runWithoutRecording(() => {{
                            deleteTextBoxById(snapshot._client_id || textId, {{ suppressHistory: true }});
                        }})
                    }};
                    mkHistory.recordAction(entry);
                }}
                return true;
            }};
            const deletePtmFromProtbox = (protboxId, meta, options = {{}}) => {{
                if (!protboxId || !meta) return false;
                const targetId = normalizeProtboxId(protboxId);
                if (!targetId) return false;
                const entries = elementGroups[targetId];
                if (!entries) return false;
                const snapshot = capturePtmSnapshot(targetId, meta);
                const shapeId = `${{meta.uniprot}}_${{meta.ptmKey}}_shape`;
                const labelId = `${{meta.uniprot}}_${{meta.ptmKey}}_label`;
                const symbolId = `${{meta.uniprot}}_${{meta.ptmKey}}_symbol`;
                [shapeId, labelId, symbolId].forEach(id => removeElementByDataId(id));
                elementGroups[targetId] = entries.filter(entry => {{
                    if (!entry?.element?.attr) return true;
                    const entryId = entry.element.attr('data-id');
                    return entryId !== shapeId && entryId !== labelId && entryId !== symbolId;
                }});
                recordPtmOverride(targetId, meta.uniprot, meta.ptmKey, {{ hidden: true }});
                updateHandleDists(targetId);
                if (attachments[targetId]) {{
                    Object.keys(attachments[targetId]).forEach(side => updateArrowPositions(targetId, side));
                }}
                Shiny?.setInputValue('delete_element', {{
                    type: 'ptm',
                    protbox_id: targetId,
                    uniprot: meta.uniprot,
                    ptm_key: meta.ptmKey,
                    timestamp: Date.now()
                }}, {{ priority: 'event' }});
                if (mkHistory && options.suppressHistory !== true && snapshot) {{
                    const entry = {{
                        kind: 'ptm-delete',
                        snapshot
                    }};
                    entry.handlers = {{
                        undo: () => mkHistory.runWithoutRecording(() => {{
                            spawnPtmForProtbox(targetId, meta.uniprot, meta.ptmKey, {{
                                recordHistory: false,
                                snapshot
                            }});
                        }}),
                        redo: () => mkHistory.runWithoutRecording(() => {{
                            deletePtmFromProtbox(targetId, meta, {{ suppressHistory: true }});
                        }})
                    }};
                    mkHistory.recordAction(entry);
                }}
                return true;
            }};
            const deleteProtboxById = (protboxId, options = {{}}) => {{
                const targetId = normalizeProtboxId(protboxId);
                if (!targetId) return false;
                const idx = protBoxes.findIndex(pb => normalizeProtboxId(pb?.protbox_id) === targetId);
                if (idx === -1) return false;
                const snapshot = cloneData(protBoxes[idx]);
                const detachedArrowPayloads = detachArrowsFromProtbox(targetId);
                (elementGroups[targetId] || []).forEach(entry => {{
                    try {{ entry.element?.remove(); }} catch (err) {{}}
                }});
                delete elementGroups[targetId];
                removeElementByDataId(targetId);
                removeElementByDataId(`${{targetId}}_label`);
                delete protboxMap[targetId];
                delete protboxSnapPoints[targetId];
                delete protboxHandleDists[targetId];
                delete currentSelected[targetId];
                protBoxes.splice(idx, 1);
                if (Array.isArray(groups)) {{
                    groups.forEach(group => {{
                        if (!group?.protbox_ids) return;
                        group.protbox_ids = group.protbox_ids.filter(id => normalizeProtboxId(id) !== targetId);
                    }});
                }}
                Shiny?.setInputValue('delete_element', {{
                    type: 'prot-box',
                    id: targetId,
                    protbox_id: targetId,
                    detached_arrows: detachedArrowPayloads,
                    timestamp: Date.now()
                }}, {{ priority: 'event' }});
                if (mkHistory && options.suppressHistory !== true && snapshot) {{
                    const entry = {{
                        kind: 'protbox-delete',
                        protbox: snapshot,
                        index: idx,
                        detached_arrows: cloneData(detachedArrowPayloads)
                    }};
                    entry.handlers = {{
                        undo: () => mkHistory.runWithoutRecording(() => {{
                            if (restoreProtboxSnapshot(entry.protbox, entry.index)) {{
                                reattachArrowPayloads(entry.detached_arrows);
                            }}
                        }}),
                        redo: () => mkHistory.runWithoutRecording(() => {{
                            deleteProtboxById(entry.protbox.protbox_id, {{ suppressHistory: true }});
                        }})
                    }};
                    mkHistory.recordAction(entry);
                }}
                return true;
            }};
            const isEditableTarget = (target) => {{
                if (!target) return false;
                const tag = target.tagName ? target.tagName.toLowerCase() : '';
                if (['input', 'textarea', 'select'].includes(tag)) return true;
                return Boolean(target.isContentEditable);
            }};
            const deleteSelectedElement = () => {{
                if (!selectedElement || !selectedType || !selectedId) return;
                const targetId = selectedId;
                const targetType = selectedType;
                const targetProtboxId = selectedProtboxId;
                if (targetType === 'prot-box') {{
                    const deletedProt = deleteProtboxById(targetId);
                    if (deletedProt) {{
                        deselectElement();
                    }}
                    return;
                }}
                deselectElement();
                let deleted = false;
                if (targetType === 'arrow' || targetType === 'arrow-start' || targetType === 'arrow-end') {{
                    const arrowId = targetType === 'arrow' ? targetId : targetId.replace(/_(start|end)$/, '');
                    deleted = deleteArrowById(arrowId);
                }} else if (targetType === 'compound') {{
                    deleted = deleteCompoundById(targetId);
                }} else if (targetType === 'text-box') {{
                    deleted = deleteTextBoxById(targetId);
                }} else if (targetType === 'ptm-shape' || targetType === 'ptm-label' || targetType === 'ptm-symbol') {{
                    const meta = parsePtmElementId(targetId);
                    if (meta && targetProtboxId) {{
                        deleted = deletePtmFromProtbox(targetProtboxId, meta);
                    }}
                }}
            }};
            document.addEventListener('keydown', e => {{
                const targetEditable = isEditableTarget(e.target);
                const keyLower = typeof e.key === 'string' ? e.key.toLowerCase() : '';
                if ((e.ctrlKey || e.metaKey) && keyLower === 'z') {{
                    if (targetEditable) return;
                    e.preventDefault();
                    if (e.shiftKey) {{
                        mkHistory?.redo();
                    }} else {{
                        mkHistory?.undo();
                    }}
                    return;
                }}
                if ((e.ctrlKey || e.metaKey) && (keyLower === 'y')) {{
                    if (targetEditable) return;
                    e.preventDefault();
                    mkHistory?.redo();
                    return;
                }}
                if (!selectedElement) return;
                if (e.key === 'Backspace') {{
                    if (isEditableTarget(e.target)) return;
                    e.preventDefault();
                    deleteSelectedElement();
                    return;
                }}
                let deltaX = 0, deltaY = 0, moveAmount = e.shiftKey ? 10 : 1;
                if (e.key === 'ArrowUp') {{ deltaY = -moveAmount; e.preventDefault(); }}
                else if (e.key === 'ArrowDown') {{ deltaY = moveAmount; e.preventDefault(); }}
                else if (e.key === 'ArrowLeft') {{ deltaX = -moveAmount; e.preventDefault(); }}
                else if (e.key === 'ArrowRight') {{ deltaX = moveAmount; e.preventDefault(); }}
                else if (e.key === 'Escape') {{ deselectElement(); e.preventDefault(); }}
                else return;
                if (deltaX !== 0 || deltaY !== 0) {{
                    if (mkHistory && trackableMoveTypes.has(selectedType)) {{
                        mkHistory.captureInstantMove({{
                            element: selectedElement,
                            type: selectedType,
                            id: selectedId,
                            protboxId: selectedProtboxId
                        }}, () => moveSelectedElement(deltaX, deltaY));
                    }} else {{
                        moveSelectedElement(deltaX, deltaY);
                    }}
                }}
            }});
            draw.node.addEventListener('click', e => {{
                if (isBackgroundTarget(e.target)) {{
                    deselectElement();
                }}
            }});
            draw.node.addEventListener('mousedown', startPanning);
            document.addEventListener('mousemove', pan);
            document.addEventListener('mouseup', endPanning);
            document.getElementById('svgCanvas').addEventListener('wheel', handleWheel, {{ passive: false }});
            document.getElementById('reset-view')?.addEventListener('click', resetView);
            const controlsContainer = document.querySelector('#svgCanvas .canvas-controls');
            if (controlsContainer) {{
                controlsContainer.style.display = '';
            }}
            const undoBtn = document.getElementById('undo-action');
            const redoBtn = document.getElementById('redo-action');
            if (mkHistory) {{
                mkHistory.attachButtons({{ undo: undoBtn, redo: redoBtn }});
            }}
            if (canvas) {{
                const focusCanvas = () => {{
                    try {{
                        canvas.focus({{ preventScroll: true }});
                    }} catch (focusErr) {{
                        try {{ canvas.focus(); }} catch (ignore) {{}}
                    }}
                }};
                if (!canvas.hasAttribute('data-mk-focus-listener')) {{
                    canvas.setAttribute('data-mk-focus-listener', '1');
                    canvas.addEventListener('mousedown', focusCanvas);
                }}
                focusCanvas();
            }}
            const textOffsetY = -14.5;
            function capturePtmSnapshot(protboxId, meta) {{
                if (!protboxId || !meta) return null;
                const baseId = `${{meta.uniprot}}_${{meta.ptmKey}}`;
                const shapeEl = ensureElementForId(`${{baseId}}_shape`);
                if (!shapeEl) return null;
                const shapePos = captureElementPosition(shapeEl, 'ptm-shape');
                if (!shapePos) return null;
                const labelEl = ensureElementForId(`${{baseId}}_label`);
                const symbolEl = ensureElementForId(`${{baseId}}_symbol`);
                const snapshot = {{
                    protbox_id: protboxId,
                    uniprot: meta.uniprot,
                    ptm_key: meta.ptmKey,
                    shape: shapePos
                }};
                if (labelEl) {{
                    const rawX = toCoordinateNumber(labelEl.x && labelEl.x());
                    const rawY = toCoordinateNumber(labelEl.y && labelEl.y());
                    snapshot.label = {{
                        x: rawX,
                        y: rawY - textOffsetY,
                        center: resolveLabelCenteringFromElement(labelEl)
                    }};
                }}
                if (symbolEl) {{
                    const rawX = toCoordinateNumber(symbolEl.x && symbolEl.x());
                    const rawY = toCoordinateNumber(symbolEl.y && symbolEl.y());
                    snapshot.symbol = {{
                        x: rawX,
                        y: rawY - textOffsetY
                    }};
                }}
                return snapshot;
            }}
            const ensureProteinRecord = (uniprot) => {{
                if (!uniprot) return null;
                if (!proteinData[uniprot] && proteinCatalog[uniprot]) {{
                    try {{
                        proteinData[uniprot] = JSON.parse(JSON.stringify(proteinCatalog[uniprot]));
                    }} catch (err) {{
                        console.warn('Failed to clone catalog entry for', uniprot, err);
                        proteinData[uniprot] = proteinCatalog[uniprot];
                    }}
                }}
                return proteinData[uniprot];
            }};

            const spawnPtmForProtbox = (protboxId, uniprot, ptmKey, options = {{}}) => {{
                const protein = ensureProteinRecord(uniprot);
                const ptm = protein?.PTMs?.[ptmKey];
                if (!protein || !ptm) return false;
                const rect = draw.find(`[data-id="${{protboxId}}"]`)[0];
                if (!rect) return false;
                const isPrimaryForProtein = isPrimaryProtbox(protboxId, uniprot);
                const entries = elementGroups[protboxId] || (elementGroups[protboxId] = []);
                const idsToRemove = [
                    `${{uniprot}}_${{ptmKey}}_shape`,
                    `${{uniprot}}_${{ptmKey}}_label`,
                    `${{uniprot}}_${{ptmKey}}_symbol`
                ];
                if (entries.length) {{
                    elementGroups[protboxId] = entries.filter(entry => {{
                        const entryId = entry?.element?.attr ? entry.element.attr('data-id') : null;
                        if (entryId && idsToRemove.includes(entryId)) {{
                            try {{ entry.element.remove(); }} catch (remErr) {{}}
                            return false;
                        }}
                        return true;
                    }});
                }}
                const workingEntries = elementGroups[protboxId] || (elementGroups[protboxId] = []);
                const pbX = rect.x();
                const pbY = rect.y();
                const pbWidth = rect.width();
                const pbHeight = rect.height();
                const spacing = settings.ptm_circle_spacing || 4;
                const snapPoints = protboxSnapPoints[protboxId];
                const preferredSnapKey = options?.preferredSnapKey || null;
                const fixedPosition = options?.fixedPosition || null;
                const forceCenterSpawn = options?.spawnInCenter === true;
                const snapshotSeed = options?.snapshot || null;
                const recordHistory = options?.recordHistory !== false;
                const existingOverride = getPtmOverride(protboxId, uniprot, ptmKey);
                const ptmMeta = {{ uniprot, ptmKey }};
                let snapChoice = null;
                if (!snapshotSeed?.shape) {{
                    if (fixedPosition && fixedPosition.x !== undefined && fixedPosition.y !== undefined) {{
                        snapChoice = {{ x: fixedPosition.x, y: fixedPosition.y, key: fixedPosition.key || null }};
                    }} else if (!preferredSnapKey && existingOverride) {{
                        const overrideX = toFiniteNumber(existingOverride.shape_x);
                        const overrideY = toFiniteNumber(existingOverride.shape_y);
                        if (overrideX !== null && overrideY !== null) {{
                            snapChoice = {{
                                x: overrideX,
                                y: overrideY,
                                key: existingOverride.ptm_position || null
                            }};
                        }}
                    }} else if (!forceCenterSpawn && preferredSnapKey && snapPoints?.[preferredSnapKey]) {{
                        snapChoice = {{ ...snapPoints[preferredSnapKey], key: preferredSnapKey }};
                    }} else if (!forceCenterSpawn && snapPoints) {{
                        const free = pickFreeSnap(protboxId);
                        if (free) snapChoice = free;
                    }}
                }}
                const centerX = pbX + pbWidth / 2;
                const centerY = pbY + pbHeight / 2;
                let newX, newY, posKey = null;
                if (snapshotSeed?.shape && snapshotSeed.shape.x !== undefined && snapshotSeed.shape.y !== undefined) {{
                    newX = snapshotSeed.shape.x;
                    newY = snapshotSeed.shape.y;
                    posKey = snapshotSeed.shape.posKey || null;
                }} else if (snapChoice) {{
                    newX = snapChoice.x;
                    newY = snapChoice.y;
                    posKey = snapChoice.key || preferredSnapKey || existingOverride?.ptm_position || null;
                }} else if (forceCenterSpawn) {{
                    newX = centerX;
                    newY = centerY;
                }} else {{
                    newX = centerX;
                    newY = pbY - spacing;
                }}
                const radius = settings.ptm_circle_radius || 5;
                const shape = (ptm.shape || 'circle').toLowerCase();
                const ptmFill = entityColor(ptm);
                let shapeObj;
                if (shape === 'circle') {{
                    shapeObj = protboxGroup.circle(radius * 2).cx(newX).cy(newY).fill(ptmFill).stroke({{ color: 'black', width: 1 }});
                }} else {{
                    shapeObj = protboxGroup.rect(radius * 2, radius * 2).move(newX - radius, newY - radius).fill(ptmFill).stroke({{ color: 'black', width: 1 }});
                }}
                shapeObj.attr({{
                    'data-id': `${{uniprot}}_${{ptmKey}}_shape`,
                    'data-type': 'ptm-shape',
                    'data-protbox-id': protboxId,
                    'data-tooltip': ptm.tooltip || '',
                    'data-tooltip-html': ptm.tooltip_html || '',
                    ...(posKey ? {{ 'data-pos-key': posKey }} : {{}})
                }});
                shapeObj.node.addEventListener('mouseenter', e => {{
                    const htmlTip = shapeObj.attr('data-tooltip-html') || '';
                    const tip = htmlTip || shapeObj.attr('data-tooltip') || '';
                    if (tip) {{
                        const pos = toContainerPosition(e);
                        setTooltipContent(tip, Boolean(htmlTip));
                        tooltip.style.display = 'block';
                        tooltip.style.left = `${{pos.x + tooltipOffsetX}}px`;
                        tooltip.style.top = `${{pos.y + tooltipOffsetY}}px`;
                    }}
                }});
                shapeObj.node.addEventListener('mousemove', e => {{
                    if (tooltip.style.display === 'block') {{
                        const pos = toContainerPosition(e);
                        tooltip.style.left = `${{pos.x + tooltipOffsetX}}px`;
                        tooltip.style.top = `${{pos.y + tooltipOffsetY}}px`;
                    }}
                }});
                shapeObj.node.addEventListener('mouseleave', () => {{
                    tooltip.style.display = 'none';
                    clearTooltipContent();
                }});
                workingEntries.push({{ element: shapeObj, offsetX: newX - pbX, offsetY: newY - pbY, jsonX: newX, jsonY: newY }});
                makeDraggable(shapeObj, 'ptm-shape', `${{uniprot}}_${{ptmKey}}_shape`, protboxId);
                bindPtmContextMenu(shapeObj, protboxId, ptmMeta);
                const resetLabel = options?.resetLabelPosition || false;
                const resetSymbol = options?.resetSymbolPosition || false;
                const resolveCoord = (value, fallback, forceReset) => (forceReset || value === undefined || value === null) ? fallback : value;
                const forceCenterPlacement = forceCenterSpawn && !existingOverride;
                const labelFallbackY = forceCenterPlacement ? (newY - 10 - textOffsetY) : newY;
                const labelSnapshot = snapshotSeed?.label;
                const symbolSnapshot = snapshotSeed?.symbol;
                const labelCenteringValue = labelSnapshot?.center || existingOverride?.label_centering || ptm.label_centering || 'center';
                const labelX = (labelSnapshot && labelSnapshot.x !== undefined)
                    ? labelSnapshot.x
                    : resolveCoord(existingOverride?.label_x ?? ptm.label_x, newX, resetLabel || forceCenterPlacement);
                const labelY = (labelSnapshot && labelSnapshot.y !== undefined)
                    ? labelSnapshot.y
                    : resolveCoord(existingOverride?.label_y ?? ptm.label_y, labelFallbackY, resetLabel || forceCenterPlacement);
                if (ptm.label) {{
                    const labelColor = Array.isArray(ptm.label_color) && ptm.label_color.length === 3 ? `rgb(${{ptm.label_color.join(',')}})` : 'black';
                    const labelCenter = labelCenteringValue.toLowerCase();
                    const textAnchor = anchorMap[labelCenter] || 'middle';
                    const labelText = protboxGroup.text(ptm.label).move(labelX, labelY + textOffsetY).font({{
                        size: settings.ptm_label_size || 10,
                        family: settings.ptm_label_font || 'Arial',
                        anchor: textAnchor,
                        leading: '1.2em'
                    }}).fill(labelColor).attr({{
                        'data-id': `${{uniprot}}_${{ptmKey}}_label`,
                        'data-type': 'ptm-label'
                    }});
                    workingEntries.push({{ element: labelText, offsetX: labelX - pbX, offsetY: (labelY + textOffsetY) - pbY, jsonX: labelX, jsonY: labelY }});
                    makeDraggable(labelText, 'ptm-label', `${{uniprot}}_${{ptmKey}}_label`, protboxId);
                    bindPtmContextMenu(labelText, protboxId, ptmMeta);
                }}
                const symbolX = (symbolSnapshot && symbolSnapshot.x !== undefined)
                    ? symbolSnapshot.x
                    : resolveCoord(existingOverride?.symbol_x ?? ptm.symbol_x, newX, resetSymbol);
                const symbolY = (symbolSnapshot && symbolSnapshot.y !== undefined)
                    ? symbolSnapshot.y
                    : resolveCoord(existingOverride?.symbol_y ?? ptm.symbol_y, newY, resetSymbol);
                if (ptm.symbol) {{
                    const symbolColor = Array.isArray(ptm.symbol_color) && ptm.symbol_color.length === 3 ? `rgb(${{ptm.symbol_color.join(',')}})` : 'black';
                    const symbolText = protboxGroup.text(ptm.symbol).move(symbolX, symbolY + textOffsetY).font({{
                        size: ptm.symbol_size || settings.ptm_label_size || 10,
                        family: ptm.symbol_font || settings.ptm_label_font || 'Arial',
                        anchor: 'middle',
                        leading: '1.2em'
                    }}).fill(symbolColor).attr({{
                        'data-id': `${{uniprot}}_${{ptmKey}}_symbol`,
                        'data-type': 'ptm-symbol',
                        'pointer-events': 'none'
                    }});
                    workingEntries.push({{ element: symbolText, offsetX: symbolX - pbX, offsetY: (symbolY + textOffsetY) - pbY, jsonX: symbolX, jsonY: symbolY }});
                }}
                if (isPrimaryForProtein) {{
                    ptm.shape_x = newX;
                    ptm.shape_y = newY;
                    ptm.ptm_position = posKey || null;
                    ptm.label_x = labelX;
                    ptm.label_y = labelY;
                    ptm.symbol_x = symbolX;
                    ptm.symbol_y = symbolY;
                }}
                recordPtmOverride(protboxId, uniprot, ptmKey, {{
                    shape_x: newX,
                    shape_y: newY,
                    ptm_position: posKey || null,
                    label_x: labelX,
                    label_y: labelY,
                    label_centering: labelCenteringValue,
                    symbol_x: symbolX,
                    symbol_y: symbolY,
                    hidden: false
                }});
                updateHandleDists(protboxId);
                Shiny?.setInputValue('ptm_spawned', {{
                    protbox_id: protboxId,
                    uniprot,
                    ptm_key: ptmKey,
                    shape_x: newX,
                    shape_y: newY,
                    ptm_position: posKey || null,
                    label_x: labelX,
                    label_y: labelY,
                    label_centering: labelCenteringValue,
                    symbol_x: symbolX,
                    symbol_y: symbolY,
                    hidden: false
                }}, {{ priority: 'event' }});
                if (mkHistory && recordHistory) {{
                    const ptmSnapshot = {{
                        protbox_id: protboxId,
                        uniprot,
                        ptm_key: ptmKey,
                        shape: {{ x: newX, y: newY, posKey: posKey || null }}
                    }};
                    if (ptm.label) {{
                        ptmSnapshot.label = {{ x: labelX, y: labelY, center: labelCenteringValue }};
                    }}
                    if (ptm.symbol) {{
                        ptmSnapshot.symbol = {{ x: symbolX, y: symbolY }};
                    }}
                    const entry = {{
                        kind: 'ptm-add',
                        snapshot: ptmSnapshot
                    }};
                    entry.handlers = {{
                        undo: () => mkHistory.runWithoutRecording(() => {{
                            deletePtmFromProtbox(protboxId, ptmMeta, {{ suppressHistory: true }});
                        }}),
                        redo: () => mkHistory.runWithoutRecording(() => {{
                            spawnPtmForProtbox(protboxId, uniprot, ptmKey, {{
                                recordHistory: false,
                                snapshot: ptmSnapshot
                            }});
                        }})
                    }};
                    mkHistory.recordAction(entry);
                }}
                return true;
            }};
            console.log('m3: about to render protBoxes, count =', Array.isArray(protBoxes) ? protBoxes.length : 0);
            try {{ var dbgEl2 = document.getElementById('debug_json'); if (dbgEl2) dbgEl2.textContent = 'Rendering protboxes: ' + (Array.isArray(protBoxes) ? protBoxes.length : 0); }} catch(e){{}}
            // Read box preview stretch factor (server provides under _box_preview.y_stretch)
            const boxPreview = data._box_preview || {{}};
            const boxYStretch = typeof boxPreview.y_stretch === 'number' ? boxPreview.y_stretch : (boxPreview.y_stretch ? Number(boxPreview.y_stretch) : 1.0);
            const renderProtbox = (pb, index) => {{ 
                try {{
                    const id = pb.protbox_id || `unknown_${{index}}`, x = pb.x || 0, y = pb.y || 0, width = pb.width || 46, height = pb.height || 17;
                    // Apply Y-stretch to the box *position* only (do not alter height)
                    const yPos = Math.round(y * boxYStretch);
                    let selectedUniprot = pb.selected_uniprot || pb.proteins?.[0] || '';
                    currentSelected[id] = selectedUniprot;
                    let protein = proteinData[selectedUniprot] || {{}};
                    const overrideMap = pb?.ptm_overrides?.[selectedUniprot] || null;
                    const isPrimaryForProtein = isPrimaryProtbox(id, selectedUniprot);
                    let label = protein.label || pb.backup_label || 'Unknown';
                    let labelColor = protein.label_color || [0, 0, 0];
                    elementGroups[id] = [];
                    const fillRgb = entityColor(protein);
                    const labelRgb = Array.isArray(labelColor) && labelColor.length === 3 && labelColor.every(c => typeof c === 'number' && c >= 0 && c <= 255) ? `rgb(${{labelColor.join(',')}})` : 'black';
                    const rect = protboxGroup.rect(width, height).move(x, yPos).fill(fillRgb).stroke({{ color: 'black', width: 1 }}).attr({{ 'data-id': id, 'data-type': 'prot-box', 'data-tooltip': (pb && pb.tooltip) || (protein && protein.tooltip) || '' }});
                    // Tooltip handlers for prot-box
                    rect.node.addEventListener('mouseenter', e => {{
                        const tip = rect.attr('data-tooltip') || '';
                        if (tip) {{
                            const pos = toContainerPosition(e);
                            setTooltipContent(tip, false);
                            tooltip.style.display = 'block';
                            tooltip.style.left = `${{pos.x + tooltipOffsetX}}px`;
                            tooltip.style.top = `${{pos.y + tooltipOffsetY}}px`;
                        }}
                    }});
                    rect.node.addEventListener('mousemove', e => {{
                        if (tooltip.style.display === 'block') {{
                            const pos = toContainerPosition(e);
                            tooltip.style.left = `${{pos.x + tooltipOffsetX}}px`;
                            tooltip.style.top = `${{pos.y + tooltipOffsetY}}px`;
                        }}
                    }});
                    rect.node.addEventListener('mouseleave', () => {{ tooltip.style.display = 'none'; clearTooltipContent(); }});
                    const text = protboxGroup.text(label).move(x + width / 2, yPos + height / 2 + textOffsetY).font({{ size: settings.prot_label_size || 12, family: settings.prot_label_font || 'Arial', anchor: 'middle', leading: '1.2em' }}).fill(labelRgb).attr({{ 'data-id': id + '_label', 'data-type': 'prot-label', 'pointer-events': 'none' }});
                    elementGroups[id].push({{ element: text, offsetX: width / 2, offsetY: height / 2, jsonX: x + width / 2, jsonY: y + height / 2 }});
                    makeDraggable(rect, 'prot-box', id, id);
                    let northHas = false, southHas = false, westHas = false, eastHas = false;
                    if (selectedUniprot && selectedUniprot in proteinData) {{
                        const proteinPtms = proteinData[selectedUniprot].PTMs || {{}};
                        const overrideMapForProt = overrideMap || {{}};
                        for (const [ptm_key, ptm] of Object.entries(proteinPtms)) {{
                            if (!ptm_key || !ptm) continue;
                            const override = overrideMapForProt[ptm_key] || null;
                            if (!isPrimaryForProtein && !override) {{
                                continue;
                            }}
                            if (override?.hidden) {{
                                continue;
                            }}
                            const ptmMeta = {{ uniprot: selectedUniprot, ptmKey: ptm_key }};
                            const effectiveShapeX = resolveCoordinate(override?.shape_x, ptm.shape_x);
                            const effectiveShapeY = resolveCoordinate(override?.shape_y, ptm.shape_y);
                            if (effectiveShapeX === null || effectiveShapeY === null) {{
                                // Leave PTMs with no stored coordinates hidden until the user spawns them.
                                continue;
                            }}
                            const ptmX = effectiveShapeX;
                            const ptmY = effectiveShapeY;
                            const shape = (ptm.shape || 'circle').toLowerCase();
                            // Apply same Y-stretch transformation to PTM positions so they align with their protbox
                            const adjPtmY = Math.round(ptmY * boxYStretch);
                            const radius = settings.ptm_circle_radius || 5;
                            const ptmFill = entityColor(ptm);
                            const shapeObj = shape === 'circle'
                                ? protboxGroup.circle(radius * 2).cx(ptmX).cy(adjPtmY).fill(ptmFill).stroke({{ color: 'black', width: 1 }})
                                : protboxGroup.rect(radius * 2, radius * 2).move(ptmX - radius, adjPtmY - radius).fill(ptmFill).stroke({{ color: 'black', width: 1 }});
                            shapeObj.attr({{
                                'data-id': `${{selectedUniprot}}_${{ptm_key}}_shape`,
                                'data-type': 'ptm-shape',
                                'data-protbox-id': id,
                                'data-tooltip': (ptm && ptm.tooltip) || '',
                                'data-tooltip-html': (ptm && ptm.tooltip_html) || ''
                            }});
                            const effectivePosKey = override?.ptm_position || ptm.ptm_position || '';
                            if (effectivePosKey) {{
                                shapeObj.attr({{ 'data-pos-key': effectivePosKey }});
                            }}
                            // Tooltip handlers for PTM shape
                            shapeObj.node.addEventListener('mouseenter', e => {{
                                const htmlTip = shapeObj.attr('data-tooltip-html') || '';
                                const tip = htmlTip || shapeObj.attr('data-tooltip') || '';
                                if (tip) {{
                                    const pos = toContainerPosition(e);
                                    setTooltipContent(tip, Boolean(htmlTip));
                                    tooltip.style.display = 'block';
                                    tooltip.style.left = `${{pos.x + tooltipOffsetX}}px`;
                                    tooltip.style.top = `${{pos.y + tooltipOffsetY}}px`;
                                }}
                            }});
                            shapeObj.node.addEventListener('mousemove', e => {{
                                if (tooltip.style.display === 'block') {{
                                    const pos = toContainerPosition(e);
                                    tooltip.style.left = `${{pos.x + tooltipOffsetX}}px`;
                                    tooltip.style.top = `${{pos.y + tooltipOffsetY}}px`;
                                }}
                            }});
                            shapeObj.node.addEventListener('mouseleave', () => {{ tooltip.style.display = 'none'; clearTooltipContent(); }});
                            elementGroups[id].push({{ element: shapeObj, offsetX: ptmX - x, offsetY: ptmY - y, jsonX: ptmX, jsonY: ptmY }});
                            makeDraggable(shapeObj, 'ptm-shape', `${{selectedUniprot}}_${{ptm_key}}_shape`, id);
                            bindPtmContextMenu(shapeObj, id, ptmMeta);
                            const ptm_pos = effectivePosKey || '';
                            if (ptm_pos.startsWith('N')) northHas = true;
                            if (ptm_pos.startsWith('S')) southHas = true;
                            if (ptm_pos.startsWith('W')) westHas = true;
                            if (ptm_pos.startsWith('E')) eastHas = true;
                            if (ptm.label) {{
                                const labelColor = Array.isArray(ptm.label_color) && ptm.label_color.length === 3 && ptm.label_color.every(c => typeof c === 'number' && c >= 0 && c <= 255) ? `rgb(${{ptm.label_color.join(',')}})` : 'black';
                                const fallbackLabelX = toFiniteNumber(ptm.label_x);
                                const fallbackLabelY = toFiniteNumber(ptm.label_y);
                                const baseLabelX = fallbackLabelX !== null ? fallbackLabelX : ptmX;
                                const baseLabelY = fallbackLabelY !== null ? fallbackLabelY : ptmY;
                                const labelX = resolveCoordinate(override?.label_x, baseLabelX) ?? baseLabelX;
                                const labelY = resolveCoordinate(override?.label_y, baseLabelY) ?? baseLabelY;
                                const adjLabelY = Math.round(labelY * boxYStretch);
                                const labelCenteringValue = override?.label_centering || ptm.label_centering || 'center';
                                const labelCenter = labelCenteringValue.toLowerCase();
                                const textAnchor = anchorMap[labelCenter] || 'middle';
                                const labelText = protboxGroup.text(ptm.label).move(labelX, adjLabelY + textOffsetY).font({{ size: settings.ptm_label_size || 10, family: settings.ptm_label_font || 'Arial', anchor: textAnchor, leading: '1.2em' }}).fill(labelColor).attr({{ 'data-id': `${{selectedUniprot}}_${{ptm_key}}_label`, 'data-type': 'ptm-label' }});
                                elementGroups[id].push({{ element: labelText, offsetX: labelX - x, offsetY: (labelY + textOffsetY) - y, jsonX: labelX, jsonY: labelY }});
                                makeDraggable(labelText, 'ptm-label', `${{selectedUniprot}}_${{ptm_key}}_label`, id);
                                bindPtmContextMenu(labelText, id, ptmMeta);
                            }}
                            if (ptm.symbol) {{
                                const fallbackSymbolX = toFiniteNumber(ptm.symbol_x);
                                const fallbackSymbolY = toFiniteNumber(ptm.symbol_y);
                                const baseSymbolX = fallbackSymbolX !== null ? fallbackSymbolX : ptmX;
                                const baseSymbolY = fallbackSymbolY !== null ? fallbackSymbolY : ptmY;
                                const symbolX = resolveCoordinate(override?.symbol_x, baseSymbolX) ?? baseSymbolX;
                                const symbolY = resolveCoordinate(override?.symbol_y, baseSymbolY) ?? baseSymbolY;
                                const adjSymbolY = Math.round(symbolY * boxYStretch);
                                const symbolColor = Array.isArray(ptm.symbol_color) && ptm.symbol_color.length === 3 && ptm.symbol_color.every(c => typeof c === 'number' && c >= 0 && c <= 255) ? `rgb(${{ptm.symbol_color.join(',')}})` : 'black';
                                const symbolText = protboxGroup.text(ptm.symbol).move(symbolX, adjSymbolY + textOffsetY).font({{ size: ptm.symbol_size || settings.ptm_label_size || 10, family: ptm.symbol_font || settings.ptm_label_font || 'Arial', anchor: 'middle', leading: '1.2em' }}).fill(symbolColor).attr({{ 'data-id': `${{selectedUniprot}}_${{ptm_key}}_symbol`, 'data-type': 'ptm-symbol', 'pointer-events': 'none' }});
                                elementGroups[id].push({{ element: symbolText, offsetX: symbolX - x, offsetY: (symbolY + textOffsetY) - y, jsonX: symbolX, jsonY: symbolY }});
                            }}
                        }}
                    }}
                                        const handleDistNorth = northHas ? 10 : 5;
                                        const handleDistSouth = southHas ? 10 : 5;
                                        const handleDistWest = westHas ? 10 : 5;
                                        const handleDistEast = eastHas ? 10 : 5;
                    elementGroups[id].push(...[ {{ element: handleGroup.line(x, yPos - handleDistNorth, x + width, yPos - handleDistNorth).stroke({{ color: 'red', width: 1 }}).attr({{ 'data-type': 'handle-line', 'data-id': id + '_north_line' }}) }}, {{ element: handleGroup.line(x, yPos + height + handleDistSouth, x + width, yPos + height + handleDistSouth).stroke({{ color: 'red', width: 1 }}).attr({{ 'data-type': 'handle-line', 'data-id': id + '_south_line' }}) }}, {{ element: handleGroup.line(x - handleDistWest, yPos, x - handleDistWest, yPos + height).stroke({{ color: 'red', width: 1 }}).attr({{ 'data-type': 'handle-line', 'data-id': id + '_west_line' }}) }}, {{ element: handleGroup.line(x + width + handleDistEast, yPos, x + width + handleDistEast, yPos + height).stroke({{ color: 'red', width: 1 }}).attr({{ 'data-type': 'handle-line', 'data-id': id + '_east_line' }}) }} ]);
                                        elementGroups[id].forEach(el => el.element.attr('data-type').startsWith('handle-') && el.element.hide());
                                        protboxHandleDists[id] = {{North: handleDistNorth,South: handleDistSouth,West: handleDistWest,East: handleDistEast}};
                                        const spacing = settings.ptm_circle_spacing || 4;
                                        // Use the stretched y-position when computing snap points / handle lines
                                        protboxSnapPoints[id] = {{'N1': {{x: x + width * 0.2, y: yPos - spacing}},'N2': {{x: x + width * 0.5, y: yPos - spacing}},'N3': {{x: x + width * 0.8, y: yPos - spacing}},'S1': {{x: x + width * 0.2, y: yPos + height + spacing}},'S2': {{x: x + width * 0.5, y: yPos + height + spacing}},'S3': {{x: x + width * 0.8, y: yPos + height + spacing}},'W1': {{x: x - spacing, y: yPos + height * 0.33 - 2}},'W2': {{x: x - spacing, y: yPos + height * 0.66 + 2}},'E1': {{x: x + width + spacing, y: yPos + height * 0.33 - 2}},'E2': {{x: x + width + spacing, y: yPos + height * 0.66 + 2}}}};
                                        for (const key in protboxSnapPoints[id]) {{
                                            const snapC = handleGroup.circle(4).cx(protboxSnapPoints[id][key].x).cy(protboxSnapPoints[id][key].y).fill('red').opacity(0.5).attr({{'data-type': 'ptm-snap-circle', 'data-pos-key': key, 'pointer-events': 'none'}});
                      elementGroups[id].push({{element: snapC}});
                      snapC.hide();
                    }}
                                        // Record the protbox's position using the stretched Y so other logic (arrows, attachments) uses the same coordinates
                                        protboxMap[id] = {{ x: x, y: yPos, width: width, height: height }};
                    rect.node.addEventListener('contextmenu', e => {{
                        e.preventDefault();
                        const existingMenu = document.querySelector('.context-menu');
                        if (existingMenu) existingMenu.remove();
                        const menu = document.createElement('div');
                        menu.className = 'context-menu';
                        menu.style.position = 'absolute';
                        const menuPos = toContainerPosition(e);
                        menu.style.left = `${{menuPos.x + menuOffsetX}}px`;
                        menu.style.top = `${{menuPos.y + menuOffsetY}}px`;
                        menu.style.backgroundColor = 'white';
                        menu.style.border = '1px solid #ccc';
                        menu.style.padding = '5px';
                        menu.style.zIndex = '1000';
                        const ul = document.createElement('ul');
                        ul.style.listStyle = 'none';
                        ul.style.margin = '0';
                        ul.style.padding = '0';
                        const liProtein = document.createElement('li');
                        liProtein.textContent = 'Protein';
                        liProtein.style.padding = '5px 10px';
                        liProtein.style.cursor = 'pointer';
                        liProtein.style.position = 'relative';
                        const submenu = document.createElement('ul');
                        submenu.style.position = 'absolute';
                        submenu.style.left = '100%';
                        submenu.style.top = '0';
                        submenu.style.backgroundColor = 'white';
                        submenu.style.border = '1px solid #ccc';
                        submenu.style.padding = '0';
                        submenu.style.margin = '0';
                        submenu.style.display = 'none';
                        submenu.style.listStyle = 'none';
                        pb.proteins.forEach(uniprot => {{
                            const prot = proteinData[uniprot] || {{label: uniprot}};
                            const subLi = document.createElement('li');
                            subLi.textContent = prot.label || uniprot;
                            subLi.style.padding = '5px 10px';
                            subLi.style.cursor = 'pointer';
                            if (uniprot === currentSelected[id]) subLi.style.fontWeight = 'bold';
                            subLi.addEventListener('click', () => {{
                                switchProtein(id, uniprot);
                                menu.remove();
                            }});
                            submenu.appendChild(subLi);
                        }});
                        liProtein.appendChild(submenu);
                        liProtein.addEventListener('mouseenter', () => {{ submenu.style.display = 'block'; }});
                        liProtein.addEventListener('mouseleave', () => {{ submenu.style.display = 'none'; }});
                        ul.appendChild(liProtein);
                        const liPTMs = document.createElement('li');
                        liPTMs.textContent = 'PTMs';
                        liPTMs.style.padding = '5px 10px';
                        liPTMs.style.cursor = 'pointer';
                        liPTMs.style.position = 'relative';
                        const ptmSubmenu = document.createElement('ul');
                        ptmSubmenu.style.position = 'absolute';
                        ptmSubmenu.style.left = '100%';
                        ptmSubmenu.style.top = '0';
                        ptmSubmenu.style.backgroundColor = 'white';
                        ptmSubmenu.style.border = '1px solid #ccc';
                        ptmSubmenu.style.padding = '0';
                        ptmSubmenu.style.margin = '0';
                        ptmSubmenu.style.display = 'none';
                        ptmSubmenu.style.listStyle = 'none';
                        const uniprot = currentSelected[id];
                        const ptms = proteinData[uniprot]?.PTMs || {{}};
                        const overrideMapForMenu = pb?.ptm_overrides?.[uniprot] || {{}};
                        Object.entries(ptms).forEach(([ptm_key, ptm]) => {{
                            const subLi = document.createElement('li');
                            subLi.style.padding = '5px 10px';
                            subLi.style.cursor = 'pointer';
                            subLi.style.display = 'flex';
                            subLi.style.alignItems = 'center';
                            subLi.style.gap = '8px';
                            const ptmPreview = document.createElement('span');
                            ptmPreview.style.display = 'inline-flex';
                            ptmPreview.style.alignItems = 'center';
                            ptmPreview.style.justifyContent = 'center';
                            ptmPreview.style.width = '18px';
                            ptmPreview.style.height = '18px';
                            const ptmShape = ptm && typeof ptm.shape === 'string' ? ptm.shape.toLowerCase() : 'circle';
                            ptmPreview.style.borderRadius = ptmShape === 'circle' ? '50%' : '3px';
                            ptmPreview.style.border = '1px solid #222';
                            ptmPreview.style.fontSize = '10px';
                            ptmPreview.style.fontWeight = 'bold';
                            ptmPreview.style.color = '#000';
                            ptmPreview.style.lineHeight = '1';
                            ptmPreview.style.flexShrink = '0';
                            try {{
                                ptmPreview.style.backgroundColor = entityColor(ptm) || '#ccc';
                            }} catch (err) {{
                                ptmPreview.style.backgroundColor = '#ccc';
                            }}
                            ptmPreview.textContent = ptm.symbol || '';
                            const ptmLabel = document.createElement('span');
                            ptmLabel.textContent = ptm.label || ptm_key;
                            subLi.appendChild(ptmPreview);
                            subLi.appendChild(ptmLabel);
                            const override = overrideMapForMenu[ptm_key] || null;
                            const isHidden = override?.hidden === true;
                            const hasCoords = toFiniteNumber(ptm.shape_x) !== null && toFiniteNumber(ptm.shape_y) !== null;
                            const isSpawned = hasCoords && !isHidden;
                            if (isSpawned) {{
                                subLi.style.backgroundColor = 'grey';
                                subLi.style.cursor = 'not-allowed';
                            }} else {{
                                subLi.style.backgroundColor = 'white';
                                subLi.addEventListener('click', () => {{
                                    if (spawnPtmForProtbox(id, uniprot, ptm_key, {{
                                        spawnInCenter: true,
                                        resetLabelPosition: true
                                    }})) {{
                                        menu.remove();
                                    }}
                                }});
                            }}
                            ptmSubmenu.appendChild(subLi);
                        }});
                        liPTMs.appendChild(ptmSubmenu);
                        liPTMs.addEventListener('mouseenter', () => {{ ptmSubmenu.style.display = 'block'; }});
                        liPTMs.addEventListener('mouseleave', () => {{ ptmSubmenu.style.display = 'none'; }});
                        ul.appendChild(liPTMs);
                        const liLinks = document.createElement('li');
                        liLinks.textContent = 'Links';
                        liLinks.style.padding = '5px 10px';
                        liLinks.style.cursor = 'pointer';
                        liLinks.style.position = 'relative';
                        const linksSubmenu = document.createElement('ul');
                        linksSubmenu.style.position = 'absolute';
                        linksSubmenu.style.left = '100%';
                        linksSubmenu.style.top = '0';
                        linksSubmenu.style.backgroundColor = 'white';
                        linksSubmenu.style.border = '1px solid #ccc';
                        linksSubmenu.style.padding = '0';
                        linksSubmenu.style.margin = '0';
                        linksSubmenu.style.display = 'none';
                        linksSubmenu.style.listStyle = 'none';
                        const createLinkMenuItem = (label, urlBuilder) => {{
                            const linkLi = document.createElement('li');
                            linkLi.textContent = label;
                            linkLi.style.padding = '5px 10px';
                            linkLi.style.cursor = 'pointer';
                            linkLi.addEventListener('click', () => {{
                                const uniprot = currentSelected[id];
                                if (uniprot) {{
                                    const safeId = encodeURIComponent(uniprot);
                                    window.open(urlBuilder(safeId), '_blank');
                                }}
                                menu.remove();
                            }});
                            return linkLi;
                        }};
                        linksSubmenu.appendChild(createLinkMenuItem('UniProt', safeId => `https://www.uniprot.org/uniprotkb/${{safeId}}/entry`));
                        linksSubmenu.appendChild(createLinkMenuItem('PSP', safeId => `https://www.phosphosite.org/uniprotAccAction?id=${{safeId}}`));
                        liLinks.appendChild(linksSubmenu);
                        liLinks.addEventListener('mouseenter', () => {{ linksSubmenu.style.display = 'block'; }});
                        liLinks.addEventListener('mouseleave', () => {{ linksSubmenu.style.display = 'none'; }});
                        ul.appendChild(liLinks);
                        menu.appendChild(ul);
                        container.appendChild(menu);
                        const removeMenu = () => {{
                            if (menu.parentNode) menu.remove();
                            document.removeEventListener('click', removeMenu);
                        }};
                        document.addEventListener('click', removeMenu);
                    }});
                }} catch (e) {{ }}
            }};
            function restoreProtboxSnapshot(snapshot, insertIndex = protBoxes.length) {{
                if (!snapshot) return null;
                const clone = cloneData(snapshot);
                if (!clone) return null;
                const safeIndex = Math.max(0, Math.min(insertIndex, protBoxes.length));
                protBoxes.splice(safeIndex, 0, clone);
                const numericId = Number(clone.protbox_id);
                if (Number.isFinite(numericId)) {{
                    protboxCounter = Math.max(protboxCounter, numericId);
                }}
                renderProtbox(clone, safeIndex);
                registerProtboxBaseGeometry(clone);
                return clone;
            }}
            const registerProtboxBaseGeometry = (pb) => {{
                if (!pb) return;
                const id = pb.protbox_id || `unknown_${{protBoxes.indexOf(pb)}}`;
                if (!id) return;
                const width = pb.width || 46;
                const height = pb.height || 17;
                protboxMap[id] = {{ x: pb.x || 0, y: pb.y || 0, width, height }};
            }};
            protBoxes.forEach((pb, index) => renderProtbox(pb, index));
            protBoxes.forEach(pb => registerProtboxBaseGeometry(pb));
            const ensureCompoundClientId = (compound, fallbackIndex = 0) => {{
                if (!compound) return null;
                if (compound._client_id) {{
                    registerCompoundAutoId(compound._client_id);
                    return compound._client_id;
                }}
                if (compound.compound_id !== undefined && compound.compound_id !== null) {{
                    compound._client_id = `compound_${{compound.compound_id}}`;
                    return compound._client_id;
                }}
                const generated = `compound_auto_${{++compoundAutoCounter}}`;
                compound._client_id = generated;
                registerCompoundAutoId(generated);
                return compound._client_id;
            }};
            const renderCompound = (compound, index = 0, options = {{}}) => {{
                if (!compound) return null;
                const id = ensureCompoundClientId(compound, index);
                const x = compound.x || 0;
                const y = compound.y || 0;
                const w = compound.width || 14;
                const h = compound.height || 14;
                const r = Math.max(w, h) / 2;
                const label = compound.label || (compound.kegg_compound || '').replace(/^cpd:/, '') || 'C?';
                const stroke = compound.border_color || 'black';
                const fill = compound.bgcolor || '#FFFFFF';
                const textColor = compound.fgcolor || '#000000';
                const g = compoundGroup.group().attr({{ 'data-id': id, 'data-type': 'compound' }});
                const circle = g.circle(r * 2).cx(x).cy(y).fill(fill).stroke({{ color: stroke, width: 1 }});
                circle.attr({{
                    'data-role': 'compound-shape',
                    'data-original-stroke': stroke,
                    'data-original-stroke-width': 1
                }});
                const labelOffset = settings.compound_label_offset || 6;
                const t = g.text(label)
                    .font({{
                        size: settings.compound_label_size || 10,
                        family: settings.compound_label_font || 'Arial',
                        anchor: 'middle',
                        leading: '1.2em'
                    }})
                    .fill(textColor);
                t.center(x, y + r + labelOffset);
                t.attr({{
                    'data-id': `${{id}}_label`,
                    'data-type': 'compound-label',
                    'pointer-events': 'none',
                    'data-original-fill': textColor
                }});
                makeDraggable(g, 'compound', id);
                if (mkHistory && options.recordHistory) {{
                    const snapshot = cloneData(compound);
                    if (snapshot) {{
                        const entry = {{
                            kind: 'compound-add',
                            snapshot,
                            index
                        }};
                        entry.handlers = {{
                            undo: () => mkHistory.runWithoutRecording(() => {{
                                deleteCompoundById(snapshot._client_id, {{ suppressHistory: true }});
                            }}),
                            redo: () => mkHistory.runWithoutRecording(() => {{
                                restoreCompoundSnapshot(snapshot, entry.index);
                            }})
                        }};
                        mkHistory.recordAction(entry);
                    }}
                }}
                return g;
            }};
            const restoreCompoundSnapshot = (snapshot, insertIndex = compounds.length) => {{
                if (!snapshot) return null;
                const clone = cloneData(snapshot);
                if (!clone) return null;
                const safeIndex = Math.max(0, Math.min(insertIndex, compounds.length));
                ensureCompoundClientId(clone, safeIndex);
                compounds.splice(safeIndex, 0, clone);
                renderCompound(clone, safeIndex, {{ recordHistory: false }});
                return clone;
            }};
            const ensureTextClientId = (textBlock, fallbackIndex = 0) => {{
                if (!textBlock) return null;
                if (textBlock._client_id) {{
                    registerTextAutoId(textBlock._client_id);
                    return textBlock._client_id;
                }}
                if (textBlock.text_id !== undefined && textBlock.text_id !== null) {{
                    textBlock._client_id = `text_${{textBlock.text_id}}`;
                    return textBlock._client_id;
                }}
                const generated = `text_auto_${{++textAutoCounter}}`;
                textBlock._client_id = generated;
                registerTextAutoId(generated);
                return textBlock._client_id;
            }};
            const renderTextBox = (textBlock, index = 0, options = {{}}) => {{
                if (!textBlock) return null;
                const id = ensureTextClientId(textBlock, index);
                const x = textBlock.x || 0;
                const y = textBlock.y || 0;
                const w = textBlock.width || 60;
                const h = textBlock.height || 20;
                const label = textBlock.label || '';
                const stroke = textBlock.border_color || 'black';
                const fill = textBlock.bgcolor || '#FFFFFF';
                const textColor = textBlock.fgcolor || '#000000';
                const g = textGroup.group().attr({{ 'data-id': id, 'data-type': 'text-box' }});
                const rect = g.rect(w, h).move(x, y).fill(fill).stroke({{ color: stroke, width: 1 }}).radius(6);
                rect.attr({{
                    'data-role': 'text-rect',
                    'data-original-stroke': stroke,
                    'data-original-stroke-width': 1
                }});
                const cx = x + w / 2;
                const cy = y + h / 2;
                const t = g.text(label)
                    .font({{
                        size: settings.textbox_label_size || 11,
                        family: settings.textbox_label_font || 'Arial',
                        anchor: 'middle',
                        leading: '1.2em'
                    }})
                    .fill(textColor);
                t.center(cx, cy);
                t.attr({{
                    'dominant-baseline': 'middle',
                    'data-id': `${{id}}_label`,
                    'data-type': 'text-label',
                    'pointer-events': 'none',
                    'data-original-fill': textColor
                }});
                makeDraggable(g, 'text-box', id);
                if (mkHistory && options.recordHistory) {{
                    const snapshot = cloneData(textBlock);
                    if (snapshot) {{
                        const entry = {{
                            kind: 'text-add',
                            snapshot,
                            index
                        }};
                        entry.handlers = {{
                            undo: () => mkHistory.runWithoutRecording(() => {{
                                deleteTextBoxById(snapshot._client_id, {{ suppressHistory: true }});
                            }}),
                            redo: () => mkHistory.runWithoutRecording(() => {{
                                restoreTextSnapshot(snapshot, entry.index);
                            }})
                        }};
                        mkHistory.recordAction(entry);
                    }}
                }}
                return g;
            }};
            const restoreTextSnapshot = (snapshot, insertIndex = textBlocks.length) => {{
                if (!snapshot) return null;
                const clone = cloneData(snapshot);
                if (!clone) return null;
                const safeIndex = Math.max(0, Math.min(insertIndex, textBlocks.length));
                ensureTextClientId(clone, safeIndex);
                textBlocks.splice(safeIndex, 0, clone);
                renderTextBox(clone, safeIndex, {{ recordHistory: false }});
                return clone;
            }};
            const createProtboxAtPosition = (uniprot, svgX, svgY, options = {{}}) => {{
                const protein = ensureProteinRecord(uniprot);
                if (!protein) {{
                    console.warn('No protein data for', uniprot);
                    return;
                }}
                const width = protein.protbox_width || protboxDefaults.width;
                const height = protein.protbox_height || protboxDefaults.height;
                const x = Math.round((svgX - width / 2));
                const displayY = svgY - height / 2;
                const y = boxYStretch ? Math.round(displayY / boxYStretch) : Math.round(displayY);
                const newId = String(++protboxCounter);
                const geneSymbol = protein.gene_symbol || protein.label || protein.backup_label || uniprot;
                const newProtbox = {{
                    protbox_id: newId,
                    proteins: [uniprot],
                    selected_uniprot: uniprot,
                    backup_label: geneSymbol,
                    x,
                    y,
                    width,
                    height,
                    ptm_overrides: {{ [uniprot]: {{}} }}
                }};
                protBoxes.push(newProtbox);
                currentSelected[newId] = uniprot;
                renderProtbox(newProtbox, protBoxes.length - 1);
                registerProtboxBaseGeometry(newProtbox);
                const ptmEntries = Object.entries(protein.PTMs || {{}});
                ptmEntries.forEach(([ptmKey], idx) => {{
                    const preferredSnapKey = prioritizedPtmPositions[idx] || null;
                    spawnPtmForProtbox(newId, uniprot, ptmKey, {{
                        preferredSnapKey,
                        resetLabelPosition: true,
                        resetSymbolPosition: true,
                        recordHistory: false
                    }});
                }});
                const proteinPayload = (() => {{
                    try {{
                        return JSON.parse(JSON.stringify(protein));
                    }} catch (err) {{
                        console.warn('Unable to clone protein payload for add_protbox', err);
                        return protein;
                    }}
                }})();
                Shiny?.setInputValue('add_protbox', {{
                    protbox: newProtbox,
                    uniprot,
                    protein: proteinPayload,
                    timestamp: Date.now()
                }}, {{ priority: 'event' }});
                if (mkHistory && options.recordHistory !== false) {{
                    const snapshot = cloneData(newProtbox);
                    if (snapshot) {{
                        const entry = {{
                            kind: 'protbox-add',
                            protbox: snapshot,
                            index: protBoxes.length - 1
                        }};
                        entry.handlers = {{
                            undo: () => mkHistory.runWithoutRecording(() => {{
                                deleteProtboxById(snapshot.protbox_id, {{ suppressHistory: true }});
                            }}),
                            redo: () => mkHistory.runWithoutRecording(() => {{
                                restoreProtboxSnapshot(snapshot, entry.index);
                            }})
                        }};
                        mkHistory.recordAction(entry);
                    }}
                }}
            }};
            const switchProtein = (protboxId, newUniprot) => {{
                const pb = protBoxes.find(p => p.protbox_id === protboxId);
                if (!pb) return;
                currentSelected[protboxId] = newUniprot;
                const rect = draw.find(`[data-id="${{protboxId}}"]`)[0];
                const text = draw.find(`[data-id="${{protboxId}}_label"]`)[0];
                const newProtein = ensureProteinRecord(newUniprot) || {{}};
                const newLabel = newProtein.label || pb.backup_label || 'Unknown';
                const newLabelColor = newProtein.label_color || [0, 0, 0];
                const fillRgb = entityColor(newProtein);
                const labelRgb = Array.isArray(newLabelColor) && newLabelColor.length === 3 ? `rgb(${{newLabelColor.join(',')}})` : 'black';
                rect.fill(fillRgb);
                text.text(newLabel).fill(labelRgb);
                const ptmTypes = ['ptm-shape', 'ptm-label', 'ptm-symbol'];
                const toRemove = elementGroups[protboxId].filter(el => ptmTypes.includes(el.element.attr('data-type')));
                toRemove.forEach(({{element}}) => element.remove());
                elementGroups[protboxId] = elementGroups[protboxId].filter(el => !ptmTypes.includes(el.element.attr('data-type')));
                let northHas = false, southHas = false, westHas = false, eastHas = false;
                const originalX = pb.x || 0;
                const originalY = pb.y || 0;
                const currentX = rect.x();
                const currentY = rect.y();
                const deltaX = currentX - originalX;
                const deltaY = currentY - originalY;
                const overrideMapForProt = pb?.ptm_overrides?.[newUniprot] || null;
                for (const [ptm_key, ptm] of Object.entries(newProtein.PTMs || {{}})) {{
                    if (!ptm_key || !ptm) continue;
                    const override = overrideMapForProt && overrideMapForProt[ptm_key] ? overrideMapForProt[ptm_key] : null;
                    if (override?.hidden) {{
                        continue;
                    }}
                    const ptmMeta = {{ uniprot: newUniprot, ptmKey: ptm_key }};
                    const resolvedShapeX = resolveCoordinate(override?.shape_x, ptm.shape_x);
                    const resolvedShapeY = resolveCoordinate(override?.shape_y, ptm.shape_y);
                    if (resolvedShapeX === null || resolvedShapeY === null) {{
                        continue;
                    }}
                    const adjusted_ptmX = resolvedShapeX + deltaX;
                    const adjusted_ptmY = resolvedShapeY + deltaY;
                    const shape = (ptm.shape || 'circle').toLowerCase();
                    const radius = settings.ptm_circle_radius || 5;
                    const ptmFill = entityColor(ptm);
                    const shapeObj = shape === 'circle'
                        ? protboxGroup.circle(radius * 2).cx(adjusted_ptmX).cy(adjusted_ptmY).fill(ptmFill).stroke({{ color: 'black', width: 1 }})
                        : protboxGroup.rect(radius * 2, radius * 2).move(adjusted_ptmX - radius, adjusted_ptmY - radius).fill(ptmFill).stroke({{ color: 'black', width: 1 }});
                    const ptm_pos = override?.ptm_position || ptm.ptm_position || '';
                    if (ptm_pos) {{
                        shapeObj.attr({{ 'data-pos-key': ptm_pos }});
                    }}
                    shapeObj.attr({{ 'data-id': `${{newUniprot}}_${{ptm_key}}_shape`, 'data-type': 'ptm-shape' }});
                    elementGroups[protboxId].push({{ element: shapeObj, offsetX: adjusted_ptmX - currentX, offsetY: adjusted_ptmY - currentY, jsonX: adjusted_ptmX, jsonY: adjusted_ptmY }});
                    makeDraggable(shapeObj, 'ptm-shape', `${{newUniprot}}_${{ptm_key}}_shape`, protboxId);
                    bindPtmContextMenu(shapeObj, protboxId, ptmMeta);
                    if (ptm_pos.startsWith('N')) northHas = true;
                    else if (ptm_pos.startsWith('S')) southHas = true;
                    else if (ptm_pos.startsWith('W')) westHas = true;
                    else if (ptm_pos.startsWith('E')) eastHas = true;
                    let shape_center_x, shape_center_y;
                    if (shape === 'circle') {{
                        shape_center_x = shapeObj.cx();
                        shape_center_y = shapeObj.cy();
                    }} else {{
                        shape_center_x = shapeObj.x() + radius;
                        shape_center_y = shapeObj.y() + radius;
                    }}
                    Shiny?.setInputValue('element_moved', {{ type: 'ptm-shape', id: `${{newUniprot}}_${{ptm_key}}_shape`, x: shape_center_x, y: shape_center_y, protbox_id: protboxId, ptm_position: ptm_pos || null }}, {{ priority: 'event' }});
                    if (ptm.label) {{
                        const fallbackLabelX = toFiniteNumber(ptm.label_x);
                        const fallbackLabelY = toFiniteNumber(ptm.label_y);
                        const baseLabelX = fallbackLabelX !== null ? fallbackLabelX : resolvedShapeX;
                        const baseLabelY = fallbackLabelY !== null ? fallbackLabelY : resolvedShapeY;
                        const resolvedLabelX = resolveCoordinate(override?.label_x, baseLabelX) ?? baseLabelX;
                        const resolvedLabelY = resolveCoordinate(override?.label_y, baseLabelY) ?? baseLabelY;
                        const adjusted_labelX = resolvedLabelX + deltaX;
                        const adjusted_labelY = resolvedLabelY + deltaY;
                        const labelColor = Array.isArray(ptm.label_color) && ptm.label_color.length === 3 ? `rgb(${{ptm.label_color.join(',')}})` : 'black';
                        const labelCenter = (override?.label_centering || ptm.label_centering || 'center').toLowerCase();
                        const textAnchor = anchorMap[labelCenter] || 'middle';
                        const labelText = protboxGroup.text(ptm.label).move(adjusted_labelX, adjusted_labelY + textOffsetY).font({{ size: settings.ptm_label_size || 10, family: settings.ptm_label_font || 'Arial', anchor: textAnchor, leading: '1.2em' }}).fill(labelColor).attr({{ 'data-id': `${{newUniprot}}_${{ptm_key}}_label`, 'data-type': 'ptm-label' }});
                        elementGroups[protboxId].push({{ element: labelText, offsetX: adjusted_labelX - currentX, offsetY: (adjusted_labelY + textOffsetY) - currentY, jsonX: adjusted_labelX, jsonY: adjusted_labelY }});
                        makeDraggable(labelText, 'ptm-label', `${{newUniprot}}_${{ptm_key}}_label`, protboxId);
                        bindPtmContextMenu(labelText, protboxId, ptmMeta);
                        Shiny?.setInputValue('element_moved', {{ type: 'ptm-label', id: `${{newUniprot}}_${{ptm_key}}_label`, x: adjusted_labelX, y: adjusted_labelY, protbox_id: protboxId }}, {{ priority: 'event' }});
                    }}
                    if (ptm.symbol) {{
                        const fallbackSymbolX = toFiniteNumber(ptm.symbol_x);
                        const fallbackSymbolY = toFiniteNumber(ptm.symbol_y);
                        const baseSymbolX = fallbackSymbolX !== null ? fallbackSymbolX : resolvedShapeX;
                        const baseSymbolY = fallbackSymbolY !== null ? fallbackSymbolY : resolvedShapeY;
                        const resolvedSymbolX = resolveCoordinate(override?.symbol_x, baseSymbolX) ?? baseSymbolX;
                        const resolvedSymbolY = resolveCoordinate(override?.symbol_y, baseSymbolY) ?? baseSymbolY;
                        const adjusted_symbolX = resolvedSymbolX + deltaX;
                        const adjusted_symbolY = resolvedSymbolY + deltaY;
                        const symbolColor = Array.isArray(ptm.symbol_color) && ptm.symbol_color.length === 3 ? `rgb(${{ptm.symbol_color.join(',')}})` : 'black';
                        const symbolText = protboxGroup.text(ptm.symbol).move(adjusted_symbolX, adjusted_symbolY + textOffsetY).font({{ size: ptm.symbol_size || settings.ptm_label_size || 10, family: ptm.symbol_font || settings.ptm_label_font || 'Arial', anchor: 'middle', leading: '1.2em' }}).fill(symbolColor).attr({{ 'data-id': `${{newUniprot}}_${{ptm_key}}_symbol`, 'data-type': 'ptm-symbol', 'pointer-events': 'none' }});
                        elementGroups[protboxId].push({{ element: symbolText, offsetX: adjusted_symbolX - currentX, offsetY: (adjusted_symbolY + textOffsetY) - currentY, jsonX: adjusted_symbolX, jsonY: adjusted_symbolY }});
                        Shiny?.setInputValue('element_moved', {{ type: 'ptm-symbol', id: `${{newUniprot}}_${{ptm_key}}_symbol`, x: adjusted_symbolX, y: adjusted_symbolY, protbox_id: protboxId }}, {{ priority: 'event' }});
                    }}
                    recordPtmOverride(protboxId, newUniprot, ptm_key, {{
                        shape_x: adjusted_ptmX,
                        shape_y: adjusted_ptmY,
                        ptm_position: ptm_pos || null,
                        label_x: (() => {{
                            const fallbackLabelXVal = toFiniteNumber(ptm.label_x);
                            const baseLabelXVal = fallbackLabelXVal !== null ? fallbackLabelXVal : resolvedShapeX;
                            return (resolveCoordinate(override?.label_x, baseLabelXVal) ?? baseLabelXVal) + deltaX;
                        }})(),
                        label_y: (() => {{
                            const fallbackLabelYVal = toFiniteNumber(ptm.label_y);
                            const baseLabelYVal = fallbackLabelYVal !== null ? fallbackLabelYVal : resolvedShapeY;
                            return (resolveCoordinate(override?.label_y, baseLabelYVal) ?? baseLabelYVal) + deltaY;
                        }})(),
                        label_centering: override?.label_centering || ptm.label_centering || 'center',
                        symbol_x: (() => {{
                            const fallbackSymbolXVal = toFiniteNumber(ptm.symbol_x);
                            const baseSymbolXVal = fallbackSymbolXVal !== null ? fallbackSymbolXVal : resolvedShapeX;
                            return (resolveCoordinate(override?.symbol_x, baseSymbolXVal) ?? baseSymbolXVal) + deltaX;
                        }})(),
                        symbol_y: (() => {{
                            const fallbackSymbolYVal = toFiniteNumber(ptm.symbol_y);
                            const baseSymbolYVal = fallbackSymbolYVal !== null ? fallbackSymbolYVal : resolvedShapeY;
                            return (resolveCoordinate(override?.symbol_y, baseSymbolYVal) ?? baseSymbolYVal) + deltaY;
                        }})()
                    }});
                }}
                protboxHandleDists[protboxId] = {{North: northHas ? 10 : 5, South: southHas ? 10 : 5, West: westHas ? 10 : 5, East: eastHas ? 10 : 5}};
                updateHandleDists(protboxId);
                if (attachments[protboxId]) {{
                    Object.keys(attachments[protboxId]).forEach(side => updateArrowPositions(protboxId, side));
                }}
                Shiny?.setInputValue('protein_switched', {{ protbox_id: protboxId, uniprot: newUniprot }}, {{ priority: 'event' }});
            }};
            
            
            // --- COMPOUNDS ---
            compounds.forEach((compound, idx) => renderCompound(compound, idx));


            // --- TEXT BLOCKS / MAP LABELS ---
            if (showTextBoxes) {{
              textBlocks.forEach((tb, idx) => renderTextBox(tb, idx));
            }}



            
            const attachedBy = {{}};
            arrows.forEach((arrow, index) => {{
                if (!arrow) return;
                const arrowId = arrowIdFromIndex(index);
                if (!arrowId) return;
                ['1', '2'].forEach(num => {{
                    const attA = arrow[`attached_arrow_${{num}}`];
                    const attE = arrow[`attached_end_${{num}}`];
                    if (attA && attE) {{
                        attachedBy[attA] = attachedBy[attA] || {{start: [], end: []}};
                        attachedBy[attA][attE].push({{slave: arrowId, slaveEnd: num === '1' ? 'start' : 'end'}});
                    }}
                }});
            }});
            arrows.forEach((arrow, index) => {{
                if (!arrow) return;
                const arrowId = arrowIdFromIndex(index);
                if (!arrowId) return;
                const id1 = arrow.protbox_id_1;
                const side1 = arrow.protbox_id_1_side;
                if (id1 && side1 && protboxMap[id1]) {{
                    attachments[id1] = attachments[id1] || {{}};
                    attachments[id1][side1] = attachments[id1][side1] || [];
                    attachments[id1][side1].push({{ type: 'start', arrow, arrowId, side: side1 }});
                }}
                const id2 = arrow.protbox_id_2;
                const side2 = arrow.protbox_id_2_side;
                if (id2 && side2 && protboxMap[id2]) {{
                    attachments[id2] = attachments[id2] || {{}};
                    attachments[id2][side2] = attachments[id2][side2] || [];
                    attachments[id2][side2].push({{ type: 'end', arrow, arrowId, side: side2 }});
                }}
            }});
            Object.keys(attachments).forEach(protboxId => {{
                Object.keys(attachments[protboxId]).forEach(side => {{
                    attachments[protboxId][side].sort((a, b) => {{
                        const targetIdA = a.type === 'start' ? a.arrow.protbox_id_2 : a.arrow.protbox_id_1;
                        const targetIdB = b.type === 'start' ? b.arrow.protbox_id_2 : b.arrow.protbox_id_1;
                        const isHorizontal = side === 'North' || side === 'South';
                        const keyA = isHorizontal ? (protboxMap[targetIdA]?.x || 0) + (protboxMap[targetIdA]?.width || 0) / 2 : (protboxMap[targetIdA]?.y || 0) + (protboxMap[targetIdA]?.height || 0) / 2;
                        const keyB = isHorizontal ? (protboxMap[targetIdB]?.x || 0) + (protboxMap[targetIdB]?.width || 0) / 2 : (protboxMap[targetIdB]?.y || 0) + (protboxMap[targetIdB]?.height || 0) / 2;
                        return keyA - keyB;
                    }});
                }});
            }});
            if (settings.show_groups && groups?.length) {{
                groups.forEach(group => {{
                    try {{
                        const protboxIds = group.protbox_ids || [];
                        if (!protboxIds.length) return;
                        let minX = Infinity, maxX = -Infinity, minY = Infinity, maxYGroup = -Infinity;
                        protboxIds.forEach(id => {{
                            const pb = protBoxes.find(p => p.protbox_id === id);
                            if (pb) {{
                                minX = Math.min(minX, pb.x || 0);
                                maxX = Math.max(maxX, (pb.x || 0) + (pb.width || 0));
                                minY = Math.min(minY, pb.y || 0);
                                maxYGroup = Math.max(maxYGroup, (pb.y || 0) + (pb.height || 0));
                            }}
                        }});
                        if (minX !== Infinity) {{
                            const padding = 5;
                            protboxGroup.rect(maxX - minX + 2 * padding, maxYGroup - minY + 2 * padding).move(minX - padding, minY - padding).fill('none').stroke({{ color: 'black', width: 1, dasharray: '5,5' }});
                        }}
                    }} catch (e) {{ }}
                }});
            }}
            const drawArrows = (arrowList) => {{
                if (!showArrows) return;
                arrowList.forEach((arrow, index) => {{
                    try {{
                        const sourceIndex = arrows.indexOf(arrow);
                        const arrowId = arrowIdFromIndex(sourceIndex);
                        if (!arrowId) return;
                        const id1 = arrow.protbox_id_1, side1 = arrow.protbox_id_1_side, box1 = protboxMap[id1];
                        const id2 = arrow.protbox_id_2, side2 = arrow.protbox_id_2_side, box2 = protboxMap[id2];
                        let x1, y1;
                        if (arrow.attached_arrow_1) {{
                            const attId = arrow.attached_arrow_1;
                            const attEnd = arrow.attached_end_1;
                            const attLine = draw.find(`[data-id="${{attId}}"]`)[0];
                            x1 = parseFloat(attLine.attr(attEnd === 'start' ? 'x1' : 'x2'));
                            y1 = parseFloat(attLine.attr(attEnd === 'start' ? 'y1' : 'y2'));
                        }} else if (id1 && side1 && box1) {{
                            const attachList1 = attachments[id1][side1] || [];
                            let frac1 = arrow.line === 'dashed_arrow' ? 0.5 : (attachList1.findIndex(att => att.arrow === arrow) + 1) / (attachList1.length + 1);
                            const width1 = box1.width || 46, height1 = box1.height || 17;
                            const handleDist1 = arrow.line === 'dashed_arrow' ? 0 : protboxHandleDists[id1]?.[side1] || 5;
                            if (side1 === 'North') {{ x1 = box1.x + frac1 * width1; y1 = box1.y - handleDist1; }}
                            else if (side1 === 'South') {{ x1 = box1.x + frac1 * width1; y1 = box1.y + height1 + handleDist1; }}
                            else if (side1 === 'West') {{ x1 = box1.x - handleDist1; y1 = box1.y + frac1 * height1; }}
                            else if (side1 === 'East') {{ x1 = box1.x + width1 + handleDist1; y1 = box1.y + frac1 * height1; }}
                        }} else {{
                            x1 = arrow.x1 || 0;
                            y1 = arrow.y1 || 0;
                        }}
                        let x2, y2;
                        if (arrow.attached_arrow_2) {{
                            const attId = arrow.attached_arrow_2;
                            const attEnd = arrow.attached_end_2;
                            const attLine = draw.find(`[data-id="${{attId}}"]`)[0];
                            x2 = parseFloat(attLine.attr(attEnd === 'start' ? 'x1' : 'x2'));
                            y2 = parseFloat(attLine.attr(attEnd === 'start' ? 'y1' : 'y2'));
                        }} else if (id2 && side2 && box2) {{
                            const attachList2 = attachments[id2][side2] || [];
                            let frac2 = arrow.line === 'dashed_arrow' ? 0.5 : (attachList2.findIndex(att => att.arrow === arrow) + 1) / (attachList2.length + 1);
                            const width2 = box2.width || 46, height2 = box2.height || 17;
                            const handleDist2 = arrow.line === 'dashed_arrow' ? 0 : protboxHandleDists[id2]?.[side2] || 5;
                            if (side2 === 'North') {{ x2 = box2.x + frac2 * width2; y2 = box2.y - handleDist2; }}
                            else if (side2 === 'South') {{ x2 = box2.x + frac2 * width2; y2 = box2.y + height2 + handleDist2; }}
                            else if (side2 === 'West') {{ x2 = box2.x - handleDist2; y2 = box2.y + frac2 * height2; }}
                            else if (side2 === 'East') {{ x2 = box2.x + width2 + handleDist2; y2 = box2.y + frac2 * height2; }}
                        }} else {{
                            x2 = arrow.x2 || 0;
                            y2 = arrow.y2 || 0;
                        }}
                        const lineType = arrow.line || 'arrow';
                        let strokeOpts = {{ color: 'black', width: 1 }};
                        if (lineType === 'dashed_arrow') {{ strokeOpts.dasharray = '5,5'; }}
                        const hitboxWidth = 10; // Adjust for desired hit area (5 pixels each side)
                        const hitbox = arrowGroup.line(x1, y1, x2, y2).stroke({{ color: 'transparent', width: hitboxWidth }}).attr({{ 'data-id': arrowId + '_hit', 'data-type': 'arrow-hitbox', 'data-arrow-id': arrowId }});
                        const line = arrowGroup.line(x1, y1, x2, y2).stroke(strokeOpts).attr({{ 'data-id': arrowId, 'data-type': 'arrow', 'data-line-type': lineType, 'pointer-events': 'none' }});
                        const dx = x2 - x1, dy = y2 - y1, angle = Math.atan2(dy, dx), arrowSize = 5;
                        let arrowHead = null;
                        if (lineType === 'arrow') {{
                            arrowHead = arrowGroup.polygon([[x2, y2],[x2 - arrowSize * Math.cos(angle + Math.PI / 6), y2 - arrowSize * Math.sin(angle + Math.PI / 6)],[x2 - arrowSize * Math.cos(angle - Math.PI / 6), y2 - arrowSize * Math.sin(angle - Math.PI / 6)]]).fill('black').stroke({{ color: 'black', width: 1 }}).attr({{ 'data-id': arrowId + '_head', 'data-type': 'arrow-head' }});
                        }} else if (lineType === 'inhibition') {{
                            let barX1, barY1, barX2, barY2;
                            if (side2 === 'North' || side2 === 'South') {{barX1 = x2 - arrowSize;barY1 = y2;barX2 = x2 + arrowSize;barY2 = y2;}} 
                            else if (side2 === 'West' || side2 === 'East') {{barX1 = x2;barY1 = y2 - arrowSize;barX2 = x2;barY2 = y2 + arrowSize;}}
                            else {{const perp = angle + Math.PI / 2;barX1 = x2 + Math.cos(perp) * arrowSize;barY1 = y2 + Math.sin(perp) * arrowSize;barX2 = x2 - Math.cos(perp) * arrowSize;barY2 = y2 - Math.sin(perp) * arrowSize;}}
                            arrowHead = arrowGroup.line(barX1, barY1, barX2, barY2).stroke({{ color: 'black', width: 2 }}).attr({{ 'data-id': arrowId + '_head', 'data-type': 'arrow-head' }});
                        }}
                        let startHandle = null;
                        if (!arrow.attached_arrow_1) {{
                            startHandle = arrowGroup.circle(6).cx(x1).cy(y1).fill('red').stroke({{ color: 'black', width: 1 }}).attr({{ 'data-id': arrowId + '_start', 'data-type': 'arrow-start' }}).hide();
                        }}
                        let endHandle = null;
                        if (!arrow.attached_arrow_2) {{
                            endHandle = arrowGroup.circle(6).cx(x2).cy(y2).fill('red').stroke({{ color: 'black', width: 1 }}).attr({{ 'data-id': arrowId + '_end', 'data-type': 'arrow-end' }}).hide();
                        }}
                        arrowHandleGroups[arrowId] = {{ startHandle, endHandle }};
                        makeDraggable(hitbox, 'arrow', arrowId);
                        if (startHandle) makeDraggable(startHandle, 'arrow-start', arrowId + '_start');
                        if (endHandle) makeDraggable(endHandle, 'arrow-end', arrowId + '_end');
                    }} catch (e) {{ }}
                }});
            }};
            if (arrows?.length) {{
                const nonLines = arrows.filter(a => a && a.line !== 'line');
                const lines = arrows.filter(a => a && a.line === 'line');
                drawArrows(nonLines);
                drawArrows(lines);
            }}
            draw.node.addEventListener('contextmenu', e => {{
                if (!isBackgroundTarget(e.target)) return;
                e.preventDefault();
                const existingMenu = document.querySelector('.context-menu');
                if (existingMenu) existingMenu.remove();
                const rect = draw.node.getBoundingClientRect();
                const screenX = e.clientX - rect.left;
                const screenY = e.clientY - rect.top;
                const svgX = viewBox.x + (screenX / rect.width) * viewBox.width;
                const svgY = viewBox.y + (screenY / rect.height) * viewBox.height;
                const menu = document.createElement('div');
                menu.className = 'context-menu';
                menu.style.position = 'absolute';
                const rootPos = toContainerPosition(e);
                menu.style.left = `${{rootPos.x + menuOffsetX}}px`;
                menu.style.top = `${{rootPos.y + menuOffsetY}}px`;
                menu.style.backgroundColor = 'white';
                menu.style.border = '1px solid #ccc';
                menu.style.padding = '5px';
                menu.style.zIndex = '1000';
                const ul = document.createElement('ul');
                ul.style.listStyle = 'none';
                ul.style.margin = '0';
                ul.style.padding = '0';
                const liAddInteraction = document.createElement('li');
                liAddInteraction.textContent = 'Add Interaction';
                liAddInteraction.style.padding = '5px 10px';
                liAddInteraction.style.cursor = 'pointer';
                liAddInteraction.style.position = 'relative';
                const interactionSubmenu = document.createElement('ul');
                interactionSubmenu.style.position = 'absolute';
                interactionSubmenu.style.left = '100%';
                interactionSubmenu.style.top = '0';
                interactionSubmenu.style.backgroundColor = 'white';
                interactionSubmenu.style.border = '1px solid #ccc';
                interactionSubmenu.style.padding = '0';
                interactionSubmenu.style.margin = '0';
                interactionSubmenu.style.display = 'none';
                interactionSubmenu.style.listStyle = 'none';
                const createInteractionItem = (label, type) => {{
                    const interactionLi = document.createElement('li');
                    interactionLi.textContent = label;
                    interactionLi.style.padding = '5px 10px';
                    interactionLi.style.cursor = 'pointer';
                    interactionLi.addEventListener('click', () => {{
                        addNewArrow(type, svgX, svgY);
                        menu.remove();
                    }});
                    interactionSubmenu.appendChild(interactionLi);
                }};
                createInteractionItem('Add Arrow', 'arrow');
                createInteractionItem('Add Line', 'line');
                createInteractionItem('Add Inhibitor Line', 'inhibition');
                liAddInteraction.appendChild(interactionSubmenu);
                liAddInteraction.addEventListener('mouseenter', () => {{ interactionSubmenu.style.display = 'block'; }});
                liAddInteraction.addEventListener('mouseleave', () => {{ interactionSubmenu.style.display = 'none'; }});
                ul.appendChild(liAddInteraction);
                const liAddProtbox = document.createElement('li');
                liAddProtbox.textContent = 'Add Protbox';
                liAddProtbox.style.padding = '5px 10px';
                liAddProtbox.style.cursor = 'pointer';
                liAddProtbox.style.position = 'relative';
                const protboxSubmenu = document.createElement('div');
                protboxSubmenu.style.position = 'absolute';
                protboxSubmenu.style.left = '100%';
                protboxSubmenu.style.top = '0';
                protboxSubmenu.style.backgroundColor = 'white';
                protboxSubmenu.style.border = '1px solid #ccc';
                protboxSubmenu.style.padding = '8px';
                protboxSubmenu.style.margin = '0';
                protboxSubmenu.style.display = 'none';
                protboxSubmenu.style.width = '260px';
                protboxSubmenu.style.boxShadow = '0 2px 6px rgba(0,0,0,0.15)';
                protboxSubmenu.classList.add('context-submenu');
                const searchInput = document.createElement('input');
                searchInput.type = 'text';
                searchInput.placeholder = 'Regex UniProt or gene symbol';
                searchInput.style.width = '100%';
                searchInput.style.boxSizing = 'border-box';
                searchInput.style.padding = '4px 6px';
                const statusRow = document.createElement('div');
                statusRow.style.fontSize = '11px';
                statusRow.style.marginTop = '6px';
                statusRow.style.color = '#666';
                const resultsList = document.createElement('ul');
                resultsList.style.listStyle = 'none';
                resultsList.style.margin = '6px 0 0 0';
                resultsList.style.padding = '0';
                resultsList.style.maxHeight = '200px';
                resultsList.style.overflowY = 'auto';
                resultsList.style.borderTop = '1px solid #eee';
                const stopEvent = (evt) => {{
                    evt.stopPropagation();
                    evt.stopImmediatePropagation();
                }};
                ['click', 'mousedown', 'mouseup', 'keydown', 'keyup', 'keypress'].forEach(evtName => {{
                    searchInput.addEventListener(evtName, stopEvent);
                }});
                protboxSubmenu.addEventListener('click', stopEvent);
                protboxSubmenu.addEventListener('mousedown', stopEvent);
                const updateProtboxResults = () => {{
                    resultsList.innerHTML = '';
                    if (!proteinSearchIndex.length) {{
                        statusRow.textContent = 'Protein catalogue unavailable.';
                        return;
                    }}
                    const query = searchInput.value.trim();
                    if (!query) {{
                        statusRow.textContent = 'Type a regex to search proteins.';
                        return;
                    }}
                    let regex = null;
                    try {{
                        regex = new RegExp(query, 'i');
                    }} catch (err) {{
                        statusRow.textContent = 'Invalid regex pattern.';
                        return;
                    }}
                    const matches = [];
                    for (const entry of proteinSearchIndex) {{
                        if (regex.test(entry.searchText)) {{
                            matches.push(entry);
                        }}
                        if (matches.length >= 40) break;
                    }}
                    if (!matches.length) {{
                        statusRow.textContent = 'No proteins matched.';
                        return;
                    }}
                    statusRow.textContent = `Showing ${{matches.length}} result${{matches.length === 1 ? '' : 's'}}`;
                    matches.forEach(entry => {{
                        const item = document.createElement('li');
                        item.style.display = 'flex';
                        item.style.alignItems = 'center';
                        item.style.gap = '6px';
                        item.style.padding = '4px 6px';
                        item.style.cursor = 'pointer';
                        const swatch = document.createElement('span');
                        const protein = searchSource[entry.uniprot];
                        swatch.style.display = 'inline-block';
                        swatch.style.width = '10px';
                        swatch.style.height = '10px';
                        swatch.style.border = '1px solid #ccc';
                        try {{
                            swatch.style.backgroundColor = entityColor(protein);
                        }} catch (err) {{
                            swatch.style.backgroundColor = '#ccc';
                        }}
                        const labelSpan = document.createElement('span');
                        labelSpan.textContent = `${{entry.uniprot}} - ${{entry.geneSymbol}}`;
                        item.appendChild(swatch);
                        item.appendChild(labelSpan);
                        item.addEventListener('mouseenter', () => {{ item.style.backgroundColor = '#eef'; }});
                        item.addEventListener('mouseleave', () => {{ item.style.backgroundColor = 'transparent'; }});
                        item.addEventListener('click', (evt) => {{
                            evt.stopPropagation();
                            createProtboxAtPosition(entry.uniprot, svgX, svgY);
                            menu.remove();
                        }});
                        resultsList.appendChild(item);
                    }});
                }};
                searchInput.addEventListener('input', updateProtboxResults);
                protboxSubmenu.appendChild(searchInput);
                protboxSubmenu.appendChild(statusRow);
                protboxSubmenu.appendChild(resultsList);
                liAddProtbox.appendChild(protboxSubmenu);
                liAddProtbox.addEventListener('mouseenter', () => {{
                    protboxSubmenu.style.display = 'block';
                    setTimeout(() => searchInput.focus(), 0);
                    updateProtboxResults();
                }});
                liAddProtbox.addEventListener('mouseleave', (evt) => {{
                    if (!liAddProtbox.contains(evt.relatedTarget)) {{
                        protboxSubmenu.style.display = 'none';
                    }}
                }});
                ul.appendChild(liAddProtbox);
                menu.appendChild(ul);
                container.appendChild(menu);
                const removeMenu = () => {{
                    if (menu.parentNode) menu.remove();
                    document.removeEventListener('click', removeMenu);
                }};
                document.addEventListener('click', removeMenu);
            }});
            const addNewArrow = (type, svgX, svgY) => {{
                const x1 = svgX - 20;
                const y1 = svgY;
                const x2 = svgX + 20;
                const y2 = svgY;
                const newArrow = {{ line: type, x1, y1, x2, y2 }};
                arrows.push(newArrow);
                drawArrows([newArrow]);
                Shiny?.setInputValue('add_arrow', {{ line: type, x1, y1, x2, y2 }}, {{ priority: 'event' }});
            }};
        }}
        // Expose initializeSvg for manual invocation from the browser console for debugging
        try {{
            window.initializeSvg = initializeSvg;
            console.log('m3: initializeSvg exposed');
        }} catch (e) {{
            console.log('m3: could not expose initializeSvg', e);
        }}
        console.log('m3: scheduling initializeSvg');
        // Increase delay slightly to allow assets to load
        setTimeout(initializeSvg, 1500);
        </script>
    """
    return ui.div(
        ui.div(
            "Instructions: Click on protein boxes, PTM shapes, PTM labels, or arrows to select them (red outline). Click arrow ends (red dots) to move one end. Drag arrow line to move entire arrow and unsnap from protboxes. Use arrow keys to move selected elements. Hold Shift for larger movements (10px). Press Escape to deselect. Scroll mouse wheel to zoom in/out. Drag background to pan. Click reset button to reset view. Hover over protein boxes or PTM shapes to see tooltips. Right-click on protein box to switch protein.",
            style="color: #666; fontSize: 12px; margin-bottom: 10px; padding: 5px; background-color: #f0f0f0; border-radius: 3px;"
        ),
        ui.HTML(f'''
            <div id="svgCanvas" style="width: {max_x}px; height: {max_y}px; background-color: white; position: relative;" tabindex="0">
                <div class="canvas-controls" style="display: none;">
                    <button id="reset-view" class="svg-control-btn" title="Reset View">
                        <i class="fa fa-search"></i>
                    </button>
                </div>
            </div>
        '''),
        ui.HTML(data_script + catalog_script + svg_js),
        class_="svg-container",
        **{"style": f"position: relative; width: {max_x + 20}px; height: {max_y + 20}px; overflow: auto; background-color: #f9f9f9; border: 1px solid #ccc;"}
    )

app_ui = ui.page_fluid(
    ui.head_content(
        ui.HTML("""
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.4/css/all.min.css">
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .svg-container {{ border: 1px solid #ccc; padding: 10px; background-color: #f9f9f9; position: relative; overflow: auto; }}
                #svgCanvas {{ background-color: white; display: block; position: relative; outline: none; }}
                #svgCanvas:focus {{ outline: 2px solid #007bff; outline-offset: -2px; }}
                #svgCanvas svg {{ width: 100%; height: 100%; display: block; position: absolute; top: 0; left: 0; }}
                .canvas-controls {{
                    position: absolute;
                    right: 10px;
                    bottom: 10px;
                    display: flex;
                    gap: 6px;
                }}
                .svg-control-btn {{
                    background-color: #fff;
                    border: 1px solid #ccc;
                    border-radius: 4px;
                    width: 32px;
                    height: 32px;
                    display: flex;
                    align-items: center;
                    justify-content: center;
                    cursor: pointer;
                    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                }}
                .svg-control-btn:disabled {{
                    opacity: 0.45;
                    cursor: default;
                }}
                .svg-control-btn:hover:not(:disabled) {{
                    background-color: #f0f0f0;
                }}
            </style>
        """)
    ),
    ui.h2("Pathway Visualization"),
    ui.output_ui("pathway_plot"),
    ui.output_text("debug_json"),
    ui.input_action_button("save_json", "Save Changes to JSON")
)


def server(input, output, session):
    json_data_reactive = reactive.Value(None)

    def _find_protbox_entry(json_data, protbox_id):
        if not json_data or protbox_id is None:
            return None
        target = str(protbox_id)
        for pb in json_data.get('protbox_data', []):
            if str(pb.get('protbox_id')) == target:
                if not isinstance(pb.get('ptm_overrides'), dict):
                    pb['ptm_overrides'] = {}
                return pb
        return None

    def _update_ptm_override(json_data, protbox_id, uniprot_id, ptm_key, updates):
        if not json_data or not updates or not uniprot_id or not ptm_key:
            return
        pb = _find_protbox_entry(json_data, protbox_id)
        if not pb:
            return
        overrides = pb.setdefault('ptm_overrides', {})
        prot_map = overrides.setdefault(uniprot_id, {})
        ptm_map = prot_map.setdefault(ptm_key, {})
        for key, value in updates.items():
            if value is None:
                continue
            ptm_map[key] = value

    def _parse_ptm_meta(element_id):
        if not element_id:
            return None, None
        parts = str(element_id).split('_')
        if len(parts) < 3:
            return None, None
        uniprot_id = parts[0]
        ptm_key = '_'.join(parts[1:-1])
        return uniprot_id, ptm_key

    def _parse_arrow_index(value):
        if value is None:
            return None
        if isinstance(value, int):
            return value
        try:
            string_value = str(value)
        except Exception:
            return None
        if string_value.startswith("arrow_"):
            suffix = string_value.split("_", 1)[1]
            if suffix.isdigit():
                return int(suffix)
        if string_value.isdigit():
            return int(string_value)
        try:
            return int(float(string_value))
        except (TypeError, ValueError):
            return None

    def _protboxes_for_uniprot(json_data, uniprot_id):
        if not json_data or not uniprot_id:
            return []
        matches = []
        for pb in json_data.get('protbox_data', []):
            proteins = pb.get('proteins') or []
            if uniprot_id in proteins:
                matches.append(pb)
        return matches

    def _is_primary_protbox(json_data, protbox_id, uniprot_id):
        if not uniprot_id:
            return True
        matches = _protboxes_for_uniprot(json_data, uniprot_id)
        if not matches:
            return True
        first_id = str(matches[0].get('protbox_id'))
        return first_id == str(protbox_id)

    @reactive.Effect
    def load_json():
        json_data_reactive.set(load_json_data())

    @output
    @render.ui
    def pathway_plot():
        json_data = json_data_reactive.get()
        if json_data is None:
            return ui.div("Error: Could not load JSON data.")
        return create_pathway_svg(json_data)

    @output
    @render.text
    def debug_json():
        json_data = json_data_reactive.get()
        if json_data is None:
            return "No JSON data loaded"
        return f"JSON keys: {list(json_data.keys())}\nProtbox count: {len(json_data.get('protbox_data', []))}"

    @reactive.Effect
    @reactive.event(input.element_moved)
    def update_positions():
        moved = input.element_moved()
        if not moved:
            return
        json_data = json_data_reactive.get()
        if not json_data:
            return
        element_type = moved.get('type')
        element_id = moved.get('id')
        x = moved.get('x')
        y = moved.get('y')
        protbox_id = moved.get('protbox_id')
        if element_type == 'prot-box':
            for pb in json_data['protbox_data']:
                if pb['protbox_id'] == element_id:
                    pb['x'] = float(x)
                    pb['y'] = float(y)
                    break
        elif element_type == 'ptm-shape':
            uniprot_id, ptm_key = _parse_ptm_meta(element_id)
            is_primary = _is_primary_protbox(json_data, protbox_id, uniprot_id)
            if uniprot_id and uniprot_id in json_data['protein_data'] and is_primary:
                ptm = json_data['protein_data'][uniprot_id]['PTMs'].get(ptm_key)
                if ptm:
                    ptm['shape_x'] = float(x)
                    ptm['shape_y'] = float(y)
            _update_ptm_override(json_data, protbox_id, uniprot_id, ptm_key, {
                'shape_x': float(x) if x is not None else None,
                'shape_y': float(y) if y is not None else None,
                'ptm_position': moved.get('ptm_position')
            })
        elif element_type == 'ptm-label':
            uniprot_id, ptm_key = _parse_ptm_meta(element_id)
            is_primary = _is_primary_protbox(json_data, protbox_id, uniprot_id)
            if uniprot_id and uniprot_id in json_data['protein_data'] and is_primary:
                ptm = json_data['protein_data'][uniprot_id]['PTMs'].get(ptm_key)
                if ptm:
                    ptm['label_x'] = float(x)
                    ptm['label_y'] = float(y)
                    if 'label_centering' in moved:
                        ptm['label_centering'] = moved['label_centering']
            _update_ptm_override(json_data, protbox_id, uniprot_id, ptm_key, {
                'label_x': float(x) if x is not None else None,
                'label_y': float(y) if y is not None else None,
                'label_centering': moved.get('label_centering')
            })
        elif element_type == 'ptm-symbol':
            uniprot_id, ptm_key = _parse_ptm_meta(element_id)
            is_primary = _is_primary_protbox(json_data, protbox_id, uniprot_id)
            if uniprot_id and uniprot_id in json_data['protein_data'] and is_primary:
                ptm = json_data['protein_data'][uniprot_id]['PTMs'].get(ptm_key)
                if ptm:
                    ptm['symbol_x'] = float(x)
                    ptm['symbol_y'] = float(y)
            _update_ptm_override(json_data, protbox_id, uniprot_id, ptm_key, {
                'symbol_x': float(x) if x is not None else None,
                'symbol_y': float(y) if y is not None else None
            })
        json_data_reactive.set(json_data)

    @reactive.Effect
    @reactive.event(input.protein_switched)
    def update_protein_selection():
        switched = input.protein_switched()
        if not switched:
            return
        json_data = json_data_reactive.get()
        if not json_data:
            return
        protbox_id = switched.get('protbox_id')
        uniprot = switched.get('uniprot')
        for pb in json_data['protbox_data']:
            if pb['protbox_id'] == protbox_id:
                pb['selected_uniprot'] = uniprot
                break
        json_data_reactive.set(json_data)

    @reactive.Effect
    @reactive.event(input.ptm_spawned)
    def spawn_ptm():
        spawned = input.ptm_spawned()
        if not spawned:
            return
        json_data = json_data_reactive.get()
        if not json_data:
            return
        uniprot = spawned.get('uniprot')
        ptm_key = spawned.get('ptm_key')
        protbox_id = spawned.get('protbox_id')
        is_primary = _is_primary_protbox(json_data, protbox_id, uniprot)
        if uniprot in json_data['protein_data'] and is_primary:
            ptms = json_data['protein_data'][uniprot]['PTMs']
            if ptm_key in ptms:
                ptm = ptms[ptm_key]
                ptm['shape_x'] = spawned['shape_x']
                ptm['shape_y'] = spawned['shape_y']
                ptm['ptm_position'] = spawned['ptm_position']
                ptm['label_x'] = spawned.get('label_x', ptm.get('label_x'))
                ptm['label_y'] = spawned.get('label_y', ptm.get('label_y'))
                if spawned.get('label_centering') is not None:
                    ptm['label_centering'] = spawned['label_centering']
                ptm['symbol_x'] = spawned.get('symbol_x', ptm.get('symbol_x'))
                ptm['symbol_y'] = spawned.get('symbol_y', ptm.get('symbol_y'))
        override_payload = {
            'shape_x': spawned.get('shape_x'),
            'shape_y': spawned.get('shape_y'),
            'ptm_position': spawned.get('ptm_position'),
            'label_x': spawned.get('label_x'),
            'label_y': spawned.get('label_y'),
            'label_centering': spawned.get('label_centering'),
            'symbol_x': spawned.get('symbol_x'),
            'symbol_y': spawned.get('symbol_y')
        }
        if 'hidden' in spawned:
            override_payload['hidden'] = bool(spawned.get('hidden'))
        _update_ptm_override(json_data, protbox_id, uniprot, ptm_key, override_payload)
        json_data_reactive.set(json_data)

    @reactive.Effect
    @reactive.event(input.add_arrow)
    def add_arrow():
        added = input.add_arrow()
        if not added:
            return
        json_data = json_data_reactive.get()
        if not json_data:
            return
        new_arrow = {
            'line': added['line'],
            'x1': added['x1'],
            'y1': added['y1'],
            'x2': added['x2'],
            'y2': added['y2']
        }
        json_data['arrows'].append(new_arrow)
        json_data_reactive.set(json_data)
    @reactive.Effect
    @reactive.event(input.add_protbox)
    def add_protbox():
        added = input.add_protbox()
        if not added:
            return
        json_data = json_data_reactive.get()
        if not json_data:
            return
        protbox = added.get('protbox')
        uniprot = added.get('uniprot')
        protein_payload = added.get('protein')
        if protbox:
            if 'ptm_overrides' not in protbox or not isinstance(protbox['ptm_overrides'], dict):
                protbox['ptm_overrides'] = {}
            json_data.setdefault('protbox_data', []).append(protbox)
        if uniprot and protein_payload:
            json_data.setdefault('protein_data', {})[uniprot] = protein_payload
        json_data_reactive.set(json_data)

    @reactive.Effect
    @reactive.event(input.delete_element)
    def delete_element():
        payload = input.delete_element()
        if not payload:
            return
        json_data = json_data_reactive.get()
        if not json_data:
            return
        element_type = payload.get('type')
        modified = False
        if element_type == 'prot-box':
            protbox_id = payload.get('protbox_id') or payload.get('id')
            if protbox_id is None:
                return
            protbox_id_str = str(protbox_id)
            protboxes = json_data.get('protbox_data', [])
            new_protboxes = [pb for pb in protboxes if str(pb.get('protbox_id')) != protbox_id_str]
            if len(new_protboxes) != len(protboxes):
                json_data['protbox_data'] = new_protboxes
                modified = True
            groups = json_data.get('groups')
            if isinstance(groups, list):
                updated_groups = []
                groups_changed = False
                for group in groups:
                    ids = group.get('protbox_ids') or []
                    filtered_ids = [gid for gid in ids if str(gid) != protbox_id_str]
                    if len(filtered_ids) != len(ids):
                        groups_changed = True
                    group['protbox_ids'] = filtered_ids
                    if filtered_ids:
                        updated_groups.append(group)
                    else:
                        groups_changed = True
                if len(updated_groups) != len(groups) or groups_changed:
                    json_data['groups'] = updated_groups
                    modified = True
            arrows = json_data.get('arrows', [])
            if isinstance(arrows, list):
                detached_entries = payload.get('detached_arrows') or []
                for entry in detached_entries:
                    arrow_index = entry.get('arrow_index')
                    if arrow_index is None:
                        arrow_index = _parse_arrow_index(entry.get('arrow_id'))
                    else:
                        arrow_index = _parse_arrow_index(arrow_index)
                    if arrow_index is None or arrow_index < 0 or arrow_index >= len(arrows):
                        continue
                    arrow = arrows[arrow_index]
                    if not isinstance(arrow, dict):
                        continue
                    end = entry.get('end')
                    if end not in {'start', 'end'}:
                        continue
                    x_val = entry.get('x')
                    y_val = entry.get('y')
                    if end == 'start':
                        arrow.pop('protbox_id_1', None)
                        arrow.pop('protbox_id_1_side', None)
                        if x_val is not None:
                            arrow['x1'] = float(x_val)
                        if y_val is not None:
                            arrow['y1'] = float(y_val)
                    else:
                        arrow.pop('protbox_id_2', None)
                        arrow.pop('protbox_id_2_side', None)
                        if x_val is not None:
                            arrow['x2'] = float(x_val)
                        if y_val is not None:
                            arrow['y2'] = float(y_val)
                    modified = True
                for arrow in arrows:
                    if not isinstance(arrow, dict):
                        continue
                    if str(arrow.get('protbox_id_1')) == protbox_id_str:
                        arrow.pop('protbox_id_1', None)
                        arrow.pop('protbox_id_1_side', None)
                        modified = True
                    if str(arrow.get('protbox_id_2')) == protbox_id_str:
                        arrow.pop('protbox_id_2', None)
                        arrow.pop('protbox_id_2_side', None)
                        modified = True
        elif element_type == 'arrow':
            arrow_idx = payload.get('arrow_index')
            if arrow_idx is None:
                arrow_idx = _parse_arrow_index(payload.get('id'))
            arrow_index = _parse_arrow_index(arrow_idx)
            if arrow_index is not None:
                arrows = json_data.get('arrows', [])
                if 0 <= arrow_index < len(arrows):
                    arrows.pop(arrow_index)
                    modified = True
        elif element_type == 'compound':
            compounds = json_data.get('compound_data')
            if isinstance(compounds, list):
                comp_id = payload.get('compound_id')
                dom_id = payload.get('id')
                def _match(comp, idx):
                    if comp_id is not None and str(comp.get('compound_id')) == str(comp_id):
                        return True
                    if dom_id:
                        generated = f"compound_{comp.get('compound_id') or idx}"
                        return generated == dom_id
                    return False
                new_compounds = [comp for idx, comp in enumerate(compounds) if not _match(comp, idx)]
                if len(new_compounds) != len(compounds):
                    json_data['compound_data'] = new_compounds
                    modified = True
        elif element_type == 'text-box':
            text_blocks = json_data.get('text_data')
            if isinstance(text_blocks, list):
                text_id = payload.get('text_id')
                dom_id = payload.get('id')
                def _match(tb, idx):
                    if text_id is not None and str(tb.get('text_id')) == str(text_id):
                        return True
                    if dom_id:
                        generated = f"text_{tb.get('text_id') or idx}"
                        return generated == dom_id
                    return False
                new_blocks = [tb for idx, tb in enumerate(text_blocks) if not _match(tb, idx)]
                if len(new_blocks) != len(text_blocks):
                    json_data['text_data'] = new_blocks
                    modified = True
        elif element_type == 'ptm':
            protbox_id = payload.get('protbox_id')
            uniprot = payload.get('uniprot')
            ptm_key = payload.get('ptm_key')
            if protbox_id and uniprot and ptm_key:
                pb = _find_protbox_entry(json_data, protbox_id)
                if pb:
                    overrides = pb.setdefault('ptm_overrides', {})
                    uni_map = overrides.setdefault(str(uniprot), {})
                    ptm_map = uni_map.setdefault(ptm_key, {})
                    ptm_map['hidden'] = True
                    modified = True
        if modified:
            json_data_reactive.set(json_data)

    @reactive.Effect
    @reactive.event(input.save_json)
    def save_json():
        json_data = json_data_reactive.get()
        if json_data:
            save_json_data(json_data)



app = App(app_ui, server)

if __name__ == "__main__":
    app.run()






