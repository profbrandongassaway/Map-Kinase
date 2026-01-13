import copy
import json
import os
from pathlib import Path
from typing import Optional
from shiny import App, ui, render, reactive

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR / "output" / "testing_file_001"
DEFAULT_JSON_PATH = Path(
    r"C:\Users\clayt\OneDrive\Desktop\MapKinase\MapKinase_WebApp\hsa04010_pathway_data_20250907_210715.json")


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


def _clone_json(payload):
    try:
        return copy.deepcopy(payload)
    except Exception:
        try:
            return json.loads(json.dumps(payload))
        except Exception:
            return payload


def _build_blank_canvas(catalog_info=None):
    base = {
        'general_data': {'settings': {'show_arrows': True, 'show_text_boxes': True}},
        'protein_data': {},
        'protbox_data': [],
        'groups': [],
        'arrows': [],
        'compound_data': [],
        'text_data': []
    }
    if catalog_info:
        try:
            base['_global_protein_catalog'] = dict(catalog_info)
        except Exception:
            base['_global_protein_catalog'] = catalog_info
    return base

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
    if isinstance(catalog_info, dict) and catalog_info.get('protein_catalog'):
        catalog_data = catalog_info.get('protein_catalog') or {}
    elif catalog_path and os.path.exists(catalog_path):
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
    full_width = bool(json_data.get('_full_width_canvas'))
    if full_width:
        # Ensure ample drawing space so objects are not clipped near the bottom of the viewer.
        min_w, min_h = 1200, 900
        max_x = max(max_x, min_w)
        max_y = max(max_y, min_h)
        container_height_px = max(max_y, min_h)
        container_width_style = "100%"
        container_height_style = f"{container_height_px}px"
        canvas_width_style = "100%"
        canvas_height_style = f"{container_height_px}px"
    else:
        container_width_style = f"{max_x + 20}px"
        container_height_style = f"{max_y + 20}px"
        canvas_width_style = f"{max_x}px"
        canvas_height_style = f"{max_y}px"
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
    (function setupGlobalErrorHooks() {{
        if (window.__mkErrorHooksInstalled) return;
        window.__mkErrorHooksInstalled = true;
        const sendErr = (label, err) => {{
            try {{
                const payload = {{
                    label,
                    message: err?.message || String(err),
                    stack: err?.stack || null,
                    time: Date.now()
                }};
                Shiny?.setInputValue('client_error', payload, {{ priority: 'event' }});
                const dbg = document.getElementById('debug_json');
                if (dbg) dbg.textContent = (dbg.textContent || '') + '\\n' + label + ': ' + payload.message;
            }} catch (e) {{ /* ignore */ }}
        }};
        window.addEventListener('error', (ev) => sendErr('window.error', ev?.error || ev));
        window.addEventListener('unhandledrejection', (ev) => sendErr('unhandledrejection', ev?.reason || ev));
    }})(); 
    var protboxMap = protboxMap || {{}};
    var attachments = attachments || {{}};
    var attachedBy = attachedBy || {{}};
    var protboxSnapPoints = protboxSnapPoints || {{}};
    var protboxHandleDists = protboxHandleDists || {{}};
    var protboxLinks = protboxLinks || {{}};
    var shiftKeyDown = typeof shiftKeyDown === 'boolean' ? shiftKeyDown : false;
    var activeDragProtboxId = typeof activeDragProtboxId === 'string' ? activeDragProtboxId : null;
    var TEXT_HIT_PADDING = typeof TEXT_HIT_PADDING === 'number' ? TEXT_HIT_PADDING : 6;
    var labelDefaults = labelDefaults || {{'N1': [-5, -5, 'right'],'N2': [0, -11, 'center'],'N3': [5, -5, 'left'],'S1': [-3, 5, 'right'],'S2': [0, 12, 'center'],'S3': [3, 5, 'left'],'W1': [-3, -2, 'right'],'W2': [-3, 2, 'right'],'E1': [3, -2, 'left'],'E2': [3, 2, 'left']}};
    var anchorMap = anchorMap || {{'left': 'start','center': 'middle','right': 'end'}};
    var mkHistory = typeof mkHistory === 'object' && mkHistory ? mkHistory : null;
    var activeSnapKey = typeof activeSnapKey === 'string' ? activeSnapKey : null;
    var textHandleGroups = textHandleGroups || {{}};
    var activeTextEditorId = typeof activeTextEditorId === 'string' ? activeTextEditorId : null;
    var isTextEditing = typeof isTextEditing === 'boolean' ? isTextEditing : false;
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
            protboxLinks = {{}};
            autoConnectSelectedProtboxes = null;
            shiftKeyDown = false;
            activeDragProtboxId = null;
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
            const defaultPtmPositionPriority = ['N1','N3','S1','S3','W1','E1','W2','E2','N2','S2'];
            const ksPtmPositionPriority = ['W1','W2','E1','E2','N1','S1','N2','S2','N3','S3'];
            let ptmPositionPriority = defaultPtmPositionPriority;
            let prioritizedPtmPositions = ptmPositionPriority.slice(0, 4);
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
              const order = ptmPositionPriority;
              for (const key of order) {{
                if (snaps[key] && !occupied.has(key)) return {{ key, x: snaps[key].x, y: snaps[key].y }};
              }}
              return null;
            }};
            let data = JSON.parse(document.getElementById('pathway-data').textContent);
            try {{
                const persistKey = data && data._bookmark_key;
                if (persistKey) {{
                    window.__mkPersisted = window.__mkPersisted || {{}};
                    const cached = window.__mkPersisted[persistKey];
                    const incomingToken = data && data._persist_token;
                    const cachedToken = cached && cached._persist_token;
                    if (cached && typeof cached === 'object' && cachedToken && incomingToken && cachedToken === incomingToken) {{
                        const incomingSettings = (data.general_data && data.general_data.settings) || {{}};
                        const cachedSettings = (cached.general_data && cached.general_data.settings) || {{}};
                        cached.general_data = cached.general_data || {{}};
                        cached.general_data.settings = Object.assign({{}}, cachedSettings, incomingSettings);
                        cached._active_fc_index = data._active_fc_index || cached._active_fc_index || 1;
                        if (data._color_preview_override) {{
                            cached._color_preview_override = data._color_preview_override;
                        }}
                        data = cached;
                    }} else if (cached && typeof cached === 'object' && !incomingToken) {{
                        data = cached;
                    }} else {{
                        window.__mkPersisted[persistKey] = data;
                    }}
                }}
            }} catch (persistErr) {{
                console.log('m3: persistence merge failed', persistErr);
            }}
            const bookmarkKey = (data && data._bookmark_key ? String(data._bookmark_key).toLowerCase() : '');
            ptmPositionPriority = bookmarkKey === 'ks' ? ksPtmPositionPriority : defaultPtmPositionPriority;
            prioritizedPtmPositions = ptmPositionPriority.slice(0, 4);
            const preview = data._kegg_preview || {{}};
            const offxVal = Number(preview.offset_x);
            const offyVal = Number(preview.offset_y);
            const offx = Number.isFinite(offxVal) ? offxVal : 0;
            const offy = Number.isFinite(offyVal) ? offyVal : 0;
            // Keep background offsets configurable; foreground layers stay unshifted so all elements align.
            const fgOffsetX = 0;
            const fgOffsetY = 0;
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
                            const offsetX = offx;
                            const offsetY = offy;
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
            const buildExportSnapshot = () => {{
                const exportProtboxes = Array.isArray(protBoxes) ? protBoxes.map(pb => {{
                    if (!pb || typeof pb !== 'object') return pb;
                    let clone = null;
                    try {{
                        clone = JSON.parse(JSON.stringify(pb));
                    }} catch (err) {{
                        clone = Object.assign({{}}, pb);
                    }}
                    const pid = normalizeProtboxId(clone.protbox_id);
                    const mapEntry = (pid && protboxMap && protboxMap[pid]) ? protboxMap[pid] : null;
                    if (mapEntry) {{
                        if (Number.isFinite(mapEntry.x)) clone.x = mapEntry.x;
                        if (Number.isFinite(mapEntry.y)) clone.y = mapEntry.y;
                        if (Number.isFinite(mapEntry.width)) clone.width = mapEntry.width;
                        if (Number.isFinite(mapEntry.height)) clone.height = mapEntry.height;
                    }}
                    const selected = (currentSelected && (currentSelected[clone.protbox_id] || currentSelected[pid])) || null;
                    if (selected) clone.selected_uniprot = selected;
                    return clone;
                }}) : protBoxes;
                const snapshot = {{
                    general_data: {{ settings }},
                    protbox_data: exportProtboxes,
                    protein_data: proteinData,
                    groups,
                    arrows,
                    compound_data: compounds,
                    text_data: textBlocks
                }};
                try {{
                    return JSON.parse(JSON.stringify(snapshot));
                }} catch (err) {{
                    return snapshot;
                }}
            }};
            const exportKey = bookmarkKey || 'default';
            window.__mkExportSnapshotMap = window.__mkExportSnapshotMap || {{}};
            window.__mkExportSnapshotMap[exportKey] = buildExportSnapshot;
            if (window.Shiny && Shiny.addCustomMessageHandler && !window.__mkExportHandlerInstalled) {{
                window.__mkExportHandlerInstalled = true;
                Shiny.addCustomMessageHandler('request_export_snapshot', function(msg) {{
                    const key = (msg && msg.prefix ? String(msg.prefix).toLowerCase() : '');
                    const snapshotFn = (window.__mkExportSnapshotMap && window.__mkExportSnapshotMap[key]) || buildExportSnapshot;
                    const payload = snapshotFn();
                    const prefix = msg && msg.prefix ? msg.prefix : '';
                    Shiny.setInputValue('export_snapshot', {{ prefix, payload, ts: Date.now() }}, {{ priority: 'event' }});
                }});
            }}
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
            const debugMode = !!settings.debug_mode;
            const activeFcIndex = Number(data._active_fc_index || 1);
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
            let figureKeyCounter = 0;
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
            const getSelectedProtboxes = () => {{
                const ids = new Set();
                selectionMap.forEach(entry => {{
                    if (!entry) return;
                    if (entry.type === 'prot-box') {{
                        ids.add(entry.id);
                    }} else if (entry.type === 'group' && typeof collectProtboxesForGroup === 'function') {{
                        collectProtboxesForGroup(entry.id).forEach(pid => ids.add(pid));
                    }}
                }});
                return Array.from(ids);
            }};
            const setArrowDashState = (arrowId, dashedOn) => {{
                const arrow = getArrowById(arrowId);
                if (!arrow) return;
                if (arrow.line === 'dashed_arrow') {{
                    arrow.line = 'arrow';
                }}
                arrow.dashed = !!dashedOn;
                updateArrowVisual(arrowId);
                translateArrowByDelta(arrowId, 0, 0);
                updateAttachedArrows(arrowId, 'start');
                updateAttachedArrows(arrowId, 'end');
            }};
            const findArrowBetweenProtboxes = (idA, idB) => {{
                const a = normalizeProtboxId(idA);
                const b = normalizeProtboxId(idB);
                if (!a || !b) return null;
                for (let i = 0; i < arrows.length; i++) {{
                    const ar = arrows[i];
                    if (!ar) continue;
                    const p1 = normalizeProtboxId(ar.protbox_id_1);
                    const p2 = normalizeProtboxId(ar.protbox_id_2);
                    if ((p1 === a && p2 === b) || (p1 === b && p2 === a)) {{
                        return arrowIdFromIndex(i);
                    }}
                }}
                return null;
            }};
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
                const str = String(elementId);
                let proto = null;
                let core = str;
                if (str.includes('__')) {{
                    const parts = str.split('__');
                    proto = parts.shift() || null;
                    core = parts.join('__');
                }}
                const parts = core.split('_');
                if (parts.length < 3) return null;
                const suffix = parts.pop();
                const uniprot = parts.shift();
                const ptmKey = parts.join('_');
                if (!uniprot || !ptmKey) return null;
                return {{ uniprot, ptmKey, suffix, protboxId: proto }};
            }};
            const ptmElementId = (protboxId, uniprot, ptmKey, suffix) => {{
                const pid = normalizeProtboxId(protboxId) || 'pb';
                return `${{pid}}__${{uniprot}}_${{ptmKey}}_${{suffix}}`;
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
            const parseFoldChange = (entity, idx = activeFcIndex) => {{
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
            const entityColor = (entity, idx = activeFcIndex) => {{
                const foldVal = parseFoldChange(entity, idx);
                if (foldVal !== null) {{
                    return gradientColorFromFold(foldVal);
                }}
                return toRgbString(entity && entity[`fc_color_${{idx}}`]);
            }};
            const trackableMoveTypes = new Set(['prot-box','ptm-shape','ptm-label','ptm-symbol','compound','text-box','figure-key','arrow','arrow-start','arrow-end']);
            const toCoordinateNumber = (value, fallback = 0) => {{
                const num = Number(value);
                return Number.isFinite(num) ? num : fallback;
            }};
            const ensureElementForId = (elementId) => {{
                if (!elementId) return null;
                return draw.find(`[data-id="${{elementId}}"]`)[0] || null;
            }};
            const applyPositionForEntry = (entry, which = 'before') => {{
                if (!entry) return;
                const pos = entry[which];
                if (!pos) return;
                const el = ensureElementForId(entry.id);
                if (!el) return;
                if (entry.type === 'prot-box') {{
                    const rect = el;
                    const curX = toCoordinateNumber(rect.x && rect.x());
                    const curY = toCoordinateNumber(rect.y && rect.y());
                    const dx = pos.x - curX;
                    const dy = pos.y - curY;
                    const handleVisit = new Set();
                    const delta = translateProtbox(entry.id, dx, dy, {{ rect, skipRectMove: false }}, handleVisit) || {{ dx, dy }};
                    propagateLinkedMove(entry.id, delta.dx, delta.dy, new Set([entry.id]), handleVisit);
                    cleanupBrokenLinks();
                    refreshGroupsForProtbox(entry.id);
                    if (entry.attached && entry.attached[which]) {{
                        reattachArrowPayloads(entry.attached[which]);
                        rebuildAttachmentsIndex();
                        updateArrowPositions(entry.id, 'North');
                        updateArrowPositions(entry.id, 'South');
                        updateArrowPositions(entry.id, 'East');
                        updateArrowPositions(entry.id, 'West');
                    }}
                    return;
                }}
                if (entry.type === 'text-box') {{
                    const tb = findTextBlockByDomId(entry.id);
                    if (tb) {{
                        tb.x = toCoordinateNumber(pos.x, tb.x);
                        tb.y = toCoordinateNumber(pos.y, tb.y);
                        applyTextBoxLayout(entry.id, tb);
                        return;
                    }}
                }}
                if (entry.type === 'ptm-shape' || entry.type === 'ptm-label' || entry.type === 'ptm-symbol') {{
                    if (entry.type === 'ptm-shape') {{
                        if (pos.posKey) {{
                            el.attr && el.attr('data-pos-key', pos.posKey);
                        }} else {{
                            try {{ el.node.removeAttribute('data-pos-key'); }} catch (ignore) {{}}
                        }}
                    }}
                }}
                if (typeof el.x === 'function' && typeof el.y === 'function') {{
                    el.move(pos.x, pos.y);
                }} else if (typeof el.cx === 'function') {{
                    el.cx(pos.x);
                    el.cy(pos.y);
                }}
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
                menu.style.padding = '2px 0';
                menu.style.fontSize = '12px';
                menu.style.zIndex = '1000';
                const removeItem = document.createElement('div');
                removeItem.textContent = 'Remove';
                removeItem.style.padding = '4px 8px';
                removeItem.style.cursor = 'pointer';
                const cleanupMenu = () => {{
                    if (menu.parentNode) menu.remove();
                    document.removeEventListener('click', cleanupMenu);
                }};
                removeItem.addEventListener('click', () => {{
                    deletePtmFromProtbox(protboxId, meta);
                    cleanupMenu();
                }});
                const isKsActive = !!document.querySelector('#bookmark_selector a.nav-link.active[data-value="ks"]');
                if (isKsActive) {{
                    const showKinItem = document.createElement('div');
                    showKinItem.textContent = "Show Substrate's Kinases";
                    showKinItem.style.padding = '4px 8px';
                    showKinItem.style.cursor = 'pointer';
                    showKinItem.addEventListener('click', () => {{
                        try {{
                            Shiny.setInputValue('ks_spawn_ptm_kinases', {{ protbox_id: protboxId, uniprot: meta.uniprot, ptm_key: meta.ptmKey, ts: Date.now() }}, {{ priority: 'event' }});
                        }} catch (e) {{}}
                        cleanupMenu();
                    }});
                    menu.appendChild(showKinItem);
                }}
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
                const dtype = target.getAttribute ? target.getAttribute('data-type') : null;
                if (dtype && dtype.startsWith('group-')) return false;
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
            const selectionMap = new Map();
            let primarySelectionKey = null;
            let activeGroupEditId = null;
            const groupMap = {{}};
            const groupVisuals = {{}};
            const groupBackgrounds = {{}};
            let groupMembership = {{}};
            let groupMembershipAll = {{}};
            let groupParents = {{}};
            let groupCounter = 0;
            let groupEditIndicator = null;
            const entryKey = (type, id) => `${{type || ''}}::${{id || ''}}`;
            const protboxGroups = (pid) => {{
                const norm = normalizeProtboxId(pid);
                return norm ? groupMembership[norm] : null;
            }};
            const protboxInAnyGroup = (pid) => {{
                const gs = protboxGroups(pid);
                return !!(gs && gs.size);
            }};
            const sharesGroup = (a, b) => {{
                const na = normalizeProtboxId(a);
                const nb = normalizeProtboxId(b);
                if (!na || !nb) return true;
                const gA = groupMembership[na];
                const gB = groupMembership[nb];
                if (!gA || !gA.size) return true;
                if (!gB || !gB.size) return false;
                for (const gid of gA) {{
                    if (gB.has(gid)) return true;
                }}
                return false;
            }};
            const groupsForEntry = (type, id) => {{
                const key = entryKey(type, id);
                return groupMembershipAll[key] || null;
            }};
            const resolveRootGroup = (groupId) => {{
                let current = `${{groupId}}`;
                const visited = new Set();
                while (current && groupParents[current] && groupParents[current].size) {{
                    if (visited.has(current)) break;
                    visited.add(current);
                    const parent = Array.from(groupParents[current])[0];
                    current = `${{parent}}`;
                }}
                return current;
            }};
            const isGroupDescendant = (gid, ancestorId) => {{
                const target = `${{ancestorId}}`;
                const stack = [`${{gid}}`];
                const visited = new Set();
                while (stack.length) {{
                    const cur = stack.pop();
                    if (!cur || visited.has(cur)) continue;
                    visited.add(cur);
                    if (cur === target) return true;
                    const parents = groupParents[cur];
                    if (parents && parents.size) {{
                        parents.forEach(p => stack.push(`${{p}}`));
                    }}
                }}
                return false;
            }};
            const allMembersForEntry = (type, id) => {{
                const gs = groupsForEntry(type, id);
                if (!gs || !gs.size) return [{{ type, id }}];
                const acc = new Set();
                gs.forEach(gid => {{
                    collectMembersForGroup(gid).forEach(m => acc.add(entryKey(m.type, m.id)));
                }});
                return Array.from(acc).map(key => {{
                    const [t, ...rest] = key.split('::');
                    return {{ type: t, id: rest.join('::') }};
                }});
            }};
            const groupBackgroundLayer = draw.group();
            const foregroundGroup = draw.group();
            const handleGroup = foregroundGroup.group();
            const shapeGroup = foregroundGroup.group(); // shapes sit beneath arrows/protboxes
            const arrowGroup = foregroundGroup.group();
            const protboxGroup = foregroundGroup.group();
            const groupOverlay = foregroundGroup.group();
            const figureKeyGroup = foregroundGroup.group();
            const compoundGroup = draw.group();  
            const textGroup = foregroundGroup.group();     
            const alignmentGuideGroup = foregroundGroup.group();
            let selectionBoxActive = false;
            let selectionBoxRect = null;
            let selectionBoxStart = null;
            let selectionBoxEnd = null;
            let selectionBoxStartScreen = null;
            let selectionBoxEndScreen = null;
            let suppressNextBackgroundClick = false;
            const normalizeBounds = (a, b) => {{
                if (!a || !b) return null;
                const left = Math.min(a.x, b.x);
                const right = Math.max(a.x, b.x);
                const top = Math.min(a.y, b.y);
                const bottom = Math.max(a.y, b.y);
                return {{ left, right, top, bottom }};
            }};
            const getScreenBounds = (el) => {{
                if (!el || !el.node || typeof el.node.getBoundingClientRect !== 'function') return null;
                const b = el.node.getBoundingClientRect();
                return normalizeBounds({{ x: b.left, y: b.top }}, {{ x: b.right, y: b.bottom }});
            }};
            const getSvgBounds = (el) => {{
                if (!el) return null;
                try {{
                    const b = el.rbox ? el.rbox(draw) : null;
                    if (b && typeof b.x === 'number') {{
                        return normalizeBounds({{ x: b.x, y: b.y }}, {{ x: b.x + b.width, y: b.y + b.height }});
                    }}
                }} catch (err) {{}}
                try {{
                    const b = el.bbox ? el.bbox() : null;
                    if (b && typeof b.x === 'number') {{
                        return normalizeBounds({{ x: b.x, y: b.y }}, {{ x: b.x + b.width, y: b.y + b.height }});
                    }}
                }} catch (err) {{}}
                return null;
            }};
            const findTextBlockByDomId = (id) => textBlocks.find(tb => tb && (tb._client_id === id || `text_${{tb.text_id}}` === id));
            const baseTextStyle = () => ({{
                fontFamily: settings?.textbox_label_font || 'Arial',
                fontSize: settings?.textbox_label_size || 14,
                color: '#000000',
                align: 'center',
                vertical: 'center',
                bold: false,
                italic: false,
                underline: false
            }});
            const resolveTextStyle = (textBlock) => {{
                const merged = Object.assign({{}}, baseTextStyle(), textBlock?.text_style || {{}});
                if (textBlock?.fgcolor) merged.color = textBlock.fgcolor;
                return merged;
            }};
            const applyStyleToEditor = (editor, style) => {{
                if (!editor || !style) return;
                editor.style.fontFamily = style.fontFamily || baseTextStyle().fontFamily;
                editor.style.fontSize = `${{style.fontSize || baseTextStyle().fontSize}}px`;
                editor.style.color = style.color || baseTextStyle().color;
                editor.style.textAlign = style.align || 'left';
                editor.style.fontWeight = style.bold ? 'bold' : 'normal';
                editor.style.fontStyle = style.italic ? 'italic' : 'normal';
                editor.style.textDecoration = style.underline ? 'underline' : 'none';
                editor.style.display = 'flex';
                editor.style.flexDirection = 'column';
                editor.style.width = '100%';
                editor.style.height = '100%';
                const vMap = {{ top: 'flex-start', center: 'center', bottom: 'flex-end' }};
                editor.style.justifyContent = vMap[style.vertical] || 'flex-start';
                const hMap = {{ left: 'flex-start', center: 'center', right: 'flex-end' }};
                editor.style.alignItems = hMap[style.align] || 'flex-start';
            }};
            const selectionInsideEditor = (editor) => {{
                const sel = window.getSelection();
                return !!(editor && sel && sel.rangeCount && editor.contains(sel.anchorNode));
            }};
            const ensureSelectionForMenu = (editor) => {{
                if (!editor) return;
                const sel = window.getSelection();
                if (!sel || sel.rangeCount === 0 || sel.isCollapsed || !editor.contains(sel.anchorNode)) {{
                    const range = document.createRange();
                    range.selectNodeContents(editor);
                    sel.removeAllRanges();
                    sel.addRange(range);
                }}
            }};
            const applyStyleToSelection = (editor, styleMap) => {{
                const sel = window.getSelection();
                if (!selectionInsideEditor(editor) || !sel || sel.isCollapsed) return false;
                const range = sel.getRangeAt(0);
                const span = document.createElement('span');
                Object.assign(span.style, styleMap || {{}});
                span.appendChild(range.extractContents());
                range.insertNode(span);
                sel.removeAllRanges();
                const newRange = document.createRange();
                newRange.selectNodeContents(span);
                sel.addRange(newRange);
                return true;
            }};
            const positionTextHandles = (id, textBlock) => {{
                const handles = textHandleGroups[id];
                if (!handles || !handles.handles) return;
                const tb = textBlock || findTextBlockByDomId(id);
                if (!tb) return;
                const {{ x = 0, y = 0, width = 60, height = 20 }} = tb;
                const dirs = Object.keys(handles.handles);
                dirs.forEach(dir => {{
                    const handle = handles.handles[dir];
                    if (!handle) return;
                    let hx = x, hy = y;
                    if (dir.includes('e')) hx = x + width;
                    if (dir.includes('s')) hy = y + height;
                    if (dir === 'n' || dir === 's') hx = x + width / 2;
                    if (dir === 'e' || dir === 'w') hy = y + height / 2;
                    handle.center(hx, hy);
                }});
            }};
            const removeTextHandles = (id) => {{
                const entry = textHandleGroups[id];
                if (!entry) return;
                Object.values(entry.handles || {{}}).forEach(h => {{
                    try {{ h.remove(); }} catch (err) {{}}
                }});
                delete textHandleGroups[id];
            }};
            function beginTextResize(evt, id, dir) {{
                const tb = findTextBlockByDomId(id);
                if (!tb) return;
                evt.preventDefault();
                evt.stopPropagation();
                if (isTextEditing && activeTextEditorId === id) {{
                    exitTextEditMode(true);
                }}
                const start = draw.point(evt.clientX, evt.clientY);
                const startGeom = {{ x: tb.x || 0, y: tb.y || 0, width: tb.width || 60, height: tb.height || 20 }};
                const minSize = 30;
                const onMove = (moveEvt) => {{
                    const pt = draw.point(moveEvt.clientX, moveEvt.clientY);
                    let dx = pt.x - start.x;
                    let dy = pt.y - start.y;
                    let {{ x, y, width, height }} = startGeom;
                    if (dir.includes('w')) {{ x = startGeom.x + dx; width = startGeom.width - dx; }}
                    if (dir.includes('n')) {{ y = startGeom.y + dy; height = startGeom.height - dy; }}
                    if (dir.includes('e')) {{ width = startGeom.width + dx; }}
                    if (dir.includes('s')) {{ height = startGeom.height + dy; }}
                    width = Math.max(minSize, width);
                    height = Math.max(minSize, height);
                    if (dir.includes('w')) {{ x = startGeom.x + (startGeom.width - width); }}
                    if (dir.includes('n')) {{ y = startGeom.y + (startGeom.height - height); }}
                    tb.x = x;
                    tb.y = y;
                    tb.width = width;
                    tb.height = height;
                    applyTextBoxLayout(id, tb);
                }};
                const onUp = () => {{
                    document.removeEventListener('mousemove', onMove);
                    document.removeEventListener('mouseup', onUp);
                    Shiny?.setInputValue('element_moved', {{ type: 'text-box', id, x: tb.x, y: tb.y, width: tb.width, height: tb.height }}, {{ priority: 'event' }});
                    persistTextBlockChange(tb, 'resize');
                }};
                document.addEventListener('mousemove', onMove);
                document.addEventListener('mouseup', onUp);
            }}
            const ensureTextHandles = (id, textBlock) => {{
                if (textHandleGroups[id]) {{
                    positionTextHandles(id, textBlock);
                    handleGroup.front();
                    Object.values(textHandleGroups[id].handles || {{}}).forEach(h => h.front());
                    return textHandleGroups[id];
                }}
                const getEditorFor = () => {{
                    const grp = draw.findOne(`[data-id="${{id}}"]`);
                    return grp ? grp.node.querySelector('.mk-text-editor') : null;
                }};
                const handleDirs = ['nw','n','ne','e','se','s','sw','w'];
                const handles = {{}};
                handleDirs.forEach(dir => {{
                    const h = handleGroup.rect(8, 8).center(textBlock?.x || 0, textBlock?.y || 0).addClass('text-resize-handle').attr({{ 'data-type': 'text-resize', 'data-handle': dir, 'data-owner': id }}).hide();
                    const cursorMap = {{ n: 'ns-resize', s: 'ns-resize', e: 'ew-resize', w: 'ew-resize', ne: 'nesw-resize', sw: 'nesw-resize', nw: 'nwse-resize', se: 'nwse-resize' }};
                    h.node.style.cursor = cursorMap[dir] || 'pointer';
                    h.on('mousedown', (evt) => beginTextResize(evt, id, dir));
                    h.node.addEventListener('contextmenu', (evt) => {{
                        evt.preventDefault();
                        evt.stopPropagation();
                        showTextBoxMenu(evt, id, getEditorFor());
                    }});
                    handles[dir] = h;
                }});
                textHandleGroups[id] = {{ handles }};
                handleGroup.front();
                Object.values(handles).forEach(h => h.front());
                positionTextHandles(id, textBlock);
                return textHandleGroups[id];
            }};
            const toggleTextHandles = (id, show = false) => {{
                const entry = textHandleGroups[id];
                if (!entry) return;
                Object.values(entry.handles || {{}}).forEach(h => {{
                    if (show) {{
                        h.show();
                        h.front();
                    }} else {{
                        h.hide();
                    }}
                }});
                if (show) handleGroup.front();
            }};
            const applyTextBoxLayout = (id, textBlock) => {{
                const group = draw.findOne(`[data-id="${{id}}"]`);
                if (!group) return;
                const rect = group.findOne('[data-role="text-rect"]');
                const fo = group.findOne('[data-role="text-fo"]');
                const hit = group.findOne('[data-role="text-hit"]');
                const outline = group.findOne('[data-role="text-outline"]');
                const w = textBlock.width || 60;
                const h = textBlock.height || 20;
                const x = textBlock.x || 0;
                const y = textBlock.y || 0;
                if (rect) rect.size(w, h).move(x, y);
                if (fo) fo.size(w, h).move(x, y);
                if (outline) outline.size(w, h).move(x, y);
                if (hit) hit.size(w + TEXT_HIT_PADDING * 2, h + TEXT_HIT_PADDING * 2).move(x - TEXT_HIT_PADDING, y - TEXT_HIT_PADDING);
                positionTextHandles(id, textBlock);
            }};
            const persistTextBlockChange = (tb, reason = 'edit') => {{
                if (!tb) return;
                const payload = {{
                    type: 'text-box',
                    id: tb._client_id,
                    text_id: tb.text_id,
                    x: tb.x,
                    y: tb.y,
                    width: tb.width,
                    height: tb.height,
                    label: tb.label,
                    html: tb.html,
                    text_style: tb.text_style,
                    bgcolor: tb.bgcolor,
                    fgcolor: tb.fgcolor,
                    border_color: tb.border_color,
                    reason,
                    timestamp: Date.now()
                }};
                Shiny?.setInputValue('text_box_changed', payload, {{ priority: 'event' }});
            }};
            const syncTextContentToBlock = (tb, editor) => {{
                if (!tb || !editor) return tb;
                tb.html = editor.innerHTML;
                tb.label = editor.textContent || '';
                tb.text_style = resolveTextStyle(tb);
                return tb;
            }};
            const enterTextEditMode = (id, opts = {{}}) => {{
                const group = draw.findOne(`[data-id="${{id}}"]`);
                if (!group) return;
                const editor = group.node.querySelector('.mk-text-editor');
                const tb = findTextBlockByDomId(id);
                if (!editor || !tb) return;
                if (activeTextEditorId && activeTextEditorId !== id) {{
                    exitTextEditMode(true);
                }}
                isTextEditing = true;
                activeTextEditorId = id;
                editor.contentEditable = true;
                editor.focus();
                const range = document.createRange();
                const sel = window.getSelection();
                if (opts.selectAll) {{
                    range.selectNodeContents(editor);
                }} else {{
                    range.selectNodeContents(editor);
                    range.collapse(false);
                }}
                sel.removeAllRanges();
                sel.addRange(range);
                applyStyleToEditor(editor, resolveTextStyle(tb));
            }};
            const exitTextEditMode = (save = true) => {{
                if (!activeTextEditorId) return;
                const group = draw.findOne(`[data-id="${{activeTextEditorId}}"]`);
                const tb = findTextBlockByDomId(activeTextEditorId);
                const editor = group?.node.querySelector('.mk-text-editor');
                if (editor) {{
                    if (save && tb) {{
                        syncTextContentToBlock(tb, editor);
                        persistTextBlockChange(tb, 'edit');
                    }}
                    editor.contentEditable = false;
                }}
                isTextEditing = false;
                activeTextEditorId = null;
            }};
            const showTextFormatMenu = (evt, id, editor) => {{
                evt.preventDefault();
                evt.stopPropagation();
                removeExistingContextMenu();
                if (!editor) return;
                const sel = window.getSelection();
                if (!selectionInsideEditor(editor) || !sel || sel.isCollapsed) return;
                const tb = findTextBlockByDomId(id);
                if (!tb) return;
                const style = resolveTextStyle(tb);
                const menu = document.createElement('div');
                menu.className = 'context-menu text-context-menu';
                menu.style.position = 'fixed';
                menu.style.left = `${{evt.clientX + menuOffsetX}}px`;
                menu.style.top = `${{evt.clientY + menuOffsetY}}px`;
                menu.style.backgroundColor = 'white';
                menu.style.border = '1px solid #ccc';
                menu.style.padding = '4px 6px';
                menu.style.fontSize = '12px';
                menu.style.zIndex = '1200';
                menu.style.minWidth = '180px';
                const applyAndPersist = (styleMap, patch = null) => {{
                    if (!applyStyleToSelection(editor, styleMap || {{}})) {{
                        tb.text_style = Object.assign(resolveTextStyle(tb), patch || styleMap || {{}});
                        applyStyleToEditor(editor, tb.text_style);
                    }}
                    if ((patch && patch.color) || (styleMap && styleMap.color)) {{
                        tb.fgcolor = patch?.color || styleMap?.color;
                    }}
                    syncTextContentToBlock(tb, editor);
                    persistTextBlockChange(tb, 'format');
                }};
                const ul = document.createElement('ul');
                ul.style.listStyle = 'none';
                ul.style.margin = '0';
                ul.style.padding = '0';
                const makeSubmenu = (parentLi) => {{
                    const sub = document.createElement('div');
                    sub.className = 'context-submenu';
                    sub.style.position = 'absolute';
                    sub.style.left = '100%';
                    sub.style.top = '0';
                    sub.style.backgroundColor = 'white';
                    sub.style.border = '1px solid #ccc';
                    sub.style.padding = '6px';
                    sub.style.margin = '0';
                    sub.style.display = 'none';
                    sub.style.minWidth = '180px';
                    sub.style.boxShadow = '0 2px 6px rgba(0,0,0,0.15)';
                    parentLi.style.position = 'relative';
                    parentLi.addEventListener('mouseenter', () => sub.style.display = 'block');
                    parentLi.addEventListener('mouseleave', () => sub.style.display = 'none');
                    parentLi.appendChild(sub);
                    return sub;
                }};
                // Font submenu with regex search
                const liFont = document.createElement('li');
                liFont.textContent = 'Font';
                liFont.style.padding = '4px 8px';
                liFont.style.cursor = 'pointer';
                const fontSub = makeSubmenu(liFont);
                const fontSearch = document.createElement('input');
                fontSearch.type = 'text';
                fontSearch.placeholder = 'Search fonts (regex)';
                fontSearch.style.width = '100%';
                fontSearch.style.boxSizing = 'border-box';
                fontSearch.style.marginBottom = '6px';
                fontSearch.addEventListener('click', (e) => {{
                    e.stopPropagation();
                }});
                fontSub.appendChild(fontSearch);
                const fontList = document.createElement('ul');
                fontList.style.listStyle = 'none';
                fontList.style.margin = '0';
                fontList.style.padding = '0';
                fontList.style.maxHeight = '200px';
                fontList.style.overflowY = 'auto';
                fontSub.appendChild(fontList);
                const fonts = ['Arial','Calibri','Helvetica','Georgia','Times New Roman','Courier New','Verdana','Tahoma','Trebuchet MS','Garamond','Palatino','Comic Sans MS'];
                const renderFonts = () => {{
                    fontList.innerHTML = '';
                    let regex = null;
                    try {{
                        const val = fontSearch.value.trim();
                        regex = val ? new RegExp(val, 'i') : null;
                    }} catch (err) {{ regex = null; }}
                    fonts.forEach(f => {{
                        if (regex && !regex.test(f)) return;
                        const item = document.createElement('li');
                        item.textContent = f;
                        item.style.padding = '4px 6px';
                        item.style.cursor = 'pointer';
                        item.style.fontFamily = f;
                        item.addEventListener('mousedown', (e) => e.stopPropagation());
                        item.addEventListener('click', (e) => {{
                            e.stopPropagation();
                            applyAndPersist({{ fontFamily: f }}, {{ fontFamily: f }});
                            fontSub.style.display = 'none';
                            menu.remove();
                        }});
                        fontList.appendChild(item);
                    }});
                }};
                fontSearch.addEventListener('input', renderFonts);
                fontSearch.addEventListener('focus', renderFonts);
                renderFonts();
                ul.appendChild(liFont);
                // Color submenu: color wheel + hex
                const liColor = document.createElement('li');
                liColor.textContent = 'Color';
                liColor.style.padding = '4px 8px';
                liColor.style.cursor = 'pointer';
                const colorSub = makeSubmenu(liColor);
                const colorRow = document.createElement('div');
                colorRow.style.display = 'flex';
                colorRow.style.alignItems = 'center';
                colorRow.style.gap = '6px';
                colorRow.style.marginBottom = '6px';
                const currentColor = (() => {{
                    try {{ return SVG.Color.toHex(style.color || '#000'); }} catch (e) {{ return '#000000'; }}
                }})();
                const colorInput = document.createElement('input');
                colorInput.type = 'color';
                colorInput.value = currentColor;
                const hexInput = document.createElement('input');
                hexInput.type = 'text';
                hexInput.value = currentColor;
                hexInput.maxLength = 7;
                hexInput.style.width = '80px';
                const applyHex = () => {{
                    let val = hexInput.value.trim();
                    if (!val.startsWith('#')) val = '#' + val;
                    if (/^#([0-9a-fA-F]{{3}}|[0-9a-fA-F]{{6}})$/.test(val)) {{
                        colorInput.value = val;
                        applyAndPersist({{ color: val }}, {{ color: val }});
                    }}
                }};
                colorInput.addEventListener('input', () => {{
                    hexInput.value = colorInput.value;
                    applyAndPersist({{ color: colorInput.value }}, {{ color: colorInput.value }});
                }});
                hexInput.addEventListener('change', applyHex);
                hexInput.addEventListener('keydown', (e) => {{
                    if (e.key === 'Enter') applyHex();
                }});
                [colorInput, hexInput].forEach(el => el.addEventListener('mousedown', (e) => e.stopPropagation()));
                colorRow.append(colorInput, hexInput);
                colorSub.appendChild(colorRow);
                const colorPalette = document.createElement('div');
                colorPalette.style.display = 'flex';
                colorPalette.style.flexWrap = 'wrap';
                colorPalette.style.gap = '6px';
                ['#000000','#333333','#666666','#999999','#FFFFFF','#AA0000','#D35400','#F1C40F','#27AE60','#2980B9','#8E44AD'].forEach(c => {{
                    const sw = document.createElement('div');
                    sw.style.width = '18px';
                    sw.style.height = '18px';
                    sw.style.border = c === '#FFFFFF' ? '1px solid #ccc' : '1px solid #666';
                    sw.style.background = c;
                    sw.style.cursor = 'pointer';
                    sw.addEventListener('mousedown', (e) => e.stopPropagation());
                    sw.addEventListener('click', (e) => {{
                        e.stopPropagation();
                        colorInput.value = c;
                        hexInput.value = c;
                        applyAndPersist({{ color: c }}, {{ color: c }});
                        colorSub.style.display = 'none';
                        menu.remove();
                    }});
                    colorPalette.appendChild(sw);
                }});
                colorSub.appendChild(colorPalette);
                ul.appendChild(liColor);
                // Size submenu: textbox + +/- buttons
                const liSize = document.createElement('li');
                liSize.textContent = 'Size';
                liSize.style.padding = '4px 8px';
                liSize.style.cursor = 'pointer';
                const sizeSub = makeSubmenu(liSize);
                const sizeRow = document.createElement('div');
                sizeRow.style.display = 'flex';
                sizeRow.style.alignItems = 'center';
                sizeRow.style.gap = '6px';
                const sizeInput = document.createElement('input');
                sizeInput.type = 'text';
                sizeInput.style.width = '70px';
                sizeInput.value = Number(style.fontSize) || baseTextStyle().fontSize || 14;
                const clampSize = (val) => Math.max(6, Math.min(200, val));
                const commitSize = (val) => {{
                    const num = Number(val);
                    if (!Number.isFinite(num)) return;
                    const clamped = clampSize(num);
                    sizeInput.value = clamped;
                    applyAndPersist({{ fontSize: `${{clamped}}px` }}, {{ fontSize: clamped }});
                }};
                const btnDec = document.createElement('button');
                btnDec.textContent = '-';
                btnDec.style.width = '32px';
                const btnInc = document.createElement('button');
                btnInc.textContent = '+';
                btnInc.style.width = '32px';
                [sizeInput, btnDec, btnInc].forEach(el => el.addEventListener('mousedown', (e) => e.stopPropagation()));
                sizeInput.addEventListener('change', () => commitSize(sizeInput.value));
                sizeInput.addEventListener('keydown', (e) => {{ if (e.key === 'Enter') commitSize(sizeInput.value); }});
                btnDec.addEventListener('click', () => commitSize((Number(sizeInput.value) || baseTextStyle().fontSize) - 1));
                btnInc.addEventListener('click', () => commitSize((Number(sizeInput.value) || baseTextStyle().fontSize) + 1));
                sizeRow.append(sizeInput, btnDec, btnInc);
                sizeSub.appendChild(sizeRow);
                ul.appendChild(liSize);
                // Highlight submenu
                const liHighlight = document.createElement('li');
                liHighlight.textContent = 'Highlight';
                liHighlight.style.padding = '4px 8px';
                liHighlight.style.cursor = 'pointer';
                const hlSub = makeSubmenu(liHighlight);
                const highlights = [
                    ['None', 'transparent'],
                    ['Yellow', '#FFF59D'],
                    ['Cyan', '#BBDEFB'],
                    ['Pink', '#F8BBD0'],
                    ['Orange', '#FFE0B2']
                ];
                highlights.forEach(([label,color]) => {{
                    const btn = document.createElement('div');
                    btn.textContent = label;
                    btn.style.padding = '4px 6px';
                    btn.style.cursor = 'pointer';
                    btn.addEventListener('mousedown', (e) => e.stopPropagation());
                    btn.addEventListener('click', (e) => {{
                        e.stopPropagation();
                        applyAndPersist({{ backgroundColor: color === 'transparent' ? 'transparent' : color }}, {{ highlight: color }});
                        hlSub.style.display = 'none';
                        menu.remove();
                    }});
                    hlSub.appendChild(btn);
                }});
                ul.appendChild(liHighlight);
                menu.appendChild(ul);
                // Inline style toggles (B / I / U)
                const toggleRow = document.createElement('div');
                toggleRow.style.display = 'flex';
                toggleRow.style.gap = '6px';
                toggleRow.style.padding = '6px 4px 0 4px';
                const makeToggle = (key, label, styleMapOn, patchKey) => {{
                    const btn = document.createElement('button');
                    btn.type = 'button';
                    btn.textContent = label;
                    btn.style.flex = '1';
                    btn.style.padding = '6px';
                    btn.style.cursor = 'pointer';
                    btn.style.border = '1px solid #ccc';
                    btn.style.background = '#f8f8f8';
                    btn.style.fontSize = '12px';
                    if (key === 'bold') btn.style.fontWeight = 'bold';
                    if (key === 'italic') btn.style.fontStyle = 'italic';
                    if (key === 'underline') btn.style.textDecoration = 'underline';
                    const currentOn = !!style[patchKey];
                    if (currentOn) btn.classList.add('active');
                    btn.addEventListener('mousedown', (e) => e.stopPropagation());
                    btn.addEventListener('click', (e) => {{
                        e.stopPropagation();
                        const nextOn = !btn.classList.contains('active');
                        btn.classList.toggle('active', nextOn);
                        const styleMap = nextOn ? styleMapOn : (key === 'underline' ? {{ textDecoration: 'none' }} : key === 'italic' ? {{ fontStyle: 'normal' }} : {{ fontWeight: 'normal' }});
                        const patch = {{ [patchKey]: nextOn }};
                        applyAndPersist(styleMap, patch);
                    }});
                    return btn;
                }};
                toggleRow.append(
                    makeToggle('bold', 'B', {{ fontWeight: 'bold' }}, 'bold'),
                    makeToggle('italic', 'I', {{ fontStyle: 'italic' }}, 'italic'),
                    makeToggle('underline', 'U', {{ textDecoration: 'underline' }}, 'underline')
                );
                menu.appendChild(toggleRow);
                document.body.appendChild(menu);
                const onDocDown = (e) => {{
                    if (!menu.contains(e.target)) {{
                        try {{ menu.remove(); }} catch (err) {{}}
                        document.removeEventListener('mousedown', onDocDown, true);
                    }}
                }};
                document.addEventListener('mousedown', onDocDown, true);
            }};
            const showTextBoxMenu = (evt, id, editor, isShape = false) => {{
                evt.preventDefault();
                evt.stopPropagation();
                removeExistingContextMenu();
                const tb = findTextBlockByDomId(id);
                if (!tb) return;
                const menu = document.createElement('div');
                menu.className = 'context-menu text-context-menu';
                menu.style.position = 'fixed';
                menu.style.left = `${{evt.clientX + menuOffsetX}}px`;
                menu.style.top = `${{evt.clientY + menuOffsetY}}px`;
                menu.style.backgroundColor = 'white';
                menu.style.border = '1px solid #ccc';
                menu.style.padding = '4px 6px';
                menu.style.fontSize = '12px';
                menu.style.zIndex = '1200';
                menu.style.minWidth = '180px';
                const updateBoxAppearance = () => {{
                    const group = draw.findOne(`[data-id="${{id}}"]`);
                    if (group) {{
                        const rect = group.findOne('[data-role="text-rect"]');
                        if (rect) rect.fill(tb.bgcolor || 'transparent');
                        const outline = group.findOne('[data-role="text-outline"]');
                        if (outline) outline.stroke({{ color: tb.border_color || 'transparent', width: tb.border_width || 1 }});
                        const ed = editor || group.node.querySelector('.mk-text-editor');
                        if (ed) {{
                            ed.style.background = tb.bgcolor || 'transparent';
                            applyStyleToEditor(ed, resolveTextStyle(tb));
                        }}
                    }}
                    applyTextBoxLayout(id, tb);
                    persistTextBlockChange(tb, 'format');
                }};
                const ul = document.createElement('ul');
                ul.style.listStyle = 'none';
                ul.style.margin = '0';
                ul.style.padding = '0';
                const makeSubmenu = (parentLi) => {{
                    const sub = document.createElement('div');
                    sub.className = 'context-submenu';
                    sub.style.position = 'absolute';
                    sub.style.left = '100%';
                    sub.style.top = '0';
                    sub.style.backgroundColor = 'white';
                    sub.style.border = '1px solid #ccc';
                    sub.style.padding = '6px';
                    sub.style.margin = '0';
                    sub.style.display = 'none';
                    sub.style.minWidth = '180px';
                    sub.style.boxShadow = '0 2px 6px rgba(0,0,0,0.15)';
                    parentLi.style.position = 'relative';
                    parentLi.addEventListener('mouseenter', () => sub.style.display = 'block');
                    parentLi.addEventListener('mouseleave', () => sub.style.display = 'none');
                    parentLi.appendChild(sub);
                    return sub;
                }};
                // Outline
                const liOutline = document.createElement('li');
                liOutline.textContent = 'Outline';
                liOutline.style.padding = '4px 8px';
                liOutline.style.cursor = 'pointer';
                const outlineSub = makeSubmenu(liOutline);
                // Outline color
                const outlineColor = document.createElement('div');
                outlineColor.textContent = 'Color';
                outlineColor.style.padding = '4px 6px';
                outlineColor.style.cursor = 'default';
                const ocRow = document.createElement('div');
                ocRow.style.display = 'flex';
                ocRow.style.alignItems = 'center';
                ocRow.style.gap = '6px';
                ocRow.style.marginTop = '4px';
                const ocInput = document.createElement('input');
                ocInput.type = 'color';
                ocInput.value = tb.border_color && tb.border_color !== 'transparent' ? tb.border_color : '#000000';
                const ocHex = document.createElement('input');
                ocHex.type = 'text';
                ocHex.maxLength = 7;
                ocHex.value = tb.border_color && tb.border_color !== 'transparent' ? tb.border_color : '#000000';
                ocHex.style.width = '80px';
                const ocNoneLabel = document.createElement('label');
                ocNoneLabel.style.display = 'flex';
                ocNoneLabel.style.alignItems = 'center';
                ocNoneLabel.style.gap = '6px';
                ocNoneLabel.style.padding = '4px 6px';
                const ocNoneRadio = document.createElement('input');
                ocNoneRadio.type = 'radio';
                ocNoneRadio.name = `outline_none_${id}`;
                ocNoneRadio.checked = !tb.border_color || tb.border_color === 'transparent';
                const ocNoneText = document.createElement('span');
                ocNoneText.textContent = 'Transparent';
                ocNoneLabel.append(ocNoneRadio, ocNoneText);
                const applyOutlineColor = (val) => {{
                    tb.border_color = val;
                    ocNoneRadio.checked = (val === 'transparent');
                    updateBoxAppearance();
                }};
                const tryOcHex = () => {{
                    let val = ocHex.value.trim();
                    if (!val.startsWith('#')) val = '#' + val;
                    if (/^#([0-9a-fA-F]{{3}}|[0-9a-fA-F]{{6}})$/.test(val)) {{
                        ocInput.value = val;
                        applyOutlineColor(val);
                    }}
                }};
                ocInput.addEventListener('input', () => {{
                    ocHex.value = ocInput.value;
                    applyOutlineColor(ocInput.value);
                    ocNoneRadio.checked = false;
                }});
                ocHex.addEventListener('change', tryOcHex);
                ocHex.addEventListener('keydown', (e) => {{ if (e.key === 'Enter') tryOcHex(); }});
                ocNoneRadio.addEventListener('mousedown', (e) => e.stopPropagation());
                ocNoneRadio.addEventListener('click', (e) => {{
                    e.stopPropagation();
                    applyOutlineColor('transparent');
                    ocHex.value = '#000000';
                    ocInput.value = '#000000';
                }});
                ocRow.append(ocInput, ocHex);
                outlineSub.appendChild(outlineColor);
                outlineSub.appendChild(ocRow);
                outlineSub.appendChild(ocNoneLabel);
                // Thickness
                const thickRow = document.createElement('div');
                thickRow.style.display = 'flex';
                thickRow.style.alignItems = 'center';
                thickRow.style.gap = '6px';
                thickRow.style.marginTop = '8px';
                const thickLabel = document.createElement('span');
                thickLabel.textContent = 'Thickness';
                thickLabel.style.fontSize = '11px';
                const thickInput = document.createElement('input');
                thickInput.type = 'number';
                thickInput.min = '0';
                thickInput.max = '20';
                thickInput.step = '0.5';
                thickInput.value = tb.border_width || 1;
                const applyThickness = () => {{
                    const val = Number(thickInput.value);
                    if (Number.isFinite(val)) {{
                        tb.border_width = Math.max(0, Math.min(20, val));
                        thickInput.value = tb.border_width;
                        updateBoxAppearance();
                    }}
                }};
                thickInput.addEventListener('change', applyThickness);
                thickInput.addEventListener('keydown', (e) => {{ if (e.key === 'Enter') applyThickness(); }});
                thickRow.append(thickLabel, thickInput);
                outlineSub.appendChild(thickRow);
                ul.appendChild(liOutline);
                // Background color
                const liBg = document.createElement('li');
                liBg.textContent = 'Color';
                liBg.style.padding = '4px 8px';
                liBg.style.cursor = 'pointer';
                const bgSub = makeSubmenu(liBg);
                const bgRow = document.createElement('div');
                bgRow.style.display = 'flex';
                bgRow.style.alignItems = 'center';
                bgRow.style.gap = '6px';
                bgRow.style.marginBottom = '6px';
                const bgInput = document.createElement('input');
                bgInput.type = 'color';
                bgInput.value = tb.bgcolor && tb.bgcolor !== 'transparent' ? tb.bgcolor : '#ffffff';
                const bgHex = document.createElement('input');
                bgHex.type = 'text';
                bgHex.maxLength = 7;
                bgHex.style.width = '80px';
                bgHex.value = tb.bgcolor && tb.bgcolor !== 'transparent' ? tb.bgcolor : '#ffffff';
                const applyBg = (val) => {{
                    tb.bgcolor = val;
                    updateBoxAppearance();
                }};
                const tryBgHex = () => {{
                    let val = bgHex.value.trim();
                    if (!val.startsWith('#')) val = '#' + val;
                    if (/^#([0-9a-fA-F]{{3}}|[0-9a-fA-F]{{6}})$/.test(val)) {{
                        bgInput.value = val;
                        applyBg(val);
                    }}
                }};
                bgInput.addEventListener('input', () => {{
                    bgHex.value = bgInput.value;
                    applyBg(bgInput.value);
                }});
                bgHex.addEventListener('change', tryBgHex);
                bgHex.addEventListener('keydown', (e) => {{ if (e.key === 'Enter') tryBgHex(); }});
                bgRow.append(bgInput, bgHex);
                const noneBg = document.createElement('label');
                noneBg.style.display = 'flex';
                noneBg.style.alignItems = 'center';
                noneBg.style.gap = '6px';
                noneBg.style.padding = '4px 6px';
                const noneRadio = document.createElement('input');
                noneRadio.type = 'radio';
                noneRadio.name = `bg_none_${{id}}`;
                noneRadio.checked = !tb.bgcolor || tb.bgcolor === 'transparent';
                noneRadio.addEventListener('mousedown', (e) => e.stopPropagation());
                noneRadio.addEventListener('click', (e) => {{
                    e.stopPropagation();
                    tb.bgcolor = 'transparent';
                    updateBoxAppearance();
                    bgSub.style.display = 'none';
                    menu.remove();
                }});
                const noneText = document.createElement('span');
                noneText.textContent = 'Transparent';
                noneBg.append(noneRadio, noneText);
                bgSub.appendChild(bgRow);
                bgSub.appendChild(noneBg);
                ul.appendChild(liBg);
                // Alignment with horizontal/vertical submenus (skip for shapes)
                if (!isShape) {{
                    const liAlign = document.createElement('li');
                    liAlign.textContent = 'Alignment';
                    liAlign.style.padding = '4px 8px';
                    liAlign.style.cursor = 'pointer';
                    const alignSub = makeSubmenu(liAlign);
                    const liH = document.createElement('div');
                    liH.textContent = 'Horizontal';
                    liH.style.padding = '4px 6px';
                    liH.style.cursor = 'pointer';
                    const hSub = makeSubmenu(liH);
                    ['left','center','right'].forEach(key => {{
                        const item = document.createElement('div');
                        item.textContent = key.charAt(0).toUpperCase() + key.slice(1);
                        item.style.padding = '4px 6px';
                        item.style.cursor = 'pointer';
                        item.addEventListener('mousedown', (e) => e.stopPropagation());
                        item.addEventListener('click', (e) => {{
                            e.stopPropagation();
                            tb.text_style = Object.assign(resolveTextStyle(tb), {{ align: key }});
                            applyStyleToEditor(editor, tb.text_style);
                            updateBoxAppearance();
                            alignSub.style.display = 'none';
                            menu.remove();
                        }});
                        hSub.appendChild(item);
                    }});
                    alignSub.appendChild(liH);
                    const liV = document.createElement('div');
                    liV.textContent = 'Vertical';
                    liV.style.padding = '4px 6px';
                    liV.style.cursor = 'pointer';
                    const vSub = makeSubmenu(liV);
                    [['top','Top'],['center','Center'],['bottom','Bottom']].forEach(([key,label]) => {{
                        const item = document.createElement('div');
                        item.textContent = label;
                        item.style.padding = '4px 6px';
                        item.style.cursor = 'pointer';
                        item.addEventListener('mousedown', (e) => e.stopPropagation());
                        item.addEventListener('click', (e) => {{
                            e.stopPropagation();
                            tb.text_style = Object.assign(resolveTextStyle(tb), {{ vertical: key }});
                            applyStyleToEditor(editor, tb.text_style);
                            updateBoxAppearance();
                            alignSub.style.display = 'none';
                            menu.remove();
                        }});
                        vSub.appendChild(item);
                    }});
                    alignSub.appendChild(liV);
                    ul.appendChild(liAlign);
                }}
                menu.appendChild(ul);
                document.body.appendChild(menu);
                const onDocDown = (e) => {{
                    if (!menu.contains(e.target)) {{
                        try {{ menu.remove(); }} catch (err) {{}}
                        document.removeEventListener('mousedown', onDocDown, true);
                    }}
                }};
                document.addEventListener('mousedown', onDocDown, true);
            }};
            // Keep the background fixed; shift all interactive layers by the inverse offset
            foregroundGroup.translate(fgOffsetX, fgOffsetY);
            groupBackgroundLayer.translate(fgOffsetX, fgOffsetY);
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
                try {{
                    const lbl = group && group.findOne ? group.findOne('[data-type="compound-label"]') : null;
                    if (lbl && lbl.font) lbl.font({{ weight: 'normal' }});
                }} catch (err) {{}}
            }};
            const updateTextBoxSelectionStyle = (group, active) => {{
                setGroupStrokeColor(group, 'text-outline', active ? 'red' : null);
                setGroupLabelColor(group, 'text-label', active ? 'red' : null);
                try {{
                    const lbl = group && group.findOne ? group.findOne('[data-type="text-label"]') : null;
                    if (lbl && lbl.font) lbl.font({{ weight: 'normal' }});
                }} catch (err) {{}}
                try {{
                    const gid = group?.attr ? group.attr('data-id') : null;
                    if (gid) {{
                        const shouldShowHandles = !!active && (!groupsForEntry('text-box', gid) || isEntryInsideActiveGroup({{ type: 'text-box', id: gid }}));
                        ensureTextHandles(gid, findTextBlockByDomId(gid));
                        toggleTextHandles(gid, shouldShowHandles);
                    }}
                }} catch (err) {{}}
            }};
            const selectionKey = (type, id) => `${{type || ''}}::${{id || ''}}`;
            const shouldShowGroupMenu = (type, id, protboxId = null) => {{
                if (selectionMap.size <= 1) return false;
                const key = selectionKey(type, id);
                const protKey = protboxId ? selectionKey('prot-box', protboxId) : null;
                return selectionMap.has(key) || (protKey && selectionMap.has(protKey));
            }};
            const isEntryInsideActiveGroup = (entry) => {{
                if (!entry || !activeGroupEditId) return false;
                if (entry.type === 'group') {{
                    return entry.id === activeGroupEditId || isGroupDescendant(entry.id, activeGroupEditId);
                }}
                const membership = groupsForEntry(entry.type, entry.id);
                if (!membership || !membership.size) return false;
                for (const gid of membership) {{
                    if (gid === activeGroupEditId || isGroupDescendant(gid, activeGroupEditId)) return true;
                }}
                return false;
            }};
            const applySelectionStyles = (entry, asPrimary = false) => {{
                if (!entry || !entry.element) return;
                const {{ element, type, id, protboxId }} = entry;
                let selectionVisual = null;
                if (type === 'compound' || type === 'figure-key') {{
                    updateCompoundSelectionStyle(element, true);
                }} else if (type === 'text-box') {{
                    updateTextBoxSelectionStyle(element, true);
                }} else if (type === 'group') {{
                    const vis = ensureGroupVisual(id);
                    vis.box.stroke({{ color: 'red', width: 1.25, dasharray: '4,2' }}).show();
                    vis.hit.show();
                    selectionVisual = vis.box;
                }} else if (type === 'arrow') {{
                    selectionVisual = draw.findOne(`[data-id="${{id}}"]`);
                    selectionVisual?.stroke({{ color: 'red', width: 1 }});
                }} else if (type === 'arrow-start' || type === 'arrow-end') {{
                    element.fill('red').stroke({{ color: 'red', width: 2 }});
                }} else if (type === 'ptm-label') {{
                    element.stroke({{ color: 'none' }});
                }} else {{
                    element.stroke({{ color: 'red', width: 1 }});
                }}
                if (asPrimary && type === 'prot-box' && id && elementGroups[id]) {{
                    elementGroups[id].forEach(el => (el.element.attr('data-type').startsWith('handle-') || el.element.attr('data-type') === 'ptm-snap-circle') && el.element.show());
                }}
                if (asPrimary && type === 'arrow') {{
                    const arrowId = id;
                    Object.values(elementGroups).forEach(group => group.forEach(el => el.element.attr('data-type').startsWith('handle-') && el.element.hide()));
                    showRelatedHandles(arrowId);
                }} else if (asPrimary && (type === 'arrow-start' || type === 'arrow-end')) {{
                    const arrowId = id.replace(/_(start|end)$/, '');
                    const arrowEl = draw.findOne(`[data-id="${{arrowId}}"]`) || element;
                    const arrowBounds = getScreenBounds(arrowEl);
                    const arrowEnds = [];
                    const endFromHandle = (h) => {{
                        const hb = getScreenBounds(h);
                        return hb ? {{ x: (hb.left + hb.right) / 2, y: (hb.top + hb.bottom) / 2 }} : null;
                    }};
                    const handleGrp = arrowHandleGroups[arrowId];
                    if (handleGrp?.startHandle) {{
                        const p = endFromHandle(handleGrp.startHandle); if (p) arrowEnds.push(p);
                    }}
                    if (handleGrp?.endHandle) {{
                        const p = endFromHandle(handleGrp.endHandle); if (p) arrowEnds.push(p);
                    }}
                    if (!arrowEnds.length && arrowBounds) {{
                        const cy = (arrowBounds.top + arrowBounds.bottom) / 2;
                        arrowEnds.push({{ x: arrowBounds.left, y: cy }});
                        arrowEnds.push({{ x: arrowBounds.right, y: cy }});
                    }}
                    const radiusPx = 50;
                    Object.values(elementGroups).forEach(group => group.forEach(el => {{
                        if (!el.element.attr('data-type').startsWith('handle-')) return;
                        if (!arrowEnds.length) {{ el.element.hide(); return; }}
                        const hb = getScreenBounds(el.element);
                        if (!hb) {{ el.element.hide(); return; }}
                        const hx = (hb.left + hb.right) / 2;
                        const hy = (hb.top + hb.bottom) / 2;
                        const within = arrowEnds.some(p => Math.hypot(hx - p.x, hy - p.y) <= radiusPx);
                        if (within) el.element.show(); else el.element.hide();
                    }}));
                    showRelatedHandles(arrowId);
                }}
                entry.visual = selectionVisual;
            }};
            const clearSelectionStyles = (entry) => {{
                if (!entry || !entry.element) return;
                const {{ element, type, id }} = entry;
                if (type === 'compound' || type === 'figure-key') {{
                    updateCompoundSelectionStyle(element, false);
                }} else if (type === 'text-box') {{
                    updateTextBoxSelectionStyle(element, false);
                }} else if (type === 'group') {{
                    const vis = groupVisuals[id];
                    if (vis) {{
                        vis.box.stroke({{ color: '#666', width: 1, dasharray: '5,5' }});
                    }}
                }} else if (type === 'arrow') {{
                    entry.visual?.stroke({{ color: 'black', width: 1 }});
                }} else if (type === 'arrow-start' || type === 'arrow-end') {{
                    element.stroke({{ color: 'black', width: 1 }});
                }} else if (type === 'ptm-label') {{
                    element.stroke({{ color: 'none' }});
                }} else {{
                    element.stroke({{ color: 'black', width: 1 }});
                }}
                if ((type === 'arrow' || type === 'arrow-start' || type === 'arrow-end') && id) {{
                    const arrowId = type === 'arrow' ? id : id.replace(/_(start|end)$/, '');
                    hideRelatedHandles(arrowId);
                    const arrow = getArrowById(arrowId);
                    if (arrow) {{
                        if (arrow.protbox_id_1 && arrow.protbox_id_1_side) updateArrowPositions(arrow.protbox_id_1, arrow.protbox_id_1_side);
                        if (arrow.protbox_id_2 && arrow.protbox_id_2_side) updateArrowPositions(arrow.protbox_id_2, arrow.protbox_id_2_side);
                    }}
                }}
                if (type === 'prot-box' && id) {{
                    const selectedGroup = elementGroups[id];
                    if (selectedGroup) {{
                        selectedGroup.forEach(el => (el.element.attr('data-type').startsWith('handle-') || el.element.attr('data-type') === 'ptm-snap-circle') && el.element.hide());
                    }}
                }}
            }};
            const deselectElement = () => {{
                selectionMap.forEach(entry => clearSelectionStyles(entry));
                selectionMap.clear();
                primarySelectionKey = null;
                Object.values(elementGroups).forEach(group => group.forEach(el => el.element.attr('data-type').startsWith('handle-') && el.element.hide()));
                Object.keys(textHandleGroups || {{}}).forEach(id => toggleTextHandles(id, false));
                selectedElement = selectedType = selectedId = selectedProtboxId = selectedVisual = null;
                updateGroupSelectionDisplay();
            }};
            const selectElement = (element, type, id, protboxId = null, options = {{}}) => {{
                let targetType = type;
                let targetId = id;
                if (type !== 'group' && !activeGroupEditId) {{
                    const gs = groupsForEntry(type, id);
                    if (gs && gs.size) {{
                        targetId = resolveRootGroup(Array.from(gs)[0]);
                        targetType = 'group';
                    }}
                }}
                const key = selectionKey(targetType, targetId);
                const additive = options.additive === true;
                const toggle = options.toggle === true;
                const alreadySelected = selectionMap.has(key);
                if (toggle && selectionMap.has(key)) {{
                    clearSelectionStyles(selectionMap.get(key));
                    selectionMap.delete(key);
                    if (primarySelectionKey === key) {{
                        primarySelectionKey = selectionMap.size ? Array.from(selectionMap.keys()).pop() : null;
                        const newPrimary = primarySelectionKey ? selectionMap.get(primarySelectionKey) : null;
                        if (newPrimary) {{
                            applySelectionStyles(newPrimary, true);
                            selectedElement = newPrimary.element;
                            selectedType = newPrimary.type;
                            selectedId = newPrimary.id;
                            selectedProtboxId = newPrimary.protboxId || null;
                            selectedVisual = newPrimary.visual || null;
                        }} else {{
                            selectedElement = selectedType = selectedId = selectedProtboxId = selectedVisual = null;
                        }}
                    }}
                    if (!selectionMap.size) {{
                        Object.values(elementGroups).forEach(group => group.forEach(el => el.element.attr('data-type').startsWith('handle-') && el.element.hide()));
                    }}
                    updateGroupSelectionDisplay();
                    return;
                }}
                if (!additive) {{
                    // Preserve existing multi-selection if clicking one of the selected items
                    if (alreadySelected && selectionMap.size > 1) {{
                        primarySelectionKey = key;
                        const current = selectionMap.get(key);
                        if (current) {{
                            applySelectionStyles(current, true);
                            selectedElement = current.element;
                            selectedType = current.type;
                            selectedId = current.id;
                            selectedProtboxId = current.protboxId || null;
                            selectedVisual = current.visual || null;
                            updateGroupSelectionDisplay();
                            return;
                        }}
                    }}
                    deselectElement();
                }}
                const entry = {{ element, type: targetType, id: targetId, protboxId: protboxId || null, visual: null }};
                selectionMap.set(key, entry);
                primarySelectionKey = key;
                selectedElement = element;
                selectedType = targetType;
                selectedId = targetId;
                selectedProtboxId = protboxId || (type === 'prot-box' ? id : (type.startsWith('ptm-') ? protboxId : null));
                selectedVisual = null;
                if (!additive && selectionMap.size > 1) {{
                    applySelectionStyles(entry, true);
                    selectedVisual = entry.visual || selectedVisual;
                    updateGroupSelectionDisplay();
                    return;
                }}
                if (activeGroupEditId && !isEntryInsideActiveGroup(entry)) {{
                    exitGroupEditMode();
                }}
                applySelectionStyles(entry, true);
                selectedVisual = entry.visual || selectedVisual;
                updateGroupSelectionDisplay();
            }};
            const ensureGroupEditIndicator = () => {{
                if (groupEditIndicator) return groupEditIndicator;
                groupEditIndicator = document.createElement('div');
                groupEditIndicator.style.position = 'absolute';
                groupEditIndicator.style.left = '10px';
                groupEditIndicator.style.top = '10px';
                groupEditIndicator.style.padding = '6px 10px';
                groupEditIndicator.style.background = 'rgba(44, 123, 229, 0.85)';
                groupEditIndicator.style.color = '#fff';
                groupEditIndicator.style.fontSize = '12px';
                groupEditIndicator.style.borderRadius = '4px';
                groupEditIndicator.style.boxShadow = '0 2px 6px rgba(0,0,0,0.2)';
                groupEditIndicator.style.zIndex = '1200';
                groupEditIndicator.style.display = 'none';
                container.appendChild(groupEditIndicator);
                return groupEditIndicator;
            }};
            const renderGroupEditIndicator = () => {{
                if (!groupEditIndicator) return;
                if (activeGroupEditId) {{
                    groupEditIndicator.textContent = `Group Edit Mode (${{
                        activeGroupEditId
                    }}) – press Esc to exit`;
                    groupEditIndicator.style.display = 'block';
                }} else {{
                    groupEditIndicator.style.display = 'none';
                }}
            }};
            function normalizeGroupEntry(raw) {{
                if (!raw || typeof raw !== 'object') return null;
                const normalized = Object.assign({{}}, raw);
                normalized.group_id = normalized.group_id || normalized.id || `${{++groupCounter}}`;
                const normalizeMemberId = (mtype, mid) => {{
                    if (mid === undefined || mid === null) return mid;
                    const text = String(mid);
                    if (mtype === 'compound') {{
                        return text.startsWith('compound_') ? text : `compound_${{text}}`;
                    }}
                    if (mtype === 'text-box') {{
                        return text.startsWith('text_') ? text : `text_${{text}}`;
                    }}
                    return mid;
                }};
                const members = Array.isArray(normalized.members) ? normalized.members.filter(Boolean).map(m => {{
                    if (typeof m === 'string' || typeof m === 'number') {{
                        return {{ type: 'prot-box', id: m }};
                    }}
                    const mtype = m.type || (m.group_id ? 'group' : 'prot-box');
                    const mid = m.id || m.protbox_id || m.group_id;
                    const fixedId = normalizeMemberId(mtype, mid);
                    return {{ type: mtype, id: fixedId }};
                }}) : [];
                if (!members.length && Array.isArray(normalized.protbox_ids)) {{
                    normalized.protbox_ids.forEach(pid => members.push({{ type: 'prot-box', id: pid }}));
                }}
                normalized.members = members;
                const numericId = Number(normalized.group_id);
                if (Number.isFinite(numericId)) {{
                    groupCounter = Math.max(groupCounter, numericId);
                }}
                const flat = new Set();
                members.forEach(m => {{
                    if (!m) return;
                    if (m.type === 'prot-box' && m.id !== undefined && m.id !== null) {{
                        const pid = normalizeProtboxId(m.id);
                        if (pid) flat.add(pid);
                    }}
                }});
                normalized.protbox_ids = Array.from(flat);
                return normalized;
            }}
function rebuildGroupIndexes() {{
                groupMembership = {{}};
                groupMembershipAll = {{}};
                groupParents = {{}};
                Object.values(groupMap).forEach(g => {{
                    (g.members || []).forEach(m => {{
                        if (!m) return;
                        const mid = m.id;
                        const mtype = m.type;
                        const keyAll = entryKey(mtype, mid);
                        if (keyAll) {{
                            groupMembershipAll[keyAll] = groupMembershipAll[keyAll] || new Set();
                            groupMembershipAll[keyAll].add(g.group_id);
                        }}
                        if (m.type === 'prot-box') {{
                            const pid = normalizeProtboxId(m.id);
                            if (!pid) return;
                            groupMembership[pid] = groupMembership[pid] || new Set();
                            groupMembership[pid].add(g.group_id);
                        }} else if (m.type === 'group') {{
                            const gid = `${{m.id}}`;
                            groupParents[gid] = groupParents[gid] || new Set();
                            groupParents[gid].add(g.group_id);
                        }}
                    }});
                }});
            }}
            function collectProtboxesForGroup(groupId, visited = new Set()) {{
                const gid = `${{groupId}}`;
                if (visited.has(gid)) return [];
                visited.add(gid);
                const group = groupMap[gid];
                if (!group) return [];
                const acc = new Set();
                (group.members || []).forEach(m => {{
                    if (!m) return;
                    if (m.type === 'prot-box') {{
                        const pid = normalizeProtboxId(m.id);
                        if (pid) acc.add(pid);
                    }} else if (m.type === 'group') {{
                        collectProtboxesForGroup(m.id, visited).forEach(pid => acc.add(pid));
                    }}
                }});
                return Array.from(acc);
            }}
            function collectMembersForGroup(groupId, visited = new Set()) {{
                const gid = `${{groupId}}`;
                if (visited.has(gid)) return [];
                visited.add(gid);
                const group = groupMap[gid];
                if (!group) return [];
                const acc = [];
                (group.members || []).forEach(m => {{
                    if (!m) return;
                    if (m.type === 'group') {{
                        collectMembersForGroup(m.id, visited).forEach(x => acc.push(x));
                    }} else {{
                        acc.push({{ type: m.type, id: m.id }});
                    }}
                }});
                return acc;
            }}
            function translateArrowByDelta(arrowId, deltaX, deltaY) {{
                const line = draw.find(`[data-id="${{arrowId}}"]`)[0];
                const hit = draw.find(`[data-id="${{arrowId}}_hit"]`)[0];
                const head = draw.find(`[data-id="${{arrowId}}_head"]`)[0];
                const startHandle = arrowHandleGroups[arrowId]?.startHandle;
                const endHandle = arrowHandleGroups[arrowId]?.endHandle;
                if (!line) return;
                let x1 = parseFloat(line.attr('x1') || 0) + deltaX;
                let y1 = parseFloat(line.attr('y1') || 0) + deltaY;
                let x2 = parseFloat(line.attr('x2') || 0) + deltaX;
                let y2 = parseFloat(line.attr('y2') || 0) + deltaY;
                line.plot(x1, y1, x2, y2);
                if (hit) hit.plot(x1, y1, x2, y2);
                if (startHandle) startHandle.cx(x1).cy(y1);
                if (endHandle) endHandle.cx(x2).cy(y2);
                const lineType = line.attr('data-line-type') || 'arrow';
                const arrowData = getArrowById(arrowId);
                if (head) {{
                    const dx = x2 - x1, dy = y2 - y1, angle = Math.atan2(dy, dx), arrowSize = 5;
                    if (lineType === 'arrow') {{
                        head.plot([[x2, y2],[x2 - arrowSize * Math.cos(angle + Math.PI / 6), y2 - arrowSize * Math.sin(angle + Math.PI / 6)],[x2 - arrowSize * Math.cos(angle - Math.PI / 6), y2 - arrowSize * Math.sin(angle - Math.PI / 6)]]);
                    }} else if (lineType === 'inhibition') {{
                        const endSide = arrowData?.protbox_id_2_side;
                        let barX1, barY1, barX2, barY2;
                        if (endSide === 'North' || endSide === 'South') {{barX1 = x2 - arrowSize; barY1 = y2; barX2 = x2 + arrowSize; barY2 = y2;}}
                        else if (endSide === 'West' || endSide === 'East') {{barX1 = x2; barY1 = y2 - arrowSize; barX2 = x2; barY2 = y2 + arrowSize;}}
                        else {{const perp = angle + Math.PI / 2; barX1 = x2 + Math.cos(perp) * arrowSize; barY1 = y2 + Math.sin(perp) * arrowSize; barX2 = x2 - Math.cos(perp) * arrowSize; barY2 = y2 - Math.sin(perp) * arrowSize;}}
                        head.plot(barX1, barY1, barX2, barY2);
                    }}
                }}
            }}
            function isProtboxInGroup(protboxId, groupId, visited = new Set()) {{
                const gid = `${{groupId}}`;
                if (visited.has(gid)) return false;
                visited.add(gid);
                const group = groupMap[gid];
                if (!group) return false;
                const pidNorm = normalizeProtboxId(protboxId);
                return (group.members || []).some(m => {{
                    if (!m) return false;
                    if (m.type === 'prot-box') {{
                        return normalizeProtboxId(m.id) === pidNorm;
                    }}
                    if (m.type === 'group') {{
                        return isProtboxInGroup(pidNorm, m.id, visited);
                    }}
                    return false;
                }});
            }}
            function computeBoundsForMembers(members) {{
                if (!members || !members.length) return null;
                let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
                members.forEach(m => {{
                    if (!m) return;
                    if (m.type === 'prot-box') {{
                        const pid = normalizeProtboxId(m.id);
                        const pb = pid ? protboxMap[pid] : null;
                        if (!pb) return;
                        minX = Math.min(minX, pb.x || 0);
                        minY = Math.min(minY, pb.y || 0);
                        maxX = Math.max(maxX, (pb.x || 0) + (pb.width || 0));
                        maxY = Math.max(maxY, (pb.y || 0) + (pb.height || 0));
                    }} else if (m.type === 'text-box' || m.type === 'compound' || m.type === 'figure-key') {{
                        const el = draw.findOne(`[data-id="${{m.id}}"]`);
                        if (!el) return;
                        const bbox = el.bbox();
                        minX = Math.min(minX, bbox.x);
                        minY = Math.min(minY, bbox.y);
                        maxX = Math.max(maxX, bbox.x2);
                        maxY = Math.max(maxY, bbox.y2);
                    }} else if (m.type === 'arrow') {{
                        const line = draw.findOne(`[data-id="${{m.id}}"]`);
                        if (!line) return;
                        const x1 = parseFloat(line.attr('x1') || 0), y1 = parseFloat(line.attr('y1') || 0);
                        const x2 = parseFloat(line.attr('x2') || 0), y2 = parseFloat(line.attr('y2') || 0);
                        minX = Math.min(minX, x1, x2);
                        minY = Math.min(minY, y1, y2);
                        maxX = Math.max(maxX, x1, x2);
                        maxY = Math.max(maxY, y1, y2);
                    }}
                }});
                if (!Number.isFinite(minX) || !Number.isFinite(minY) || !Number.isFinite(maxX) || !Number.isFinite(maxY)) return null;
                return {{
                    x: minX,
                    y: minY,
                    width: (maxX - minX),
                    height: (maxY - minY)
                }};
            }}
            function computeGroupBounds(groupId, padding = 8) {{
                const members = collectMembersForGroup(groupId);
                const base = computeBoundsForMembers(members);
                if (!base) return null;
                return {{
                    x: base.x - padding,
                    y: base.y - padding,
                    width: base.width + padding * 2,
                    height: base.height + padding * 2
                }};
            }}
            function computeGroupBackgroundBounds(groupId, padding = 10, anchorId = null) {{
                const membersAll = collectMembersForGroup(groupId);
                if (!membersAll.length) return null;
                let members = membersAll;
                if (anchorId) {{
                    const anchorNorm = normalizeProtboxId(anchorId);
                    members = membersAll.filter(m => !(m.type === 'prot-box' && normalizeProtboxId(m.id) === anchorNorm));
                }}
                let base = computeBoundsForMembers(members);
                if (!base) base = computeBoundsForMembers(membersAll);
                if (!base) return null;
                let left = base.x - padding;
                let right = base.x + base.width + padding;
                let top = base.y - padding;
                let bottom = base.y + base.height + padding;
                if (anchorId) {{
                    const anchorNorm = normalizeProtboxId(anchorId);
                    const pb = anchorNorm ? protboxMap[anchorNorm] : null;
                    if (pb) {{
                        const cx = (pb.x || 0) + (pb.width || 0) / 2;
                        const cy = (pb.y || 0) + (pb.height || 0) / 2;
                        const distLeft = Math.abs(cx - left);
                        const distRight = Math.abs(right - cx);
                        const distTop = Math.abs(cy - top);
                        const distBottom = Math.abs(bottom - cy);
                        const minDist = Math.min(distLeft, distRight, distTop, distBottom);
                        const inside = cx >= left && cx <= right && cy >= top && cy <= bottom;
                        if (!inside || minDist <= padding) {{
                            if (minDist === distLeft) left = cx;
                            else if (minDist === distRight) right = cx;
                            else if (minDist === distTop) top = cy;
                            else bottom = cy;
                        }}
                    }}
                }}
                return {{
                    x: left,
                    y: top,
                    width: right - left,
                    height: bottom - top
                }};
            }}
            function ensureGroupBackground(groupId) {{
                if (groupBackgrounds[groupId]) return groupBackgrounds[groupId];
                const box = groupBackgroundLayer.rect(10, 10)
                    .fill('#e0e0e0')
                    .stroke({{ color: '#888', width: 1 }})
                    .attr({{ 'data-type': 'group-background', 'data-group-id': groupId }});
                box.node.style.pointerEvents = 'none';
                groupBackgrounds[groupId] = {{ box }};
                return groupBackgrounds[groupId];
            }}
            function ensureGroupVisual(groupId) {{
                if (groupVisuals[groupId]) return groupVisuals[groupId];
                const box = groupOverlay.rect(10, 10).fill('none').stroke({{ color: '#666', width: 1, dasharray: '5,5' }}).attr({{ 'data-type': 'group-box', 'data-group-id': groupId }}).hide();
                const hit = groupOverlay.rect(10, 10).fill({{ color: '#fff', opacity: 0.001 }}).stroke({{ color: 'transparent', width: 0 }}).attr({{ 'data-type': 'group-hit', 'data-group-id': groupId }}).hide();
                const label = groupOverlay.text('Group Edit').font({{ size: 10, family: 'Arial', anchor: 'start' }}).fill('#2c7be5').attr({{ 'data-type': 'group-label', 'data-group-id': groupId }}).hide();
                hit.node.addEventListener('dblclick', evt => {{
                    evt.preventDefault();
                    evt.stopPropagation();
                    selectElement(hit, 'group', groupId, null);
                    enterGroupEditMode(groupId);
                }});
                hit.node.addEventListener('mousedown', evt => evt.stopPropagation());
                hit.node.addEventListener('click', evt => evt.stopPropagation());
                makeDraggable(hit, 'group', groupId);
                groupVisuals[groupId] = {{ box, hit, label }};
                return groupVisuals[groupId];
            }}
            function refreshGroupVisual(groupId) {{
                const group = groupMap[groupId] || {{}};
                const bounds = computeGroupBounds(groupId);
                const vis = ensureGroupVisual(groupId);
                const showBackground = !!group.show_box;
                if (showBackground) {{
                    const padRaw = Number(group.box_padding);
                    const pad = Number.isFinite(padRaw) ? padRaw : 10;
                    const radiusRaw = Number(group.box_radius);
                    const radius = Number.isFinite(radiusRaw) ? radiusRaw : 8;
                    const bgBounds = computeGroupBackgroundBounds(groupId, pad, group.anchor_member);
                    const bg = ensureGroupBackground(groupId);
                    if (!bgBounds) {{
                        bg.box.hide();
                    }} else {{
                        bg.box.width(bgBounds.width).height(bgBounds.height).move(bgBounds.x, bgBounds.y);
                        bg.box.radius(radius);
                        bg.box.show();
                    }}
                }} else if (groupBackgrounds[groupId]) {{
                    groupBackgrounds[groupId].box.hide();
                }}
                if (!bounds) {{
                    vis.box.hide();
                    vis.hit.hide();
                    vis.label.hide();
                    return;
                }}
                vis.box.width(bounds.width).height(bounds.height).move(bounds.x, bounds.y);
                vis.hit.width(bounds.width).height(bounds.height).move(bounds.x, bounds.y);
                vis.label.move(bounds.x + 4, bounds.y - 10);
            }}
            document.addEventListener('dblclick', (e) => {{
                if (!activeGroupEditId) return;
                if (!isBackgroundTarget(e.target)) return;
                exitGroupEditMode();
            }});
            function groupHasSelectedMember(groupId) {{
                let hasSelected = false;
                selectionMap.forEach(entry => {{
                    if (hasSelected || !entry) return;
                    if (entry.type === 'group') {{
                        if (`${{entry.id}}` === `${{groupId}}`) {{
                            hasSelected = true;
                        }}
                        return;
                    }}
                    const pbId = entry.type === 'prot-box' ? entry.id : entry.protboxId;
                    if (pbId && isProtboxInGroup(pbId, groupId)) {{
                        hasSelected = true;
                    }}
                }});
                return hasSelected;
            }}
            function updateGroupSelectionDisplay() {{
                Object.keys(groupMap).forEach(gid => {{
                    const vis = ensureGroupVisual(gid);
                    const isSelected = selectionMap.has(selectionKey('group', gid));
                    const isEdit = activeGroupEditId === gid;
                    const show = isSelected || groupHasSelectedMember(gid) || isEdit;
                    if (!show) {{
                        vis.box.hide();
                        vis.hit.show();
                        vis.hit.attr('pointer-events', 'all');
                        vis.label.hide();
                        return;
                    }}
                    vis.box.show();
                    vis.hit.show();
                    if (isEdit) {{
                        vis.box.stroke({{ color: '#2c7be5', width: 1.5, dasharray: '4,2' }});
                        vis.hit.attr('pointer-events', 'none');
                        vis.label.text('Group Edit Mode').show();
                    }} else {{
                        vis.box.stroke({{ color: isSelected ? 'red' : '#666', width: 1, dasharray: '5,5' }});
                        vis.hit.attr('pointer-events', 'all');
                        vis.label.hide();
                    }}
                }});
                renderGroupEditIndicator();
            }}
            function refreshGroupsForProtbox(protboxId) {{
                const pid = normalizeProtboxId(protboxId);
                const queued = pid && groupMembership[pid] ? Array.from(groupMembership[pid]) : [];
                const visited = new Set();
                while (queued.length) {{
                    const gid = queued.pop();
                    if (!gid || visited.has(gid)) continue;
                    visited.add(gid);
                    refreshGroupVisual(gid);
                    if (groupParents[gid]) {{
                        groupParents[gid].forEach(parent => queued.push(parent));
                    }}
                }}
                updateGroupSelectionDisplay();
            }}
            function syncGroupsToServer() {{
                const payload = Object.values(groupMap).map(g => {{
                    const protIds = collectProtboxesForGroup(g.group_id);
                    return {{
                        group_id: g.group_id,
                        members: g.members || [],
                        protbox_ids: protIds,
                        show_box: g.show_box || false,
                        box_padding: g.box_padding,
                        box_radius: g.box_radius,
                        anchor_member: g.anchor_member
                    }};
                }});
                Shiny?.setInputValue('groups_changed', {{ groups: payload, timestamp: Date.now() }}, {{ priority: 'event' }});
            }}
            function initializeGroupState() {{
                Object.keys(groupVisuals).forEach(k => {{
                    const vis = groupVisuals[k];
                    try {{ vis.box.remove(); vis.hit.remove(); vis.label?.remove(); }} catch (err) {{}}
                    delete groupVisuals[k];
                }});
                Object.keys(groupBackgrounds).forEach(k => {{
                    const vis = groupBackgrounds[k];
                    try {{ vis.box.remove(); }} catch (err) {{}}
                    delete groupBackgrounds[k];
                }});
                Object.keys(groupMap).forEach(k => delete groupMap[k]);
                groupCounter = 0;
                if (Array.isArray(groups)) {{
                    groups.forEach(raw => {{
                        const norm = normalizeGroupEntry(raw);
                        if (norm) {{
                            groupMap[norm.group_id] = norm;
                        }}
                    }});
                }}
                Object.keys(groupMap).forEach(gid => {{
                    groupMap[gid].protbox_ids = collectProtboxesForGroup(gid);
                    const hasMembers = Array.isArray(groupMap[gid].members) && groupMap[gid].members.length;
                    if (!groupMap[gid].protbox_ids.length && !hasMembers) {{
                        delete groupMap[gid];
                    }}
                }});
                if (Array.isArray(groups)) {{
                    groups.length = 0;
                    Object.values(groupMap).forEach(g => groups.push(g));
                }}
                rebuildGroupIndexes();
                Object.keys(groupMap).forEach(gid => refreshGroupVisual(gid));
                updateGroupSelectionDisplay();
            }}
            function createGroupFromSelection() {{
                const members = [];
                const seen = new Set();
                selectionMap.forEach(entry => {{
                    if (!entry) return;
                    if (entry.type === 'prot-box') {{
                        const pid = normalizeProtboxId(entry.id);
                        if (!pid || seen.has('pb:' + pid)) return;
                        seen.add('pb:' + pid);
                        members.push({{ type: 'prot-box', id: pid }});
                    }} else if (entry.type === 'group') {{
                        const gid = `${{entry.id}}`;
                        if (seen.has('g:' + gid)) return;
                        seen.add('g:' + gid);
                        members.push({{ type: 'group', id: gid }});
                    }} else if (entry.type === 'arrow' || entry.type === 'text-box' || entry.type === 'compound' || entry.type === 'figure-key') {{
                        const key = entryKey(entry.type, entry.id);
                        if (seen.has(key)) return;
                        seen.add(key);
                        members.push({{ type: entry.type, id: entry.id }});
                    }}
                }});
                if (members.length < 2) return;
                const newId = `${{++groupCounter}}`;
                const protIds = new Set();
                members.forEach(m => {{
                    if (m.type === 'prot-box') {{
                        const pid = normalizeProtboxId(m.id);
                        if (pid) protIds.add(pid);
                    }} else if (m.type === 'group') {{
                        collectProtboxesForGroup(m.id).forEach(pid => protIds.add(pid));
                    }}
                }});
                const newGroup = {{ group_id: newId, members, protbox_ids: Array.from(protIds) }};
                groups?.push(newGroup);
                groupMap[newId] = newGroup;
                rebuildGroupIndexes();
                refreshGroupVisual(newId);
                updateGroupSelectionDisplay();
                syncGroupsToServer();
                deselectElement();
                const vis = ensureGroupVisual(newId);
                selectElement(vis.hit || vis.box, 'group', newId, null);
            }}
            function ungroupGroup(groupId) {{
                const gid = `${{groupId}}`;
                if (!groupMap[gid]) return;
                delete groupMap[gid];
                if (Array.isArray(groups)) {{
                    for (let i = groups.length - 1; i >= 0; i--) {{
                        if (`${{groups[i]?.group_id}}` === gid) {{
                            groups.splice(i, 1);
                        }} else if (Array.isArray(groups[i]?.members)) {{
                            groups[i].members = groups[i].members.filter(m => !(m && m.type === 'group' && `${{m.id}}` === gid));
                        }}
                    }}
                }}
                const vis = groupVisuals[gid];
                if (vis) {{
                    try {{ vis.box.remove(); vis.hit.remove(); vis.label?.remove(); }} catch (err) {{}}
                    delete groupVisuals[gid];
                }}
                const bg = groupBackgrounds[gid];
                if (bg) {{
                    try {{ bg.box.remove(); }} catch (err) {{}}
                    delete groupBackgrounds[gid];
                }}
                if (activeGroupEditId === gid) {{
                    activeGroupEditId = null;
                }}
                rebuildGroupIndexes();
                Object.keys(groupMap).forEach(id => refreshGroupVisual(id));
                updateGroupSelectionDisplay();
                syncGroupsToServer();
                deselectElement();
            }}
            function adjustGroupStacking(groupId, direction = 'front') {{
                const ids = collectProtboxesForGroup(groupId);
                ids.forEach(pid => {{
                    const rect = draw.findOne(`[data-id="${{pid}}"]`);
                    const entries = elementGroups[pid] || [];
                    entries.forEach(ent => {{
                        if (typeof ent.element?.[direction] === 'function') {{
                            ent.element[direction]();
                        }}
                    }});
                    if (rect && typeof rect[direction] === 'function') rect[direction]();
                }});
                const vis = groupVisuals[groupId];
                if (vis) {{
                    ['box', 'hit', 'label'].forEach(key => {{
                        const el = vis[key];
                        if (el && typeof el[direction] === 'function') el[direction]();
                    }});
                }}
            }}
            function moveGroupMembers(groupId, deltaX, deltaY, options = {{}}) {{
                const members = collectMembersForGroup(groupId);
                const visitedHandles = new Set();
                members.forEach(m => {{
                    if (!m) return;
                    if (m.type === 'prot-box') {{
                        const pid = normalizeProtboxId(m.id);
                        const rect = draw.findOne(`[data-id="${{pid}}"]`);
                        const delta = translateProtbox(pid, deltaX, deltaY, {{ rect, skipRectMove: false }}, visitedHandles) || {{ dx: deltaX, dy: deltaY }};
                        propagateLinkedMove(pid, delta.dx, delta.dy, new Set([pid]), visitedHandles);
                        if (options.snap !== false && !(shiftKeyDown && activeDragProtboxId === pid)) {{
                            maybeSnapProtbox(pid, rect);
                        }}
                    }} else if (m.type === 'text-box') {{
                        const tb = findTextBlockByDomId(m.id);
                        if (tb) {{
                            tb.x = toCoordinateNumber((tb.x || 0) + deltaX, 0);
                            tb.y = toCoordinateNumber((tb.y || 0) + deltaY, 0);
                            applyTextBoxLayout(m.id, tb);
                            Shiny?.setInputValue('element_moved', {{ type: 'text-box', id: m.id, x: tb.x, y: tb.y, width: tb.width, height: tb.height }}, {{ priority: 'event' }});
                        }} else {{
                            const grp = draw.findOne(`[data-id="${{m.id}}"]`);
                            if (grp) grp.dmove(deltaX, deltaY);
                        }}
                        positionTextHandles(m.id);
                    }} else if (m.type === 'compound' || m.type === 'figure-key') {{
                        const grp = draw.findOne(`[data-id="${{m.id}}"]`);
                        if (grp) grp.dmove(deltaX, deltaY);
                    }} else if (m.type === 'arrow') {{
                        translateArrowByDelta(m.id, deltaX, deltaY);
                    }}
                }});
                collectProtboxesForGroup(groupId).forEach(pid => refreshGroupsForProtbox(pid));
            }}
            function resolveGroupingMeta(target) {{
                let node = target;
                while (node && node !== draw.node) {{
                    if (node.getAttribute) {{
                        const dtype = node.getAttribute('data-type');
                        const gid = node.getAttribute('data-group-id');
                        const did = node.getAttribute('data-id');
                        const arrowOwner = node.getAttribute('data-arrow-id');
                        if (dtype === 'group-box' || dtype === 'group-hit') {{
                            return {{ type: 'group', id: gid }};
                        }}
                        if (dtype === 'prot-box') {{
                            return {{ type: 'prot-box', id: did }};
                        }}
                        if (dtype === 'text-box') {{
                            return {{ type: 'text-box', id: did }};
                        }}
                        if (dtype === 'compound') {{
                            return {{ type: 'compound', id: did }};
                        }}
                        if (dtype === 'figure-key') {{
                            return {{ type: 'figure-key', id: did }};
                        }}
                        if (dtype === 'arrow' || dtype === 'arrow-hitbox' || dtype === 'arrow-head') {{
                            const aid = arrowOwner || did;
                            if (aid) {{
                                return {{ type: 'arrow', id: aid }};
                            }}
                        }}
                        if (dtype === 'arrow-start' || dtype === 'arrow-end') {{
                            const aid = (arrowOwner || did || '').replace(/_(start|end)$/, '');
                            if (aid) {{
                                return {{ type: 'arrow', id: aid }};
                            }}
                        }}
                        if (dtype === 'ptm-shape' || dtype === 'ptm-label' || dtype === 'ptm-symbol') {{
                            return {{ type: dtype, id: did, protboxId: node.getAttribute('data-protbox-id') }};
                        }}
                    }}
                    node = node.parentNode;
                }}
                return null;
            }}
            function showGroupingMenu(evt, options) {{
                const existing = document.querySelector('.group-context-menu');
                if (existing) existing.remove();
                const menu = document.createElement('div');
                menu.className = 'group-context-menu';
                menu.style.position = 'absolute';
                const pos = toContainerPosition(evt);
                menu.style.left = `${{pos.x + menuOffsetX}}px`;
                menu.style.top = `${{pos.y + menuOffsetY}}px`;
                menu.style.background = '#fff';
                menu.style.border = '1px solid #ccc';
                menu.style.boxShadow = '0 2px 8px rgba(0,0,0,0.15)';
                menu.style.padding = '2px 0';
                menu.style.minWidth = '120px';
                menu.style.zIndex = '1200';
                const addItem = (label, handler) => {{
                    const item = document.createElement('div');
                    item.textContent = label;
                    item.style.padding = '6px 10px';
                    item.style.cursor = 'pointer';
                    item.addEventListener('mouseenter', () => item.style.backgroundColor = '#eef');
                    item.addEventListener('mouseleave', () => item.style.backgroundColor = 'transparent');
                    item.addEventListener('click', () => {{
                        try {{ handler(); }} catch (err) {{}}
                        menu.remove();
                    }});
                    menu.appendChild(item);
                }};
                if (options.mode === 'create') {{
                    const selProtIds = getSelectedProtboxes();
                    if (selProtIds.length > 1) {{
                        addItem('Add Interactions', () => {{
                            try {{
                                const made = autoConnectSelectedProtboxes();
                                if (!made) console.log && console.log('mk: add interactions produced no arrows');
                            }} catch (err) {{
                                console.warn && console.warn('mk: add interactions failed', err);
                            }}
                        }});
                    }}
                    const selectedProtIds = [];
                    selectionMap.forEach(entry => {{
                        if (entry?.type === 'prot-box') selectedProtIds.push(entry.id);
                    }});
                    const canAddInteraction = selectedProtIds.length === 2;
                    addItem('Group', () => createGroupFromSelection());
                    if (canAddInteraction) {{
                        const liInteraction = document.createElement('div');
                        liInteraction.textContent = 'Add Interaction';
                        liInteraction.style.padding = '6px 10px';
                        liInteraction.style.cursor = 'pointer';
                        liInteraction.style.position = 'relative';
                        const sub = document.createElement('div');
                        sub.style.display = 'none';
                        sub.style.flexDirection = 'column';
                        sub.style.position = 'absolute';
                        sub.style.left = '100%';
                        sub.style.top = '0';
                        sub.style.background = '#fff';
                        sub.style.border = '1px solid #ccc';
                        sub.style.boxShadow = '0 2px 8px rgba(0,0,0,0.15)';
                        ['arrow','inhibition','line'].forEach(t => {{
                            const opt = document.createElement('div');
                            opt.textContent = t === 'arrow' ? 'Arrow' : (t === 'line' ? 'Line' : 'Inhibitor');
                            opt.style.padding = '6px 10px';
                            opt.style.cursor = 'pointer';
                            opt.addEventListener('mouseenter', () => opt.style.backgroundColor = '#eef');
                            opt.addEventListener('mouseleave', () => opt.style.backgroundColor = 'transparent');
                            opt.addEventListener('click', () => {{
                                const targetId = options?.targetProtboxId || selectedProtIds[0];
                                const otherId = selectedProtIds.find(pid => pid !== targetId) || selectedProtIds[0];
                                createInteractionBetweenProtboxes(targetId, otherId, t);
                                menu.remove();
                            }});
                            sub.appendChild(opt);
                        }});
                        liInteraction.addEventListener('mouseenter', () => {{
                            liInteraction.style.backgroundColor = '#eef';
                            sub.style.display = 'flex';
                        }});
                        liInteraction.addEventListener('mouseleave', () => {{
                            liInteraction.style.backgroundColor = 'transparent';
                            sub.style.display = 'none';
                        }});
                        liInteraction.appendChild(sub);
                        menu.appendChild(liInteraction);
                    }}
                }} else if (options.mode === 'group') {{
                    addItem('Ungroup', () => ungroupGroup(options.groupId));
                    addItem('Bring to Front', () => adjustGroupStacking(options.groupId, 'front'));
                    addItem('Send to Back', () => adjustGroupStacking(options.groupId, 'back'));
                }}
                container.appendChild(menu);
                const cleanup = () => menu.remove();
                setTimeout(() => document.addEventListener('click', cleanup, {{ once: true }}), 0);
            }}
            const showArrowMenu = (evt, arrowId) => {{
                const existing = document.querySelector('.context-menu');
                if (existing) existing.remove();
                evt.preventDefault();
                evt.stopPropagation();
                const arrow = getArrowById(arrowId);
                if (!arrow) return;
                const menu = document.createElement('div');
                menu.className = 'context-menu';
                menu.style.position = 'absolute';
                const pos = toContainerPosition(evt);
                menu.style.left = `${{pos.x + menuOffsetX}}px`;
                menu.style.top = `${{pos.y + menuOffsetY}}px`;
                menu.style.background = '#fff';
                menu.style.border = '1px solid #ccc';
                menu.style.boxShadow = '0 2px 8px rgba(0,0,0,0.15)';
                menu.style.padding = '2px 0';
                menu.style.minWidth = '140px';
                menu.style.zIndex = '1200';
                const addItem = (label, handler, disabled = false) => {{
                    const item = document.createElement('div');
                    item.textContent = label;
                    item.style.padding = '4px 8px';
                    item.style.fontSize = '12px';
                    item.style.cursor = disabled ? 'default' : 'pointer';
                    item.style.color = disabled ? '#999' : '#000';
                    if (!disabled) {{
                        item.addEventListener('mouseenter', () => item.style.backgroundColor = '#eef');
                        item.addEventListener('mouseleave', () => item.style.backgroundColor = 'transparent');
                        item.addEventListener('click', () => {{
                            try {{ handler(); }} catch (err) {{}}
                            menu.remove();
                        }});
                    }}
                    menu.appendChild(item);
                }};
                const flipArrow = () => {{
                    const line = draw.findOne(`[data-id="${{arrowId}}"]`);
                    const hit = draw.findOne(`[data-id="${{arrowId}}_hit"]`);
                    if (line) {{
                        const x1 = line.attr('x1'), y1 = line.attr('y1'), x2 = line.attr('x2'), y2 = line.attr('y2');
                        line.attr({{ x1: x2, y1: y2, x2: x1, y2: y1 }});
                        if (hit) hit.attr({{ x1: x2, y1: y2, x2: x1, y2: y1 }});
                    }}
                    const swapFields = (obj, a, b) => {{ const tmp = obj[a]; obj[a] = obj[b]; obj[b] = tmp; }};
                    swapFields(arrow, 'x1', 'x2');
                    swapFields(arrow, 'y1', 'y2');
                    swapFields(arrow, 'protbox_id_1', 'protbox_id_2');
                    swapFields(arrow, 'protbox_id_1_side', 'protbox_id_2_side');
                    swapFields(arrow, 'attached_arrow_1', 'attached_arrow_2');
                    swapFields(arrow, 'attached_end_1', 'attached_end_2');
                    updateArrowVisual(arrowId);
                    rebuildAttachmentsIndex();
                    updateAttachedArrows(arrowId, 'start');
                    updateAttachedArrows(arrowId, 'end');
                }};
                const changeArrowType = (newType) => {{
                    arrow.line = newType;
                    updateArrowVisual(arrowId);
                    translateArrowByDelta(arrowId, 0, 0);
                    rebuildAttachmentsIndex();
                    updateAttachedArrows(arrowId, 'start');
                    updateAttachedArrows(arrowId, 'end');
                }};
                const toggleDash = () => {{
                    arrow.dashed = !arrow.dashed;
                    updateArrowVisual(arrowId);
                }};
                const disconnectArrow = () => {{
                    const line = draw.findOne(`[data-id="${{arrowId}}"]`);
                    const hit = draw.findOne(`[data-id="${{arrowId}}_hit"]`);
                    const x1 = line ? parseFloat(line.attr('x1') || 0) : toFiniteNumber(arrow.x1) || 0;
                    const y1 = line ? parseFloat(line.attr('y1') || 0) : toFiniteNumber(arrow.y1) || 0;
                    const x2 = line ? parseFloat(line.attr('x2') || 0) : toFiniteNumber(arrow.x2) || 0;
                    const y2 = line ? parseFloat(line.attr('y2') || 0) : toFiniteNumber(arrow.y2) || 0;
                    arrow.x1 = x1; arrow.y1 = y1; arrow.x2 = x2; arrow.y2 = y2;
                    delete arrow.protbox_id_1; delete arrow.protbox_id_1_side;
                    delete arrow.protbox_id_2; delete arrow.protbox_id_2_side;
                    cleanupProtboxAttachmentsForArrow(arrowId);
                    if (line) line.plot(x1, y1, x2, y2);
                    if (hit) hit.plot(x1, y1, x2, y2);
                    ensureArrowHandle(arrowId, 'start');
                    ensureArrowHandle(arrowId, 'end');
                    rebuildAttachmentsIndex();
                    updateAttachedArrows(arrowId, 'start');
                    updateAttachedArrows(arrowId, 'end');
                }};
                addItem('Flip Interaction', flipArrow);
                const liType = document.createElement('div');
                liType.textContent = 'Change Type';
                liType.style.padding = '6px 10px';
                liType.style.cursor = 'pointer';
                liType.style.position = 'relative';
                liType.addEventListener('mouseenter', () => {{
                    liType.style.backgroundColor = '#eef';
                    typeSub.style.display = 'flex';
                }});
                liType.addEventListener('mouseleave', () => {{
                    liType.style.backgroundColor = 'transparent';
                    typeSub.style.display = 'none';
                }});
                const typeSub = document.createElement('div');
                typeSub.style.display = 'none';
                typeSub.style.flexDirection = 'column';
                typeSub.style.position = 'absolute';
                typeSub.style.left = '100%';
                typeSub.style.top = '0';
                typeSub.style.background = '#fff';
                typeSub.style.border = '1px solid #ccc';
                typeSub.style.boxShadow = '0 2px 8px rgba(0,0,0,0.15)';
                typeSub.style.minWidth = '140px';
                typeSub.style.zIndex = '1201';
                ['arrow','line','inhibition'].forEach(t => {{
                    const opt = document.createElement('div');
                    opt.textContent = t === 'arrow' ? 'Arrow' : (t === 'line' ? 'Line' : 'Inhibitor');
                    opt.style.padding = '6px 10px';
                    opt.style.cursor = t === arrow.line ? 'default' : 'pointer';
                    opt.style.color = t === arrow.line ? '#999' : '#000';
                    if (t !== arrow.line) {{
                        opt.addEventListener('mouseenter', () => opt.style.backgroundColor = '#eef');
                        opt.addEventListener('mouseleave', () => opt.style.backgroundColor = 'transparent');
                        opt.addEventListener('click', () => {{ changeArrowType(t); menu.remove(); }});
                    }}
                    typeSub.appendChild(opt);
                }});
                liType.appendChild(typeSub);
                menu.appendChild(liType);
                addItem('Disconnect', disconnectArrow);
                container.appendChild(menu);
                setTimeout(() => document.addEventListener('click', () => menu.remove(), {{ once: true }}), 0);
            }};
            function enterGroupEditMode(groupId) {{
                activeGroupEditId = `${{groupId}}`;
                ensureGroupEditIndicator();
                renderGroupEditIndicator();
                updateGroupSelectionDisplay();
            }}
            function exitGroupEditMode() {{
                if (!activeGroupEditId) return;
                activeGroupEditId = null;
                renderGroupEditIndicator();
                updateGroupSelectionDisplay();
            }}
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
                if (!isBackgroundTarget(e.target)) return;
                if (e.ctrlKey || e.metaKey) {{
                    e.preventDefault();
                    const pt = draw.point(e.clientX, e.clientY);
                    selectionBoxActive = true;
                    selectionBoxStart = pt;
                    selectionBoxEnd = pt;
                    selectionBoxStartScreen = {{ x: e.clientX, y: e.clientY }};
                    selectionBoxEndScreen = {{ x: e.clientX, y: e.clientY }};
                    if (selectionBoxRect) selectionBoxRect.remove();
                    selectionBoxRect = handleGroup.rect(0, 0).move(pt.x, pt.y).fill({{ color: '#4da3ff', opacity: 0.12 }}).stroke({{ color: '#4da3ff', width: 1, dasharray: '4,2' }});
                    document.getElementById('svgCanvas').style.cursor = 'crosshair';
                    return;
                }}
                e.preventDefault();
                isPanning = true;
                startPanX = e.clientX;
                startPanY = e.clientY;
                startViewBoxX = viewBox.x;
                startViewBoxY = viewBox.y;
                document.getElementById('svgCanvas').style.cursor = 'grab';
            }};
            const pan = e => {{
                if (selectionBoxActive && selectionBoxRect && selectionBoxStart) {{
                    e.preventDefault();
                    const pt = draw.point(e.clientX, e.clientY);
                    selectionBoxEnd = pt;
                    selectionBoxEndScreen = {{ x: e.clientX, y: e.clientY }};
                    const x = Math.min(selectionBoxStart.x, pt.x);
                    const y = Math.min(selectionBoxStart.y, pt.y);
                    const w = Math.abs(pt.x - selectionBoxStart.x);
                    const h = Math.abs(pt.y - selectionBoxStart.y);
                    selectionBoxRect.size(w, h).move(x, y);
                    return;
                }}
                if (!isPanning) return;
                e.preventDefault();
                viewBox.x = startViewBoxX - (e.clientX - startPanX) / zoomLevel;
                viewBox.y = startViewBoxY - (e.clientY - startPanY) / zoomLevel;
                draw.viewbox(viewBox.x, viewBox.y, viewBox.width, viewBox.height);
            }};
            const endPanning = (e) => {{
                if (selectionBoxActive) {{
                    selectionBoxActive = false;
                    document.getElementById('svgCanvas').style.cursor = 'default';
                    if (selectionBoxRect) {{
                        const latestPt = e ? draw.point(e.clientX, e.clientY) : null;
                        selectionBoxEnd = latestPt || selectionBoxEnd || selectionBoxStart;
                        selectionBoxEndScreen = e ? {{ x: e.clientX, y: e.clientY }} : selectionBoxEndScreen;
                        const selectionBoundsScreen = normalizeBounds(selectionBoxStartScreen, selectionBoxEndScreen) || getScreenBounds(selectionBoxRect);
                        const selectionBoundsSvg = normalizeBounds(selectionBoxStart, selectionBoxEnd);
                        selectionBoxRect.remove();
                        selectionBoxRect = null;
                        selectionBoxStart = null;
                        selectionBoxEnd = null;
                        selectionBoxStartScreen = null;
                        selectionBoxEndScreen = null;
                        const selectedEntries = [];
                        const hitScreen = (box) => {{
                            if (!box || !selectionBoundsScreen) return false;
                            return box.right >= selectionBoundsScreen.left && box.left <= selectionBoundsScreen.right && box.bottom >= selectionBoundsScreen.top && box.top <= selectionBoundsScreen.bottom;
                        }};
                        Object.keys(protboxMap).forEach(pid => {{
                            const el = draw.findOne(`[data-id=\"${{pid}}\"]`);
                            if (!el) return;
                            let hit = false;
                            if (selectionBoundsSvg) {{
                                const boxSvg = getSvgBounds(el);
                                if (boxSvg) {{
                                    const px1 = boxSvg.left;
                                    const py1 = boxSvg.top;
                                    const px2 = boxSvg.right;
                                    const py2 = boxSvg.bottom;
                                    hit = px2 >= selectionBoundsSvg.left && px1 <= selectionBoundsSvg.right && py2 >= selectionBoundsSvg.top && py1 <= selectionBoundsSvg.bottom;
                                }}
                            }}
                            if (!hit) {{
                                const box = getScreenBounds(el);
                                hit = hitScreen(box);
                            }}
                            if (hit) selectedEntries.push({{ type: 'prot-box', id: pid, protboxId: pid }});
                        }});
                        // Text boxes (including shapes rendered as text-box with shape_type)
                        const textEls = draw.find('[data-type=\"text-box\"]') || [];
                        textEls.forEach(el => {{
                            const tid = el.attr('data-id');
                            if (!tid) return;
                            const box = getScreenBounds(el);
                            if (hitScreen(box)) selectedEntries.push({{ type: 'text-box', id: tid }});
                        }});
                        // Compounds (treated as shapes)
                        const compoundEls = draw.find('[data-type=\"compound\"]') || [];
                        compoundEls.forEach(el => {{
                            const cid = el.attr('data-id');
                            if (!cid) return;
                            const box = getScreenBounds(el);
                            if (hitScreen(box)) selectedEntries.push({{ type: 'compound', id: cid }});
                        }});
                        deselectElement();
                        selectedEntries.forEach((entry, idx) => {{
                            const target = draw.findOne(`[data-id=\"${{entry.id}}\"]`);
                            if (target) selectElement(target, entry.type, entry.id, entry.protboxId || null, {{ additive: idx > 0 }});
                        }});
                suppressNextBackgroundClick = true;
            }}
            return;
        }}
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
                    let arrowElement = draw.find(`[data-id="${{arrowId}}"]`)[0];
                    let arrowHitbox = draw.find(`[data-id="${{arrowId}}_hit"]`)[0];
                    let arrowHead = draw.find(`[data-id="${{arrowId}}_head"]`)[0];
                    if (!arrowElement) {{
                        // Draw on demand if the element has not been created yet
                        drawArrows([arrow], {{ force: true }});
                        arrowElement = draw.find(`[data-id="${{arrowId}}"]`)[0];
                        arrowHitbox = draw.find(`[data-id="${{arrowId}}_hit"]`)[0];
                        arrowHead = draw.find(`[data-id="${{arrowId}}_head"]`)[0];
                        if (!arrowElement) {{
                            console.warn && console.warn('mk: missing arrow element for', arrowId);
                            return;
                        }}
                    }}
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
                    if (arrowHitbox && arrowHitbox.plot) arrowHitbox.plot(x1, y1, x2, y2);
                    const startHandle = arrowHandleGroups[arrowId]?.startHandle;
                    const endHandle = arrowHandleGroups[arrowId]?.endHandle;
                    if (startHandle) startHandle.cx(x1).cy(y1);
                    if (endHandle) endHandle.cx(x2).cy(y2);
                    if (arrowHead && arrowHead.plot) {{
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
            const PROTBOX_SNAP_DISTANCE = 2;
            const PROTBOX_UNSNAP_OFFSET = 5;
            const ALIGNMENT_TOLERANCE = 2;
            const alignmentGuides = [];
            const clearAlignmentGuides = () => {{
                while (alignmentGuides.length) {{
                    const g = alignmentGuides.pop();
                    try {{ g.remove(); }} catch (err) {{}}
                }}
            }};
            const addAlignmentGuide = (x1, y1, x2, y2) => {{
                const line = alignmentGuideGroup.line(x1, y1, x2, y2).stroke({{ color: '#ff4d6a', width: 1, dasharray: '2,2' }});
                alignmentGuides.push(line);
                return line;
            }};
            const hasVerticalLineOfSight = (xLine, y1, y2, ignore = new Set()) => {{
                const minY = Math.min(y1, y2), maxY = Math.max(y1, y2);
                for (const pid in protboxMap) {{
                    if (ignore.has(pid)) continue;
                    const box = protboxMap[pid];
                    if (!box) continue;
                    const bx2 = (box.x || 0) + (box.width || 0);
                    const by2 = (box.y || 0) + (box.height || 0);
                    if (xLine >= (box.x || 0) && xLine <= bx2 && by2 > minY && (box.y || 0) < maxY) {{
                        return false;
                    }}
                }}
                return true;
            }};
            const hasHorizontalLineOfSight = (yLine, x1, x2, ignore = new Set()) => {{
                const minX = Math.min(x1, x2), maxX = Math.max(x1, x2);
                for (const pid in protboxMap) {{
                    if (ignore.has(pid)) continue;
                    const box = protboxMap[pid];
                    if (!box) continue;
                    const bx2 = (box.x || 0) + (box.width || 0);
                    const by2 = (box.y || 0) + (box.height || 0);
                    if (yLine >= (box.y || 0) && yLine <= by2 && bx2 > minX && (box.x || 0) < maxX) {{
                        return false;
                    }}
                }}
                return true;
            }};
            const linkProtboxes = (idA, idB, axis, sideA = null, sideB = null) => {{
                if (!idA || !idB || idA === idB) return;
                protboxLinks[idA] = protboxLinks[idA] || {{}};
                protboxLinks[idB] = protboxLinks[idB] || {{}};
                protboxLinks[idA][idB] = {{ axis: axis || null, side: sideA, otherSide: sideB }};
                protboxLinks[idB][idA] = {{ axis: axis || null, side: sideB, otherSide: sideA }};
            }};
            const unlinkPair = (idA, idB) => {{
                if (!idA || !idB) return;
                if (protboxLinks[idA]) {{
                    delete protboxLinks[idA][idB];
                    if (!Object.keys(protboxLinks[idA]).length) delete protboxLinks[idA];
                }}
                if (protboxLinks[idB]) {{
                    delete protboxLinks[idB][idA];
                    if (!Object.keys(protboxLinks[idB]).length) delete protboxLinks[idB];
                }}
            }};
            const pickSideFromTheta = (thetaDeg) => {{
                if (thetaDeg >= -45 && thetaDeg < 45) return 'East';
                if (thetaDeg >= 45 && thetaDeg < 135) return 'South';
                if (thetaDeg >= -135 && thetaDeg < -45) return 'North';
                return 'West';
            }};
            const oppositeSide = (side) => {{
                const map = {{ East: 'West', West: 'East', North: 'South', South: 'North' }};
                return map[side] || side;
            }};
            const protboxAnchorPoint = (pid, side) => {{
                const pb = protboxMap[pid];
                if (!pb) return {{ x: 0, y: 0 }};
                const dists = protboxHandleDists[pid] || {{ North: 5, South: 5, West: 5, East: 5 }};
                const x = pb.x || 0, y = pb.y || 0, w = pb.width || 0, h = pb.height || 0;
                switch (side) {{
                    case 'North': return {{ x: x + w / 2, y: y - dists.North }};
                    case 'South': return {{ x: x + w / 2, y: y + h + dists.South }};
                    case 'West': return {{ x: x - dists.West, y: y + h / 2 }};
                    case 'East': return {{ x: x + w + dists.East, y: y + h / 2 }};
                    default: return {{ x: x + w / 2, y: y + h / 2 }};
                }}
            }};
            const createInteractionBetweenProtboxes = (idA, idB, type = 'arrow', options = {{}}) => {{
                const a = protboxMap[idA], b = protboxMap[idB];
                if (!a || !b) return;
                const acx = (a.x || 0) + (a.width || 0) / 2;
                const acy = (a.y || 0) + (a.height || 0) / 2;
                const bcx = (b.x || 0) + (b.width || 0) / 2;
                const bcy = (b.y || 0) + (b.height || 0) / 2;
                const theta = Math.atan2(bcy - acy, bcx - acx) * 180 / Math.PI;
                const sideA = pickSideFromTheta(theta);
                const sideB = oppositeSide(sideA);
                const p1 = protboxAnchorPoint(idA, sideA);
                const p2 = protboxAnchorPoint(idB, sideB);
                const newArrow = {{
                    line: type,
                    dashed: options.dashed === true,
                    protbox_id_1: idA,
                    protbox_id_1_side: sideA,
                    protbox_id_2: idB,
                    protbox_id_2_side: sideB,
                    x1: p1.x, y1: p1.y, x2: p2.x, y2: p2.y
                }};
                arrows.push(newArrow);
                drawArrows([newArrow], {{ force: true }});
                rebuildAttachmentsIndex();
                updateArrowPositions(idA, sideA);
                updateArrowPositions(idB, sideB);
                updateAttachedArrows(arrowIdFromIndex(arrows.length - 1), 'start');
                updateAttachedArrows(arrowIdFromIndex(arrows.length - 1), 'end');
                Shiny?.setInputValue('add_arrow', {{
                    line: type,
                    protbox_id_1: idA,
                    protbox_id_1_side: sideA,
                    protbox_id_2: idB,
                    protbox_id_2_side: sideB,
                    x1: p1.x, y1: p1.y, x2: p2.x, y2: p2.y
                }}, {{ priority: 'event' }});
                return arrowIdFromIndex(arrows.length - 1);
            }};
            const findAlignmentSnap = (movingId, movingRect) => {{
                if (!movingId || !movingRect || !shiftKeyDown) return null;
                const ignore = new Set([`${{movingId}}`]);
                const moving = {{
                    x: movingRect.x,
                    y: movingRect.y,
                    width: movingRect.width,
                    height: movingRect.height,
                    cx: movingRect.x + movingRect.width / 2,
                    cy: movingRect.y + movingRect.height / 2
                }};
                let bestX = null;
                let bestY = null;
                Object.keys(protboxMap).forEach(otherId => {{
                    if (otherId === `${{movingId}}`) return;
                    if (!sharesGroup(movingId, otherId)) return;
                    const other = protboxMap[otherId];
                    if (!other) return;
                    const otherGeom = {{
                        x: other.x || 0,
                        y: other.y || 0,
                        width: other.width || 0,
                        height: other.height || 0,
                    }};
                    otherGeom.cx = otherGeom.x + otherGeom.width / 2;
                    otherGeom.cy = otherGeom.y + otherGeom.height / 2;
                    const localIgnore = new Set(ignore);
                    localIgnore.add(`${{otherId}}`);
                    const guideX = (moving.cx + otherGeom.cx) / 2;
                    const xTargets = [
                        {{ movingVal: moving.x, otherVal: otherGeom.x, axis: 'x', y1: moving.cy, y2: otherGeom.cy, guideX }},
                        {{ movingVal: moving.x + moving.width, otherVal: otherGeom.x + otherGeom.width, axis: 'x', y1: moving.cy, y2: otherGeom.cy, guideX }},
                        {{ movingVal: moving.cx, otherVal: otherGeom.cx, axis: 'x', y1: moving.cy, y2: otherGeom.cy, guideX }},
                    ];
                    const guideY = (moving.cy + otherGeom.cy) / 2;
                    const yTargets = [
                        {{ movingVal: moving.y, otherVal: otherGeom.y, axis: 'y', x1: moving.cx, x2: otherGeom.cx, guideY }},
                        {{ movingVal: moving.y + moving.height, otherVal: otherGeom.y + otherGeom.height, axis: 'y', x1: moving.cx, x2: otherGeom.cx, guideY }},
                        {{ movingVal: moving.cy, otherVal: otherGeom.cy, axis: 'y', x1: moving.cx, x2: otherGeom.cx, guideY }},
                    ];
                    xTargets.forEach(t => {{
                        const gap = Math.abs(t.movingVal - t.otherVal);
                        if (gap > ALIGNMENT_TOLERANCE) return;
                        if (!hasVerticalLineOfSight(t.otherVal, t.y1, t.y2, localIgnore)) return;
                        if (!bestX || gap < bestX.gap) {{
                            bestX = {{ gap, delta: t.otherVal - t.movingVal, x: t.guideX || t.otherVal, y1: t.y1, y2: t.y2 }};
                        }}
                    }});
                    yTargets.forEach(t => {{
                        const gap = Math.abs(t.movingVal - t.otherVal);
                        if (gap > ALIGNMENT_TOLERANCE) return;
                        if (!hasHorizontalLineOfSight(t.otherVal, t.x1, t.x2, localIgnore)) return;
                        if (!bestY || gap < bestY.gap) {{
                            bestY = {{ gap, delta: t.otherVal - t.movingVal, y: t.guideY || t.otherVal, x1: t.x1, x2: t.x2 }};
                        }}
                    }});
                }});
                if (!bestX && !bestY) return null;
                const guides = [];
                if (bestX) guides.push({{ x1: bestX.x, y1: bestX.y1, x2: bestX.x, y2: bestX.y2 }});
                if (bestY) guides.push({{ x1: bestY.x1, y1: bestY.y, x2: bestY.x2, y2: bestY.y }});
                return {{
                    dx: bestX ? bestX.delta : 0,
                    dy: bestY ? bestY.delta : 0,
                    guides
                }};
            }};
            const isSideFree = (protboxId, side, allowPeer = null) => {{
                if (!protboxId || !side) return true;
                const links = protboxLinks[protboxId];
                if (!links) return true;
                for (const peer in links) {{
                    if (allowPeer && peer === allowPeer) continue;
                    const meta = links[peer];
                    if (!meta) continue;
                    // Only block the side this protbox is using for that link
                    if (meta.side === side) return false;
                }}
                return true;
            }};
            const linkStillValid = (idA, idB, meta) => {{
                const a = protboxMap[idA];
                const b = protboxMap[idB];
                if (!a || !b || !meta) return false;
                const sideA = meta.side;
                const sideB = meta.otherSide;
                if (!sideA || !sideB) return false;
                const tol = PROTBOX_SNAP_DISTANCE + 0.5;
                const centerTol =  (sideA === 'East' || sideA === 'West') 
                    ? (a.height || 0)/2 + (b.height || 0)/2 + tol 
                    : (a.width || 0)/2 + (b.width || 0)/2 + tol;
                const aCx = a.x + (a.width || 0)/2;
                const aCy = a.y + (a.height || 0)/2;
                const bCx = b.x + (b.width || 0)/2;
                const bCy = b.y + (b.height || 0)/2;
                if (sideA === 'East' && sideB === 'West') {{
                    if (Math.abs((a.x + (a.width||0)) - b.x) > tol) return false;
                    if (Math.abs(aCy - bCy) > centerTol) return false;
                    return true;
                }}
                if (sideA === 'West' && sideB === 'East') {{
                    if (Math.abs(a.x - (b.x + (b.width||0))) > tol) return false;
                    if (Math.abs(aCy - bCy) > centerTol) return false;
                    return true;
                }}
                if (sideA === 'South' && sideB === 'North') {{
                    if (Math.abs((a.y + (a.height||0)) - b.y) > tol) return false;
                    if (Math.abs(aCx - bCx) > centerTol) return false;
                    return true;
                }}
                if (sideA === 'North' && sideB === 'South') {{
                    if (Math.abs(a.y - (b.y + (b.height||0))) > tol) return false;
                    if (Math.abs(aCx - bCx) > centerTol) return false;
                    return true;
                }}
                return false;
            }};
            const cleanupBrokenLinks = () => {{
                const removals = [];
                Object.keys(protboxLinks).forEach(idA => {{
                    const peers = protboxLinks[idA];
                    if (!peers) return;
                    Object.keys(peers).forEach(idB => {{
                        const meta = peers[idB];
                        if (!linkStillValid(idA, idB, meta)) {{
                            removals.push([idA, idB]);
                        }}
                    }});
                }});
                removals.forEach(pair => unlinkPair(pair[0], pair[1]));
            }};
            const unlinkProtbox = (protboxId) => {{
                const links = protboxLinks[protboxId];
                const axes = links ? Object.values(links).map(l => (l && l.axis)).filter(Boolean) : [];
                if (!links) return null;
                Object.keys(links).forEach(peerId => {{
                    if (protboxLinks[peerId]) {{
                        delete protboxLinks[peerId][protboxId];
                        if (!Object.keys(protboxLinks[peerId]).length) delete protboxLinks[peerId];
                    }}
                }});
                delete protboxLinks[protboxId];
                return axes.length ? axes[0] : null;
            }};
            const translateProtbox = (protboxId, deltaX, deltaY, options = {{}}, visitedHandles = new Set()) => {{
                if (!protboxId) return null;
                const rect = options.rect || draw.findOne(`[data-id="${{protboxId}}"]`);
                if (!rect) return null;
                const startX = rect.x();
                const startY = rect.y();
                if (!options.skipRectMove) {{
                    rect.dmove(deltaX, deltaY);
                }}
                const actualDX = rect.x() - startX;
                const actualDY = rect.y() - startY;
                const group = elementGroups[protboxId];
                if (group) {{
                    group.forEach(assocElement => {{
                        const assocType = assocElement.element.attr('data-type');
                        if (assocType === 'prot-box') return;
                        if (assocType === 'handle-line') {{
                            if (visitedHandles.has(assocElement.element)) return;
                            visitedHandles.add(assocElement.element);
                            assocElement.element.plot(
                                parseFloat(assocElement.element.attr('x1') || 0) + actualDX,
                                parseFloat(assocElement.element.attr('y1') || 0) + actualDY,
                                parseFloat(assocElement.element.attr('x2') || 0) + actualDX,
                                parseFloat(assocElement.element.attr('y2') || 0) + actualDY
                            );
                        }} else {{
                            assocElement.element.dmove(actualDX, actualDY);
                        }}
                    }});
                }}
                if (protboxId in attachments) {{
                    Object.keys(attachments[protboxId]).forEach(side => updateArrowPositions(protboxId, side));
                }}
                const spacing = settings.ptm_circle_spacing || 4;
                const width = rect.width();
                const height = rect.height();
                const pbX = rect.x();
                const pbY = rect.y();
                protboxSnapPoints[protboxId] = {{'N1': {{x: pbX + width * 0.2, y: pbY - spacing}},'N2': {{x: pbX + width * 0.5, y: pbY - spacing}},'N3': {{x: pbX + width * 0.8, y: pbY - spacing}},'S1': {{x: pbX + width * 0.2, y: pbY + height + spacing}},'S2': {{x: pbX + width * 0.5, y: pbY + height + spacing}},'S3': {{x: pbX + width * 0.8, y: pbY + height + spacing}},'W1': {{x: pbX - spacing, y: pbY + height * 0.33 - 2}},'W2': {{x: pbX - spacing, y: pbY + height * 0.66 + 2}},'E1': {{x: pbX + width + spacing, y: pbY + height * 0.33 - 2}},'E2': {{x: pbX + width + spacing, y: pbY + height * 0.66 + 2}}}};
                if (group) {{
                    group.filter(el => el.element.attr('data-type') === 'ptm-snap-circle').forEach(snapC => {{
                        const key = snapC.element.attr('data-pos-key');
                        snapC.element.cx(protboxSnapPoints[protboxId][key].x).cy(protboxSnapPoints[protboxId][key].y);
                    }});
                }}
                if (protboxMap[protboxId]) {{
                    protboxMap[protboxId].x = pbX;
                    protboxMap[protboxId].y = pbY;
                }}
                refreshGroupsForProtbox(protboxId);
                return {{ dx: actualDX, dy: actualDY }};
            }};
            const propagateLinkedMove = (anchorId, deltaX, deltaY, visited = new Set(), visitedHandles = new Set()) => {{
                if (!anchorId || visited.has(anchorId)) return;
                visited.add(anchorId);
                const links = protboxLinks[anchorId];
                if (!links) return;
                Object.keys(links).forEach(peerId => {{
                    const rect = draw.findOne(`[data-id="${{peerId}}"]`);
                    const delta = translateProtbox(peerId, deltaX, deltaY, {{ rect }}, visitedHandles) || {{ dx: deltaX, dy: deltaY }};
                    propagateLinkedMove(peerId, delta.dx, delta.dy, visited, visitedHandles);
                }});
            }};
            const maybeSnapProtbox = (movingId, rectRef = null) => {{
                cleanupBrokenLinks();
                const rect = rectRef || draw.findOne(`[data-id="${{movingId}}"]`);
                if (!movingId || !rect) return;
                const movingWidth = rect.width();
                const movingHeight = rect.height();
                const movingX = rect.x();
                const movingY = rect.y();
                const movingCx = movingX + movingWidth / 2;
                const movingCy = movingY + movingHeight / 2;
                let best = null;
                const considerCandidate = (gap, align, axis, otherId, movingCx, movingCy, otherCx, otherCy) => {{
                    if (gap > PROTBOX_SNAP_DISTANCE) return;
                    const centerPenalty = (align === 'left-of' || align === 'right-of') ? Math.abs(movingCy - otherCy) : Math.abs(movingCx - otherCx);
                    const score = gap + centerPenalty * 0.25; // bias toward the nearest, better-aligned box
                    if (!best || score < best.score) {{
                        best = {{ otherId, axis, dist: gap, align, score }};
                    }}
                }};
                Object.keys(protboxMap).forEach(otherId => {{
                    if (otherId === movingId) return;
                    if (!sharesGroup(movingId, otherId)) return;
                    const other = protboxMap[otherId];
                    if (!other) return;
                    const otherWidth = other.width || 0;
                    const otherHeight = other.height || 0;
                    const otherCx = other.x + otherWidth / 2;
                    const otherCy = other.y + otherHeight / 2;
                    // Horizontal adjacency (moving left/right of other)
                    const horizOverlap = Math.abs(movingCy - otherCy) <= (movingHeight + otherHeight) / 2 + PROTBOX_SNAP_DISTANCE;
                    if (horizOverlap) {{
                        const distLeft = Math.abs((movingX + movingWidth) - other.x); // moving right edge to other left
                        const distRight = Math.abs(movingX - (other.x + otherWidth)); // moving left edge to other right
                        considerCandidate(distLeft, 'left-of', 'y', otherId, movingCx, movingCy, otherCx, otherCy);
                        considerCandidate(distRight, 'right-of', 'y', otherId, movingCx, movingCy, otherCx, otherCy);
                    }}
                    // Vertical adjacency (stacked)
                    const vertOverlap = Math.abs(movingCx - otherCx) <= (movingWidth + otherWidth) / 2 + PROTBOX_SNAP_DISTANCE;
                    if (vertOverlap) {{
                        const distAbove = Math.abs((movingY + movingHeight) - other.y); // moving bottom to other top
                        const distBelow = Math.abs(movingY - (other.y + otherHeight)); // moving top to other bottom
                        considerCandidate(distAbove, 'above', 'x', otherId, movingCx, movingCy, otherCx, otherCy);
                        considerCandidate(distBelow, 'below', 'x', otherId, movingCx, movingCy, otherCx, otherCy);
                    }}
                }});
                if (!best) return;
                let extraDX = 0, extraDY = 0;
                const target = protboxMap[best.otherId];
                const targetWidth = target.width || 0;
                const targetHeight = target.height || 0;
                let movingSide = null, targetSide = null;
                if (best.align === 'left-of') {{
                    extraDX = target.x - (movingX + movingWidth);
                    extraDY = (target.y + targetHeight / 2) - movingCy;
                    movingSide = 'East';
                    targetSide = 'West';
                }} else if (best.align === 'right-of') {{
                    extraDX = (target.x + targetWidth) - movingX;
                    extraDY = (target.y + targetHeight / 2) - movingCy;
                    movingSide = 'West';
                    targetSide = 'East';
                }} else if (best.align === 'above') {{
                    extraDY = target.y - (movingY + movingHeight);
                    extraDX = (target.x + targetWidth / 2) - movingCx;
                    movingSide = 'South';
                    targetSide = 'North';
                }} else if (best.align === 'below') {{
                    extraDY = (target.y + targetHeight) - movingY;
                    extraDX = (target.x + targetWidth / 2) - movingCx;
                    movingSide = 'North';
                    targetSide = 'South';
                }} else if (best.axis === 'x') {{
                    const targetCx = target.x + targetWidth / 2;
                    extraDX = targetCx - movingWidth / 2 - movingX;
                    movingSide = 'West';
                    targetSide = 'East';
                }} else {{
                    const targetCy = target.y + targetHeight / 2;
                    extraDY = targetCy - movingHeight / 2 - movingY;
                    movingSide = 'North';
                    targetSide = 'South';
                }}
                // block snapping if either side is already occupied by a different protbox
                if (!isSideFree(movingId, movingSide, best.otherId) || !isSideFree(best.otherId, targetSide, movingId)) {{
                    return;
                }}
                if (extraDX || extraDY) {{
                    const handleVisit = new Set();
                    const delta = translateProtbox(movingId, extraDX, extraDY, {{ rect, skipRectMove: false }}, handleVisit) || {{ dx: extraDX, dy: extraDY }};
                    propagateLinkedMove(movingId, delta.dx, delta.dy, new Set([movingId]), handleVisit);
                }}
                linkProtboxes(movingId, best.otherId, best.axis, movingSide, targetSide);
            }};
            const unsnapAndNudge = (protboxId) => {{
                if (!protboxId) return;
                const rect = draw.findOne(`[data-id="${{protboxId}}"]`);
                if (!rect) return;
                const axis = unlinkProtbox(protboxId);
                if (!axis) return;
                const dx = axis === 'x' ? PROTBOX_UNSNAP_OFFSET : 0;
                const dy = axis === 'y' ? PROTBOX_UNSNAP_OFFSET : 0;
                const mover = () => {{
                    const delta = translateProtbox(protboxId, dx, dy, {{ rect }}) || {{ dx, dy }};
                    propagateLinkedMove(protboxId, delta.dx, delta.dy);
                }};
                if (mkHistory) {{
                    mkHistory.captureInstantMove({{ type: 'prot-box', id: protboxId, protboxId }}, mover);
                }} else {{
                    mover();
                }}
            }};
            const moveSelectedElement = (deltaX, deltaY, pointerX = null, pointerY = null) => {{
                if (!selectedElement) return;
                const memberGroups = !activeGroupEditId ? groupsForEntry(selectedType, selectedId) : null;
                const isGroupDrag = selectionMap.size > 1 || (selectedType === 'group') || (memberGroups && memberGroups.size);
                if (!(selectedType === 'prot-box' && shiftKeyDown && activeDragProtboxId)) {{
                    clearAlignmentGuides();
                }}
                if (selectionMap.size > 1) {{
                    selectionMap.forEach(entry => {{
                        if (!entry || !entry.element) return;
                        if (entry.type === 'group') {{
                            moveGroupMembers(entry.id, deltaX, deltaY, {{ snap: false }});
                        }} else if (entry.type === 'prot-box') {{
                            const handleVisit = new Set();
                            const delta = translateProtbox(entry.id, deltaX, deltaY, {{ rect: entry.element, skipRectMove: false }}, handleVisit) || {{ dx: deltaX, dy: deltaY }};
                            propagateLinkedMove(entry.id, delta.dx, delta.dy, new Set([entry.id]), handleVisit);
                            cleanupBrokenLinks();
                            if (!isGroupDrag && !(shiftKeyDown && activeDragProtboxId === entry.id)) {{
                                maybeSnapProtbox(entry.id, entry.element);
                            }}
                            refreshGroupsForProtbox(entry.id);
                        }} else if (entry.type === 'text-box') {{
                            const tb = findTextBlockByDomId(entry.id);
                            if (tb) {{
                                tb.x = toCoordinateNumber((tb.x || 0) + deltaX, 0);
                                tb.y = toCoordinateNumber((tb.y || 0) + deltaY, 0);
                                applyTextBoxLayout(entry.id, tb);
                                Shiny?.setInputValue('element_moved', {{ type: 'text-box', id: entry.id, x: tb.x, y: tb.y, width: tb.width, height: tb.height }}, {{ priority: 'event' }});
                            }} else {{
                                entry.element.dmove(deltaX, deltaY);
                                positionTextHandles(entry.id);
                            }}
                        }} else {{
                            entry.element.dmove(deltaX, deltaY);
                        }}
                    }});
                    return;
                }}
                const isCircle = selectedType === 'ptm-shape' || selectedType === 'arrow-start' || selectedType === 'arrow-end' || (selectedElement && selectedElement.type === 'circle');
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
                        if (pointerX != null && pointerY != null) {{
                            newX = pointerX;
                            newY = pointerY;
                        }}
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
                if (selectedType === 'group') {{
                    moveGroupMembers(selectedId, deltaX, deltaY, {{ snap: false }});
                    return;
                }}
                if (selectedType === 'prot-box' && shiftKeyDown && activeDragProtboxId === selectedProtboxId) {{
                    const w = selectedElement.width ? selectedElement.width() : 0;
                    const h = selectedElement.height ? selectedElement.height() : 0;
                    const movingRect = {{ x: newX, y: newY, width: w, height: h }};
                    const snap = findAlignmentSnap(selectedProtboxId, movingRect);
                    clearAlignmentGuides();
                    if (snap) {{
                        deltaX += snap.dx;
                        deltaY += snap.dy;
                        newX += snap.dx;
                        newY += snap.dy;
                        (snap.guides || []).forEach(g => addAlignmentGuide(g.x1, g.y1, g.x2, g.y2));
                    }}
                }}
                if (selectedType === 'prot-box' && protboxInAnyGroup(selectedId) && !activeGroupEditId) {{
                    const members = allMembersForEntry('prot-box', selectedId);
                    const visitedHandles = new Set();
                    members.forEach(m => {{
                    if (!m) return;
                    if (m.type === 'prot-box') {{
                        const pid = normalizeProtboxId(m.id);
                        const rect = draw.findOne(`[data-id="${{pid}}"]`);
                        const delta = translateProtbox(pid, deltaX, deltaY, {{ rect, skipRectMove: false }}, visitedHandles) || {{ dx: deltaX, dy: deltaY }};
                        propagateLinkedMove(pid, delta.dx, delta.dy, new Set([pid]), visitedHandles);
                    }} else if (m.type === 'text-box') {{
                        const tb = findTextBlockByDomId(m.id);
                        if (tb) {{
                            tb.x = toCoordinateNumber((tb.x || 0) + deltaX, 0);
                            tb.y = toCoordinateNumber((tb.y || 0) + deltaY, 0);
                                applyTextBoxLayout(m.id, tb);
                                Shiny?.setInputValue('element_moved', {{ type: 'text-box', id: m.id, x: tb.x, y: tb.y, width: tb.width, height: tb.height }}, {{ priority: 'event' }});
                            }} else {{
                                const grp = draw.findOne(`[data-id="${{m.id}}"]`);
                                if (grp) grp.dmove(deltaX, deltaY);
                                positionTextHandles(m.id);
                            }}
                        }} else if (m.type === 'compound' || m.type === 'figure-key') {{
                            const grp = draw.findOne(`[data-id="${{m.id}}"]`);
                            if (grp) grp.dmove(deltaX, deltaY);
                        }} else if (m.type === 'arrow') {{
                            translateArrowByDelta(m.id, deltaX, deltaY);
                        }}
                    }});
                    members.forEach(m => {{
                        if (m.type === 'prot-box') refreshGroupsForProtbox(normalizeProtboxId(m.id));
                    }});
                    return;
                }}
                if (selectedType === 'text-box') {{
                    const tb = findTextBlockByDomId(selectedId);
                    if (tb) {{
                        tb.x = toCoordinateNumber((tb.x || 0) + deltaX, 0);
                        tb.y = toCoordinateNumber((tb.y || 0) + deltaY, 0);
                        applyTextBoxLayout(selectedId, tb);
                        Shiny?.setInputValue('element_moved', {{ type: 'text-box', id: selectedId, x: tb.x, y: tb.y, width: tb.width, height: tb.height }}, {{ priority: 'event' }});
                    }}
                    return;
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
                        const baseLineType = lineType.startsWith('dashed_') ? lineType.replace(/^dashed_/, '') : lineType;
                        const dx = x2 - x1, dy = y2 - y1, angle = Math.atan2(dy, dx), arrowSize = 5;
                        if (baseLineType === 'arrow') {{
                            arrowHead.plot([[x2, y2],[x2 - arrowSize * Math.cos(angle + Math.PI / 6), y2 - arrowSize * Math.sin(angle + Math.PI / 6)],[x2 - arrowSize * Math.cos(angle - Math.PI / 6), y2 - arrowSize * Math.sin(angle - Math.PI / 6)]]);
                        }} else if (baseLineType === 'inhibition') {{
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
                                const baseLineType = lineType.startsWith('dashed_') ? lineType.replace(/^dashed_/, '') : lineType;
                                if (baseLineType === 'arrow') {{
                                    masterHead.plot([[mx2, my2],[mx2 - arrowSize * Math.cos(angle + Math.PI / 6), my2 - arrowSize * Math.sin(angle + Math.PI / 6)],[mx2 - arrowSize * Math.cos(angle - Math.PI / 6), my2 - arrowSize * Math.sin(angle - Math.PI / 6)]]);
                                }} else if (baseLineType === 'inhibition') {{
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
                        if (selectedType !== 'prot-box') {{
                            selectedElement.move(newX, newY);
                        }}
                    }}
                if (selectedType === 'prot-box' && selectedProtboxId) {{
                    const handleVisit = new Set();
                    const delta = translateProtbox(selectedProtboxId, deltaX, deltaY, {{ rect: selectedElement, skipRectMove: false }}, handleVisit) || {{ dx: deltaX, dy: deltaY }};
                    propagateLinkedMove(selectedProtboxId, delta.dx, delta.dy, new Set([selectedProtboxId]), handleVisit);
                    cleanupBrokenLinks();
                    const snapAllowed = !(shiftKeyDown && activeDragProtboxId === selectedProtboxId);
                    if (snapAllowed) {{
                        maybeSnapProtbox(selectedProtboxId, selectedElement);
                    }}
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
                    const arrowData = context?.arrowData || getArrowById(arrowId);
                    if (arrowElement) {{
                        selectedElement?.cx(newX).cy(newY); // keep handle with cursor/snap
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
                            const baseLineType = lineType.startsWith('dashed_') ? lineType.replace(/^dashed_/, '') : lineType;
                            if (baseLineType === 'arrow') {{
                                arrowHead.plot([[x2, y2],[x2 - arrowSize * Math.cos(angle + Math.PI / 6), y2 - arrowSize * Math.sin(angle + Math.PI / 6)],[x2 - arrowSize * Math.cos(angle - Math.PI / 6), y2 - arrowSize * Math.sin(angle - Math.PI / 6)]]);
                            }} else if (baseLineType === 'inhibition') {{
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
                let multiMoveBefore = null;
                element.node.style.cursor = 'pointer';
                element.node.addEventListener('mousedown', e => {{
                    if (e.button === 2) {{
                        return;
                    }}
                    e.preventDefault();
                    e.stopPropagation();
                    if (type === 'text-box' && isTextEditing && activeTextEditorId === id) {{
                        return;
                    }}
                    if (type === 'text-box' && isTextEditing && activeTextEditorId && activeTextEditorId !== id) {{
                        exitTextEditMode(true);
                    }}
                    if (type !== 'text-box' && isTextEditing) {{
                        exitTextEditMode(true);
                    }}
                    isDragging = true;
                    if (type === 'prot-box') {{
                        activeDragProtboxId = protboxId || id;
                    }}
                    const additiveSelect = e.ctrlKey || e.metaKey;
                    let selectTarget = element;
                    let selectType = type;
                    let selectId = id;
                    let selectProtbox = protboxId;
                    if (additiveSelect && (type === 'ptm-shape' || type === 'ptm-label' || type === 'ptm-symbol') && protboxId) {{
                        selectType = 'prot-box';
                        selectId = protboxId;
                        selectProtbox = protboxId;
                        const pbRect = draw.findOne(`[data-id="${{protboxId}}"]`);
                        if (pbRect) selectTarget = pbRect;
                    }} else if (type === 'text-box') {{
                        const grp = draw.findOne('[data-id="' + id + '"]');
                        if (grp) selectTarget = grp;
                    }}
                    selectElement(selectTarget, selectType, selectId, selectProtbox, {{ additive: additiveSelect, toggle: additiveSelect }});
                    if (selectType === 'text-box') {{
                        ensureTextHandles(selectId, findTextBlockByDomId(selectId));
                        toggleTextHandles(selectId, true);
                    }}
                    if (mkHistory && selectionMap.size > 1) {{
                        multiMoveBefore = [];
                        selectionMap.forEach(entry => {{
                            if (!entry || !trackableMoveTypes.has(entry.type)) return;
                            const el = ensureElementForId(entry.id);
                            const before = captureElementPosition(el, entry.type);
                            if (!before) return;
                            const snap = {{
                                id: entry.id,
                                type: entry.type,
                                protboxId: entry.protboxId || (entry.type === 'prot-box' ? entry.id : null),
                                before
                            }};
                            if (snap.protboxId && snap.type === 'prot-box') {{
                                snap.attached = {{ before: captureAttachedArrowSnapshots(snap.protboxId) }};
                            }}
                            multiMoveBefore.push(snap);
                        }});
                    }} else {{
                        multiMoveBefore = null;
                    }}
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
                    if (mkHistory && multiMoveBefore && selectionMap.size > 1) {{
                        const items = [];
                        multiMoveBefore.forEach(prev => {{
                            const el = ensureElementForId(prev.id);
                            const after = captureElementPosition(el, prev.type);
                            if (!after || !positionsDiffer(prev.before, after)) return;
                            const item = {{
                                id: prev.id,
                                type: prev.type,
                                protboxId: prev.protboxId || null,
                                before: prev.before,
                                after
                            }};
                            if (prev.attached && prev.attached.before && prev.protboxId && prev.type === 'prot-box') {{
                                item.attached = {{
                                    before: prev.attached.before,
                                    after: captureAttachedArrowSnapshots(prev.protboxId)
                                }};
                            }}
                            items.push(item);
                        }});
                        if (items.length) {{
                            const entry = {{
                                kind: 'multi-move',
                                items
                            }};
                            entry.handlers = {{
                                undo: () => mkHistory.runWithoutRecording(() => {{
                                    items.slice().reverse().forEach(it => applyPositionForEntry(it, 'before'));
                                }}),
                                redo: () => mkHistory.runWithoutRecording(() => {{
                                    items.forEach(it => applyPositionForEntry(it, 'after'));
                                }})
                            }};
                            mkHistory.recordAction(entry);
                        }}
                        multiMoveBefore = null;
                    }}
                    if (type === 'prot-box') {{
                        activeDragProtboxId = null;
                        clearAlignmentGuides();
                    }}
                    if (mkHistory && trackableMoveTypes.has(type)) {{
                        mkHistory.finalizeMoveSession();
                    }}
                    let currentX, currentY;
                    if (selectedType === 'arrow-start' || selectedType === 'arrow-end') {{
                        const arrowId = selectedId.replace(/_(start|end)$/, '');
                        const end = selectedType.replace('arrow-', '');
                        const handle = selectedType === 'arrow-start' ? arrowHandleGroups[arrowId]?.startHandle : arrowHandleGroups[arrowId]?.endHandle;
                        const arrow = getArrowById(arrowId);
                        const snapped = updateArrowAttachments(arrowId, end, handle.cx(), handle.cy());
                        if (snapped.isSnapped) {{
                            if (arrow) {{
                                if (snapped.protbox_id) {{
                                    arrow.protbox_id = snapped.protbox_id;
                                    arrow.side = snapped.side;
                                }} else if (snapped.arrow_id) {{
                                    arrow.attached_arrow = snapped.arrow_id;
                                    arrow.attached_end = snapped.arrow_end;
                                }}
                            }}
                        }}
                        updateAttachedArrows(arrowId, end);
                        if (handle) {{
                            Shiny?.setInputValue('arrow_moved', {{ id: arrowId, end, x: handle.cx(), y: handle.cy() }}, {{ priority: 'event' }});
                        }}
                        const lineEl = draw.findOne(`[data-id="${{arrowId}}"]`);
                        if (arrow && lineEl) {{
                            arrow.x1 = parseFloat(lineEl.attr('x1') || handle.cx());
                            arrow.y1 = parseFloat(lineEl.attr('y1') || handle.cy());
                            arrow.x2 = parseFloat(lineEl.attr('x2') || handle.cx());
                            arrow.y2 = parseFloat(lineEl.attr('y2') || handle.cy());
                        }}
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
                }} else if (selectedType === 'group') {{
                    const ids = collectProtboxesForGroup(selectedId);
                    ids.forEach(pid => {{
                        const rect = draw.findOne(`[data-id="${{pid}}"]`);
                        if (!rect) return;
                        const payload = {{ type: 'prot-box', id: pid, x: parseFloat(rect.x() || 0), y: parseFloat(rect.y() || 0), protbox_id: pid }};
                        Shiny?.setInputValue('element_moved', payload, {{ priority: 'event' }});
                        refreshGroupsForProtbox(pid);
                    }});
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
                    if (selectedType === 'text-box') {{
                        const tb = findTextBlockByDomId(selectedId);
                        if (tb) {{
                            currentX = toCoordinateNumber(tb.x, currentX);
                            reportedY = toCoordinateNumber(tb.y, reportedY);
                            eventPayload.x = currentX;
                            eventPayload.y = reportedY;
                            if (tb.width !== undefined) eventPayload.width = tb.width;
                            if (tb.height !== undefined) eventPayload.height = tb.height;
                        }}
                    }}
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
            const updateArrowVisual = (arrowId) => {{
                const arrow = getArrowById(arrowId);
                const line = draw.findOne(`[data-id="${{arrowId}}"]`);
                const hit = draw.findOne(`[data-id="${{arrowId}}_hit"]`);
                let head = draw.findOne(`[data-id="${{arrowId}}_head"]`);
                if (!arrow || !line) return;
                const type = arrow.line || 'arrow';
                const baseType = type.startsWith('dashed_') ? type.replace(/^dashed_/, '') : type;
                const strokeOpts = {{ color: 'black', width: 1 }};
                if (type.startsWith('dashed') || arrow.dashed) strokeOpts.dasharray = '5,5'; else strokeOpts.dasharray = null;
                line.attr('data-line-type', type);
                line.stroke(strokeOpts);
                if (strokeOpts.dasharray) {{
                    line.attr('stroke-dasharray', strokeOpts.dasharray);
                }} else {{
                    line.attr('stroke-dasharray', null);
                }}
                if (hit) hit.stroke({{ color: 'transparent', width: parseFloat(hit.attr('stroke-width') || 22) }});
                if (head) {{ try {{ head.remove(); }} catch (err) {{}} head = null; }}
                const x1 = parseFloat(line.attr('x1') || 0);
                const y1 = parseFloat(line.attr('y1') || 0);
                const x2 = parseFloat(line.attr('x2') || 0);
                const y2 = parseFloat(line.attr('y2') || 0);
                const dx = x2 - x1, dy = y2 - y1, angle = Math.atan2(dy, dx), arrowSize = 5;
                if (baseType === 'arrow') {{
                    head = arrowGroup.polygon([[x2, y2],[x2 - arrowSize * Math.cos(angle + Math.PI / 6), y2 - arrowSize * Math.sin(angle + Math.PI / 6)],[x2 - arrowSize * Math.cos(angle - Math.PI / 6), y2 - arrowSize * Math.sin(angle - Math.PI / 6)]]).fill('black').stroke({{ color: 'black', width: 1 }}).attr({{ 'data-id': `${{arrowId}}_head`, 'data-type': 'arrow-head' }});
                }} else if (baseType === 'inhibition') {{
                    head = arrowGroup.line(x2 - arrowSize, y2, x2 + arrowSize, y2).stroke({{ color: 'black', width: 2 }}).attr({{ 'data-id': `${{arrowId}}_head`, 'data-type': 'arrow-head' }});
                }}
                return head;
            }};
            const rebuildAttachmentsIndex = () => {{
                attachments = {{}};
                arrows.forEach((arrow, index) => {{
                    if (!arrow) return;
                    const arrowId = arrowIdFromIndex(index);
                    if (!arrowId) return;
                    const pushAtt = (pid, side, type) => {{
                        if (!pid || !side || !protboxMap[pid]) return;
                        attachments[pid] = attachments[pid] || {{}};
                        attachments[pid][side] = attachments[pid][side] || [];
                        attachments[pid][side].push({{ type, arrow, arrowId, side }});
                    }};
                    pushAtt(arrow.protbox_id_1, arrow.protbox_id_1_side, 'start');
                    pushAtt(arrow.protbox_id_2, arrow.protbox_id_2_side, 'end');
                }});
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
            const deleteFigureKeyById = (keyId) => {{
                if (!keyId) return false;
                removeElementByDataId(keyId);
                if (elementGroups[keyId]) {{
                    delete elementGroups[keyId];
                }}
                return true;
            }};
            const deleteTextBoxById = (textId, options = {{}}) => {{
                if (!textId) return false;
                removeElementByDataId(textId);
                removeElementByDataId(`${{textId}}_label`);
                removeTextHandles(textId);
                if (activeTextEditorId === textId) {{
                    exitTextEditMode(false);
                }}
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
                        if (Array.isArray(group?.protbox_ids)) {{
                            group.protbox_ids = group.protbox_ids.filter(id => normalizeProtboxId(id) !== targetId);
                        }}
                        if (Array.isArray(group?.members)) {{
                            group.members = group.members.filter(m => !(m && m.type === 'prot-box' && normalizeProtboxId(m.id) === targetId));
                        }}
                    }});
                }}
                if (Object.keys(groupMap).length) {{
                    const toRemove = [];
                    Object.keys(groupMap).forEach(gid => {{
                        const g = groupMap[gid];
                        if (!g) return;
                        if (Array.isArray(g.members)) {{
                            g.members = g.members.filter(m => !(m && m.type === 'prot-box' && normalizeProtboxId(m.id) === targetId));
                        }}
                        g.protbox_ids = (g.protbox_ids || []).filter(pid => normalizeProtboxId(pid) !== targetId);
                        if (!g.members || !g.members.length) {{
                            toRemove.push(gid);
                        }}
                    }});
                    toRemove.forEach(gid => ungroupGroup(gid));
                    rebuildGroupIndexes();
                    Object.keys(groupMap).forEach(gid => refreshGroupVisual(gid));
                    if (Array.isArray(groups)) {{
                        groups.length = 0;
                        Object.values(groupMap).forEach(g => groups.push(g));
                    }}
                    updateGroupSelectionDisplay();
                    syncGroupsToServer();
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
                }} else if (targetType === 'group') {{
                    ungroupGroup(targetId);
                    return;
                }}
                deselectElement();
                let deleted = false;
                if (targetType === 'arrow' || targetType === 'arrow-start' || targetType === 'arrow-end') {{
                    const arrowId = targetType === 'arrow' ? targetId : targetId.replace(/_(start|end)$/, '');
                    deleted = deleteArrowById(arrowId);
                }} else if (targetType === 'compound') {{
                    deleted = deleteCompoundById(targetId);
                }} else if (targetType === 'figure-key') {{
                    deleted = deleteFigureKeyById(targetId);
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
                if (e.key === 'Shift') {{
                    shiftKeyDown = true;
                }}
                const targetEditable = isEditableTarget(e.target);
                if (isTextEditing && activeTextEditorId) {{
                    if (e.key === 'Escape') {{
                        e.preventDefault();
                        exitTextEditMode(true);
                        deselectElement();
                    }}
                    return;
                }}
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
                const selectedProtIds = getSelectedProtboxes();
                const selectedArrowId = (selectedType === 'arrow' || selectedType === 'arrow-hitbox')
                    ? selectedId
                    : ((selectedType === 'arrow-start' || selectedType === 'arrow-end') ? selectedId.replace(/_(start|end)$/, '') : null);
                const applyToArrowId = (arrowId, action) => {{
                    if (!arrowId) return null;
                    const arrow = getArrowById(arrowId);
                    if (!arrow) return null;
                    const result = action(arrowId, arrow);
                    return result === undefined ? arrowId : result;
                }};
                const toggleDashForArrowId = (arrowId) => applyToArrowId(arrowId, () => setArrowDashState(arrowId, !getArrowById(arrowId)?.dashed));
                const flipArrowId = (arrowId) => applyToArrowId(arrowId, (aid, ar) => {{
                    const line = draw.findOne(`[data-id=\"${{aid}}\"]`);
                    const hit = draw.findOne(`[data-id=\"${{aid}}_hit\"]`);
                    if (line) {{
                        const x1 = line.attr('x1'), y1 = line.attr('y1'), x2 = line.attr('x2'), y2 = line.attr('y2');
                        line.attr({{ x1: x2, y1: y2, x2: x1, y2: y1 }});
                        if (hit) hit.attr({{ x1: x2, y1: y2, x2: x1, y2: y1 }});
                    }}
                    const swapFields = (obj, a, b) => {{ const tmp = obj[a]; obj[a] = obj[b]; obj[b] = tmp; }};
                    swapFields(ar, 'x1', 'x2'); swapFields(ar, 'y1', 'y2');
                    swapFields(ar, 'protbox_id_1', 'protbox_id_2');
                    swapFields(ar, 'protbox_id_1_side', 'protbox_id_2_side');
                    swapFields(ar, 'attached_arrow_1', 'attached_arrow_2');
                    swapFields(ar, 'attached_end_1', 'attached_end_2');
                    updateArrowVisual(aid);
                    rebuildAttachmentsIndex();
                    updateAttachedArrows(aid, 'start');
                    updateAttachedArrows(aid, 'end');
                }});
                const cycleTypeForArrowId = (arrowId) => applyToArrowId(arrowId, (aid, ar) => {{
                    const order = ['arrow', 'inhibition', 'line'];
                    const cur = ar.line === 'dashed_arrow' ? 'arrow' : (ar.line || 'arrow');
                    const idx = order.indexOf(cur);
                    const next = order[(idx + 1) % order.length];
                    ar.line = next;
                    updateArrowVisual(aid);
                    translateArrowByDelta(aid, 0, 0);
                    rebuildAttachmentsIndex();
                    updateAttachedArrows(aid, 'start');
                    updateAttachedArrows(aid, 'end');
                }});
                autoConnectSelectedProtboxes = () => {{
                    const protIds = getSelectedProtboxes();
                    if (!Array.isArray(protIds) || protIds.length < 2) return 0;
                    let created = 0;
                    const createdEntries = [];
                    const nodes = protIds.map(id => {{
                        const pb = protboxMap[id];
                        if (!pb) return null;
                        const w = pb.width || 0, h = pb.height || 0;
                        return {{ id, x: pb.x || 0, y: pb.y || 0, h, w, cx: (pb.x || 0) + w / 2, cy: (pb.y || 0) + h / 2 }};
                    }}).filter(Boolean);
                    if (!nodes.length) return 0;
                    const deg = {{}};
                    const connSet = new Set();
                    const addConnKey = (a, b) => connSet.add(a < b ? `${{a}}|${{b}}` : `${{b}}|${{a}}`);
                    const hasConn = (a, b) => connSet.has(a < b ? `${{a}}|${{b}}` : `${{b}}|${{a}}`) || !!findArrowBetweenProtboxes(a, b);
                    (arrows || []).forEach(ar => {{
                        if (!ar) return;
                        const a = ar.protbox_id_1, b = ar.protbox_id_2;
                        if (a && protIds.includes(a)) deg[a] = (deg[a] || 0) + 1;
                        if (b && protIds.includes(b)) deg[b] = (deg[b] || 0) + 1;
                        if (protIds.includes(a) && protIds.includes(b)) addConnKey(a, b);
                    }});
                    const recordConn = (fromId, toId) => {{
                        addConnKey(fromId, toId);
                        deg[fromId] = (deg[fromId] || 0) + 1;
                        deg[toId] = (deg[toId] || 0) + 1;
                    }};
                    const canAttach = (node) => (deg[node.id] || 0) < 2;
                    const tryConnect = (from, to) => {{
                        if (!from || !to) return false;
                        if (to.cx <= from.cx - 5) return false;
                        if (!canAttach(from) || !canAttach(to)) return false;
                        if (hasConn(from.id, to.id)) return false;
                        const newId = createInteractionBetweenProtboxes(from.id, to.id, 'arrow');
                        if (newId) {{
                            const idx = parseArrowIndex(newId);
                            const snap = idx !== null ? cloneData(arrows[idx]) : null;
                            if (snap) createdEntries.push({{ id: newId, idx, snapshot: snap }});
                        }}
                        recordConn(from.id, to.id);
                        created += 1;
                        return true;
                    }};
                    // Cluster into rows by vertical proximity
                    const rowThreshold = 30;
                    const rows = [];
                    [...nodes].sort((a, b) => a.cy - b.cy).forEach(n => {{
                        const row = rows.find(r => Math.abs(r.centerY - n.cy) <= rowThreshold);
                        if (row) {{
                            row.nodes.push(n);
                            row.centerY = row.nodes.reduce((s, v) => s + v.cy, 0) / row.nodes.length;
                        }} else {{
                            rows.push({{ nodes: [n], centerY: n.cy }});
                        }}
                    }});
                    rows.forEach(r => r.nodes.sort((a, b) => a.cx - b.cx || a.cy - b.cy));
                    // Connect within each row only
                    rows.forEach(r => {{
                        const chain = r.nodes;
                        for (let i = 0; i < chain.length - 1; i++) {{
                            tryConnect(chain[i], chain[i + 1]);
                        }}
                        // Ensure each node has at least one connection in its row
                        chain.forEach((n, idx) => {{
                            if ((deg[n.id] || 0) === 0) {{
                                const right = idx < chain.length - 1 ? chain[idx + 1] : null;
                                const left = idx > 0 ? chain[idx - 1] : null;
                                if (right && tryConnect(n, right)) return;
                                if (left) tryConnect(left, n);
                            }}
                        }});
                        // Optional second connection in-row (still max degree 2)
                        chain.forEach((n, idx) => {{
                            if (!canAttach(n)) return;
                            const right = idx < chain.length - 1 ? chain[idx + 1] : null;
                            const left = idx > 0 ? chain[idx - 1] : null;
                            if (right && tryConnect(n, right)) return;
                            if (left) tryConnect(left, n);
                        }});
                    }});
                    try {{
                        console.log('mk: autoConnect created', created, 'arrows for', protIds);
                        Shiny?.setInputValue('auto_connect_result', {{ count: created, protboxes: protIds }}, {{ priority: 'event' }});
                        if (mkHistory && createdEntries.length) {{
                            const entry = {{
                                kind: 'auto-connect',
                                arrows: cloneData(createdEntries)
                            }};
                            entry.handlers = {{
                                undo: () => mkHistory.runWithoutRecording(() => {{
                                    createdEntries.slice().reverse().forEach(e => deleteArrowById(e.id, {{ silent: true }}));
                                    rebuildAttachmentsIndex();
                                }}),
                                redo: () => mkHistory.runWithoutRecording(() => {{
                                    createdEntries.forEach(e => {{
                                        if (typeof e.idx === 'number' && e.idx >= arrows.length) {{
                                            while (arrows.length <= e.idx) arrows.push(null);
                                        }}
                                        const clone = cloneData(e.snapshot);
                                        arrows[e.idx] = clone;
                                        drawArrows([clone], {{ force: true }});
                                        rebuildAttachmentsIndex();
                                        if (clone.protbox_id_1 && clone.protbox_id_1_side) updateArrowPositions(clone.protbox_id_1, clone.protbox_id_1_side);
                                        if (clone.protbox_id_2 && clone.protbox_id_2_side) updateArrowPositions(clone.protbox_id_2, clone.protbox_id_2_side);
                                        const aid = arrowIdFromIndex(e.idx);
                                        updateAttachedArrows(aid, 'start');
                                        updateAttachedArrows(aid, 'end');
                                    }});
                                }})
                            }};
                            mkHistory.recordAction(entry);
                        }}
                    }} catch (err) {{}}
                    return created;
                }};
                const ensureArrowBetweenSelected = (type, opts = {{}}) => {{
                    if (selectedProtIds.length !== 2) return null;
                    const [a, b] = selectedProtIds;
                    let arrowId = findArrowBetweenProtboxes(a, b);
                    if (!arrowId) {{
                        createInteractionBetweenProtboxes(a, b, type, opts);
                        arrowId = arrowIdFromIndex(arrows.length - 1);
                    }} else {{
                        const arrow = getArrowById(arrowId);
                        if (arrow) {{
                            arrow.line = type;
                            if (opts.dashed !== undefined) arrow.dashed = opts.dashed;
                            updateArrowVisual(arrowId);
                            translateArrowByDelta(arrowId, 0, 0);
                            rebuildAttachmentsIndex();
                            updateAttachedArrows(arrowId, 'start');
                            updateAttachedArrows(arrowId, 'end');
                        }}
                    }}
                    return arrowId;
                }};
                if (e.ctrlKey && !targetEditable) {{
                    if (keyLower === 'a' && selectedProtIds.length >= 2) {{ e.preventDefault(); autoConnectSelectedProtboxes(); return; }}
                    if (keyLower === 'i' && selectedProtIds.length === 2) {{ e.preventDefault(); ensureArrowBetweenSelected('inhibition'); return; }}
                    if (keyLower === 'l' && selectedProtIds.length === 2) {{ e.preventDefault(); ensureArrowBetweenSelected('line'); return; }}
                    if (keyLower === 'd') {{
                        if (selectedProtIds.length === 2) {{
                            e.preventDefault();
                            const aid = findArrowBetweenProtboxes(selectedProtIds[0], selectedProtIds[1]);
                            if (aid) toggleDashForArrowId(aid);
                            return;
                        }} else if (selectedArrowId) {{
                            e.preventDefault();
                            toggleDashForArrowId(selectedArrowId);
                            return;
                        }}
                    }}
                    if (keyLower === 'f') {{
                        if (selectedProtIds.length === 2) {{
                            e.preventDefault();
                            const aid = findArrowBetweenProtboxes(selectedProtIds[0], selectedProtIds[1]);
                            if (aid) flipArrowId(aid);
                            return;
                        }} else if (selectedArrowId) {{
                            e.preventDefault();
                            flipArrowId(selectedArrowId);
                            return;
                        }}
                    }}
                    if (keyLower === 's') {{
                        if (selectedProtIds.length === 2) {{
                            e.preventDefault();
                            const aid = findArrowBetweenProtboxes(selectedProtIds[0], selectedProtIds[1]) || ensureArrowBetweenSelected('arrow');
                            if (aid) cycleTypeForArrowId(aid);
                            return;
                        }} else if (selectedArrowId) {{
                            e.preventDefault();
                            cycleTypeForArrowId(selectedArrowId);
                            return;
                        }}
                    }}
                }}
                if (e.key === 'Escape') {{
                    e.preventDefault();
                    exitGroupEditMode();
                    deselectElement();
                    return;
                }}
                if (!selectedElement) return;
                if (e.key === 'Backspace' || e.key === 'Delete') {{
                    const isShapeSelection = selectedType === 'text-box' && findTextBlockByDomId(selectedId)?.shape_type;
                    if (!isShapeSelection && isEditableTarget(e.target)) return;
                    e.preventDefault();
                    deleteSelectedElement();
                    return;
                }}
                if (selectedType === 'text-box' && e.key === 'Enter') {{
                    if (isEditableTarget(e.target)) return;
                    e.preventDefault();
                    enterTextEditMode(selectedId, {{ selectAll: false }});
                    return;
                }}
                let deltaX = 0, deltaY = 0, moveAmount = e.shiftKey ? 10 : 1;
                if (e.key === 'ArrowUp') {{ deltaY = -moveAmount; e.preventDefault(); }}
                else if (e.key === 'ArrowDown') {{ deltaY = moveAmount; e.preventDefault(); }}
                else if (e.key === 'ArrowLeft') {{ deltaX = -moveAmount; e.preventDefault(); }}
                else if (e.key === 'ArrowRight') {{ deltaX = moveAmount; e.preventDefault(); }}
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
            document.addEventListener('keyup', e => {{
                if (e.key === 'Shift') {{
                    shiftKeyDown = false;
                    clearAlignmentGuides();
                }}
            }});
            draw.node.addEventListener('click', e => {{
                if (suppressNextBackgroundClick) {{
                    suppressNextBackgroundClick = false;
                    return;
                }}
                if (isBackgroundTarget(e.target)) {{
                    if (isTextEditing) {{
                        exitTextEditMode(true);
                    }}
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
                const shapeIdToRemove = ptmElementId(protboxId, uniprot, ptmKey, 'shape');
                const labelIdToRemove = ptmElementId(protboxId, uniprot, ptmKey, 'label');
                const symbolIdToRemove = ptmElementId(protboxId, uniprot, ptmKey, 'symbol');
                const idsToRemove = [
                    shapeIdToRemove,
                    labelIdToRemove,
                    symbolIdToRemove
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
                const shapeId = ptmElementId(protboxId, uniprot, ptmKey, 'shape');
                shapeObj.attr({{
                    'data-id': shapeId,
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
                makeDraggable(shapeObj, 'ptm-shape', shapeId, protboxId);
                bindPtmContextMenu(shapeObj, protboxId, ptmMeta);
                const resetLabel = options?.resetLabelPosition || false;
                const resetSymbol = options?.resetSymbolPosition || false;
                const resolveCoord = (value, fallback, forceReset) => (forceReset || value === undefined || value === null) ? fallback : value;
                const forceCenterPlacement = forceCenterSpawn && !existingOverride;
                const labelPosKey = (posKey || ptm.ptm_position || '');
                const ld = labelDefaults[labelPosKey] || [];
                const defaultLabelX = newX + (typeof ld[0] === 'number' ? ld[0] : 0);
                const defaultLabelY = newY + (typeof ld[1] === 'number' ? ld[1] : 0);
                const labelFallbackY = forceCenterPlacement ? (newY - 10 - textOffsetY) : defaultLabelY;
                const labelSnapshot = snapshotSeed?.label;
                const symbolSnapshot = snapshotSeed?.symbol;
                const labelCenteringValue = labelSnapshot?.center || existingOverride?.label_centering || ptm.label_centering || (ld[2] || 'center');
                const labelX = (labelSnapshot && labelSnapshot.x !== undefined)
                    ? labelSnapshot.x
                    : resolveCoord(existingOverride?.label_x ?? ptm.label_x, defaultLabelX, resetLabel || forceCenterPlacement);
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
                        'data-id': ptmElementId(protboxId, uniprot, ptmKey, 'label'),
                        'data-type': 'ptm-label',
                        'data-protbox-id': protboxId
                    }});
                workingEntries.push({{ element: labelText, offsetX: labelX - pbX, offsetY: (labelY + textOffsetY) - pbY, jsonX: labelX, jsonY: labelY }});
                makeDraggable(labelText, 'ptm-label', ptmElementId(protboxId, uniprot, ptmKey, 'label'), protboxId);
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
                    'data-id': ptmElementId(protboxId, uniprot, ptmKey, 'symbol'),
                    'data-type': 'ptm-symbol',
                    'pointer-events': 'none',
                    'data-protbox-id': protboxId
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
                    let label = debugMode ? String(id) : (protein.label || pb.backup_label || 'Unknown');
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
                                'data-id': ptmElementId(id, selectedUniprot, ptm_key, 'shape'),
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
                            makeDraggable(shapeObj, 'ptm-shape', ptmElementId(id, selectedUniprot, ptm_key, 'shape'), id);
                            bindPtmContextMenu(shapeObj, id, ptmMeta);
                            const ptm_pos = effectivePosKey || '';
                            if (ptm_pos.startsWith('N')) northHas = true;
                            if (ptm_pos.startsWith('S')) southHas = true;
                            if (ptm_pos.startsWith('W')) westHas = true;
                            if (ptm_pos.startsWith('E')) eastHas = true;
                            if (ptm.label) {{
                                const labelColor = Array.isArray(ptm.label_color) && ptm.label_color.length === 3 && ptm.label_color.every(c => typeof c === 'number' && c >= 0 && c <= 255) ? `rgb(${{ptm.label_color.join(',')}})` : 'black';
                                const ld = labelDefaults[ptm_pos] || [];
                                const fallbackLabelX = toFiniteNumber(ptm.label_x);
                                const fallbackLabelY = toFiniteNumber(ptm.label_y);
                                const baseLabelX = fallbackLabelX !== null ? fallbackLabelX : ptmX + (typeof ld[0] === 'number' ? ld[0] : 0);
                                const baseLabelY = fallbackLabelY !== null ? fallbackLabelY : ptmY + (typeof ld[1] === 'number' ? ld[1] : 0);
                                const labelX = resolveCoordinate(override?.label_x, baseLabelX) ?? baseLabelX;
                                const labelY = resolveCoordinate(override?.label_y, baseLabelY) ?? baseLabelY;
                                const adjLabelY = Math.round(labelY * boxYStretch);
                                const labelCenteringValue = override?.label_centering || ptm.label_centering || (ld[2] || 'center');
                                const labelCenter = labelCenteringValue.toLowerCase();
                                const textAnchor = anchorMap[labelCenter] || 'middle';
                                const labelText = protboxGroup.text(ptm.label).move(labelX, adjLabelY + textOffsetY).font({{ size: settings.ptm_label_size || 10, family: settings.ptm_label_font || 'Arial', anchor: textAnchor, leading: '1.2em' }}).fill(labelColor).attr({{ 'data-id': ptmElementId(id, selectedUniprot, ptm_key, 'label'), 'data-type': 'ptm-label', 'data-protbox-id': id }});
                                elementGroups[id].push({{ element: labelText, offsetX: labelX - x, offsetY: (labelY + textOffsetY) - y, jsonX: labelX, jsonY: labelY }});
                                makeDraggable(labelText, 'ptm-label', ptmElementId(id, selectedUniprot, ptm_key, 'label'), id);
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
                                const symbolText = protboxGroup.text(ptm.symbol).move(symbolX, adjSymbolY + textOffsetY).font({{ size: ptm.symbol_size || settings.ptm_label_size || 10, family: ptm.symbol_font || settings.ptm_label_font || 'Arial', anchor: 'middle', leading: '1.2em' }}).fill(symbolColor).attr({{ 'data-id': ptmElementId(id, selectedUniprot, ptm_key, 'symbol'), 'data-type': 'ptm-symbol', 'pointer-events': 'none', 'data-protbox-id': id }});
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
                        if (shouldShowGroupMenu('prot-box', id)) {{
                            showGroupingMenu(e, {{ mode: 'create' }});
                            return;
                        }}
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
                        menu.style.padding = '2px 0';
                menu.style.fontSize = '12px';
                        menu.style.zIndex = '1000';
                        const ul = document.createElement('ul');
                        ul.style.listStyle = 'none';
                        ul.style.margin = '0';
                        ul.style.padding = '0';
                        const liProtein = document.createElement('li');
                        liProtein.textContent = 'Protein';
                        liProtein.style.padding = '4px 8px';
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
                            subLi.style.padding = '4px 8px';
                            subLi.style.cursor = 'pointer';
                            subLi.style.display = 'flex';
                            subLi.style.alignItems = 'center';
                            subLi.style.gap = '8px';
                            const swatch = document.createElement('span');
                            swatch.style.display = 'inline-block';
                            swatch.style.width = '12px';
                            swatch.style.height = '12px';
                            swatch.style.border = '1px solid #ccc';
                            swatch.style.flexShrink = '0';
                            try {{
                                swatch.style.backgroundColor = entityColor(prot) || '#ccc';
                            }} catch (err) {{
                                swatch.style.backgroundColor = '#ccc';
                            }}
                            const labelSpan = document.createElement('span');
                            labelSpan.textContent = prot.label || uniprot;
                            subLi.appendChild(swatch);
                            subLi.appendChild(labelSpan);
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
                        liPTMs.style.padding = '4px 8px';
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
                            subLi.style.padding = '4px 8px';
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
                        liLinks.style.padding = '4px 8px';
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
                            linkLi.style.padding = '4px 8px';
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
                const spacing = settings.ptm_circle_spacing || 4;
                const pbX = pb.x || 0;
                const pbY = pb.y || 0;
                protboxSnapPoints[id] = {{
                    'N1': {{ x: pbX + width * 0.2, y: pbY - spacing }},
                    'N2': {{ x: pbX + width * 0.5, y: pbY - spacing }},
                    'N3': {{ x: pbX + width * 0.8, y: pbY - spacing }},
                    'S1': {{ x: pbX + width * 0.2, y: pbY + height + spacing }},
                    'S2': {{ x: pbX + width * 0.5, y: pbY + height + spacing }},
                    'S3': {{ x: pbX + width * 0.8, y: pbY + height + spacing }},
                    'W1': {{ x: pbX - spacing, y: pbY + height * 0.33 - 2 }},
                    'W2': {{ x: pbX - spacing, y: pbY + height * 0.66 + 2 }},
                    'E1': {{ x: pbX + width + spacing, y: pbY + height * 0.33 - 2 }},
                    'E2': {{ x: pbX + width + spacing, y: pbY + height * 0.66 + 2 }},
                }};
            }};
            protBoxes.forEach((pb, index) => {{
                renderProtbox(pb, index);
                registerProtboxBaseGeometry(pb);
            }});
            const autoSpawnPtmsForProtbox = (pb) => {{
                if (!pb || !pb.protbox_id) return;
                const proteins = Array.isArray(pb.proteins) ? pb.proteins : [];
                // Only spawn PTMs for the protein currently shown in the protbox
                const selectedUni = pb.selected_uniprot || proteins[0];
                if (!selectedUni) return;
                const maxPtm = Number(settings.ptm_max_display || 0);
                const snapOrder = ptmPositionPriority;
                const prot = proteinData && proteinData[selectedUni];
                const ptms = prot && prot.PTMs;
                if (!ptms || typeof ptms !== 'object') return;
                const entries = Object.entries(ptms);
                if (!entries.length) return;
                const limit = maxPtm > 0 ? Math.min(maxPtm, entries.length) : entries.length;
                let idx = 0;
                for (let i = 0; i < limit; i++) {{
                    const [ptmKey] = entries[i];
                    const shapeId = ptmElementId(pb.protbox_id, selectedUni, ptmKey, 'shape');
                    if (draw.find(`[data-id="${{shapeId}}"]`).length) {{
                        continue;
                    }}
                    const posKey = snapOrder[idx % snapOrder.length];
                    idx += 1;
                    spawnPtmForProtbox(pb.protbox_id, selectedUni, ptmKey, {{
                        preferredSnapKey: posKey,
                        spawnInCenter: false,
                        recordHistory: false
                    }});
                }}
            }};
            protBoxes.forEach(autoSpawnPtmsForProtbox);
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
                textBlock.width = toCoordinateNumber(textBlock.width, 160);
                textBlock.height = toCoordinateNumber(textBlock.height, 80);
                textBlock.x = toCoordinateNumber(textBlock.x, 0);
                textBlock.y = toCoordinateNumber(textBlock.y, 0);
                textBlock.bgcolor = textBlock.bgcolor || 'transparent';
                textBlock.fgcolor = textBlock.fgcolor || '#000000';
                textBlock.border_color = textBlock.border_color || '#000000';
                textBlock.border_width = textBlock.border_width || 1;
                const style = resolveTextStyle(textBlock);
                textBlock.text_style = style;
                const isShape = !!textBlock.shape_type;
                const htmlContent = (textBlock.html !== undefined) ? textBlock.html : (textBlock.label || '');
                const parentGroup = isShape ? shapeGroup : textGroup;
                const g = parentGroup.group().attr({{ 'data-id': id, 'data-type': 'text-box' }});
                const hitRect = g.rect(textBlock.width + TEXT_HIT_PADDING * 2, textBlock.height + TEXT_HIT_PADDING * 2)
                    .move(textBlock.x - TEXT_HIT_PADDING, textBlock.y - TEXT_HIT_PADDING)
                    .fill({{ color: '#000', opacity: 0.001 }})
                    .stroke({{ width: 0, color: 'transparent' }});
                hitRect.back().attr({{ 'data-role': 'text-hit' }});
                hitRect.node.style.pointerEvents = 'all';
                hitRect.back();
                const isTransparentShape = isShape && (!textBlock.bgcolor || textBlock.bgcolor === 'transparent' || textBlock.bgcolor === 'none');
                const outlineRadius = (() => {{
                    if (textBlock.shape_type === 'rounded') return 12;
                    if (textBlock.shape_type === 'circle') return Math.max(textBlock.width, textBlock.height) / 2;
                    return 4;
                }})();
                const outlineRect = g.rect(textBlock.width, textBlock.height)
                    .move(textBlock.x, textBlock.y)
                    .fill('transparent')
                    .stroke({{ color: textBlock.border_color, width: textBlock.border_width || 1 }})
                    .radius(outlineRadius)
                    .attr({{ 'data-role': 'text-outline', 'data-original-stroke': textBlock.border_color || 'black', 'data-original-stroke-width': textBlock.border_width || 1 }});
                const rectRadius = outlineRadius;
                const rect = g.rect(textBlock.width, textBlock.height).move(textBlock.x, textBlock.y).fill(textBlock.bgcolor).stroke({{ color: 'transparent', width: 0 }}).radius(rectRadius);
                rect.attr({{
                    'data-role': 'text-rect',
                    'data-original-stroke': textBlock.border_color || 'black',
                    'data-original-stroke-width': textBlock.border_width || 1
                }});
                const foreign = g.foreignObject(textBlock.width, textBlock.height).move(textBlock.x, textBlock.y).attr({{ 'data-role': 'text-fo' }});
                const editor = document.createElement('div');
                editor.className = 'mk-text-editor';
                editor.setAttribute('data-text-id', id);
                editor.style.background = isShape ? 'transparent' : textBlock.bgcolor;
                editor.innerHTML = isShape ? '' : (htmlContent || '');
                applyStyleToEditor(editor, style);
                editor.contentEditable = false;
                if (isTransparentShape) {{
                    hitRect.node.style.pointerEvents = 'none';
                    rect.node.style.pointerEvents = 'none';
                    foreign.node.style.pointerEvents = 'none';
                }}
                if (isShape) {{
                    editor.style.display = 'none';
                    foreign.node.style.display = 'none';
                    // Render shape label directly on SVG if provided
                    if (textBlock.label) {{
                        const alignMap = {{ left: 'start', center: 'middle', right: 'end' }};
                        const verticalMap = {{ top: 'top', center: 'middle', bottom: 'bottom' }};
                        const styleAlign = (textBlock.text_style && textBlock.text_style.align) ? textBlock.text_style.align.toLowerCase() : 'center';
                        const styleVertical = (textBlock.text_style && textBlock.text_style.vertical) ? textBlock.text_style.vertical.toLowerCase() : 'center';
                        const anchor = alignMap[styleAlign] || 'middle';
                        let tx = textBlock.x + textBlock.width / 2;
                        if (anchor === 'start') tx = textBlock.x + 4;
                        if (anchor === 'end') tx = textBlock.x + textBlock.width - 4;
                        const fontSize = textBlock.text_style?.fontSize || 12;
                        let ty = textBlock.y + textBlock.height / 2;
                        if (styleVertical === 'top') ty = textBlock.y + fontSize / 2 + 2;
                        if (styleVertical === 'bottom') ty = textBlock.y + textBlock.height - fontSize / 2 - 2;
                        const lbl = g.text(textBlock.label)
                            .font({{
                                size: fontSize,
                                family: textBlock.text_style?.fontFamily || baseTextStyle().fontFamily,
                                anchor,
                                weight: textBlock.text_style?.bold ? 'bold' : 'normal',
                                style: textBlock.text_style?.italic ? 'italic' : 'normal',
                                leading: '1.2em'
                            }})
                            .fill(textBlock.fgcolor || '#000000');
                        lbl.center(tx, ty);
                        lbl.attr({{ 'data-role': 'shape-label', 'pointer-events': 'all' }});
                        // Allow independent dragging of the shape label
                        makeDraggable(lbl, 'shape-label', `${{id}}_shape_label`, id);
                    }}
                }}
                textBlock.label = editor.textContent || textBlock.label || '';
                foreign.node.appendChild(editor);
                const focusEditor = (selectAll = false) => {{
                    selectElement(g, 'text-box', id);
                    enterTextEditMode(id, {{ selectAll }});
                }};
                if (!isShape) {{
                    editor.addEventListener('mousedown', evt => {{
                        selectTextBox(evt);
                    }});
                    editor.addEventListener('dblclick', evt => {{
                        evt.stopPropagation();
                        focusEditor(true);
                    }});
                    editor.addEventListener('contextmenu', evt => {{
                        selectTextBox(evt);
                        const sel = window.getSelection();
                        const hasSelection = selectionInsideEditor(editor) && sel && sel.rangeCount > 0 && !sel.isCollapsed;
                        if (selectionMap.size > 1) {{
                            evt.preventDefault();
                            evt.stopPropagation();
                            showGroupingMenu(evt, {{ mode: 'create' }});
                            return;
                        }}
                        if (hasSelection) {{
                            showTextFormatMenu(evt, id, editor);
                        }} else {{
                            showTextBoxMenu(evt, id, editor, isShape);
                        }}
                    }});
                    g.node.addEventListener('dblclick', () => focusEditor(true));
                    rect.node.addEventListener('dblclick', () => focusEditor(true));
                }}
                const selectTextBox = (evt = null) => {{
                    const additive = !!(evt && ((evt.ctrlKey || evt.metaKey) || evt.button === 2));
                    const toggle = !!(evt && (evt.ctrlKey || evt.metaKey));
                    selectElement(g, 'text-box', id, null, {{ additive, toggle }});
                }};
                const showHandles = () => {{
                    ensureTextHandles(id, textBlock);
                    toggleTextHandles(id, true);
                }};
                rect.on('mousedown', e => {{
                    if (isTextEditing && activeTextEditorId === id) return;
                    selectTextBox(e);
                    showHandles();
                }});
                outlineRect.on('mousedown', e => {{
                    if (isTextEditing && activeTextEditorId === id) return;
                    selectTextBox(e);
                    showHandles();
                }});
                foreign.on('mousedown', e => {{
                    if (isTextEditing && activeTextEditorId === id) return;
                    selectTextBox(e);
                    showHandles();
                }});
                hitRect.on('mousedown', e => {{
                    if (isTextEditing && activeTextEditorId === id) return;
                    selectTextBox(e);
                    showHandles();
                }});
                const maybeShowGroupMenu = (evt) => {{
                    const key = selectionKey('text-box', id);
                    if (selectionMap.has(key) && selectionMap.size > 1) {{
                        evt.preventDefault();
                        evt.stopPropagation();
                        showGroupingMenu(evt, {{ mode: 'create' }});
                        return true;
                    }}
                    return false;
                }};
                rect.node.addEventListener('contextmenu', e => {{
                    selectTextBox(e);
                    if (!maybeShowGroupMenu(e)) {{
                        showTextBoxMenu(e, id, editor, isShape);
                    }}
                }});
                outlineRect.node.addEventListener('contextmenu', e => {{
                    selectTextBox(e);
                    if (!maybeShowGroupMenu(e)) {{
                        showTextBoxMenu(e, id, editor, isShape);
                    }}
                }});
                hitRect.node.addEventListener('contextmenu', e => {{
                    selectTextBox(e);
                    if (!maybeShowGroupMenu(e)) {{
                        showTextBoxMenu(e, id, editor, isShape);
                    }}
                }});
                if (!isTransparentShape && !(isShape && textBlock.is_background)) {{
                    makeDraggable(hitRect, 'text-box', id);
                    makeDraggable(rect, 'text-box', id);
                    makeDraggable(foreign, 'text-box', id);
                }}
                if (!(isShape && textBlock.is_background)) {{
                    makeDraggable(g, 'text-box', id);
                    makeDraggable(outlineRect, 'text-box', id);
                }} else {{
                    hitRect.node.style.pointerEvents = 'none';
                    rect.node.style.pointerEvents = 'none';
                    foreign.node.style.pointerEvents = 'none';
                    outlineRect.node.style.pointerEvents = 'none';
                }}
                ensureTextHandles(id, textBlock);
                toggleTextHandles(id, false);
                if (options.startEditing) {{
                    setTimeout(() => focusEditor(true), 0);
                }}
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
            const createFigureKey = (orientation = 'vertical', svgX = 0, svgY = 0) => {{
                const id = `figure_key_${{++figureKeyCounter}}`;
                const isHorizontal = orientation === 'horizontal';
                const width = isHorizontal ? 160 : 26;
                const height = isHorizontal ? 26 : 160;
                const posColor = toRgbString(gradientConfig.posColor);
                const negColor = toRgbString(gradientConfig.negColor);
                const white = 'rgb(255,255,255)';
                const grad = draw.gradient('linear', add => {{
                    if (isHorizontal) {{
                        add.stop(0, negColor);
                        add.stop(0.5, white);
                        add.stop(1, posColor);
                        add.from(0, 0).to(1, 0);
                    }} else {{
                        add.stop(0, posColor);
                        add.stop(0.5, white);
                        add.stop(1, negColor);
                        add.from(0, 0).to(0, 1);
                    }}
                }});
                const g = figureKeyGroup.group().attr({{ 'data-id': id, 'data-type': 'figure-key' }});
                const rect = g.rect(width, height).move(svgX, svgY).fill(grad).stroke({{ color: 'black', width: 1 }});
                rect.attr({{
                    'data-type': 'compound-shape',
                    'data-original-stroke': 'black',
                    'data-original-stroke-width': 1
                }});
                const fontSize = 11;
                if (isHorizontal) {{
                    const labelY = svgY + height + fontSize + 2;
                    const negLabel = g.text(String(gradientConfig.maxNeg)).font({{ size: fontSize, family: 'Arial', anchor: 'start' }}).fill('#000');
                    negLabel.attr({{ 'data-type': 'compound-label', 'data-original-fill': '#000' }});
                    negLabel.move(svgX, labelY);
                    const zeroLabel = g.text('0').font({{ size: fontSize, family: 'Arial', anchor: 'middle' }}).fill('#000');
                    zeroLabel.attr({{ 'data-type': 'compound-label', 'data-original-fill': '#000' }});
                    zeroLabel.center(svgX + width / 2, labelY + fontSize / 2);
                    const posLabel = g.text(String(gradientConfig.maxPos)).font({{ size: fontSize, family: 'Arial', anchor: 'end' }}).fill('#000');
                    posLabel.attr({{ 'data-type': 'compound-label', 'data-original-fill': '#000' }});
                    posLabel.move(svgX + width - 2, labelY);
                }} else {{
                    const labelX = svgX + width + 6;
                    const posLabel = g.text(String(gradientConfig.maxPos)).font({{ size: fontSize, family: 'Arial', anchor: 'start' }}).fill('#000');
                    posLabel.attr({{ 'data-type': 'compound-label', 'data-original-fill': '#000' }});
                    posLabel.move(labelX, svgY - 2);
                    const zeroLabel = g.text('0').font({{ size: fontSize, family: 'Arial', anchor: 'start' }}).fill('#000');
                    zeroLabel.attr({{ 'data-type': 'compound-label', 'data-original-fill': '#000' }});
                    zeroLabel.move(labelX, svgY + height / 2 - fontSize);
                    const negLabel = g.text(String(gradientConfig.maxNeg)).font({{ size: fontSize, family: 'Arial', anchor: 'start' }}).fill('#000');
                    negLabel.attr({{ 'data-type': 'compound-label', 'data-original-fill': '#000' }});
                    negLabel.move(labelX, svgY + height - fontSize - 2);
                }}
                makeDraggable(g, 'figure-key', id);
                elementGroups[id] = [{{ element: g, offsetX: 0, offsetY: 0 }}];
                return g;
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
                const newLabel = debugMode ? String(protboxId) : (newProtein.label || pb.backup_label || 'Unknown');
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
                    const shapeId = ptmElementId(protboxId, newUniprot, ptm_key, 'shape');
                    shapeObj.attr({{ 'data-id': shapeId, 'data-type': 'ptm-shape' }});
                    elementGroups[protboxId].push({{ element: shapeObj, offsetX: adjusted_ptmX - currentX, offsetY: adjusted_ptmY - currentY, jsonX: adjusted_ptmX, jsonY: adjusted_ptmY }});
                    makeDraggable(shapeObj, 'ptm-shape', shapeId, protboxId);
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
                    Shiny?.setInputValue('element_moved', {{ type: 'ptm-shape', id: shapeId, x: shape_center_x, y: shape_center_y, protbox_id: protboxId, ptm_position: ptm_pos || null }}, {{ priority: 'event' }});
                    if (ptm.label) {{
                        const ld = labelDefaults[ptm_pos] || [];
                        const fallbackLabelX = toFiniteNumber(ptm.label_x);
                        const fallbackLabelY = toFiniteNumber(ptm.label_y);
                        const baseLabelX = fallbackLabelX !== null ? fallbackLabelX : resolvedShapeX + (typeof ld[0] === 'number' ? ld[0] : 0);
                        const baseLabelY = fallbackLabelY !== null ? fallbackLabelY : resolvedShapeY + (typeof ld[1] === 'number' ? ld[1] : 0);
                        const resolvedLabelX = resolveCoordinate(override?.label_x, baseLabelX) ?? baseLabelX;
                        const resolvedLabelY = resolveCoordinate(override?.label_y, baseLabelY) ?? baseLabelY;
                        const adjusted_labelX = resolvedLabelX + deltaX;
                        const adjusted_labelY = resolvedLabelY + deltaY;
                        const labelColor = Array.isArray(ptm.label_color) && ptm.label_color.length === 3 ? `rgb(${{ptm.label_color.join(',')}})` : 'black';
                        const labelCenter = (override?.label_centering || ptm.label_centering || (ld[2] || 'center')).toLowerCase();
                        const textAnchor = anchorMap[labelCenter] || 'middle';
                        const labelId = ptmElementId(protboxId, newUniprot, ptm_key, 'label');
                        const labelText = protboxGroup.text(ptm.label).move(adjusted_labelX, adjusted_labelY + textOffsetY).font({{ size: settings.ptm_label_size || 10, family: settings.ptm_label_font || 'Arial', anchor: textAnchor, leading: '1.2em' }}).fill(labelColor).attr({{ 'data-id': labelId, 'data-type': 'ptm-label', 'data-protbox-id': protboxId }});
                        elementGroups[protboxId].push({{ element: labelText, offsetX: adjusted_labelX - currentX, offsetY: (adjusted_labelY + textOffsetY) - currentY, jsonX: adjusted_labelX, jsonY: adjusted_labelY }});
                        makeDraggable(labelText, 'ptm-label', labelId, protboxId);
                        bindPtmContextMenu(labelText, protboxId, ptmMeta);
                        Shiny?.setInputValue('element_moved', {{ type: 'ptm-label', id: labelId, x: adjusted_labelX, y: adjusted_labelY, protbox_id: protboxId }}, {{ priority: 'event' }});
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
                        const symbolId = ptmElementId(protboxId, newUniprot, ptm_key, 'symbol');
                        const symbolText = protboxGroup.text(ptm.symbol).move(adjusted_symbolX, adjusted_symbolY + textOffsetY).font({{ size: ptm.symbol_size || settings.ptm_label_size || 10, family: ptm.symbol_font || settings.ptm_label_font || 'Arial', anchor: 'middle', leading: '1.2em' }}).fill(symbolColor).attr({{ 'data-id': symbolId, 'data-type': 'ptm-symbol', 'pointer-events': 'none', 'data-protbox-id': protboxId }});
                        elementGroups[protboxId].push({{ element: symbolText, offsetX: adjusted_symbolX - currentX, offsetY: (adjusted_symbolY + textOffsetY) - currentY, jsonX: adjusted_symbolX, jsonY: adjusted_symbolY }});
                        Shiny?.setInputValue('element_moved', {{ type: 'ptm-symbol', id: symbolId, x: adjusted_symbolX, y: adjusted_symbolY, protbox_id: protboxId }}, {{ priority: 'event' }});
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
              const shapes = textBlocks.filter(tb => tb && tb.shape_type);
              const nonShapes = textBlocks.filter(tb => !tb || !tb.shape_type);
              shapes.sort((a, b) => {{
                const areaA = (Number(a.width) || 0) * (Number(a.height) || 0);
                const areaB = (Number(b.width) || 0) * (Number(b.height) || 0);
                return areaB - areaA; // larger first => lower z
              }});
              const orderedTextBlocks = [...shapes, ...nonShapes];
              orderedTextBlocks.forEach((tb, idx) => {{
                try {{
                  renderTextBox(tb, idx);
                }} catch (err) {{
                  console.error('renderTextBox failed', err, tb);
                  const dbg = document.getElementById('debug_json');
                  if (dbg) dbg.textContent = 'Client error rendering text box: ' + (err?.message || err);
                }}
              }});
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
            initializeGroupState();
            const drawArrows = (arrowList, options = {{}}) => {{
                const force = options.force === true;
                if (!showArrows && !force) return;
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
                        const baseLineType = lineType.startsWith('dashed_') ? lineType.replace(/^dashed_/, '') : lineType;
                        let strokeOpts = {{ color: 'black', width: 1 }};
                        if (lineType.startsWith('dashed') || arrow.dashed) {{ strokeOpts.dasharray = '5,5'; }}
                        const hitboxWidth = 22; // make hit area larger (approx +6px padding each side)
                        const ctrlPoints = Array.isArray(arrow.control_points) ? arrow.control_points : [];
                        let hitbox = null;
                        let line = null;
                        if (ctrlPoints.length > 0) {{
                            const cp = ctrlPoints[0];
                            const pathStr = `M ${{x1}} ${{y1}} Q ${{cp.x}} ${{cp.y}} ${{x2}} ${{y2}}`;
                            hitbox = arrowGroup.path(pathStr).fill('none').stroke({{ color: 'transparent', width: hitboxWidth }}).attr({{ 'data-id': arrowId + '_hit', 'data-type': 'arrow-hitbox', 'data-arrow-id': arrowId }});
                            line = arrowGroup.path(pathStr).fill('none').stroke(strokeOpts).attr({{ 'data-id': arrowId, 'data-type': 'arrow', 'data-line-type': lineType, 'pointer-events': 'none' }});
                        }} else {{
                            hitbox = arrowGroup.line(x1, y1, x2, y2).stroke({{ color: 'transparent', width: hitboxWidth }}).attr({{ 'data-id': arrowId + '_hit', 'data-type': 'arrow-hitbox', 'data-arrow-id': arrowId }});
                            line = arrowGroup.line(x1, y1, x2, y2).stroke(strokeOpts).attr({{ 'data-id': arrowId, 'data-type': 'arrow', 'data-line-type': lineType, 'pointer-events': 'none' }});
                        }}
                        const dx = x2 - x1, dy = y2 - y1, angle = Math.atan2(dy, dx), arrowSize = 5;
                        let arrowHead = null;
                        if (baseLineType === 'arrow') {{
                            arrowHead = arrowGroup.polygon([[x2, y2],[x2 - arrowSize * Math.cos(angle + Math.PI / 6), y2 - arrowSize * Math.sin(angle + Math.PI / 6)],[x2 - arrowSize * Math.cos(angle - Math.PI / 6), y2 - arrowSize * Math.sin(angle - Math.PI / 6)]]).fill('black').stroke({{ color: 'black', width: 1 }}).attr({{ 'data-id': arrowId + '_head', 'data-type': 'arrow-head' }});
                        }} else if (baseLineType === 'inhibition') {{
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
                        const openArrowMenu = (evt) => {{
                            const additive = evt && (evt.ctrlKey || evt.metaKey);
                            selectElement(hitbox, 'arrow', arrowId, null, {{ additive, toggle: additive }});
                            if (shouldShowGroupMenu('arrow', arrowId)) {{
                                showGroupingMenu(evt, {{ mode: 'create' }});
                                return;
                            }}
                            showArrowMenu(evt, arrowId);
                        }};
                        hitbox.node.addEventListener('contextmenu', openArrowMenu);
                        if (line && line.node) line.node.addEventListener('contextmenu', openArrowMenu);
                        if (arrowHead && arrowHead.node) arrowHead.node.addEventListener('contextmenu', openArrowMenu);
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
                const meta = resolveGroupingMeta(e.target);
                if (meta) {{
                    let metaKey = null;
                    if (meta.type === 'group') {{
                        metaKey = selectionKey('group', meta.id);
                    }} else if (meta.type === 'prot-box') {{
                        metaKey = selectionKey('prot-box', meta.id);
                    }} else if (meta.protboxId) {{
                        metaKey = selectionKey('prot-box', meta.protboxId);
                    }} else if (meta.type) {{
                        metaKey = selectionKey(meta.type, meta.id);
                    }}
                    const isSelected = metaKey ? selectionMap.has(metaKey) : false;
                    if (meta.type === 'group' && isSelected) {{
                        e.preventDefault();
                        e.stopPropagation();
                        showGroupingMenu(e, {{ mode: 'group', groupId: meta.id }});
                        return;
                    }}
                    if (selectionMap.size > 1 && isSelected) {{
                        e.preventDefault();
                        e.stopPropagation();
                        const opt = {{ mode: 'create' }};
                        if (meta.type === 'prot-box') opt.targetProtboxId = meta.id;
                        showGroupingMenu(e, opt);
                        return;
                    }}
                }}
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
                menu.style.padding = '2px 0';
                menu.style.fontSize = '12px';
                menu.style.zIndex = '1000';
                const ul = document.createElement('ul');
                ul.style.listStyle = 'none';
                ul.style.margin = '0';
                ul.style.padding = '0';
                const liAddInteraction = document.createElement('li');
                liAddInteraction.textContent = 'Add Arrow';
                liAddInteraction.style.padding = '4px 8px';
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
                const liAddText = document.createElement('li');
                liAddText.textContent = 'Add Text Box';
                liAddText.style.padding = '4px 8px';
                liAddText.style.cursor = 'pointer';
                liAddText.addEventListener('click', () => {{
                    addNewTextBox(svgX, svgY);
                    menu.remove();
                }});
                const createInteractionItem = (iconType, type, options = {{}}) => {{
                    const dashed = options.dashed === true;
                    const label = options.label || (type.charAt(0).toUpperCase() + type.slice(1));
                    const interactionLi = document.createElement('li');
                    interactionLi.style.padding = '4px 8px';
                    interactionLi.style.cursor = 'pointer';
                    interactionLi.style.display = 'flex';
                    interactionLi.style.alignItems = 'center';
                    interactionLi.title = label;
                    const svgNS = 'http://www.w3.org/2000/svg';
                    const icon = document.createElementNS(svgNS, 'svg');
                    icon.setAttribute('width', '56');
                    icon.setAttribute('height', '18');
                    icon.setAttribute('viewBox', '0 0 56 18');
                    icon.style.display = 'block';
                    const line = document.createElementNS(svgNS, 'line');
                    line.setAttribute('x1', '6');
                    line.setAttribute('y1', '9');
                    line.setAttribute('x2', '50');
                    line.setAttribute('y2', '9');
                    line.setAttribute('stroke', '#222');
                    line.setAttribute('stroke-width', '2');
                    line.setAttribute('stroke-linecap', 'round');
                    if (dashed) line.setAttribute('stroke-dasharray', '5,4');
                    icon.appendChild(line);
                    if (iconType === 'arrow') {{
                        const poly = document.createElementNS(svgNS, 'polygon');
                        poly.setAttribute('points', '44,4 52,9 44,14');
                        poly.setAttribute('fill', '#222');
                        icon.appendChild(poly);
                    }} else if (iconType === 'inhibition') {{
                        const bar = document.createElementNS(svgNS, 'line');
                        bar.setAttribute('x1', '48');
                        bar.setAttribute('y1', '5');
                        bar.setAttribute('x2', '48');
                        bar.setAttribute('y2', '13');
                        bar.setAttribute('stroke', '#222');
                        bar.setAttribute('stroke-width', '2.5');
                        if (dashed) bar.setAttribute('stroke-dasharray', '5,4');
                        icon.appendChild(bar);
                    }}
                    interactionLi.appendChild(icon);
                    interactionLi.addEventListener('click', () => {{
                        addNewArrow(type, svgX, svgY);
                        menu.remove();
                    }});
                    interactionSubmenu.appendChild(interactionLi);
                }};
                createInteractionItem('arrow', 'arrow');
                createInteractionItem('inhibition', 'inhibition');
                createInteractionItem('line', 'line');
                createInteractionItem('arrow', 'dashed_arrow', {{ dashed: true, label: 'Dashed Arrow' }});
                createInteractionItem('inhibition', 'dashed_inhibition', {{ dashed: true, label: 'Dashed Inhibitor Line' }});
                createInteractionItem('line', 'dashed_line', {{ dashed: true, label: 'Dashed Line' }});
                liAddInteraction.appendChild(interactionSubmenu);
                liAddInteraction.addEventListener('mouseenter', () => {{ interactionSubmenu.style.display = 'block'; }});
                liAddInteraction.addEventListener('mouseleave', () => {{ interactionSubmenu.style.display = 'none'; }});
                const liAddShape = document.createElement('li');
                liAddShape.textContent = 'Add Shape';
                liAddShape.style.padding = '4px 8px';
                liAddShape.style.cursor = 'pointer';
                liAddShape.style.position = 'relative';
                const shapeSub = document.createElement('ul');
                shapeSub.style.position = 'absolute';
                shapeSub.style.left = '100%';
                shapeSub.style.top = '0';
                shapeSub.style.backgroundColor = 'white';
                shapeSub.style.border = '1px solid #ccc';
                shapeSub.style.padding = '4px';
                shapeSub.style.margin = '0';
                shapeSub.style.display = 'none';
                shapeSub.style.listStyle = 'none';
                const makeShapeItem = (type, radius) => {{
                    const li = document.createElement('li');
                    li.style.padding = '4px 6px';
                    li.style.cursor = 'pointer';
                    const icon = document.createElement('div');
                    icon.style.width = type === 'circle' ? '22px' : '24px';
                    icon.style.height = type === 'circle' ? '22px' : '18px';
                    icon.style.border = '2px solid #444';
                    icon.style.background = '#f5f5f5';
                    icon.style.borderRadius = radius;
                    li.appendChild(icon);
                    li.addEventListener('click', () => {{
                        addNewShape(type, svgX, svgY);
                        menu.remove();
                    }});
                    shapeSub.appendChild(li);
                }};
                makeShapeItem('square', '2px');
                makeShapeItem('rounded', '10px');
                makeShapeItem('circle', '50%');
                liAddShape.appendChild(shapeSub);
                liAddShape.addEventListener('mouseenter', () => {{ shapeSub.style.display = 'block'; }});
                liAddShape.addEventListener('mouseleave', () => {{ shapeSub.style.display = 'none'; }});
                ul.appendChild(liAddShape);
                const liAddProtbox = document.createElement('li');
                liAddProtbox.textContent = 'Add Protbox';
                liAddProtbox.style.padding = '4px 8px';
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
                const liFigureKey = document.createElement('li');
                liFigureKey.textContent = 'Add Legend';
                liFigureKey.style.padding = '4px 8px';
                liFigureKey.style.cursor = 'pointer';
                liFigureKey.style.position = 'relative';
                const figureKeySubmenu = document.createElement('ul');
                figureKeySubmenu.style.position = 'absolute';
                figureKeySubmenu.style.left = '100%';
                figureKeySubmenu.style.top = '0';
                figureKeySubmenu.style.backgroundColor = 'white';
                figureKeySubmenu.style.border = '1px solid #ccc';
                figureKeySubmenu.style.padding = '0';
                figureKeySubmenu.style.margin = '0';
                figureKeySubmenu.style.display = 'none';
                figureKeySubmenu.style.listStyle = 'none';
                const addKeyOption = (label, orientation) => {{
                    const li = document.createElement('li');
                    li.textContent = label;
                    li.style.padding = '4px 8px';
                    li.style.cursor = 'pointer';
                    li.addEventListener('click', () => {{
                        createFigureKey(orientation, svgX, svgY);
                        menu.remove();
                    }});
                    figureKeySubmenu.appendChild(li);
                }};
                addKeyOption('Horizontal', 'horizontal');
                addKeyOption('Vertical', 'vertical');
                liFigureKey.appendChild(figureKeySubmenu);
                liFigureKey.addEventListener('mouseenter', () => {{ figureKeySubmenu.style.display = 'block'; }});
                liFigureKey.addEventListener('mouseleave', () => {{ figureKeySubmenu.style.display = 'none'; }});
                // Append items in desired order
                ul.appendChild(liAddProtbox);
                ul.appendChild(liAddInteraction);
                ul.appendChild(liAddShape);
                ul.appendChild(liAddText);
                ul.appendChild(liFigureKey);
                menu.appendChild(ul);
                container.appendChild(menu);
                const removeMenu = () => {{
                    if (menu.parentNode) menu.remove();
                    document.removeEventListener('click', removeMenu);
                }};
                document.addEventListener('click', removeMenu);
            }});
            const addNewTextBox = (svgX, svgY) => {{
                const tb = {{
                    text_id: Date.now(),
                    x: svgX,
                    y: svgY,
                    width: 180,
                    height: 80,
                    bgcolor: 'transparent',
                    fgcolor: '#000000',
                    border_color: 'transparent',
                    label: '',
                    html: '',
                    text_style: baseTextStyle()
                }};
                ensureTextClientId(tb, textBlocks.length);
                textBlocks.push(tb);
                renderTextBox(tb, textBlocks.length - 1, {{ recordHistory: true, startEditing: true }});
                Shiny?.setInputValue('add_text_box', {{ text_block: tb }}, {{ priority: 'event' }});
            }};
            const addNewShape = (shapeType, svgX, svgY) => {{
                const tb = {{
                    text_id: Date.now(),
                    x: svgX,
                    y: svgY,
                    width: 120,
                    height: 80,
                    bgcolor: '#f5f5f5',
                    fgcolor: '#000000',
                    border_color: '#000000',
                    border_width: 1,
                    label: '',
                    html: '',
                    shape_type: shapeType,
                    text_style: baseTextStyle()
                }};
                ensureTextClientId(tb, textBlocks.length);
                textBlocks.push(tb);
                renderTextBox(tb, textBlocks.length - 1, {{ recordHistory: true, startEditing: false }});
                Shiny?.setInputValue('add_text_box', {{ text_block: tb }}, {{ priority: 'event' }});
            }};
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
        window.__mkSafeInit = window.__mkSafeInit || function() {{
            try {{
                initializeSvg();
            }} catch (err) {{
                console.error('m3: initializeSvg failed', err);
                const dbg = document.getElementById('debug_json');
                if (dbg) dbg.textContent = 'Client init error: ' + (err?.message || err);
                Shiny?.setInputValue('client_error', {{ label: 'initializeSvg', message: err?.message || String(err), stack: err?.stack || null, time: Date.now() }}, {{ priority: 'event' }});
            }}
        }};
        // Increase delay slightly to allow assets to load
        setTimeout(window.__mkSafeInit, 1500);
        </script>
    """
    return ui.div(
        ui.div(
            "",
            style="color: #666; fontSize: 12px; margin-bottom: 10px; padding: 5px; background-color: #f0f0f0; border-radius: 3px;"
        ),
        ui.HTML(f'''
            <div id="svgCanvas" style="width: {canvas_width_style}; height: {canvas_height_style}; background-color: white; position: relative; max-width: 100%;" tabindex="0">
                <div class="canvas-controls" style="display: none;">
                    <button id="reset-view" class="svg-control-btn" title="Reset View">
                        <i class="fa fa-search"></i>
                    </button>
                </div>
            </div>
        '''),
        ui.HTML(data_script + catalog_script + svg_js),
        class_="svg-container",
        **{"style": f"position: relative; width: {container_width_style}; height: {container_height_style}; overflow: auto; background-color: #f9f9f9; border: 1px solid #ccc; max-width: 100%;"}
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
                .mk-text-editor {{
                    width: 100%;
                    height: 100%;
                    box-sizing: border-box;
                    padding: 6px;
                    white-space: pre-wrap;
                    word-break: break-word;
                    outline: none;
                    overflow: hidden;
                    color: #000;
                    font-family: Arial, sans-serif;
                    font-size: 14px;
                }}
                .mk-text-editor:focus {{
                    outline: 2px solid #2c7be5;
                }}
                .text-context-menu {{
                    background: #fff;
                    border: 1px solid #ccc;
                    box-shadow: 0 4px 10px rgba(0,0,0,0.12);
                    padding: 6px;
                    border-radius: 4px;
                    display: grid;
                    gap: 4px;
                    min-width: 200px;
                    z-index: 1200;
                    font-size: 12px;
                }}
                .text-context-menu label {{
                    display: flex;
                    align-items: center;
                    gap: 6px;
                    font-size: 11px;
                    color: #444;
                }}
                .text-context-menu input[type="color"],
                .text-context-menu select {{
                    flex: 1;
                }}
                .text-context-menu .text-toggle-row {{
                    display: flex;
                    gap: 4px;
                }}
                .text-context-menu .text-toggle-row button {{
                    flex: 1;
                    padding: 4px 6px;
                    border: 1px solid #ccc;
                    background: #f8f8f8;
                    cursor: pointer;
                    border-radius: 4px;
                    font-size: 12px;
                }}
                .text-context-menu {{
                    overflow: visible;
                }}
                .text-context-menu .text-toggle-row button.active {{
                    background: #e6f0ff;
                    border-color: #2c7be5;
                    color: #174ea6;
                }}
                /* Compact context menus */
                .context-menu, .group-context-menu {{
                    font-family: Arial, sans-serif;
                    font-size: 12px;
                    padding: 2px 0 !important;
                }}
                .context-menu ul, .group-context-menu ul {{
                    padding: 0 !important;
                    margin: 0 !important;
                }}
                .context-menu li, .group-context-menu li,
                .context-menu div, .group-context-menu div {{
                    padding: 4px 8px !important;
                    line-height: 1.2;
                }}
                .text-menu-item {{
                    position: relative;
                    width: 100%;
                }}
                .text-menu-label {{
                    width: 100%;
                    text-align: left;
                    padding: 5px 8px;
                    border: 1px solid #ccc;
                    background: #f8f8f8;
                    cursor: pointer;
                    border-radius: 4px;
                    font-size: 11px;
                }}
                .text-menu-label:hover {{
                    background: #e8f0ff;
                    border-color: #bcd0f7;
                }}
                .text-submenu {{
                    position: absolute;
                    left: 100%;
                    top: 0;
                    background: #fff;
                    border: 1px solid #ccc;
                    box-shadow: 0 3px 8px rgba(0,0,0,0.12);
                    padding: 6px;
                    border-radius: 4px;
                    display: none;
                    z-index: 4001;
                    min-width: 140px;
                }}
                .text-menu-item.open > .text-submenu {{
                    display: block;
                }}
                .text-submenu button {{
                    display: block;
                    width: 100%;
                    padding: 4px 6px;
                    margin: 2px 0;
                    border: 1px solid #ddd;
                    background: #f8f8f8;
                    cursor: pointer;
                    font-size: 11px;
                    text-align: left;
                }}
                .text-submenu button:hover {{
                    background: #e8f0ff;
                    border-color: #bcd0f7;
                }}
                .text-submenu .color-row {{
                    display: flex;
                    align-items: center;
                    gap: 6px;
                    margin-bottom: 6px;
                }}
                .text-submenu .color-row input[type="text"] {{
                    width: 90px;
                    padding: 3px 4px;
                    border: 1px solid #ccc;
                    border-radius: 3px;
                    font-size: 11px;
                }}
                .text-resize-handle {{
                    fill: #fff;
                    stroke: #ff4d4f;
                    stroke-width: 1;
                    cursor: pointer;
                }}
            </style>
        """)
    ),
    ui.navset_tab(
        ui.nav_panel("Pathway", ui.div(), value="pathway"),
        ui.nav_panel("Import Data", ui.div("Import data page (coming soon)."), value="import"),
        ui.nav_panel("Custom", ui.div(), value="custom"),
        id="page_selector",
        selected="pathway",
    ),
    ui.panel_conditional(
        "input.page_selector != 'import'",
        ui.TagList(
            ui.h2(ui.output_text("page_title")),
            ui.output_ui("pathway_plot"),
            ui.output_text("debug_json"),
            ui.input_action_button("save_json", "Save Changes to JSON")
        )
    )
)


def server(input, output, session):
    json_data_reactive = reactive.Value(None)
    active_page = reactive.Value("pathway")
    pathway_json = reactive.Value(None)
    custom_json = reactive.Value(None)

    def _extract_catalog_info(data):
        if not data:
            return None
        catalog = data.get('_global_protein_catalog')
        if catalog:
            return catalog
        env_path = os.environ.get('GLOBAL_PROTEIN_CATALOG_PATH')
        if env_path:
            return {'path': env_path}
        return None

    def _set_active_json(data):
        json_data_reactive.set(data)
        page = active_page.get()
        if page == "custom":
            custom_json.set(data)
        elif page == "pathway":
            pathway_json.set(data)

    def _ensure_custom_state():
        data = custom_json.get()
        if data is None:
            catalog_info = _extract_catalog_info(pathway_json.get() or json_data_reactive.get())
            data = _build_blank_canvas(catalog_info)
            custom_json.set(data)
        _set_active_json(_clone_json(data))

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
        core = str(element_id)
        if "__" in core:
            core = core.split("__", 1)[1]
        parts = core.split('_')
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
        if pathway_json.get() is not None:
            return
        loaded = load_json_data()
        cloned = _clone_json(loaded) if loaded else None
        pathway_json.set(cloned)
        if active_page.get() == "pathway":
            _set_active_json(_clone_json(cloned) if cloned is not None else None)

    @reactive.Effect
    def sync_page_selection():
        selected = input.page_selector() or "pathway"
        if selected == active_page.get():
            return
        active_page.set(selected)
        if selected == "custom":
            _ensure_custom_state()
        elif selected == "pathway":
            data = pathway_json.get()
            if data is None:
                data = load_json_data()
                pathway_json.set(_clone_json(data) if data else None)
            _set_active_json(_clone_json(data))
        else:
            json_data_reactive.set(None)

    @output
    @render.text
    def page_title():
        page = active_page.get()
        if page == "custom":
            return "Custom Pathway Builder"
        if page == "import":
            return "Import Data"
        return "Pathway Visualization"

    @output
    @render.ui
    def pathway_plot():
        json_data = json_data_reactive.get()
        if json_data is None:
            page = active_page.get()
            if page == "import":
                return ui.div("Import data coming soon.")
            return ui.div("Error: Could not load JSON data.")
        return create_pathway_svg(json_data)

    @output
    @render.text
    def debug_json():
        json_data = json_data_reactive.get()
        if json_data is None:
            return "No JSON data loaded"
        page = active_page.get()
        return f"Page: {page} | JSON keys: {list(json_data.keys())}\nProtbox count: {len(json_data.get('protbox_data', []))}"

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
        elif element_type == 'text-box':
            blocks = json_data.get('text_data') or []
            target = None
            for tb in blocks:
                if tb is None:
                    continue
                if str(tb.get('text_id')) == str(element_id) or f"text_{tb.get('text_id')}" == str(element_id) or str(tb.get('_client_id')) == str(element_id) or str(tb.get('id')) == str(element_id):
                    target = tb
                    break
            if target is not None:
                if x is not None:
                    try:
                        target['x'] = float(x)
                    except (TypeError, ValueError):
                        pass
                if y is not None:
                    try:
                        target['y'] = float(y)
                    except (TypeError, ValueError):
                        pass
                if moved.get('width') is not None:
                    try:
                        target['width'] = float(moved.get('width'))
                    except (TypeError, ValueError):
                        pass
                if moved.get('height') is not None:
                    try:
                        target['height'] = float(moved.get('height'))
                    except (TypeError, ValueError):
                        pass
                if 'html' in moved:
                    target['html'] = moved.get('html')
                if 'label' in moved:
                    target['label'] = moved.get('label') or ''
                if 'text_style' in moved and isinstance(moved.get('text_style'), dict):
                    target['text_style'] = moved.get('text_style')
                if 'fgcolor' in moved:
                    target['fgcolor'] = moved.get('fgcolor')
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
        _set_active_json(json_data)

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
        _set_active_json(json_data)

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
        _set_active_json(json_data)

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
        _set_active_json(json_data)
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
        _set_active_json(json_data)

    @reactive.Effect
    @reactive.event(input.add_text_box)
    def add_text_box():
        added = input.add_text_box()
        if not added:
            return
        json_data = json_data_reactive.get()
        if not json_data:
            return
        text_block = added.get('text_block') if isinstance(added, dict) else None
        if text_block:
            json_data.setdefault('text_data', []).append(text_block)
            _set_active_json(json_data)

    @reactive.Effect
    @reactive.event(input.text_box_changed)
    def update_text_box():
        payload = input.text_box_changed()
        if not payload:
            return
        json_data = json_data_reactive.get()
        if not json_data:
            return
        text_blocks = json_data.setdefault('text_data', [])
        text_id = payload.get('text_id')
        dom_id = payload.get('id')
        target = None
        for tb in text_blocks:
            if not isinstance(tb, dict):
                continue
            if text_id is not None and str(tb.get('text_id')) == str(text_id):
                target = tb
                break
            if dom_id:
                generated = f"text_{tb.get('text_id')}" if tb.get('text_id') is not None else tb.get('_client_id')
                if str(generated) == str(dom_id) or str(tb.get('_client_id')) == str(dom_id) or str(tb.get('id')) == str(dom_id):
                    target = tb
                    break
        if not target:
            return
        modified = False
        for key in ('x', 'y', 'width', 'height'):
            if payload.get(key) is not None:
                try:
                    target[key] = float(payload.get(key))
                    modified = True
                except (TypeError, ValueError):
                    pass
        for key in ('label', 'html', 'bgcolor', 'fgcolor', 'border_color'):
            if key in payload:
                target[key] = payload.get(key)
                modified = True
        if isinstance(payload.get('text_style'), dict):
            target['text_style'] = payload['text_style']
            modified = True
        if modified:
            _set_active_json(json_data)

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
                    members = group.get('members')
                    if isinstance(members, list):
                        filtered_members = [m for m in members if not (isinstance(m, dict) and m.get('type') == 'prot-box' and str(m.get('id')) == protbox_id_str)]
                        if len(filtered_members) != len(members):
                            groups_changed = True
                        group['members'] = filtered_members
                    if filtered_ids or (group.get('members')):
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
            _set_active_json(json_data)

    @reactive.Effect
    @reactive.event(input.groups_changed)
    def sync_groups():
        payload = input.groups_changed()
        if not payload:
            return
        groups_payload = payload.get('groups')
        if not isinstance(groups_payload, list):
            return
        json_data = json_data_reactive.get()
        if not json_data:
            return
        cleaned = []
        for group in groups_payload:
            if not isinstance(group, dict):
                continue
            clean_entry = {}
            if group.get('group_id') is not None:
                clean_entry['group_id'] = str(group.get('group_id'))
            if 'show_box' in group:
                clean_entry['show_box'] = bool(group.get('show_box'))
            if group.get('box_padding') is not None:
                try:
                    clean_entry['box_padding'] = float(group.get('box_padding'))
                except (TypeError, ValueError):
                    pass
            if group.get('box_radius') is not None:
                try:
                    clean_entry['box_radius'] = float(group.get('box_radius'))
                except (TypeError, ValueError):
                    pass
            if group.get('anchor_member') is not None:
                clean_entry['anchor_member'] = str(group.get('anchor_member'))
            members_payload = group.get('members')
            if isinstance(members_payload, list):
                clean_members = []
                for member in members_payload:
                    if not isinstance(member, dict):
                        continue
                    mtype = member.get('type')
                    mid = member.get('id') or member.get('group_id') or member.get('protbox_id')
                    if mtype in ('group', 'prot-box', 'compound', 'text-box', 'figure-key', 'arrow') and mid is not None:
                        clean_members.append({'type': mtype, 'id': mid})
                clean_entry['members'] = clean_members
            prot_ids = group.get('protbox_ids')
            if isinstance(prot_ids, list):
                clean_entry['protbox_ids'] = [str(pid) for pid in prot_ids]
            cleaned.append(clean_entry)
        json_data['groups'] = cleaned
        _set_active_json(json_data)

    @reactive.Effect
    @reactive.event(input.save_json)
    def save_json():
        json_data = json_data_reactive.get()
        if json_data:
            save_json_data(json_data)



app = App(app_ui, server)

if __name__ == "__main__":
    app.run()

