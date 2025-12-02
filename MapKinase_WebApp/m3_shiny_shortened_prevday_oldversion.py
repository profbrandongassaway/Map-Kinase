import json
import os
from shiny import App, ui, render, reactive

json_file_path = r"C:\Users\clayt\OneDrive - Brigham Young University\Desktop\MapKinase\MapKinase_WebApp\hsa04010_pathway_data_20250907_210715.json"
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
    if not protbox_data:
        max_x, max_y = 800, 600
    else:
        max_x = max(pb.get('x', 0) + pb.get('width', 0) for pb in protbox_data if pb.get('x') is not None and pb.get('width') is not None) + 50
        max_y = max(pb.get('y', 0) + pb.get('height', 0) for pb in protbox_data if pb.get('y') is not None and pb.get('height') is not None) + 50
        max_x = max(max_x, 800)
        max_y = max(max_y, 600)
    # Embed the full JSON (including any 'kegg_bg_image' and preview settings) for the client
    data_script = f"""<script type="application/json" id="pathway-data">{json.dumps(json_data)}</script>"""
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
    var protboxMap = protboxMap || {{}};
    var attachments = attachments || {{}};
    var attachedBy = attachedBy || {{}};
    var protboxSnapPoints = protboxSnapPoints || {{}};
    var protboxHandleDists = protboxHandleDists || {{}};
    var labelDefaults = labelDefaults || {{'N1': [-17, -10, 'right'],'N2': [-6, -12, 'center'],'N3': [5, -10, 'left'],'S1': [-18, 4, 'right'],'S2': [0, 10, 'center'],'S3': [6, 2, 'left'],'W1': [-18, -6, 'right'],'W2': [-18, 2, 'right'],'E1': [5, -6, 'left'],'E2': [5, 2, 'left']}};
    var anchorMap = anchorMap || {{'left': 'start','center': 'middle','right': 'end'}};
    // Pick the first unoccupied PTM snap point for a protbox
    var ptmSnapRadius = 12;
    var currentSelected = currentSelected || {{}};
    function initializeSvg() {{
            const container = document.querySelector('.svg-container');
            if (!container || !document.getElementById('svgCanvas')) return;
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
            Object.assign(tooltip.style, {{position: 'absolute', backgroundColor: 'rgba(0,0,0,0.8)', color: 'white',padding: '5px 10px', borderRadius: '4px', fontSize: '12px', pointerEvents: 'none',display: 'none', maxWidth: '300px', whiteSpace: 'normal', wordWrap: 'break-word'}});
            if (!document.getElementById('mk-tooltip')) container.appendChild(tooltip);
            const elementGroups = {{}};
            const foregroundGroup = draw.group();
            const handleGroup = foregroundGroup.group();
            const arrowGroup = foregroundGroup.group();
            const protboxGroup = foregroundGroup.group();
            const compoundGroup = draw.group();  
            const textGroup = foregroundGroup.group();     
            // Keep the background fixed; shift all interactive layers by the inverse offset
            foregroundGroup.translate(fgOffsetX, fgOffsetY);

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
                const arrow = arrows.find(a => `arrow_${{arrows.indexOf(a)}}` === arrowId);
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
                const arrow = arrows.find(a => `arrow_${{arrows.indexOf(a)}}` === arrowId);
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
            const selectElement = (element, type, id, protboxId = null) => {{
                if (selectedElement) {{
                    if (selectedType === 'arrow') {{
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
                    if (selectedType === 'arrow') {{
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
                        const arrow = arrows.find(a => `arrow_${{arrows.indexOf(a)}}` === arrowId);
                        if (arrow) {{
                            if (arrow.protbox_id_1 && arrow.protbox_id_1_side) updateArrowPositions(arrow.protbox_id_1, arrow.protbox_id_1_side);
                            if (arrow.protbox_id_2 && arrow.protbox_id_2_side) updateArrowPositions(arrow.protbox_id_2, arrow.protbox_id_2_side);
                        }}
                    }}
                    if (selectedType === 'prot-box' && selectedId) {{
                        elementGroups[selectedId].forEach(el => (el.element.attr('data-type').startsWith('handle-') || el.element.attr('data-type') === 'ptm-snap-circle') && el.element.hide());
                    }}
                }}
                if (selectedProtboxId) {{
                    elementGroups[selectedProtboxId].forEach(el => el.element.attr('data-type') === 'ptm-snap-circle' && el.element.hide());
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
                if (e.target === draw.node) {{
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
                const arrow = arrows.find(a => `arrow_${{arrows.indexOf(a)}}` === arrowId);
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
                const newX = currentX + deltaX, newY = currentY + deltaY;
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
                            const arrow = arrows.find(a => `arrow_${{arrows.indexOf(a)}}` === arrowId);
                            const endSide = arrow?.protbox_id_2_side;
                            let barX1, barY1, barX2, barY2;
                            if (endSide === 'North' || endSide === 'South') {{barX1 = x2 - arrowSize;barY1 = y2;barX2 = x2 + arrowSize;barY2 = y2;}} 
                            else if (endSide === 'West' || endSide === 'East') {{barX1 = x2;barY1 = y2 - arrowSize;barX2 = x2;barY2 = y2 + arrowSize;}} 
                            else {{const perp = angle + Math.PI / 2;barX1 = x2 + Math.cos(perp) * arrowSize;barY1 = y2 + Math.sin(perp) * arrowSize;barX2 = x2 - Math.cos(perp) * arrowSize;barY2 = y2 - Math.sin(perp) * arrowSize;}}
                            arrowHead.plot(barX1, barY1, barX2, barY2);
                        }}
                    }}
                    const arrow = arrows.find(a => `arrow_${{arrows.indexOf(a)}}` === arrowId);
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
                                    const masterArrow = arrows.find(a => `arrow_${{arrows.indexOf(a)}}` === masterId);
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
                    const arrowId = selectedId.replace(/_(start|end)$/, '');
                    const arrowElement = draw.find(`[data-id="${{arrowId}}"]`)[0];
                    const arrowHitbox = draw.find(`[data-id="${{arrowId}}_hit"]`)[0];
                    const arrowHead = draw.find(`[data-id="${{arrowId}}_head"]`)[0];
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
                                const arrow = arrows.find(a => `arrow_${{arrows.indexOf(a)}}` === arrowId);
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
                    if (type === 'arrow') {{
                        const arrowId = id;
                        const arrow = arrows.find(a => `arrow_${{arrows.indexOf(a)}}` === arrowId);
                        if (arrow) {{
                            if (attachedBy[arrowId]) {{
                                ['start', 'end'].forEach(endType => {{
                                    if (attachedBy[arrowId][endType]) {{
                                        attachedBy[arrowId][endType].forEach(slaveInfo => {{
                                            const sArrow = arrows.find(a => `arrow_${{arrows.indexOf(a)}}` === slaveInfo.slave);
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
                        const arrow = arrows.find(a => `arrow_${{arrows.indexOf(a)}}` === arrowId);
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
                        const isCircle = (selectedType === 'ptm-shape' && selectedElement.type === 'circle');
                        currentX = parseFloat(isCircle ? selectedElement.cx() : selectedElement.x() || 0);
                        currentY = parseFloat(isCircle ? selectedElement.cy() : selectedElement.y() || 0);
                        Shiny?.setInputValue('element_moved', {{ type: selectedType, id: selectedId, x: currentX, y: currentY }}, {{ priority: 'event' }});
                        if (selectedType === 'ptm-shape') {{
                            updateHandleDists(selectedProtboxId);
                            if (selectedProtboxId in attachments) {{
                                Object.keys(attachments[selectedProtboxId]).forEach(side => updateArrowPositions(selectedProtboxId, side));
                            }}
                        }}
                        activeSnapKey = null;
                    }}
                }};
                document.addEventListener('mousemove', dragging);
                document.addEventListener('mouseup', endDragging);
            }};
            document.addEventListener('keydown', e => {{
                if (!selectedElement) return;
                let deltaX = 0, deltaY = 0, moveAmount = e.shiftKey ? 10 : 1;
                if (e.key === 'ArrowUp') {{ deltaY = -moveAmount; e.preventDefault(); }}
                else if (e.key === 'ArrowDown') {{ deltaY = moveAmount; e.preventDefault(); }}
                else if (e.key === 'ArrowLeft') {{ deltaX = -moveAmount; e.preventDefault(); }}
                else if (e.key === 'ArrowRight') {{ deltaX = moveAmount; e.preventDefault(); }}
                else if (e.key === 'Escape') {{ deselectElement(); e.preventDefault(); }}
                else return;
                if (deltaX !== 0 || deltaY !== 0) moveSelectedElement(deltaX, deltaY);
            }});
            draw.node.addEventListener('click', e => e.target === draw.node && deselectElement());
            draw.node.addEventListener('mousedown', startPanning);
            document.addEventListener('mousemove', pan);
            document.addEventListener('mouseup', endPanning);
            document.getElementById('svgCanvas').addEventListener('wheel', handleWheel, {{ passive: false }});
            document.getElementById('reset-view')?.addEventListener('click', resetView);
            const textOffsetY = -14.5;
            console.log('m3: about to render protBoxes, count =', Array.isArray(protBoxes) ? protBoxes.length : 0);
            try {{ var dbgEl2 = document.getElementById('debug_json'); if (dbgEl2) dbgEl2.textContent = 'Rendering protboxes: ' + (Array.isArray(protBoxes) ? protBoxes.length : 0); }} catch(e){{}}
            // Read box preview stretch factor (server provides under _box_preview.y_stretch)
            const boxPreview = data._box_preview || {{}};
            const boxYStretch = typeof boxPreview.y_stretch === 'number' ? boxPreview.y_stretch : (boxPreview.y_stretch ? Number(boxPreview.y_stretch) : 1.0);
            protBoxes.forEach((pb, index) => {{ 
                try {{
                    const id = pb.protbox_id || `unknown_${{index}}`, x = pb.x || 0, y = pb.y || 0, width = pb.width || 46, height = pb.height || 17;
                    // Apply Y-stretch to the box *position* only (do not alter height)
                    const yPos = Math.round(y * boxYStretch);
                    let selectedUniprot = pb.selected_uniprot || pb.proteins?.[0] || '';
                    currentSelected[id] = selectedUniprot;
                    let protein = proteinData[selectedUniprot] || {{}};
                    let label = protein.label || pb.backup_label || 'Unknown';
                    let fcColor = protein.fc_color_1 || [128, 128, 128];
                    let labelColor = protein.label_color || [0, 0, 0];
                    elementGroups[id] = [];
                    const fillRgb = Array.isArray(fcColor) && fcColor.length === 3 && fcColor.every(c => typeof c === 'number' && c >= 0 && c <= 255) ? `rgb(${{fcColor.join(',')}})` : 'rgb(128,128,128)';
                    const labelRgb = Array.isArray(labelColor) && labelColor.length === 3 && labelColor.every(c => typeof c === 'number' && c >= 0 && c <= 255) ? `rgb(${{labelColor.join(',')}})` : 'black';
                    const rect = protboxGroup.rect(width, height).move(x, yPos).fill(fillRgb).stroke({{ color: 'black', width: 1 }}).attr({{ 'data-id': id, 'data-type': 'prot-box', 'data-tooltip': (pb && pb.tooltip) || (protein && protein.tooltip) || '' }});
                    // Tooltip handlers for prot-box
                    rect.node.addEventListener('mouseenter', e => {{
                        const tip = rect.attr('data-tooltip') || '';
                        if (tip) {{ tooltip.textContent = tip; tooltip.style.display = 'block'; tooltip.style.left = (e.pageX + tooltipOffsetX) + 'px'; tooltip.style.top = (e.pageY + tooltipOffsetY) + 'px'; }}
                    }});
                    rect.node.addEventListener('mousemove', e => {{ if (tooltip.style.display === 'block') {{ tooltip.style.left = (e.pageX + tooltipOffsetX) + 'px'; tooltip.style.top = (e.pageY + tooltipOffsetY) + 'px'; }} }});
                    rect.node.addEventListener('mouseleave', () => {{ tooltip.style.display = 'none'; tooltip.textContent = ''; }});
                    const text = protboxGroup.text(label).move(x + width / 2, yPos + height / 2 + textOffsetY).font({{ size: settings.prot_label_size || 12, family: settings.prot_label_font || 'Arial', anchor: 'middle', leading: '1.2em' }}).fill(labelRgb).attr({{ 'data-id': id + '_label', 'data-type': 'prot-label', 'pointer-events': 'none' }});
                    elementGroups[id].push({{ element: text, offsetX: width / 2, offsetY: height / 2, jsonX: x + width / 2, jsonY: y + height / 2 }});
                    makeDraggable(rect, 'prot-box', id, id);
                    let northHas = false, southHas = false, westHas = false, eastHas = false;
                    if (selectedUniprot && selectedUniprot in proteinData) {{
                        for (const [ptm_key, ptm] of Object.entries(proteinData[selectedUniprot].PTMs || {{}})) {{
                            if (!ptm_key) continue;
                            if (ptm.shape_x === undefined || ptm.shape_y === undefined) continue;
                            const ptmX = ptm.shape_x || 0, ptmY = ptm.shape_y || 0, shape = (ptm.shape || 'circle').toLowerCase();
                            // Apply same Y-stretch transformation to PTM positions so they align with their protbox
                            const adjPtmY = Math.round(ptmY * boxYStretch);
                            const radius = settings.ptm_circle_radius || 5;
                            const fcColor = Array.isArray(ptm.fc_color_1) && ptm.fc_color_1.length === 3 && ptm.fc_color_1.every(c => typeof c === 'number' && c >= 0 && c <= 255) ? `rgb(${{ptm.fc_color_1.join(',')}})` : 'rgb(128,128,128)';
                                                        const shapeObj = shape === 'circle' ? protboxGroup.circle(radius * 2).cx(ptmX).cy(adjPtmY).fill(fcColor).stroke({{ color: 'black', width: 1 }}) : protboxGroup.rect(radius * 2, radius * 2).move(ptmX - radius, adjPtmY - radius).fill(fcColor).stroke({{ color: 'black', width: 1 }});
                                                        shapeObj.attr({{ 
                                                            'data-id': `${{selectedUniprot}}_${{ptm_key}}_shape`,
                                                            'data-type': 'ptm-shape',
                                                            'data-protbox-id': id,
                                                            'data-tooltip': (ptm && ptm.tooltip) || ''
                                                        }});
                                                        // Tooltip handlers for PTM shape
                                                        shapeObj.node.addEventListener('mouseenter', e => {{
                                                                const tip = shapeObj.attr('data-tooltip') || '';
                                                                if (tip) {{ tooltip.textContent = tip; tooltip.style.display = 'block'; tooltip.style.left = (e.pageX + tooltipOffsetX) + 'px'; tooltip.style.top = (e.pageY + tooltipOffsetY) + 'px'; }}
                                                        }});
                                                        shapeObj.node.addEventListener('mousemove', e => {{ if (tooltip.style.display === 'block') {{ tooltip.style.left = (e.pageX + tooltipOffsetX) + 'px'; tooltip.style.top = (e.pageY + tooltipOffsetY) + 'px'; }} }});
                                                        shapeObj.node.addEventListener('mouseleave', () => {{ tooltip.style.display = 'none'; tooltip.textContent = ''; }});
                            if (ptm.ptm_position) {{
                              shapeObj.attr({{ 'data-pos-key': ptm.ptm_position }});
                            }}

                            elementGroups[id].push({{ element: shapeObj, offsetX: ptmX - x, offsetY: ptmY - y, jsonX: ptmX, jsonY: ptmY }});
                            makeDraggable(shapeObj, 'ptm-shape', `${{selectedUniprot}}_${{ptm_key}}_shape`, id);
                            const ptm_pos = ptm.ptm_position || '';
                            if (ptm_pos.startsWith('N')) northHas = true;
                            if (ptm_pos.startsWith('S')) southHas = true;
                            if (ptm_pos.startsWith('W')) westHas = true;
                            if (ptm_pos.startsWith('E')) eastHas = true;
                            if (ptm.label) {{
                                const labelColor = Array.isArray(ptm.label_color) && ptm.label_color.length === 3 && ptm.label_color.every(c => typeof c === 'number' && c >= 0 && c <= 255) ? `rgb(${{ptm.label_color.join(',')}})` : 'black';
                                const labelX = ptm.label_x || ptmX, labelY = ptm.label_y || ptmY;
                                const adjLabelY = Math.round(labelY * boxYStretch);
                                const textAnchor = anchorMap[(ptm.label_centering || 'center').toLowerCase()] || 'middle';
                                const labelText = protboxGroup.text(ptm.label).move(labelX, adjLabelY + textOffsetY).font({{ size: settings.ptm_label_size || 10, family: settings.ptm_label_font || 'Arial', anchor: textAnchor, leading: '1.2em' }}).fill(labelColor).attr({{ 'data-id': `${{selectedUniprot}}_${{ptm_key}}_label`, 'data-type': 'ptm-label' }});
                                elementGroups[id].push({{ element: labelText, offsetX: labelX - x, offsetY: labelY + textOffsetY - y, jsonX: labelX, jsonY: labelY }});
                                makeDraggable(labelText, 'ptm-label', `${{selectedUniprot}}_${{ptm_key}}_label`, id);
                            }}
                            if (ptm.symbol) {{
                                const symbolX = ptm.symbol_x || ptmX, symbolY = ptm.symbol_y || ptmY;
                                const adjSymbolY = Math.round(symbolY * boxYStretch);
                                const symbolColor = Array.isArray(ptm.symbol_color) && ptm.symbol_color.length === 3 && ptm.symbol_color.every(c => typeof c === 'number' && c >= 0 && c <= 255) ? `rgb(${{ptm.symbol_color.join(',')}})` : 'black';
                                const symbolText = protboxGroup.text(ptm.symbol).move(symbolX, adjSymbolY + textOffsetY).font({{ size: ptm.symbol_size || settings.ptm_label_size || 10, family: ptm.symbol_font || settings.ptm_label_font || 'Arial', anchor: 'middle', leading: '1.2em' }}).fill(symbolColor).attr({{ 'data-id': `${{selectedUniprot}}_${{ptm_key}}_symbol`, 'data-type': 'ptm-symbol', 'pointer-events': 'none' }});
                                elementGroups[id].push({{ element: symbolText, offsetX: symbolX - x, offsetY: symbolY + textOffsetY - y, jsonX: symbolX, jsonY: symbolY }});
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
                        menu.style.left = `${{e.pageX + menuOffsetX}}px`;
                        menu.style.top = `${{e.pageY + menuOffsetY}}px`;
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
                        Object.entries(ptms).forEach(([ptm_key, ptm]) => {{
                            const subLi = document.createElement('li');
                            subLi.textContent = ptm.label || ptm_key;
                            subLi.style.padding = '5px 10px';
                            subLi.style.cursor = 'pointer';
                            const isSpawned = ptm.shape_x !== undefined && ptm.shape_y !== undefined;
                            if (isSpawned) {{
                                subLi.style.backgroundColor = 'grey';
                            }} else {{
                                subLi.style.backgroundColor = 'white';
                                subLi.addEventListener('click', () => {{
                                    const rect = draw.find(`[data-id="${{id}}"]`)[0];
                                    if (!rect) return;
                                    const pbX = rect.x();
                                    const pbY = rect.y();
                                    const pbWidth = rect.width();
                                    const pbHeight = rect.height();
                                    const spacing = settings.ptm_circle_spacing || 4;
                                    // Prefer an open snap point around this protbox
                                    let snapChoice = pickFreeSnap(id);
                                    let newX, newY, posKey = null;
                                    if (snapChoice) {{
                                        newX = snapChoice.x;
                                        newY = snapChoice.y;
                                        posKey = snapChoice.key;
                                    }} else {{
                                        // fallback: top-center
                                        newX = pbX + pbWidth / 2;
                                        newY = pbY - spacing;
                                    }}
                                    const radius = settings.ptm_circle_radius || 5;
                                    const shape = (ptm.shape || 'circle').toLowerCase();
                                    const fcColor = Array.isArray(ptm.fc_color_1) && ptm.fc_color_1.length === 3 ? `rgb(${{ptm.fc_color_1.join(',')}})` : 'rgb(128,128,128)';
                                    let shapeObj;
                                    if (shape === 'circle') {{
                                        shapeObj = protboxGroup.circle(radius * 2).cx(newX).cy(newY).fill(fcColor).stroke({{ color: 'black', width: 1 }});
                                    }} else {{
                                        shapeObj = protboxGroup.rect(radius * 2, radius * 2).move(newX - radius, newY - radius).fill(fcColor).stroke({{ color: 'black', width: 1 }});
                                    }}
                                    shapeObj.attr({{
                                      'data-id': `${{uniprot}}_${{ptm_key}}_shape`,
                                      'data-type': 'ptm-shape',
                                      'data-protbox-id': id,
                                      'data-tooltip': (ptm && ptm.tooltip) || '',
                                      ...(posKey ? {{ 'data-pos-key': posKey }} : {{}})
                                    }});
                                    // Tooltip handlers for spawned PTM
                                    shapeObj.node.addEventListener('mouseenter', e => {{
                                        const tip = shapeObj.attr('data-tooltip') || '';
                                        if (tip) {{ tooltip.textContent = tip; tooltip.style.display = 'block'; tooltip.style.left = (e.pageX + 10) + 'px'; tooltip.style.top = (e.pageY + 10) + 'px'; }}
                                    }});
                                    shapeObj.node.addEventListener('mousemove', e => {{ if (tooltip.style.display === 'block') {{ tooltip.style.left = (e.pageX + 10) + 'px'; tooltip.style.top = (e.pageY + 10) + 'px'; }} }});
                                    shapeObj.node.addEventListener('mouseleave', () => {{ tooltip.style.display = 'none'; tooltip.textContent = ''; }});
                                    elementGroups[id].push({{ element: shapeObj, offsetX: newX - pbX, offsetY: newY - pbY, jsonX: newX, jsonY: newY }});
                                    makeDraggable(shapeObj, 'ptm-shape', `${{uniprot}}_${{ptm_key}}_shape`, id);
                                    protboxHandleDists[id].North = 10;
                                    updateHandleDists(id);
                                    let labelX = ptm.label_x || newX;
                                    let labelY = ptm.label_y || newY;
                                    if (ptm.label) {{
                                        const labelColor = Array.isArray(ptm.label_color) && ptm.label_color.length === 3 ? `rgb(${{ptm.label_color.join(',')}})` : 'black';
                                        const textAnchor = anchorMap[(ptm.label_centering || 'center').toLowerCase()] || 'middle';
                                        const labelText = protboxGroup.text(ptm.label).move(labelX, labelY + textOffsetY).font({{ size: settings.ptm_label_size || 10, family: settings.ptm_label_font || 'Arial', anchor: textAnchor, leading: '1.2em' }}).fill(labelColor).attr({{ 'data-id': `${{uniprot}}_${{ptm_key}}_label`, 'data-type': 'ptm-label' }});
                                        elementGroups[id].push({{ element: labelText, offsetX: labelX - pbX, offsetY: (labelY + textOffsetY) - pbY, jsonX: labelX, jsonY: labelY }});
                                        makeDraggable(labelText, 'ptm-label', `${{uniprot}}_${{ptm_key}}_label`, id);
                                    }}
                                    let symbolX = ptm.symbol_x || newX;
                                    let symbolY = ptm.symbol_y || newY;
                                    if (ptm.symbol) {{
                                        const symbolColor = Array.isArray(ptm.symbol_color) && ptm.symbol_color.length === 3 ? `rgb(${{ptm.symbol_color.join(',')}})` : 'black';
                                        const symbolText = protboxGroup.text(ptm.symbol).move(symbolX, symbolY + textOffsetY).font({{ size: ptm.symbol_size || settings.ptm_label_size || 10, family: ptm.symbol_font || settings.ptm_label_font || 'Arial', anchor: 'middle', leading: '1.2em' }}).fill(symbolColor).attr({{ 'data-id': `${{uniprot}}_${{ptm_key}}_symbol`, 'data-type': 'ptm-symbol', 'pointer-events': 'none' }});
                                        elementGroups[id].push({{ element: symbolText, offsetX: symbolX - pbX, offsetY: (symbolY + textOffsetY) - pbY, jsonX: symbolX, jsonY: symbolY }});
                                    }}
                                    ptm.shape_x = newX;
                                    ptm.shape_y = newY;
                                    ptm.ptm_position = 'N2';
                                    ptm.label_x = labelX;
                                    ptm.label_y = labelY;
                                    ptm.symbol_x = symbolX;
                                    ptm.symbol_y = symbolY;
                                    Shiny?.setInputValue('ptm_spawned', {{ protbox_id: id, uniprot, ptm_key, shape_x: newX, shape_y: newY, ptm_position: posKey || null, ptm_position: 'N2', label_x: labelX, label_y: labelY, symbol_x: symbolX, symbol_y: symbolY }}, {{ priority: 'event' }});
                                    menu.remove();
                                }});
                            }}
                            ptmSubmenu.appendChild(subLi);
                        }});
                        liPTMs.appendChild(ptmSubmenu);
                        liPTMs.addEventListener('mouseenter', () => {{ ptmSubmenu.style.display = 'block'; }});
                        liPTMs.addEventListener('mouseleave', () => {{ ptmSubmenu.style.display = 'none'; }});
                        ul.appendChild(liPTMs);
                        const liUniprot = document.createElement('li');
                        liUniprot.textContent = 'UniProt';
                        liUniprot.style.padding = '5px 10px';
                        liUniprot.style.cursor = 'pointer';
                        liUniprot.addEventListener('click', () => {{
                            const uniprot = currentSelected[id];
                            if (uniprot) {{
                                window.open(`https://www.uniprot.org/uniprotkb/${{uniprot}}/entry`, '_blank');
                            }}
                            menu.remove();
                        }});
                        ul.appendChild(liUniprot);
                        menu.appendChild(ul);
                        container.appendChild(menu);
                        const removeMenu = () => {{
                            if (menu.parentNode) menu.remove();
                            document.removeEventListener('click', removeMenu);
                        }};
                        document.addEventListener('click', removeMenu);
                    }});
                }} catch (e) {{ }}
            }});
            const switchProtein = (protboxId, newUniprot) => {{
                const pb = protBoxes.find(p => p.protbox_id === protboxId);
                if (!pb) return;
                currentSelected[protboxId] = newUniprot;
                const rect = draw.find(`[data-id="${{protboxId}}"]`)[0];
                const text = draw.find(`[data-id="${{protboxId}}_label"]`)[0];
                const newProtein = proteinData[newUniprot] || {{}};
                const newLabel = newProtein.label || pb.backup_label || 'Unknown';
                const newFcColor = newProtein.fc_color_1 || [128, 128, 128];
                const newLabelColor = newProtein.label_color || [0, 0, 0];
                const fillRgb = Array.isArray(newFcColor) && newFcColor.length === 3 ? `rgb(${{newFcColor.join(',')}})` : 'rgb(128,128,128)';
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
                for (const [ptm_key, ptm] of Object.entries(newProtein.PTMs || {{}})) {{
                    if (!ptm_key) continue;
                    if (ptm.shape_x === undefined || ptm.shape_y === undefined) continue;
                    const adjusted_ptmX = (ptm.shape_x || 0) + deltaX;
                    const adjusted_ptmY = (ptm.shape_y || 0) + deltaY;
                    const shape = (ptm.shape || 'circle').toLowerCase();
                    const radius = settings.ptm_circle_radius || 5;
                    const fcColor = Array.isArray(ptm.fc_color_1) && ptm.fc_color_1.length === 3 ? `rgb(${{ptm.fc_color_1.join(',')}})` : 'rgb(128,128,128)';
                    const shapeObj = shape === 'circle' ? protboxGroup.circle(radius * 2).cx(adjusted_ptmX).cy(adjusted_ptmY).fill(fcColor).stroke({{ color: 'black', width: 1 }}) : protboxGroup.rect(radius * 2, radius * 2).move(adjusted_ptmX - radius, adjusted_ptmY - radius).fill(fcColor).stroke({{ color: 'black', width: 1 }});
                    shapeObj.attr({{ 'data-id': `${{newUniprot}}_${{ptm_key}}_shape`, 'data-type': 'ptm-shape' }});
                    elementGroups[protboxId].push({{ element: shapeObj, offsetX: adjusted_ptmX - currentX, offsetY: adjusted_ptmY - currentY, jsonX: adjusted_ptmX, jsonY: adjusted_ptmY }});
                    makeDraggable(shapeObj, 'ptm-shape', `${{newUniprot}}_${{ptm_key}}_shape`, protboxId);
                    const ptm_pos = ptm.ptm_position || '';
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
                    Shiny?.setInputValue('element_moved', {{ type: 'ptm-shape', id: `${{newUniprot}}_${{ptm_key}}_shape`, x: shape_center_x, y: shape_center_y }}, {{ priority: 'event' }});
                    if (ptm.label) {{
                        const adjusted_labelX = (ptm.label_x ? ptm.label_x + deltaX : adjusted_ptmX);
                        const adjusted_labelY = (ptm.label_y ? ptm.label_y + deltaY : adjusted_ptmY);
                        const labelColor = Array.isArray(ptm.label_color) && ptm.label_color.length === 3 ? `rgb(${{ptm.label_color.join(',')}})` : 'black';
                        const textAnchor = anchorMap[(ptm.label_centering || 'center').toLowerCase()] || 'middle';
                        const labelText = protboxGroup.text(ptm.label).move(adjusted_labelX, adjusted_labelY + textOffsetY).font({{ size: settings.ptm_label_size || 10, family: settings.ptm_label_font || 'Arial', anchor: textAnchor, leading: '1.2em' }}).fill(labelColor).attr({{ 'data-id': `${{newUniprot}}_${{ptm_key}}_label`, 'data-type': 'ptm-label' }});
                        elementGroups[protboxId].push({{ element: labelText, offsetX: adjusted_labelX - currentX, offsetY: (adjusted_labelY + textOffsetY) - currentY, jsonX: adjusted_labelX, jsonY: adjusted_labelY }});
                        makeDraggable(labelText, 'ptm-label', `${{newUniprot}}_${{ptm_key}}_label`, protboxId);
                        Shiny?.setInputValue('element_moved', {{ type: 'ptm-label', id: `${{newUniprot}}_${{ptm_key}}_label`, x: labelText.x(), y: labelText.y() - textOffsetY }}, {{ priority: 'event' }});
                    }}
                    if (ptm.symbol) {{
                        const adjusted_symbolX = (ptm.symbol_x ? ptm.symbol_x + deltaX : adjusted_ptmX);
                        const adjusted_symbolY = (ptm.symbol_y ? ptm.symbol_y + deltaY : adjusted_ptmY);
                        const symbolColor = Array.isArray(ptm.symbol_color) && ptm.symbol_color.length === 3 ? `rgb(${{ptm.symbol_color.join(',')}})` : 'black';
                        const symbolText = protboxGroup.text(ptm.symbol).move(adjusted_symbolX, adjusted_symbolY + textOffsetY).font({{ size: ptm.symbol_size || settings.ptm_label_size || 10, family: ptm.symbol_font || settings.ptm_label_font || 'Arial', anchor: 'middle', leading: '1.2em' }}).fill(symbolColor).attr({{ 'data-id': `${{newUniprot}}_${{ptm_key}}_symbol`, 'data-type': 'ptm-symbol', 'pointer-events': 'none' }});
                        elementGroups[protboxId].push({{ element: symbolText, offsetX: adjusted_symbolX - currentX, offsetY: (adjusted_symbolY + textOffsetY) - currentY, jsonX: adjusted_symbolX, jsonY: adjusted_symbolY }});
                        Shiny?.setInputValue('element_moved', {{ type: 'ptm-symbol', id: `${{newUniprot}}_${{ptm_key}}_symbol`, x: symbolText.x(), y: symbolText.y() - textOffsetY }}, {{ priority: 'event' }});
                    }}
                }}
                protboxHandleDists[protboxId] = {{North: northHas ? 10 : 5, South: southHas ? 10 : 5, West: westHas ? 10 : 5, East: eastHas ? 10 : 5}};
                updateHandleDists(protboxId);
                if (attachments[protboxId]) {{
                    Object.keys(attachments[protboxId]).forEach(side => updateArrowPositions(protboxId, side));
                }}
                Shiny?.setInputValue('protein_switched', {{ protbox_id: protboxId, uniprot: newUniprot }}, {{ priority: 'event' }});
            }};
            protBoxes.forEach(pb => {{
                const id = pb.protbox_id || `unknown_${{protBoxes.indexOf(pb)}}`;
                if (id) {{
                    protboxMap[id] = {{x: pb.x || 0,y: pb.y || 0,width: pb.width || 46,height: pb.height || 17}};
                }}
            }});
            
            // --- COMPOUNDS ---
            compounds.forEach((c, idx) => {{
              const id = c.compound_id ? 'compound_' + c.compound_id : 'compound_' + idx;
              const x = c.x || 0, y = c.y || 0;
              const w = c.width || 14, h = c.height || 14;
              const r = Math.max(w, h) / 2;
              const label = c.label || (c.kegg_compound || '').replace(/^cpd:/, '') || 'C?';
            
              const stroke = 'black';
              const fill = c.bgcolor || '#FFFFFF';
              const textColor = c.fgcolor || '#000000';
            
              const g = compoundGroup.group().attr({{ 'data-id': id, 'data-type': 'compound' }});
            
              // Circle centered at (x, y)
              g.circle(r * 2).cx(x).cy(y).fill(fill).stroke({{ color: stroke, width: 1 }});
            
              // Label centered horizontally, a few px below the circle
              const labelOffset = (settings.compound_label_offset || 6); // tweak in settings if you want
              const t = g.text(label)
                .font({{
                  size: (settings.compound_label_size || 10),
                  family: settings.compound_label_font || 'Arial',
                  anchor: 'middle',
                  leading: '1.2em'
                }})
                .fill(textColor);
            
              t.center(x, y + r + labelOffset);
              t.attr({{
                'data-id': id + '_label',
                'data-type': 'compound-label',
                'pointer-events': 'none'
              }});
            
              makeDraggable(g, 'compound', id);
            }});

            
            // --- TEXT BLOCKS / MAP LABELS ---
            textBlocks.forEach((tb, idx) => {{
              const id = tb.text_id ? 'text_' + tb.text_id : 'text_' + idx;
              const x = tb.x || 0, y = tb.y || 0;
              const w = tb.width || 60, h = tb.height || 20;
              const label = tb.label || '';
              const stroke = 'black';
              const fill = tb.bgcolor || '#FFFFFF';
              const textColor = tb.fgcolor || '#000000';
            
              const g = textGroup.group().attr({{ 'data-id': id, 'data-type': 'text-box' }});
              g.rect(w, h).move(x, y).fill(fill).stroke({{ color: stroke, width: 1 }}).radius(6);
            
              const cx = x + w / 2;
              const cy = y + h / 2;
            
              const t = g.text(label)
                .font({{
                  size: (settings.textbox_label_size || 11),
                  family: settings.textbox_label_font || 'Arial',
                  anchor: 'middle',
                  leading: '1.2em'
                }})
                .fill(textColor);
            
              // Center both horizontally and vertically
              t.center(cx, cy);
              t.attr({{
                'dominant-baseline': 'middle',     // vertical centering
                'data-id': id + '_label',
                'data-type': 'text-label',
                'pointer-events': 'none'
              }});
            
              makeDraggable(g, 'text-box', id);
            }});



            
            const attachedBy = {{}};
            arrows.forEach((arrow, index) => {{
                const arrowId = `arrow_${{index}}`;
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
                const arrowId = `arrow_${{index}}`;
                const id1 = arrow.protbox_id_1;
                const side1 = arrow.protbox_id_1_side;
                if (id1 && side1 && protboxMap[id1]) {{
                    attachments[id1] = attachments[id1] || {{}};
                    attachments[id1][side1] = attachments[id1][side1] || [];
                    attachments[id1][side1].push({{ type: 'start', arrow, arrowId }});
                }}
                const id2 = arrow.protbox_id_2;
                const side2 = arrow.protbox_id_2_side;
                if (id2 && side2 && protboxMap[id2]) {{
                    attachments[id2] = attachments[id2] || {{}};
                    attachments[id2][side2] = attachments[id2][side2] || [];
                    attachments[id2][side2].push({{ type: 'end', arrow, arrowId }});
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
                arrowList.forEach((arrow, index) => {{
                    try {{
                        const arrowId = `arrow_${{arrows.indexOf(arrow)}}`;
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
                const nonLines = arrows.filter(a => a.line !== 'line');
                const lines = arrows.filter(a => a.line === 'line');
                drawArrows(nonLines);
                drawArrows(lines);
            }}
            draw.node.addEventListener('contextmenu', e => {{
                if (e.target !== draw.node) return;
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
                menu.style.left = `${{e.pageX + menuOffsetX}}px`;
                menu.style.top = `${{e.pageY + menuOffsetY}}px`;
                menu.style.backgroundColor = 'white';
                menu.style.border = '1px solid #ccc';
                menu.style.padding = '5px';
                menu.style.zIndex = '1000';
                const ul = document.createElement('ul');
                ul.style.listStyle = 'none';
                ul.style.margin = '0';
                ul.style.padding = '0';
                const addArrowLi = document.createElement('li');
                addArrowLi.textContent = 'Add Arrow';
                addArrowLi.style.padding = '5px 10px';
                addArrowLi.style.cursor = 'pointer';
                addArrowLi.addEventListener('click', () => {{
                    addNewArrow('arrow', svgX, svgY);
                    menu.remove();
                }});
                ul.appendChild(addArrowLi);
                const addInhibitorLi = document.createElement('li');
                addInhibitorLi.textContent = 'Add Inhibitor Line';
                addInhibitorLi.style.padding = '5px 10px';
                addInhibitorLi.style.cursor = 'pointer';
                addInhibitorLi.addEventListener('click', () => {{
                    addNewArrow('inhibition', svgX, svgY);
                    menu.remove();
                }});
                ul.appendChild(addInhibitorLi);
                const addLineLi = document.createElement('li');
                addLineLi.textContent = 'Add Line';
                addLineLi.style.padding = '5px 10px';
                addLineLi.style.cursor = 'pointer';
                addLineLi.addEventListener('click', () => {{
                    addNewArrow('line', svgX, svgY);
                    menu.remove();
                }});
                ul.appendChild(addLineLi);
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
            </div>
            <button id="reset-view" title="Reset View" style="position: absolute; bottom: 10px; right: 10px; padding: 5px; font-size: 16px; cursor: pointer;">
                <i class="fa fa-search"></i>
            </button>
        '''),
        ui.HTML(data_script + svg_js),
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
                #reset-view {{ background-color: #fff; border: 1px solid #ccc; border-radius: 3px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
                #reset-view:hover {{ background-color: #f0f0f0; }}
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
        if element_type == 'prot-box':
            for pb in json_data['protbox_data']:
                if pb['protbox_id'] == element_id:
                    pb['x'] = float(x)
                    pb['y'] = float(y)
                    break
        elif element_type == 'ptm-shape':
            uniprot_id, ptm_key, _ = element_id.split('_', 2)
            if uniprot_id in json_data['protein_data']:
                ptm = json_data['protein_data'][uniprot_id]['PTMs'].get(ptm_key)
                if ptm:
                    ptm['shape_x'] = float(x)
                    ptm['shape_y'] = float(y)
        elif element_type == 'ptm-label':
            uniprot_id, ptm_key, _ = element_id.split('_', 2)
            if uniprot_id in json_data['protein_data']:
                ptm = json_data['protein_data'][uniprot_id]['PTMs'].get(ptm_key)
                if ptm:
                    ptm['label_x'] = float(x)
                    ptm['label_y'] = float(y)
                    if 'label_centering' in moved:
                        ptm['label_centering'] = moved['label_centering']
        elif element_type == 'ptm-symbol':
            uniprot_id, ptm_key, _ = element_id.split('_', 2)
            if uniprot_id in json_data['protein_data']:
                ptm = json_data['protein_data'][uniprot_id]['PTMs'].get(ptm_key)
                if ptm:
                    ptm['symbol_x'] = float(x)
                    ptm['symbol_y'] = float(y)
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
        if uniprot in json_data['protein_data']:
            ptms = json_data['protein_data'][uniprot]['PTMs']
            if ptm_key in ptms:
                ptm = ptms[ptm_key]
                ptm['shape_x'] = spawned['shape_x']
                ptm['shape_y'] = spawned['shape_y']
                ptm['ptm_position'] = spawned['ptm_position']
                ptm['label_x'] = spawned.get('label_x', ptm.get('label_x'))
                ptm['label_y'] = spawned.get('label_y', ptm.get('label_y'))
                ptm['symbol_x'] = spawned.get('symbol_x', ptm.get('symbol_x'))
                ptm['symbol_y'] = spawned.get('symbol_y', ptm.get('symbol_y'))
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
    @reactive.event(input.save_json)
    def save_json():
        json_data = json_data_reactive.get()
        if json_data:
            save_json_data(json_data)



app = App(app_ui, server)

if __name__ == "__main__":
    app.run()

