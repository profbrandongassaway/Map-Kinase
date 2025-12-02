<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Pathway Visualization</title>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.4/css/all.min.css">
    <script src="https://cdnjs.cloudflare.com/ajax/libs/svg.js/3.2.0/svg.min.js"></script>
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f9f9f9;
        }
        .svg-container {
            border: 1px solid #ccc;
            padding: 10px;
            background-color: white;
            position: relative;
            overflow: auto;
            width: 820px;
            height: 620px;
            display: none; /* Hidden until data loads */
        }
        #svgCanvas {
            width: 100%;
            height: 100%;
            outline: none;
        }
        #svgCanvas:focus {
            outline: 2px solid #007bff;
            outline-offset: -2px;
        }
        #reset-view {
            position: absolute;
            bottom: 10px;
            right: 10px;
            padding: 5px;
            font-size: 16px;
            cursor: pointer;
            background-color: #fff;
            border: 1px solid #ccc;
            border-radius: 3px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        #reset-view:hover {
            background-color: #f0f0f0;
        }
        .instructions {
            color: #666;
            font-size: 12px;
            margin-bottom: 10px;
            padding: 5px;
            background-color: #f0f0f0;
            border-radius: 3px;
        }
        #loading, #error {
            display: none;
            color: #d9534f;
            margin-top: 10px;
        }
        #loading.show, #error.show {
            display: block;
        }
    </style>
</head>
<body>
    <h2>Pathway Visualization</h2>
    <div class="instructions">
        Instructions: Click on protein boxes to select them (red outline). Use arrow keys to move selected elements. Hold Shift for larger movements (10px). Press Escape to deselect. Scroll mouse wheel to zoom in/out. Drag background to pan. Click reset button to reset view.
    </div>
    <div class="svg-container" id="svgCanvas" tabindex="0"></div>
    <button id="reset-view" title="Reset View"><i class="fa fa-search"></i></button>
    <div id="loading">Loading pathway data...</div>
    <div id="error"></div>

    <script>
        let svg;
        let draw;
        let zoomLevel = 1;
        const minZoom = 0.5;
        const maxZoom = 2;
        let viewBox = { x: 0, y: 0, width: 800, height: 600 };
        let isPanning = false;
        let startPanX = 0;
        let startPanY = 0;
        let startViewBoxX = 0;
        let startViewBoxY = 0;
        let selectedElement = null;
        let selectedId = null;

        // Show loading state
        document.getElementById("loading").classList.add("show");

        // Load JSON data
        function loadPathwayData() {
            return fetch("/data/pathway_data.json")
                .then(response => {
                    if (!response.ok) throw new Error(`HTTP error! status: ${response.status}`);
                    return response.json();
                })
                .catch(error => {
                    console.error("Error loading JSON:", error);
                    // Fallback test data
                    return {
                        protbox_data: [
                            { "protbox_id": "pb1", "x": 100, "y": 100, "width": 46, "height": 17, "proteins": ["P12345"], "backup_label": "Protein1" }
                        ],
                        protein_data: {
                            "P12345": { "label": "Protein1", "fc_color_1": [128, 128, 128], "label_color": [0, 0, 0] }
                        },
                        groups: [],
                        arrows: [],
                        general_data: {
                            settings: {
                                "prot_label_size": 12,
                                "prot_label_font": "Arial",
                                "ptm_circle_radius": 5,
                                "ptm_label_size": 10,
                                "ptm_label_font": "Arial"
                            }
                        }
                    };
                });
        }

        // Initialize SVG
        function initializeSvg() {
            console.log("Debug: Initializing SVG");
            const canvas = document.getElementById("svgCanvas");
            if (!canvas) {
                console.error("Debug: SVG canvas not found");
                return;
            }
            svg = SVG().addTo("#svgCanvas").size("100%", "100%");
            draw = SVG.Svg(canvas);
            svg.viewbox(0, 0, viewBox.width, viewBox.height);
            console.log(`Debug: SVG initialized with size ${viewBox.width} x ${viewBox.height}`);

            loadPathwayData().then(data => {
                renderPathway(data);
                setupInteractivity();
                document.getElementById("loading").classList.remove("show");
                document.querySelector(".svg-container").style.display = "block";
            }).catch(error => {
                document.getElementById("error").textContent = `Error: Failed to load pathway data. ${error.message}. Using fallback data.`;
                document.getElementById("error").classList.add("show");
                document.getElementById("loading").classList.remove("show");
                renderPathway({
                    protbox_data: [{ "protbox_id": "pb1", "x": 100, "y": 100, "width": 46, "height": 17, "proteins": ["P12345"], "backup_label": "Protein1" }],
                    protein_data: { "P12345": { "label": "Protein1", "fc_color_1": [128, 128, 128], "label_color": [0, 0, 0] } },
                    groups: [], arrows: [], general_data: { settings: { "prot_label_size": 12, "prot_label_font": "Arial", "ptm_circle_radius": 5, "ptm_label_size": 10, "ptm_label_font": "Arial" } }
                });
                setupInteractivity();
                document.querySelector(".svg-container").style.display = "block";
            });
        }

        // Render pathway
        function renderPathway(data) {
            const protboxData = data.protbox_data || [];
            protboxData.forEach(pb => {
                const rect = draw.rect(pb.width, pb.height)
                    .move(pb.x, pb.y)
                    .fill('none')
                    .stroke({ width: 1, color: '#000' })
                    .attr('id', pb.protbox_id);
                if (pb.proteins && pb.proteins.length > 0) {
                    draw.text(pb.backup_label || pb.proteins[0])
                        .move(pb.x + 2, pb.y + 2)
                        .font({ size: 12, family: 'Arial' })
                        .fill('#000');
                }
            });
        }

        // Setup interactivity
        function setupInteractivity() {
            const canvas = document.getElementById("svgCanvas");
            canvas.addEventListener("wheel", handleZoom);
            canvas.addEventListener("mousedown", handlePanStart);
            canvas.addEventListener("mousemove", handlePanMove);
            canvas.addEventListener("mouseup", handlePanEnd);
            canvas.addEventListener("mouseleave", handlePanEnd);
            document.addEventListener("keydown", handleKeyPress);
            document.getElementById("reset-view").addEventListener("click", resetView);

            canvas.addEventListener("click", function(e) {
                const rect = e.target.closest("rect");
                if (rect) {
                    if (selectedElement) selectedElement.stroke({ color: '#000', width: 1 });
                    selectedElement = rect;
                    selectedId = rect.getAttribute("id");
                    rect.stroke({ color: 'red', width: 2 });
                    console.log(`Selected ${selectedId} at ${e.clientX}, ${e.clientY}`);
                } else if (e.target === canvas && selectedElement) {
                    selectedElement.stroke({ color: '#000', width: 1 });
                    selectedElement = null;
                    selectedId = null;
                }
            });
        }

        // Zoom handler
        function handleZoom(e) {
            e.preventDefault();
            const zoomFactor = e.deltaY < 0 ? 1.1 : 0.9;
            zoomLevel = Math.min(maxZoom, Math.max(minZoom, zoomLevel * zoomFactor));
            svg.viewbox(viewBox.x, viewBox.y, viewBox.width / zoomLevel, viewBox.height / zoomLevel);
        }

        // Pan handlers
        function handlePanStart(e) {
            if (e.button === 0 && !e.target.closest("rect")) {
                isPanning = true;
                startPanX = e.clientX;
                startPanY = e.clientY;
                startViewBoxX = viewBox.x;
                startViewBoxY = viewBox.y;
            }
        }

        function handlePanMove(e) {
            if (isPanning) {
                const dx = (e.clientX - startPanX) / zoomLevel;
                const dy = (e.clientY - startPanY) / zoomLevel;
                viewBox.x = startViewBoxX - dx;
                viewBox.y = startViewBoxY - dy;
                svg.viewbox(viewBox.x, viewBox.y, viewBox.width / zoomLevel, viewBox.height / zoomLevel);
            }
        }

        function handlePanEnd() {
            isPanning = false;
        }

        // Key press handler
        function handleKeyPress(e) {
            if (selectedElement) {
                const step = e.shiftKey ? 10 : 1;
                let newX = parseFloat(selectedElement.attr("x"));
                let newY = parseFloat(selectedElement.attr("y"));
                switch (e.key) {
                    case "ArrowUp": newY -= step; break;
                    case "ArrowDown": newY += step; break;
                    case "ArrowLeft": newX -= step; break;
                    case "ArrowRight": newX += step; break;
                    case "Escape":
                        selectedElement.stroke({ color: '#000', width: 1 });
                        selectedElement = null;
                        selectedId = null;
                        break;
                }
                if (newX !== parseFloat(selectedElement.attr("x")) || newY !== parseFloat(selectedElement.attr("y"))) {
                    selectedElement.move(newX, newY);
                }
            }
        }

        // Reset view
        function resetView() {
            zoomLevel = 1;
            viewBox.x = 0;
            viewBox.y = 0;
            svg.viewbox(0, 0, viewBox.width, viewBox.height);
        }

        // Initialize when DOM is ready
        document.addEventListener("DOMContentLoaded", initializeSvg);
    </script>
</body>
</html>