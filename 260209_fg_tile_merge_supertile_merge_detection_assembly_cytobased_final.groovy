import static qupath.lib.gui.scripting.QPEx.*

import qupath.lib.objects.PathObjects
import qupath.lib.roi.ROIs
import qupath.lib.roi.RoiTools
import qupath.ext.biop.cellpose.Cellpose2D
import qupath.lib.analysis.features.ObjectMeasurements
import qupath.lib.regions.ImagePlane

import java.awt.geom.Area

// =======================
// PARAMETERS
// =======================
double tileSizeMicrons = 500
double overlapMicrons  = 30
double downsample = 1.0
double superTileSizeMicrons = 3730 // --> 8 tiles = 1 super tile

// =======================
// IMAGE INFO
// =======================
def imageData = getCurrentImageData()
def server = imageData.getServer()
def cal = server.getPixelCalibration()
double pxSize = cal.getAveragedPixelSizeMicrons()

// Global containers
def cytos     = []
def nucs_omni = []
def nucs_qp   = []

// =======================
// STAGE 0 — FOREGROUND DETECTION (FULL IMAGE, NO MERGE)
// =======================
println "STAGE 0 — Foreground detection (full image, unmerged)"

resetSelection()
createFullImageAnnotation(true)

runPlugin(
    'qupath.imagej.detect.cells.WatershedCellDetection',
    '{"detectionImage":"AF488",' +
    '"requestedPixelSizeMicrons":0.5,' +
    '"backgroundRadiusMicrons":5.0,' +
    '"sigmaMicrons":3.0,' +
    '"minAreaMicrons":1.0,' +
    '"maxAreaMicrons":1000.0,' +
    '"threshold":40.0,' +
    '"watershedPostProcess":false,' +
    '"cellExpansionMicrons":8.0}'
)

// Convert each detection to its OWN foreground annotation
def fgAnnotations = []
getDetectionObjects().each { d ->
    fgAnnotations << PathObjects.createAnnotationObject(
        d.getROI(),
        getPathClass("Foreground")
    )
}

removeDetections()
addObjects(fgAnnotations)
fireHierarchyUpdate()

println "✔ Created ${fgAnnotations.size()} foreground annotations"


// =======================
// STAGE 1 — CREATE TILES
// =======================
println "STAGE 1 — Creating tiles"

clearDetections()
def width  = server.getWidth()
def height = server.getHeight()

double tilePx    = tileSizeMicrons / pxSize
double overlapPx = overlapMicrons / pxSize
double step      = tilePx - overlapPx

def tiles = []

for (double y = 0; y < height; y += step) {
    for (double x = 0; x < width; x += step) {
        double w = Math.min(tilePx, width - x)
        double h = Math.min(tilePx, height - y)
        def roi = ROIs.createRectangleROI(x, y, w, h, null)
        def tileAnno = PathObjects.createAnnotationObject(roi, getPathClass("Tile"))
        addObject(tileAnno)
        tiles << tileAnno
    }
}
fireHierarchyUpdate()
println "Created ${tiles.size()} tiles"

// =======================
// STAGE 2 — MERGE FOREGROUNDS PER TILE
// =======================
println "STAGE 2 — Merge foregrounds per tile"

def allForegrounds = getAnnotationObjects()
    .findAll { it.getPathClass()?.getName() == "Foreground" }

def mergedTileAnnotations = []

tiles.eachWithIndex { tile, idx ->
    println "Processing tile ${idx+1}/${tiles.size()}"

    Area tileArea = new Area()

    allForegrounds.each { fg ->
        def inter = RoiTools.intersection(fg.getROI(), tile.getROI())
        if (inter != null && inter.getArea() > 0)
            tileArea.add(new Area(inter.getShape()))
    }

    if (tileArea.isEmpty())
        return

    def plane = tile.getROI().getImagePlane()
    def mergedROI = ROIs.createAreaROI(tileArea, plane)

    def fgTile = PathObjects.createAnnotationObject(
        mergedROI,
        getPathClass("ForegroundTile")
    )

    addObject(fgTile)
    mergedTileAnnotations << fgTile
}

fireHierarchyUpdate()
println "✔ Created ${mergedTileAnnotations.size()} per-tile foregrounds"


// =======================
// STAGE 3 — Build supertiles
// =======================
println "STAGE 3 — Hierarchical merge (super-tiles)"

double superTilePx = superTileSizeMicrons / pxSize
def superMergedForegrounds = []

// Counters for tracking
int superTileIndex = 0
int numY = (int) Math.ceil(height / superTilePx)
int numX = (int) Math.ceil(width / superTilePx)

for (int yi = 0; yi < numY; yi++) {
    for (int xi = 0; xi < numX; xi++) {

        double x = xi * superTilePx
        double y = yi * superTilePx
        double w = Math.min(superTilePx, width - x)
        double h = Math.min(superTilePx, height - y)

        def superROI = ROIs.createRectangleROI(x, y, w, h, ImagePlane.getDefaultPlane())
        Area superArea = new Area()

        mergedTileAnnotations.each { fg ->
            def inter = RoiTools.intersection(superROI, fg.getROI())
            if (inter != null && inter.getArea() > 0) {
                superArea.add(new Area(inter.getShape()))
            }
        }

        if (!superArea.isEmpty()) {
            def roi = ROIs.createAreaROI(superArea, superROI.getImagePlane())
            superMergedForegrounds << roi
        }

        superTileIndex++
        println "Processed supertile ${superTileIndex} at (x=${(int)x}, y=${(int)y}), width=${(int)w}, height=${(int)h}"
    }
}

println "Created ${superMergedForegrounds.size()} merged super-tile foreground ROIs"

// =======================
// STAGE 4 — CREATE SUPER-TILE ANNOTATIONS
// =======================
println "STAGE 4 — Creating super-tile foreground annotations"

def superTileFGAnnos = []

superMergedForegrounds.eachWithIndex { roi, i ->
    def anno = PathObjects.createAnnotationObject(
        roi,
        getPathClass("ForegroundSuperTile")
    )
    addObject(anno)
    superTileFGAnnos << anno
    println "Processing supertile ${i+1}/${superMergedForegrounds.size()}"
}

fireHierarchyUpdate()
println "Created ${superTileFGAnnos.size()} super-tile foreground annotations"


// =======================
// STAGE 5 — Cell detection per supertile
// =======================
println "STAGE 5 — Cell detection per supertile"

superTileFGAnnos.eachWithIndex { fgAnno, idx ->
    println "Running Cellpose on supertile ${idx+1}/${superTileFGAnnos.size()}"

    resetSelection()
    selectObjects(fgAnno)


    def cytoModel = Cellpose2D.builder('bact_fluor_omni')
        .pixelSize(0.3)
        .setOverlap(30)
        .channels('DAPI', 'AF488')
        .tileSize(4096)
        .cellposeChannels(2, 1)
        .useOmnipose()
        .build()

    cytoModel.detectObjects(imageData, getSelectedObjects())

    getDetectionObjects().each { d ->
        d.setPathClass(getPathClass('Cytoplasm'))
        cytos << d
    }

    clearDetections()
    System.gc()
}
removeObjects(mergedTileAnnotations, false)
removeObjects(superTileFGAnnos, false)


// =======================
// STAGE 6 — Nuclei detection
// =======================
println "STAGE 6 — Nuclei detection"
resetSelection()
createFullImageAnnotation(true)

runPlugin(
    'qupath.imagej.detect.cells.WatershedCellDetection',
    '{"detectionImage":"DAPI",' +
    '"requestedPixelSizeMicrons":0.5,' +
    '"backgroundRadiusMicrons":2.0,' +
    '"sigmaMicrons":1.5,' +
    '"minAreaMicrons":0.5,' +
    '"maxAreaMicrons":170.0,' +
    '"threshold":50.0,' +
    '"watershedPostProcess":true,' +
    '"cellExpansionMicrons":0.0}'
)

nucs_qp.clear()
getDetectionObjects().each { d ->
    d.setPathClass(getPathClass('Nucleus_QuPath'))
    nucs_qp << d
}
clearDetections()
println "✔ Nuclei detected: ${nucs_qp.size()}"

// =======================
// STAGE 7 — Cell assembly
// =======================
println "STAGE 7 — Cell assembly (foreground-gated)"
addObjects(cytos)  // make sure cytos contains all Cellpose cytoplasms

if (nucs_qp.isEmpty()) {
    println "⚠ No nuclei detected — cannot assemble cells!"
} else {

    // Use cytoplasms as foregrounds for gating
    def foregrounds = cytos
    println "✔ Using ${foregrounds.size()} cytoplasm objects as foregrounds for gating"

    // Index nuclei by cytoplasm
    def nucToCyto = [:]  // map each nucleus → its cytoplasm

    nucs_qp.each { nuc ->
        def nROI = nuc.getROI()
        if (nROI == null) return

        def bestCyto = null
        double bestOverlap = 0
        foregrounds.each { cyto ->
            def inter = RoiTools.intersection(nROI, cyto.getROI())
            if (inter != null && inter.getArea() > bestOverlap) {
                bestOverlap = inter.getArea()
                bestCyto = cyto
            }
        }
        if (bestOverlap > 0) nucToCyto[nuc] = bestCyto
    }

    // Group nuclei by cytoplasm
    def cytoToNucs = [:].withDefault { [] }
    nucToCyto.each { nuc, cyto -> cytoToNucs[cyto] << nuc }

    // Create cells
    def newCells = []

    // Cytoplasm-based cells
    cytoToNucs.each { cyto, nucleiInCyto ->
        def cell = PathObjects.createCellObject(
            cyto.getROI(),
            nucleiInCyto[0].getROI(),
            getPathClass(nucleiInCyto.size() > 1 ? 'Cell_Binucleated' : 'Cell_Cyto'),
            null
        )
        nucleiInCyto.each { cell.addChildObject(it) }
        newCells << cell
    }

    // Nucleus-only cells (no overlapping cytoplasm)
    nucs_qp.each { nuc ->
        if (!nucToCyto.containsKey(nuc)) {
            def cell = PathObjects.createCellObject(nuc.getROI(), nuc.getROI(), getPathClass('Cell_NucleusOnly'), null)
            cell.addChildObject(nuc)
            newCells << cell
        }
    }

    // Cleanup
    removeObjects(cytos, false)
    removeObjects(nucs_qp, false)

    // Add cells & measure
    addObjects(newCells)
    fireHierarchyUpdate()
    println "✅ DONE — Created ${newCells.size()} cells (foreground-gated, Stage 7)"
}


