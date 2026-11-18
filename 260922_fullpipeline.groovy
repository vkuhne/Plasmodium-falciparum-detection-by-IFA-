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
double superTileSizeMicrons = 3790 // --> 4 tiles = 1 super tile

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
// STAGE 1 — Build supertiles
// =======================
println "STAGE 1 — Create supertiles"

double superTilePx = superTileSizeMicrons / pxSize

def superTileAnnos = []

int superTileIndex = 0
int numY = (int)Math.ceil(server.getHeight() / superTilePx)
int numX = (int)Math.ceil(server.getWidth() / superTilePx)

for (int yi = 0; yi < numY; yi++) {
    for (int xi = 0; xi < numX; xi++) {

        double x = xi * superTilePx
        double y = yi * superTilePx
        double w = Math.min(superTilePx, server.getWidth() - x)
        double h = Math.min(superTilePx, server.getHeight() - y)

        def roi = ROIs.createRectangleROI(x, y, w, h, ImagePlane.getDefaultPlane())
        def anno = PathObjects.createAnnotationObject(roi)

        superTileAnnos << anno

        superTileIndex++
        println "Created supertile ${superTileIndex} at (${(int)x}, ${(int)y})"
    }
}

addObjects(superTileAnnos)

println "✔ Created ${superTileAnnos.size()} supertiles"

// =======================
// STAGE 2 — Cell detection per supertile
// =======================
println "STAGE 2 — Cellpose detection per supertile"

def cytoModel = Cellpose2D.builder('bact_fluor_omni')
    .pixelSize(0.3)
    .setOverlap(30)
    .channels('DAPI', 'AF488')
    .tileSize(4096)
    .cellposeChannels(2, 1)
    .useOmnipose()
    .build()

superTileAnnos.eachWithIndex { tileAnno, idx ->

    println "Running Cellpose on supertile ${idx+1}/${superTileAnnos.size()}"

    resetSelection()
    selectObjects(tileAnno)

    cytoModel.detectObjects(imageData, getSelectedObjects())

    getDetectionObjects().each { d ->
        d.setPathClass(getPathClass('Cytoplasm'))
        cytos << d
    }

    clearDetections()

    // Prevent RAM buildup
    System.gc()
}

println "✔ Cytoplasms detected: ${cytos.size()}"

// remove supertiles after use
removeObjects(superTileAnnos, false)


// =======================
// STAGE 3 — Nuclei detection
// =======================
println "STAGE 3 — Nuclei detection"
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
// STAGE 4 — Cell assembly (spatially indexed)
// =======================
println "STAGE 4 — Cell assembly (spatially indexed)"

// Safety checks
if (cytos.isEmpty()) {
    println "❌ No cytoplasms found — cannot assemble cells"
    return
}
if (nucs_qp.isEmpty()) {
    println "❌ No nuclei found — cannot assemble cells"
    return
}

println "✔ Cytoplasms: ${cytos.size()}"
println "✔ Nuclei: ${nucs_qp.size()}"

// ----------------------------------------------------
// 1. Build nucleus spatial index (grid-based)
// ----------------------------------------------------
double gridSize = 50.0  // microns-equivalent in pixels (tune if needed)
def nucGrid = [:].withDefault { [] }

nucs_qp.each { nuc ->
    def r = nuc.getROI()
    if (r == null) return

    int gx = (int) Math.floor(r.getBoundsX() / gridSize)
    int gy = (int) Math.floor(r.getBoundsY() / gridSize)
    nucGrid[[gx, gy]] << nuc
}

println "✔ Built nucleus spatial grid (${nucGrid.size()} cells)"


// ----------------------------------------------------
// 2. Assign nuclei → cytoplasm (fast spatial pruning)
// ----------------------------------------------------
def nucToCyto = [:]

cytos.eachWithIndex { cyto, idx ->
    def cROI = cyto.getROI()
    if (cROI == null) return

    double cx0 = cROI.getBoundsX()
    double cy0 = cROI.getBoundsY()
    double cx1 = cx0 + cROI.getBoundsWidth()
    double cy1 = cy0 + cROI.getBoundsHeight()

    int gx0 = (int) Math.floor(cx0 / gridSize)
    int gy0 = (int) Math.floor(cy0 / gridSize)
    int gx1 = (int) Math.floor(cx1 / gridSize)
    int gy1 = (int) Math.floor(cy1 / gridSize)

    for (int gx = gx0; gx <= gx1; gx++) {
        for (int gy = gy0; gy <= gy1; gy++) {
            nucGrid[[gx, gy]].each { nuc ->
                if (nucToCyto.containsKey(nuc)) return

                def nROI = nuc.getROI()
                if (nROI == null) return

                // Fast bounding-box reject
                double nx0 = nROI.getBoundsX()
                double ny0 = nROI.getBoundsY()
                double nx1 = nx0 + nROI.getBoundsWidth()
                double ny1 = ny0 + nROI.getBoundsHeight()

                if (nx1 < cx0 || nx0 > cx1 || ny1 < cy0 || ny0 > cy1) return

                // Precise geometry check
                def inter = RoiTools.intersection(nROI, cROI)
                if (inter != null && inter.getArea() > 0) {
                    nucToCyto[nuc] = cyto
                }
            }
        }
    }

    if (idx % 1000 == 0 && idx > 0)
        println "  Processed ${idx}/${cytos.size()} cytoplasms"
}

println "✔ Assigned ${nucToCyto.size()} nuclei to cytoplasms"


// ----------------------------------------------------
// 3. Group nuclei by cytoplasm
// ----------------------------------------------------
def cytoToNucs = [:].withDefault { [] }
nucToCyto.each { nuc, cyto ->
    cytoToNucs[cyto] << nuc
}


// ----------------------------------------------------
// 4. Create cells
// ----------------------------------------------------
def newCells = []

// Cytoplasm-based cells
cytoToNucs.each { cyto, nuclei ->
    def cell = PathObjects.createCellObject(
        cyto.getROI(),
        nuclei[0].getROI(),
        getPathClass(nuclei.size() > 1 ? 'Cell_Binucleated' : 'Cell_Cyto'),
        null
    )
    nuclei.each { cell.addChildObject(it) }
    newCells << cell
}

// Nucleus-only cells
nucs_qp.each { nuc ->
    if (!nucToCyto.containsKey(nuc)) {
        def cell = PathObjects.createCellObject(
            nuc.getROI(),
            nuc.getROI(),
            getPathClass('Cell_NucleusOnly'),
            null
        )
        cell.addChildObject(nuc)
        newCells << cell
    }
}

println "✔ Cells created: ${newCells.size()}"


// ----------------------------------------------------
// 5. Cleanup + add cells
// ----------------------------------------------------
removeObjects(cytos, false)
removeObjects(nucs_qp, false)
addObjects(newCells)
// Measure intensity & shape
def measurements = ObjectMeasurements.Measurements.values() as List
def compartments = ObjectMeasurements.Compartments.values() as List
def shapes = ObjectMeasurements.ShapeFeatures.values() as List
newCells.each { cell ->
    ObjectMeasurements.addIntensityMeasurements(server, cell, downsample, measurements, compartments)
    ObjectMeasurements.addCellShapeMeasurements(cell, cal, shapes)
}
fireHierarchyUpdate()

println "✅ DONE — Stage 4 completed successfully"

runObjectClassifier("260917_classifierA")