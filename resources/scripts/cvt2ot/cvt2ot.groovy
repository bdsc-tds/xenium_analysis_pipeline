System.setProperty("java.awt.headless", "true")

import qupath.lib.images.writers.ome.OMEPyramidWriter
import qupath.lib.images.servers.bioformats.BioFormatsServerBuilder
import loci.formats.ImageReader

// Setup output directory
def outputDir = null
if (args != null && args.size() >= 2 && args[0] == "output") {
    outputDir = args[1]
} else {
    try {
        outputDir = buildFilePath(PROJECT_BASE_DIR, "Exported_Regions")
    } catch (Exception e) {
        print "No output directory specified in args\n"
        return
    }
}

mkdirs(outputDir)
print "Output directory: ${outputDir}\n"

def project = getProject()

// Identify images to process
def imagesToProcess = []
if (project != null) {
    imagesToProcess = project.getImageList()
} else {
    // If no project, get the currently open image
    def currentImageData = getCurrentImageData()
    if (currentImageData != null) {
        imagesToProcess = [currentImageData]
    } else {
        print "No project and no image open. Nothing to do!"
        return
    }
}

// Loop through every image in the project
for (item in imagesToProcess) {

    def server = item.getServer()
    def uri = server.getURIs()[0]
    
    // Create a BioFormats server builder to access all series
    def builder = new BioFormatsServerBuilder()
    int nSeries = 1
    try {
        def reader = new ImageReader()
        def path = new File(uri).getAbsolutePath()
        reader.setId(path)
        nSeries = reader.getSeriesCount()
        reader.close()
        print "Detected ${nSeries} series in ${path}\n"
    } catch (Exception e) {
        print "Could not determine series count using Bio-Formats ImageReader: ${e.getMessage()}\n"
    }
    
    int seriesIndex = 0
    while(seriesIndex < nSeries) {
        try {
            def sceneServer = null
            if (seriesIndex == 0) {
                sceneServer = server
            } else {
                // Build a server for the specific series index using BioFormats args
                def bfArgs = ["--series", seriesIndex.toString()]
                sceneServer = builder.buildServer(uri, bfArgs as String[])
            }
            
            if (sceneServer == null) {
                break
            }
            
            print "Series ${seriesIndex}: ${sceneServer.getWidth()}x${sceneServer.getHeight()} pixels\n"
            
            // Skip small utility images (Labels/Thumbnails)
            if (sceneServer.getWidth() < 1000 || sceneServer.getHeight() < 1000) {
                print "Series ${seriesIndex}: Skipping (Utility image)\n"
                // Only close if it's NOT the primary server
                if (seriesIndex > 0) sceneServer.close()
                seriesIndex++
                continue
            }

            // Build a robust filename using the image URI (original file) and metadata
            def metaName = sceneServer.getMetadata().getName()
            def uriStr = uri.toString()

            // Extract original filename from URI (strip any path)
            def origFilename = uriStr.replaceAll('^.*[\\/]', '')
            try {
                origFilename = java.net.URLDecoder.decode(origFilename, 'UTF-8')
            } catch(Exception e) {
                // ignore decoding errors and use raw value
            }

            // Base name without extension (from original file)
            def baseName = origFilename.replaceAll(/\.[^.]+$/, '')

            // Decide whether to include a region/series suffix
            def includeSuffix = (nSeries > 1) || (metaName != null && metaName != '' && metaName != origFilename)

            def suffix = "series${seriesIndex}"
            if (metaName != null && metaName != '') {
                suffix = metaName.replaceAll(/\.[^.]+$/, '')
            }

            // Sanitize parts: replace whitespace and unsafe chars with underscores
            def sanitize = { s ->
                if (s == null) return ''
                s = s.toString()
                s = s.replaceAll("\\s+", "_")
                s = s.replaceAll("[^A-Za-z0-9_\\-\\.]", "_")
                return s
            }

            def cleanBase = sanitize(baseName)
            def cleanSuffix = sanitize(suffix)

            // If multiple series are present, always keep a disambiguating suffix.
            // This avoids overwriting outputs when metadata names are identical.
            if (cleanSuffix == '' || cleanSuffix == cleanBase) {
                if (nSeries > 1) {
                    cleanSuffix = "series${seriesIndex}"
                } else {
                    includeSuffix = false
                }
            }

            def finalName = includeSuffix ? "${cleanBase}_${cleanSuffix}.ome.tif" : "${cleanBase}.ome.tif"
            
            // Define the output path
            def outputPath = buildFilePath(outputDir, finalName)
            
            print "Starting export for: ${finalName}\n"
            
            // Build the writer with pyramid requirements
            new OMEPyramidWriter.Builder(sceneServer)
                .downsamples(1, 2, 4, 8, 16, 32) // Pyramid scale 2
                .tileSize(1024)
                .compression(OMEPyramidWriter.CompressionType.ZLIB)
                .parallelize(true)
                .build()
                .writePyramid(outputPath)
                
            print "Successfully exported: ${outputPath}\n"
            
            // Only close if it's NOT the primary server
            if (seriesIndex > 0) sceneServer.close()
            seriesIndex++
            
        } catch (Exception e) {
            def msg = e.getMessage() ?: ""
            if (msg.contains("Invalid series") || msg.contains("index=")) {
                print "Finished processing available series for ${uri}\n"
                break
            } else {
                print "ERROR: Failed to process series ${seriesIndex}. Error: ${msg}\n"
                e.printStackTrace()
                seriesIndex++
            }
        }
    }
}

print "All regions processed.\n"

