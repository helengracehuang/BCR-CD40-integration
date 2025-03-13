//-----------------------------------------------------------
// IMAGEJ MACRO: Count DRAQ7+ (dying) cells vs. H2B+ total
// Based on overlaying nucleus ROIs on top of a DRAQ7 mask
// Similar structure and commands to your original script
//-----------------------------------------------------------

//-----------------------------------------
// USER PARAMETERS (tune for your data)
//-----------------------------------------
// Only use 'var' for global variables! 
// (Avoid declaring them inside macros/functions in ImageJ.)
var homeDir = "/Users/helenhuang/Documents/1st Year PhD/Hoffmann Lab/Microscopy/VK014/";
var subDir  = "VK014_CD40HighBCRHigh_24h/"; // CHANGE to your folder
var gridIter = 4;            // CHANGE every time
var label    = "cellID";     // "cellID" or "morph+ID"

// Channel assignments in your hyperstack (original image):
CH_DEATH   = 1; // e.g. DRAQ7 (APC)
CH_NUCLEUS = 2; // e.g. H2B (mCherry)

// Analyze Particles settings for nuclei
var MIN_SIZE = 50;   // example, in pixels
var MAX_SIZE = 1e5;  // large upper bound
var MIN_CIRC = 0.10;
var MAX_CIRC = 1.00;

// Border exclusion in pixels
var NPIXELS  = 25;

// Auto-threshold method
var THRESHOLD_STYLE = "Auto Threshold";
var DC_AUTO_METHOD  = "Moments";  // For dead cell channel
var NUC_AUTO_METHOD = "Moments";     // For nucleus channel

// We'll decide a "DeadMean" threshold to call it "dead"
var DEAD_MEAN_THRESHOLD = 2; // Adjust to your data

//-----------------------------------------
// Check if the output file already exists
//-----------------------------------------
resultsFile = homeDir + subDir + "Results_death" + gridIter + ".csv";
if (File.exists(resultsFile)) {
    waitForUser("Step 0: Warning! Are you sure?",
    "This file already exists. Click 'Cancel' to abort, or 'OK' to proceed.");
}
// Delete previous segmentation
nRois = roiManager("count");
for (i = (nRois-1); i > 0; i--) {
	roiManager("select", i);
	roiManager("Delete");
}
//-----------------------------------------
// Initialize measurement + ROI defaults
//-----------------------------------------
run("Set Measurements...", "area mean modal centroid median display redirect=None decimal=3");
run("Roi Defaults...", "color=pink stroke=1 group=0 name=AOI");
run("Roi Defaults...", "color=green stroke=1 group=1 name=nucleus");
run("Clear Results");

//-----------------------------------------
// STEP 0: Verify the "AOI" ROI exists
//-----------------------------------------
waitForUser("Step 0: ROI Setup",
"Please box an area for analysis, name it 'AOI', and add it to the ROI Manager.\n" +
"Click 'OK' once you have added the 'AOI' ROI.");

if (roiManager("count") == 0) {
   waitForUser("Step 0: No ROIs found!",
   "Please add an 'AOI' ROI to the ROI Manager, then click 'OK'.");
}

var nRois = roiManager("count");
var aoiIndex = -1;
for (i = 0; i < nRois; i++) {
    name = RoiManager.getName(i);
    if (name == "AOI") {
        aoiIndex = i;
        break;
    }
}
if (aoiIndex < 0) {
    exit("ROI named 'AOI' not found in ROI Manager. Aborting.");
}
print("AOI index = " + aoiIndex);

//-----------------------------------------
// STEP 1: Duplicate the DRAQ7 channel (CH_DEATH)
//         and create a DeadCellMask
//-----------------------------------------
selectWindow("Composite");   // The multi-channel image should be named "Composite"
roiManager("select", aoiIndex);

// Duplicate channel for DRAQ7
run("Duplicate...", "title=AOI_C1 duplicate channels=" + CH_DEATH);
selectWindow("AOI_C1");
run("8-bit");
rename("DeadCell_Raw");

// Let user adjust threshold if desired
waitForUser("Step 1: Dead cell channel thresholding",
"Adjust threshold parameters as needed, then click OK to proceed.");

run(THRESHOLD_STYLE, "method=" + DC_AUTO_METHOD + " ignore_black white");
run("Convert to Mask");
rename("DeadCellMask");

// A bit of morphological cleanup
run("Open");  // Erode then Dilate
run("Erode");
run("Fill Holes");

//-----------------------------------------
// STEP 2: Duplicate the H2B channel (CH_NUCLEUS)
//         for nucleus segmentation
//-----------------------------------------
selectWindow("Composite");
roiManager("select", aoiIndex);
run("Duplicate...", "title=AOI_C2 duplicate channels=" + CH_NUCLEUS +"-"+ CH_DEATH);
selectWindow("Composite");
run("Duplicate...", "title=Nucleus_Raw duplicate channels=" + CH_NUCLEUS);
selectWindow("Nucleus_Raw");
run("8-bit");

//-----------------------------------------
// STEP 3: Create a BlackMask to exclude edges
//-----------------------------------------
selectWindow("Nucleus_Raw");
run("Duplicate...", "title=BlackRegion duplicate");
selectWindow("BlackRegion");
setThreshold(0, 0);
run("Convert to Mask");
rename("BlackMask");

// Expand black mask for border exclusion
run("Duplicate...", "title=BlackMask_expanded");
selectWindow("BlackMask_expanded");
for (d = 0; d < NPIXELS; d++) {
    run("Dilate");
}
// Invert the original so that white = valid region
selectWindow("BlackMask");
run("Invert");

//-----------------------------------------
// STEP 4: Nucleus segmentation
//-----------------------------------------
selectWindow("Nucleus_Raw");
waitForUser("Step 4: Nucleus thresholding",
"Adjust threshold if needed, then click OK to proceed.");

run(THRESHOLD_STYLE, "method=" + NUC_AUTO_METHOD + " ignore_black white");
run("Convert to Mask");
run("Fill Holes");
//run("Watershed");

// Analyze Particles for nuclei
run("Analyze Particles...",
    "size=" + MIN_SIZE + "-" + MAX_SIZE +
    " pixel circularity=" + MIN_CIRC + "-" + MAX_CIRC +
    " exclude add");

nucCountStart = roiManager("count");
nNUC = 0;

// Remove nuclei overlapping black edges
for (i = nucCountStart - 1; i >= 1; i--) {
    selectWindow("BlackMask_expanded");
    roiManager("select", i);
    roiManager("Measure");
    GridMed = getResult("Median", 0);
    run("Clear Results");

    if (GridMed > 128) {
        roiManager("Delete");
    } else {
        nNUC++;
        newName = "nuc," + nNUC;
        roiManager("rename", newName);
        RoiManager.setGroup(1);
    }
}

waitForUser("Step 4: Nucleus segmentation completed",
nNUC + " nuclei found. Inspect them, then press OK to continue.");

//-----------------------------------------
// STEP 5: Determine which nuclei are DRAQ7+
//         by measuring in DeadCellMask
//-----------------------------------------
var dyingCount = 0;
for (j = 1; j < nNUC; j++) { // excluding AOI
    // Nucleus index in manager
    nucIndex = j;  // These are all nucleus ROIs now
    selectWindow("DeadCellMask");
    roiManager("select", nucIndex);
    roiManager("Measure");
    DeadMean = getResult("Mean", 0);
    run("Clear Results");

    // Decide if nucleus is "dead" by Mean intensity
    if (DeadMean > DEAD_MEAN_THRESHOLD) {
        dyingCount++;
        oldName = RoiManager.getName(nucIndex);
        newName = "dead_" + oldName;
        roiManager("rename", newName);
    }
}
fractionDying = dyingCount / nNUC;

//-----------------------------------------
// STEP 6: Output and save results
//-----------------------------------------
print("Total nuclei (H2B+): " + nNUC);
print("DRAQ7+ (dying) nuclei: " + dyingCount);
print("Fraction dying: " + fractionDying);

// Put them in the Results table row 0
setResult("TotalCells", 0, nNUC);
setResult("DyingCells", 0, dyingCount);
setResult("FractionDying", 0, fractionDying);

selectWindow("AOI_C2");
roiManager("Show All");
waitForUser("Step 6: Save results",
"Check everything, then click OK to save results and images.");

// If file exists, increment gridIter
resultsFile = homeDir + subDir + "Results_death" + gridIter + ".csv";
if (File.exists(resultsFile)) {
    waitForUser("Warning!",
    "This file exists. Click 'Cancel' to abort, or 'OK' to increment gridIter.");
    gridIter = gridIter + 1;
    resultsFile = homeDir + subDir + "Results_death" + gridIter + ".csv";
}

// Save Results table
saveAs("Results", resultsFile);

// Save ROI set
roiManager("Save", homeDir + subDir + "RoiSet_death" + gridIter + ".zip");

// Save final images
selectWindow("DeadCellMask");
saveAs("tiff", homeDir + subDir + "DeadCellMask_death" + gridIter + ".tif");
selectWindow("Nucleus_Raw");
saveAs("tiff", homeDir + subDir + "Nucleus_death" + gridIter + ".tif");
selectWindow("AOI_C2");
saveAs("tiff", homeDir + subDir + "Scene_death" + gridIter + ".tif");

// Cleanup
selectWindow("BlackMask"); close();
selectWindow("BlackMask_expanded"); close();
selectWindow("DeadCellMask_death" + gridIter + ".tif"); close();
selectWindow("Nucleus_death" + gridIter + ".tif"); close();
selectWindow("Scene_death" + gridIter + ".tif"); close();

print("Done! " + nNUC + " total nuclei; " + dyingCount +
      " are DRAQ7+. Fraction: " + fractionDying);
roiManager("Deselect");
roiManager("Delete");
/*
--------------------------------------------------------------------------------
HOW THIS SCRIPT WORKS:

1) Duplicate the DRAQ7 channel to create "DeadCell_Raw", then threshold it => 
   "DeadCellMask". (We do not run Analyze Particles on DRAQ7; we just need the mask.)

2) Duplicate the H2B channel => "Nucleus_Raw", threshold, and use Analyze Particles 
   to segment individual nuclei as "nuc,X". Remove those on black edges.

3) For each nucleus ROI, measure the Mean intensity in "DeadCellMask". 
   If it's above DEAD_MEAN_THRESHOLD, label it as "dead_nuc,X".

4) Count how many are above threshold => fraction of dying cells.

5) Save the Results table, ROI set, and images. 
   Adjust morphological steps, thresholds, and sizes as needed.
--------------------------------------------------------------------------------
*/