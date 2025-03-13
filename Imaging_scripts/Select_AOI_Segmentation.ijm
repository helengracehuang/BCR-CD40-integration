//-----------------------------------------
// USER PARAMETERS (tune for your data)
//-----------------------------------------
// Only use 'var' for global variables! Do not use them inside macro or 
// function code blocks. Using 'var' in a macro or f=unction may cause it to fail.
// ImageJ runs all the global variable declaration before executing other code.

var homeDir = "/Users/helenhuang/Documents/1st Year PhD/Hoffmann Lab/Microscopy/VK014/";
var subDir = "VK014_CD40High_24h/"; // CHANGE!!!!! to condition folder
var gridIter = 5; // CHANGE every time!!!!!
var label = "cellID";// "cellID" or both "morph+ID" for final Results table

// Channels in your hyperstack (original image):
CH_DEATH       = 1; // channel for death marker, e.g. DRAQ7 (red)
CH_WHOLECELL_1 = 1 +CH_DEATH; // channel for whole cell marker, e.g. cRel (cyan)
CH_NUCLEUS     = 2 +CH_DEATH; // channel for nucleus marker, e.g. H2B (magenta)
CH_WHOLECELL_2 = 3 +CH_DEATH; // 2nd channel for whole cell marker, e.g. RelA (yellow)
var includeH2BforWcChannel = true; // include H2B channel for whole cell segmentation
var includeDeathChannel = true; // CHANGE!!!!! if DRAQ7-APC channel for dead cells

// Analyze Particles settings for cells/nuclei
var MIN_SIZE = 50;    // example, in pixels
var MAX_SIZE = 1e4;    // large upper bound for later timepoint
var MIN_CIRC = 0.50;   
var MAX_CIRC = 1.00;   // perfectly round

var NPIXELS = 25;  // (radius of biggest cell) border exclusion pixels

// Auto-threshold method you prefer (e.g., Otsu, Default, etc.)
var THRESHOLD_STYLE = "Auto Threshold"; // "Auto Threshold" or "Auto Local Threshold"
var DC_AUTO_METHOD = "Moments"; // Dead cell auto thresholding
var WC_AUTO_METHOD = "Otsu"; // "Otsu" for 24h, "Moments" for C-h B-h 72h
var NUC_AUTO_METHOD = "IsoData"; // "Moments" for most, "Otsu" for some (extreme H2B signal), 
//"Shanbhag" or "Intermodes" for much later timepoint

// Morphological processing
var WC_EROSION = 2;   // number of erosions after watershed for whole cell
var NUC_EROSION = 1;   // number of erosions after watershed for nucleus

// Check if file iter exists yet
resultsFile = homeDir+subDir+"Results_auto"+gridIter+".csv";
if(File.exists(resultsFile)){
	waitForUser("Step 0: Warning! Are you sure?", 
	"This file exists. Click 'cancel' to abort now, or click 'OK' to proceed.");
}

run("Set Measurements...", "area mean modal centroid median display redirect=None decimal=3");
run("Roi Defaults...", "color=pink stroke=1 group=0 name=AOI"); // AOI & BG
run("Roi Defaults...", "color=green stroke=1 group=1 name=wholeCell"); // whole cell
run("Roi Defaults...", "color=red stroke=1 group=2 name=nucleus"); // nucleus
run("Roi Defaults...", "color=blue stroke=1 group=3 name=cytoplasm"); // cytoplasm
run("Clear Results");
//-----------------------------------------
// STEP 0: Check that the "AOI" ROI exists
//-----------------------------------------
waitForUser("Step 0: Please box an area (ROI) for analysis and rename it 'AOI'.\n" + 
"Try to select scenes that are not overly dense. Click 'OK' to proceed.");

if (roiManager("count") == 0) {
   waitForUser("Step 0: No ROIs in ROI Manager. Please add 'AOI' first, then click 'OK'.");
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
   exit("ROI named 'AOI' not found in ROI Manager.");
}
print("The index of the AOI is "+aoiIndex);
//-----------------------------------------
// STEP 1: Duplicate the AOI region 
//         and build a dead cell mask (optional)
//-----------------------------------------
selectWindow("Composite");
roiManager("select", aoiIndex);

run("Duplicate...", 
    "title=AOI_C567 duplicate channels="+CH_WHOLECELL_1+"-"+CH_WHOLECELL_2);
delta = CH_WHOLECELL_1;
CH_WHOLECELL_1 = CH_WHOLECELL_1-delta+1;
CH_NUCLEUS = CH_NUCLEUS-delta+1;
CH_WHOLECELL_2 = CH_WHOLECELL_2-delta+1; // relative channel index

if (includeDeathChannel) { // if DRAQ7 was taken
	// duplicate the AOI region of APC channel
	selectWindow("Composite");
	roiManager("select", aoiIndex);
	run("Duplicate...", 
    "title=AOI_C1 duplicate channels="+CH_DEATH);
    selectWindow("AOI_C1");
    
	run("Duplicate...", "title=DeadCell_Raw duplicate channels="+CH_DEATH);
	selectWindow("DeadCell_Raw");
	run("8-bit");
	waitForUser("Step 1: Dead cell segmentation beginning", 
	"You might want to try all thresholding algorithms before proceeding.");
	
	// segmentation to identify dead cells
	run(THRESHOLD_STYLE, "method="+DC_AUTO_METHOD+" ignore_black white");
	run("Convert to Mask");
	rename("DeadCellMask");
	run("Open"); // Erosion + dilation: smoothes objects and removes isolated pixels.
	run("Erode");
	run("Fill Holes");
	
	waitForUser("Step 1: Dead cell segmentation completed", 
	"Please inspect DeadCellMask before proceeding.");
	selectWindow("AOI_C1"); close();
}

//-----------------------------------------
// STEP 2: Build a "WholeCell_Raw" image 
//-----------------------------------------
// Switch to the AOI-cropped hyperstack
selectWindow("AOI_C567");
run("Duplicate...", "title=Temp_C5 duplicate channels="+CH_WHOLECELL_1);
selectWindow("AOI_C567");
run("Duplicate...", "title=Temp_C6 duplicate channels="+CH_NUCLEUS);
selectWindow("AOI_C567");
run("Duplicate...", "title=Temp_C7 duplicate channels="+CH_WHOLECELL_2);

// Sum channels (Image Calculator => Max)
run("Image Calculator...", "image1=Temp_C5 operation=Max image2=Temp_C7 create");
rename("WholeCell_Temp");
if (includeH2BforWcChannel) { // include H2B channel for WC segmentation
	run("Image Calculator...", "image1=WholeCell_Temp operation=Max image2=Temp_C6 create");
	selectWindow("WholeCell_Temp"); close();
}
rename("WholeCell_Raw");

// Cleanup
selectWindow("Temp_C5"); close();
selectWindow("Temp_C6"); close();
selectWindow("Temp_C7"); close();
roiManager("Show None"); roiManager("Show All");

//-----------------------------------------
// STEP 3: Make a mask of the "black edges"
//         so we know which area to exclude
//-----------------------------------------
// We'll threshold any zero-intensity pixels in "WholeCell_Raw"
// to identify that "dark region." Then we invert to get a mask 
// of the "non-zero" region, so we can track that as well.

selectWindow("WholeCell_Raw");
run("Duplicate...", "title=BlackRegion duplicate");
selectWindow("BlackRegion");
run("8-bit");
setThreshold(0, 0); // specifically pick intensity=0
run("Convert to Mask");
rename("BlackMask"); // white = black(0) pixels in WholeCell_Raw

// Create an expanded mask of the edge
run("Duplicate...", "title=BlackMask_expanded");
selectWindow("BlackMask_expanded");

// Morphologically dilate NPIXELS times
for (d = 0; d < NPIXELS; d++) {
    run("Dilate");
}

// Invert if you want the black region to be black in the mask
// but for intersection checks, we actually want the black region
// as a 'white' object:
selectWindow("BlackMask");
run("Invert");

//-----------------------------------------
// STEP 4: Auto-threshold the "whole cell"
//-----------------------------------------
selectWindow("WholeCell_Raw");

run("8-bit");
//run("Gaussian Blur...", "sigma=2"); // example pre-filter
//run("Sharpen");
waitForUser("Step 4: Whole cell segmentation beginning", 
"You might want to try all thresholding algorithms (shortcut 7) before proceeding.");

selectWindow("WholeCell_Raw");
run(THRESHOLD_STYLE, "method="+WC_AUTO_METHOD+" ignore_black white");
run("Convert to Mask");

//run("Open"); // Erosion + dilation: smoothes objects and removes isolated pixels.
//run("Close-"); // Dilation + erosion: smoothes objects and fills in small holes.
run("Fill Holes");
run("Watershed");
for (d = 0; d < WC_EROSION; d++) {
    run("Erode");
}

// Analyze particles: going over all the objects
run("Analyze Particles...", 
    "size="+MIN_SIZE+"-"+MAX_SIZE+
    " pixel circularity="+MIN_CIRC+"-"+MAX_CIRC+
    " exclude add"); 

var wcCountStart = roiManager("count");
var nWC = 0;

// Remove cells that fall on the border of each well
// Rename good ROIs -> "wc1", "wc2", ...
// Move from later to earlier indices because deletion changes later indices
for (i = wcCountStart - 1; i >= 1; i--) {
	// Step 1: switch to black mask image
	selectWindow("BlackMask_expanded");
	// Step 2: read the median from the Results table
	roiManager("Select", i);
	roiManager("Measure");
	GridMed = getResult("Median", 0);
	run("Clear Results"); // keep table clean
	print("Grid Mask median value: "+GridMed);
	
	if (includeDeathChannel) {
		selectWindow("DeadCellMask");
		roiManager("Select", i);
		roiManager("Measure");
		DeathMean = getResult("Mean", 0);
		run("Clear Results"); // keep table clean
		print("Cell Death Mask median value: "+DeathMean);
	}

	// Step 3: decide if it's inside the grid or death mask
	if (GridMed > 128 || DeathMean > 64) {
		// white = 256, black = 0
		// centroid is on/inside edges or a quarter overlaps with death mask (white) → exclude
		roiManager("Delete");
		print("Skipping ROI index " + i + " due to black-edge overlap.");
	} else {
        // It's valid, so rename it "wcX" and store index
        nWC++;
        name = "wc," + nWC;
        roiManager("rename", name);
        RoiManager.setGroup(1); // whole cell in group 1
	}
}
roiManager("Show None"); roiManager("Show All");

waitForUser("Step 4: Whole cell segmentation completed", 
nWC+" whole cell ROIs passed the filters.\n"+
"Please inspect before clicking 'OK' to segment the nuclei.");
selectWindow("WholeCell_Raw"); close();

//-----------------------------------------
// STEP 5: Auto-threshold the nucleus channel
//-----------------------------------------
selectWindow("AOI_C567");
run("Duplicate...", "title=Nucleus_Raw duplicate channels="+CH_NUCLEUS);
selectWindow("Nucleus_Raw");

run("8-bit");
//run("Gaussian Blur...", "sigma=1.5"); // tweak as needed
waitForUser("Step 5: Nuclear segmentation beginning", 
"You might want to try all thresholding algorithms (shortcut 7) before proceeding.");

selectWindow("Nucleus_Raw");
run(THRESHOLD_STYLE, "method="+NUC_AUTO_METHOD+" ignore_black white");
run("Convert to Mask");

run("Fill Holes"); 
run("Watershed"); // Try not to run watershed because some activated B cell nucleus is cleaved
for (d = 0; d < NUC_EROSION; d++) {
    run("Erode");
}

run("Analyze Particles...", 
    "size="+MIN_SIZE+"-"+MAX_SIZE+
    " pixel circularity=0.00-1.00"+
    " exclude add");

var nucCountStart = roiManager("count");
var nNUC = 0;

// Remove nuclei that fall on the border of each well
// Rename good ROIs -> "nuc1", "nuc2", ...
for (i = nucCountStart - 1; i >= 1+nWC; i--) {
	// Step 1: switch to black mask image
	selectWindow("BlackMask_expanded");
	
	// Step 2: read the median from the Results table
	roiManager("Select", i);
	roiManager("Measure");
	GridMed = getResult("Median", 0);
	run("Clear Results"); // keep table clean
	print("Grid Mask median value: "+GridMed);
	
	if (includeDeathChannel) {
		selectWindow("DeadCellMask");
		roiManager("Select", i);
		roiManager("Measure");
		DeathMean = getResult("Mean", 0);
		run("Clear Results"); // keep table clean
		print("Cell Death Mask median value: "+DeathMean);
	}

	// Step 3: decide if it's inside the mask
	if (GridMed > 128 || DeathMean > 32) {
		// centroid is on/inside black edges → exclude
		roiManager("Delete");
		print("Skipping ROI index " + i + " due to mask overlap.");
	} else {
        // It's valid, so rename it "nucX" and store index
        nNUC++;
        name = "nuc," + nNUC;
        roiManager("rename", name);
        RoiManager.setGroup(2); // nucleus in group 2
	}
}
roiManager("Show None"); roiManager("Show All");
selectWindow("BlackMask"); close();
selectWindow("BlackMask_expanded"); close();

waitForUser("Step 5: Nuclear segmentation completed", 
nNUC+" nuclear ROIs passed the filters. Please inspect before proceeding.");
selectWindow("Nucleus_Raw"); close();

//-----------------------------------------
// STEP 6: Match each nucleus to the best 
//         valid "whole-cell"
//-----------------------------------------
// AOI: 1, nWC, nNUC

var badNUC = 0;
var badNucIndexes = newArray();
var goodWC = 0;

for (j = 0; j < nNUC; j++) {
    nucIndex = 1 + nWC + j;
    roiManager("select", nucIndex);
    nucName = RoiManager.getName(nucIndex);
    roiManager("Measure");
    nucArea = getResult("Area", 0);
    run("Clear Results");

    // We'll find which wc it overlaps most
    bestWC = -1;
    bestOverlap = 0;

    for (k = 0; k < nWC; k++) {
        wcIndex = 1 + k;
        roiManager("select", wcIndex);
        wcName = RoiManager.getName(wcIndex);
        overlapArea = 0;

        // Intersect the two and measure overlap if there is any
        roiManager("Select", newArray(wcIndex,nucIndex));
        roiManager("AND");
        if(selectionType()!=-1){ // if there is overlap
        	roiManager("Add");
        	roiManager("Select", roiManager("count")-1);
        	roiManager("Measure");
        	overlapArea = getResult("Area", 0);
        	run("Clear Results");
        	roiManager("Delete");
        }
        percentOverlap = overlapArea / nucArea; // overlap at least 30% with wc
        if (percentOverlap > 0.3 && overlapArea > bestOverlap) {
            bestOverlap = overlapArea;
            bestWC = wcIndex;
        }
    }
//    bestWCs = Array.concat(bestWCs, bestWC);
//    bestOverlaps = Array.concat(bestOverlaps, bestOverlap);

    if (bestWC >= 0) {
        // Found a match with some overlap, get the wc info
        roiManager("select", bestWC);
        wcName = RoiManager.getName(bestWC);
        if (!startsWith(wcName, "T")) {
        	cellID = substring(wcName, 3); // after 'wc,'
        	newWcName = "Twc," + cellID;
        	roiManager("rename", newWcName);
        	goodWC++;
        } else { // some whole cells are already renamed
        	cellID = substring(wcName, 4); // after 'Twc,'
        }
        
        // Rename nucleus to match cell number
        roiManager("select", nucIndex);
        newNucName = "Tnuc," + cellID; // consistent cell ID with wc
        roiManager("rename", newNucName);
        
    } else {
    	badNUC++;
    	badNucIndexes = Array.concat(badNucIndexes, nucIndex);
        print("No valid wc for nucleus " + nucName + "; skipping.");
    }
}
roiManager("Show None"); roiManager("Show All");

nRois = roiManager("count");
waitForUser("Step 6: Cell pairing completed", 
"There were "+nRois+" ROIs prior to matching wc and nuc pairs. \n"+
"Those without prefix 'T' are unmatched and will be deleted, leaving "+goodWC+" pairs.\n"+
"Please inspect them before clicking 'OK' to delete them.");
// Delete all the nucleus and whole cell with no overlapping partner
nRois = roiManager("count");
for (i = (nRois-1); i > 0; i--) {
	if (i > (roiManager("count")-1)) { // somehow there's an issue if the last ROI was deleted
		print(i+" > total ROI count, skip and continue.");
		continue;
	}
	roiManager("select", i);
	roiName = RoiManager.getName(i);
	if (!startsWith(roiName, "T")) {
		roiManager("select", i);
		roiManager("Delete");
		print("Deleting ROI "+i+", named "+roiName+". "+roiManager("count")+" ROIs left.");
	}
}
roiManager("Show None"); roiManager("Show All");

//-----------------------------------------
// STEP 7: Combine nuclei & compute cytoplasm
//-----------------------------------------
roiManager("Sort");
var sameCell = newArray(); // tracks ROI index
var deleteROI = newArray(); // tracks ROI with repetitive CellID
var goodCells = 0;
var prevCellID = -1;

// find wc that includes multiple nucs & combine them into 1 ROI
for (j = 0; j < (nNUC-badNUC); j++) {
    nucIndex = 1 + j;
    roiManager("select", nucIndex);
    nucName = RoiManager.getName(nucIndex);
    cellID = parseInt(substring(nucName, 5)); // after 'Tnuc,';
    
    if (cellID == prevCellID) { // same ID as the previous nucleus
    	sameCell = Array.concat(sameCell, nucIndex);
    } else { // not the same cell as previous
    	if (sameCell.length > 1) { 
    		// if multiple nucs, combine ROIs from the list before moving on
    		roiManager("Select", sameCell);
    		roiManager("Combine");
			roiManager("Add");
			roiManager("Select", roiManager("count") - 1);
			roiManager("rename", "Tnuc,"+prevCellID);
			RoiManager.setGroup(2); // nucleus in group 2
			deleteROI = Array.concat(deleteROI, sameCell);
			print(sameCell.length+" nucs associated with wc"+prevCellID);
			print("Delete list expanded to "+deleteROI.length);
    	}
    	sameCell = newArray(); // reset the list
    	sameCell = Array.concat(sameCell, nucIndex);
    	goodCells++; // NOTE: goodCells != cellID because some cells were removed
    	prevCellID = cellID;
    }
    if (j == (nNUC-badNUC-1)) { // Edge case: final cell won't have the next cell
    	if (sameCell.length > 1) { 
    		// if multiple nucs, combine ROIs from the list before moving on
    		roiManager("Select", sameCell);
    		roiManager("Combine");
			roiManager("Add");
			roiManager("Select", roiManager("count") - 1);
			roiManager("rename", "Tnuc,"+cellID);
			RoiManager.setGroup(2); // nucleus in group 2
			deleteROI = Array.concat(deleteROI, sameCell);
			print(sameCell.length+" nucs associated with wc"+prevCellID);
			print("Delete list expanded to "+deleteROI.length);
    	}
    }
}
// delete the original nucs
for (i = (deleteROI.length-1); i >= 0; i--) {
	roiManager("Select", deleteROI[i]);
	currRoiName = RoiManager.getName(deleteROI[i]);
	print("Deleting ROI "+deleteROI[i]+", name "+currRoiName);
	roiManager("Delete");
}
roiManager("Show None"); roiManager("Show All");
roiManager("Sort");

waitForUser("Step 7: Nuclei combination completed", 
"Inspect the "+goodCells+" remaining whole cell and nucleus pairs.\n"+
"Click 'OK' to proceed with enlarging whole cell based on associated nucleus.");

if ((2*goodCells+1) != roiManager("count")) {
	waitForUser("Step 7: Warning!", 
	"There is an mismatch in the number of whole cell and nuclei!");
}

goodCells = (roiManager("count")-1)/2;
var veryGoodCells = 0;
// Sometimes whole cell does not include nucleus (e.g. donuts), 
// thus expand wc based on OR gate with nucleus
for (j = 0; j < goodCells; j++) {
	nucIndex = j+1;
    wcIndex = j+1+goodCells;
    
    roiManager("Select", newArray(wcIndex,nucIndex));
    roiManager("OR");
    roiManager("Add");
    
    wcIndex = roiManager("count") - 1;
    roiManager("Select", wcIndex);
    roiManager("Measure");
    wcArea = getResult("Area", 0);
    run("Clear Results");
    roiManager("Select", nucIndex);
    roiManager("Measure");
    nucArea = getResult("Area", 0);
    run("Clear Results");
    if (wcArea > nucArea) { // avoid scenario where nucleus >= whole cell
    	roiManager("Select", wcIndex);
    	newWcName = "Fwc," + veryGoodCells; // final wc, reassign cell ID due to removed cells
    	roiManager("rename", newWcName);
    	RoiManager.setGroup(1); // whole cell in group 1
    
    	roiManager("Select", nucIndex);
    	newNucName = "Fnuc," + veryGoodCells; // final nuc, reassign cell ID due to removed cells
    	roiManager("rename", newNucName); // final nucleus, also reassign ID
    	RoiManager.setGroup(2); // nucleus in group 2
    	
    	veryGoodCells++;
    }
}
roiManager("Show None"); roiManager("Show All");

// Delete all the temporary whole cells & nuclei
for (i = roiManager("count")-1; i > 0; i--) {
	roiManager("select", i);
	roiName = RoiManager.getName(i);
	if (!startsWith(roiName, "F")) {
	roiManager("Delete");
	}
}
selectWindow("AOI_C567");
roiManager("Show None"); roiManager("Show All");

waitForUser("Step 7: Cytoplasm computation beginning", 
"Inspect the "+veryGoodCells+" remaining whole cell and nucleus pairs.\n"+
"[Delete a pair] if they contain multiple cells or clearly dead cell.\n"+
"Shortcut 8 to fit spline, 0 to scale size, and Freehand selection tool to move around.\n"+
"Click 'OK' to proceed with cytoplasm computing.");

// Recalculate veryGoodCells
nROIs = roiManager("count")-1; // excluded the AOI
if ((nROIs % 2) == 0) { // even number
	veryGoodCells = nROIs / 2;
}
else {
	waitForUser("You did not delete wc and nuc ROIs in pairs! Please try again.");
	veryGoodCells = (roiManager("count")-1)/2;
}

// cells should already be in order
// compute cytoplasm
for (j = 0; j < veryGoodCells; j++) {
	nucIndex = j+1;
    wcIndex = j+1+veryGoodCells;
    roiManager("select", wcIndex);
    RoiManager.setGroup(1); // whole cell in group 1
	roiName = RoiManager.getName(wcIndex);
	cellID = parseInt(substring(roiName, 4)); // After 'Fwc,'
	
	roiManager("select", nucIndex);
    RoiManager.setGroup(2); // nucleus in group 2
    
	// Cytoplasm = wc - nucleus
    // XOR with nucleus => difference
    roiManager("Select", newArray(wcIndex,nucIndex));
    roiManager("XOR");
	if(selectionType()!=-1){ // if there is difference
        roiManager("Add");
        cytIndex = roiManager("count") - 1;
     	roiManager("Select", cytIndex);
     	roiManager("rename", "Fcyt,"+cellID);
     	RoiManager.setGroup(3); // cytoplasm in group 3
	} // scenarios where nucleus >= whole cell was avoided, so shouldn't be else cases
}
roiManager("Show None"); roiManager("Show All");

//-----------------------------------------
// STEP 8: Fluorescence intensity measurement
//-----------------------------------------
waitForUser("Step 8: Fluorescence intensity measurement beginning", 
"Please double check all the ROIs.\n"+ 
"!!!! Add an ROI NOW for background measurement !!!!\n"+
"Click 'OK' to proceed with fluorescence intensity measurements");

// Check if background ROI was drawn
roiName = RoiManager.getName(roiManager("count")-1);
if (substring(roiName, 0, 2) != "BG") { // if final ROI does not start with 'BG'
	waitForUser("Warning!!!", 
"Please add an ROI for background measurement with name starting with 'BG'.\n");
}

// Measure statistics of the 3 fluorescence channels in each ROI
selectWindow("AOI_C567");
//selectWindow("Scene_auto2.tif");
roiManager("Show All");
roiManager("multi-measure measure_all");
// remove excessive predixes from the label, leaving only ROI names
for (i=0; i<nResults; i++) {
    oldLabel = getResultLabel(i);
    delimiter = indexOf(oldLabel, ":");
    newLabel = substring(oldLabel, delimiter+1);
    delimiter = indexOf(newLabel, ":");
    newLabel = substring(newLabel, 0, delimiter);
    if (label=="cellID") {
    	delimiter = indexOf(newLabel, ","); // remove the prefix before , (Fnuc, Fwc, Fcyt)
    	newLabel = substring(newLabel, delimiter+1); // since Group contains that info
    }
    setResult("Label", i, newLabel);
  }

resultsFile = homeDir+subDir+"Results_auto"+gridIter+".csv";
if(File.exists(resultsFile)){
	waitForUser("Warning! Are you sure?", 
	"This file exists. Click 'cancel' to abort, or click 'OK' to write to +1 iter.");
	gridIter = gridIter+1;
	resultsFile = homeDir+subDir+"Results_auto"+gridIter+".csv";
}
saveAs("Results", homeDir+subDir+"Results_auto"+gridIter+".csv"); // save results table
roiManager("Save", homeDir+subDir+"RoiSet_auto"+gridIter+".zip"); // save ROI selections
selectWindow("AOI_C567");
saveAs("tiff", homeDir+subDir+"Scene_auto"+gridIter+".tif");
saveAs("Jpeg", homeDir+subDir+"Scene_auto"+gridIter+".jpeg");

if (includeDeathChannel) { // if DRAQ7 was taken
	selectWindow("DeadCellMask"); close();
}
