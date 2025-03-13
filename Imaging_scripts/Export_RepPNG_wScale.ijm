var homeDir = "/Users/helenhuang/Documents/1st Year PhD/Hoffmann Lab/Microscopy/VK014/";
var subDir = "VK014_CD40Low_7h/"; // CHANGE!!!!! to condition folder
// BCRHigh

timepoint = 1; //0: 0h, 1: 7h, 2: 24h, 3: 48h, 4: 72h
//tif_name = "Scene_auto5.tif";
tif_name = "Composite";

cRelmins = newArray(952, 924, 921, 477, 515);
cRelmaxes = newArray(3125, 2924, 4489, 1178, 1543);

H2Bmins = newArray(559, 540, 601, 424, 439);
H2Bmaxes = newArray(6641, 5088, 5437, 1060, 1236);

RelAmins = newArray(976, 908, 937, 408, 410);
RelAmaxes = newArray(6683, 7668, 9417, 528, 572);

run("Set Scale...", "distance=100 known=12.897 unit=μm");
selectImage(tif_name);
run("Duplicate...", "title=AOI duplicate");
selectImage("AOI");
run("Scale Bar...", "width=10 height=10 location=[Lower Left] horizontal bold overlay");
//run("Scale Bar...", "width=10 height=10 horizontal bold overlay");
run("Split Channels");

selectImage("C1-AOI");
setMinAndMax(cRelmins[timepoint], cRelmaxes[timepoint]);

selectImage("C2-AOI");
setMinAndMax(H2Bmins[timepoint], H2Bmaxes[timepoint]);

selectImage("C3-AOI");
setMinAndMax(RelAmins[timepoint], RelAmaxes[timepoint]);

selectImage("C1-AOI"); saveAs("PNG", homeDir+subDir+"C1-AOI3.png"); close();
selectImage("C2-AOI"); saveAs("PNG", homeDir+subDir+"C2-AOI3.png"); close();
selectImage("C3-AOI"); saveAs("PNG", homeDir+subDir+"C3-AOI3.png"); close();