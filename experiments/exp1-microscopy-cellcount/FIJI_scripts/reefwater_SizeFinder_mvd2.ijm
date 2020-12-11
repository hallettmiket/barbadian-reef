filePath = File.openDialog("Locate your file:");
print(filePath);
//C:\Users\Admin\Desktop\Shawn_reefwater_20200723\Shawn_reefwater_20200723.mvd2

fileArray = split(filePath, "\\");
dirSpot = "C:";
for (x = 1;x<lengthOf(fileArray)-1;x++){
	dirSpot = dirSpot + "\\" + fileArray[x];
}
print(dirSpot);
//fList=getFileList(fPath);
//fPref=split(fList[0],"."); //these two lines allow the renaming of output files
	//Array.print(fPref);
	
File.makeDirectory(dirSpot+"/output/");
run("Clear Results"); 
updateResults();
run("Bio-Formats Macro Extensions");
run("Set Measurements...", "area mean shape redirect=None decimal=3");
setBackgroundColor(0,0,0);
setForegroundColor(255,255,255);

//run("Bio-Formats", "open="+filePath+" autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");

Ext.setId(filePath)
Ext.getSeriesCount(seriesCount);
//close();

for(f = 0; f<seriesCount;f++){
		run("Clear Results"); 
	updateResults();
	run("Bio-Formats", "open="+filePath+" autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+f+1);
	roiManager("Reset");
	fName = getTitle();
	fName = split(fName , "-");
	fName= fName[1];
	print(fName);
	//fName = fileArray[lengthOf(fileArray)-1];
	//fName=substring(fName,0,lengthOf(fName)-5);
	run("Set Scale...", "distance=1 known=1 pixel=0.102 unit=micron");
	run("Z Project...", "projection=[Max Intensity]");
	run("Duplicate...", " ");
	run("Subtract Background...", "rolling=50");
	run("Green");
	run("Gaussian Blur...", "sigma=1");
	setAutoThreshold("Triangle dark");
	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("Analyze Particles...", "add");
	close();
	roiManager("Show None");
	roiManager("Show All");
	waitForUser("Check it Out!");
	roiManager("multi-measure append");
	saveAs("Results",dirSpot+"/output/"+fName+"_results.csv");

	run("Close All");
}


