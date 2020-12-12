fPath=getDirectory("Choose a Directory");
fList=getFileList(fPath);
//fPref=split(fList[0],"."); //these two lines allow the renaming of output files
	//Array.print(fPref);
	
File.makeDirectory(fPath+"output/");
run("Clear Results"); 
updateResults();
run("Bio-Formats Macro Extensions");
run("Set Measurements...", "area mean shape redirect=None decimal=3");
setBackgroundColor(0,0,0);
setForegroundColor(255,255,255);

for(f = 0; f<lengthOf(fList);f++){
	if (fList[f] != "output/"){
		run("Bio-Formats", "open="+fPath+fList[f]+" autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		roiManager("Reset");
		fName = fList[f];
		fName=substring(fName,0,lengthOf(fName)-4);
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
		saveAs("Results",fPath+"output/"+fName+"_results.csv");
		run("Clear Results"); 
		updateResults();
		run("Close All");
	}
}

