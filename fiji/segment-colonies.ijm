// WEKA Segmentation of Trypanosoma Colonies
// former file name "Trypanosomen Segmentierung WEKA-4 eng.ijm"
// 24-06-17
// Magdalena Schüttler




relative = getDirectory("file");		// folder, that contains the script file	// allows for relative pathing

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// -> the default entries in the dialog box can be changed here:
// file paths that include ä, ö, ü can cause errors during saving, these should be avoided

// -> path to the Input file:
defaultinput = relative + "sample-image.tif";
// -> path to the classifier file:
defaultclassifier = relative + "classifier.model";
// -> path to the Output folder:
defaultoutput = relative + "output\\";

// -> sigma parameter for Gaussian Blur (does not need to be changed usually):
sigma = 1.5;
//////////////////////////////////////////////////////////////////////////////////////////////////////////





// DIALOG WINDOW
Dialog.createNonBlocking("WEKA Segmentation");
Dialog.setInsets(0, 0, 0);
Dialog.addMessage("Segmentation of Trypanosoma colonies using Trainable WEKA Segmentation");

Dialog.setInsets(15, 0, 0);
Dialog.addFile("Choose Input file:", defaultinput);
Dialog.setInsets(0, 60, 0);
Dialog.addMessage("Input file is a non-cropped raw image stack in Tiff format.");

Dialog.setInsets(15, 0, 0);
Dialog.addNumber("Number of horizontal plates:", 4);							// nhor
Dialog.addNumber("Number of vertical plates:", 2);								// nvert

Dialog.setInsets(15, 0, 0);
methodsArray = newArray("automatically", "manually");
Dialog.addChoice("Crop method:", methodsArray, methodsArray[0]);				// method
Dialog.setInsets(0, 0, 0);
regtimeArray = newArray("after Cropping", "after Segmentation");
Dialog.addChoice("Registration order", regtimeArray);							// regtime

Dialog.setInsets(15, 0, 0);
Dialog.addFile("Choose classifier file:", defaultclassifier);

Dialog.setInsets(15, 0, 0);
Dialog.addDirectory("Choose empty Output folder:", defaultoutput);

Dialog.show()

// get file paths and create subfolders
inputFile = Dialog.getString();
classifier = Dialog.getString();
outputDir = Dialog.getString();

if (!File.exists(outputDir)) {
	File.makeDirectory(outputDir);
}

resultsDir = outputDir + "results\\";
if (!File.exists(resultsDir)){
	File.makeDirectory(resultsDir);
}

tempDir = outputDir + "temp\\";
if (!File.exists(tempDir)){
	File.makeDirectory(tempDir);
}

cropDir = tempDir + "crop\\";
if (!File.exists(cropDir)) {
	File.makeDirectory(cropDir);
}

regDir = tempDir + "reg\\";
if (!File.exists(regDir)) {
	File.makeDirectory(regDir);
}

splitDir = tempDir + "split\\";
if (!File.exists(splitDir)) {
	File.makeDirectory(splitDir);
}

splitSegDir = tempDir + "split seg\\";
if (!File.exists(splitSegDir)) {
	File.makeDirectory(splitSegDir);
}

// get variables
nhor = Dialog.getNumber();
nvert = Dialog.getNumber();
ncol = nhor*nvert;
method = Dialog.getChoice();
regtime = Dialog.getChoice();
checkZ = Dialog.getCheckbox();




// PREPROCESSING
// open stack
open(inputFile);
title = getTitle();
print("\\Clear");
print("Input file: ", title);
print("Number of colonies: ", ncol);

if(bitDepth()!=8) {
	run("8-bit");
}

Stack.getDimensions(sw, sh, sc, ss, sf);

// reverse stack
if (ss == 3) {
	Stack.swap(1, 3);
}
if (ss == 4) {
	Stack.swap(1, 4);
	Stack.swap(2, 3);
}
if (ss != 3 && ss != 4) {
	waitForUser("Error", "The current stack could not be reversed automatically. Please reverse it manually using the Stack Sorter: Image>Stacks>Tools>Stack Sorter>Reverse \nThen click \"OK\".");
}




if (method == "automatically") {
	// define dimensions for cropping
	selw = sw/nhor;
	selh = sh/nvert;
	selhred = selh*0.8; // needed to cut part of the edge of the plate for its removal
	
	print(" ");
	print("Cropping colonies...");
	
	// crop and save all colonies
	colcount = 1;
	for (i=0; i<nvert; i++) {
		for (j=0; j<nhor; j++){
			run("Specify...", "width="+selw+" height="+selh+" x="+(j*selw)+" y="+(i*selh)+" slice=1");
			run("Duplicate...", "duplicate");
			titleCrop = getTitle();
			saveAs("Tiff", cropDir + "Colony-"+colcount);
			close();
			colcount++;
		}
	}
	close();
	print("\\Update:Cropping completed.");
}




if (method == "manually") {	
	print(" ");
	print("Cropping colonies manually...");
	
	// crop and save all colonies
	colcount = 1;
	for (i=0; i<nvert; i++) {
		for (j=0; j<nhor; j++){
			setTool("rectangle");
			waitForUser("Crop Colonies", "Please mark colony " + colcount + ". \nThe edge of the plate should not be included.");
			
			run("Duplicate...", "duplicate");
			titleCrop = getTitle();
			saveAs("Tiff", cropDir + "Colony-"+colcount);
			close();
			colcount++;
		}
	}
	close();
	print("\\Update:Cropping completed.");
}

if (regtime == "after Cropping") {	
	// REGISTRATION with StackReg and PREPROCESSING after cropping
	print("Registration, splitting and preprocessing of colonies...");
	list = getFileList(cropDir);
	colcount = 1;
	for (i=0; i<ncol; i++) {
		open(cropDir + "/" + list[i]);
		
		// register colonies
		setOption("ScaleConversions", true);
		run("StackReg", "transformation=[Rigid Body]");
		saveAs("Tiff", regDir + "/Colony-"+colcount+ " reg");
		
		// split stack and apply Gaussian Blur
		for (j=0; j<ss; j++) {
			setSlice((j+1));
			run("Duplicate...", "use");
			run("Gaussian Blur...", "sigma="+sigma);
			saveAs("Tiff", splitDir + "Colony-"+colcount+"-"+(j+1));
			close();
		}
		colcount++;
		close();
	}
	print("\\Update:Registration, splitting and preprocessing completed.");
} else {
	// PREPROCESSING
	print("Splitting and preprocessing of colonies...");
	list = getFileList(cropDir);
	colcount = 1;
	for (i=0; i<ncol; i++) {
		open(cropDir + "/" + list[i]);
			
		// split stack and apply Gaussian Blur
		for (j=0; j<ss; j++) {
			setSlice((j+1));
			run("Duplicate...", "use");
			run("Gaussian Blur...", "sigma=1.5");
			saveAs("Tiff", splitDir + "Colony-"+colcount+"-"+(j+1));
			close();
		}
		colcount++;
		close();
	}
}




// SEGMENTATION with Trainable WEKA Segmentation and POSTPROCESSING
print(" ");
print("Log-Output of Trainable WEKA Segmentation: -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  ");

list = getFileList(splitDir);
open(splitDir + "/" + list[0]); // open sample image, does not matter which
run("Trainable Weka Segmentation");
while(getInfo("window.title") != "Trainable Weka Segmentation v3.3.4") {
	wait(100);
}
selectWindow("Trainable Weka Segmentation v3.3.4");
call("trainableSegmentation.Weka_Segmentation.loadClassifier", classifier)

for (i=0; i<list.length; i++) {
	// segment colonies
	selectWindow("Trainable Weka Segmentation v3.3.4");
	close("\\Others");
	selectWindow("Trainable Weka Segmentation v3.3.4");
	call("trainableSegmentation.Weka_Segmentation.applyClassifier", splitDir, "/"+list[i], "showResults=true", "storeResults=false", "probabilityMaps=false", ""); 
	
	// postprocessing
	selectWindow(list[i]);
	titleOld = File.nameWithoutExtension();
	close();
	
	selectWindow("Classification result");
	run("Duplicate...", "use");
	titleDup = getTitle();
	
	selectWindow("Classification result");
	run("Close");
	
	selectWindow(titleDup);
	
	setThreshold(0, 0, "raw"); 
	setOption("BlackBackground", false);
	run("Convert to Mask");
	
	saveAs("Tiff", splitSegDir + titleOld + " seg");
}

close("*");
print(" -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  ");
print(" ");
print("Segmentation and postprocessing is completed.");




// RECOMBINING colonies into stacks and POSTPROCESSING
print("Recombining colonies...");
list = getFileList(splitSegDir);
filecounter = 0;
zprojerror = 0;
for (i=0; i<ncol; i++) {
	// recombine colonies
	for (j=0; j<ss; j++) {
		open(splitSegDir + list[filecounter]);
		filecounter++;
		
		// postprocessing
		// removes edge of plate for automatic cropping
		if (method == "automatically") { 			
			// remove largest object, which is usually the edge of the plate
			
			t = getTitle();
			run("Remove Largest Region");
			close(t);
		}
		
		t = getTitle();
		run("Analyze Particles...", "  show=Masks exclude");
		close(t);
		t = getTitle();
		run("Keep Largest Region");
		close(t);
		t = getTitle();
		run("Fill Holes (Binary/Gray)");
		close(t);
	}
	
	run("Images to Stack", "name=[Colony-" + (i+1) + " seg] use");
	title = getTitle();
	
	
	if (regtime == "after Segmentation") {	
		// REGISTRATION with StackReg
		print("Registration of colonies...");
		list = getFileList(splitSegDir);
		colcount = 1;
		for (k=0; k<ncol; k++) {
			open(splitSegDir + list[i]);
			
			// register colonies
			setOption("ScaleConversions", true);
			run("StackReg", "transformation=[Rigid Body]");
			saveAs("Tiff", regDir + "Colony-"+colcount+ " reg");
			
			colcount++;
			close();
		}
		print("\\Update:Registration completed.");
	}
	
	saveAs("Tiff", resultsDir + title);
	close("*");
}
print("\\Update:Recombining colonies completed.");

print(" ");
print("Application of the script is completed..");