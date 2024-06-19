// Z Projections for Trypanosoma Colonies
// former file name "Z-Projektionen für Trypanosomen-1 eng.ijm"
// 24-06-17
// Magdalena Schüttler




if (nImages != 0) {
	waitForUser("Z-project colonies", "The images that are currently opened will be closed. \nIf necessary, please save them now.");
	close("*");
	waitForUser("Z-project colonies", "Please open all segmentations that should be \ncombined into the Z projection and click \"OK\".");
} else {
	waitForUser("Z-project colonies", "Please open all segmentations that should be \ncombined into the Z projection and click \"OK\".");
}



n = nImages;

for (i = 0; i < n; i++) {
	title = getTitle();
	selectWindow(title);
	run("Z Project...", "projection=[Average Intensity]");
	zDir = File.getParent(getDirectory("file")) + "\\z-projections\\";
	if (!File.exists(zDir)) {
		File.makeDirectory(zDir);
	}
	title2 = File.getNameWithoutExtension(title);
	saveAs("Tiff", zDir + title2 + "-Z-projection");
	titleZ = getTitle();
	close(title);
	close(titleZ);
}