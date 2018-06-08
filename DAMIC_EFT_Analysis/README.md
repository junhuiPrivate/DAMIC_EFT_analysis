Author: Jun Liao (junhui.private@gmail.com)
This is a set of example script showing setting 90% CL upper limit for 10 GeV/c^2 WIMP with EFT O1.
The same script framework can be used for setting limits for other WIMP masses and operators.
The scripts are based on the tutorial twiki page: https://twiki.cern.ch/twiki/bin/view/RooStats/RooStatsTutorialsJune2013. 
CountingModel.C describes background and signal models. Running this script can generate CountingModel.root which will be used in limit setting script, SimpleHypoTestInv.C. The scripts, BkgLinearPdf.cxx and OOne10GeVPdf.cxx, are the background and signal models, respectively.

To run it:
(1). Running this command to generate CountingModel.root, "rl CountingModel.C".
(2). Running this command to get the scanned p-values of POI and detailed statistical results, "rl SimpleHypoTestInv.C".
