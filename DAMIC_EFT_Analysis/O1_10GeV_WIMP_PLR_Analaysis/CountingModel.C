/// This script has been modified to read data from a tree/histogram in a separate root file.

#include "TTree.h"
#include "TH1.h"

using namespace RooFit;
using namespace RooStats;

void CountingModel(  int nobs = 53,           // number of observed events
                     double n_bkg_pdf = 50,           // number of background events
                     double sigmab = 7 )   // relative uncertainty in b
{
   RooWorkspace w("w");

  RooClassFactory::makePdf("OOne10GeVPdf","cof") ;
#ifdef __CINT__
  gROOT->ProcessLineSync(".! cp BackupOOne10GeVPdf.cxx OOne10GeVPdf.cxx") ;
  gROOT->ProcessLineSync(".! cp BackupOOne10GeVPdf.h OOne10GeVPdf.h") ;
  gROOT->ProcessLineSync(".x OOne10GeVPdf.cxx+") ;
  gROOT->ProcessLineSync(".L OOne10GeVPdf.cxx+") ;
#endif

   w.factory("OOne10GeVPdf:n_sig_pdf(cof[5.0E-2,5.0E-3,0.9E-1])");
   w.factory("sum:nexp(n_sig_pdf, n_bkg_pdf[53,0,100])");
///Use the # of events to build Poisson distribution for signal and Gaussian for background
// Poisson of (n | s+b)
   w.factory("Poisson:pdf(nobs[0,100],nexp)");///don't change the maxiumum value of 100 to an arbitary big one.
   w.factory("Gaussian:constraint(b0[0,100],n_bkg_pdf,sigmab[8])"); /// 7 ~= sqrt(53).
// make Poisson model * Gaussian constraint
   w.factory("PROD:model(pdf,constraint)");

   
   ///~~~~~~~  Begining of creating the ModelConfig object and set "b0" as a global observable ~~~///

   w.var("b0")->setVal(n_bkg_pdf);   
   w.var("b0")->setConstant(true); // needed for being treated as global observables
   w.var("sigmab")->setVal(sigmab); 
   
   ModelConfig mc("ModelConfig",&w);/// declare the name and reference of ModelConfig as "mc" and "ModelConfig".
   mc.SetPdf(*w.pdf("model")); 
   mc.SetParametersOfInterest(*w.var("cof")); 
   ///mc.SetParametersOfInterest(*w.var("n_sig_pdf")); 
   mc.SetObservables(*w.var("nobs"));  
   ///w.defineSet("nuisanceParams","a,b,c");
   ///mc.SetNuisanceParameters(*w.var("nuisanceParams"));
   mc.SetNuisanceParameters(*w.var("n_bkg_pdf"));

   // these are needed for the hypothesis tests
   mc.SetSnapshot(*w.var("cof"));
   // need now to set the global observable
   mc.SetGlobalObservables(*w.var("b0"));

   mc.Print();
   // import model in the workspace 
   w.import(mc);



   ///~~~~~~~  End of creating the ModelConfig object and set "b0" as a global observable ~~~///

  // import a tree
  TFile *treefile = TFile::Open("total54PRD2016PaperEvents.root"); /// total 53 events.
  TTree* tree = (TTree*) treefile->Get("ntuple");

  /// make dataset with the number of ovserved events.
  RooDataSet ds("ds","ds",*w.var("nobs"),Import(*tree)) ; ///refered from rf102_dataimport.C
  w.var("nobs")->setVal(53);/// for 53 observed events in the tree
  ds.add(*w.var("nobs") );
 
   /// Draw the histogram imported.
   /// Above three lines should be commented before uncommenting following lines to plot the data.
/*  
  RooRealVar x("x","x",1,0.21,10) ; 
  RooDataSet ds("ds","ds",*w.var("x"),Import(*tree)) ;
  ds.add( *w.var("x"));
  RooPlot* frame3 = x.frame(Title("Recorded energy in CCDs")) ;
  ds.plotOn(frame3) ;
  mc.GetPdf().plotOn(frame3) ;

  TCanvas* c = new TCanvas("Recorded energy","Recorded energy",800,800) ;
   gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.4) ; frame3->Draw() ;
*/ 


  w.import(ds);
   w.Print();

   TString fileName = "CountingModel.root"; 

   // write workspace in the file (recreate file if already existing)
   w.writeToFile(fileName, true);

}
