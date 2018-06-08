using namespace RooStats;
using namespace RooFit;
#include "TMath.h"


/// plot histograms for all detailed outputs
bool PlotAllDetailedOutput(HypoTestResult * r, double poiValue, bool nullOutput) {
   RooDataSet * data = 0; 
   if (nullOutput) 
      data = r->GetNullDetailedOutput();
   else 
      data = r->GetAltDetailedOutput();

   if (!data) { 
      Error("PlotDetailedOutput","Detailed output is not present in result file");
      return false;
   }
   
   ///loop on the data variables 

   const RooArgSet * vars = data->get(); 

   TCanvas * c2 = new TCanvas(); 
   c2->SetTitle(TString::Format("Output for poi = %f",poiValue) );
   int n = vars->getSize(); 
   int ny = TMath::CeilNint(TMath::Sqrt(n));
   int nx = TMath::CeilNint(double(n)/ny);
   c2->Divide( nx,ny);

   TFile * f = new TFile("detailedOutputHistograms.root","UPDATE");
   
   TString type = (nullOutput) ? "null" : "alt"; 
   
   RooArgList lvars(*vars);
   for (int i = 0; i < n; ++i) { 
      RooAbsArg & arg = lvars[i]; 
      TString hname = TString::Format("%s/poi_%f/h_%s",type.Data(),poiValue,arg.GetName() ); 
      RooAbsRealLValue * var = dynamic_cast<RooAbsRealLValue*>(&arg);
      TH1 * h1 = new TH1D(hname,hname,100,1,0); 
      for (int j = 0; j < data->numEntries(); ++j) {
         const RooArgSet * ev = data->get(j); 
         h1->Fill(ev->getRealValue(arg.GetName() ) );
      }
      c2->cd(i+1);
      h1->Draw();
      h1->Write();
   }
   gPad->Update();

   ///get also the fit information
   RooAbsData * fitData = r->GetFitInfo(); 
   if (!fitData) return false; 
   const RooArgSet * fitvars = fitData->get(0);

  ///loop on all variable with null or all 

   RooArgList fitVarList(*fitvars);
   if (fitVarList.getSize() < 1) return false; 
   
   ///make a label histogram with the fit data information
   TString hname1 = TString::Format("fitData_%s_poi_%f",type.Data(),poiValue ); 
   TH1D * hf = new TH1D(hname1,hname1,fitvars->getSize()/2,0,1);

   for (int i = 0; i < fitVarList.getSize(); ++i) { 
      RooAbsArg & arg = fitVarList[i];
      TString name = arg.GetName(); 
      if (nullOutput && name.Contains("fitNull") ) { 
         name.ReplaceAll("fitNull_","");
      }
      if (nullOutput && name.Contains("fitAlt") ) { 
         name.ReplaceAll("fitAlt_","");
      }
      double value = ( (RooRealVar &) arg).getVal();
      double error = ( (RooRealVar &) arg).getError();
      int bin = hf->Fill(name,value);
      hf->SetBinError(bin,error);
      hf->SetMinimum(-3); hf->SetMaximum(3);         
   }
   hf->Write();


   f->Close();
   return true;


}





void SimpleHypoTestInv( const char* infile =  "CountingModel.root", 
                        const char* workspaceName = "w",
                        const char* modelConfigName = "ModelConfig",
                        ///const char* bmodelConfigName = "BModelConfig",
                        const char* dataName = "ds" )
{

   ///ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(3);
   ///RooStats::AsymptoticCalculator::SetPrintLevel(1);  

  /////////////////////////////////////////////////////////////
  // First part is just to access the workspace file 
  ////////////////////////////////////////////////////////////
  RooClassFactory::makePdf("OOne10GeVPdf","x,cof") ;
#ifdef __CINT__
  gROOT->ProcessLineSync(".! cp BackupOOne10GeVPdf.cxx OOne10GeVPdf.cxx") ;
  gROOT->ProcessLineSync(".! cp BackupOOne10GeVPdf.h OOne10GeVPdf.h") ;
  gROOT->ProcessLineSync(".x OOne10GeVPdf.cxx+") ;
  gROOT->ProcessLineSync(".L OOne10GeVPdf.cxx+") ;
#endif

  // open input file 
  TFile *file = TFile::Open(infile);
  if (!file) return;

  // get the workspace out of the file
  RooWorkspace* w = (RooWorkspace*) file->Get(workspaceName);


  // get the modelConfig out of the file
  RooStats::ModelConfig* mc = (RooStats::ModelConfig*) w->obj(modelConfigName);
  ///RooStats::ModelConfig* bmc = (RooStats::ModelConfig*) w->obj(bmodelConfigName);

  // get the data out of the file
  RooAbsData* data = w->data(dataName);

  ModelConfig*  sbModel = (RooStats::ModelConfig*) w->obj(modelConfigName);
  RooRealVar* poi = (RooRealVar*) sbModel->GetParametersOfInterest()->first();
  /// background model created from s+b model by setting nsig = 0.
  ModelConfig * bModel = (ModelConfig*) sbModel->Clone();
  bModel->SetName(TString(sbModel->GetName())+TString("_with_poi_0"));      
  poi->setVal(0);
  bModel->SetSnapshot( *poi  );
  
  FrequentistCalculator  fc(*data, *bModel, *sbModel);
  ///fc.SetToys(3000,1500);  
  fc.SetToys(1000,1000);  
  ///fc.SetToys(100,50);  
  ///fc.SetToys(5,5);  

  // asymptotic calculator
   ///AsymptoticCalculator  ac(*data, *bModel, *sbModel);
   ///ac.SetOneSided(true);  // for one-side tests (limits)
  //  ac->SetQTilde(true);
  /// AsymptoticCalculator::SetPrintLevel(-1);


  // create hypotest inverter 
  // passing the desired calculator 
  ///HypoTestInverter calc(ac);    // for asymptotic 
  HypoTestInverter calc(fc);  // for frequentist

  // set confidence level (e.g. 95% upper limits)
  calc.SetConfidenceLevel(0.90);
  
  // for CLS
  bool useCLs = true;
  calc.UseCLs(useCLs);
  calc.SetVerbose(false);
  ///calc.SetVerbose(2);

  // configure ToyMC Samler (needed only for frequentit calculator)
  ToyMCSampler *toymcs = (ToyMCSampler*)calc.GetHypoTestCalculator()->GetTestStatSampler();
  
  // profile likelihood test statistics 
  ProfileLikelihoodTestStat profll(*sbModel->GetPdf());
  // for CLs (bounded intervals) use one-sided profile likelihood
  if (useCLs) profll.SetOneSided(true);

  // ratio of profile likelihood - need to pass snapshot for the alt 
  // RatioOfProfiledLikelihoodsTestStat ropl(*sbModel->GetPdf(), *bModel->GetPdf(), bModel->GetSnapshot());
   
  // set the test statistic to use 
  toymcs->SetTestStatistic(&profll);

  // if the pdf is not extended (e.g. in the Poisson model) 
  // we need to set the number of events
  if (!sbModel->GetPdf()->canBeExtended())
     toymcs->SetNEventsPerToy(1);


  int npoints = 10;  // for a tree, it's the number of bins.
  ///int npoints = 565;  // for a histogram, it's number of entries or bins ?
  // min and max (better to choose smaller intervals)
  double poimin = poi->getMin();
  double poimax = poi->getMax();
  ///double poimin = 0., poimax = 10.; /// poimin is the min of x to be scanned, similarly to poimax.

  std::cout << "Doing a fixed scan  in interval : " << poimin << " , " << poimax << std::endl;
  calc.SetFixedScan(npoints,poimin,poimax);
  ///calc.SetAutoScan();
 
  ///to enable detailed ouput and should be put before "calc.GetInterval()".
  profll.EnableDetailedOutput();

  ///std::cout << std::endl << "test flag 6" << std::endl;
  // run inverter, only one line, :-).
  HypoTestInverterResult * r = calc.GetInterval();
  ///std::cout << std::endl << "test flag 7" << std::endl;

  ///sbModel.GetPdf().fitTo(data); 
  double upperLimit = r->UpperLimit();

  std::cout << "The computed upper limit is: " << upperLimit << std::endl;
  
  // compute expected limit
  std::cout << "Expected upper limits, using the B (alternate) model : " << std::endl;
  std::cout << " expected limit (median) " << r->GetExpectedUpperLimit(0) << std::endl;
  std::cout << " expected limit (-1 sig) " << r->GetExpectedUpperLimit(-1) << std::endl;
  std::cout << " expected limit (+1 sig) " << r->GetExpectedUpperLimit(1) << std::endl;
  

  // plot now the result of the scan 

  HypoTestInverterPlot *plot = new HypoTestInverterPlot("HTI_Result_Plot","HypoTest Scan Result",r);

  // plot in a new canvas with style
  TCanvas * c1 = new TCanvas("HypoTestInverter Scan"); 
  c1->SetLogy(false);

  plot->Draw("CLb 2CL");  // plot also CLb and CLs+b 
  ///plot->Draw("CLb");  // plot also CLb and CLs+b 
  //plot->Draw("OBS");  *// plot only observed p-value


  // plot also in a new canvas the test statistics distributions 
  
  // plot test statistics distributions for the two hypothesis
  // when distribution is generated (case of FrequentistCalculators)

///  /*
  const int n = r->ArraySize();
  if (n> 0 &&  r->GetResult(0)->GetNullDistribution() ) { 
     TCanvas * c2 = new TCanvas("Test Statistic Distributions","",2);
     if (n > 1) {
        int ny = TMath::CeilNint( sqrt(n) );
        int nx = TMath::CeilNint(double(n)/ny);
        c2->Divide( nx,ny);
     }
     for (int i=0; i<n; i++) {
        if (n > 1) c2->cd(i+1);
        SamplingDistPlot * pl = plot->MakeTestStatPlot(i);
        pl->SetLogYaxis(true);
        pl->Draw();
     }
  }

  ///loop on the points to save test parameters into a file.
   for (int i = 0; i < r->ArraySize(); ++i) {

      HypoTestResult * r0 = r->GetResult(i); 

      bool ok = false; 
      ///bool nullOutput = true; 
      bool nullOutput = false; 
      ok = PlotAllDetailedOutput(r0, r->GetXValue(i),nullOutput);          

      if (!ok) return;

   }

///  */
}
