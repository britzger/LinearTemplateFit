//#include "include/LTF/LTF.h"
#include "LTF/LTF.h"
//#include "LTF_Tools.cxx"
//#include "plot_LTF1D.cxx"
#include "LTF_ROOTTools.h"

#include <TH1D.h>
#include <TGraphErrors.h>

// __________________________________________________________________________________ //
// main
int example_LTF_gaus() {
   using namespace std;

   gStyle->SetOptStat(0);
   gSystem->Load("libLTF.so");
   TH1D::AddDirectory(false);
   TH1::SetDefaultSumw2(true);

   // --- binning of the histograms
   //const vector<double> bins{162, 163, 164, 165, 166,167, 168,169, 170, 171, 172,173, 174, 175, 176, 177, 178};
   const vector<double> bins{169,170,171,172,173, 174, 175, 176, 177, 178,179,180,181,182,183};

   const double sigmadata      = 6.2;
   const double meandata       = 170.2;
   const int nEventsData       = 500;
   const int seeddata          = 4322;

   // ------------------------------------------------ //
   // ---  make templates
   // ------------------------------------------------ //
   map<double,TH1D*> templates;
   { // make templates
      const int nEventsTemplates = 40000;
      const double sigma = 6;
      const vector<double> reference_values{169, 169.5, 170, 170.5, 171, 171.5, 172}; // template reference points
      //const vector<double> reference_values{170, 172}; // template reference points
      int seed = 1234;
      for ( double mean : reference_values ) {
         templates[mean] = LTF_ROOTTools::MakeHistogram(nEventsTemplates, seed++, mean , sigma ,bins);
         templates[mean]->SetTitle(Form("m = %.1f",mean));
         templates[mean]->Scale(double(nEventsData)/nEventsTemplates);
      }
   }
   
   // ------------------------------------------------ //
   // ---  generate pseudo-data
   // ------------------------------------------------ //
   TH1D* data = LTF_ROOTTools::MakeHistogram(nEventsData, seeddata, meandata , sigmadata ,bins);
   data->SetLineColor(kBlack);
   data->SetMarkerSize(1.8);
   


   // ------------------------------------------------ //
   // ---  Do linear template fit
   // ------------------------------------------------ //
   // --- instantiate LTF object
   LTF ltf;
   ltf.SetGamma(vector<double>{1});
   ltf.UseNuisanceParameters(false);
   ltf.UseLogNormalUncertainties(false);

   // --- initialize templates
   for ( auto [MM,hist] : templates ) {
      ltf.AddTemplate(MM,  bins.size()-1,  hist->GetArray()+1 ); // set template
      ltf.AddTemplateErrorSquared("statY", MM , bins.size()-1, hist->GetSumw2()->GetArray()+1, 0.); // set template error dY
   }
   
   // --- initialize data
   ltf.SetData( bins.size()-1,data->GetArray()+1);
   ltf.AddUncorrelatedErrorSquared("stat.", bins.size()-1, data->GetSumw2()->GetArray()+1);
   // for ( const auto& s : shiftsnuisance ) ltf.AddError("",N,s->GetArray()+1,1.);
   // for ( const auto& s : shifts ) ltf.AddError("",N,s->GetArray()+1,1.);
   // for ( const auto& s : shifts ) ltf.AddError("",N,s->GetArray()+1,0.5,LTF::Uncertainty::External);
   

   LTF::LiTeFit fit = ltf.DoLiTeFit();
   fit.PrintFull();
   //fit.DoIterativeFitNewton(6,0.6,2,1);
   //fit.DoIterativeFitTaylor();
   //fit.PrintFull();
   LTF_ROOTTools::plotLiTeFit(fit,bins);
   
   return 0;

}


int main() {
   return example_LTF_gaus();
}
