/*
  Copyright (c) 2021, D. Britzger, Max-Planck-Institute for Physics, Munich, Germany

  Permission is hereby granted, free of charge, to any person obtaining
  a copy of this software and associated documentation files (the
  "Software"), to deal in the Software without restriction, including
  without limitation the rights to use, copy, modify, merge, publish,
  distribute, sublicense, and/or sell copies of the Software, and to
  permit persons to whom the Software is furnished to do so, subject to
  the following conditions:

  The above copyright notice and this permission notice shall be
  included in all copies or substantial portions of the Software.

  The Software is provided "as is", without warranty of any kind,
  express or implied, including but not limited to the warranties of
  merchantability, fitness for a particular purpose and
  noninfringement. In no event shall the authors or copyright holders be
  liable for any claim, damages or other liability, whether in an action
  of contract, tort or otherwise, arising from, out of or in connection
  with the Software or the use or other dealings in the Software.
*/
// -------------------------------------------------------------------- //

#include <iostream>

#if defined __WITH_ROOT__ || defined __CLING__

#include "LTF/LTF.h"
#include "LTF/LTF_ROOTTools.h"
#include <TSystem.h>

void PrintAsciiTable(const map<double,TH1D*>&, TH1D* data);

// __________________________________________________________________________________ //
// main
int example1_LTF_gaus() {
   using namespace std;

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
      const int nEventsTemplates = 40000;
      const double sigma = 6;
      const vector<double> reference_values{169, 169.5, 170, 170.5, 171, 171.5, 172}; // template reference points
   { // make templates
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

   
   PrintAsciiTable(templates,data);

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


//! ------------------------------------------------------------------------ //
//! --- write templates, and data, to ascii file
void PrintAsciiTable(const map<double,TH1D*>& templates, TH1D* data){
   cout<<endl;
   printf(" %11s %11s",Form("Data"),Form("Stat"));
   for ( auto [mean,hist] : templates ) 
      printf(" %11s %11s",Form("Tmpl_%5.2f",mean),Form("Stat_%5.2f",mean));
   cout<<endl;
   for ( int i=1; i<=data->GetNbinsX() ;i++ ) {
      printf(" %11.4f %11.4f",data->GetBinContent(i),data->GetBinError(i));
   for ( auto [mean,hist] : templates ) 
         printf(" %11.4f %11.4f",hist->GetBinContent(i),hist->GetBinError(i));
      cout<<endl;
   }
   cout<<endl;
}

#endif // end WITH_ROOT



//! ------------------------------------------------------------------------ //
//! main function
int main() {
#if defined __WITH_ROOT__ || defined __CLING__
   return example1_LTF_gaus();
#else
   //! printout if compiled without ROOT
   std::cout<<"This example is working only if ROOT is available. Otherwise, please see 'example1_LTF_gaus_NoROOT'."<<std::endl;
   return 0;
#endif
}
