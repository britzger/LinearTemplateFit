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
#include <TH1D.h>
#include <TSystem.h>
#include <string>

using namespace std;


void PrintAsciiTable(const map<vector<double>,TH1D*>& templates, TH1D* data);


// __________________________________________________________________________________ //
//! main
int example2_LTF_gaus2D() {

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
   map<vector<double>,TH1D*> templates;

   { // make templates
      const int nEventsTemplates = 40000;
      //const double sigma = 6;
      const vector<double> reference_values1{169.5, 170, 170.5, 171}; // template reference points
      const vector<double> reference_values2{5.8, 6.0, 6.2, 6.4}; // template reference points
      const set<double>    central_values{170,170.5,6.0,6.2 }; // drop 'extreme' 2D points
      //const vector<double> reference_values{170, 172}; // template reference points
      int seed = 1234;
      for ( double mean : reference_values1 ) {
         for ( double sigm : reference_values2 ) {
            if ( central_values.find(mean) != central_values.end() || central_values.find(sigm) != central_values.end() ) {
               templates[{mean,sigm}] = LTF_ROOTTools::MakeHistogram(nEventsTemplates, seed++, mean , sigm ,bins);
               templates[{mean,sigm}]->SetTitle(Form("m = %.1f",mean));
               templates[{mean,sigm}]->Scale(double(nEventsData)/nEventsTemplates);
            }
         }
      }
   }
   
   // ------------------------------------------------ //
   // ---  generate pseudo-data
   // ------------------------------------------------ //
   for ( int ii = 0 ; ii<1 ; ii++ ) 
   {

      TH1D* data = LTF_ROOTTools::MakeHistogram(nEventsData, seeddata+ii, meandata , sigmadata ,bins);

      PrintAsciiTable(templates,data);

   // ------------------------------------------------ //
   // ---  Do linear template fit
   // ------------------------------------------------ //

   // --- instantiate LTF object
   LTF ltf;
   ltf.SetGamma(vector<double>{1,1});
   ltf.UseNuisanceParameters(true);
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

   TH1D* normerr = (TH1D*)data->Clone("normerr"); 
   normerr->Scale(0.1);
   // ltf.AddError("normaelisation", bins.size()-1,normerr->GetArray()+1 , 1., LTF::Uncertainty::Constrained);
   // ltf.AddError("nordmalisation", bins.size()-1,normerr->GetArray()+1 , 1., LTF::Uncertainty::Constrained);
   // ltf.AddError("nfdormalisation", bins.size()-1,normerr->GetArray()+1 , 1., LTF::Uncertainty::Constrained);
   // ltf.AddError("normalisation", bins.size()-1,normerr->GetArray()+1 , 1., LTF::Uncertainty::Constrained);


   LTF::LiTeFit fit2 = ltf.DoIterativeFitNewton(2,2,1);
   //LTF::LiTeFit fit2 = ltf.DoLiTeFit();
   fit2.PrintFull();

   // fit2.DoIterativeFitNewton() ;

   // LTF::LiTeFit fit2 = ltf.DoLiTeFit();
   // fit2.PrintFull();
   // fit2.DoIterativeFitNewton() ;

   // LTF::LiTeFit fit2 = ltf.DoLiTeFit();
   // fit2.DoIterativeFit(6,0.6,2,0);
   // fit2.PrintFull();

   // fit2.DoIterativeFit(6,0.6,2,1);
   //fit2.PrintFull();


   // LTF::LiTeFit fit = ltf.DoLiTeFit();
   // //fit.PrintFull();
   // Eigen::VectorXd ahat = fit.ahat;
   // for ( int kk = 0 ; kk<4 ; kk++ ) {
   //    //ahat = (fit.SolveLinearTemplateFit(2,1,ahat) - ahat)*0.6+ahat;
   //    fit.DoLiTeFit(2,1,ahat);
   //    ahat = (fit.ahat - ahat)*0.6+ahat;
   //    cout<<"kk: "<<kk<<endl<<ahat<<endl;
   // }

   // if ( ii==0 ) {
   //    fit.PrintFull();
   //    //
   LTF_ROOTTools::plotLiTeFit_2D(fit2,bins);
   // }
   // else {
   //    fit.PrintShort();
   // }

   }

   
   return 0;
  
}



//! ------------------------------------------------------------------------ //
//! --- write templates, and data, to ascii file
void PrintAsciiTable(const map<vector<double>,TH1D*>& templates, TH1D* data){
   cout<<endl;
   printf(" %11s %11s",Form("Data"),Form("Stat"));
   for ( auto [refs,hist] : templates )
      printf(" %11s %11s",Form("T_%5.1f_%3.1f",refs[0],refs[1]),
             Form("S_%5.1f_%3.1f",refs[0],refs[1]));
   cout<<endl;
   for ( int i=1; i<=data->GetNbinsX() ;i++ ) {
      printf(" %11.4f %11.4f",data->GetBinContent(i),data->GetBinError(i));
      for ( auto [refs,hist] : templates )
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
   return example2_LTF_gaus2D();
#else
   //! printout if compiled without ROOT
   std::cout<<"This example is working only if ROOT is available. Otherwise, please see 'example2_LTF_gaus2D_NoROOT'."<<std::endl;
   return 0;
#endif
}

