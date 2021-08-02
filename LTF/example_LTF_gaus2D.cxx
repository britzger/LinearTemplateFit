//#include "include/LTF/LTF.h"
#include "LTF/LTF.h"
#include "LTF_ROOTTools.h"
#include <TH1D.h>
#include <string>

using namespace std;



// __________________________________________________________________________________ //
//! main
int example_LTF_gaus2D() {

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


int main() {
   return example_LTF_gaus2D();
}
