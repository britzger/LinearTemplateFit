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

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>
#include <LTF/LTF_ROOTTools.h>
#include <LTF/LTF.h>

void PrintAsciiTable(const map<double,TH1D*>&, TH1D* data);


//! ------------------------------------------------------------------------ //
//! main function
#ifndef __CLING__

int example_ATLAS_topmass();

int main(int ,const char **) {

   gROOT->SetBatch();

   return example_ATLAS_topmass();
}
#endif

int example_ATLAS_topmass() {
   using namespace std;

#ifdef __CLING__
   TH1D::AddDirectory(false);
#endif
   TH1::SetDefaultSumw2(true);

   map<double,TH1D*> templates;

   const TString histname     = "m_bl";
   const int     iRebin       = 1;
   const TString histnamedata = "unfolding_mbl_selected_NOSYS";
   const int     iRebinData   = 1;
   const TString pseudodatafile     = "/ptmp/mpp/jhessler/LTF/LinearTemplateFit/LTF_ROOT/examples/data/Ana_S3beta_Cluster__mtop_170_H_1246.root";
   const TString datafile = "/ptmp/mpp/jhessler/LTF/LinearTemplateFit/LTF_ROOT/examples/data/unfolding_SR_Whad_Final_l_Whad_particle_TUnfoldStandalone_OptionA_data_nonClosureAlternative.root";
   const TString datauncfile  = "/ptmp/mpp/jhessler/LTF/LinearTemplateFit/LTF_ROOT/examples/data/unfolding_SR_Whad_Final_l_Whad_particle_TUnfoldStandalone_OptionA_data_nonClosureAlternative.root";
   const TString migmaname    = "h2_mupt2_4p";

   templates[165] = TFile::Open("/ptmp/mpp/jhessler/LTF/LinearTemplateFit/LTF_ROOT/examples/data/Ana_S3beta_Cluster__mtop_165_H_1246.root")->Get<TH1D>(histname);
   templates[170] = TFile::Open("/ptmp/mpp/jhessler/LTF/LinearTemplateFit/LTF_ROOT/examples/data/Ana_S3beta_Cluster__mtop_170_H_1246.root")->Get<TH1D>(histname);
   templates[175] = TFile::Open("/ptmp/mpp/jhessler/LTF/LinearTemplateFit/LTF_ROOT/examples/data/Ana_S3beta_Cluster__mtop_175_H_1246.root")->Get<TH1D>(histname);

   for ( auto [MM,hist] : templates ) hist->Rebin(iRebin);
   TH1D* data      = TFile::Open(pseudodatafile)->Get<TH1D>(histname); // pseudo data, use Sherpa 3 with m_t = 170 GeV for now
   data -> Rebin(iRebinData);

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
      ltf.AddTemplate(MM,  hist->GetNbinsX(),  hist->GetArray()+1 ); // set template
      hist->Sumw2();
      ltf.AddTemplateErrorSquared("statY", MM , hist->GetNbinsX(), hist->GetSumw2()->GetArray()+1, 0.); // set template error dY
      cout<<"MM: "<<MM<<"\tnBins: "<<hist->GetNbinsX()<<endl;
   }
   
   // --- initialize data
   ltf.SetData( data->GetNbinsX(), data->GetArray()+1);
   cout<<"data:     "<<"\tnBins: "<<data->GetNbinsX()<<endl;
   // --- initialize data uncertainties
   data->Sumw2(); // only for now
   TH1D* data_ref = TFile::Open(datauncfile)->Get<TH1D>(histnamedata);
   data_ref->Print("All");
   //TFile *file = TFile::Open(datauncfile);

   std::unique_ptr<TFile> file(TFile::Open(datauncfile));

   if (!file || file->IsOpen() == kFALSE) {
      std::cerr << "Error: Couldn't open the file!" << std::endl;
      return 1;
   }

   // Get the list of keys in the file (this represents all objects)
   TIter next(file->GetListOfKeys());
   TObject *obj;
   TKey* key;
   // Loop over all keys (objects) in the file
   while ((key = (TKey*) next())) {
      obj = file->Get<TObject>(key->GetName());
      if (obj->InheritsFrom("TH1D")) {
         std::unique_ptr<TH1D> hist(dynamic_cast<TH1D*>(obj));
         if (hist) {
            string title = hist->GetTitle();
            if ((title.find("unfolding_error_mbl_selected_direct_envelope_") != std::string::npos) && (title.find("__1up") != std::string::npos)){
               cout<<title<<endl;
               //hist->Add(data_ref, -1);
               vector<double> tmp;
               double* entries = hist->GetArray();
               for (int i=1; i< hist->GetNbinsX(); i++) {
                  cout<<"Added uncertainty: "<<entries[i]<<" bin content "<<hist->GetBinContent(i)<<" nbins "<<hist->GetNbinsX()<<endl;
                  tmp.push_back(hist->GetBinContent(i));
               }
               ltf.AddErrorRelative(title, tmp);
               //ltf.AddCorrelatedError(title, tmp);
            }
            else {
               continue;
            }
         }
      }
   }

   //ltf.AddUncorrelatedErrorSquared("stat.", data->GetNbinsX(), data->GetSumw2()->GetArray()+1);

   // for ( const auto& s : shiftsnuisance ) ltf.AddError("",N,s->GetArray()+1,1.);
   // for ( const auto& s : shifts ) ltf.AddError("",N,s->GetArray()+1,1.);
   // for ( const auto& s : shifts ) ltf.AddError("",N,s->GetArray()+1,0.5,LTF::Uncertainty::External);
   

   LTF::LiTeFit fit = ltf.DoLiTeFit();
   fit.PrintFull();
   //fit.DoIterativeFitNewton(6,0.6,2,1);
   //fit.DoIterativeFitTaylor();
   //fit.PrintFull();
   //LTF_ROOTTools::plotLiTeFit(fit,bins);
   
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
