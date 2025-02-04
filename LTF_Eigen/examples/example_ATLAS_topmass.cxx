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
#include <TH2D.h>

void PrintAsciiTable(const map<double,TH1D*>&, TH1D* data);

vector<vector<double > > TH2D_to_vecvec(TH2D* hist2D) {
   const int n = hist2D->GetNbinsX() - 1;
   vector<vector<double > > vecvec(n);
   for ( size_t i = 0 ; i<n ; i++ ) {
      vecvec[i].resize(n);
      for ( size_t j = 0 ; j<n ; j++ ) {
         vecvec[i][j] = hist2D->GetBinContent(i+1,j+1);
      }
   }
   return vecvec;
}



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
   const TString histnamedata = "unfolding_mbl_selected_NOSYS";
   //const TString histname     = "m_bw";
   //const TString histnamedata = "unfolding_mbwhad_selected_NOSYS";
   //const TString histname     = "m_wbbl";
   //const TString histnamedata = "unfolding_mwhadbbl_NOSYS";
   //const TString histname     = "m_minimax";
   //const TString histnamedata = "unfolding_minimax_whadbbl_NOSYS";
   const int     iRebin       = 1;
   const int     iRebinData   = 1;
   const TString pseudodatafile     = "/ptmp/mpp/jhessler/LTF/LinearTemplateFit/LTF_ROOT/examples/data/output/Ana_S3beta_Cluster_H_mtop_170_1248.root";
   const TString datafile = "/ptmp/mpp/jhessler/LTF/LinearTemplateFit/LTF_ROOT/examples/data/unfolding_SR_Whad_Final_l_Whad_particle_TUnfoldStandalone_OptionA_data_nonClosureAlternative.root";
   const TString datauncfile  = "/ptmp/mpp/jhessler/LTF/LinearTemplateFit/LTF_ROOT/examples/data/unfolding_SR_Whad_Final_l_Whad_particle_TUnfoldStandalone_OptionA_data_nonClosureAlternative.root";
   const TString migmaname    = "h2_mupt2_4p";
   
   double scale = TFile::Open(datafile)->Get<TH1D>(histnamedata)->Integral();
   cout<<"Integral "<<TFile::Open(datafile)->Get<TH1D>(histnamedata)->Integral()<<endl;
   templates[155] = TFile::Open("/ptmp/mpp/jhessler/LTF/LinearTemplateFit/LTF_ROOT/examples/data/output/Ana_S3beta_Cluster_H_mtop_155_1258.root")->Get<TH1D>(histname);
   templates[160] = TFile::Open("/ptmp/mpp/jhessler/LTF/LinearTemplateFit/LTF_ROOT/examples/data/output/Ana_S3beta_Cluster_H_mtop_160_1256.root")->Get<TH1D>(histname);
   templates[165] = TFile::Open("/ptmp/mpp/jhessler/LTF/LinearTemplateFit/LTF_ROOT/examples/data/output/Ana_S3beta_Cluster_H_mtop_165_1246.root")->Get<TH1D>(histname);
   templates[170] = TFile::Open("/ptmp/mpp/jhessler/LTF/LinearTemplateFit/LTF_ROOT/examples/data/output/Ana_S3beta_Cluster_H_mtop_170_1248.root")->Get<TH1D>(histname);
   templates[175] = TFile::Open("/ptmp/mpp/jhessler/LTF/LinearTemplateFit/LTF_ROOT/examples/data/output/Ana_S3beta_Cluster_H_mtop_175_1250.root")->Get<TH1D>(histname);
   templates[180] = TFile::Open("/ptmp/mpp/jhessler/LTF/LinearTemplateFit/LTF_ROOT/examples/data/output/Ana_S3beta_Cluster_H_mtop_180_1252.root")->Get<TH1D>(histname);
   templates[185] = TFile::Open("/ptmp/mpp/jhessler/LTF/LinearTemplateFit/LTF_ROOT/examples/data/output/Ana_S3beta_Cluster_H_mtop_185_1254.root")->Get<TH1D>(histname);


   for ( auto [MM,hist] : templates ) {
      hist->Rebin(iRebin);
      //hist->Scale(scale);
   }
   TH1D* data      = TFile::Open(pseudodatafile)->Get<TH1D>(histname); // pseudo data, use Sherpa 3 with m_t = 170 GeV for now
   data -> Rebin(iRebinData);
   //data->Scale(scale);
   data->SetLineColor(kBlack);
   data->SetMarkerSize(1.8);
   
   PrintAsciiTable(templates,data);

   // ------------------------------------------------ //
   // ---  Do linear template fit
   // ------------------------------------------------ //
   // --- instantiate LTF object
   LTF ltf;
   ltf.SetGamma(vector<double>{1});
   ltf.UseNuisanceParameters(true);
   ltf.UseLogNormalUncertainties(false);

   // --- initialize templates
   int test_tmp = 1;
   for ( auto [MM,hist] : templates ) {
      ltf.AddTemplate(MM,  hist->GetNbinsX(),  hist->GetArray()+1 ); // set template
      hist->Sumw2();
      ltf.AddTemplateErrorSquared("statY", MM , hist->GetNbinsX(), hist->GetSumw2()->GetArray()+1, 0.); // set template error dY
      cout<<"MM: "<<MM<<"\tnBins: "<<hist->GetNbinsX()<<endl;
     
   }

   // In case you need to rebin by hand (e.g. to remove the last (overflow) bin
   //double bin_edges[] = {0,   40,  80,  120, 160, 200, 280, 360, 480, 640, 960};
   //TH1D* tmp_data = new TH1D("h1", "h1", 10, bin_edges);
   //for ( int i = 1; i <= tmp_data->GetNbinsX(); i++ ) {
   //   tmp_data->SetBinContent(i, data->GetBinContent(i));
   //   tmp_data->SetBinError(i, std::sqrt(data->GetBinContent(i)));
   //}

   // --- initialize data
   //for ( int i = 1; i <= data->GetNbinsX(); i++ ) data->SetBinError(i, std::sqrt(data->GetBinContent(i)));
   ltf.SetData( data->GetNbinsX(), data->GetArray()+1);
   ltf.AddUncorrelatedErrorSquared("stat.", data->GetNbinsX(), data->GetSumw2()->GetArray()+1);
   
   // --- initialize data uncertainties
   data->Sumw2(); // only for now
   TH1D* data_ref = TFile::Open(datauncfile)->Get<TH1D>(histnamedata);
   data_ref->Print("All");

  
   std::unique_ptr<TFile> file(TFile::Open(datauncfile));

   if (!file || file->IsOpen() == kFALSE) {
      std::cerr << "Error: Couldn't open the file!" << std::endl;
      return 1;
   }
   
   TString histnameCovStat(histnamedata); // "unfolding_mbl_selected_NOSYS"  ->  unfolding_covariance_matrix_ptl1_covariance_STAT_DATA                                                                      
   histnameCovStat.ReplaceAll("unfolding_","unfolding_covariance_matrix_");
   histnameCovStat.ReplaceAll("_NOSYS","_covariance_STAT_DATA");
   TH2D* cov_stat_dat = TFile::Open(datauncfile)->Get<TH2D>(histnameCovStat);
   if ( !cov_stat_dat ) { cerr<<"Could not find covariance matrix " << histnameCovStat <<endl; exit(1);}
   //vector<vector<double > > vecCov2 = TH2D_to_vecvec(cov_stat_dat);
   //ltf.AddError("stat.", vecCov2); //Johannes this is where the covariance matrix was passed

   // cov_stat_dat->Print("all");                                                                                                                                                                           
   // cout<<"data_ref"<<endl;                                                                                                                                                                               
   // data_ref->Print("all");                                                                                                                                                                               
   // cout<<"data"<<endl;                                                                                                                                                                                   
   // data->Print("all");                                                                                                                                                                                   
   // cout<<"Ndata points in templates: "<<templates[170]->GetNbinsX()<<endl;                                                                                                                               
   // cout<<"Template: " <<endl;                                                                                                                                                                            
   // templates[170]->Print("all");         

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
            //if ((title.find("unfolding_error_mbwhad_selected_direct_envelope_") != std::string::npos) && (title.find("__1up") != std::string::npos)){
            //if ((title.find("unfolding_error_mwhadbbl_direct_envelope_") != std::string::npos) && (title.find("__1up") != std::string::npos)){
            //if ((title.find("unfolding_error_minimax_whadbbl_direct_envelope_") != std::string::npos) && (title.find("__1up") != std::string::npos)){
               cout<<title<<endl;
               vector<double> tmp;
               double* entries = hist->GetArray();
               for (int i=1; i< hist->GetNbinsX(); i++) {
                  tmp.push_back(hist->GetBinContent(i));
               }
               double corr = 1;
               if ( title.find("_STAT_")!= std::string::npos) corr = 0;
               if ( title.find("_TOTAL_")!= std::string::npos) ltf.AddErrorRelative(title, tmp, corr, LTF::Uncertainty::External);
               else  ltf.AddErrorRelative(title, tmp, corr, LTF::Uncertainty::Constrained);
               //ltf.AddCorrelatedError(title, tmp);
            }
            else {
               continue;
            }
         }
      }
   }
   ltf.AddUncorrelatedErrorSquared("stat.", data->GetNbinsX(), data->GetSumw2()->GetArray()+1);

   // for ( const auto& s : shiftsnuisance ) ltf.AddError("",N,s->GetArray()+1,1.);
   // for ( const auto& s : shifts ) ltf.AddError("",N,s->GetArray()+1,1.);
   // for ( const auto& s : shifts ) ltf.AddError("",N,s->GetArray()+1,0.5,LTF::Uncertainty::External);
   

   LTF::LiTeFit fit = ltf.DoLiTeFit();
   fit.PrintFull();


   //fit.DoIterativeFitNewton(6,0.6,2,1);
   //fit.DoIterativeFitTaylor();
   //fit.PrintFull();
   //const vector<double> bins{0, 40.0,  60.0,  80.0,  100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 220.0 };
   vector<double> bins{};
   for(int i = 1; i <= data->GetNbinsX(); i++) {
      bins.push_back(data->GetBinLowEdge(i));
      cout<<"Bin "<<i<<" "<<data->GetBinLowEdge(i)<<endl;
   }
   bins.push_back(data->GetXaxis()->GetBinUpEdge(data->GetNbinsX()));

   double* bins1 = data->GetArray()+1;
   cout<<bins1[0]<<" "<<bins1[-1]<<endl;
//for ( auto number : bins1 ) std::cout<<number<<std::endl;
   LTF_ROOTTools::plotLiTeFit(fit, bins,"1/#sigma d#sigma/dx","","m_{bl}");
   
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

