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
   const int var_name_index = 1;
   const vector<string> var_name = {"mbl_selected", "mbwhad_selected", "mwhadbbl", "minimax_whadbbl", "dRbl_selected", "dRbwhad_selected", "ptl1", "ptb1", "mwhad", "rapiditywhad"};
   const vector<string> var_name_short = {"m_bl", "m_bw", "m_wbbl", "m_minimax", "dr_bl", "dr_bw", "pT_lep1", "pT_bjet1", "m_whad", "y_whad"};
   const vector<TString> fit_vars = {"m_bl", "m_bw"};

   if ( !(var_name_index<var_name.size())  ) { cerr<<"Make sure to set var_name_index to a correct value!"<<endl; exit(1);}

   const TString histname     = var_name_short[var_name_index];
   const TString histnamedata = "unfolding_"+var_name[var_name_index]+"_NOSYS";
   const int     iRebin       = 1;
   const int     iRebinData   = 1;
   const TString pseudodatafile     = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_170_1248.root";
   const TString datafile = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/unfolding_SR_Whad_Final_l_Whad_particle_TUnfoldStandalone_OptionA_data_nonClosureAlternative.root";
   //const TString datauncfile  = "/ptmp/mpp/jhessler/LTF/LinearTemplateFit/LTF_ROOT/examples/data/unfolding_SR_Whad_Final_l_Whad_particle_TUnfoldStandalone_OptionA_data_nonClosureAlternative.root";
   const TString migmaname    = "h2_mupt2_4p";
   
   //double scale = TFile::Open(datafile)->Get<TH1D>(histnamedata)->Integral();
   //cout<<"Integral "<<TFile::Open(datafile)->Get<TH1D>(histnamedata)->Integral()<<endl;

   int bins_number = 0;
   for ( auto& tmp: fit_vars ) {
     TH1D* tmp_data = TFile::Open(pseudodatafile)->Get<TH1D>(tmp);
     bins_number += tmp_data->GetNbinsX();
     tmp_data->Clear();
   }
   TH1D* combined_data = new TH1D("combined_data", "combined_data", bins_number, 0, bins_number);
   int bin_offset = 0;
   for ( auto& tmp: fit_vars ) {
     TH1D* tmp_data = TFile::Open(pseudodatafile)->Get<TH1D>(tmp);
     for ( int i = 1; i <= tmp_data->GetNbinsX(); i++ ) {
       combined_data->SetBinContent(i+bin_offset, tmp_data->GetBinContent(i));
     }
     bin_offset += tmp_data->GetNbinsX();
   }
   
   TH1D* data = TFile::Open(pseudodatafile)->Get<TH1D>(histname); // pseudo data, use Sherpa 3 with m_t = 170 GeV for now

   double binning[9] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0};
   //TH1D* data_tmp      = TFile::Open(pseudodatafile)->Get<TH1D>(histname); // pseudo data, use Sherpa 3 with m_t = 170 GeV for now
   //TH1D* data      = new TH1D("data", "data", 8, binning);
   //for ( int i =1; i <= data->GetNbinsX(); i++) {
   //   data->SetBinContent(i, data_tmp->GetBinContent(i));
   //   data->SetBinError(i, data_tmp->GetBinError(i));
   //}
   data -> Rebin(iRebinData);
   //data->Scale(scale);
   data->SetLineColor(kBlack);
   data->SetMarkerSize(1.8);

   if ( var_name_index == 4 || var_name_index == 5 ) {
      
      TH1D* template_1_tmp = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_155_1258.root")->Get<TH1D>(histname);
      TH1D* template_2_tmp = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_160_1256.root")->Get<TH1D>(histname);
      TH1D* template_3_tmp = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_165_1246.root")->Get<TH1D>(histname);
      TH1D* template_4_tmp = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_170_1248.root")->Get<TH1D>(histname);
      TH1D* template_5_tmp = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_175_1250.root")->Get<TH1D>(histname);
      TH1D* template_6_tmp = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_180_1252.root")->Get<TH1D>(histname);
      TH1D* template_7_tmp = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_185_1254.root")->Get<TH1D>(histname);
      
      TH1D* template_1   = new TH1D("template_155", "template_155", 8, binning);
      TH1D* template_2   = new TH1D("template_160", "template_160", 8, binning);
      TH1D* template_3   = new TH1D("template_165", "template_165", 8, binning);
      TH1D* template_4   = new TH1D("template_170", "template_170", 8, binning);
      TH1D* template_5   = new TH1D("template_175", "template_175", 8, binning);
      TH1D* template_6   = new TH1D("template_180", "template_180", 8, binning);
      TH1D* template_7   = new TH1D("template_185", "template_185", 8, binning);
      for ( int i =1; i <= data->GetNbinsX(); i++) {
         template_1->SetBinContent(i, template_1_tmp->GetBinContent(i));
         template_1->SetBinError(  i, template_1_tmp->GetBinError(i));
         template_2->SetBinContent(i, template_2_tmp->GetBinContent(i));
         template_2->SetBinError(  i, template_2_tmp->GetBinError(i));
         template_3->SetBinContent(i, template_3_tmp->GetBinContent(i));
         template_3->SetBinError(  i, template_3_tmp->GetBinError(i));
         template_4->SetBinContent(i, template_4_tmp->GetBinContent(i));
         template_4->SetBinError(  i, template_4_tmp->GetBinError(i));
         template_5->SetBinContent(i, template_5_tmp->GetBinContent(i));
         template_5->SetBinError(  i, template_5_tmp->GetBinError(i));
         template_6->SetBinContent(i, template_6_tmp->GetBinContent(i));
         template_6->SetBinError(  i, template_6_tmp->GetBinError(i));
         template_7->SetBinContent(i, template_7_tmp->GetBinContent(i));
         template_7->SetBinError(  i, template_7_tmp->GetBinError(i));
      }
      
      templates[155] = template_1;
      templates[160] = template_2;
      templates[165] = template_3;
      templates[170] = template_4;
      templates[175] = template_5;
      templates[180] = template_6;
      templates[185] = template_7;
   }
   else {
     TH1D* combined_template_155 = new TH1D("combined_template_155", "combined_template_155", bins_number, 0, bins_number);
     TH1D* combined_template_160 = new TH1D("combined_template_160", "combined_template_160", bins_number, 0, bins_number);
     TH1D* combined_template_165 = new TH1D("combined_template_165", "combined_template_165", bins_number, 0, bins_number);
     TH1D* combined_template_170 = new TH1D("combined_template_170", "combined_template_170", bins_number, 0, bins_number);
     TH1D* combined_template_175 = new TH1D("combined_template_175", "combined_template_175", bins_number, 0, bins_number);
     TH1D* combined_template_180 = new TH1D("combined_template_180", "combined_template_180", bins_number, 0, bins_number);
     TH1D* combined_template_185 = new TH1D("combined_template_185", "combined_template_185", bins_number, 0, bins_number);
     int bin_offset = 0;
     for ( auto& tmp: fit_vars ) {
       TH1D* h_tmp_155 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_155_1258.root")->Get<TH1D>(tmp);
       TH1D* h_tmp_160 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_160_1256.root")->Get<TH1D>(tmp);
       TH1D* h_tmp_165 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_165_1246.root")->Get<TH1D>(tmp);
       TH1D* h_tmp_170 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_170_1248.root")->Get<TH1D>(tmp);
       TH1D* h_tmp_175 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_175_1250.root")->Get<TH1D>(tmp);
       TH1D* h_tmp_180 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_180_1252.root")->Get<TH1D>(tmp);
       TH1D* h_tmp_185 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_185_1254.root")->Get<TH1D>(tmp);
       for ( int i = 1; i<= h_tmp_155->GetNbinsX(); i++ ) {
	 combined_template_155->SetBinContent(i+bin_offset, h_tmp_155->GetBinContent(i));
	 combined_template_160->SetBinContent(i+bin_offset, h_tmp_160->GetBinContent(i));
	 combined_template_165->SetBinContent(i+bin_offset, h_tmp_165->GetBinContent(i));
	 combined_template_170->SetBinContent(i+bin_offset, h_tmp_170->GetBinContent(i));
	 combined_template_175->SetBinContent(i+bin_offset, h_tmp_175->GetBinContent(i));
	 combined_template_180->SetBinContent(i+bin_offset, h_tmp_180->GetBinContent(i));
	 combined_template_185->SetBinContent(i+bin_offset, h_tmp_185->GetBinContent(i));
       }
       bin_offset =+ h_tmp_155->GetNbinsX();
     }
     templates[155] = combined_template_155;
     templates[160] = combined_template_160;
     templates[165] = combined_template_165;
     templates[170] = combined_template_170;
     templates[175] = combined_template_175;
     templates[180] = combined_template_180;
     templates[185] = combined_template_185;
     
   }
   exit(0);

   for ( auto [MM,hist] : templates ) {
      hist->Rebin(iRebin);
      //hist->Scale(scale);
   }
   
   PrintAsciiTable(templates,data);

   // ------------------------------------------------ //
   // ---  Do linear template fit
   // ------------------------------------------------ //
   // --- instantiate LTF object
   LTF ltf;
   ltf.SetGamma(vector<double>{1});
   ltf.UseNuisanceParameters(true);
   ltf.UseLogNormalUncertainties(true);

   // --- initialize templates
   int test_tmp = 1;
   for ( auto [MM,hist] : templates ) {
      ltf.AddTemplate(MM,  hist->GetNbinsX(),  hist->GetArray()+1 ); // set template
      //for(int i = 1; i <= hist->GetNbinsX(); i++) {hist->SetBinError(i, std::sqrt(hist->GetBinContent(i)));}
      ltf.AddTemplateErrorSquared("statY", MM , hist->GetNbinsX(), hist->GetSumw2()->GetArray()+1, 0.); // set template error dY //Johannes add the correct value to the template root files!
      cout<<"MM: "<<MM<<"\tnBins: "<<hist->GetNbinsX()<<endl;
   }


   // --- initialize data
   ltf.SetData( data->GetNbinsX(), data->GetArray()+1);
   //ltf.AddUncorrelatedErrorSquared("stat.", data->GetNbinsX(), data->GetSumw2()->GetArray()+1);
   
   std::unique_ptr<TFile> file(TFile::Open(datafile));

   if (!file || file->IsOpen() == kFALSE) {
      std::cerr << "Error: Couldn't open the file!" << std::endl;
      return 1;
   }

   bool passDataCovMatrix = true;
   if ( passDataCovMatrix ) {
      TString histnameCovStat(histnamedata); // "unfolding_mbl_selected_NOSYS"  ->  unfolding_covariance_matrix_ptl1_covariance_STAT_DATA
      histnameCovStat.ReplaceAll("unfolding_","unfolding_covariance_matrix_");
      histnameCovStat.ReplaceAll("_NOSYS","_covariance_STAT_DATA");
      TH2D* cov_stat_dat = TFile::Open(datafile)->Get<TH2D>(histnameCovStat);
      if ( !cov_stat_dat ) { cerr<<"Could not find covariance matrix " << histnameCovStat <<endl; exit(1);}
      else cout<<"Found covarinace matrix "<<histnameCovStat<<endl;
   
      vector<vector<double > > vecCov2 = TH2D_to_vecvec(cov_stat_dat);
      for ( auto& tmp_vec: vecCov2 ){
         for ( auto& tmp: tmp_vec ) cout<<tmp<<"\t";
         cout<<endl;
      }
      ltf.AddErrorRelative("stat.", vecCov2);
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
            if ((title.find("unfolding_error_"+var_name[var_name_index]+"_direct_envelope_") != std::string::npos) && (title.find("__1up") != std::string::npos)){
               cout<<title<<endl;
               vector<double> tmp;
               double* entries = hist->GetArray();
               for (int i=1; i< hist->GetNbinsX(); i++) {
                  tmp.push_back(hist->GetBinContent(i));
               }
               double corr = 1;
               if ( title.find("_STAT_DATA")!= std::string::npos ) {
                  if ( !passDataCovMatrix ) ltf.AddErrorRelative("stat.", tmp, 0.0, LTF::Uncertainty::Constrained);
               }
               else if ( title.find("_STAT_MC")!= std::string::npos )                ltf.AddErrorRelative(title, tmp, 0.0, LTF::Uncertainty::Constrained);
               else if ( title.find("_FULL_SYS_SUM_")!= std::string::npos)           ltf.AddErrorRelative(title, tmp, corr, LTF::Uncertainty::External);
               else if ( title.find("_TOTAL_SYSONLY__1up")!= std::string::npos)      ltf.AddErrorRelative(title, tmp, corr, LTF::Uncertainty::External);
               else if ( title.find("_TOTAL__1up")!= std::string::npos)              ltf.AddErrorRelative(title, tmp, corr, LTF::Uncertainty::External);
               else if ( title.find("_TOTAL_NO_DR_DS__1up")!= std::string::npos)     ltf.AddErrorRelative(title, tmp, corr, LTF::Uncertainty::External);
               else  ltf.AddErrorRelative(title, tmp, corr, LTF::Uncertainty::Constrained);
            }
            else {
               continue;
            }
         }
      }
   }

   LTF::LiTeFit fit = ltf.DoLiTeFit();
   fit.PrintFull();

   //fit.DoIterativeFitNewton(6,0.6,2,1);
   //fit.DoIterativeFitTaylor();
   //fit.PrintFull();
   vector<double> bins{};
   for(int i = 1; i <= data->GetNbinsX(); i++) {
      bins.push_back(data->GetBinLowEdge(i));
   }
   bins.push_back(data->GetXaxis()->GetBinUpEdge(data->GetNbinsX()));

   double* bins1 = data->GetArray()+1;
   LTF_ROOTTools::plotLiTeFit(fit, bins,"1/#sigma d#sigma/dx","",var_name_index);
   
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

