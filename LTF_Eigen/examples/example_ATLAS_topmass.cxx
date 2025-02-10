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
   //const int var_name_index = 0;
   //const vector<string> var_name = {"mbl_selected", "mbwhad_selected", "mwhadbbl", "minimax_whadbbl", "dRbl_selected", "dRbwhad_selected", "ptl1", "ptb1", "mwhad", "rapiditywhad"};
   //const vector<string> var_name_short = {"m_bl", "m_bw", "m_wbbl", "m_minimax", "dr_bl", "dr_bw", "pT_lep1", "pT_bjet1", "m_whad", "y_whad"};
   const vector<TString> fit_vars = {"mbl_selected", "mbwhad_selected"};
   const vector<TString> fit_vars_short = {"m_bl", "m_bw"};
   
   const int     iRebin       = 1;
   const int     iRebinData   = 1;
   const TString pseudodatafile     = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_170_1248.root";
   const TString datafile = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/unfolding_SR_Whad_Final_l_Whad_particle_TUnfoldStandalone_OptionA_data_nonClosureAlternative.root";
   
   //double scale = TFile::Open(datafile)->Get<TH1D>(histnamedata)->Integral();
   //cout<<"Integral "<<TFile::Open(datafile)->Get<TH1D>(histnamedata)->Integral()<<endl;

   int bins_number = 0;
   for ( auto& tmp: fit_vars_short ) {
     TH1D* tmp_data = TFile::Open(pseudodatafile)->Get<TH1D>(tmp);
     bins_number += tmp_data->GetNbinsX();
     tmp_data->Clear();
   }
   TH1D* combined_data = new TH1D("combined_data", "combined_data", bins_number, 0, bins_number);
   int bin_offset = 0;
   for ( auto& tmp: fit_vars_short ) {
     TH1D* tmp_data = TFile::Open(pseudodatafile)->Get<TH1D>(tmp);
     for ( int i = 1; i <= tmp_data->GetNbinsX(); i++ ) {
       combined_data->SetBinContent(i+bin_offset, tmp_data->GetBinContent(i));
     }
     bin_offset += tmp_data->GetNbinsX();
   }
   combined_data -> Rebin(iRebinData);
   //combined_data->Scale(scale);
   combined_data->SetLineColor(kBlack);
   combined_data->SetMarkerSize(1.8);

   /*   
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
   else */
   {
     TH1D* combined_template_155 = new TH1D("combined_template_155", "combined_template_155", bins_number, 0, bins_number);
     TH1D* combined_template_160 = new TH1D("combined_template_160", "combined_template_160", bins_number, 0, bins_number);
     TH1D* combined_template_165 = new TH1D("combined_template_165", "combined_template_165", bins_number, 0, bins_number);
     TH1D* combined_template_170 = new TH1D("combined_template_170", "combined_template_170", bins_number, 0, bins_number);
     TH1D* combined_template_175 = new TH1D("combined_template_175", "combined_template_175", bins_number, 0, bins_number);
     TH1D* combined_template_180 = new TH1D("combined_template_180", "combined_template_180", bins_number, 0, bins_number);
     TH1D* combined_template_185 = new TH1D("combined_template_185", "combined_template_185", bins_number, 0, bins_number);
     int bin_offset = 0;
     for ( auto& tmp: fit_vars_short ) {
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

   for ( auto [MM,hist] : templates ) {
      hist->Rebin(iRebin);
      //hist->Scale(scale);
   }
   
   PrintAsciiTable(templates,combined_data);

   // ------------------------------------------------ //
   // --- List of uncertainties
   // ------------------------------------------------ //

   vector<string> uncertainties = {"EG_RESOLUTION_ALL", "EG_SCALE_ALL", "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR", "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR",
     "EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR", "EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR", "EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR", "MUON_SAGITTA_DATASTAT",
     "MUON_SAGITTA_RESBIAS", "MUON_EFF_BADMUON_SYS", "MUON_EFF_ISO_STAT", "MUON_EFF_ISO_SYS", "MUON_EFF_RECO_STAT", "MUON_EFF_RECO_SYS", "MUON_EFF_TTVA_STAT",
     "MUON_EFF_TTVA_SYS", "MUON_EFF_TrigStatUncertainty", "MUON_EFF_TrigSystUncertainty", // Lepton uncertainties
     "JET_EffectiveNP_Detector1", "JET_EffectiveNP_Detector2", "JET_EffectiveNP_Mixed1", "JET_EffectiveNP_Mixed2", "JET_EffectiveNP_Mixed3", "JET_EffectiveNP_Modelling1",
     "JET_EffectiveNP_Modelling2", "JET_EffectiveNP_Modelling3", "JET_EffectiveNP_Modelling4", "JET_EffectiveNP_Statistical1", "JET_EffectiveNP_Statistical2",
     "JET_EffectiveNP_Statistical3", "JET_EffectiveNP_Statistical4", "JET_EffectiveNP_Statistical5", "JET_EffectiveNP_Statistical6", "JET_EtaIntercalibration_Modelling",
     "JET_EtaIntercalibration_NonClosure_2018data", "JET_EtaIntercalibration_NonClosure_highE", "JET_EtaIntercalibration_NonClosure_negEta",
     "JET_EtaIntercalibration_NonClosure_posEta", "JET_EtaIntercalibration_TotalStat", "JET_Flavor_Composition_prop", "JET_Flavor_Response_prop",
     "JET_Pileup_OffsetMu", "JET_Pileup_OffsetNPV", "JET_Pileup_PtTerm", "JET_Pileup_RhoTopology", "JET_PunchThrough_MC16", // JES uncertainty
     "JET_JER_DataVsMC_MC16_PseudoData", "JET_JER_EffectiveNP_1_PseudoData", "JET_JER_EffectiveNP_2_PseudoData", "JET_JER_EffectiveNP_3_PseudoData",
     "JET_JER_EffectiveNP_4_PseudoData", "JET_JER_EffectiveNP_5_PseudoData", "JET_JER_EffectiveNP_6_PseudoData", "JET_JER_EffectiveNP_7_PseudoData",
     "JET_JER_EffectiveNP_8_PseudoData", "JET_JER_EffectiveNP_9_PseudoData", "JET_JER_EffectiveNP_10_PseudoData", "JET_JER_EffectiveNP_11_PseudoData",
     "JET_JER_EffectiveNP_12restTerm_PseudoData", // JER uncertainty
     "JET_JvtEfficiency", "PRW_DATASF", "MET_SoftTrk_ResoPara", "MET_SoftTrk_ResoPerp", "MET_SoftTrk_Scale", "WMASS_VAR_signal", //MET+JVT+PileUp
     "FT_EFF_Eigen_B_0", "FT_EFF_Eigen_B_1", "FT_EFF_Eigen_B_2", "FT_EFF_Eigen_B_3", "FT_EFF_Eigen_B_4", "FT_EFF_Eigen_B_5", "FT_EFF_Eigen_B_6", "FT_EFF_Eigen_B_7",
     "FT_EFF_Eigen_B_8", "FT_EFF_Eigen_C_0", "FT_EFF_Eigen_C_1", "FT_EFF_Eigen_C_2", "FT_EFF_Eigen_C_3", "FT_EFF_Eigen_Light_0", "FT_EFF_Eigen_Light_1",
     "FT_EFF_Eigen_Light_2", "FT_EFF_Eigen_Light_3", "FT_EFF_extrapolation", "FT_EFF_extrapolation_from_charm", // b-tagging uncertainties
     "THEORY_CROSS_SECTION_signal", "THEORY_SHOWERING_HERWIG7_signal", "THEORY_SCALE_FACTORISATION_signal", "THEORY_SCALE_RENORMALISATION_signal",
     "THEORY_ISR_signal", "THEORY_FSR_signal", "THEORY_HDAMP_signal", "THEORY_PTHARD_signal", "THEORY_TOPRECOIL_signal", "THEORY_TOP_MASS_signal",
     "THEORY_PDF4LHC_VARIATION_signal", "THEORY_DR_DS_signal", // modelling uncertainties
     "THEORY_CROSS_SECTION_Wjets", "THEORY_SCALE_COMBINED_Wjets", "THEORY_PDF4LHC_VARIATION_Wjets", "THEORY_EWK_Wjets", "THEORY_SCALE_COMBINED_multiboson_noW",
     "THEORY_PDF4LHC_VARIATION_multiboson_noW", "THEORY_EWK_multiboson_noW", "THEORY_CROSS_SECTION_other_top_noWt", "FAKES_Electron", "FAKES_Muon", // bkgd uncertainties
     "LUMINOSITY" // lumi uncertainty
   };

   vector<string> statistical_uncertainties = {"STAT_DATA", "STAT_MC"};
   
   vector<string> external_uncertainties = {"FULL_SYS_SUM", "TOTAL_SYSONLY", "TOTAL", "TOTAL_NO_DR_DS"};

   // ------------------------------------------------ //
   // ---  Do linear template fit
   // ------------------------------------------------ //
   // --- instantiate LTF object
   LTF ltf;
   ltf.SetGamma(vector<double>{1});
   ltf.UseNuisanceParameters(true);
   ltf.UseLogNormalUncertainties(true);

   // --- initialize templates
   for ( auto [MM,hist] : templates ) {
      ltf.AddTemplate(MM,  hist->GetNbinsX(),  hist->GetArray()+1 ); // set template
      //for(int i = 1; i <= hist->GetNbinsX(); i++) {hist->SetBinError(i, std::sqrt(hist->GetBinContent(i)));}
      ltf.AddTemplateErrorSquared("statY", MM , hist->GetNbinsX(), hist->GetSumw2()->GetArray()+1, 0.); // set template error dY //Johannes add the correct value to the template root files!
      cout<<"MM: "<<MM<<"\tnBins: "<<hist->GetNbinsX()<<endl;
   }

   // --- initialize data
   ltf.SetData( combined_data->GetNbinsX(), combined_data->GetArray()+1);
   //ltf.AddUncorrelatedErrorSquared("stat.", data->GetNbinsX(), data->GetSumw2()->GetArray()+1);
   
   std::unique_ptr<TFile> file(TFile::Open(datafile));

   if (!file || file->IsOpen() == kFALSE) {
      std::cerr << "Error: Couldn't open the file!" << std::endl;
      return 1;
   }
   /*
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
   */

   // Systematical uncertainties
   for ( auto& uncertainty: uncertainties ) {
     vector<double> combined_error;
     for ( auto& fit_variable: fit_vars ) {
       TH1D* hist = file->Get<TH1D>("unfolding_error_"+fit_variable+"_direct_envelope_"+uncertainty+"__1up");
       for (int i=1; i< hist->GetNbinsX(); i++) {
	 combined_error.push_back(hist->GetBinContent(i));
       }
     }
     double corr = 1;
     if ( combined_error.size() > 0 ) ltf.AddErrorRelative(uncertainty, combined_error, corr, LTF::Uncertainty::Constrained);
     combined_error.clear();
   }
   // Stat uncertainties
   for ( auto& uncertainty: statistical_uncertainties ) {
     vector<double> combined_error;
     for ( auto& fit_variable: fit_vars ) {
       TH1D* hist = file->Get<TH1D>("unfolding_error_"+fit_variable+"_direct_envelope_"+uncertainty+"__1up");
       for (int i=1; i< hist->GetNbinsX(); i++) {
         combined_error.push_back(hist->GetBinContent(i));
       }
     }
     double corr = 0;
     if ( combined_error.size() > 0 ) ltf.AddErrorRelative(uncertainty, combined_error, corr, LTF::Uncertainty::Constrained);
     combined_error.clear();
   }

   // External uncertainties
   for ( auto& uncertainty: external_uncertainties ) {
     vector<double> combined_error;
     for ( auto& fit_variable: fit_vars ) {
       TH1D* hist = file->Get<TH1D>("unfolding_error_"+fit_variable+"_direct_envelope_"+uncertainty+"__1up");
       for (int i=1; i< hist->GetNbinsX(); i++) {
         combined_error.push_back(hist->GetBinContent(i));
       }
     }
     double corr = 1;
     if ( combined_error.size() > 0 ) ltf.AddErrorRelative(uncertainty, combined_error, corr, LTF::Uncertainty::External);
     combined_error.clear();
   }

   LTF::LiTeFit fit = ltf.DoLiTeFit();
   fit.PrintFull();

   //fit.DoIterativeFitNewton(6,0.6,2,1);
   //fit.DoIterativeFitTaylor();
   //fit.PrintFull();
   vector<double> bins{};
   for(int i = 1; i <= combined_data->GetNbinsX(); i++) {
      bins.push_back(combined_data->GetBinLowEdge(i));
   }
   bins.push_back(combined_data->GetXaxis()->GetBinUpEdge(combined_data->GetNbinsX()));

   LTF_ROOTTools::plotLiTeFit(fit, bins,"1/#sigma d#sigma/dx", "m_{bl}\t m_{bW}","m_{t} [GeV]");
   
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

