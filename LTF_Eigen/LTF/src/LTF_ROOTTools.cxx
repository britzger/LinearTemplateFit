//#include "include/LTF/LTF.h"
#include "LTF/LTF.h"
#include "LTF/LTF_ROOTTools.h"

#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMatrixDfwd.h>
#include <TMatrixD.h>
#include <TMatrixTSym.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TSystem.h>
#include <TFile.h>
#include <string>
#include <fstream>
#include <TRandom3.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TF1.h>
#include <TF2.h>
#include <TLine.h>
#include <TRandom3.h>
#include <TGraph2DErrors.h>
#include <TGraph2D.h>
#include <set>

using namespace std;

// __________________________________________________________________________________ //
//!
//! read_input_table2()
//!
//! read input data table from file 
//! 
std::map < std::string, std::vector<double> > LTF_ROOTTools::read_input_table2(std::string filename, int ncol ) 
{
   std::map < std::string, std::vector<double> > ret;
   std::vector<std::string> cols;
   // open file for reading                                           
   std::ifstream istrm(filename.c_str(), std::ios::binary);
   if (!istrm.is_open()) {
      std::cout << "failed to open " << filename << std::endl;
      exit(1);
   }
   for ( int c=0 ; c<ncol ; c++ ) {
      std::string colname;
      istrm >> colname;
      ret[colname] = vector<double>();
      cols.push_back(colname);
   }
   while ( istrm.good()) {
      double value;
      for ( int c=0 ; c<ncol ; c++ ) {
         istrm >> value;
         if ( !istrm.good() ) break;
         std::string colname = cols[c];
         ret[colname].push_back(value);
      }
   }
   cout<<"Info.  [read_input_table]  Read "<<ret.size()<<" rows."<<endl;
   return ret;
}



// __________________________________________________________________________________ //
//! 
//!  MakeTGraph
//!
//!  make a TGraph for plotting purposes
//!
TGraphErrors* LTF_ROOTTools::MakeTGraph(const Eigen::VectorXd& xvalues, int ibin, const Eigen::MatrixXd& Y,
                         const std::map<std::string,Eigen::MatrixXd >& VSysY ) 
{

   TGraphErrors* graph = new TGraphErrors();
   for ( int i = 0 ; i<xvalues.size() ; i++ ) {
      double yvalue = Y(ibin,i);
      graph->SetPoint(i,xvalues(i),yvalue);
      double ey2 = 0;
      for ( auto [name,Vy] : VSysY ) {
         ey2 += pow(Vy(ibin,i),2);
      }
      graph->SetPointError(i,0,sqrt(ey2));
   }
   return graph;
}


// __________________________________________________________________________________ //
//! 
//!  MakeTGraph2D
//!
//!  make a TGraph2D for plotting purposes
//!
TGraph2DErrors* LTF_ROOTTools::MakeTGraph2D(const Eigen::VectorXd& xvalues, const Eigen::VectorXd& yvalues,
                                            int ibin, const Eigen::MatrixXd& Y,
                                            const std::map<std::string,Eigen::MatrixXd >& VSysY ) 
{

   TGraph2DErrors* graph = new TGraph2DErrors();
   for ( int i = 0 ; i<xvalues.size() ; i++ ) {
      double zvalue = Y(ibin,i);
      graph->SetPoint(i,xvalues(i),yvalues(i),zvalue);
      double ey2 = 0;
      for ( auto [name,Vy] : VSysY ) {
         ey2 += pow(Vy(ibin,i),2);
      }
      graph->SetPointError(i,0,0,sqrt(ey2));
   }
   return graph;
}


// __________________________________________________________________________________ //
//!
//! MakeHistogram
//!
//! make a histogram and fill it with random events according to a gauss
//! distribution around M
//!
TH1D* LTF_ROOTTools::MakeHistogram(int nEvents, int seed, double mean, double sigma, vector<double> bins ) 
{
   TRandom3 rn(seed);
   
   TH1D* hist = new TH1D("hist","hist",bins.size()-1, &bins[0]);
   hist->Sumw2();
   for ( int i=0 ; i<nEvents; i++ ) {
      hist->Fill( rn.Gaus(mean,sigma), 1 );
   }
   return hist;
}




// __________________________________________________________________________________ //
//!
//!  MakeHistogram
//!
//!  make a histogram from an Eigen::Vector for plotting purposes
//!
TH1D* LTF_ROOTTools::MakeHistogram(const Eigen::VectorXd& values, vector<double> bins, const std::vector<std::pair<std::string,Eigen::MatrixXd > >& V ) 
{
   TH1D* hist = bins.empty() ?
      new TH1D("hist","hist",values.size(),0,values.size() ) :
      new TH1D("hist","hist",bins.size()-1, &bins[0]);
   if ( bins.size() && int(values.size()+1) != int(bins.size()) ) {cout<<"ERROR! binning and number of entries does not fit!"<<endl;exit(1);}
   for ( size_t i = 0 ;i<bins.size()-1; i++ ) {
      hist->SetBinContent(i+1, values(i));
      if ( V.size() ) {
         double e2 = 0;
         for ( auto apair : V ) {
            e2 += apair.second(i,i);
         }
         hist->SetBinError(i+1, sqrt(e2) );
      }
   }
   return hist;
}



// __________________________________________________________________________________ //
//!
//!
//!  Plot a LiTeFit object using ROOT
//!
//!  The binning needs to be provided to the plotting function,
//!  since this is not included in LTF::LiTeFit
//! 
void LTF_ROOTTools::plotLiTeFit(const LTF::LiTeFit& fit, const vector<double>& bins, 
				const char*   ps_name,
				const string& yaxistitle,
				const string& xaxistitle,
				const string& referencename) 
{
   gStyle->SetOptStat(0);
   gSystem->mkdir("plots");
   auto& M = fit.M;
   
   // sanity check
   if ( M.cols() != 2 ) {cout<<"Error! only 1-dim plotting implemented."<<endl;exit(1);}
   Eigen::VectorXd reference_values = M.col(1);
   //TH1D* hist = new TH1D("hist","hist",bins.size()-1, &bins[0]);
   
   //const vector<string> var_name_latex = {"m_{bl}", "m_{bW}", "m_{Wbbl}", "m_{Wb,bl}^{minimax}", "#Delta R(b,l)", "#Delta R(b,W)", "p_{T}(l)",  "p_{T}(b_{1})", "m(W_{had})", "y(W_{had})"};

   map<double,TH1D*> templates;
   for ( int iref = 0 ; iref<reference_values.size() ; iref++ ) {
      //double ref = reference_values(iref);
      templates[iref] = MakeHistogram(fit.Y.col(iref),bins);
   }

   TH1D* data    = MakeHistogram(fit.Dt,bins,fit.Vs);
   TH1D* TheoFit = MakeHistogram(fit.TheoFit,bins);
   
   TCanvas c1("c1","LTF plots",800,800);
   c1.SetRightMargin(0.05);
   c1.SetLeftMargin(0.15);
   c1.SetTopMargin(0.08);
   // c1.SetRightMargin(0.02);

   //const char* ps_name = fit.GetLogNormal() ?
   //   "plots/LTFlog_plots.ps" :
   //   "plots/LTF_plots.ps";
   c1.Print( (string(ps_name)+"[").c_str() );
   // ---------------------------------------------- //
   // main plot
   // ---------------------------------------------- //
   for ( int iref = 0 ; iref<reference_values.size() ; iref++ ) {
      templates[iref]->SetLineWidth(2);
      if ( iref == 0 ) {
         templates[0]->SetTitle(("Linear Template Fit;"+xaxistitle+";"+yaxistitle).c_str());
         templates[0]->SetLineColor(kRed+1);
         if ( templates[0]->GetMaximum()>0 )templates[0]->SetMinimum(0);
         if ( !fit.GetLogNormal() ) 
            templates[0]->SetMaximum(templates[0]->GetMaximum()*1.7);
         // else 
         //    templates[0]->SetMaximum(templates[0]->GetMaximum()*3.); 
         templates[0]->SetLineWidth(3);
         templates[0]->DrawClone("hist");
      }
      else if ( iref==reference_values.size()-1) {
         if ( templates[0]->GetMaximum()>0 ) templates[iref]->SetFillColorAlpha(kBlue,0.15);
         templates[iref]->SetLineColor(kBlue+2);
         templates[iref]->SetLineWidth(3);
         templates[iref]->Draw("histsame");
      }
      else {
         templates[iref]->SetLineColor(iref+2);
         templates[iref]->SetLineWidth(2);
         templates[iref]->Draw("histsame");
      }
   }   
   if ( templates[0]->GetMaximum()>0 )templates[0]->SetFillColorAlpha(kRed,0.15);
   templates[0]->Draw("histsame");

   data->SetMarkerStyle(20);
   data->SetMarkerSize(1.4);
   data->SetLineColor(kBlack);
   data->Draw("e0same");

   TheoFit->SetLineColor(923);
   TheoFit->SetLineWidth(4);
   TheoFit->SetLineStyle(2);
   TheoFit->Draw("histsame");

   TLegend legend(0.18,0.70,0.94,0.92,"","NDC");
   legend.SetNColumns(3);
   legend.SetFillStyle(0);
   legend.SetBorderSize(0);
   legend.AddEntry(data,"Data","E0P");
   for ( int iref = 0 ; iref<reference_values.size() ; iref++ ) {
      legend.AddEntry(templates[iref],Form("Template #alpha=%6.2f",reference_values[iref]),"FL");
   }
   legend.AddEntry(TheoFit,"Estimated best model","L");
   legend.Draw();

   c1.Print(ps_name);
   //c1.Print("plots/LTF_plot.pdf");


   // ---------------------------------------------- //
   // print relative size of all errors
   // ---------------------------------------------- //
   {
      c1.Clear();
      c1.SetLogy(0);
      c1.SetLeftMargin(0.2);
      
      vector<string> lepton_uncertainties = {"EG_RESOLUTION_ALL",
                                             "EG_SCALE_ALL",
                                             "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR",
                                             "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR",
                                             "EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR",
                                             "EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR",
                                             "EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR",
                                             "MUON_SAGITTA_DATASTAT",
                                             "MUON_SAGITTA_RESBIAS",
                                             "MUON_EFF_BADMUON_SYS",
                                             "MUON_EFF_ISO_STAT",
                                             "MUON_EFF_ISO_SYS",
                                             "MUON_EFF_RECO_STAT",
                                             "MUON_EFF_RECO_SYS",
                                             "MUON_EFF_TTVA_STAT",
                                             "MUON_EFF_TTVA_SYS",
                                             "MUON_EFF_TrigStatUncertainty",
                                             "MUON_EFF_TrigSystUncertainty"};
      
      vector<string> jes_uncertainties = {"JET_EffectiveNP_Detector1",
                                          "JET_EffectiveNP_Detector2",
                                          "JET_EffectiveNP_Mixed1",
                                          "JET_EffectiveNP_Mixed2",
                                          "JET_EffectiveNP_Mixed3",
                                          "JET_EffectiveNP_Modelling1",
                                          "JET_EffectiveNP_Modelling2",
                                          "JET_EffectiveNP_Modelling3",
                                          "JET_EffectiveNP_Modelling4",
                                          "JET_EffectiveNP_Statistical1",
                                          "JET_EffectiveNP_Statistical2",
                                          "JET_EffectiveNP_Statistical3",
                                          "JET_EffectiveNP_Statistical4",
                                          "JET_EffectiveNP_Statistical5",
                                          "JET_EffectiveNP_Statistical6",
                                          "JET_EtaIntercalibration_Modelling",
                                          "JET_EtaIntercalibration_NonClosure_2018data",
                                          "JET_EtaIntercalibration_NonClosure_highE",
                                          "JET_EtaIntercalibration_NonClosure_negEta",
                                          "JET_EtaIntercalibration_NonClosure_posEta",
                                          "JET_EtaIntercalibration_TotalStat",
                                          "JET_Flavor_Composition_prop",
                                          "JET_Flavor_Response_prop",
                                          "JET_Pileup_OffsetMu",
                                          "JET_Pileup_OffsetNPV",
                                          "JET_Pileup_PtTerm",
                                          "JET_Pileup_RhoTopology",
                                          "JET_PunchThrough_MC16"};


      vector<string> jer_uncertainties = {"JET_JER_DataVsMC_MC16_PseudoData",
                                          "JET_JER_EffectiveNP_1_PseudoData",
                                          "JET_JER_EffectiveNP_2_PseudoData",
                                          "JET_JER_EffectiveNP_3_PseudoData",
                                          "JET_JER_EffectiveNP_4_PseudoData",
                                          "JET_JER_EffectiveNP_5_PseudoData",
                                          "JET_JER_EffectiveNP_6_PseudoData",
                                          "JET_JER_EffectiveNP_7_PseudoData",
                                          "JET_JER_EffectiveNP_8_PseudoData",
                                          "JET_JER_EffectiveNP_9_PseudoData",
                                          "JET_JER_EffectiveNP_10_PseudoData",
                                          "JET_JER_EffectiveNP_11_PseudoData",
                                          "JET_JER_EffectiveNP_12restTerm_PseudoData"};

      vector<string> jvt_pileup_met_uncertainties = {"JET_JvtEfficiency",
                                                     "PRW_DATASF",
                                                     "MET_SoftTrk_ResoPara",
                                                     "MET_SoftTrk_ResoPerp",
                                                     "MET_SoftTrk_Scale",
                                                     "WMASS_VAR_signal"};


      vector<string> b_tagging_uncertainties = {"FT_EFF_Eigen_B_0",
                                                "FT_EFF_Eigen_B_1",
                                                "FT_EFF_Eigen_B_2",
                                                "FT_EFF_Eigen_B_3",
                                                "FT_EFF_Eigen_B_4",
                                                "FT_EFF_Eigen_B_5",
                                                "FT_EFF_Eigen_B_6",
                                                "FT_EFF_Eigen_B_7",
                                                "FT_EFF_Eigen_B_8",
                                                "FT_EFF_Eigen_C_0",
                                                "FT_EFF_Eigen_C_1",
                                                "FT_EFF_Eigen_C_2",
                                                "FT_EFF_Eigen_C_3",
                                                "FT_EFF_Eigen_Light_0",
                                                "FT_EFF_Eigen_Light_1",
                                                "FT_EFF_Eigen_Light_2",
                                                "FT_EFF_Eigen_Light_3",
                                                "FT_EFF_extrapolation",
                                                "FT_EFF_extrapolation_from_charm"};

      // Theory uncertainties                                                                                                                                                                               
      vector<string> modelling_uncertainties = {"THEORY_CROSS_SECTION_signal",
                                                "THEORY_SHOWERING_HERWIG7_signal",
                                                "THEORY_SCALE_FACTORISATION_signal",
                                                "THEORY_SCALE_RENORMALISATION_signal",
                                                "THEORY_ISR_signal",
                                                "THEORY_FSR_signal",
                                                "THEORY_HDAMP_signal",
                                                "THEORY_PTHARD_signal",
                                                "THEORY_TOPRECOIL_signal",
                                                "THEORY_TOP_MASS_signal",
                                                "THEORY_PDF4LHC_VARIATION_signal",
                                                "THEORY_DR_DS_signal"};


      vector<string> bkgd_uncertainties = {"THEORY_CROSS_SECTION_Wjets",
                                           "THEORY_SCALE_COMBINED_Wjets",
                                           "THEORY_PDF4LHC_VARIATION_Wjets",
                                           "THEORY_EWK_Wjets",
                                           "THEORY_SCALE_COMBINED_multiboson_noW",
                                           "THEORY_PDF4LHC_VARIATION_multiboson_noW",
                                           "THEORY_EWK_multiboson_noW",
                                           "THEORY_CROSS_SECTION_other_top_noWt",
                                           "FAKES_Electron",
                                           "FAKES_Muon"};

      vector<string> other_uncertainties = {"STAT_DATA",
                                            "STAT_MC",
                                            "LUMINOSITY"};

      std::map<string, double> error_summary;

      error_summary.insert({"Lepton",           makeErrorPlot(c1, ps_name, "Lepton uncertainties", fit, lepton_uncertainties)});
      error_summary.insert({"JES",              makeErrorPlot(c1, ps_name, "JES uncertainties", fit, jes_uncertainties)});
      error_summary.insert({"JER",              makeErrorPlot(c1, ps_name, "JER uncertainties", fit, jer_uncertainties)});
      error_summary.insert({"JVT+PileUp+MET",   makeErrorPlot(c1, ps_name, "JVT+PileUp+MET uncertainties", fit, jvt_pileup_met_uncertainties)});
      error_summary.insert({"b-tagging",        makeErrorPlot(c1, ps_name, "b-tagging uncertainties", fit, b_tagging_uncertainties)});
      error_summary.insert({"theory modelling", makeErrorPlot(c1, ps_name, "modelling uncertainties", fit, modelling_uncertainties)});
      error_summary.insert({"bkgd modelling",   makeErrorPlot(c1, ps_name, "background uncertainties", fit, bkgd_uncertainties)});
      error_summary.insert({"Stat.+Lumi",       makeErrorPlot(c1, ps_name, "statistical uncertainties", fit, other_uncertainties)});

      TH1D* h  = new TH1D("Full error breakdown", "Full error breakdown", error_summary.size()+1, 0, error_summary.size()+1);
      double sum_error_sq = 0;
      for( auto& tmp_err: error_summary ) {
         h->Fill(tmp_err.first.c_str(), tmp_err.second);
         sum_error_sq += pow(tmp_err.second,2);
      }

      h->SetBinContent(h->GetNbinsX(), std::sqrt(sum_error_sq));
      h->GetXaxis()->SetBinLabel(h->GetNbinsX(), "Total unc.");
      h->SetBarWidth(0.85);
      h->GetYaxis()->SetLabelSize(0.03);
      h->GetYaxis()->SetTitle("Uncertainty [GeV]");
      h->GetXaxis()->SetLabelSize(0.02);
      h->GetXaxis()->SetTickLength(0);
      h->GetXaxis()->LabelsOption("<");
      h->Draw("hbar");

      c1.Print(ps_name);
      c1.Clear();
      
   }

   // ---------------------------------------------- //
   // print linear-functions in every bin
   // ---------------------------------------------- //
   c1.Clear();
   c1.SetRightMargin(0.02);
   c1.SetTopMargin(0.02);
   c1.SetLeftMargin(0.16);
   c1.SetBottomMargin(0.16);

   gPad->SetTicky(1);

   gStyle->SetLabelSize(0.05,"XYZ");
   gStyle->SetTitleSize(0.05,"XYZ");
   gStyle->SetTitleOffset(1.1,"X");
   gStyle->SetTitleOffset(1.6,"Y");
   gStyle->SetMarkerSize(2);


   //map < string, vector<double> > input_table = read_input_table2("data/CMS_data.txt",32);

   for ( int ibin = 0 ; ibin<fit.Dt.size() ; ibin++ ) {
   //for ( int ibin = 0 ; ibin<0 ; ibin++ ) {
   //for ( int ibin = 0 ; ibin<6 ; ibin++ ) {
      //TGraphErrors* data  = MakeTGraph(fit.ahat.row(0),fit.Dt.row(ibin));
      TGraphErrors* gdata = new TGraphErrors();
      gdata->SetPoint(0,fit.ahat(0),fit.Dt(ibin));
      gdata->SetPointError(0,0,data->GetBinError(ibin+1));
      gdata->SetMarkerStyle(20);

      TGraphErrors* graph = MakeTGraph(fit.M.col(1),ibin,fit.Y,fit.SysY);

      graph->Fit("pol1","Q");
      TF1* pol1w = (TF1*)graph->GetFunction("pol1")->Clone("pol1w");
      pol1w->SetLineColor(922);
      pol1w->SetLineStyle(3);
      pol1w->SetLineWidth(2);

      TF1* f2= new TF1("f2","[0]+[1]*pow(x,2)", -FLT_MIN,FLT_MAX );
      graph->Fit(f2,"QW");
      TF1* f1= new TF1("f1","[0]+[1]*pow(x,[2])", -FLT_MIN,FLT_MAX );
      

      graph->Fit("pol1","QW"); // "W": Ignore all point errors when fitting a TGraphErrors 
      bool UseLTFOutput = true;
      if ( UseLTFOutput ) {
         Eigen::VectorXd ltfpol1param =(fit.Y*fit.Mc().transpose()).row(ibin);
         //cout<<"M+*Y:"<<endl<< ltfpol1param <<endl;
         graph->GetFunction("pol1")->SetParameter(0,ltfpol1param(0));
         graph->GetFunction("pol1")->SetParameter(1,ltfpol1param(1));
         f1->SetParameter(0,ltfpol1param(0));
         f1->SetParameter(1,ltfpol1param(1));
         f1->SetParameter(2,fit.Gamma[0]);
      }
      if ( fit.Gamma[0]!=1 )   cout<<"Warning! Plotting with gamma factor !=1 not correctly implmeneted!"<<endl;

      graph->GetFunction("pol1")->SetLineWidth(3);
      
      graph->Fit("pol2","QW");
      graph->GetFunction("pol2")->SetLineWidth(3);

      graph->SetMarkerStyle(47);
      graph->SetMarkerColor(kRed+2);
      graph->SetLineColor(kRed+2);
      if ( fit.GetLogNormal() ) 
         graph->SetTitle((";"+referencename+";log("+yaxistitle+")").c_str()); // log(value/unit)
      else
         graph->SetTitle((";"+referencename+";"+yaxistitle).c_str());


      TF1* f1log = NULL;
      if ( fit.GetLogNormal() ) {
         //f1log = new TF1("pol1log","[0]+[1]*exp(x)", -FLT_MIN,FLT_MAX );
         f1log = new TF1("pol1log","log([0]+[1]*pow(x,1))", -FLT_MIN,FLT_MAX );

         TGraph* gexp = new TGraph();
         for ( int i = 0 ; i<graph->GetN() ; i++ ) 
            gexp->SetPoint(i,graph->GetX()[i], exp(graph->GetY()[i]));
         gexp->Fit("pol1","QW");

         Eigen::VectorXd ltfpol1param =(fit.Y*fit.Mc().transpose()).row(ibin);
         if ( graph->GetFunction("pol1") )
            graph->GetFunction("pol1")->SetParameter(0,ltfpol1param(0));

         f1log->SetParameter(0, gexp->GetFunction("pol1")->GetParameter(0));
         f1log->SetParameter(1, gexp->GetFunction("pol1")->GetParameter(1));
         if ( fit.Gamma[0]!=1 )   cout<<"Warning! Plotting with gamma factor !=1 not correctly implmeneted!"<<endl;

         // //f1 = new TF1("pol1","[0]+[1]*x", -FLT_MIN,FLT_MAX );
         // graph->Fit(f1log,"QW"); // "W": Ignore all point errors when fitting a TGraphErrors 
         f1log->SetLineColor(kBlue+1);
         f1log->SetLineStyle(7);
         f1log->SetLineWidth(2);
      }


      if ( gdata->GetY()[0] > 0 )
         graph->SetMinimum(0);

      graph->SetMaximum( max(gdata->GetY()[0],max(graph->GetY()[0],graph->GetY()[graph->GetN()-1]))*1.2);
      graph->GetYaxis()->SetRangeUser(0,20);
      //graph->SetMinimum( -20 );
      graph->Draw("APE0");
      if ( fit.GetLogNormal() )
         f1log->Draw("Lsame");
      else
         pol1w->Draw("Lsame");
      //graph->GetFunction("pol1")->Draw("Lsame");
    if ( reference_values.size()+1<= 6 ) 
      graph->GetHistogram()->GetXaxis()->SetNdivisions(graph->GetN()+1);
    else
       graph->GetHistogram()->GetXaxis()->SetNdivisions(int(graph->GetN()/2)+1+200);
      //f2->Draw("same");

      gdata->Draw("PE0");
      
      if ( ibin==0 ) {
         double xmin = fit.GetLogNormal() ? 0.36 : 0.45;
         TLegend legend(xmin,0.19,0.96,0.47,"","NDC");
         //legend.SetNColumns(3);
         legend.SetFillStyle(0);
         legend.SetBorderSize(0);
         legend.SetTextSize(0.05);
         legend.AddEntry(data,"Data","E0P");
         legend.AddEntry(graph,"Templates","PE0");
         if ( fit.GetLogNormal() ) {
            legend.AddEntry(graph->GetFunction("pol1"),"Linear log(model)","L");
            legend.AddEntry(f1log,"#scale[0.9]{Linearized model #scale[0.7]{(unused)}}","L");
         }
         else {
            legend.AddEntry(graph->GetFunction("pol1"),"Linearized model","L");
            legend.AddEntry(pol1w,"Weighted fit #scale[0.7]{(unused)}","L");
         }
         legend.DrawClone();
      }

      TLatex text;
      text.SetNDC();
      text.SetTextFont(42);
      text.SetTextAlign(11);
      text.SetTextSize(0.05);
      
      text.SetTextSize(0.04);
      //text.DrawLatex(0.20,0.93,Form("%3.1f_{ }<_{ }|y|_{ }<_{ }%3.1f",input_table["ylow"][ibin],input_table["yhigh"][ibin]));
      TString infotext = "_{ }<_{ }" + xaxistitle + "_{ }<_{ }";
      infotext.Prepend(Form("%3.0f", bins[ibin]));
      infotext.Append(Form("%3.0f", bins[ibin+1]));
      infotext.Append("_{}");
      text.DrawLatex(0.20,0.93, infotext);

      c1.Print(ps_name);
      //if ( fit.GetLogNormal() )
      //   c1.Print( Form("plots/LTFlog_bin_%02d.pdf",ibin));
      //else
      //   c1.Print( Form("plots/LTF_bin_%02d.pdf",ibin));
      
   }
   

   // ---------------------------------------------- //
   //   chisq plot
   // ---------------------------------------------- //
    TGraph* gChi2 = new TGraph();
    int ndf = (fit.Dt.rows()-(fit.M.cols()-1));
    for ( int itmpl = 0 ; itmpl<fit.chisq_y.size() ; itmpl++ )  {
       gChi2->SetPoint(itmpl, reference_values[itmpl], fit.chisq_y(itmpl)/ndf);
    }
    gChi2->Fit("pol2","QW");
    gChi2->SetMarkerStyle(20);
    
    TGraph* gChi2LTF = new TGraph();
    gChi2LTF->SetPoint(0, fit.ahat(0), fit.chisq/ndf);
    gChi2LTF->SetMarkerStyle(29);
    gChi2LTF->SetMarkerSize(3.1);
    gChi2LTF->SetMarkerColor(kViolet+2);

    TGraph* gChi2chk = new TGraph();
    gChi2chk->SetPoint(0, fit.achk(0), fit.achk_chisq/ndf);
    gChi2chk->SetMarkerSize(2.2);
    gChi2chk->SetMarkerStyle(24);
    gChi2chk->SetMarkerColor(kRed);
    
    //gChi2LTF->Print("all");
    
    //gChi2->SetTitle(";#alpha_{0} [unit];#chi^{2}/ndf");
    gChi2->SetTitle((";"+referencename+";#chi^{2}/ndf").c_str());
    gChi2->SetMinimum(0.0);
    gChi2->SetMaximum(2.5); 
    gChi2->SetMinimum(0.75); // use for CMS jet fits
    gChi2->SetMaximum(1.00); // use for CMS jet fits 

    gChi2->Draw("ap");
    if ( reference_values.size()+1<= 8 ) 
       gChi2->GetHistogram()->SetNdivisions(reference_values.size()+1+300,"X");
    else
       gChi2->GetHistogram()->SetNdivisions(int(reference_values.size()/2)+1+200,"X");
    
    TLine line;
    line.SetLineColor(920);
    line.SetLineStyle(3);
    line.DrawLine( 
       gChi2->GetHistogram()->GetXaxis()->GetXmin(),1.,
       gChi2->GetHistogram()->GetXaxis()->GetXmax(), 1);

    line.DrawLine( 
       fit.achk_chisq/ndf+1,1.,
       fit.achk_chisq/ndf+1,1.);

    gChi2chk->Draw("Psame");
    gChi2LTF->Draw("Psame");

    {
       TLegend legend(0.18,0.70,0.66,0.97,"","NDC");
       //legend.SetNColumns(3);
       legend.SetFillStyle(0);
       legend.SetBorderSize(0);
       legend.SetTextSize(0.045);
       legend.AddEntry(gChi2LTF,"#hat#chi^{2} of the Linear Template Fit","P");
       legend.AddEntry(gChi2,   "#chi^{2}_{#font[12]{j}} of the individual templates","P");
       legend.AddEntry(gChi2->GetFunction("pol2"),"Parabola","L");
       legend.AddEntry(gChi2chk,"Minimum of #chi^{2} parabola #scale[0.8]{(#check#chi^{2})}","P"); //  (#check#chi^{2})
       legend.DrawClone();
    }
       
    c1.Print(ps_name);
    //c1.Print( "plots/LTF_chi2.pdf");

    c1.Print( (string(ps_name)+"]").c_str() );
}


// __________________________________________________________________________________ //
//!
//!
//!  Plot a LiTeFit object using ROOT
//!
//!  The binning needs to be provided to the plotting function,
//!  since this is not included in LTF::LiTeFit
//! 
void LTF_ROOTTools::plotLiTeFit_2D(const LTF::LiTeFit& fit, const vector<double> bins )
{

   gStyle->SetOptStat(0);
   gSystem->mkdir("plots");
   auto& M = fit.M;
   if ( M.cols() != 3 ) {cout<<"Error! only 2-dim plotting implemented."<<endl;exit(1);}
   Eigen::VectorXd reference_values1 = M.col(1);
   Eigen::VectorXd reference_values2 = M.col(2);

   map<double,TH1D*> templates;
   set<double> r1,r2;
   for ( int iref = 0 ; iref<reference_values1.size() ; iref++ ) {
      //double ref = reference_values1(iref);
      templates[iref] = MakeHistogram(fit.Y.col(iref),bins);
      r1.insert(reference_values1(iref));
      r2.insert(reference_values2(iref));
   }
   
   if ( r1.size()<=1 ) {
      cout<<"ERROR! at least two distinct reference points for dimension 1 must be given"<<endl;
      exit(1);
   }
   if ( r2.size()<=1 ) {
      cout<<"ERROR! at least two distinct reference points for dimension 2 must be given"<<endl;
      exit(1);
   }

   double xmin = 2.* (*r1.begin())  - (*(++r1.begin()))*0.9999;
   double xmax = 2.* (*r1.rbegin()) - (*(++r1.rbegin()))*1.00001;
   double ymin = 2.* (*r2.begin())  - (*(++r2.begin()))*0.9999;
   double ymax = 2.* (*r2.rbegin()) - (*(++r2.rbegin()))*1.00001;
   
   TH1D* data    = MakeHistogram(fit.Dt,bins,fit.Vs);
   TH1D* TheoFit = MakeHistogram(fit.TheoFit,bins);
   
   TCanvas c1("c1","LTF plots",800,800);
   c1.SetRightMargin(0.05);
   c1.SetLeftMargin(0.15);
   c1.SetTopMargin(0.08);
   // c1.SetRightMargin(0.02);

   const char* ps_name = "LTF2D_plots.ps";
   c1.Print( (string(ps_name)+"[").c_str() );

   // ---------------------------------------------- //
   // main plot
   // ---------------------------------------------- //
   for ( int iref = 0 ; iref<reference_values1.size() ; iref++ ) {
      templates[iref]->SetLineWidth(2);
      if ( iref == 0 ) {
         templates[0]->SetTitle("Linear Template Fit;Observable [unit]; Value [unit]");
         templates[0]->SetLineColor(kRed+2);
         templates[0]->SetMinimum(0);
         templates[0]->SetMaximum(templates[0]->GetMaximum()*1.7);
         templates[0]->SetLineWidth(3);
         templates[0]->DrawClone("hist");
      }
      else if ( iref==reference_values1.size()-1) {
         templates[iref]->SetFillColorAlpha(kBlue,0.15);
         templates[iref]->SetLineColor(kBlue+2);
         templates[iref]->SetLineWidth(3);
         templates[iref]->Draw("histsame");
      }
      else {
         int color = iref+1;
         if (color >= 10 ) color=(color-10)*2+28;
         templates[iref]->SetLineColor(color);
         //templates[iref]->SetFillColorAlpha(color,0.08);
         templates[iref]->Draw("histsame");
      }
   }   
   templates[0]->SetFillColorAlpha(kRed,0.15);
   templates[0]->Draw("histsame");

   data->SetMarkerStyle(20);
   data->SetMarkerSize(1.4);
   data->SetLineColor(kBlack);
   data->Draw("e0same");

   TheoFit->SetLineColor(923);
   TheoFit->SetLineWidth(4);
   TheoFit->SetLineStyle(2);
   TheoFit->Draw("histsame");

   TLegend legend(0.18,0.70,0.94,0.92,"","NDC");
   legend.SetNColumns(3);
   legend.SetFillStyle(0);
   legend.SetBorderSize(0);
   legend.AddEntry(data,"Data","E0P");
   for ( int iref = 0 ; iref<reference_values1.size() ; iref++ ) {
      legend.AddEntry(templates[iref],Form("Tpl. #alpha_{0}=%5.1f, #alpha_{1}=%3.1f",reference_values1[iref],reference_values2[iref]),"FL");
   }
   legend.AddEntry(TheoFit,"Estimated best model","L");
   legend.Draw();

   c1.Print(ps_name);
   c1.Print("plots/LTF2D_plot.pdf");


   // ---------------------------------------------- //
   // print linear-functions in every bin
   // ---------------------------------------------- //
   c1.Clear();
   c1.SetRightMargin(0.02);
   c1.SetTopMargin(0.02);

   c1.SetLeftMargin(0.12);
   c1.SetBottomMargin(0.12);

   for ( int ibin = 0 ; ibin<fit.Dt.size() ; ibin++ ) {
      //TGraphErrors* data  = MakeTGraph(fit.ahat.row(0),fit.Dt.row(ibin));
      TGraph2DErrors* gdata = new TGraph2DErrors();
      gdata->SetPoint(0, fit.ahat(0), fit.ahat(1), fit.Dt(ibin));
      gdata->SetPointError(0,0,0,data->GetBinError(ibin+1));
      gdata->SetMarkerStyle(20);
      gdata->SetMarkerSize(2);


      TGraph2DErrors* graph = MakeTGraph2D(fit.M.col(1),fit.M.col(2),ibin,fit.Y,fit.SysY);
      TF2* f2 = new TF2(Form("Linear Template Fit (2-dim., bin %d)",ibin),"[0]+[1]*x+[2]*y",
                        xmin,xmax,ymin,ymax
         );
      graph->Fit(f2,"QW");



      TGraph2D* gtheo = new TGraph2D();
      gtheo->SetPoint(0, fit.ahat(0), fit.ahat(1), f2->Eval(fit.ahat(0), fit.ahat(1)));
      gtheo->SetMarkerStyle(20);
      gtheo->SetMarkerSize(0.4);
      gtheo->SetMarkerColor(kBlack);
      for ( int it = 0 ; it<fit.M.col(1).size(); it++ ) {
         gtheo->SetPoint(it+1, reference_values1(it), reference_values2(it),
                         f2->Eval(reference_values1(it), reference_values2(it)));
      }

      f2->SetLineColor(922);
      f2->SetLineStyle(3);
      f2->SetLineWidth(1);
      f2->SetMinimum(0);
      f2->Draw("");
      f2->GetHistogram()->GetXaxis()->SetNdivisions(r1.size()+1);
      f2->GetHistogram()->GetYaxis()->SetNdivisions(r2.size()+1);
      f2->GetHistogram()->GetXaxis()->SetTitleOffset(2.2);
      f2->GetHistogram()->GetYaxis()->SetTitleOffset(2.2);
      f2->GetHistogram()->GetZaxis()->SetTitleOffset(1.8);
      f2->GetHistogram()->SetTitle(";Reference value 0 (#alpha_{0}) [unit]  ;Reference value 1 (#alpha_{1}) [unit];Value [unit]");

      graph->SetMarkerStyle(25);
      //graph->SetMarkerSize(2);
      graph->SetMarkerColor(kRed+2);
      graph->SetLineColor(kRed+2);
      graph->SetTitle(";Reference value (#alpha) [unit];Value [unit]");

      //graph->SetMinimum(0);
      //graph->SetMaximum( max(gdata->GetY()[0],max(graph->GetY()[0],graph->GetY()[graph->GetN()-1]))*1.2);
      //graph->Draw("ap");
      //pol1w->Draw("Lsame");
      //graph->GetFunction("pol1")->Draw("Lsame");
      //graph->GetHistogram()->GetXaxis()->SetNdivisions(graph->GetN()+1);

      f2->Draw("surf1");
      gtheo->Draw("P0 same");
      graph->Draw("err p0 same");
      gdata->Draw("err P0 same");


      TF2* f2leg = (TF2*)f2->Clone("f2leg");
      f2leg->SetLineWidth(3);
      gtheo->SetMarkerStyle(24);
      gtheo->SetMarkerSize(0.4);
      gdata->SetMarkerStyle(24);
      graph->SetMarkerStyle(24);
      if ( ibin==1 ) {
         TLegend legend(0.42,0.18,0.80,0.43,"","NDC");
         //legend.SetNColumns(3);
         legend.SetFillStyle(0);
         legend.SetBorderSize(0);
         legend.SetTextSize(0.03);
         legend.AddEntry(graph,"Templates","PE0");
         legend.AddEntry(f2leg,"Linearized model","L");
         legend.AddEntry(gtheo,"Projections onto model","P");
         legend.AddEntry(gdata,"Data","E0P");
         legend.DrawClone();
      }

      TLatex text;
      text.SetNDC();
      text.SetTextFont(42);
      text.SetTextAlign(13);
      text.SetTextSize(0.05);
      //text.DrawLatex(0.75,0.25,Form("Bin %d",ibin));
      text.DrawLatex(0.02,0.97,Form("Bin %d",ibin));

      c1.Print(ps_name);
      c1.Print( Form("plots/LTF2D_bin_%02d.pdf",ibin));
      
   }
   
   c1.Print( (string(ps_name)+"]").c_str() );
   
}






// __________________________________________________________________________________ //
//!
//!
//!  Plot a LiTeFit object using ROOT
//!
//!  The binning needs to be provided to the plotting function,
//!  since this is not included in LTF::LiTeFit
//! 
void LTF_ROOTTools::plotLiTeFitPol2Test(const LTF::LiTeFit& fit, const vector<double>& bins, 
                 const string& yaxistitle,
                 const string& referencename,
                 const string& observablename) 
{
   gStyle->SetOptStat(0);
   gSystem->mkdir("plots");
   auto& M = fit.M;

   // sanity check
   if ( M.cols() != 2 ) {cout<<"Error! only 1-dim plotting implemented."<<endl;exit(1);}
   Eigen::VectorXd reference_values = M.col(1);
   
   map<double,TH1D*> templates;
   for ( int iref = 0 ; iref<reference_values.size() ; iref++ ) {
      templates[iref] = MakeHistogram(fit.Y.col(iref),bins);
   }

   TH1D* data    = MakeHistogram(fit.Dt,bins,fit.Vs);
   TH1D* TheoFit = MakeHistogram(fit.TheoFit,bins);
   
   TCanvas c1("c1","LTF plots",800,800);
   c1.SetRightMargin(0.05);
   c1.SetLeftMargin(0.15);
   c1.SetTopMargin(0.08);

   const char* ps_name = fit.GetLogNormal() ?
      "LTFlog_plots.ps" :
      "LTF_plots.ps";
   c1.Print( (string(ps_name)+"[").c_str() );

   // ---------------------------------------------- //
   // main plot
   // ---------------------------------------------- //
   for ( int iref = 0 ; iref<reference_values.size() ; iref++ ) {
      templates[iref]->SetLineWidth(2);
      if ( iref == 0 ) {
         templates[0]->SetTitle(("Quadratic Template Fit;"+observablename+";"+yaxistitle).c_str());
         templates[0]->SetLineColor(kRed+1);
         if ( templates[0]->GetMaximum()>0 )templates[0]->SetMinimum(0);
         if ( !fit.GetLogNormal() ) 
            templates[0]->SetMaximum(templates[0]->GetMaximum()*1.7);
         // else 
         //    templates[0]->SetMaximum(templates[0]->GetMaximum()*3.); 
         templates[0]->SetLineWidth(3);
         templates[0]->DrawClone("hist");
      }
      else if ( iref==reference_values.size()-1) {
         if ( templates[0]->GetMaximum()>0 ) templates[iref]->SetFillColorAlpha(kBlue,0.15);
         templates[iref]->SetLineColor(kBlue+2);
         templates[iref]->SetLineWidth(3);
         templates[iref]->Draw("histsame");
      }
      else {
         templates[iref]->SetLineColor(iref+2);
         templates[iref]->SetLineWidth(2);
         templates[iref]->Draw("histsame");
      }
   }   
   if ( templates[0]->GetMaximum()>0 )templates[0]->SetFillColorAlpha(kRed,0.15);
   templates[0]->Draw("histsame");

   data->SetMarkerStyle(20);
   data->SetMarkerSize(1.4);
   data->SetLineColor(kBlack);
   data->Draw("e0same");

   TheoFit->SetLineColor(923);
   TheoFit->SetLineWidth(4);
   TheoFit->SetLineStyle(2);
   TheoFit->Draw("histsame");

   TLegend legend(0.18,0.70,0.94,0.92,"","NDC");
   legend.SetNColumns(3);
   legend.SetFillStyle(0);
   legend.SetBorderSize(0);
   legend.AddEntry(data,"Data","E0P");
   for ( int iref = 0 ; iref<reference_values.size() ; iref++ ) {
      legend.AddEntry(templates[iref],Form("Template #alpha=%6.2f",reference_values[iref]),"FL");
   }
   legend.AddEntry(TheoFit,"Estimated best model","L");
   legend.Draw();

   c1.Print(ps_name);
   c1.Print("plots/LTF_plot.pdf");


   // ---------------------------------------------- //
   // print linear-functions in every bin
   // ---------------------------------------------- //
   c1.Clear();
   c1.SetRightMargin(0.02);
   c1.SetTopMargin(0.02);
   c1.SetLeftMargin(0.16);
   c1.SetBottomMargin(0.16);

   gPad->SetTicky(1);

   gStyle->SetLabelSize(0.05,"XYZ");
   gStyle->SetTitleSize(0.05,"XYZ");
   gStyle->SetTitleOffset(1.1,"X");
   gStyle->SetTitleOffset(1.6,"Y");
   gStyle->SetMarkerSize(2);


   map < string, vector<double> > input_table = read_input_table2("data/CMS_data.txt",32);

   for ( int ibin = 0 ; ibin<fit.Dt.size() ; ibin++ ) {
      //TGraphErrors* data  = MakeTGraph(fit.ahat.row(0),fit.Dt.row(ibin));
      TGraphErrors* gdata = new TGraphErrors();
      gdata->SetPoint(0,fit.ahat(0),fit.Dt(ibin));
      gdata->SetPointError(0,0,data->GetBinError(ibin+1));
      gdata->SetMarkerStyle(20);

      TGraphErrors* graph = MakeTGraph(fit.M.col(1),ibin,fit.Y,fit.SysY);
      graph->Fit("pol1","QW"); // in this function it is "QW" (unweighted
      TF1* pol1 = (TF1*)graph->GetFunction("pol1")->Clone("pol1");
      pol1->SetLineColor(kRed);
      //pol1w->SetLineStyle(3);
      pol1->SetLineWidth(2);

      TF1* f2= new TF1("f2","[0]+[1]*pow(x,2)", -FLT_MIN,FLT_MAX );
      graph->Fit(f2,"QW");
      TF1* f1= new TF1("f1","[0]+[1]*pow(x,[2])", -FLT_MIN,FLT_MAX );
      

      graph->Fit("pol1","QW"); // "W": Ignore all point errors when fitting a TGraphErrors 
      bool UseLTFOutput = true;
      if ( UseLTFOutput ) {
         Eigen::VectorXd ltfpol1param =(fit.Y*fit.Mc().transpose()).row(ibin);
         //cout<<"M+*Y:"<<endl<< ltfpol1param <<endl;
         graph->GetFunction("pol1")->SetParameter(0,ltfpol1param(0));
         graph->GetFunction("pol1")->SetParameter(1,ltfpol1param(1));
         f1->SetParameter(0,ltfpol1param(0));
         f1->SetParameter(1,ltfpol1param(1));
         f1->SetParameter(2,fit.Gamma[0]);
      }
      if ( fit.Gamma[0]!=1 )   cout<<"Warning! Plotting with gamma factor !=1 not correctly implmeneted!"<<endl;

      graph->GetFunction("pol1")->SetLineWidth(3);
      
      graph->Fit("pol2","QW");
      graph->GetFunction("pol2")->SetLineWidth(83);

      graph->SetMarkerStyle(47);
      graph->SetMarkerColor(kRed+2);
      graph->SetLineColor(kRed+2);
      if ( fit.GetLogNormal() ) 
         graph->SetTitle((";"+referencename+";log("+yaxistitle+")").c_str()); // log(value/unit)
      else
         graph->SetTitle((";"+referencename+";"+yaxistitle).c_str());


      TF1* f1log = NULL;
      if ( fit.GetLogNormal() ) {
         //f1log = new TF1("pol1log","[0]+[1]*exp(x)", -FLT_MIN,FLT_MAX );
         f1log = new TF1("pol1log","log([0]+[1]*pow(x,1))", -FLT_MIN,FLT_MAX );

         TGraph* gexp = new TGraph();
         for ( int i = 0 ; i<graph->GetN() ; i++ ) 
            gexp->SetPoint(i,graph->GetX()[i], exp(graph->GetY()[i]));
         gexp->Fit("pol1","QW");

         Eigen::VectorXd ltfpol1param =(fit.Y*fit.Mc().transpose()).row(ibin);
         graph->GetFunction("pol1")->SetParameter(0,ltfpol1param(0));

         f1log->SetParameter(0, gexp->GetFunction("pol1")->GetParameter(0));
         f1log->SetParameter(1, gexp->GetFunction("pol1")->GetParameter(1));
         if ( fit.Gamma[0]!=1 )   cout<<"Warning! Plotting with gamma factor !=1 not correctly implmeneted!"<<endl;

         // //f1 = new TF1("pol1","[0]+[1]*x", -FLT_MIN,FLT_MAX );
         // graph->Fit(f1log,"QW"); // "W": Ignore all point errors when fitting a TGraphErrors 
         f1log->SetLineColor(kBlue+1);
         f1log->SetLineStyle(7);
         f1log->SetLineWidth(2);
      }


      if ( gdata->GetY()[0] > 0 )
         graph->SetMinimum(0);

      graph->SetMaximum( max(gdata->GetY()[0],max(graph->GetY()[0],graph->GetY()[graph->GetN()-1]))*1.2);
      graph->Fit("pol2","QW");
      graph->GetFunction("pol2")->SetLineColor(kBlue);;

      TF1* tang = new TF1("tang","[0]+[1]*x", -FLT_MIN,FLT_MAX );
      double f0  = graph->GetFunction("pol2")->Eval(gdata->GetX()[0]);
      double fp0 = graph->GetFunction("pol2")->GetParameter(1) + 2.*graph->GetFunction("pol2")->GetParameter(2)*gdata->GetX()[0];
      tang->SetParameter(0, f0 - fp0*gdata->GetX()[0] );
      tang->SetParameter(1, fp0);
      tang->SetLineWidth(2);
      tang->SetLineColor(kTeal+4);
      tang->SetLineStyle(7);
      
      pol1->SetLineStyle(3);
      graph->Draw("ap");
      if ( fit.GetLogNormal() )
         f1log->Draw("Lsame");
      else
         pol1->Draw("Lsame");
      //graph->GetFunction("pol1")->Draw("Lsame");
      if ( reference_values.size()+1<= 8 ) 
         graph->GetHistogram()->GetXaxis()->SetNdivisions(graph->GetN()+1);
      else
         graph->GetHistogram()->GetXaxis()->SetNdivisions(int(graph->GetN()/2)+1+200);
      //f2->Draw("same");

      tang->Draw("Lsame");
      gdata->Draw("P");
      
      if ( ibin==0 ) {
         //double xmin = fit.GetLogNormal() ? 0.36 : 0.45;
         //TLegend legend(xmin,0.19,0.96,0.47,"","NDC");
         TLegend legend(0.19,0.19,0.56,0.47,"","NDC"); 
         //legend.SetNColumns(3);
         legend.SetFillStyle(0);
         legend.SetBorderSize(0);
         legend.SetTextSize(0.045);
         legend.AddEntry(data,"Data","E0P");
         legend.AddEntry(graph,"Templates","PE0");
         if ( fit.GetLogNormal() ) {
            legend.AddEntry(graph->GetFunction("pol1"),"Linear log(model)","L");
            legend.AddEntry(f1log,"#scale[0.9]{Linearized model #scale[0.7]{(unused)}}","L");
         }
         else {
            legend.AddEntry(graph->GetFunction("pol2"),"Second-order model","L");
            legend.AddEntry(tang,"Linear model","L");
            legend.AddEntry(pol1,"Linear Template Fit","L");
         }
         legend.DrawClone();
      }

      TLatex text;
      text.SetNDC();
      text.SetTextFont(42);
      text.SetTextAlign(11);
      text.SetTextSize(0.045);
      //text.DrawLatex(0.20,0.20,Form("Bin %d",ibin));

      text.SetTextAlign(31);
      text.DrawLatex(0.955,0.20,Form("Bin %d",ibin));

      // text.DrawLatex(0.20,0.30,"CMS inclusive jets");
      // text.SetTextSize(0.04);
      // text.DrawLatex(0.20,0.25,Form("%3.1f_{ }<_{ }|y|_{ }<_{ }%3.1f",input_table["ylow"][ibin],input_table["yhigh"][ibin]));
      // text.DrawLatex(0.20,0.20,Form("%3.0f_{ }<_{ }p_{T}_{ }<_{ }%3.0f_{ }GeV",input_table["ptlow"][ibin],input_table["pthigh"][ibin]));
      
      // text.SetTextSize(0.04);
      // text.DrawLatex(0.20,0.93,Form("%3.1f_{ }<_{ }|y|_{ }<_{ }%3.1f",input_table["ylow"][ibin],input_table["yhigh"][ibin]));
      // text.DrawLatex(0.20,0.88,Form("%3.0f_{ }<_{ }p_{T}_{ }<_{ }%3.0f_{ }GeV",input_table["ptlow"][ibin],input_table["pthigh"][ibin]));

      cout<<input_table["ylow"][ibin]<<"\t"
          <<input_table["yhigh"][ibin]<<"\t"
          <<input_table["ptlow"][ibin]<<"\t"
          <<input_table["pthigh"][ibin]<<endl;

      c1.Print(ps_name);
      if ( fit.GetLogNormal() )
         c1.Print( Form("plots/LTFlog_bin_%02d.pdf",ibin));
      else
         c1.Print( Form("plots/LTF_bin_%02d.pdf",ibin));
      
   }
   

   // ---------------------------------------------- //
   //   chisq plot
   // ---------------------------------------------- //
    TGraph* gChi2 = new TGraph();
    int ndf = (fit.Dt.rows()-(fit.M.cols()-1));
    for ( int itmpl = 0 ; itmpl<fit.chisq_y.size() ; itmpl++ )  {
       gChi2->SetPoint(itmpl, reference_values[itmpl], fit.chisq_y(itmpl)/ndf);
    }
    gChi2->Fit("pol2","QW");
    gChi2->SetMarkerStyle(20);
    
    TGraph* gChi2LTF = new TGraph();
    gChi2LTF->SetPoint(0, fit.ahat(0), fit.chisq/ndf);
    gChi2LTF->SetMarkerStyle(29);
    gChi2LTF->SetMarkerSize(3.1);
    gChi2LTF->SetMarkerColor(kViolet+2);

    TGraph* gChi2chk = new TGraph();
    gChi2chk->SetPoint(0, fit.achk(0), fit.achk_chisq/ndf);
    gChi2chk->SetMarkerSize(2.2);
    gChi2chk->SetMarkerStyle(24);
    gChi2chk->SetMarkerColor(kRed);
    
    //gChi2LTF->Print("all");
    
    //gChi2->SetTitle(";#alpha_{0} [unit];#chi^{2}/ndf");
    gChi2->SetTitle((";"+referencename+";#chi^{2}/ndf").c_str());
    gChi2->SetMinimum(0.);
    gChi2->SetMaximum(20.2);
    gChi2->Draw("apc");
    if ( reference_values.size()+1<= 8 ) 
       gChi2->GetHistogram()->SetNdivisions(reference_values.size()+1+500,"X");
    else
       gChi2->GetHistogram()->SetNdivisions(int(reference_values.size()/2)+1+500,"X");
    
    TLine line;
    line.SetLineColor(920);
    line.SetLineStyle(3);
    line.DrawLine( 
       gChi2->GetHistogram()->GetXaxis()->GetXmin(),1.,
       gChi2->GetHistogram()->GetXaxis()->GetXmax(), 1);

    line.DrawLine( 
       fit.achk_chisq/ndf+1,1.,
       fit.achk_chisq/ndf+1,1.);

    gChi2chk->Draw("Psame");
    gChi2LTF->Draw("Psame");

    {
       TLegend legend(0.18,0.75,0.56,0.97,"","NDC");
       //legend.SetNColumns(3);
       legend.SetFillStyle(0);
       legend.SetBorderSize(0);
       legend.SetTextSize(0.045);
       legend.AddEntry(gChi2LTF,"#hat#chi^{2} of the Quadratic Template Fit","P");
       legend.AddEntry(gChi2,   "#chi^{2}_{#font[12]{j}} of the individual templates","PL");
       legend.AddEntry(gChi2->GetFunction("pol2"),"Parabola","L");
       //legend.AddEntry(gChi2chk,"Minimum of #chi^{2} parabola #scale[0.8]{(#check#chi^{2})}","P"); //  (#check#chi^{2})
       legend.DrawClone();
    }
       
    c1.Print(ps_name);
    c1.Print( "plots/LTF_chi2.pdf");

    c1.Print( (string(ps_name)+"]").c_str() );
   
}


double LTF_ROOTTools::makeErrorPlot(TCanvas& c, const char* ps_name, const char* title, const LTF::LiTeFit& fit, const vector<string> &uncertainties)
{
   bool useNuisanceParameter = true;
   int nPar = 1; //M.cols()-1; 
   double sum_error = 0; // this needs to be a vector in the case of more than 1 parameter
   if ( !useNuisanceParameter ) {
     for ( int i = 0 ; i<nPar ; i++ ) {
         TH1D* h  = new TH1D(title, title, uncertainties.size()+1, 0, uncertainties.size()+1);
         
         for ( const string &source: uncertainties ) {
            double error = 0;
            //if (source.find("stat.")!= std::string::npos ) error = std::sqrt(fabs(fit.Vsource.find(source)->second(1,1)));
            //else error = std::sqrt(fabs(fit.Vsource.find("unfolding_error_"+variable+"_direct_envelope_"+source+"__1up")->second(i,i)));
	    
            h->Fill(source.c_str(), error);
            sum_error += pow(error,2);
         }
         h->SetBinContent(h->GetNbinsX(), std::sqrt(sum_error));
         h->GetXaxis()->SetBinLabel(h->GetNbinsX(), "Total unc.");
         h->SetBarWidth(0.85);
         h->GetYaxis()->SetLabelSize(0.03);
         h->GetYaxis()->SetTitle("Uncertainty [GeV]");
         h->GetXaxis()->SetLabelSize(0.02);
         h->GetXaxis()->SetTickLength(0);
         h->Draw("hbar");
      }
   }
   else {
      for ( int i = 0 ; i<nPar ; i++ ) {
         c.Divide(2,1);
         c.cd(1);
         TH1D* h  = new TH1D(title, title, uncertainties.size()+1, 0, uncertainties.size()+1);

	 for ( auto& [name,val]: fit.Vsource ) cout<<name<<" "<<val(1,1)<<endl;
	 
	 for ( const string &source: uncertainties ) {
	   double error = 0;
	   if (source.find("STAT_DATA")!= std::string::npos ) error = std::sqrt(fabs(fit.Vsource.find(source)->second(1,1)));
	   else if (source.find("STAT_MC")!= std::string::npos ) error = std::sqrt(fabs(fit.Vsource.find(source)->second(1,1)));
	   else	   error = fabs(fit.DeltaSys.find(source)->second(i));

	   h->Fill(source.c_str(), error);
	   sum_error += pow(error,2);
	 }
	 
	 h->SetBinContent(h->GetNbinsX(), std::sqrt(sum_error));
         h->GetXaxis()->SetBinLabel(h->GetNbinsX(), "Total unc.");
         h->SetBarOffset(0.1);
         h->SetBarWidth(0.8);
         h->GetYaxis()->SetLabelSize(0.03);
         h->GetYaxis()->SetTitle("Uncertainty [GeV]");
         h->GetXaxis()->SetLabelSize(0.02);
         h->GetXaxis()->SetTickLength(0);
         h->Draw("hbar");
         
         c.cd(2)->SetLeftMargin(0.15);
         TH1D* h1  = new TH1D("NP", "NP", uncertainties.size()+1, 0, uncertainties.size()+1);
         TGraphErrors* g = new TGraphErrors(uncertainties.size());
         for ( long unsigned int j = 0; j < uncertainties.size(); j++ ) {
	    h1->Fill(uncertainties[j].c_str(), fit.map_nuisance.find(uncertainties[j])->second.first);
            h1->SetBinError(j,fit.map_nuisance.find(uncertainties[j])->second.second);
            g->SetPoint(j, fit.map_nuisance.find(uncertainties[j])->second.first, j+0.5);
            g->SetPointError(j, fit.map_nuisance.find(uncertainties[j])->second.second, 0);
         }
         h1->SetBinContent(h1->GetNbinsX(), 0.);
         h1->GetXaxis()->SetBinLabel(h1->GetNbinsX(), "");
         gStyle->SetHistMinimumZero();
         h1->SetBarOffset(0.95);
         h1->SetBarWidth(0);
         h1->SetLineColor(10);
         h1->SetFillColor(10);
         h1->SetLineColorAlpha(10,0);

         h1->GetYaxis()->SetLabelSize(0.03);
         h1->GetYaxis()->SetTitle("Nuisance parameter");
         h1->GetXaxis()->SetLabelSize(0.02);
         h1->GetXaxis()->SetTickLength(0);
         g->SetMarkerSize(1);
         g->SetMarkerStyle(20);
         g->SetMarkerColor(kBlue);
         h1->Draw("hbar e");
         g->Draw("same PE");
         
         TLine* l = new TLine(0, 0, 0, h1->GetNbinsX());
         l->SetLineStyle(3);
         l->Draw("same");

         c.Print(ps_name);
         c.Clear();
         h->Delete();
         h1->Delete();
         g->Delete();
         l->Delete();
      }
   }
   return std::sqrt(sum_error);
}

