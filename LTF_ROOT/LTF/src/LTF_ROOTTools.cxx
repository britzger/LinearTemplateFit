#include "LTF/LTF.h"
#include "LTF/LTF_ROOTTools.h"

#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMatrixD.h>
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
TGraphErrors* LTF_ROOTTools::MakeTGraph(const TVectorD& xvalues, int ibin, const TMatrixD& Y,
                         const std::map<std::string,TMatrixD >& VSysY ) 
{

   TGraphErrors* graph = new TGraphErrors();
   for ( int i = 0 ; i<xvalues.GetNoElements() ; i++ ) {
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
TGraph2DErrors* LTF_ROOTTools::MakeTGraph2D(const TVectorD& xvalues, const TVectorD& yvalues,
                                            int ibin, const TMatrixD& Y,
                                            const std::map<std::string,TMatrixD >& VSysY ) 
{

   TGraph2DErrors* graph = new TGraph2DErrors();
   for ( int i = 0 ; i<xvalues.GetNoElements() ; i++ ) {
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
TH1D* LTF_ROOTTools::MakeHistogram(const TVectorD& values, vector<double> bins, const std::vector<std::pair<std::string,TMatrixDSym > >& V ) 
{
   TH1D* hist = bins.empty() ?
      new TH1D("hist","hist",values.GetNoElements(),0,values.GetNoElements() ) :
      new TH1D("hist","hist",bins.size()-1, &bins[0]);
   if ( bins.size() && int(values.GetNoElements()+1) != int(bins.size()) ) {cout<<"ERROR! binning and number of entries does not fit!"<<endl;exit(1);}
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
                 const string& yaxistitle,
                 const string& referencename,
                 const string& observablename ) 
{
   gStyle->SetOptStat(0);
   gSystem->mkdir("plots");
   auto& M = fit.M;
   
   // sanity check
   if ( M.GetNcols() != 2 ) {cout<<"Error! only 1-dim plotting implemented."<<endl;exit(1);}
   TVectorD reference_values = TMatrixDColumn_const(M,1);
   //TH1D* hist = new TH1D("hist","hist",bins.size()-1, &bins[0]);
   
   map<double,TH1D*> templates;
   for ( int iref = 0 ; iref<reference_values.GetNoElements() ; iref++ ) {
      //double ref = reference_values(iref);
      templates[iref] = MakeHistogram(TMatrixDColumn_const(fit.Y,iref),bins);
   }

   TH1D* data    = MakeHistogram(fit.Dt,bins,fit.Vs);
   TH1D* TheoFit = MakeHistogram(fit.TheoFit,bins);
   
   TCanvas c1("c1","LTF plots",800,800);
   c1.SetRightMargin(0.05);
   c1.SetLeftMargin(0.15);
   c1.SetTopMargin(0.08);
   // c1.SetRightMargin(0.02);

   const char* ps_name = fit.GetLogNormal() ?
      "plots/LTFlog_plots.ps" :
      "plots/LTF_plots.ps";
   c1.Print( (string(ps_name)+"[").c_str() );
   bool doLogPlot = true;
   if ( doLogPlot) c1.cd()->SetLogy();
   // ---------------------------------------------- //
   // main plot
   // ---------------------------------------------- //
   for ( int iref = 0 ; iref<reference_values.GetNoElements() ; iref++ ) {
      templates[iref]->SetLineWidth(2);
      if ( iref == 0 ) {
         templates[0]->SetTitle(("Linear Template Fit;"+observablename+";"+yaxistitle).c_str());
         templates[0]->SetLineColor(kRed+1);
         if ( templates[0]->GetMaximum()>0 && doLogPlot ){
            templates[0]->SetMinimum(0.0000001); //was 0 johannes
            templates[0]->SetMaximum(templates[0]->GetMaximum()*60);
         }
         else if ( templates[0]->GetMaximum()>0 && !doLogPlot ){
            templates[0]->SetMinimum(0);
            if ( !fit.GetLogNormal() ) 
               templates[0]->SetMaximum(templates[0]->GetMaximum()*1.7);
         }

         // else 
         //    templates[0]->SetMaximum(templates[0]->GetMaximum()*3.); 
         templates[0]->SetLineWidth(3);
         templates[0]->DrawClone("hist");
      }
      else if ( iref==reference_values.GetNoElements()-1) {
         //if ( templates[0]->GetMaximum()>0 ) templates[iref]->SetFillColorAlpha(kBlue,0.15); //johannes
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
   //if ( templates[0]->GetMaximum()>0 )templates[0]->SetFillColorAlpha(kRed,0.15);//johannes
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
   for ( int iref = 0 ; iref<reference_values.GetNoElements() ; iref++ ) {
      legend.AddEntry(templates[iref],Form("Template #alpha=%6.2f",reference_values[iref]),"FL");
   }
   legend.AddEntry(TheoFit,"Estimated best model","L");
   legend.Draw();

   c1.Print(ps_name);
   c1.Print("plots/LTF_plot.pdf");

   // ---------------------------------------------- //
   // print relative size of all errors //Johannes
   // ---------------------------------------------- //
   c1.Clear();
   c1.SetLogy(0);
   // set margins?

   // map for all histos: Make one histo for each summary
   
   TH1D* unc_single  = new TH1D("Uncertainty_single","Uncertainty_single",fit.Vsource.size(), 0, fit.Vsource.size());
   //TH1D* unc_summray = new TH1D("Uncertainty_summary","Uncertainty_summary",m_unc_summary.size(), 0, m_unc_summary.size());

   //std::map<std::string, std::double> m_unc_summary{};
   string lepton_uncertainties[] = {"EG_RESOLUTION_ALL", 
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
   
   TH1D* unc_single  = new TH1D("lepton_uncertainty","lepton_uncertainty",lepton_uncertainties.size(), 0, lepton_uncertainties.size());


   for ( const string &source: lepton_uncertainties ) {
      double test = std::sqrt(fit.Vsource.find("unfolding_error_mbl_selected_direct_envelope_"+source+"__1up")->second(0,0));
      cout<<source<<" "<<test<<endl;
   }


   int nPar = fit.M.GetNcols()-1;
   //for ( int i = 0 ; i<nPar ; i++ ) {
      int bin = 1;
      for ( auto& [name,V] : fit.Vsource ) {
         //printf("  myprintout  +/- % 8.6f (%s)\n", std::sqrt(V(i,i)), name.c_str());
         std::string prefix = "unfolding_error_mbl_selected_direct_envelope_";
         std::string suffix = "__1up";
         std::string title = name.substr(prefix.length());
         title.substr(0, title.length()-suffix.length());
         cout<<title<<endl;

         unc_single->SetBinContent(bin, std::sqrt(V(0,0))); //std::sqrt(V(i,i)));
         unc_single->GetXaxis()->SetBinLabel(bin, title.c_str());
         bin++;
      }
      double error = fit.Vsource.find("unfolding_error_mbl_selected_direct_envelope_WMASS_VAR_signal__1up")->second(0,0);
      cout<<"WMASS_VAR_signal__1up "<<std::sqrt(error)<<endl;

   unc_single->Draw("B"); //Draw car chart
   c1.Print(ps_name);


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

   const TMatrixD Mc_T = fit.Mc().T();
   for ( int ibin = 0 ; ibin<fit.Dt.GetNoElements() ; ibin++ ) {
   //for ( int ibin = 0 ; ibin<0 ; ibin++ ) {
   //for ( int ibin = 0 ; ibin<6 ; ibin++ ) {
      //TGraphErrors* data  = MakeTGraph(fit.ahat.row(0),fit.Dt.row(ibin));
      TGraphErrors* gdata = new TGraphErrors();
      gdata->SetPoint(0,fit.ahat(0),fit.Dt(ibin));
      gdata->SetPointError(0,0,data->GetBinError(ibin+1));
      gdata->SetMarkerStyle(20);

      TGraphErrors* graph = MakeTGraph(TMatrixDColumn_const(fit.M,1),ibin,fit.Y,fit.SysY);
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
         const TVectorD ltfpol1param =TMatrixDRow_const(fit.Y*Mc_T,ibin);
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

         const TVectorD ltfpol1param = TMatrixDRow_const(fit.Y*Mc_T,ibin);
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

      graph->Draw("ap");
      if ( fit.GetLogNormal() )
         f1log->Draw("Lsame");
      else
         pol1w->Draw("Lsame");
      //graph->GetFunction("pol1")->Draw("Lsame");
    if ( reference_values.GetNoElements()+1<= 6 ) 
      graph->GetHistogram()->GetXaxis()->SetNdivisions(graph->GetN()+1);
    else
       graph->GetHistogram()->GetXaxis()->SetNdivisions(int(graph->GetN()/2)+1+200);
      //f2->Draw("same");

      gdata->Draw("P");
      
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
      //text.DrawLatex(0.20,0.20,Form("Bin %d",ibin));
      // text.DrawLatex(0.20,0.30,"CMS inclusive jets");
      // text.SetTextSize(0.04);
      // text.DrawLatex(0.20,0.25,Form("%3.1f_{ }<_{ }|y|_{ }<_{ }%3.1f",input_table["ylow"][ibin],input_table["yhigh"][ibin]));
      // text.DrawLatex(0.20,0.20,Form("%3.0f_{ }<_{ }p_{T}_{ }<_{ }%3.0f_{ }GeV",input_table["ptlow"][ibin],input_table["pthigh"][ibin]));
      
      text.SetTextSize(0.04);
      //text.DrawLatex(0.20,0.93,Form("%3.1f_{ }<_{ }|y|_{ }<_{ }%3.1f",input_table["ylow"][ibin],input_table["yhigh"][ibin]));
      TString infotext = "_{ }<_{ }" + observablename + "_{ }<_{ }";
      infotext.Prepend(Form("%3.0f", bins[ibin]));
      infotext.Append(Form("%3.0f", bins[ibin+1]));
      infotext.Append("_{}");
      text.DrawLatex(0.20,0.93, infotext);

/*
      cout<<input_table["ylow"][ibin]<<"\t"
          <<input_table["yhigh"][ibin]<<"\t"
          <<input_table["ptlow"][ibin]<<"\t"
          <<input_table["pthigh"][ibin]<<endl;
*/

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
    int ndf = (fit.Dt.GetNrows()-(fit.M.GetNcols()-1));
    for ( int itmpl = 0 ; itmpl<fit.chisq_y.GetNoElements() ; itmpl++ )  {
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
    if ( reference_values.GetNoElements()+1<= 8 ) 
       gChi2->GetHistogram()->SetNdivisions(reference_values.GetNoElements()+1+300,"X");
    else
       gChi2->GetHistogram()->SetNdivisions(int(reference_values.GetNoElements()/2)+1+200,"X");
    
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
    c1.Print( "plots/LTF_chi2.pdf");

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
   if ( M.GetNcols() != 3 ) {cout<<"Error! only 2-dim plotting implemented."<<endl;exit(1);}
   TVectorD reference_values1 = TMatrixDColumn_const(M,1);
   TVectorD reference_values2 = TMatrixDColumn_const(M,2);

   map<double,TH1D*> templates;
   set<double> r1,r2;
   for ( int iref = 0 ; iref<reference_values1.GetNoElements() ; iref++ ) {
      //double ref = reference_values1(iref);
      templates[iref] = MakeHistogram(TMatrixDColumn_const(fit.Y,iref),bins);
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
   for ( int iref = 0 ; iref<reference_values1.GetNoElements() ; iref++ ) {
      templates[iref]->SetLineWidth(2);
      if ( iref == 0 ) {
         templates[0]->SetTitle("Linear Template Fit;Observable [unit]; Value [unit]");
         templates[0]->SetLineColor(kRed+2);
         templates[0]->SetMinimum(0);
         templates[0]->SetMaximum(templates[0]->GetMaximum()*1.7);
         templates[0]->SetLineWidth(3);
         templates[0]->DrawClone("hist");
      }
      else if ( iref==reference_values1.GetNoElements()-1) {
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
   for ( int iref = 0 ; iref<reference_values1.GetNoElements() ; iref++ ) {
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

   for ( int ibin = 0 ; ibin<fit.Dt.GetNoElements() ; ibin++ ) {
      //TGraphErrors* data  = MakeTGraph(fit.ahat.row(0),fit.Dt.row(ibin));
      TGraph2DErrors* gdata = new TGraph2DErrors();
      gdata->SetPoint(0, fit.ahat(0), fit.ahat(1), fit.Dt(ibin));
      gdata->SetPointError(0,0,0,data->GetBinError(ibin+1));
      gdata->SetMarkerStyle(20);
      gdata->SetMarkerSize(2);


      TGraph2DErrors* graph = MakeTGraph2D(TMatrixDColumn_const(fit.M,1),TMatrixDColumn_const(fit.M,2),ibin,fit.Y,fit.SysY);
      TF2* f2 = new TF2(Form("Linear Template Fit (2-dim., bin %d)",ibin),"[0]+[1]*x+[2]*y",
                        xmin,xmax,ymin,ymax
         );
      graph->Fit(f2,"QW");



      TGraph2D* gtheo = new TGraph2D();
      gtheo->SetPoint(0, fit.ahat(0), fit.ahat(1), f2->Eval(fit.ahat(0), fit.ahat(1)));
      gtheo->SetMarkerStyle(20);
      gtheo->SetMarkerSize(0.4);
      gtheo->SetMarkerColor(kBlack);
      for ( int it = 0 ; it<TVectorD(TMatrixDColumn_const(fit.M,1)).GetNoElements(); it++ ) {
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
   if ( M.GetNcols() != 2 ) {cout<<"Error! only 1-dim plotting implemented."<<endl;exit(1);}
   TVectorD reference_values = TMatrixDColumn_const(M,1);
   //TH1D* hist = new TH1D("hist","hist",bins.size()-1, &bins[0]);
   
   map<double,TH1D*> templates;
   for ( int iref = 0 ; iref<reference_values.GetNoElements() ; iref++ ) {
      //double ref = reference_values(iref);
      templates[iref] = MakeHistogram(TMatrixDColumn_const(fit.Y,iref),bins);
   }

   TH1D* data    = MakeHistogram(fit.Dt,bins,fit.Vs);
   TH1D* TheoFit = MakeHistogram(fit.TheoFit,bins);
   
   TCanvas c1("c1","LTF plots",800,800);
   c1.SetRightMargin(0.05);
   c1.SetLeftMargin(0.15);
   c1.SetTopMargin(0.08);
   // c1.SetRightMargin(0.02);

   const char* ps_name = fit.GetLogNormal() ?
      "LTFlog_plots.ps" :
      "LTF_plots.ps";
   c1.Print( (string(ps_name)+"[").c_str() );

   // ---------------------------------------------- //
   // main plot
   // ---------------------------------------------- //
   for ( int iref = 0 ; iref<reference_values.GetNoElements() ; iref++ ) {
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
      else if ( iref==reference_values.GetNoElements()-1) {
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
   for ( int iref = 0 ; iref<reference_values.GetNoElements() ; iref++ ) {
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


   //map < string, vector<double> > input_table = read_input_table2("data/CMS_data.txt",32);

   for ( int ibin = 0 ; ibin<fit.Dt.GetNoElements() ; ibin++ ) {
   //for ( int ibin = 0 ; ibin<0 ; ibin++ ) {
   //for ( int ibin = 0 ; ibin<6 ; ibin++ ) {
      //TGraphErrors* data  = MakeTGraph(fit.ahat.row(0),fit.Dt.row(ibin));
      TGraphErrors* gdata = new TGraphErrors();
      gdata->SetPoint(0,fit.ahat(0),fit.Dt(ibin));
      gdata->SetPointError(0,0,data->GetBinError(ibin+1));
      gdata->SetMarkerStyle(20);

      TGraphErrors* graph = MakeTGraph(TMatrixDColumn_const(fit.M,1),ibin,fit.Y,fit.SysY);
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
         const TVectorD ltfpol1param =TMatrixDRow_const(fit.Y*fit.Mc().T(),ibin);
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

         TVectorD ltfpol1param = TMatrixDRow_const(fit.Y*fit.Mc().T(),ibin);
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
      if ( reference_values.GetNoElements()+1<= 8 ) 
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

      
      // cout<<input_table["ylow"][ibin]<<"\t"
      //     <<input_table["yhigh"][ibin]<<"\t"
      //     <<input_table["ptlow"][ibin]<<"\t"
      //     <<input_table["pthigh"][ibin]<<endl;

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
    int ndf = (fit.Dt.GetNrows()-(fit.M.GetNcols()-1));
    for ( int itmpl = 0 ; itmpl<fit.chisq_y.GetNoElements() ; itmpl++ )  {
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
    if ( reference_values.GetNoElements()+1<= 8 ) 
       gChi2->GetHistogram()->SetNdivisions(reference_values.GetNoElements()+1+500,"X");
    else
       gChi2->GetHistogram()->SetNdivisions(int(reference_values.GetNoElements()/2)+1+500,"X");
    
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


