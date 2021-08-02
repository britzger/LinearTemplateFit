//#include "include/LTF/LTF.h"
#include "LTF/LTF.h"

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
#include <TLine.h>
#include <TRandom3.h>

using namespace std;

class LTF_ROOTTools {
public:

// __________________________________________________________________________________ //
//!
//! read_input_table2()
//!
//! read input data table from file 
//!
static
std::map < std::string, std::vector<double> > read_input_table2(std::string filename, int ncol ) {
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
//!  MakeHistogramPlot
//!
//!  make a histogram and fill it 
//!  with random events according to a gauss
//!  distribution around M
//!
static
TH1D* MakeHistogramPlot(int nEvents, int seed, double mean, double sigma, vector<double> bins ) {
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
//!  MakeHistogramPlot
//!  make a histogram from an Eigen::Vector for plotting purposes
//!
static
TH1D* MakeHistogramPlot(const Eigen::VectorXd& values, vector<double> bins ={}, const std::vector<std::pair<std::string,Eigen::MatrixXd > >& V = {} ) {
   TH1D* hist = bins.empty() ?
      new TH1D("hist","hist",values.size(),0,values.size() ) :
      new TH1D("hist","hist",bins.size()-1, &bins[0]);
   if ( bins.size() && values.size() != bins.size()-1 ) {cout<<"ERROR! binning and number of entries does not fit!"<<endl;exit(1);}
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
//!  Make a TGraph for plotting purposes
//!
static
TGraphErrors* MakeTGraphPlot(const Eigen::VectorXd& xvalues, int ibin, const Eigen::MatrixXd& Y,
                         const std::map<std::string,Eigen::MatrixXd >& VSysY = {}
   ) {
   TGraphErrors* graph = new TGraphErrors();
   for ( size_t i = 0 ; i<xvalues.size() ; i++ ) {
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
//!  MakeTGraph
//!
//!  make a TGraph for plotting purposes
//!
static
TGraphErrors* MakeTGraph(const Eigen::VectorXd& xvalues, int ibin, const Eigen::MatrixXd& Y,
                         const std::map<std::string,Eigen::MatrixXd >& VSysY = {}
   ) {

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
//! MakeHistogram
//!
//! make a histogram and fill it with random events according to a gauss
//! distribution around M
//!
static
TH1D* MakeHistogram(int nEvents, int seed, double mean, double sigma, vector<double> bins ) {
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
static
TH1D* MakeHistogram(const Eigen::VectorXd& values, vector<double> bins ={}, const std::vector<std::pair<std::string,Eigen::MatrixXd > >& V = {} ) {
   TH1D* hist = bins.empty() ?
      new TH1D("hist","hist",values.size(),0,values.size() ) :
      new TH1D("hist","hist",bins.size()-1, &bins[0]);
   if ( bins.size() && values.size() != bins.size()-1 ) {cout<<"ERROR! binning and number of entries does not fit!"<<endl;exit(1);}
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
static
void plotLiTeFit(const LTF::LiTeFit& fit, const vector<double>& bins, 
                 const string& yaxistitle    = "value [unit]",
                 const string& referencename = "Reference value (#alpha) [unit]",
                 const string& observablename = "Observable [unit]"
   ) {
   gStyle->SetOptStat(0);
   gSystem->mkdir("plots");
   auto& M = fit.M;
   
   // sanity check
   if ( M.cols() != 2 ) {cout<<"Error! only 1-dim plotting implemented."<<endl;exit(1);}
   Eigen::VectorXd reference_values = M.col(1);
   //TH1D* hist = new TH1D("hist","hist",bins.size()-1, &bins[0]);
   
   map<double,TH1D*> templates;
   for ( int iref = 0 ; iref<reference_values.size() ; iref++ ) {
      //double ref = reference_values(iref);
      templates[iref] = MakeHistogramPlot(fit.Y.col(iref),bins);
   }

   TH1D* data    = MakeHistogramPlot(fit.Dt,bins,fit.Vs);
   TH1D* TheoFit = MakeHistogramPlot(fit.TheoFit,bins);
   
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
   for ( int iref = 0 ; iref<reference_values.size() ; iref++ ) {
      templates[iref]->SetLineWidth(2);
      if ( iref == 0 ) {
         templates[0]->SetTitle(("Linear Template Fit;"+observablename+";"+yaxistitle).c_str());
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
   //for ( int ibin = 0 ; ibin<0 ; ibin++ ) {
   //for ( int ibin = 0 ; ibin<6 ; ibin++ ) {
      //TGraphErrors* data  = MakeTGraph(fit.ahat.row(0),fit.Dt.row(ibin));
      TGraphErrors* gdata = new TGraphErrors();
      gdata->SetPoint(0,fit.ahat(0),fit.Dt(ibin));
      gdata->SetPointError(0,0,data->GetBinError(ibin+1));
      gdata->SetMarkerStyle(20);

      TGraphErrors* graph = MakeTGraphPlot(fit.M.col(1),ibin,fit.Y,fit.SysY);
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
      graph->SetMarkerStyle(47);
      graph->SetMarkerColor(kRed+2);
      graph->SetLineColor(kRed+2);
      if ( fit.GetLogNormal() ) 
         graph->SetTitle((";"+referencename+";log("+yaxistitle+")").c_str()); // log(value/unit)
      else
         graph->SetTitle((";"+referencename+";"+yaxistitle+")").c_str());


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
      graph->Draw("ap");
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

      gdata->Draw("P");
      gdata->Print("all");
      
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
      text.DrawLatex(0.20,0.93,Form("%3.1f_{ }<_{ }|y|_{ }<_{ }%3.1f",input_table["ylow"][ibin],input_table["yhigh"][ibin]));
      text.DrawLatex(0.20,0.88,Form("%3.0f_{ }<_{ }p_{T}_{ }<_{ }%3.0f_{ }GeV",input_table["ptlow"][ibin],input_table["pthigh"][ibin]));

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
    gChi2->SetMinimum(0.75);
    gChi2->SetMaximum(1.2);
    gChi2->Draw("ap");
    if ( reference_values.size()+1<= 6 ) 
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
       TLegend legend(0.18,0.70,0.66,0.97,"","NDC");
       //legend.SetNColumns(3);
       legend.SetFillStyle(0);
       legend.SetBorderSize(0);
       legend.SetTextSize(0.045);
       legend.AddEntry(gChi2LTF,"#hat#chi^{2} of the linear template fit","P");
       legend.AddEntry(gChi2,   "#chi^{2}_{#font[12]{j}} of the individual templates","P");
       legend.AddEntry(gChi2->GetFunction("pol2"),"Parabola","L");
       legend.AddEntry(gChi2chk,"Minimum of #chi^{2} parabola #scale[0.8]{(#check#chi^{2})}","P"); //  (#check#chi^{2})
       legend.DrawClone();
    }
       
    c1.Print(ps_name);
    c1.Print( "plots/LTF_chi2.pdf");

    c1.Print( (string(ps_name)+"]").c_str() );
   
}


};
