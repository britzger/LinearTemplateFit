
#ifdef __WITH_ROOT__
#include "plot_LTF1D.cxx"
#endif
//#include "LTF.cxx"
#include <iostream>
#include <fstream>
#include "LTF/LTF.h"

std::vector<std::vector<double > > read_correlations(const std::string& filename, int ncol) { 
   // read correlation files
   std::vector<std::vector<double > > ret(ncol);
   std::ifstream istrm(filename.c_str(), std::ios::binary);
   if (!istrm.is_open()) {
      std::cout << "failed to open " << filename << std::endl;
      exit(1);
   }
   for ( int icol=0 ; icol<ncol ; icol++ ) {
      double value;
      for ( int i=0 ; i<ncol ; i++ ) {
         istrm >> value;
         ret[icol].push_back(value);
         // if ( icol==i && value!= 1 ) {
         //    std::cout<<"ERROR! diagonal element must be 1, but is: "<<value<<endl;
         //    exit(1);
         // }
      }
      //std::cout<<"read: "<< ret[icol].size()<<" values in row "<<icol<<std::endl;
   }
   return ret;
}

std::vector<std::vector<double > > corr_to_cov( const std::vector<std::vector<double > >& corr, const std::vector<double >& percenterr, const std::vector<double >& data) { 
   std::vector<double > err(data.size());
   for ( int i = 0 ; i<err.size() ; i++ ) 
      err[i] = percenterr[i]/100.*data[i];
   
   std::vector<std::vector<double > > cov(err.size());
   for ( int i = 0 ; i<err.size() ; i++ ) {
      cov[i].resize(err.size());
      for ( int j = 0 ; j<err.size() ; j++ ) {
         cov[i][j] = corr[i][j]*err[i]*err[j];
         //if ( corr[i][j] != 0 ) cout<<i<<", "<<j<<"\t"<<corr[i][j]<<"\t"<<percenterr[i]<<"\t"<<data[i]<<"\te: "<<err[i]<<"\tcov: "<<cov[i][j]<<endl;
      }
   }

   for ( int i = 0 ; i<err.size() ; i++ ) {
      for ( int j = 0 ; j<i ; j++ ) {
         if ( fabs((cov[i][j] / cov[j][i]-1.) > 1e-3) )
            cout<<i<<", "<<j<<"\t"<<corr[i][j]<<"\t"<<percenterr[i]<<"\t"<<data[i]<<"\te: "<<err[i]<<"\tcov: "<<cov[i][j]<<"\t ratio ij/ji: "<<(cov[i][j] / cov[j][i]-1.)<<endl;
      }
   }

   return cov;
}



std::map < std::string, std::vector<double> > read_input_table(std::string filename, int ncol ) {
   // read input data table from file
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
      //cout<<colname<<endl;
      ret[colname] = vector<double>();
      cols.push_back(colname);
   }
   while ( istrm.good()) {
      double value;
      for ( int c=0 ; c<ncol ; c++ ) {
         istrm >> value;
         if ( !istrm.good() ) break;
         //cout<<value<<"  ";
         std::string colname = cols[c];
         ret[colname].push_back(value);
      }      
      //cout<<endl;
   }
   cout<<"Info.  [read_input_table]  Read "<<ret.size()<<" rows."<<endl;
   return ret;
}

// ----------------------------------------------------------------------------------------------------
//
// Example: determine the value of the strong coupling constant from 
// CMS inclusive jet data arXiv:1212.6660 using NLO QCD predictions.
///
// ----------------------------------------------------------------------------------------------------
int example_CMSinclusivejets_NN30_BRSSWpaper() {

#ifdef __CLING__
   gSystem->Load("libLTF.so");
   TH1D::AddDirectory(false);
#endif
   //gSystem->Load("LTF_cxx.so");
   using namespace std;
   cout<<"Start."<<endl;

   // --- for 'CMS'-type fits
   //map < string, vector<double> > templates_CT10nlo = read_input_table("CMS_CT10nlo_all.txt",16); string PDF="CT10"; 
   //map < string, vector<double> > templates_CT10nlo = read_input_table("CMS_MSTWnlo_all.txt",16); string PDF="MSTW"; // almost ok, if rescaling, and sign of PDF uncertainty is droped

   // --- for 'common'-type and 'H1'-type fits 
   //map < string, vector<double> > templates_CT10nlo = read_input_table("CT10nlo_fixedPDF.txt",16); string PDF="CT10";
   //map < string, vector<double> > templates_CT10nlo = read_input_table("MSTWnlo_fixedPDF.txt",16);  string PDF="MSTW";
   //map < string, vector<double> > templates_CT10nlo = read_input_table("CT14nlo_fixedPDF.txt",16);  string PDF="CT14"; // ok, wie bei BRSSW
   map < string, vector<double> > templates_CT10nlo = read_input_table("data/NN30nlo_fixedPDF.txt",16);  string PDF="NN30"; // ok, wie bei BRSSW. Central result

   // --- data
   map < string, vector<double> > input_table = read_input_table("data/CMS_data.txt",32);
   std::vector<std::vector<double > > corr =read_correlations("data/CMS_correlations.txt",133);
   vector<double> data    = input_table["Sigma"];

   // apply np_corr and ewk_cor to NLO predictions
   vector<double> np_cor  = input_table["np_cor"];
   vector<double> ewk_cor = input_table["ewk_cor"];
   for ( auto& [name,v] : templates_CT10nlo ) {
      for ( int i = 0 ; i<data.size(); i++ ) {
         v[i] *= np_cor[i] * ewk_cor[i];
      }
   }
   

   // --- linear template fit
   LTF ltf_CMS;
   ltf_CMS.SetData( data );

   // some settings
   ltf_CMS.SetGamma(vector<double>{1});
   ltf_CMS.UseNuisanceParameters(false);
   ltf_CMS.UseLogNormalUncertainties(true);

   
   // ---- set templates
   ltf_CMS.AddTemplate( 0.112 , templates_CT10nlo[PDF+"nlo_0112"] );
   ltf_CMS.AddTemplate( 0.113 , templates_CT10nlo[PDF+"nlo_0113"] );
   ltf_CMS.AddTemplate( 0.114 , templates_CT10nlo[PDF+"nlo_0114"] );
   ltf_CMS.AddTemplate( 0.115 , templates_CT10nlo[PDF+"nlo_0115"] );
   ltf_CMS.AddTemplate( 0.116 , templates_CT10nlo[PDF+"nlo_0116"] );
   ltf_CMS.AddTemplate( 0.117 , templates_CT10nlo[PDF+"nlo_0117"] );
   ltf_CMS.AddTemplate( 0.118 , templates_CT10nlo[PDF+"nlo_0118"] );
   ltf_CMS.AddTemplate( 0.119 , templates_CT10nlo[PDF+"nlo_0119"] );
   ltf_CMS.AddTemplate( 0.120 , templates_CT10nlo[PDF+"nlo_0120"] );
   ltf_CMS.AddTemplate( 0.121 , templates_CT10nlo[PDF+"nlo_0121"] );

   // ltf_CMS.AddTemplate( 0.122 , templates_CT10nlo[PDF+"nlo_0122"] );
   // ltf_CMS.AddTemplate( 0.123 , templates_CT10nlo[PDF+"nlo_0123"] );
   // ltf_CMS.AddTemplate( 0.124 , templates_CT10nlo[PDF+"nlo_0124"] );
   // ltf_CMS.AddTemplate( 0.125 , templates_CT10nlo[PDF+"nlo_0125"] );
   // ltf_CMS.AddTemplate( 0.126 , templates_CT10nlo[PDF+"nlo_0126"] );
   // ltf_CMS.AddTemplate( 0.127 , templates_CT10nlo[PDF+"nlo_0127"] );

   

   // ---- set uncertainties
   //ltf_CMS.AddErrorPercent( "stat",  input_table["stat"],  0.  );
   auto cov  = corr_to_cov(corr, input_table["stat"], data);
   ltf_CMS.AddError( "stat",  cov  );
   ltf_CMS.AddErrorPercent( "uncor"       , input_table["uncor"      ], 0.  );
   ltf_CMS.AddErrorPercent( "npcorerr"    , input_table["npcorerr"   ], 1.  ); // CMS not included in chi2 (LTF::Uncertainty::External), BRSSW included
   ltf_CMS.AddErrorPercent( "lumi"        , input_table["CMS_lumierr"], 1.  ); //, LTF::Uncertainty::External );
   ltf_CMS.AddErrorPercent( "unfolding"   , input_table["unf_err"    ], 1.  );
   ltf_CMS.AddErrorPercent( "JEC0"        , input_table["JEC0"       ], 1.  );
   ltf_CMS.AddErrorPercent( "JEC1"        , input_table["JEC1"       ], 1.  );
   ltf_CMS.AddErrorPercent( "JEC3"        , input_table["JEC3"       ], 1.  );
   ltf_CMS.AddErrorPercent( "JEC4"        , input_table["JEC4"       ], 1.  );
   ltf_CMS.AddErrorPercent( "JEC5"        , input_table["JEC5"       ], 1.  );
   ltf_CMS.AddErrorPercent( "JEC6"        , input_table["JEC6"       ], 1.  );
   ltf_CMS.AddErrorPercent( "JEC9"        , input_table["JEC9"       ], 1.  );
   ltf_CMS.AddErrorPercent( "JEC8"        , input_table["JEC8"       ], 1.  );
   ltf_CMS.AddErrorPercent( "JEC11"       , input_table["JEC11"      ], 1.  );
   ltf_CMS.AddErrorPercent( "JEC12"       , input_table["JEC12"      ], 1.  );
   ltf_CMS.AddErrorPercent( "JEC13"       , input_table["JEC13"      ], 1.  );
   ltf_CMS.AddErrorPercent( "JEC14"       , input_table["JEC14"      ], 1.  );
   ltf_CMS.AddErrorPercent( "JEC15"       , input_table["JEC15"      ], 1.  );
   ltf_CMS.AddErrorPercent( "JEC2a"       , input_table["JEC2a"      ], 1.  );
   ltf_CMS.AddErrorPercent( "JEC2b"       , input_table["JEC2b"      ], 1.  );
   ltf_CMS.AddErrorPercent( "JEC2c"       , input_table["JEC2c"      ], 1.  );
   ltf_CMS.AddErrorPercent( "JEC2d"       , input_table["JEC2d"      ], 1.  );
   ltf_CMS.AddErrorPercent( "JEC2e"       , input_table["JEC2e"      ], 1.  );
   ltf_CMS.AddErrorPercent( "JEC7"        , input_table["JEC7"       ], 1.  ); // zero
   ltf_CMS.AddErrorPercent( "JEC10"       , input_table["JEC10"      ], 1.  ); // zero
   
   if ( PDF=="NN30" || PDF=="NN31") {
      auto covPDF_rel = read_correlations("data/NNPDF30_nlo_as_0118.PDFuncertainties.rel.txt",133);
      ltf_CMS.AddErrorRelative("NNPDF", covPDF_rel );
   }
   else {
      map < string, vector<double> > pdferrors;
      if ( PDF == "CT10" ) pdferrors = read_input_table("data/log.CT10.errorAbsolut.txt",26);
      if ( PDF == "MSTW" ) pdferrors = read_input_table("data/log.MSTW2008nlo68cl.errorFabsAbsolut.txt",20);
      //if ( PDF == "MSTW" ) pdferrors = read_input_table("log.MSTW2008nlo68cl.errorAbsolut.txt",20);
      //map < string, vector<double> > pdferrors = read_input_table("log.CT10.errorPercent.txt",26);
      //map < string, vector<double> > pdferrors = read_input_table("log.CT10.errorFabsAbsolut.txt",26);
      //map < string, vector<double> > pdferrors = read_input_table("log.CT10.errorAbsolut.txt",26);
      //map < string, vector<double> > pdferrors = read_input_table("log.MMHT2014nlo68clas118.errorAbsolut.txt",25);
      //map < string, vector<double> > pdferrors = read_input_table("log.MSTW2008nlo68cl.errorPercent.txt",20);
      //map < string, vector<double> > pdferrors = read_input_table("log.MSTW2008nlo68cl.errorAbsolut.txt",20);
      //map < string, vector<double> > pdferrors = read_input_table("log.CT14.errorPercent.txt",28);
      for ( int i = 1; i<= pdferrors.size() ; i++ ) {
         //string ename = Form("PDF_%02d",i);
         char buffer [20];
         sprintf(buffer,"PDF_%02d",i);
         string ename = buffer;
         // //for ( int b = 0 ; b<pdferrors[ename].size() ;b++ ) pdferrors[ename][b] *=0.5;
         if ( PDF=="CT10" || PDF=="CT14") { // 68%CL
            for ( int b = 0 ; b<pdferrors[ename].size() ;b++ ) pdferrors[ename][b] /= 1.645;
            //for ( int b = 0 ; b<pdferrors[ename].size() ;b++ ) pdferrors[ename][b] /= templates_CT10nlo[PDF+"nlo_0118"][b] / data[b];

         }
         // fabs??
         //for ( int b = 0 ; b<pdferrors[ename].size() ;b++ ) pdferrors[ename][b] = fabs(pdferrors[ename][b]);

         //ltf_CMS.AddErrorPercent( ename      , pdferrors[ename], 1.);
     
         ltf_CMS.AddError( ename      , pdferrors[ename], 1.);

         // ltf_CMS.AddErrorPercent( ename      , pdferrors[ename], 1.);
      }
   }

   // ---- do the fit !
   LTF::LiTeFit fit = ltf_CMS.DoLiTeFit();
   fit.PrintFull();
   fit.DoIterativeFitNewton() ;
   if ( !fit.GetLogNormal() ) { // repeat fit with 'expected' uncertainties
      //fit.RescaleInputErrors(fit.CalculatePrediction(0.117));
      //fit.DoIterativeFitNewton() ;
      //fit.DoLiTeFit();
   }
   fit.PrintFull();
   
   //return 0;

   //fit.RescaleInputErrors(templates_CT10nlo[PDF+"nlo_0118"]);
   //fit.RescaleInputErrors(fit.CalculatePrediction(0.1185));
   //fit.DoLiTeFit();
   //fit.PrintFull();

   // fit.RescaleInputErrors(fit.CalculatePrediction(0.116));
   // fit.DoLiTeFit();
   // fit.PrintFull();
   
   // fit.DoIterativeFitNewton() ;
   // fit.PrintFull();


   vector<double> bins;
   for ( int i = 0; i<data.size()+1; i++ ) bins.push_back(i);
   cout<<"data.size: "<<data.size()<<"\t"<<bins.size()<<endl;

#ifdef __WITH_ROOT__
   plotLiTeFit(fit,bins,"d^{2}#sigma_{jet}/dp_{T}d|y| [pb/GeV]","Reference value:  #alpha_{s}(M_{Z})","Bin");//, const vector<double> bins)
#endif
   return 0;


   // fit.DoIterativeFitNewton() ;
   // fit.PrintFull();
   
   // for ( int b = 0 ; b<data.size() ;b++ ) {
   //    cout<<"Ratio theo/data, bin "<<b<<"\t"<<templates_CT10nlo[PDF+"nlo_0117"][b]/data[b]<<"\tCalcPred: "<< fit.CalculatePrediction(0.117)(b)/data[b] <<endl;
   // }
   //return 0;

   // iteratively rescale errors to avoid bias.
   if ( !fit.GetLogNormal() ) { 
      for ( int i=0; i<3; i++ ) {
         for ( int b = 0 ; b<data.size() ;b++ ) {
            cout<<"Ratio theo/data, bin "<<b<<"\t"<<fit.TheoFit(b)/data[b]<<endl;
         }

         fit.RescaleInputErrors(fit.TheoFit);
         fit.DoLiTeFit();
      }
      fit.PrintFull();
   }

   return 0;
}


int main() {
   return example_CMSinclusivejets_NN30_BRSSWpaper();
}
