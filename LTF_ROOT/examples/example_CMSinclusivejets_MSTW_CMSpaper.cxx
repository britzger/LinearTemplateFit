/*
  Copyright 2021, D. Britzger, Max-Planck-Institute for Physics, Munich, Germany

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

#include <LTF/LTF_ROOTTools.h>
#include <LTF/LTF_Tools.h>
#include <LTF/LTF.h>


// -------------------------------------------------------------------- //
//!
//! example_CMSinclusivejets_MSTW_CMSpaper()
//!
//! Example application of the Linear Template Fit
//!
//! Determine the value of the strong coupling constant from
//! CMS inclusive jet data arXiv:1212.6660 using NLO QCD predictions,
//! similar to Eur.Phys.J.C 75 (2015) 288 [arXiv:1410.6765]
//!
//! Data published in Phys.Rev.D 87 (2013) 112002, PRD 87 (2013) 119902 (erratum) [arXiv:1212.6660]
//!
// -------------------------------------------------------------------- //

//! ------------------------------------------------------------------------ //
//! main function
#ifndef __CLING__

int example_CMSinclusivejets_MSTW_CMSpaper();

int main(int ,const char **) {

   gROOT->SetBatch();

   return example_CMSinclusivejets_MSTW_CMSpaper();
}
#endif

int example_CMSinclusivejets_MSTW_CMSpaper() {
   using namespace std;
#ifdef __CLING__
   TH1D::AddDirectory(false);
#endif

   // --- for 'CMS'-type fits  (similar to EPJ C 75 (2015) 288 [arXiv:1410.6765])
   map < string, vector<double> > templates = LTF_Tools::read_input_table("data/CMS_MSTWnlo_all.txt",16); string PDF="MSTW"; // CMS-type NLO prediction


   // --- data
   map < string, vector<double> > input_table = LTF_Tools::read_input_table("data/CMS_data.txt",32);
   std::vector<std::vector<double > > corr    = LTF_Tools::read_correlations("data/CMS_correlations.txt",133,true);
   vector<double> data                        = input_table["Sigma"];

   // apply np_corr and ewk_cor to NLO predictions
   vector<double> np_cor  = input_table["np_cor"];
   vector<double> ewk_cor = input_table["ewk_cor"];
   for ( auto& [name,v] : templates ) {
      for ( unsigned int i = 0 ; i<data.size(); i++ ) {
         v[i] *= np_cor[i] * ewk_cor[i];
      }
   }
   

   // --- linear template fit
   LTF ltf_CMS;
   ltf_CMS.SetData( data );

   // some settings
   ltf_CMS.SetGamma(vector<double>{1});
   ltf_CMS.UseNuisanceParameters(false);
   ltf_CMS.UseLogNormalUncertainties(false);

   
   // ---- set templates
   ltf_CMS.AddTemplate( 0.112 , templates[PDF+"nlo_0112"] );
   ltf_CMS.AddTemplate( 0.113 , templates[PDF+"nlo_0113"] );
   ltf_CMS.AddTemplate( 0.114 , templates[PDF+"nlo_0114"] );
   ltf_CMS.AddTemplate( 0.115 , templates[PDF+"nlo_0115"] );
   ltf_CMS.AddTemplate( 0.116 , templates[PDF+"nlo_0116"] );
   ltf_CMS.AddTemplate( 0.117 , templates[PDF+"nlo_0117"] );
   ltf_CMS.AddTemplate( 0.118 , templates[PDF+"nlo_0118"] );
   ltf_CMS.AddTemplate( 0.119 , templates[PDF+"nlo_0119"] );
   ltf_CMS.AddTemplate( 0.120 , templates[PDF+"nlo_0120"] );
   ltf_CMS.AddTemplate( 0.121 , templates[PDF+"nlo_0121"] );
   // unused templates
   // ltf_CMS.AddTemplate( 0.122 , templates[PDF+"nlo_0122"] );
   // ltf_CMS.AddTemplate( 0.123 , templates[PDF+"nlo_0123"] );
   // ltf_CMS.AddTemplate( 0.124 , templates[PDF+"nlo_0124"] );
   // ltf_CMS.AddTemplate( 0.125 , templates[PDF+"nlo_0125"] );
   // ltf_CMS.AddTemplate( 0.126 , templates[PDF+"nlo_0126"] );
   // ltf_CMS.AddTemplate( 0.127 , templates[PDF+"nlo_0127"] );

   

   // --- rescale uncertainties to avoid d'Agostini bias
   //     use NLO-prediction with 0.116
   string rescale = PDF+"nlo_0119";
   if( PDF=="MSTW" ) rescale = PDF+"nlo_0116";
   for ( unsigned int i = 0; i<data.size()+1; i++ ) {
      // additive uncertainties: do not rescale!
      //input_table["stat"       ][i]*= templates[rescale][i] / data[i];  // do not rescale this 'additive' uncertaintiy
      //input_table["uncor"      ][i]*= templates[rescale][i] / data[i];  // do not rescale this 'additive' uncertaintiy
      //input_table["npcorerr"   ][i]*= templates[rescale][i] / data[i];  // do not rescale this 'additive' uncertaintiy
      // multiplicative uncertainties
      // "The JES, unfolding, and luminosity uncertainties are treated as multiplicative" [paper]
      input_table["CMS_lumierr"][i]*= templates[rescale][i] / data[i];
      input_table["unf_err"    ][i]*= templates[rescale][i] / data[i];
      input_table["JEC0"       ][i]*= templates[rescale][i] / data[i];
      input_table["JEC1"       ][i]*= templates[rescale][i] / data[i];
      input_table["JEC3"       ][i]*= templates[rescale][i] / data[i];
      input_table["JEC4"       ][i]*= templates[rescale][i] / data[i];
      input_table["JEC5"       ][i]*= templates[rescale][i] / data[i];
      input_table["JEC6"       ][i]*= templates[rescale][i] / data[i];
      input_table["JEC7"       ][i]*= templates[rescale][i] / data[i];
      input_table["JEC8"       ][i]*= templates[rescale][i] / data[i];
      input_table["JEC9"       ][i]*= templates[rescale][i] / data[i];
      input_table["JEC10"      ][i]*= templates[rescale][i] / data[i];
      input_table["JEC11"      ][i]*= templates[rescale][i] / data[i];
      input_table["JEC12"      ][i]*= templates[rescale][i] / data[i];
      input_table["JEC13"      ][i]*= templates[rescale][i] / data[i];
      input_table["JEC14"      ][i]*= templates[rescale][i] / data[i];
      input_table["JEC15"      ][i]*= templates[rescale][i] / data[i];
      input_table["JEC2a"      ][i]*= templates[rescale][i] / data[i];
      input_table["JEC2b"      ][i]*= templates[rescale][i] / data[i];
      input_table["JEC2c"      ][i]*= templates[rescale][i] / data[i];
      input_table["JEC2d"      ][i]*= templates[rescale][i] / data[i];
      input_table["JEC2e"      ][i]*= templates[rescale][i] / data[i];
   }

   // --- set uncertainties
   auto cov  = LTF_Tools::corr_to_cov(corr, input_table["stat"], data);
   ltf_CMS.AddError( "stat",  cov  );
   ltf_CMS.AddErrorPercent( "uncor"       , input_table["uncor"      ], 0.  );
   ltf_CMS.AddErrorPercent( "npcorerr"    , input_table["npcorerr"   ], 1. ,  LTF::Uncertainty::External ); // not included in chi2, but uncertainty source is still propagated linearly
   ltf_CMS.AddErrorPercent( "lumi"        , input_table["CMS_lumierr"], 1.  );
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
   ltf_CMS.AddErrorPercent( "JEC7"        , input_table["JEC7"       ], 1.  );
   ltf_CMS.AddErrorPercent( "JEC10"       , input_table["JEC10"      ], 1.  );
   
   // --- PDF uncertainties
   map < string, vector<double> > pdferrors;
   if ( PDF == "MSTW" ) pdferrors = LTF_Tools::read_input_table("data/MSTW2008nlo68cl.errorFabsAbsolut.txt",20);
   for ( unsigned int i = 1; i<= pdferrors.size() ; i++ ) {
      //string ename = Form("PDF_%02d",i);
      char buffer [20];
      sprintf(buffer,"PDF_%02d",i);
      string ename = buffer;
      ltf_CMS.AddError( ename      , pdferrors[ename], 1.);
   }


   // ---- do the fit !
   LTF::LiTeFit fit = ltf_CMS.DoLiTeFit();
   fit.PrintFull();

   // --- Do iterative fit with 2nd order representation of the model
   fit.DoIterativeFitNewton() ;
   fit.PrintFull();
   

   // --- plot the results if ROOT is available
   vector<double> bins;
   for ( unsigned int i = 0; i<data.size()+1; i++ ) bins.push_back(i); // dummy binning for plotting-function
   LTF_ROOTTools::plotLiTeFit(fit,bins,"d^{2}#sigma_{jet}/dp_{T}d|y| [pb/GeV]","#alpha_{s}(m_{Z})","Bin");

   return 0; 
}
