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

#include <string>
#include <set>
#include "LTF/LTF.h"
#include "LTF/LTF_Tools.h"
#if defined __WITH_ROOT__ || defined __CLING__
#include <TSystem.h>
#endif

using namespace std;

// __________________________________________________________________________________ //
//! main
int example2_LTF_gaus2D_NoROOT() {

#if defined __WITH_ROOT__ || defined __CLING__
   gSystem->Load("libLTF.so");
#endif

   // ------------------------------------------------ //
   // ---  Do linear template fit   
   LTF ltf;
   ltf.SetGamma(vector<double>{1});
   ltf.UseNuisanceParameters(false);
   ltf.UseLogNormalUncertainties(false);

   // --- binning of the histograms
   //const vector<double> bins{169,170,171,172,173, 174, 175, 176, 177, 178,179,180,181,182,183};

   // --- read and set templates
   const vector<double> reference_values1{169.5, 170, 170.5, 171}; // template reference points
   const vector<double> reference_values2{5.8, 6.0, 6.2, 6.4}; // template reference points
   const set<double>    central_values{170,170.5,6.0,6.2 }; // drop 'extreme' 2D points   

   map < string, vector<double> > templates   = LTF_Tools::read_input_table("data/example2_LTF_gaus2D.txt",26);
   for ( double mean : reference_values1 ) {
      for ( double sigm : reference_values2 ) {
         if ( central_values.find(mean) != central_values.end() || central_values.find(sigm) != central_values.end() ) {
            // set template
            char buffer [20];
            sprintf(buffer,"T_%.1f_%.1f",mean,sigm);
            const vector<double>& tmpl = templates[string(buffer)];
            ltf.AddTemplate( {mean,sigm},  tmpl ); // set template

            // set template uncertainty (optional)
            sprintf(buffer,"S_%.1f_%.1f",mean,sigm);
            const vector<double>& stat = templates[string(buffer)];
            ltf.AddTemplateError("statY", {mean,sigm} , stat, 0.); // set template error dY
         }
      }
   }
   

   // --- initialize data
   const vector<double>& data = templates["Data"];
   const vector<double>& stat = templates["Stat"];
   ltf.SetData( data);
   ltf.AddError("stat.", stat,0.);


   LTF::LiTeFit fit2 = ltf.DoLiTeFit();
   //LTF::LiTeFit fit2 = ltf.DoQuadraticTemplateFit();
   fit2.PrintFull();

   return 0;
  
}


//! ------------------------------------------------------------------------ //
//! main function
int main() {
   return example2_LTF_gaus2D_NoROOT();
}

