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

#include "LTF/LTF.h"
#include "LTF_Tools.h"
#if defined __WITH_ROOT__ || defined __CLING__
#include <TSystem.h>
#endif

// __________________________________________________________________________________ //
//! 
//! example_LTF_gaus_NoROOT
//!
//! Example for the linear template fit
//! The templates are generated as a gauss-distribution and stored in the file 
//! "data/example_LTF_gaus.txt"
//! The templates are generated with ROOT::TRandom3. See example_LTF_gaus.cxx
//! for more details.
//!
//! More details are given in the paper.
//!
int example1_LTF_gaus_NoROOT() {
   using namespace std;
#if defined __WITH_ROOT__ || defined __CLING__
   gSystem->Load("libLTF.so");
#endif

   // ------------------------------------------------ //
   // --- read templates and data
   const vector<double> reference_values{169, 169.5, 170, 170.5, 171, 171.5, 172}; // template reference points
   map < string, vector<double> > templates   = LTF_Tools::read_input_table("data/example_LTF_gaus.txt",16);

   // ------------------------------------------------ //
   // ---  Do linear template fit
   LTF ltf;
   ltf.SetGamma(vector<double>{1});
   ltf.UseNuisanceParameters(false);
   ltf.UseLogNormalUncertainties(false);

   // --- read and set templates
   for ( double ref : reference_values ) {
      char buffer [20];
      sprintf(buffer,"Tmpl_%.2f",ref);
      const vector<double>& tmpl = templates[string(buffer)];
      ltf.AddTemplate(ref,  tmpl ); // set template
      sprintf(buffer,"Stat_%.2f",ref);
      const vector<double>& stat = templates[string(buffer)];
      ltf.AddTemplateError("statY", ref , stat, 0.); // set template error dY
   }
   
   // --- set data
   const vector<double>& data = templates["Data"];
   const vector<double>& stat = templates["Stat"];
   ltf.SetData( data);
   ltf.AddError("stat.", stat,0.);
   

   // --- perform the linear template fit and print results
   LTF::LiTeFit fit = ltf.DoLiTeFit();
   fit.PrintFull();

   return 0;
}

//! ------------------------------------------------------------------------ //
//! main function
int main() {
   return example1_LTF_gaus_NoROOT();
}

