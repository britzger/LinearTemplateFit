//#include "include/LTF/LTF.h"
#include "LTF.h"

#include <TH1D.h>
#include <TGraphErrors.h>
#include <TGraph2DErrors.h>
#include <set>

using namespace std;

namespace LTF_ROOTTools {

// __________________________________________________________________________________ //
//!
//! read_input_table2()
//!
//! read input data table from file 
//!
std::map < std::string, std::vector<double> > read_input_table2(std::string filename, int ncol );



// __________________________________________________________________________________ //
//! 
//!  MakeTGraph
//!
//!  make a TGraph for plotting purposes
//!
TGraphErrors* MakeTGraph(const TVectorD& xvalues, int ibin, const TMatrixD& Y,
                         const std::map<std::string,TMatrixD >& VSysY = {} );


// __________________________________________________________________________________ //
//! 
//!  MakeTGraph2D
//!
//!  make a TGraph2D for plotting purposes
//!
TGraph2DErrors* MakeTGraph2D(const TVectorD& xvalues, const TVectorD& yvalues,
                             int ibin, const TMatrixD& Y,
                             const std::map<std::string,TMatrixD >& VSysY = {} );


// __________________________________________________________________________________ //
//!
//! MakeHistogram
//!
//! make a histogram and fill it with random events according to a gauss
//! distribution around M
//!
TH1D* MakeHistogram(int nEvents, int seed, double mean, double sigma, vector<double> bins );



// __________________________________________________________________________________ //
//!
//!  MakeHistogram
//!
//!  make a histogram from an Eigen::Vector for plotting purposes
//!
TH1D* MakeHistogram(const TVectorD& values, vector<double> bins ={}, const std::vector<std::pair<std::string,TMatrixD > >& V = {} );


// __________________________________________________________________________________ //
//!
//!
//!  Plot a LiTeFit object using ROOT
//!
//!  The binning needs to be provided to the plotting function,
//!  since this is not included in LTF::LiTeFit
//! 
void plotLiTeFit(const LTF::LiTeFit& fit, const vector<double>& bins, 
                 const string& yaxistitle    = "value [unit]",
                 const string& referencename = "Reference value (#alpha) [unit]",
                 const string& observablename = "Observable [unit]");


// __________________________________________________________________________________ //
//!
//!
//!  Plot a LiTeFit object using ROOT
//!
//!  The binning needs to be provided to the plotting function,
//!  since this is not included in LTF::LiTeFit
//! 
void plotLiTeFit_2D(const LTF::LiTeFit& fit, const vector<double> bins );


// __________________________________________________________________________________ //
//!
//!
//!  Plot a LiTeFit object using ROOT
//!
//!  The binning needs to be provided to the plotting function,
//!  since this is not included in LTF::LiTeFit
//! 
void plotLiTeFitPol2Test(const LTF::LiTeFit& fit, const vector<double>& bins, 
                 const string& yaxistitle    = "value [unit]",
                 const string& referencename = "Reference value (#alpha) [unit]",
                 const string& observablename = "Observable [unit]");


}
