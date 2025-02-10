//#include "include/LTF/LTF.h"
#include "LTF.h"

#include <TH1D.h>
#include <TGraphErrors.h>
#include <TGraph2DErrors.h>
#include <set>
#include <TCanvas.h>

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
std::map < std::string, std::vector<double> > read_input_table2(std::string filename, int ncol );



// __________________________________________________________________________________ //
//! 
//!  MakeTGraph
//!
//!  make a TGraph for plotting purposes
//!
static
TGraphErrors* MakeTGraph(const Eigen::VectorXd& xvalues, int ibin, const Eigen::MatrixXd& Y,
                         const std::map<std::string,Eigen::MatrixXd >& VSysY = {} );


// __________________________________________________________________________________ //
//! 
//!  MakeTGraph2D
//!
//!  make a TGraph2D for plotting purposes
//!
static
TGraph2DErrors* MakeTGraph2D(const Eigen::VectorXd& xvalues, const Eigen::VectorXd& yvalues,
                             int ibin, const Eigen::MatrixXd& Y,
                             const std::map<std::string,Eigen::MatrixXd >& VSysY = {} );


// __________________________________________________________________________________ //
//!
//! MakeHistogram
//!
//! make a histogram and fill it with random events according to a gauss
//! distribution around M
//!
static
TH1D* MakeHistogram(int nEvents, int seed, double mean, double sigma, vector<double> bins );



// __________________________________________________________________________________ //
//!
//!  MakeHistogram
//!
//!  make a histogram from an Eigen::Vector for plotting purposes
//!
static
TH1D* MakeHistogram(const Eigen::VectorXd& values, vector<double> bins ={}, const std::vector<std::pair<std::string,Eigen::MatrixXd > >& V = {} );


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
		 const string& xaxistitle    = "value [unit]",
		 const string& referencename = "Reference value (#alpha) [unit]");
//const string& observablename = "Observable [unit]");


// __________________________________________________________________________________ //
//!
//!
//!  Plot a LiTeFit object using ROOT
//!
//!  The binning needs to be provided to the plotting function,
//!  since this is not included in LTF::LiTeFit
//! 
static
void plotLiTeFit_2D(const LTF::LiTeFit& fit, const vector<double> bins );


// __________________________________________________________________________________ //
//!
//!
//!  Plot a LiTeFit object using ROOT
//!
//!  The binning needs to be provided to the plotting function,
//!  since this is not included in LTF::LiTeFit
//! 
static
void plotLiTeFitPol2Test(const LTF::LiTeFit& fit, const vector<double>& bins, 
                 const string& yaxistitle    = "value [unit]",
                 const string& referencename = "Reference value (#alpha) [unit]",
                 const string& observablename = "Observable [unit]");

static
   double makeErrorPlot(TCanvas&, const string&, const char*, const LTF::LiTeFit& fit, const vector<string>&);

};
