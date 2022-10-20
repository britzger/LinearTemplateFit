
#include "LTF/LTF.h"

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>

class LTF_Tools {
public:

// _____________________________________________________________________________________ //
//! 
//!  read_correlations
//!
//!  Read a matrix from a simple ascii file
//!
static
   std::vector<std::vector<double > > read_correlations(std::string filename, int ncol, bool DoCorrTest=false);


// _____________________________________________________________________________________ //
//! 
//!  corr_to_cov
//!
//!  calculate a covariance matrix from a correlation matrix and uncertainties
//!
static
std::vector<std::vector<double > > corr_to_cov( const std::vector<std::vector<double > >& corr, const std::vector<double >& percenterr, const std::vector<double >& data);


// _____________________________________________________________________________________ //
//! 
//!  read_input_table
//!
//! read input data table from file and return it as a map
//!
static
std::map < std::string, std::vector<double> > read_input_table(std::string filename, int ncol );

};
