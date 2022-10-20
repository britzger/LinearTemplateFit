#include "LTF/LTF.h"
#include "LTF/LTF_Tools.h"

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>


// _____________________________________________________________________________________ //
//! 
//!  read_correlations
//!
//!  Read a matrix from a simple ascii file
//!
std::vector<std::vector<double > > LTF_Tools::read_correlations(std::string filename, int ncol, bool DoCorrTest) { 
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
         if ( DoCorrTest && icol==i && value!= 1 ) {
            std::cout<<"ERROR! diagonal element must be 1, but is: "<<value<<endl;
            exit(1);
         }
      }
      //std::cout<<"read: "<< ret[icol].size()<<" values in row "<<icol<<std::endl;
   }
   return ret;
}



// _____________________________________________________________________________________ //
//! 
//!  corr_to_cov
//!
//!  calculate a covariance matrix from a correlation matrix and uncertainties
//!
std::vector<std::vector<double > > LTF_Tools::corr_to_cov( const std::vector<std::vector<double > >& corr, const std::vector<double >& percenterr, const std::vector<double >& data) { 
   std::vector<double > err(data.size());
   for ( size_t i = 0 ; i<err.size() ; i++ ) 
      err[i] = percenterr[i]/100.*data[i];
   
   std::vector<std::vector<double > > cov(err.size());
   for ( size_t i = 0 ; i<err.size() ; i++ ) {
      cov[i].resize(err.size());
      for ( size_t j = 0 ; j<err.size() ; j++ ) {
         cov[i][j] = corr[i][j]*err[i]*err[j];
         //if ( corr[i][j] != 0 ) cout<<i<<", "<<j<<"\t"<<corr[i][j]<<"\t"<<percenterr[i]<<"\t"<<data[i]<<"\te: "<<err[i]<<"\tcov: "<<cov[i][j]<<endl;
      }
   }

   for ( size_t i = 0 ; i<err.size() ; i++ ) {
      for ( size_t j = 0 ; j<i ; j++ ) {
         if ( fabs((cov[i][j] / cov[j][i]-1.) > 1e-3) )
            cout<<i<<", "<<j<<"\t"<<corr[i][j]<<"\t"<<percenterr[i]<<"\t"<<data[i]<<"\te: "<<err[i]<<"\tcov: "<<cov[i][j]<<"\t ratio ij/ji: "<<(cov[i][j] / cov[j][i]-1.)<<endl;
      }
   }

   return cov;
}


// _____________________________________________________________________________________ //
//! 
//!  read_input_table
//!
//! read input data table from file and return it as a map
//!
std::map < std::string, std::vector<double> > LTF_Tools::read_input_table(std::string filename, int ncol ) {
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

