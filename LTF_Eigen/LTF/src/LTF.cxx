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
/*
 * LTF
 * The Linear template fit, implemented using the Eigen library.
 *
 */
// -------------------------------------------------------------------- //

#include "LTF/LTF.h"
#include "Eigen/Sparse"
#include "Eigen/Dense"
#include <string>
#include <map>
#include <unordered_map>
#include <iostream>
#include <memory>

#include "Eigen/src/Core/util/StaticAssert.h"

using namespace std;

// -------------------------------------------------------------------- //
//                     Eigen template instantiations
// -------------------------------------------------------------------- //
//template class Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
//template class Eigen::Matrix< long double, Eigen::Dynamic, 1              >;

// -------------------------------------------------------------------- //
//                     Implementations
// -------------------------------------------------------------------- //
void LTF::AddError(std::string name, size_t n, const double* error, double corr, LTF::Uncertainty constr) {
   // add a new error source to LTF
   if ( name=="" ) name = "Sys_" + to_string(fSys.size()+fVs.size()+fSysExt.size()+fVsExt.size());
   if ( corr<0 || corr>1 ) { cout<<"Error! Correlation coefficient must be between [0,1]."<<endl; exit(1);}
   if ( corr!= 1 && constr==LTF::Uncertainty::Unconstrained ) { cout<<"Warning! only fully correlated errors can be treated as 'unconstrained'."<<endl;}
   if ( n != size_t(fDt.rows()) ) { 
      std::cerr<<"Error in AddError()."
               <<" Number of values ("<<n<<") not consistent with number of data points ("<<fDt.rows()<<")."<<endl;
      exit(1);
   }
   if (corr<1) {
      auto uname =name;
      if ( corr!= 0 ) uname += "_u";
      std::vector<double> vect(error, error + n); // copy to non-const
      for ( double& el : vect ) el *= el * ( 1-corr );
      AddUncorrelatedErrorSquared(uname,n,&vect[0],constr);
   }
   if ( corr>0 ) {
      //fSys[name] = Eigen::Map<Eigen::VectorXd>(&vect[0], n);
      if ( corr!= 1 ) name += "_c";
      std::vector<double> vect(error, error + n); // copy to non-const
      for ( double& el : vect ) el *= sqrt(corr);
      AddCorrelatedError(name,vect,constr);
   }
}


// -------------------------------------------------------------------- //
void LTF::AddError(const std::string& name, const std::vector<std::vector<double> >& cov, LTF::Uncertainty constr) {
   // add a new error source to LTF
   if ( cov.size() != size_t(fDt.rows()) ) { 
      std::cerr<<"Error in AddError()."
               <<" Number of error values ("<<cov.size()<<") not consistent with number of data points ("<<fDt.rows()<<")."<<endl;
      exit(1);
   }
   vector<double> vect;
   for ( auto v1 : cov ) 
      for ( auto v2 : v1 ) 
         vect.push_back(v2);
   if ( constr == LTF::Uncertainty::External )
      fVsExt[name] = Eigen::Map<Eigen::MatrixXd>(&vect[0], cov.size(),cov.size());
   else {
      fVs.push_back(make_pair(name,Eigen::Map<Eigen::MatrixXd>(&vect[0], cov.size(),cov.size())));
      if ( constr == LTF::Uncertainty::Unconstrained )
         cout<<"Warning! Uncorrelated error sources or covariance matrices cannot be 'unconstrained'. using 'constrained' instead.."<<endl; 
   }
}


// -------------------------------------------------------------------- //
void LTF::AddErrorRelative(const std::string& name, const std::vector<std::vector<double> >& cov_rel, LTF::Uncertainty constr) {
   // add a new error source to LTF
   if ( cov_rel.size() != size_t(fDt.rows()) ) { 
      std::cerr<<"Error in AddError()."
               <<" Number of error values ("<<cov_rel.size()<<") not consistent with number of data points ("<<fDt.rows()<<")."<<endl;
      exit(1);
   }
   std::vector<std::vector<double> > cov = cov_rel;
   for ( size_t i = 0 ; i<cov.size() ;i++ ){
      for ( size_t j = 0 ; j<cov.size() ; j++ ){
         cov[i][j] *= fDt(i)*fDt(j);
      }
   }
   AddError(name,cov,constr);
}


// -------------------------------------------------------------------- //
void LTF::SetCorrSys(const std::string& name, double c) {
   //set correlation coefficient for syst. uncertainties
   if ( c!=1 && c!= 0 ) {
      cerr<<"Error in LTF::SetCorrSys! Correlation coefficient must be 0 or 1."<<endl; 
      exit(1);
   }
   if ( fcorrSys.count(name) && fcorrSys[name] != c ) {
      cout<<"Warning in LTF::SetCorrSys! Correlation coefficient for uncertainty '"
          <<name<<"' already initialized with a different value ("<<fcorrSys[name]<<")."<<endl;
   }
   fcorrSys[name] = c;
}

// -------------------------------------------------------------------- //
void LTF::AddTemplateError(const std::string& name, const vector<double>& mvals, size_t n, const double* error, double corr) {
   // specify an uncertainty of the templates

   // 1) find the reference values in 'Y'
   int k = -1;
   for ( auto [mvs,yy] : fTmplDistN ) {
      k++;
      if ( mvs == mvals ) break;
   }
   if ( k==-1 ) {
      cerr<<"Error in LTF::AddTemplateError! No template for the given reference values were found. Ignoring call!"<<endl; 
      return;
   }
   // 2) add the uncertainty to the member variables
   SetCorrSys(name,corr);
   std::vector<double> vect(error, error + n); // copy to non-const
   fSysY[name][mvals] = Eigen::Map<Eigen::VectorXd>(&vect[0],n);
}




// -------------------------------------------------------------------- //
void LTF::AddTemplate(const vector<double>& mvals, const vector<double>& dist) {
   // check if template is already available
   for ( auto r : fTmplDistN ) {
      if ( mvals.size() != r.first.size() ) { 
         std::cerr<<"Error. Different number of reference points given. Exiting."<<std::endl; 
         exit(1);
      }
      bool allsame = true;
      for ( size_t i = 0 ; i<mvals.size() ; i++ )  allsame &= (mvals[i]==r.first[i]);
      if ( allsame ) 
         std::cout<<"Warning! Template reference value[0] (mval="<<mvals[0]<<") already specified. Replacing previous template."<<std::endl;
   }
   std::vector<double> tmpl = dist; // copy to non-const for eigen...
   //fTmplDistN[mvals] = Eigen::Map<Eigen::VectorXd>(&tmpl[0], dist.size());
   fTmplDistN.push_back(std::make_pair(mvals,Eigen::Map<Eigen::VectorXd>(&tmpl[0], dist.size())));
}


// -------------------------------------------------------------------- //
void LTF::SetResponseMatrix(const std::vector<std::vector<double> >& A) {
   // set a new response matrix
   // If vector A is empty, then delete the previous setting.
   if ( A.empty() ) {
      std::cout<<"Info.    Response matrix has zero bins. Deleting response matrix."<<std::endl;
      fA = Eigen::MatrixXd();
   }
   else {
      if ( fDt.rows() == 0 || fTmplDistN.empty()) {
         std::cerr<<"Error in SetResponseMatrix. Data or templates not yet available. "
                  <<"Please first set data and template distributions."<<endl;
         exit(1);
      }
      if ( A.size() != size_t(fDt.rows()) ) {
         std::cerr<<"Error. size of response matrix does not fit size of data"
                  <<". A.size: "<<A.size()
                  <<", Data.size: "<<fDt.rows()
                  <<std::endl;
         exit(1);
      }
      if ( A[0].size() != size_t(fTmplDistN.begin()->second.size()) ) {
         std::cerr<<"Error. size of response matrix does not fit size of templates."<<std::endl;
         exit(1);
      }
      fA = LTF::Std_to_Mat(A);
   }
}


// -------------------------------------------------------------------- //
void LTF::SetLiTeFitInput() {
   //!< \brief Apply fit range, according to input parameters specified by SetFitRange()
   //!< Initialize 'input' parameters of output object fLTF.
   // --- input paramters of fLTF 
  if ( fGamma.size() ) fGamma.resize(fTmplDistN.begin()->first.size(),1);// resize, if Gamma was set before number of fit parameters was known.
   fLTF.Gamma       = fGamma;
   fLTF.Dt          = fDt;              
   fLTF.M           = Calc_M();
   fLTF.Y           = Calc_Y();
   fLTF.X           = Calc_X();
   fLTF.A           = fA;
   fLTF.Vs          = fVs;  // Covariance matrices
   fLTF.Sys         = fSys; // systematic uncertainty
   fLTF.VsExt       = fVsExt;  // Covariance matrices, treated as external
   fLTF.SysExt      = fSysExt; // systematic uncertainty
   fLTF.SysIsConstr = fSysIsConstr;//
   fLTF.SysY        = Calc_dY(fSysY); // uncertainty of the templates,  dY
   fLTF.SysA        = fSysA; // uncertainty of response matrix A, dA
   fLTF.corrSys     = fcorrSys; // correlation coefficients of dY, dA
   fLTF.y0new       = fy0New; // new nominal template
   fLTF.a0new       = fa0New; // reference value of new nominal template

   // request log-normal uncertainties
   if ( fUseLog ) 
      fLTF.ApplyLogNormalDistribution();
   
   // rescale errors, if requested
   if ( fErrorScale.rows() )
      fLTF.RescaleInputErrors(fErrorScale);

   // --- adapt fit range
   if ( fBinNN==0 ) { // --- sanity checks
      fBinNN=N();
      return;  
   }

   if ( fBinMin+fBinNN > fLTF.Dt.rows() ) {
      cerr<<"Error! More bins are requested by SetFitRange than available. BinMin="<<fBinMin
          <<", nBins="<<fBinNN<<", nBinsData="<<fLTF.Dt.rows()<<endl;
      exit(1);
   }

   // no changes to M, Mc
   std::cout<<"Info.    Apply fit range. First bin "<<fBinMin<<", last bin: "<<fBinMin+fBinNN-1<<std::endl;
   fLTF.Dt    = fLTF.Dt.segment(fBinMin,fBinNN).eval(); 
   fLTF.Y     = fLTF.Y.block(fBinMin,0,fBinNN,fLTF.M.rows()).eval(); //(don't forget to call 'eval()'!
   for ( auto& [name,V] : fLTF.Vs )    V = V.block(fBinMin,fBinMin,fBinNN,fBinNN).eval();
   for ( auto& [name,s] : fLTF.Sys )   s = s.segment(fBinMin,fBinNN).eval();
   for ( auto& [name,V] : fLTF.VsExt ) V = V.block(fBinMin,fBinMin,fBinNN,fBinNN).eval();
   for ( auto& [name,s] : fLTF.SysExt) s = s.segment(fBinMin,fBinNN).eval();
   for ( auto& [name,V] : fLTF.SysY )  V = V.block(fBinMin,0,fBinNN,V.cols()).eval();
   if ( fLTF.SysA.size() ) { // warning. not implemented
      cerr<<"Selection of fit range not implemented for uncertainties on migration matrix."<<endl;
      cout<<"Please provide this (trivial) feature."<<endl;
      exit(1);
   }
   //for ( auto& [name,V] : fLTF.SysA )   V = V.block(fBinMin,fBinMin,fBinNN,fBinNN).eval();
   
}



// -------------------------------------------------------------------- //
size_t LTF::N() const {
   // return number of points. 
   // This value is not valid after applying apply Fit Range
   size_t n = fDt.size(); 
   if ( n==0 ) {
      if ( fTmplDistN.size() ) {
         n = fTmplDistN.begin()->second.size();
      }
   }
   if ( n==0 ) 
      cout<<"Warning! Number of points is zero."<<endl;
   return n;
}


// -------------------------------------------------------------------- //
Eigen::MatrixXd LTF::Calc_M() const {
   // Calculate matrix M from template reference points
   Eigen::MatrixXd M(fTmplDistN.size(),fTmplDistN.begin()->first.size()+1);
   for ( size_t i = 0 ; i<size_t(M.rows()) ; i++ )  M(i,0) = 1;
   int i = 0;
   for ( auto p : fTmplDistN ) {
      int j  = 1;
      for ( auto a : p.first ) M(i,j++) = a;
      i++;
   }
   return M;
}

// -------------------------------------------------------------------- //
Eigen::MatrixXd LTF::Calc_Y() const {
   // Calculate matrix Y, representing the templates
   Eigen::MatrixXd Y;
   if ( fA.rows() > 0  && fTmplDistN.size()> 0 )  
      Y = Eigen::MatrixXd::Zero(fA.rows(),fTmplDistN.size());
   else if ( N()>0 && fTmplDistN.size()>0 )
      Y = Eigen::MatrixXd::Zero(N(),fTmplDistN.size());
   else {
      std::cout<<"Error! Cannot initialize matrix Y. Nr. of templates: "<<fTmplDistN.size()<<std::endl; 
      exit(1);
   }
   int i = 0;
   for ( auto v : fTmplDistN )  {
      if ( v.second.rows() != Y.col(i).rows() ) {
         std::cout<<"Error! Size of (template) vector is incompatible with size of template matrix Y."<<std::endl; 
         std::cout<<"       Y.rows: "<<Y.col(i).rows()<<"\t, y.rows: "<<v.second.rows()<<endl;
         exit(1);
      }
      if ( fA.rows() > 0 )  Y.col(i++) = fA * v.second;
      else                  Y.col(i++) =      v.second;
   }
   return Y;
}


// -------------------------------------------------------------------- //
Eigen::MatrixXd LTF::Calc_X() const {
   // Calculate matrix X, representing the templates
   if ( fA.rows() == 0 )  // no response matrix, therefore X==Y
      return Eigen::MatrixXd();
   Eigen::MatrixXd X = Eigen::MatrixXd(fA.cols(),fTmplDistN.size());
   int i = 0;
   for ( auto [MM,y] : fTmplDistN )  {
      X.col(i++) =  y;
   }
   return X;
}



// -------------------------------------------------------------------- //
std::map<std::string,Eigen::MatrixXd > 
LTF::Calc_dY(std::map<std::string,std::map<vector<double>,Eigen::VectorXd> >& InSysY) const {
   // calcluate template uncertainties for LTF::LiTeFit
   std::map<std::string,Eigen::MatrixXd > retSysY; // return value
   for ( auto [name,mvalvec] : InSysY ) {
      // the number of rows is a bit trick here, if A is used: (Y is actually A*Y, but also A may not be set yet)
      int nn = fTmplDistN.begin()->second.rows();
      retSysY[name] = Eigen::MatrixXd::Zero(nn, fTmplDistN.size()); // (nn,fTmplDistN.size())
      int i = 0;
      for ( auto [mvs,yy] : fTmplDistN )  { // matrix Y is ordered like fTmplDistN
         if ( mvalvec.count(mvs) ) {
            retSysY[name].col(i) = mvalvec[mvs]; // if A is used: sysY are sysX (but note: Y=A*X)
         }
         i++; // next row
      }
   }
   return retSysY;
}


// -------------------------------------------------------------------- //
//!  \brief Run the linear template fit and return the LiTeFit object with 
//!  all results.
const LTF::LiTeFit& LTF::DoLiTeFit() {
   // run the Linear template fit.
   // --- set LTF input and apply fit range
   SetLiTeFitInput();
   fLTF.DoLiTeFit();
   return fLTF;
}

// -------------------------------------------------------------------- //
//!  \brief Run the template fit using the iterative newton algorithm
//!         and a 2nd order model approximation
const LTF::LiTeFit& LTF::DoIterativeFitNewton(int nIter, int nPol, int nInfrc) {
   SetLiTeFitInput();
   fLTF.DoIterativeFitNewton(nIter,nPol,nInfrc);
   return fLTF;
}


// -------------------------------------------------------------------- //
//              Implementations for class LTF::LiTeFit
// -------------------------------------------------------------------- //
//! calculate matrix called M^+
//! M^+ = ( M^T  W  M )^-1 M^T  W
//! with W=1, since we use unweighted linear regression
//! 
//! exponential Gamma-factor is applied for each fit-parameter
//!
//! The function returns the g-inverse for n-th-order regression:
//!    powmax = 1:   linear regression
//!    powmax = 2:   second order polynomial
//!    etc....
//! powmax=1 is the default
//! 
Eigen::MatrixXd LTF::LiTeFit::Mc(int powmax, int interference) const { 
   int nfitpar = this->M.cols()-1;
   if ( this->Gamma.size() ) { // print-out
      for ( int j = 0 ; j<nfitpar ;j++ ) {
         if ( this->Gamma[j] == 0. ) { // info-message
            cout<<"Error.   Gamma factor for parameter "<<j<<" cannot be zero!"<<endl;
            exit(1); }
         if ( this->Gamma[j] != 1. )  // info-message
            cout<<"Info.    Applying Gamma factor for "<<j
                <<"'s fit parameter. Gamma="<<this->Gamma[j]<<endl;
      }
   }
   int ncols = 1 + nfitpar*powmax + interference*((nfitpar*nfitpar-nfitpar)/2);
   //cout<<"LTF::Mc(). nfitpar="<<nfitpar<<"\tpowmax="<<powmax<<"\tintfrce="<<interference<<"\tncols="<<ncols<<endl;
   Eigen::MatrixXd Mtmp = Eigen::MatrixXd::Zero( M.rows(), ncols );
   int icol = 0;
   for ( int i = 0 ; i<M.rows() ;i++ )  Mtmp(i,icol) = 1; // unity's in first row
   icol++;
   for ( int ipow = 1 ; ipow <= powmax ; ipow++ ) {// power of alpha: alpha^1, alpha^2, ...
      for ( int j = 1 ; j<M.cols() ;j++ ) { // each fit parameter
         double gamma = this->Gamma.size() ? this->Gamma[j-1] : 1;
         for ( int i = 0 ; i<M.rows() ;i++ ) {
            Mtmp(i,icol) = pow(M(i,j),gamma*ipow);
         }
         icol++;
      }
   }

   // interference terms
   if ( interference != 0 && interference != 1 ) {
      cout<<"Error in LTF::LiTeFit::Mc()! interference terms only up to n=1 is implemented."<<endl; exit(1);
   }
   if ( interference > 0 ) { // should become a for-loop for all interference terms.
      for ( int j0 = 0 ; j0<nfitpar ; j0++ ) {
         for ( int j1 = j0+1 ; j1<nfitpar ; j1++ ) {
            double gamma0 = this->Gamma.size() ? this->Gamma[j0] : 1;
            double gamma1 = this->Gamma.size() ? this->Gamma[j1] : 1;
            for ( int i = 0 ; i<M.rows() ;i++ ) {
               Mtmp(i,icol) = pow(M(i,1+j0),gamma0) * pow(M(i,1+j1),gamma1);
            }
            icol++;
         }
      }
   }      

   // // apply gamma factor
   // if ( this->Gamma.size() ) {
   //    for ( int j = 1 ; j<M.cols() ;j++ ) {
   //       if ( this->Gamma[j-1] != 1. )  // info-message
   //          cout<<"Info.    Applying Gamma factor for "<<j-1
   //              <<"'s fit parameter. Gamma="<<this->Gamma[j-1]<<endl;
   //       for ( int i = 0 ; i<M.rows() ;i++ ) { 
   //          Mtmp(i,j) = pow(Mtmp(i,j),this->Gamma[j-1]);
   //       }
   //    }
   // }
   Eigen::MatrixXd MT = Mtmp.transpose().eval();
   // if ( (MT*Mtmp).determinant() < 1.e-8 ) {
   //    cout<<"Warning! Determinant of (M^T*M) is very small! Problem may not be solvable."<<endl;
   //    cout<<"         Possibly, not enough templates were provided."<<endl;
   //    cout<<"Printing matrix M (reference point) "<<endl<<M<<endl;
   //    cout<<"Printing matrix (M^T*M):"<<endl<<(MT*Mtmp)<<endl<<endl;
   // }
   return (MT*Mtmp).inverse()*MT;
}


// -------------------------------------------------------------------- //
//!< \brief calculated approximate matrix M^+ from non-linear matrix M^+  
//!< Return matrix Mc that is obtained from first-order approximation
//!< of non-linear Mc
//!< aexp specifies the expansion point
Eigen::MatrixXd LTF::LiTeFit::McApprox(int mPolN, int mOrdInfrc, const Eigen::VectorXd& aexp0) const {
   // --- vector mbar, and matrix Mtilde
   const int nPar    = M.cols()-1;
   Eigen::MatrixXd MXYtmp  = this->Mc(mPolN,mOrdInfrc); // Matrix M^+, using higher-order polynoms and interference terms
   Eigen::MatrixXd MXYtil = MXYtmp.block(1,0,nPar,MXYtmp.cols()).transpose(); // matrix: Mtilde
   Eigen::VectorXd mXYbar = MXYtmp.row(0).transpose(); // vector: mbar
   
   // --- linearly approximation for non-linear model (if mPolN>1)
   //     n-parameter fit 2nd-order approx. with and w/o interference
   if ( mPolN == 2 ) {
      Eigen::VectorXd aexp = aexp0;
      // --- apply gamma factor first to expansion values
      for ( size_t j = 0 ; j<this->Gamma.size() ;j++ ) 
         aexp(j) = pow( aexp(j), this->Gamma[j] ); 

      for ( int ia =0  ; ia<nPar; ia++ ) {
         mXYbar         -=    MXYtmp.row(1+nPar+ia).transpose() * aexp(ia,0) * aexp(ia,0);
         MXYtil.col(ia) += 2.*MXYtmp.row(1+nPar+ia).transpose() * aexp(ia,0);
      } 
      // interference terms
      if ( mOrdInfrc==1 ) {
         int icol=1+nPar*nPar;
         for ( int j0 = 0 ; j0<nPar ; j0++ ) {
            for ( int j1 = j0+1 ; j1<nPar ; j1++ ) {
               //cout<<"Evaluating column "<<icol<<" for interference terms."<<endl;
               mXYbar         -=  1.*MXYtmp.row(icol).transpose() * aexp(j0,0) * aexp(j1,0); // constant term from interference aexp0*aexp1
               MXYtil.col(j0) +=     MXYtmp.row(icol).transpose() * aexp(j1,0); // aexp1*aexp
               MXYtil.col(j1) +=     MXYtmp.row(icol).transpose() * aexp(j0,0); // aexp0*a1
               icol++;
            }
         }
      }      
   }
   // build matrix Mc again:
   Eigen::MatrixXd McReturn = Eigen::MatrixXd::Zero(1+nPar,mXYbar.rows());
   McReturn.row(0) = mXYbar.transpose();
   McReturn.block(1,0,nPar,MXYtmp.cols()) = MXYtil.transpose();
   return McReturn;
}


// -------------------------------------------------------------------- //
//! \brief Calculate inverse (covariance) matrix
//! Input are all covariance matrices in V
//! Note: at least one (uncorrelated) uncertainty must be 
//! present in the fit
Eigen::MatrixXd LTF::LiTeFit::W() const {
   // Calculate inverse covariance matrix W
   if ( Vs.size() == 0 ) { 
      std::cout<<"Error! No uncorrelated or stat. uncertainty specified. Please use AddError()."<<std::endl;
      exit(1);
   }
   // Eigen::MatrixXd Vsum = Eigen::MatrixXd::Zero(Vs.begin()->second.rows(),Vs.begin()->second.cols());
   // for ( auto& [name,V] : Vs ) Vsum += V;
   Eigen::MatrixXd Vsum = LTF::VSum(Vs);
   // if ( Vsum.determinant() < 1.e-8 ) {
   //    cout<<"Warning! Determinant of matrix V is very small! "<<endl;
   //    cout<<"Printing matrix Vsum"<<endl<<Vsum<<endl<<endl;
   // }
   return Vsum.inverse().eval();
}


// -------------------------------------------------------------------- //
void LTF::LiTeFit::PrintShort() const {
   for ( int i = 0 ; i<(int(M.cols())-1) ; i++ ) {
      std::cout<<"  LTF. Result["<<i<<"]:   "<<ahat(i)<<"  +/- " 
               << ahat_errorFit(i)<<" (fit)  +/- " << ahat_errorExt(i)<< " (ext)"
               <<std::endl;
   }
}

// -------------------------------------------------------------------- //
void LTF::LiTeFit::PrintFull() const {
   // print all results
   std::cout<<std::endl;
   std::cout<<"------------------------------------------------------------------------"<<std::endl;
   std::cout<<"        Result of the "<<(M.cols()-1)<<"-dimensional linear template fit"<<std::endl;
   std::cout<<"        ```````````````````````````````````````````````"<<std::endl;
   PrintShort();
   int nPar = M.cols()-1;
   std::cout<<std::endl;
   for ( int i = 0 ; i<nPar ; i++ ) {
      printf("  LTF. Result[%d]:   %f\n",i,ahat(i));
      for ( auto& [name,V] : Vsource )         printf("                             +/- % 8.6f (%s)\n", std::sqrt(V(i,i)), name.c_str());
      for ( auto& [name,v] : DeltaSys )        printf("                             +/- % 8.6f (%s)\n", v(i), name.c_str());
      std::cout<<std::endl; 
      for ( auto& [name,V] : VsourceExt )      printf("                             +/- % 8.6f (%s)\n", std::sqrt(V(i,i)), name.c_str());
      for ( auto& [name,v] : DeltaSysExt )     printf("                             +/- % 8.6f (%s)\n", v(i), name.c_str());
      for ( auto& [name,s] : DeltaSysY )       printf("                             +/- % 8.6f (%s)\n", s(i), name.c_str());
      for ( auto& [name,s] : DeltaSysA )       printf("                             +/- % 8.6f (%s)\n", s(i), name.c_str());
   }
   Eigen::MatrixXd V = VFit();
   if ( ahat.rows() > nPar ) {
      std::cout<<std::endl;
      std::cout<<"  Nuisance parameters       ";
      printf("          % 5.3f  +/-  %5.3f  (%s)\n",ahat(nPar),sqrt(V(nPar,nPar)),Sys[0].first.c_str() );
      for ( int i = nPar+1 ; i<int(ahat.rows()) ; i++ )
         printf("                                      % 5.3f  +/-  %5.3f  (%s)\n",ahat(i),sqrt(V(i,i)),Sys[i-nPar].first.c_str() );
   }
   if ( ahat.rows() > 1 ) {
      std::cout<<std::endl;
      streamsize tmpprec = std::cout.precision();
      std::cout<<"  Correlation matrix "<<std::endl;
      std::cout.precision(3);
      std::cout<<Cov_to_Cor(V)<<std::endl;
      std::cout.precision(tmpprec);
   }
   std::cout<<std::endl;
   std::cout<<"  Chi^2                                "<< chisq << "  +/- "<<chisq_error<<std::endl;
   std::cout<<"  Chi^2/ndf                            "<< ((Dt.rows()-nPar)!=0 ? chisq / (Dt.rows()-nPar) : 0 ) <<std::endl;
   std::cout<<"  Chi^2 for each template              "<< chisq_y.transpose() <<endl;
   //std::cout<<"  Partial Chi^2:                  ";
   int kk=0;
   for ( auto [name,c] : chisq_part ) {
      if ( (kk++ %4) == 0  )  // line break
         cout<<"  Partial chi^2                    ";
      printf("   %5.2f (%-10s",c,(name+")").c_str());
      //cout<<"   "<<c<<" ("<<name<<")";
      if ( kk%4 == 0 ) std::cout<<std::endl;
   }
   std::cout<<std::endl;
   

   // --- pol2-fit to chi2-'parabola
   std::cout<<std::endl;
   std::cout<<"  Alternative result from pol2-fit to chi^2 values of the templates"<<std::endl;
   for ( size_t i = 0 ; i<size_t(achk.rows()) ; i++ ) {
      std::cout<<"     Result(pol2-chi2)["<<i<<"]              "<<achk(i)<<"   +/- "<<achk_errorFit(i)<<" (fit)"<<endl;
   }
   std::cout<<"     Chi^2(min)                        "<<achk_chisq<<endl;


   // --- linearized second-order model approximation 
   std::cout<<std::endl;
   std::cout<<"  Result from linearized second-order model approximation"<<std::endl;
   // for ( int i = 0 ; i< nPar ; i++ )
   //    std::cout<<"     Result(pol2)["<<i<<"]                "<<ahat21(i)<<endl;
   // if ( ahat21.rows() > nPar ) {
   //    std::cout<<"     Nuisance parameters ";
   //    printf("           % 5.3f   (%s)\n",ahat21(0),Sys[0].first.c_str() );
   //    for ( int i = nPar+1 ; i<int(ahat.rows()) ; i++ )
   //       printf("                                    % 5.3f   (%s)\n",ahat21(i),Sys[i-nPar].first.c_str() );
   // }
   for ( int i = 0 ; i< nPar ; i++ ) 
      std::cout<<"     EDM (Taylor) ["<<i<<"]                  "<<ahat21(i)-ahat(i)<<endl;
   std::cout<<endl;

   // --- EDM (Newton distance)
   std::cout<<"  Expected distance to minimum"<<std::endl;
   for ( int i = 0 ; i< nPar ; i++ )
      std::cout<<"     EDM (Newton) ["<<i<<"]                  "<<EDM21(i)<<endl;
   std::cout<<"     EDM Delta-chi^2                   "<<EDM21_chisq<<endl;
   // one could also print nuisance parameters here...
   std::cout<<endl;

   std::cout<<std::endl;
   std::cout<<"  'Uncertainties' on M from second-order model approximation"<<std::endl;
   for ( int i = 0 ; i< nPar ; i++ ) 
      std::cout<<"     sigma_m["<<i<<"]                       "<<delta_m21(i)<<endl;
   for ( int i = nPar ; i<int(ahat.rows()) ; i++ )
      printf("     sigma_m[%-24s %4.2f\n",(Sys[i-nPar].first+"]").c_str(),delta_m21(i));
   
   // --- settings
   std::cout<<std::endl;
   printf("  Settings\n");
   printf("     Number of data points             %ld\n",Dt.rows());
   printf("     Number of templates               %ld\n",Y.cols());
   printf("     Number of fit parameters          %ld\n",(M.cols()-1));
   for ( int i = 1 ; i<M.cols() ; i++ )
      std::cout<<"     Template reference values ["<<i-1<<"]     " << M.col(i).transpose() <<endl;
   printf("     Number of nuisance param.         %ld\n",ahat.rows()-(M.cols()-1));
   int i=0;
   for ( auto& g : Gamma )printf("     Gamma[%d]                          %.3f\n",i++,g);
   printf("     LogNormal                         %s\n",(LogNormal?"true":"false"));
   printf("     Response matrix                   %s\n",(A.rows()?"true":"false"));
   std::cout<<std::endl;
   std::cout<<"  please cite: (to be published)."<<std::endl;
   std::cout<<"------------------------------------------------------------------------"<<std::endl;
   std::cout<<std::endl;
}


// -------------------------------------------------------------------- //
void LTF::LiTeFit::ApplyGamma(std::map<std::string,Eigen::VectorXd>& Delta) const {
   // apply gamma factor to the uncertainties
   for ( auto& [name,s] : Delta ) {    // systematic shifts
      for ( size_t i = 0 ; i<Gamma.size() ;i++ ) { // apply gamma factor
         s(i) = s(i) / ( pow(ahat[i],Gamma[i]-1) * Gamma[i] );
      }
   }
}

// -------------------------------------------------------------------- //
void LTF::LiTeFit::ApplyGamma(std::map<std::string,Eigen::MatrixXd>& VMat) const {
   // apply gamma factor to the uncertainties
   for ( auto& [name,V] : VMat ) {    // covariance matrices
      for ( size_t i = 0 ; i<Gamma.size() ;i++ ) { // apply gamma factor
         for ( size_t j = 0 ; j<Gamma.size() ;j++ ) { // apply gamma factor
            V(i,j) = V(i,j ) / ( pow(ahat[i],Gamma[i]-1) * Gamma[i] ) / ( pow(ahat[i],Gamma[i]-1) * Gamma[i] );
         }
      }
   }
}


// -------------------------------------------------------------------- //
//! Modify input paramters in order to obtain relative uncertainties
void LTF::LiTeFit::ApplyLogNormalDistribution() {
   // Use LogNormal distributed uncertainties, i.e. normal-distibuted
   // *relative*  uncertainties
   if ( LogNormal ) {
      std::cerr<<"Error in ApplyLogNormalDistribution()! "
               <<"Relative uncertainties were already calculated. "
               <<"Please re-initialize LiTeFit first."<<std::endl;
      exit(1);
   }
   
   LogNormal = true;
   std::cout<<"Info.    Use log-normal probability density functions, i.e. relative uncertainties."<<std::endl;

   // reference scale of the uncertainties
   Eigen::VectorXd ref = ErrorScale.rows() > 0 ? ErrorScale : Dt;

   // absolute -> relative uncertainties
   for ( int i = 0 ; i<ref.rows(); i++ ) {
      if ( ref(i) <= 0 ) { // sanity check
         if ( ref(i) == 0 ) std::cerr<<"Error! Data value cannot be zero when working with relative uncertainties."<<std::endl;
         if ( ref(i) <  0 ) std::cerr<<"Error! Data value cannot be negative when working with relative uncertainties (not implemented)."<<std::endl;
         exit(1);
      }
      for ( int j = 0 ; j<ref.rows(); j++ )  {
         for ( auto& [name,V] : this->Vs  )    V(i,j) = V(i,j) / ref(i) / ref(j);
         for ( auto& [name,V] : this->VsExt  ) V(i,j) = V(i,j) / ref(i) / ref(j);
      }
      for ( auto& [name,s] : this->Sys )     s(i) = s(i)/ref(i);
      for ( auto& [name,s] : this->SysExt )  s(i) = s(i)/ref(i);
   }
   
   for ( auto& [name,dA] : this->SysA ) {   
      if ( this->corrSys[name] == 1. ) {
         this->SysAX.push_back(std::make_pair(name,(A+dA) * X - (A*X)));
      }
      else if ( this->corrSys[name] == 0. ) { // generate a AX matrix for each element of A
         for ( int i = 0 ; i<dA.rows() ;i++ ) {
            for ( int t = 0 ; t<dA.cols() ;t++ ) {
               if ( dA(i,t) != 0 ) {
                  auto AdA = A;
                  AdA(i,t) += dA(i,t);
                  this->SysAX.push_back(make_pair(name,(AdA) * X - (A*X)));
               }
            }
         }
      }
      else {
         std::cerr<<"Error!  correlation must be 1 or 0 for dA."<<std::endl; 
         exit(1);
      }
      this->SysA.erase(name);
   }

   // template uncertainties
   for ( auto& [name,dY] : this->SysY ) {   
      for ( int i = 0 ; i<dY.rows() ; i++ ) {
         for ( int j = 0 ; j<dY.cols() ; j++ ) {
            if ( A.rows() > 0 ) dY(i,j) = dY(i,j) / X(i,j);
            else                dY(i,j) = dY(i,j) / Y(i,j);
         }
      }
   }
   for ( auto& [name,dAX] : this->SysAX ) {   
      for ( int i = 0 ; i<dAX.rows() ; i++ ) {
         for ( int j = 0 ; j<dAX.cols() ; j++ ) {
            dAX(i,j) = dAX(i,j) / Y(i,j);
         }
      }
   }

   // last: data and Y
   for ( int i = 0 ; i<ref.rows(); i++ ) {
      Dt(i) = log(Dt(i));
      for ( int j = 0 ; j<Y.cols(); j++ ) {
         if ( Y(i,j) <= 0 ) { // sanity check
            if ( Y(i,j) == 0 ) std::cerr<<"Error! Template value cannot be zero when working with relative uncertainties."<<std::endl;
            if ( Y(i,j)  < 0 ) std::cerr<<"Error! Template value cannot be negative when working with relative uncertainties (possible, but not implemented)."<<std::endl;
            exit(1);
         }
         Y(i,j)  = log(Y(i,j));
      }
      if ( y0new.rows() > 0 )
	 y0new(i) = log(y0new(i));
   }
}


// -------------------------------------------------------------------- //
//! \brief Rescale uncertainties. See LTF::RescaleErrors() for more details
void LTF::LiTeFit::RescaleInputErrors(const Eigen::VectorXd& NewScale ) {
   if ( LogNormal ) { // sanity
      std::cout<<"Warning. Error rescaling is not reasonable when working with relative uncertainties (LogNormal). Skipping...."<<std::endl;
      return;
   }
   if ( NewScale.rows() != Dt.rows() ) {
      std::cerr<<"Error!   LTF::LiTeFit::RescaleInputErrors(). Incompatible size of the input vector for rescaling. Size="
               <<NewScale.rows()<<", but data has: "<<Dt.rows()<<". Exiting."<<std::endl;
      exit(1);      
   }

   // --- set present scale
   if ( ErrorScale.rows() == 0 ) { // assuming data as default
      ErrorScale = Dt;
   }
   const auto& ref = ErrorScale;

   // // --- rescale
   // for ( int i = 0 ; i<ref.rows(); i++ ) {
   //    cout<<"Rescaling uncertainties in bin  "<<i<<", using factor: "<<NewScale(i)/ref(i)<<endl;
   //    for ( int j = 0 ; j<ref.rows(); j++ )  {
   //       for ( auto& [name,V] : this->Vs  )    V(i,j) = V(i,j) / ref(i) / ref(j) * NewScale(i) * NewScale(j);
   //       for ( auto& [name,V] : this->VsExt  ) V(i,j) = V(i,j) / ref(i) / ref(j) * NewScale(i) * NewScale(j);
   //    }
   //    for ( auto& [name,s] : this->Sys )     s(i) = s(i)/ref(i) * NewScale(i);
   //    for ( auto& [name,s] : this->SysExt )  s(i) = s(i)/ref(i) * NewScale(i);
   // }


   // --- rescale
   for ( auto& [name,V] : this->Vs  )  {
      if ( name.find("stat") != string::npos || name.find("Stat") !=string::npos  || name.find("ncor") !=string::npos ) { // skip statistical uncertainty
         cout<<"Info. LiTeFit::RescaleInputErrors. Skipping uncertianty source '"<<name<<"' because it is a statistical uncertainty."<<endl;
         continue;
      }
      for ( int i = 0 ; i<ref.rows(); i++ ) {
         for ( int j = 0 ; j<ref.rows(); j++ )  {
            V(i,j) = V(i,j) / ref(i) / ref(j) * NewScale(i) * NewScale(j);
         }
      }
   }

   for ( int i = 0 ; i<ref.rows(); i++ ) {
      cout<<"Rescaling uncertainties in bin  "<<i<<", using factor: "<<NewScale(i)/ref(i)<<endl;
      for ( int j = 0 ; j<ref.rows(); j++ )  {
         for ( auto& [name,V] : this->VsExt  ) V(i,j) = V(i,j) / ref(i) / ref(j) * NewScale(i) * NewScale(j);
      }
      for ( auto& [name,s] : this->Sys )     s(i) = s(i)/ref(i) * NewScale(i);
      for ( auto& [name,s] : this->SysExt )  s(i) = s(i)/ref(i) * NewScale(i);
   }

   
   // store new scale:
   ErrorScale = NewScale;
}



// -------------------------------------------------------------------- //
//!  \brief Compute the value of the next Newton step, when using higher-order model approximation
Eigen::VectorXd LTF::LiTeFit::ComputeNewtonEstimator(int mPolN, int nInfrc, const Eigen::VectorXd& ain) /*const*/ {
   
   // --- check input
   if ( mPolN !=1  && mPolN !=2 ) { 
      cout<<"Error. only polynomials with power 1 and 2 are implemented in this function."<<endl; 
      exit(1);
   }
   if ( nInfrc !=0  && nInfrc !=1 ) { 
      cout<<"Error. only interference terms up to first order are implemented in this function."<<endl; 
      exit(1);
   }

   // local parameters
   const int nPar    = M.cols()-1;
   const int nSys    = this->Sys.size();
   const int nParEps = nPar + this->Sys.size(); //nPar+S.cols();
   Eigen::VectorXd aexp = ain;
   for ( int k = 0 ; k<nPar ; k++ ) {
      double gamma = this->Gamma.size() ? this->Gamma[k] : 1;  // apply gamma factor first to ahat
      aexp(k) = pow(ain(k),gamma);
   }

   // sanity check
   if ( mPolN > 1 && (aexp.rows() != nPar && aexp.rows() != nParEps) ) {
      cout<<"Error. aexp is needed as expansion point for linear approximation."<<endl;
      cout<<"aexp.cols: "<<aexp.rows()<<"\t# parameters: "<<nPar<<",\t and including nuisance param.: "<<nParEps<<endl;
      exit(1);
   }


   // calculate matrix of uncertainies: W, S and E
   Eigen::MatrixXd W  = this->W();
   Eigen::MatrixXd E  = Eigen::MatrixXd::Zero(nParEps,nParEps);
   Eigen::MatrixXd Es = Eigen::MatrixXd::Zero(this->Sys.size(),this->Sys.size());
   Eigen::MatrixXd S(Dt.rows(),this->Sys.size());
   {
      int j=0;
      for ( const auto& [name,s] : this->Sys ) {
         S.col(j) = s;
         if ( this->SysIsConstr.at(name) == LTF::Uncertainty::Constrained )  {
            E(j+nPar,j+nPar) = 1;
            Es(j,j)   = 1;
         }
         j++;
      }
   }

   // regression matrix
   Eigen::MatrixXd McN = this->Mc(mPolN,nInfrc);

   // --- residuum
   Eigen::VectorXd a0pow = this->a0pow(mPolN,nInfrc,aexp);
   Eigen::VectorXd r0    = Dt - Y*McN.transpose()*a0pow; // epsilons=0

   Eigen::VectorXd eps;
   if ( nSys>0 ) {
      if ( aexp.rows() == nPar+nSys ) 
         eps = aexp.tail(nSys);
      else {
         eps = +1.*(S.transpose()*W*S + Es).inverse()*S.transpose()*W*(Dt-Y*McN.transpose()*a0pow);
         aexp.resize(nPar+nSys);
         aexp.tail(nSys) = eps;
      }
      //r0  = r0; // eps=0 as expansion point
      r0  = r0 - S*eps;
   }

   // --- first derivative of model
   vector<Eigen::VectorXd> r1(nPar+nSys); // first derivatives of residuum (without Y)
   for ( int k1 = 0 ; k1<nPar ; k1++ )
      r1[k1]  = -1. * McN.row(1+k1)     .transpose();               // init & first order terms
   for ( int k1 = 0, kinfrc=0 ; mPolN>1 && k1<nPar ; k1++ ) {        // only for mPolN>1
      r1[k1] += -2. * McN.row(1+nPar+k1).transpose() * aexp(k1);    // second order terms
      for ( int k2 = k1+1 ; nInfrc && k2<nPar ; k2++ ) {
         r1[k1] += -1. * McN.row(1+2*nPar+kinfrc).transpose() * aexp(k2); // interference terms
         r1[k2] += -1. * McN.row(1+2*nPar+kinfrc).transpose() * aexp(k1); // interference terms
         kinfrc++;
      }
   }
   for ( int k = 0 ; k<nPar ; k++ )  // times Y
      r1[k] = Y*r1[k];
   for ( int l = 0 ; l<nSys ; l++ ) { // nuisance parameters
      r1[nPar+l] = -1.*S.col(l);
   }

   // ---- vector with first derivatives at ahat
   Eigen::VectorXd Grad(nPar+nSys);
   for ( int k = 0 ; k<nPar+nSys ; k++ ) {
      Grad(k) = 2*r1[k].transpose() * W * r0;
   }
   for ( int l = 0 ; l<nSys ; l++ ) {
      Grad(nPar+l) += 2*eps(l); // comment, if eps=0 as expansion point
   }
   //cout<<" # Grad:"<<endl<<Grad<<endl;
         
   // --- second derivative of model   
   // all second derivatives with nuisance parameters are zero.
   vector<vector<Eigen::VectorXd > > r2(nPar+nSys); // second derivatives of residuum (without Y)
   for ( int k = 0 ; k<nPar+nSys ; k++ ) { // initialize r2
      r2[k].resize(nPar+nSys);
      for ( int k2 = 0 ; k2<nPar+nSys ; k2++ ) {
         r2[k][k2] = Eigen::VectorXd::Zero(Dt.rows());
      }
   }
   for ( int k = 0, kinfrc=0 ; mPolN>1 &&  k<nPar ; k++ ) { // compute second deriv's (mPolN>1)
      r2[k][k] = -2.*Y*McN.row(1+nPar+k).transpose(); // times Y
      for ( int k2 = k+1 ; nInfrc && k2<nPar ; k2++ ) {
         r2[k][k2] = -1.*Y*McN.row(1+2*nPar+kinfrc).transpose(); // infrc,  // times Y
         r2[k2][k] = r2[k][k2];
         kinfrc++;
      }
   }
   
   // --- Hesse matrix
   Eigen::MatrixXd Hesse(nPar+nSys,nPar+nSys);
   Eigen::MatrixXd HesseApprox(nPar+nSys,nPar+nSys);
   for ( int k1 = 0 ; k1<nPar+nSys ; k1++ ) { 
      for ( int k2 = k1 ; k2<nPar+nSys ; k2++ ) {
         Hesse(k1,k2)       = ( 2.*(r2[k1][k2]).transpose() * W * r0
                                + 2.*(r1[k1]).transpose() * W * (r1[k2]))(0,0);
         HesseApprox(k1,k2) =  2.*((r1[k1]).transpose() * W * r1[k2])(0,0);
         if ( k1!=k2 ) Hesse(k2,k1)       = Hesse(k1,k2);
         if ( k1!=k2 ) HesseApprox(k2,k1) = HesseApprox(k1,k2);
      }
   }
   for ( int l = 0 ; l<nSys ; l++ ) {
      Hesse(nPar+l,nPar+l)       += 2. ;
      HesseApprox(nPar+l,nPar+l) += 2.;
   }   
   //cout<<"Hesse: "<<endl<<Hesse<<endl; 
   
   // --- calculate final results
   Eigen::MatrixXd diff = -1.*Hesse.inverse()*Grad;
   Eigen::MatrixXd diffApprox = -1.*HesseApprox.inverse()*Grad; // this is equivalent to the taylor approximation
   double          chi2aprox = (Grad.transpose()*diff + 0.5*diff.transpose()*Hesse*diff)(0,0);
   if ( nInfrc==1 ) EDM21_chisq = chi2aprox; // keep it for printing

   Eigen::VectorXd ahatnew = aexp+diff;
   // --- apply 'inverse' gamma factor 
   for ( size_t j = 0 ; j<this->Gamma.size() ;j++ ) 
      ahatnew(j) = pow( ahatnew(j), 1./this->Gamma[j] ); 

   return ahatnew;
   
} 



// -------------------------------------------------------------------- //
//! \brief Solve the (linear) template fit
//!  The equations of the linear template fit are applied
//!  and the best estimator(s) are determined.           
//!  With default arguments, the linear template fit     
//!  is applied.                                         
//!  Non-linear approximations of the model are          
//!  performed for mPolN > 1, while only mPolN=2 is      
//!  implemented (quadratic function)                    
//!  In addition interference terms can be considered with mOrdInfrc=1     
Eigen::VectorXd LTF::LiTeFit::ComputeLinearTemplateFit ( int mPolN, int mOrdInfrc,  const Eigen::VectorXd& aexp ) const {
   // solve the linear template fit
   // when using mPolN>0 and/or mOrdInrc, the
   // model is approximated with a quadratic function
   // andthe approximated linearly at  the expansion point aexp.
   
   // --- check input
   if ( mPolN !=1  && mPolN !=2 ) { 
      cout<<"Error. only polynomials with power 1 and 2 are implemented in this function."<<endl; 
      exit(1);
   }
   if ( mOrdInfrc !=0  && mOrdInfrc !=1 ) { 
      cout<<"Error. only interference terms up to first order are implemented in this function."<<endl; 
      exit(1);
   }

   // local parameters
   const int nPar    = M.cols()-1;
   const int nParEps = nPar + this->Sys.size(); //nPar+S.cols();
   // sanity check
   if ( mPolN > 1 && (aexp.rows() != nPar && aexp.rows() != nParEps) ) {
      cout<<"Error. aexp is needed as expansion point for linear approximation."<<endl;
      cout<<"aexp.cols: "<<aexp.rows()<<"\t# parameters: "<<nPar<<",\t and including nuisance param.: "<<nParEps<<endl;
      exit(1);
   }

   Eigen::MatrixXd McApprox = this->McApprox(mPolN, mOrdInfrc, aexp);
   Eigen::VectorXd mXYbar   = McApprox.row(0).transpose();
   Eigen::MatrixXd MXYtil   = McApprox.bottomRows(McApprox.rows()-1).transpose(); // ytil(Dt.rows(),nPar)


   // calculate matrix of uncertainies: W, S and E
   Eigen::MatrixXd W = this->W();
   Eigen::MatrixXd E = Eigen::MatrixXd::Zero(nParEps,nParEps);
   Eigen::MatrixXd S(Dt.rows(),this->Sys.size());
   int j=0;
   for ( const auto& [name,s] : this->Sys ) {
      S.col(j) = s;
      if ( this->SysIsConstr.at(name) == LTF::Uncertainty::Constrained )  
         E(j+nPar,j+nPar) = 1;
      j++;
   }

   // useful matrix (Y*Mtilde S)
   Eigen::MatrixXd YtSXY = Eigen::MatrixXd::Zero(Dt.rows(),nParEps); // 'n x 1'  matrix. elements: ytilde
   YtSXY.block(0,0,Dt.rows(),nPar) = Y*MXYtil;
   YtSXY.block(0,nPar,Dt.rows(),S.cols()) = S;
      

   // --- do the template fit!
   Eigen::MatrixXd ahatnew = ((YtSXY).transpose()*W*YtSXY + E).inverse() * (YtSXY).transpose() * W * (Dt-Y*mXYbar);

      // --- apply 'inverse' gamma factor 
   for ( size_t j = 0 ; j<this->Gamma.size() ;j++ ) 
      ahatnew(j) = pow( ahatnew(j), 1./this->Gamma[j] ); 

   return ahatnew;
}


// -------------------------------------------------------------------- //
//! \brief Calculate the prediction as it is used by the fit for a given parameter value
Eigen::VectorXd LTF::LiTeFit::CalculatePrediction(Eigen::VectorXd v ) const {
   if ( v.rows() < M.cols()-1 ) {
      std::cerr<<"Error in CalculatePrediction(). Number of given values ("<<v.rows()
               <<") must be equivalent to the number of fit parameters (here: "<<M.cols()-1
               <<")."<<std::endl;
      exit(1);
   }
   if ( v.rows() > M.cols()-1 ) {
      std::cout<<"Warning in CalculatePrediction(). Too many values given. Ignoring unrelevant values"<<std::endl;
   }
   for ( size_t k = 0 ; k<Gamma.size() ;k++ )  // apply gamma factor first to ahat
      v(k) = pow( v(k), this->Gamma[k] ); 
   
   Eigen::VectorXd prediction = ybar + Ytil*v;
   if ( LogNormal ) { 
      for ( int i=0 ; i<prediction.rows() ; i++ ) {
         prediction(i) = exp(prediction(i));
      }
   }
   return prediction;
}


// -------------------------------------------------------------------- //
//! \brief Compute vector a to be multiplied with Mc
//! \param nPow     maximum power of polynomial model approximation (nPow=1,2)
//! \param nInfrc   maximum power of interference terms (nInfrc=0,1) (only applicable for nPow>=2)
//! \param ahat     input values for parameters alpha
//!
//!  Return vector may be similar to [and includes the gamma factors]
//!  a0pow = ( 1  a0  a1   a0^2  a1^2  a1*a2 )  
//!
Eigen::VectorXd LTF::LiTeFit::a0pow(int nPow, int nInfrc, const Eigen::VectorXd& ahat) const {
   // check input
   if ( nPow != 1 && nPow!=2 ) {
      cout<<"Error. Invalid input in LTF::LiTeFit::a0pow. nPow="<<nPow<<endl;
      exit(1);
   }
   if ( nInfrc != 0 && nInfrc!=1 ){
      cout<<"Error. Invalid input in LTF::LiTeFit::a0pow. nInfrc="<<nInfrc<<endl;
      exit(1);
   }

   int nPar = this->M.cols()-1;
   int nSize = 1+nPow*nPar;
   if ( nInfrc==1 ) nSize += (nPar*nPar-nPar)/2;
   Eigen::VectorXd a0pow(nSize);
   a0pow(0) = 1.;
   for ( int k = 0, kinfrc = 0 ; k<nPar ; k++ ) {
      double gamma = this->Gamma.size() ? this->Gamma[k] : 1;  // apply gamma factor first to ahat
      a0pow(1+k) = pow(ahat(k),gamma);
      if ( nPow==2 ) a0pow(1+nPar+k) = pow(ahat(k),gamma*2);
      for ( int k2 = k+1 ; nInfrc && k2<nPar ; k2++ ) {
         double gamma2 = this->Gamma.size() ? this->Gamma[k2] : 1;  // apply gamma factor first to ahat
         a0pow(1+2*nPar) = pow(ahat(k),gamma)*pow(ahat(k2),gamma2);
         kinfrc++;
      }
   }
   return a0pow;
}


// -------------------------------------------------------------------- //
//! \brief Run the linear template fit iteratively with non-linear model approximation
//! \brief using the Newton algorithm
//! 
//! \param nIter       Number of iterations beyond first linear template fit
//! \param nPol        polynomial of model approximation (default =2, 1: linear template fit)
//! \param nInfrc      polynomial of interference terms for nPol>=2 (default =1)
//!
double LTF::LiTeFit::DoIterativeFitNewton(int nIter, int nPol, int nInfrc) {
   Eigen::VectorXd ahat0 = ComputeLinearTemplateFit();
   for ( int kk = 0 ; kk<nIter ; kk++ ) {
      ahat0  = ComputeNewtonEstimator(2, 1, ahat0);
   }
   return DoLiTeFit(nPol, nInfrc, ahat0);
}
   

// -------------------------------------------------------------------- //
//! \brief Run the linear template fit iteratively with non-linear model approximation
//! \brief using a taylor expansion of the non-linear model
//! 
//! \param nIter       Number of iterations beyond first linear template fit
//! \param StepSize    Learning rate, 0<StepSize<=1
//! \param nPol        polynomial of model approximation (default =2, 1: linear template fit)
//! \param nInfrc      polynomial of interference terms for nPol>=2 (default =1)
//!
double LTF::LiTeFit::DoIterativeFitTaylor(int nIter, double StepSize, int nPol, int nInfrc) {
   // run the linear template fit

   Eigen::VectorXd ahat0 = ComputeLinearTemplateFit();
   for ( int kk = 0 ; kk<nIter-1 ; kk++ ) {
      ahat0 = ahat0 + (ComputeLinearTemplateFit(nPol,nInfrc,ahat0)-ahat)*StepSize;
   }
   return DoLiTeFit(nPol, nInfrc, ahat);
}

// -------------------------------------------------------------------- //
//! \brief Run linear template fit (LiTeFit,LTF) with previously set input parameters
//! This function fills all member variables of LiTeFit.
//! All results are collected in the member variables
//! For a printout of the results, call PrintFull() or PrintShort()
double LTF::LiTeFit::DoLiTeFit() {
   // run the linear template fit
   return this->DoLiTeFit(1,0,Eigen::VectorXd());
}


// -------------------------------------------------------------------- //
//! \brief Run linear template fit (LiTeFit,LTF) with previously set input parameters
//! This function fills all member variables of LiTeFit.
//! All results are collected in the member variables
//! For a printout of the results, call PrintFull() or PrintShort()
double LTF::LiTeFit::DoLiTeFit(int mPolN, int mOrdInfrc,  const Eigen::VectorXd& aexp) {
   // run the linear template fit
   // in this function, only few sanity checks are performed
   // since we expect that all member variables are initialized
   // with reasonable values (e.g. matrix-sizes, etc...)

   //Eigen::MatrixXd  Mc = this->Mc(1,0) ; // linear template fit: use Mc(1,0);
   Eigen::MatrixXd  Mc = this->McApprox(mPolN,mOrdInfrc,aexp) ; // get linearised (non-linear) model
   Eigen::MatrixXd  W  = this->W()  ;
   
   // --- useful quantities
   const int nPar           = M.cols()-1;
   const int nParEps        = nPar + this->Sys.size();
   // Eigen::MatrixXd Mtil     = Mc.bottomRows(Mc.rows()-1).transpose(); // ytil(Dt.rows(),nPar)
   // Eigen::MatrixXd Ytil     = Y * Mtil; // ytil(Dt.rows(),nPar)
   Eigen::MatrixXd Mtil     = Mc.bottomRows(Mc.rows()-1).transpose(); // ytil(Dt.rows(),nPar)
   this->Ytil     = Y * Mtil; // ytil(Dt.rows(),nPar)
   this->ybar     = Y * Mc.row(0).transpose();
   // get ybar: if there is a new nominal template, then recalcluate Ytil accordingly.
   if ( y0new.rows() == Dt.rows() && a0new.rows() == Mc.rows()-1 ) {
      cout<<"old ybar: "<<endl;
      cout<<ybar<<endl;
      cout<<"Info. considering a new nominal template, but keeping the derivatives from the template matrix."<<endl;
      this->ybar = y0new - Ytil * a0new;
      cout<<"new ybar: "<<endl;
      cout<<ybar<<endl;
   }
   else if (  y0new.rows() !=0 || a0new.rows() != 0 ) {
      cerr<<"Error! A new nominal template is provided, but the size of the vector is incompatible with other input parameters."<<endl;
      cout<<"y0new.rows():: "<<y0new.rows()<<endl;
      cout<<"y0new.cols():: "<<y0new.cols()<<endl;
      cout<<" a0new.rows(): "<<a0new.rows()<<endl;
      cout<<" a0new.cols(): "<<a0new.cols()<<endl;
      cout<<"Mc.cols():     "<<Mc.cols()<<endl;
      cout<<"Mc.rows():     "<<Mc.rows()<<endl;
      exit(1);
   }
   Eigen::VectorXd d_ybar   = Dt - ybar;
   
   // construct matrix YtS 
   Eigen::MatrixXd E        = Eigen::MatrixXd::Zero(nParEps,nParEps);
   Eigen::MatrixXd Es       = Eigen::MatrixXd::Zero(this->Sys.size(),this->Sys.size());
   Eigen::MatrixXd YtS  (Dt.rows(),nParEps); // 'n x 1'  matrix. elements: ytilde
   Eigen::MatrixXd S    (Dt.rows(),this->Sys.size());// useful below
   int i=0; int j=0;
   for ( ; i<nPar ; i++ ) YtS.col(i)   = Ytil.col(i);
   for ( const auto& [name,s] : this->Sys ) {
      YtS.col(i) = s;
      S.col(j) = s;
      if ( this->SysIsConstr[name] == LTF::Uncertainty::Constrained ) {
         E(i,i)    = 1;
         Es(j,j)   = 1;
      }
      i++;
      j++;
   }

   // --- calculate result for mt and nuisance parameters epsilon
   Eigen::MatrixXd YtSTW        = YtS.transpose()*W;
   Eigen::MatrixXd D = YtSTW*YtS + E;   // 'discriminant' matrix
   if ( D.determinant() < 1.e-8 ) {
      cout<<"Warning! Determinant of matrix D is very small! "<<endl;
      cout<<"Printing matrix D"<<endl<<D<<endl<<endl;
      cout<<"Printing matrix M^T "<<endl<<Mc<<endl<<endl;
   }
   Eigen::MatrixXd Dinv         = D.inverse();
   this->F           = Dinv * (YtSTW); // the LTF-master-matrix
   this->ahat        = F * d_ybar; // this is the final fit-result for the multivariate LTF !
   this->a0          = Gamma.empty() ? ahat(0,0) : pow( ahat(0,0), 1./this->Gamma[0] ) ; // this is the very final fit-result for the one-parameter fit !


   cout<<"LTF::DoLiTeFit(). Fit done. Result:  "<<this->a0<< std::flush;

   
   // ---- cross check, alternative implementation
   /* {
      Eigen::MatrixXd YM      = Y * Mc.transpose();
      YM.conservativeResize(YM.rows(), nParEps+1); // n * nParEps+1, add syst. shifts
      Eigen::MatrixXd E2      = Eigen::MatrixXd::Zero(nParEps+1,nParEps+1);
      for ( size_t i = 0 ; i<this->Sys.size() ; i++ ) { // add syst. shifts
      YM.col(i+nPar+1) = this->Sys[i].second;
      if ( SysIsConstr[this->Sys[i].first] == LTF::Uncertainty::Constrained ) E2(i+nPar+1,i+nPar+1)=1;
      }
      Eigen::MatrixXd YMWYM = YM.transpose()*W*YM;
      YMWYM.row(0)  = Eigen::VectorXd::Zero(nPar+1);  // patch first row for constant term ybarWd
      YMWYM(0,0)    = (ybar.transpose()*W*Dt)(0,0); // patch first elem.
      cout<<"Result: "<<endl<<((YMWYM+E2).inverse()*YM.transpose()*W*Dt).bottomRows(nParEps)<<endl; // 0th element is '1' and could serve a s a cross check
      } 
   */


   // --- best-fit theoretical prediction  (note: ahat is still with gamma!)
   this->TheoFit   = ybar + Ytil*ahat.topRows(M.cols()-1);
   if ( LogNormal ) { 
      for ( int i=0 ; i<TheoFit.rows() ; i++ ) TheoFit(i) = exp(TheoFit(i));
   }

   // --- Error propagation systematic errors
   for ( const auto& [name,V] : this->Vs )     // uncorr./cov uncertainties
      Vsource[name]       = F * V * F.transpose();
   for ( const auto& [name,s] : this->Sys )     // corr. uncertainties
      DeltaSys[name]      = F * s;
   for ( const auto& [name,V] : this->VsExt )  // external uncorr./cov uncertainties
      VsourceExt[name]    = F * V * F.transpose();
   for ( const auto& [name,s] : this->SysExt )     // corr. uncertainties
      DeltaSysExt[name]   = F * s;

   // ---- propagate template uncertainties
   if ( this->SysY.size() ) {
      // --- partial derivatives
      const auto& dYtmp = SysY.begin()->second; // just for the size, which is Y(i*j) or X(t*j)
      map<int,map<int,Eigen::MatrixXd> > dYijs;
      for ( int i = 0 ; i<dYtmp.rows() ; i++ ) {
         for ( int j = 0 ; j<dYtmp.cols() ; j++ ) {
            // note: for performance critical applications the following calculation could be
            // improved by considering the many '0's only at the end
            Eigen::MatrixXd d1ij =  Eigen::MatrixXd::Zero(dYtmp.rows(),dYtmp.cols()); // derivative. i.e. 1_(ij)
            d1ij(i,j) = 1.;
            if ( A.rows() > 0 )  { // with A
               if ( LogNormal ) { // with A & LogNormal
                  d1ij(i,j) = X(i,j);
                  d1ij = (A*d1ij).eval();
                  for ( int k = 0 ; k<d1ij.rows(); k ++ ) {
                     for ( int l = 0 ; l<d1ij.cols(); l ++ ) {
                        if ( d1ij(k,l) ) d1ij(k,l) *= 1. / exp(Y(k,l));
                     }
                  }      
               }         
               else { // with A
                  d1ij = (A*d1ij).eval();
               }
            }
            auto d1M = d1ij*Mtil; // note: this can be improved since 1*Mt selects only a single row of Mt
            Eigen::MatrixXd d1M0 = Eigen::MatrixXd::Zero(YtS.rows(),YtS.cols());
            for ( int iii=0; iii<nPar ; iii++ ) d1M0.col(iii) = d1M.col(iii); // note: this can be improved since 1*Mt selects only a single row of Mtf
            // we use: YtSTW  & Dinv from above
            Eigen::MatrixXd dYij_1_2 = d1M0.transpose() * W * d_ybar; // 1st & 2nd term
            Eigen::MatrixXd dYij_3   = YtSTW * d1ij * Mc.row(0).transpose(); // 3rd term 
            Eigen::MatrixXd dYij_4   = (YtSTW * d1M0 + d1M0.transpose() * W * YtS )  * ahat; //4th term (from D^-1)
            // partial derivative:
            dYijs[i][j] = Dinv * ( dYij_1_2 - dYij_3 - dYij_4 ) ;
         }
      }

      // -- propagate template uncertainty (corr. & uncorr.)
      for ( const auto& [name,dY] : this->SysY ) {   
         DeltaSysY[name] =  Eigen::VectorXd::Zero(ahat.rows());
         double corr = corrSys[name];
         for ( int i = 0 ; i<dY.rows() ; i++ ) {
            for ( int j = 0 ; j<dY.cols() ; j++ ) {
               auto dYijDYij = dYijs[i][j].col(0) * dY(i,j); // error propagation
               if ( corr == 1 )   // correlated
                  DeltaSysY[name]  +=  dYijDYij; 
               else {             // uncorrelated (apply sqrt below!)
                  for ( int kk = 0 ; kk<nPar ; kk++ ) {
                     DeltaSysY[name](kk)  +=  dYijDYij(kk,0) * dYijDYij(kk,0); // coefficient wise pow(2);
                  }
               }
            }
         }
         if ( corr == 0 ) 
            DeltaSysY[name] =  DeltaSysY[name].cwiseSqrt().eval();
      }
   }
   
   // ---- propagate response matrix uncertainties
   if ( this->SysA.size() ) { // only if !LogNormal
      // --- partial derivatives
      map<int,map<int,Eigen::MatrixXd> > dAits;
      int N_i = SysA.begin()->second.rows();
      int N_t = SysA.begin()->second.cols();
      for ( int i = 0 ; i<N_i ; i++ ) {
         for ( int t = 0 ; t<N_t ; t++ ) {
            Eigen::MatrixXd d1it =  Eigen::MatrixXd::Zero(N_i, N_t); // derivative. i.e. 1_(it)
            d1it(i,t) = 1.;
            // note: for performance critical applications the following calculation could be
            // improved by considering the many '0's only at the end and constructing directly d1XM and dD/dA
            auto d1XM = d1it*X*Mtil;
            Eigen::MatrixXd d1XM0 = Eigen::MatrixXd::Zero(YtS.rows(),YtS.cols());
            for ( int iii=0; iii<nPar ; iii++ ) d1XM0.col(iii)   = d1XM.col(iii); // note: this can be improved since 1*Mt selects only a single row of Mtf
            // we use: YtSTW from above
            Eigen::MatrixXd dAit_1_2 = d1XM0.transpose() * W * d_ybar; // 1st & 2nd term
            Eigen::MatrixXd dAit_3   = YtSTW * d1it * X * Mc.row(0).transpose(); // 3rd term 
            Eigen::MatrixXd dAit_4   =  ( YtSTW * d1XM0 +  d1XM0.transpose() * W * YtS) * ahat; //4th term (from D^-1)
            // partial derivative:
            dAits[i][t] = Dinv * ( dAit_1_2 - dAit_3 - dAit_4 ) ;
         }
      }

      // -- propagate uncertainties of response matrix (corr. & uncorr.)
      for ( const auto& [name,dA] : this->SysA ) {   
         DeltaSysA[name] =  Eigen::VectorXd::Zero(ahat.rows());
         double corr = corrSys[name];
         for ( int i = 0 ; i<N_i ; i++ ) {
            for ( int t = 0 ; t<N_t ; t++ ) {
               auto dAitDAit = dAits[i][t].col(0) * dA(i,t); // error propagation
               if ( corr == 1 )   // correlated
                  DeltaSysA[name]  +=  dAitDAit; 
               else {             // uncorrelated (apply sqrt below!)
                  for ( int kk = 0 ; kk<nPar ; kk++ ) {
                     DeltaSysA[name](kk)  +=  dAitDAit(kk,0) * dAitDAit(kk,0); // coefficient wise pow(2);
                  }
               }
            }
         }
         if ( corr == 0 ) 
            DeltaSysA[name] =  DeltaSysA[name].cwiseSqrt().eval();
      }
      //for ( auto [name,err] : DeltaSysY ) cout<<name<<"\t"<<err<<endl;
   }

   // ---- propagate matrix uncertainties when using log-normal
   if ( this->SysAX.size() ) { // if log-normal
      // --- partial derivatives
      auto dYtmp = SysAX.begin()->second; // just for the size, which is Y(i*j) or X(t*j)
      map<int,map<int,Eigen::MatrixXd> > dYijs;
      for ( int i = 0 ; i<dYtmp.rows() ; i++ ) {
         for ( int j = 0 ; j<dYtmp.cols() ; j++ ) {
            // note: for performance critical applications the following calculation could be
            // improved by considering the many '0's only at the end
            Eigen::MatrixXd d1ij =  Eigen::MatrixXd::Zero(dYtmp.rows(),dYtmp.cols()); // derivative. i.e. 1_(ij)
            d1ij(i,j) = 1.;
            //if ( A.rows() > 0 ) d1ij = (A*d1ij).eval();
            auto d1M = d1ij*Mtil; // note: this can be improved since 1*Mt selects only a single row of Mt
            Eigen::MatrixXd d1M0 = Eigen::MatrixXd::Zero(YtS.rows(),YtS.cols());
            for ( int iii=0; iii<nPar ; iii++ ) d1M0.col(iii) = d1M.col(iii); // note: this can be improved since 1*Mt selects only a single row of Mtf
            // we use: YtSTW  & Dinv from above
            Eigen::MatrixXd dYij_1_2 = d1M0.transpose() * W * d_ybar; // 1st & 2nd term
            Eigen::MatrixXd dYij_3   = YtSTW * d1ij * Mc.row(0).transpose(); // 3rd term 
            Eigen::MatrixXd dYij_4   = (YtSTW * d1M0 + d1M0.transpose() * W * YtS )  * ahat; //4th term (from D^-1)
            // partial derivative:
            dYijs[i][j] = Dinv * ( dYij_1_2 - dYij_3 - dYij_4 ) ;
         }
      }

      // -- propagate template uncertainty (corr. & uncorr.)
      // note in SysAX each uncorr. uncertainties 'name' may be present multiple times
      for ( const auto& [name,dY] : this->SysAX ) {   
         DeltaSysA[name] =  Eigen::VectorXd::Zero(ahat.rows()); // initialize first here!
      }
      for ( const auto& [name,dY] : this->SysAX ) {   
         double corr = corrSys[name];
         Eigen::VectorXd dSysAtmp = Eigen::VectorXd::Zero(ahat.rows());
         for ( int i = 0 ; i<Y.rows() ; i++ ) {
            for ( int j = 0 ; j<Y.cols() ; j++ ) {
               auto dYijDYij = dYijs[i][j].col(0) * dY(i,j); // error propagation
               if ( corr == 1 )   // correlated
                  DeltaSysA[name]  +=  dYijDYij; 
               else   // built uncorr. sum below, and apply final sqrt below!
                  dSysAtmp += dYijDYij;
            }
         }
         if ( corr == 0 ) {
            for ( int kk = 0 ; kk<nPar ; kk++ ) {
               DeltaSysA[name](kk)  +=  dSysAtmp(kk,0) * dSysAtmp(kk,0); // coefficient wise pow(2);
            }
         }
      }
      for ( auto[ name, ds ] : DeltaSysA ) { // only once at the end!
         if ( corrSys[name] == 0 ) 
            DeltaSysA[name] =  DeltaSysA[name].cwiseSqrt().eval();
      }
      //for ( auto [name,err] : DeltaSysY ) cout<<name<<"\t"<<err<<endl;
   }
   // --- end of error propagation !
   //  ---------------------------------------------------------------- //


   //  ---------------------------------------------------------------- //
   // correct for gamma factor(s)
   //  ---------------------------------------------------------------- //
   Eigen::VectorXd ahatgamma = ahat; // keep it
   for ( size_t k = 0 ; k<Gamma.size() ;k++ )  // apply gamma factor first to ahat
      this->ahat(k) = pow( ahat(k), 1./this->Gamma[k] ); 
   ApplyGamma(Vsource);
   ApplyGamma(VsourceExt);
   ApplyGamma(DeltaSys);
   ApplyGamma(DeltaSysExt);
   ApplyGamma(DeltaSysY);
   ApplyGamma(DeltaSysA);


   // --- calculate prediction (note: ahat is now without gamma!)
   //TheoFit = CalculatePrediction(ahat);

   //  ---------------------------------------------------------------- //
   // --- calculate total uncertainties for convenience
   //  ---------------------------------------------------------------- //
   ahat_errorFit  = Eigen::VectorXd::Zero(nPar);
   for ( int i = 0 ; i<nPar ; i++ ) {
      for ( const auto& s : this->Vsource  )    this->ahat_errorFit(i) += s.second(i,i);
      for ( const auto& s : this->DeltaSys )    this->ahat_errorFit(i) += s.second(i)*s.second(i);
   }
   ahat_errorExt  = Eigen::VectorXd::Zero(nPar);   
   for ( int i = 0 ; i<nPar ; i++ ) {
      for ( const auto& s : this->VsourceExt  ) this->ahat_errorExt(i) += s.second(i,i);
      for ( const auto& s : this->DeltaSysExt ) this->ahat_errorExt(i) += s.second(i)*s.second(i);
   }
   ahat_errorTmpl = Eigen::VectorXd::Zero(nPar);
   for ( int i = 0 ; i<nPar ; i++ ) {
      for ( const auto& [name,dsY] : this->DeltaSysY  ) this->ahat_errorTmpl(i) += dsY(i)*dsY(i);
   }
   ahat_errorRespMat = Eigen::VectorXd::Zero(nPar);
   for ( int i = 0 ; i<nPar ; i++ ) {
      for ( const auto& [name,dsA] : this->DeltaSysA  ) this->ahat_errorRespMat(i) += dsA(i)*dsA(i);
   }
   
   this->ahat_errorFit   = ahat_errorFit.cwiseSqrt().eval();
   this->ahat_errorExt   = ahat_errorExt.cwiseSqrt().eval();
   this->ahat_errorTmpl  = ahat_errorTmpl.cwiseSqrt().eval();
   this->ahat_errorRespMat = ahat_errorRespMat.cwiseSqrt().eval();
   this->a0_errorFit      = ahat_errorFit(0);
   this->a0_errorExt      = ahat_errorExt(0);
   this->a0_errorTmpl     = ahat_errorTmpl(0);
   this->a0_errorRespMat  = ahat_errorRespMat(0);

   //  ---------------------------------------------------------------- //
   // ---  chi^2
   //  ---------------------------------------------------------------- //
   //chisq = -42;//((Dt - Ymbar - Ymtil*pow(fLTF.mt,fDelta)).transpose() * W * (Dt - Ymbar - Ymtil*pow(fLTF.mt,fDelta)))(0,0);      
   Eigen::VectorXd aGamma = Eigen::VectorXd::Zero(Gamma.size());
   for ( size_t i = 0 ; i<Gamma.size() ; i++ )  aGamma(i) = pow(ahat(i),Gamma[i]);
   Eigen::VectorXd resid = Dt - ybar - Ytil * aGamma;
   double sumeps2 = 0;
   for (  int i = M.cols()-1 ; i <ahat.rows() ; i++ ) {
      sumeps2 += ahat(i) * ahat(i); // ahat(1) is the nuisance parameter epsilon
      resid = (resid - YtS.col(i)*ahat(i)).eval(); // it is minus s*epsilon here!
   }
   this->chisq = resid.transpose() * W * resid + sumeps2; // Gamma ?
   // --- uncertainty of chi2
   Eigen::VectorXd dChidd = Eigen::VectorXd(Dt.size()); //uncertainty propagation to chi2 
   chisq_error = 0;
   for ( int i=0 ; i<dChidd.rows() ;i++ ) dChidd(i) = 2.*(W.row(i) * resid)(0,0); // = (W.row(i) * resid + (resid).transpose() * W.col(i))(0,0);
   //for ( const auto& [n,V] : this->Vs  )   chisq_error += (dChidd.transpose()*V*dChidd)(0,0) ; 
   chisq_error += (dChidd.transpose()*LTF::VSum(Vs)*dChidd)(0,0) ; 
   for ( const auto& [n,s] : this->Sys )   chisq_error += pow((dChidd.transpose()*s)(0,0),2) ; 
   chisq_error = sqrt(chisq_error);

   //  ---------------------------------------------------------------- //
   // --- partial chi^2 values
   //  ---------------------------------------------------------------- //
   for ( const auto& [name,V] : this->Vs )     // uncorr./cov uncertainties
      chisq_part[name]       = resid.transpose() * W * V * W * resid ;
   for ( int i = nPar ; i<int(ahat.rows()) ; i++ )  // corr. uncertainties (this->Sys)
      chisq_part[Sys[i-nPar].first] = ahat(i) * ahat(i);
   // -- include also  external syst. Note: the fit-chi2 will be preserved only with above's uncertainties
   for ( const auto& [name,V] : this->VsExt )     // external uncertainties
      chisq_part[name]       = resid.transpose() * W * V * W * resid ;
   for ( const auto& [name,s] : this->SysExt )     // external shifts
      chisq_part[name]       = pow( (s.transpose() * W * resid)(0,0),2) ;

   // --- (pseudo-)nuisance paramters for external uncertainties
   if ( this->SysExt.size() ) {
      Eigen::MatrixXd Vtot = LTF::VSum(Vs,Sys);
      Eigen::MatrixXd Vinv = Vtot.inverse();
      for ( const auto& [name,s] : this->SysExt )  {
         double n = (s.transpose()*W*resid)(0,0);
         double dn = sqrt((s.transpose()*Vinv*s)(0,0));
         nuisance_SysExt[name] = make_pair(n,dn);
      }
   }



   //  ---------------------------------------------------------------- //
   //  --- chi^2 for each template
   //  ---------------------------------------------------------------- //
   this->chisq_y = Eigen::VectorXd::Zero(Y.cols());
   Eigen::MatrixXd Fc2s = (S.transpose() * W * S + Es).inverse() * S.transpose() * W;
   for ( int k=0 ; k<Y.cols(); k++ ) {
      Eigen::VectorXd eps = Fc2s * (Dt - Y.col(k));
      Eigen::VectorXd rr = Dt - Y.col(k);
      for ( int is=0 ; is<S.cols(); is++ ) {
         rr = (rr - S.col(is) * eps(is) ).eval();
         chisq_y(k) += eps(is)*eps(is);
      }
      chisq_y(k) = (rr.transpose() * W * rr)(0);
   }


   //  ---------------------------------------------------------------- //
   // --- fit chi2-parabola
   //  ---------------------------------------------------------------- //
   // (note: not tested for multi-dimensional fits!)
   Eigen::MatrixXd M2tmp(chisq_y.rows(), 1+nPar*2); // quadratic regression
   for ( int i = 0 ; i<chisq_y.rows() ;i++ ) {
      M2tmp(i,0) = 1.;
      for ( int j = 0 ; j<nPar ;j++ ) {
         M2tmp(i,1+j*2) = M(i,1+j);
         M2tmp(i,2+j*2) = M(i,1+j)*M(i,1+j);
      }
   }
   auto MT = M2tmp.transpose();
   Eigen::MatrixXd Mc2 = (MT*M2tmp).inverse()*MT;
   Eigen::VectorXd abc = Mc2 * chisq_y; // do regression !
   // collect results...
   achk          = Eigen::VectorXd::Zero(nPar); // initialize result
   achk_errorFit = Eigen::VectorXd::Zero(nPar); // initialize
   achk_chisq    = abc(0);
   for ( int j = 0 ; j<nPar ;j++ ) {
      double d    = -2.*abc(2+2*j);
      achk(j)     = abc(1+2*j) / d;
      achk_chisq += abc(1+2*j) * achk(j) + abc(2+2*j) * achk(j)*achk(j);
   }
   const double delta_chi2 = 1.; // get uncertainty from delta-chi2 criterion
   for ( int j = 0 ; j<nPar ;j++ ) {
      double d    = -2.*abc(2+2*j);
      double stmp = abc(1+2*j)*abc(1+2*j) - 4.*abc(2+2*j) * ( abc(0) - (achk_chisq+delta_chi2) ) ;
      // double Dachk2 = (-abc(1+2*j) - sqrt(stmp) ) / (2.*abc(2+2*j))  -  achk(j);
      //achk_errorFit(j) = fabs(  (-abc(1+2*j) + sqrt(stmp) ) / (2.*abc(2+2*j))  - achk(j) );
      if ( stmp>0) achk_errorFit(j) = fabs( sqrt(stmp)/d );
   }

   //  ---------------------------------------------------------------- //
   // --- non-linear approximations
   //  ---------------------------------------------------------------- //
   this->ahat21 = ComputeLinearTemplateFit(2,1,ahat ); // second order approximation incl. interference terms
   //this->ahat20 = ComputeLinearTemplateFit(2,0,ahat );
   //this->ahat   = ComputeLinearTemplateFit(); // linear template fit
   //this->EDM20  = ComputeNewtonEstimator(2, 0, ahat) - ahat;
   this->EDM21  = ComputeNewtonEstimator(2, 1, ahat) - ahat;

   //  ---------------------------------------------------------------- //
   // --- error propagation of m
   //  ---------------------------------------------------------------- //
   delta_m21 = Eigen::VectorXd::Zero(ahat.rows());
   //if ( nPar == 1 ) 
   { // validated only for 1-parameter fit. Resulting uncertainties in n-parameter fits appear a bit large
      Eigen::MatrixXd McApprox = this->McApprox(2,1, ahat);
      Eigen::MatrixXd dmbar2 = McApprox.row(0).transpose()-Mc.row(0).transpose(); // mbar2 - mbar1
      Eigen::MatrixXd dMtil2 = McApprox.bottomRows(McApprox.rows()-1).transpose()-Mtil;                  // Mtilde2 - Mtilde
      for ( int j = 0 ; j<dmbar2.rows() ; j++ ) {
         //mbar
         Eigen::VectorXd j0 = Eigen::VectorXd::Zero(dmbar2.rows());
         j0(j) = 1.;
         delta_m21 +=  (-F*Y*j0) * dmbar2(j,0);
         //mtilde
         for ( int k = 0 ; k<nPar ; k++ ) {
            // Eigen::MatrixXd Yj2 = Eigen::MatrixXd::Zero(Y.rows(),nPar); // Yj2 = Y*j1
            // Yj2.col(k) += Y.col(j);
            Eigen::MatrixXd YtSXY = Eigen::MatrixXd::Zero(Dt.rows(),nParEps); // 'n x 1'  matrix. elements: ytilde
            YtSXY.col(k) += Y.col(j);
            YtSXY.block(0,nPar,Dt.rows(),S.cols()) = S;
            delta_m21 += (Dinv*( (YtSXY).transpose()*W*(d_ybar) - ((YtSXY).transpose()*W*YtS + YtSTW*YtSXY) * ahatgamma )) * dMtil2(j,k);
         }
      }
      // apply gamma factor
      for ( size_t k = 0 ; k<Gamma.size() ;k++ ) 
         delta_m21(k) = delta_m21(k) / ( pow(ahat[k],Gamma[k]-1) * Gamma[k] );
   }

   //  ---------------------------------------------------------------- //
   // --- printout
   //  ---------------------------------------------------------------- //
   std::cout
      //<<"LTF::DoLiTeFit().     done. Result:  "<<mt
      <<"  +/- " << a0_errorFit  << " (fit)"
      <<"  +/- " << a0_errorExt  << " (ext)"
      <<"  +/- " << a0_errorTmpl << " (tmpl)"
      <<"  +/- " << a0_errorRespMat << " (MatA)"
      <<"   [chi^2=" << chisq<<"]"
      <<std::endl;

   return chisq;
}



// -------------------------------------------------------------------- //
//! \brief Set correlated error
//! \param name    Name of the error source
//! \param error   size of the error (or systematic shift), taking care for the sign of the shifts
//! \param constr  treat uncertainty constrained or unconstrained in fit, or only propagate it after the fit ('External')
//!
//! Note: at least one uncorrelated (or covariance) error must be specified.
//! Function used internally for initializing member variables
void LTF::AddCorrelatedError(const std::string& name, vector<double> error, LTF::Uncertainty constr) {
   if ( constr != LTF::Uncertainty::External) {
      if ( fUseNuisance ) { // use nuisance parameter 
         fSys.push_back(make_pair(name, Eigen::Map<Eigen::VectorXd>(&error[0], error.size()) ));
         fSysIsConstr[name] = constr;
      }
      else { // use covariance matrix instead
         size_t n = error.size();
         std::vector<std::vector<double> > cov(n);
         for ( size_t i = 0 ; i<n ;i++ ) {
            cov[i].resize(n);
            for ( size_t j = 0 ; j<n ;j++ ) {
               cov[i][j] = error[i]*error[j];
            }
         }
         AddError(name, cov, constr);         
      }
   } 
   else 
      fSysExt[name] = Eigen::Map<Eigen::VectorXd>(&error[0], error.size());
}



// -------------------------------------------------------------------- //
//! \brief calculate covariance matrix from an uncertainty
//static 
Eigen::MatrixXd LTF::GetV_Delta(const Eigen::VectorXd& d, double corr) { 
   if ( corr == 1 ) return d * d.transpose();
   else if ( corr == 0 ) {
      Eigen::MatrixXd ret = Eigen::MatrixXd::Zero(d.rows(),d.rows());
      for ( size_t i = 0 ; i<size_t(d.rows()) ; i++ ) ret(i,i) = d(i)*d(i);
      return ret;
   }
   else {
      cout<<"Error. corr must be 1 or 0. corr="<<corr<<endl; exit(1);
   }
}


// -------------------------------------------------------------------- //
//! \brief Calculate covariance matrix as sum of individual matrices and systematic shifts
//! \param Vs        Covariance matrices, i.e. uncorr. and/or cov. uncertainties
//! \param DeltaSys  correlated systematic uncertainty
//static 
Eigen::MatrixXd LTF::VSum( const std::map<std::string,Eigen::MatrixXd>& Vs, const std::map<std::string,Eigen::VectorXd>& DeltaSys ) {
   if ( Vs.empty() && DeltaSys.empty() ) std::cerr << "LTF::VSum(). No input given!"<<endl;
   int n = Vs.size() ? Vs.begin()->second.rows() : Vs.begin()->second.rows();
   Eigen::MatrixXd ret = Eigen::MatrixXd::Zero( n, n);
   for ( const auto& V : Vs )       ret += V.second;
   for ( const auto& v : DeltaSys ) ret += LTF::GetV_Delta(v.second,1);
   return ret;
}



// -------------------------------------------------------------------- //
//! \brief Calculate covariance matrix as sum of individual matrices and systematic shifts
//! \param Vs        Covariance matrices, i.e. uncorr. and/or cov. uncertainties
//! \param Sys       Correlated systematic uncertainties
//static 
Eigen::MatrixXd LTF::VSum( const std::vector<std::pair<std::string,Eigen::MatrixXd > >& Vs, const std::vector<std::pair<std::string,Eigen::VectorXd > >& Sys ) {
   int n = Vs.size() ? Vs.begin()->second.rows() : 0;
   Eigen::MatrixXd ret = Eigen::MatrixXd::Zero( n, n);
   if ( Vs.empty() ) {
      std::cerr << "LTF::VSum(). No input given!"<<endl;
      return ret;
   }
   for ( const auto& V : Vs )       ret += V.second;
   for ( const auto& v : Sys )      ret += LTF::GetV_Delta(v.second,1);
   return ret;
}


// -------------------------------------------------------------------- //
//! \brief Calculate correlation matrix from given covariacne matrix V
//static 
Eigen::MatrixXd LTF::Cov_to_Cor(const Eigen::MatrixXd& V) {
   if ( V.rows() != V.cols() ) {
      std::cerr<<"Error! LTF::Cov_to_Cor(). This is not a diagonal matrix."<<std::endl;
      exit(1);
   }
   size_t n = V.cols();
   Eigen::MatrixXd Rho = Eigen::MatrixXd::Zero(n,n);
   std::vector<double> delta;
   for ( size_t i = 0 ; i<n ; i++ ) delta.push_back(sqrt(V(i,i)));
   for ( size_t i = 0 ; i<n ; i++ ) {
      for ( size_t j = 0 ; j<n ; j++ )  {
         Rho(i,j) = V(i,j)/(delta[i]*delta[j]);
      }
   }
   return Rho;
}


// -------------------------------------------------------------------- //
//! \brief Calculate an Eigen-matrix from std::vector<std::vector<>>. T=double
//static 
Eigen::MatrixXd LTF::Std_to_Mat(const std::vector<std::vector<double> >& A) {
   int nj = A.size() > 0 ? A[0].size() : 0;
   Eigen::MatrixXd ret = Eigen::MatrixXd(A.size(), nj);
   for ( size_t i = 0 ; i<A.size(); i++ ) 
      for ( size_t j = 0 ; j<A[i].size(); j++ ) 
         ret(i,j) = A[i][j];
   return ret;
}

