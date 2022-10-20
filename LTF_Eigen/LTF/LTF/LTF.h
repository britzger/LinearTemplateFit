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
#ifndef LTF_h
#define LTF_h

// --------------------------------------------------------------- //
/*
 *
 * LTF
 *
 * The Linear Template Fit implemented using the Eigen library
 *
 */
// --------------------------------------------------------------- //

#include "Eigen/Sparse"
#include "Eigen/Dense"
#include <string>
#include <map>
#include <unordered_map>
#include <iostream>
#include <memory>

using namespace std;

/**
 * \class LTF
 *
 * \brief Linear template fit
 * \brief implemented with the linear algebra package eigen
 *
 */
class LTF {

public:
   enum class Uncertainty { Constrained=1, Unconstrained=0, External=-1};
   
public:

   // ------------------------------------------- //
   //! \class LTF::LiTeFit
   //! \brief Small helper class to collect all results, 
   //! All input parameters are set by class LTF
   //! The user can access all input and output quantities
   class LiTeFit {
   public:

      LiTeFit() {}
      ~LiTeFit(){}
      double DoLiTeFit();                                                        //!< \brief Run linear template fit (LiTeFit,LTF) with previously set input parameters
      double DoLiTeFit(int mPolN, int mOrdInfrc, const Eigen::VectorXd& ahat);   //!< \brief Run non-linear template fit with previously set input parameters
      double DoQuadraticTemplateFit(int nIter=3) {return DoIterativeFitNewton(nIter);};        //!< \brief Run a quadratic template fit with nIter iterations
      double DoIterativeFitNewton(int nIter=3, int nPol=2, int nInfrc=1);        //!< \brief Run an non-linear template fit with iterative improvements
      double DoIterativeFitTaylor(int nIter=4, double StepSize=0.6, int nPol=2, int nInfrc=1);   //!< \brief Run an non-linear template fit with iterative improvements
      void ApplyLogNormalDistribution();                                         //!< \brief Assume underlying log-normal distributed probability density functions.
      void RescaleInputErrors(std::vector<double> values ){                      //!< \brief Rescale uncertainties. See LTF::RescaleErrors() for more details
         RescaleInputErrors(Eigen::Map<Eigen::VectorXd>(&values[0],values.size()));
      }
      void RescaleInputErrors(const Eigen::VectorXd& values );                   //!< \brief Rescale uncertainties. See LTF::RescaleErrors() for more details

      // --- input, as used in LiTeFit
      Eigen::VectorXd Dt;                                                        //!< data distribution
      Eigen::MatrixXd M ;                                                        //!< Matrix M (reference values)
      Eigen::MatrixXd Y ;                                                        //!< Matrix Y (templates)
      Eigen::MatrixXd X ;                                                        //!< Matrix X (templates) [only required for error propagation]
      Eigen::MatrixXd A ;                                                        //!< Response matrix (optional) [only used for (relative) error propagation]
      std::vector<std::pair<std::string,Eigen::MatrixXd > > Vs    ;              //!< (uncor./Cov) uncertainties
      std::vector<std::pair<std::string,Eigen::VectorXd > > Sys   ;              //!< (correlated) systematic uncertainties 
      std::map<std::string,Eigen::MatrixXd >                VsExt ;              //!< external (uncor./Cov) uncertainties
      std::map<std::string,Eigen::VectorXd >                SysExt;              //!< external (correlated) systematic uncertainties
      std::map<std::string,LTF::Uncertainty>                SysIsConstr;         //!< uncertainty specification for Sys
      std::map<std::string,Eigen::MatrixXd >                SysY;                //!< uncertainty of the templates
      std::map<std::string,Eigen::MatrixXd >                SysA;                //!< uncertainty of the response matrix
      std::map<std::string,double>                          corrSys;             //!< Correlation coefficient of the systematic uncertainties
      vector<double> Gamma;                                                      //!< Use non-linear Gamma-factor

      // --- output
      double a0 = 0;                                                             //!< Fit result
      double a0_errorFit = 0;                                                    //!< uncertainty from all uncertainties included in fit (no 'external' uncertainties)
      double a0_errorExt  = 0;                                                   //!< uncertainty from 'external' uncertainties
      double a0_errorTmpl = 0;                                                   //!< uncertainty from templates
      double a0_errorRespMat = 0;                                                //!< uncertainty from response matrix
      double chisq = 0;                                                          //!< chi^2
      double chisq_error = 0;                                                    //!< uncertainty of chi^2
      Eigen::MatrixXd F;                                                         //!< Projection matrix of LiTeFit
      Eigen::VectorXd chisq_y;                                                   //!< chi^2 for each template
      std::map<std::string,double> chisq_part;                                   //!< 'partial' chi^2 for each error source
      Eigen::VectorXd ahat;                                                      //!< results (a0, a1, ..., epsilon1, epsilon2,...)
      Eigen::VectorXd ahat_errorFit;                                             //!< Uncertainties included in fit
      Eigen::VectorXd ahat_errorExt;                                             //!< Uncertainties not included in fit
      Eigen::VectorXd ahat_errorTmpl;                                            //!< Uncertainties from templates
      Eigen::VectorXd ahat_errorRespMat;                                         //!< Uncertainties from response matrix
      Eigen::VectorXd TheoFit;                                                   //!< best-fit theory predictions
      std::map<std::string,Eigen::MatrixXd> Vsource;                             //!< Covariance matrices from uncorr. and cov.  uncertainties 
      std::map<std::string,Eigen::VectorXd> DeltaSys;                            //!< Correlated systematic uncertainties
      std::map<std::string,Eigen::MatrixXd> VsourceExt;                          //!< Covariance matrices from 'external' uncorr. and cov.  uncertainties
      std::map<std::string,Eigen::VectorXd> DeltaSysExt;                         //!< 'external' correlated systematic uncertainty
      std::map<std::string,Eigen::VectorXd> DeltaSysY;                           //!< systematic uncertainty from the templates
      std::map<std::string,Eigen::VectorXd> DeltaSysA;                           //!< systematic uncertainty from the response matrix
      std::map<std::string,std::pair<double,double> > nuisance_SysExt;           //!< nuisance parameters (and their uncertainties) for external syst. uncertainties

      // --- result from alternative pol2-fit to chi2-values of the templates
      Eigen::VectorXd achk;                                                      //!< result from polynomial fit to chi2's of templates
      Eigen::VectorXd achk_errorFit;                                             //!< result from polynomial fit to chi2's of templates
      double  achk_chisq;                                                        //!< chisq_min from polynomial fit to chi2's of templates
      
      // --- result from linearized second-order model approximation
      Eigen::VectorXd ahat21;                                                    //!< result from linearized second-order model approximation  (incl. interference)
      Eigen::VectorXd EDM21;                                                     //!< expected distance to minimum (EDM), calculating the newton step for a non-linear model
      double          EDM21_chisq;                                               //!< expected distance to minimum (EDM), calculating the newton step for a non-linear model
      Eigen::VectorXd delta_m21;                                                 //!< uncertainty on M from linearized second-order model  (incl. interference)
      
      // ---- printout, and useful getters
      void PrintShort() const;                                                   //!< \brief Print results (short version)
      void PrintFull()  const;                                                   //!< \brief Print all results. Do not print covariance or correlation matrices though.
      Eigen::MatrixXd Mc(int nPow=1, int nInfrc=0) const;                        //!< \brief calculated matrix M^+ from matrix M
      Eigen::VectorXd a0pow(int nPow, int nInfrc, const Eigen::VectorXd& ahat) const; //!< \brief calculated vector (a0,a0^2,infrc) to be multiplied with Mc
      Eigen::MatrixXd McApprox(int powmax, int intfrc, const Eigen::VectorXd& ahat) const;  //!< \brief calculated approximate matrix M^+ from non-linear matrix M^+
      Eigen::MatrixXd W() const;                                                 //!< \brief calculate inverse covariance matrix from V's
      Eigen::MatrixXd VFit() const { return LTF::VSum(Vsource,DeltaSys); }       //!< \brief calculate covariance matrix of all uncertainties included in fit
      Eigen::MatrixXd VExt() const { return LTF::VSum(VsourceExt,DeltaSysExt); } //!< \brief calculate covariance matrix of all uncertainties included in fit
      bool GetLogNormal() const { return LogNormal; }                            //!< \brief return boolean if fit was performed with relative uncertainties
      Eigen::VectorXd CalculatePrediction(double v ) const {                     //!< \brief Calculate the prediction as it is used by the fit for a given parameter value
         return CalculatePrediction(std::vector<double>{v});}
      Eigen::VectorXd CalculatePrediction(std::vector<double> v ) const {        //!< \brief Calculate the prediction as it is used by the fit for a given parameter value
         return CalculatePrediction(Eigen::Map<Eigen::VectorXd>(&v[0],v.size()));}
      Eigen::VectorXd CalculatePrediction(Eigen::VectorXd v ) const;             //!< \brief Calculate the prediction as it is used by the fit for a given parameter value

      Eigen::VectorXd ComputeLinearTemplateFit(int mPolN = 1, int mOrdInfrc = 0, const Eigen::VectorXd& ahat = Eigen::VectorXd() ) const;  //!  \brief Solve (linear) template fit
      Eigen::VectorXd ComputeNewtonEstimator(int mPolN = 1, int mOrdInfrc = 0, const Eigen::VectorXd& ahat = Eigen::VectorXd() ) /*const*/;  //!  \brief Compute the value of the next Newton step, when using a higher-order model approximation

   protected:
      void ApplyGamma(std::map<std::string,Eigen::VectorXd>& Delta) const;       //!< \brief apply gamma factor to the uncertainties
      void ApplyGamma(std::map<std::string,Eigen::MatrixXd>& VMat) const;        //!< \brief apply gamma factor to the uncertainties
      bool LogNormal  = false;                                                   //!< Use LogNormal distributed uncertainties, i.e. normal-distibuted relative  uncertainties
      std::vector<std::pair<std::string,Eigen::MatrixXd> >     SysAX;            //!< uncertainty of the response matrix when using log-normal
      Eigen::VectorXd ErrorScale;                                                //!< Reference values for error rescaling (default: data)
      Eigen::VectorXd ybar;
      Eigen::MatrixXd Ytil;

   };  // end class LTF::LiTeFitResult
   // ------------------------------------------- //

public:
   LTF() {}
   ~LTF() {}

   // ------------------------------------------- //
   // --- set data
   // ------------------------------------------- //
   //!< \brief Set data distribution
   void SetData(const vector<double>& data) {  
      std::vector<double> tmp = data; // copy to no-const for eigen...
      fDt = Eigen::Map<Eigen::VectorXd>(&tmp[0], tmp.size());
   }
   //!< \brief Set data distribution
   void SetData(size_t n, const double* data) {
      SetData(vector<double>(data,data+n)); 
   }


   // ------------------------------------------- //
   // --- add a new template
   // ------------------------------------------- //
   //! \brief Add new template for a given reference point
   //! \param mval  reference point(s) of the template distribution
   //! \param n     number of bins
   //! \param dist  template distribution
   void AddTemplate(double mval, const vector<double>& dist) {
      AddTemplate(vector<double>{mval},dist);
   }
   //! \brief Add new template for a given reference point
   void AddTemplate(double mval, size_t n , const double* dist) {
      //! Add a new template distribution for a 
      //! given reference values mval
      AddTemplate(mval,std::vector<double>(dist, dist+n));
   }
   //!< \brief Add new template for a given reference point
   void AddTemplate(const vector<double> mvals, size_t n , const double* dist) {
      AddTemplate(mvals,std::vector<double>(dist, dist+n));
   }
   //!< \brief Add new template for a given reference point
   void AddTemplate(const vector<double>& mvals, const vector<double>& dist);


   // ------------------------------------------- //
   // --- set error
   // ------------------------------------------- //
   //! \brief Add a new error source to LTF
   //! \param name    Name of the error source
   //! \param n       number of bins
   //! \param error   size of the error
   //! \param corr    size of correlation (0=uncorrelated, 1=correlated, (0-1) split-up)
   //! \param constr  treat uncertainty constrained or unconstrained in fit, or only propagate it after the fit ('External')
   void AddError(std::string name, size_t n, const double* error, double corr=0., LTF::Uncertainty constr = LTF::Uncertainty::Constrained); 

   //! \brief  Add a new error source to LTF
   void AddError(const std::string& name, const std::vector<double>& error, double corr=0., LTF::Uncertainty constr = LTF::Uncertainty::Constrained) {
      AddError(name,error.size(),&error[0],corr,constr);
   }

   //! \brief Set covariance matrix as error source
   //! \param name  Name of the error source
   //! \param cov   covariance matrix (full-matrix notation, symmetric)
   //! \param constr  treat uncertainty constrained or unconstrained in fit, or only propagate it after the fit ('External')
   //! 
   //! note: covariance matrix is always 'constrained' in fit.
   //! note: Function used interally for initializing member variables
   void AddError(const std::string& name, const std::vector<std::vector<double> >& cov, LTF::Uncertainty constr = LTF::Uncertainty::Constrained);

   //! \brief Set covariance matrix as error source
   //! \param name  Name of the error source
   //! \param cov   covariance matrix (full-matrix notation, symmetric, with relative uncertainties)
   //! \param constr  treat uncertainty constrained or unconstrained in fit, or only propagate it after the fit ('External')
   //! 
   //! note: covariance matrix is always 'constrained' in fit.
   //! note: internally, the relative uncertainties are first multiplied by data, and later on again divided by data.
   void AddErrorRelative(const std::string& name, const std::vector<std::vector<double> >& cov_rel, LTF::Uncertainty constr = LTF::Uncertainty::Constrained);

   //! \brief  Add a new error source to LTF with relative values
   //! \param error   size of the relative error
   //! Pass a relative uncertainty to LTF
   //! Note: This setter function only passes the uncertainty
   //! values as relative uncertainties to LTF, whereas the
   //! probability density function in the linear template fit
   //! can be chosen with LTF::UseLogNormalUncertainties()
   void AddErrorRelative(std::string name, size_t n, const double* error, double corr=0., LTF::Uncertainty constr = LTF::Uncertainty::Constrained) {
      AddErrorRelative(name,std::vector<double>(error,error+n),corr,constr);
   }



   //! \brief  Add a new error source to LTF with relative values
   void AddErrorRelative(const std::string& name, std::vector<double> error, double corr=0., LTF::Uncertainty constr = LTF::Uncertainty::Constrained) {
      if ( fDt.rows() == 0 || size_t(fDt.rows()) != error.size() ) {
         std::cerr<<"Error!  Data distribution not set or its size is inconsistent with input uncertainty vector."<<std::endl;
         exit(1);
      }
      for ( size_t i=0 ; i<error.size( ); i++ ) error[i] *= fDt(i);
      AddError(name,error,corr,constr);
   }

   //! \brief  Add a new error source to LTF
   //! \param error   size of the relative error in percent (%)
   void AddErrorPercent(std::string name, size_t n, const double* error, double corr=0., LTF::Uncertainty constr = LTF::Uncertainty::Constrained) {
      AddErrorPercent(name,std::vector<double>(error,error+n),corr,constr);
   }

   //! \brief  Add a new error source to LTF with relative values
   void AddErrorPercent(const std::string& name, std::vector<double> error, double corr=0., LTF::Uncertainty constr = LTF::Uncertainty::Constrained) {
      if ( fDt.rows() == 0 || size_t(fDt.rows()) != error.size() ) {
         std::cerr<<"Error!  Data distribution not set or its size is inconsistent with input uncertainty vector."<<std::endl;
         exit(1);
      }
      for ( size_t i=0 ; i<error.size( ); i++ ) error[i] *= fDt(i) / 100.;
      AddError(name,error,corr,constr);
   }

   //! \brief Set uncorrelated error
   //! \param name    Name of the error source
   //! \param n       number of bins
   //! \param error2  errors squared
   //! \param constr  treat uncertainty constrained or unconstrained in fit, or only propagate it after the fit ('External')
   //!
   //! Note: at least one uncorrelated (or covariance) error must be specified.
   void AddUncorrelatedErrorSquared(const std::string& name, size_t n, const double* error2, LTF::Uncertainty constr = LTF::Uncertainty::Constrained) {
      if ( constr == LTF::Uncertainty::Unconstrained ) {
         std::cerr<<"Error!  An uncorrelated error cannot be 'unconstrained'. Error name: "<<name<<std::endl; 
         exit(1);
      }
      std::vector<std::vector<double> > cov(n);
      for ( size_t i = 0 ; i<n ;i++ ) {
         cov[i].resize(n);
         cov[i][i] = error2[i];
      }
      AddError(name, cov, constr);
   }

   //! \brief Set correlated error
   //! \param name    Name of the error source
   //! \param error   size of the error (or systematic shift), taking care for the sign of the shifts
   //! \param constr  treat uncertainty constrained or unconstrained in fit, or only propagate it after the fit ('External')
   //!
   //! Note: at least one uncorrelated (or covariance) error must be specified.
   //! Function used internally for initializing member variables
   void AddCorrelatedError(const std::string& name, vector<double> error, LTF::Uncertainty constr = LTF::Uncertainty::Constrained);
   
   // ------------------------------------------- //
   //! \brief Add a new error source of the templates to LTF
   //! Note: a repeated call with an already existent name will
   //!       overwrite a previous uncertainty, or alternatively sets/adds
   //!       a value if the uncertianty was unspecified or zero before.
   //! \param name    Name of the error source
   //! \param mvals   reference point(s) of the template distribution
   //! \param n       number of bins
   //! \param error   size of the error
   //! \param corr    size of correlation (0=uncorrelated, 1=correlated, (0-1) split-up [not implemented])
   void AddTemplateError(const std::string& name, const std::vector<double>& mvals, size_t n, const double* error, double corr=0.);

   //! \brief Add a new error source of the templates to LTF
   void AddTemplateError(const std::string& name, double mval, size_t n, const double* error, double corr=0.) {
      AddTemplateError(name, std::vector<double>{mval}, n, error, corr);
   }

   //! \brief Add a new error source of the templates to LTF
   void AddTemplateError(const std::string& name, const std::vector<double>& mvals, const std::vector<double>& error, double corr=0.) {
      AddTemplateError(name,mvals,error.size(),&error[0],corr);
   }
   //! \brief Add a new error source of the templates to LTF
   void AddTemplateError(const std::string& name, double mval, const std::vector<double>& error, double corr=0.) {
      AddTemplateError(name,std::vector<double>{mval},error.size(),&error[0],corr);
   }

   //! \brief Add a new error source of the templates to LTF
   void AddTemplateErrorSquared(const std::string& name, const std::vector<double>& mvals, size_t n, const double* error2, double corr=0.) {
      AddTemplateErrorSquared(name, mvals, std::vector<double>(error2, error2+n), corr);
   }
   //! \brief Add a new error source of the templates to LTF
   void AddTemplateErrorSquared(const std::string& name, double mval, size_t n, const double* error2, double corr=0.) {
      AddTemplateErrorSquared(name, std::vector<double>{mval}, std::vector<double>(error2, error2+n), corr);
   }

   //! \brief Add a new error source of the templates to LTF
   void AddTemplateErrorSquared(const std::string& name, const std::vector<double>& mvals, const std::vector<double>& error2, double corr=0.) {
      std::vector<double> vect = error2;
      for ( double& e : vect ) e = std::sqrt(e);
      AddTemplateError(name, mvals, vect.size(), &vect[0], corr);
   }

   
   // ------------------------------------------- //
   // --- set a (detector) response matrix [optional]
   // ------------------------------------------- //
   //! \brief set a (detector) response matrix (optional)
   //! A detector response matrix can represent 
   //! resolution and acceptance effects of the experimental
   //! device. The response matrix A 'maps' the predictions
   //! to the detector level:  y=Ax
   //! A response matrix furthermore may account for different
   //! binnings of the templates and the data, or it may
   //! represent background contributions, etc...
   void SetResponseMatrix(const std::vector<std::vector<double> >& A);

   //! \brief set a (detector) response matrix (optional)
   void SetResponseMatrix(const Eigen::MatrixXd& A) { fA = A; }


   // ------------------------------------------- //
   //! \brief Add a new error source of the templates to LTF
   //! Note: a repeated call with an already existent name will
   //!       overwrite a previous uncertainty, or alternatively sets/adds
   //!       a value if the uncertianty was unspecified or zero before.
   //! \param name    Name of the error source
   //! \param error   size of the error of the elements of A
   //! \param corr    size of correlation (0=uncorrelated, 1=correlated, (0-1) split-up [not implemented])
   void AddResponseMatrixError(const std::string& name, const std::vector<std::vector<double> >& errorA, double corr=0.) {
      AddResponseMatrixError(name, LTF::Std_to_Mat(errorA), corr);
   }
   //! \brief Add a new error source of the templates to LTF
   void AddResponseMatrixError(const std::string& name, const Eigen::MatrixXd& errorA, double corr=0.) {
      SetCorrSys(name,corr);
      fSysA[name] = errorA;
   }
   //! \brief Add a new error source of the templates to LTF
   void AddResponseMatrixErrorSquared(const std::string& name, std::vector<std::vector<double> > error2A, double corr=0.) {
      for ( auto& v : error2A )
         for ( double& e2 : v ) e2 = sqrt(e2);
      AddResponseMatrixError(name, LTF::Std_to_Mat(error2A), corr);
   }
   //! \brief Add a new error source of the templates to LTF
   void AddResponseMatrixErrorSquared(const std::string& name, const Eigen::MatrixXd& error2A, double corr=0.) {
      AddResponseMatrixError(name, error2A.cwiseSqrt(), corr);
   }

   // ------------------------------------------- //
   // --- set exponential factor(s) gamma
   // ------------------------------------------- //
   //! \brief set exponential factor(s) gamma
   //! \param gamma  
   //! \param binmax  last bin that is included in the fit
   //! Note: better use: SetGamma({gamma1,gamma2,...,});
   void SetGamma(double gamma = 1) {
      if ( fTmplDistN.empty() ) { 
         std::cerr<<"Warning in SetGamma. Templates should be set first!"<<std::endl;
         fGamma = vector<double>(50,gamma);
      }
      else  {
         if ( fTmplDistN.begin()->first.size() > 1 ) 
            std::cout<<"Warning."
                     <<" SetGamma(double) is ambigous when performing a multivariate LTF."
                     <<" Better use SetGamma(vector<double>)."<<std::endl;
         
         fGamma = vector<double>(fTmplDistN.begin()->first.size(),gamma);
      }
   }
   //! \brief set exponential factor(s) gamma
   void SetGamma(const vector<double>& gamma) {
      fGamma = gamma;
   }

   // ------------------------------------------- //
   // --- fit range
   // ------------------------------------------- //
   //! \brief set fit range, i.e. exclude bins
   //! \param binmin  first bin to be included in the fit
   //! \param binmax  last bin that is included in the fit
   void SetFitRange(int binmin, int binmax) {
      fBinMin = binmin;
      fBinNN  = binmax -binmin+1;
      if ( fBinNN<=0 ) {
         cout<<"Error. No bins left for fit."<<endl;
         exit(1);
      }
   }

   // ------------------------------------------- //
   // --- perform linear template fit
   // ------------------------------------------- //
   //! \brief     Perform linear template fit (LiTeFit,LTF)
   //! \returns   LiTeFit-object with all relevant quantities
   const LiTeFit& DoLiTeFit();
   const LiTeFit& DoIterativeFitNewton(int nIter=3, int nPol=2, int nInfrc=1);

   //!< \brief Run a quadratic template fit with nIter iterations
   //! \returns   LiTeFit-object with all relevant quantities
   const LiTeFit& DoQuadraticTemplateFit(int nIter=3) {return DoIterativeFitNewton(nIter);};        

   // ------------------------------------------- //
   // --- Choose normal or log-normal uncertainties 
   // ------------------------------------------- //
   //! \brief     Select normal distributed or log-normal distributed uncertainties
   //! note: log-normal distributed uncertainties are equivalent
   //! to normal-distributed *relative* uncertainties and can
   //! be considered as a good approximation to a poisson distribution
   void UseLogNormalUncertainties(bool LogNormal = true) {
      fUseLog = LogNormal;
   }
   //! \brief get boolean for relative uncertainties (LogNormal)
   bool GetLogNormalUncertainties() const {
      return fUseLog;
   }

   // ------------------------------------------- //
   // --- use nuisance parameter or covariance matrices
   // ------------------------------------------- //
   //!   \brief Consider correlated systematic uncertainties as 'shifts' with nuisance
   //!   \brief parameter or as a covariance matrix. 
   //!   default: true (use nuisance parameters). 
   //!   Note: this setting does *not* apply to already present uncertainties, but
   //!   only to correlated uncertainties set in the following.
   void UseNuisanceParameters(bool UseNuisance=true) {
      fUseNuisance = UseNuisance;
   }


   // ------------------------------------------- //
   // --- rescale (absolute) uncertainties
   // ------------------------------------------- //
   //!  \brief Rescale uncertainties
   //!  All uncertainties are stored interally as absolute uncertainties. We assume that
   //!  these are proportional to the data. In oreder to avoid (or reduce) a bias, one
   //!  can re-scale all uncertainties, e.g. with the predictions (from a previous fit, or
   //!  from a selected template.) The same function is availalbe for LTF::LiTeFit
   //!  and a previous LiTeFit can be repeated with rescaled uncertainties.
   //!
   //!  Note: this option cannot be applied when using relative uncertainties in the chi2.
   //!
   //!  Note: all errors passed to LTF are expected to scale with data, or as relative uncertainties.
   //!  This method then 're-scales' the size of the absolute uncertainties.
   //!
   //!  Note: template and response-matrix uncertainties are omitted.
   //!  Note: also statistical and uncorrelated uncertainties are rescaled!
   void RescaleInputErrors(std::vector<double> values ) {
      if ( values.size() != size_t(fDt.rows()) ) {
         std::cerr<<"Error in LTF::RescaleInputErrors. Size of input vector not compatible with size of data vector."<<std::endl;
         exit(1);
      }
      fErrorScale = Eigen::Map<Eigen::VectorXd>(&values[0],values.size());
   }


// ------------------------------------------- //
// --- members
protected:
   Eigen::VectorXd fDt;                                                       //!< data distribution
   std::vector<std::pair<vector<double>,Eigen::VectorXd> > fTmplDistN;        //!< templates
   Eigen::MatrixXd                                       fA;                  //!< the response matrix (opional)
   std::vector<std::pair<std::string,Eigen::MatrixXd > > fVs;                 //!< Covariance matrices
   std::vector<std::pair<std::string,Eigen::VectorXd > > fSys;                //!< correlated (syst.) uncertainty
   std::map<std::string,LTF::Uncertainty>                fSysIsConstr;        //!< uncertainty specification for fSys
   std::map<std::string,Eigen::MatrixXd >                fVsExt;              //!< Covariance matrices, treated as external
   std::map<std::string,Eigen::VectorXd >                fSysExt;             //!< correlated uncertainties, treated as external
   std::map<std::string,std::map<vector<double>,Eigen::VectorXd> > fSysY;     //!< uncertainty of the templates (uncorrelated)
   std::map<std::string,Eigen::MatrixXd>                 fSysA;               //!< uncertainty of the response matrix (uncorrelated)
   std::map<std::string,double>                          fcorrSys;            //!< Correlation coefficient of the systematic uncertainties (dA,dY)
                                                                           
   Eigen::VectorXd fErrorScale;                                               //!< Reference values for error rescaling (default: data)
   std::vector<double> fGamma;                                                //!< exponential parameter (a+b*x^\gamma)
   bool   fUseLog = false;                                                    //!< Use relative uncertainties
   bool fUseNuisance = true;                                                  //!< Nuisance parameter or a covariance matrices
   LiTeFit fLTF;                                                              //!< result
                                                                              
   int fBinMin = 0;                                                           //!< Fit range
   int fBinNN  = 0;                                                           //!< Fit range

// ------------------------------------------- //
// --- protected methods
protected:
   size_t N() const;                                                          //!< \brief return number of points. This value is not valid after applying apply Fit Range
   Eigen::MatrixXd Calc_M() const;                                            //!< \brief Calculate matrix M from template reference points
   Eigen::MatrixXd Calc_Y() const;                                            //!< \brief Calculate matrix Y, representing the templates
   Eigen::MatrixXd Calc_X() const;                                            //!< \brief Calculate matrix X, representing the templates
   std::map<std::string,Eigen::MatrixXd > Calc_dY(
      std::map<std::string,std::map<vector<double>,Eigen::VectorXd> >& fSysY) const; //!< \brief Calculate template uncertainties dY
   void SetLiTeFitInput();                                                    //!< \brief Set the input to LiTeFit
   void SetCorrSys(const std::string& name, double c);                        //!< \brief set correlation coefficient for syst. uncertainties
   

// ------------------------------------------- //
// ---- static helper functions
public:

   //! \brief calculate uncertainty from a covariance matrix
   static Eigen::VectorXd GetDelta_V(const Eigen::MatrixXd& V) { 
      return V.diagonal().cwiseSqrt();
   }

   //! \brief calculate covariance matrix from an uncertainty
   static Eigen::MatrixXd GetV_Delta(const Eigen::VectorXd& d, double corr);

   //! \brief calculate covariance matrices from (systematic) shifts
   template< class T >
   static std::map<std::string,Eigen::MatrixXd> GetV_DeltaSys(const T& DeltaSys) {     
      std::map<std::string,Eigen::MatrixXd> ret;
      for ( auto& s : DeltaSys ) ret[s.first] = LTF::GetV_Delta(s.second,1);
      return ret; 
   }

   //! \brief calculate size of uncertainty from covariance matrix
   template< class T >
   static std::map<std::string,Eigen::VectorXd> GetDelta_Vsource(const T& Vsource) {  
      std::map<std::string,Eigen::VectorXd> ret;
      for ( auto& s : Vsource ) ret[s.first] = LTF::GetDelta_V(s.second);
      return ret; 
   }

   //! \brief Calculate covariance matrix as sum of individual matrices and systematic shifts
   //! \param Vs        Covariance matrices, i.e. uncorr. and/or cov. uncertainties
   //! \param DeltaSys  correlated systematic uncertainty
   static Eigen::MatrixXd VSum( const std::map<std::string,Eigen::MatrixXd>& Vs, const std::map<std::string,Eigen::VectorXd>& DeltaSys );

   //! \brief Calculate covariance matrix as sum of individual matrices and systematic shifts
   //! \param Vs        Covariance matrices, i.e. uncorr. and/or cov. uncertainties
   //! \param Sys       correlated systematic uncertainties
   static Eigen::MatrixXd VSum( const std::vector<std::pair<std::string,Eigen::MatrixXd > >& Vs, const std::vector<std::pair<std::string,Eigen::VectorXd > >& Sys = std::vector<std::pair<std::string,Eigen::VectorXd > >() );

   //! \brief Calculate correlation matrix from given covariacne matrix V
   static Eigen::MatrixXd Cov_to_Cor(const Eigen::MatrixXd& V);

   //! \brief Calculate an Eigen-matrix from std::vector<std::vector<>>. T=double
   static Eigen::MatrixXd Std_to_Mat(const std::vector<std::vector<double> >& A);
   
};

#endif
