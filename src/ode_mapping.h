# ifndef ODE_MAPPING_H_
# define ODE_MAPPING_H_

#include <RcppArmadillo.h>

class mapping {
public:
  unsigned int d, nrow, ncol;

private:
  // The mapping function: R^d -> R^(nrow-x-ncol)
  virtual arma::mat map(arma::vec x) = 0;



  /* ------------------------------
   * Interpolation helper function
   * ------------------------------*/

  // Interpolation of two points x1 at time t1 and x2 at time t2, target time tm
  arma::vec interpolate(const double t1, const arma::vec x1, const double t2, const arma::vec x2, const double tm);



  /* ------------------------------------
   * Simpson integration helper functions
   * ------------------------------------*/

  // int(t1, t2, map(x(s))) claimed by in y1 * (t2 - t1)
  // Convention: y2 becomes the new y1, x1 is the old x2
  void simpson(arma::vec &x1, arma::mat &y1, arma::vec x2, arma::mat &y2);

  // int(t1, t2, map(x(s))) stored in y1 based on the intermediates xm at tm. If no intermediates, it just calls above simpson.
  void simpson(const double t1, arma::vec &x1, arma::mat &y1, const double t2, arma::vec x2, arma::mat &y2, arma::vec tm, arma::mat xm);


public:
  mapping(unsigned int d_, unsigned int nrow_, unsigned int ncol_);

  // Interpolation of time points xin (stored column wise) at time tin, to time points tout.
  // Assumes tin and tout are non-decreasing!! and length(tin) and ncol(xin) matches
  // returns interpolated value of x at tout, stored in columns (uses linear interpolate, unless oob)
  arma::mat interpolate(const arma::vec tin, const arma::vec tout, const arma::mat xin);

  // Simpson integration
  // returns {int(tout(i), tout(i + 1), map(x(s)))}_i for i = 0,...,end-1, stacked on top of each other
  // length(tin) and ncol(xin) matches and length(tout) and ncol(xout) matches
  arma::mat simpson(arma::vec tin, arma::vec tout, arma::mat xin, arma::mat xout);

};



class power : public mapping {
  arma::mat A;

public:
  power(arma::mat A_);
  arma::mat map(arma::vec x);
};



class fracpower : public mapping {
  arma::mat A, B;

public:
  fracpower(arma::mat A_, arma::mat B_);
  arma::mat map(arma::vec x);
};



# endif
