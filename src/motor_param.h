# ifndef MOTOR_PARAM_H_
# define MOTOR_PARAM_H_

#include <RcppArmadillo.h>
#include "motor_ctrl.h"
#include "motor_reg.h"

class param {
public:
  // Use default copy constructor, but create virtual Clone function for cloning derived classes correctly
  param (arma::vec position_, ctrl c_, reg* preg_, bool fixed_, bool screen_, double lambda_factor_, unsigned int nrow_);

  // Full position has length p_full, the rest have length p
  arma::vec position_full;
  arma::uvec active;
  arma::vec position, differential, gnhda;  // Position, differential at that and gauss newton hessian diag approx

  // Sizes
  unsigned int p_full, p;

  bool fixed, screen;   // Should this parameter remain fixed? Should it be screened? Has it converged?
  ctrl c;
  reg* preg;
  double lambda_factor;       // Multiplied on lambda
  double pen;
  void set_lambda(double lambda_);  // Gets multiplied


  // For using matrix-type semantics (i.e., if one wants to )
  unsigned int nrow;    // (if not 0, then position can be cast as matrix with nrow rows by calling get_position)
  arma::uvec coord;     // The coordinates of the actives, assuming all 0-columns have been removed
  void set_coord(void);
  unsigned int r_sub;   // Number of columns with non zeros
  arma::mat get_position(void); // returns position as matrix with nrow rows
  arma::mat get_dparam(const arma::vec v);  // Returns differential of get_position() * v wrt to active param

  // Constructor
  param (Rcpp::List param_);
  param (Rcpp::List reg_struct, arma::vec x, unsigned int nrow_);
  ~param();

  // Run proximal operator in given direction
  arma::vec proximal(arma::vec position_, arma::vec direction);

  // Steins degrees of freedom
  double dfs(arma::mat X);

  // Virtual restore/reduce for the system (reset_active(void) restores full system)
  void set_active(arma::uvec active_);
  void reset_active(void);
};

arma::uvec get_span(const arma::vec &v);



# endif
