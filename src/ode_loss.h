# ifndef ODE_LOSS_H_
# define ODE_LOSS_H_

#include <RcppArmadillo.h>
#include "ode_solver.h"
#include "motor_reg.h"
#include "motor_fun.h"


/* ode_loss class.
 *  Defines the loss function to be passed to optim-class
 */
class ode_loss : public fun {
private:
  arma::mat w_nan_guard (arma::mat w_, arma::mat y_);
  arma::mat y_nan_guard (arma::mat y_);

public:
  arma::vec time;
  unsigned int n_obs;
  unsigned int n_series;

protected:
  std::vector<sc>* pvsc_full;

public:
  std::vector<sc> vsc;  // Will have length one less than vparam
  likeli* plikeli;

  arma::mat y;  // Data (samples stored in columns),
  bool w_exist;
  arma::mat w;  // weights in loss function, same format as data

  field* pfield;      // ODE should be public, so that parameters at will and one can reduce the system (from within optim class)

  arma::uvec exact_diff;  // Logicals, length as vparam, indicating if differentials are exact

private:
  std::vector<arma::vec> unscaled_pos;    // Has length one less than vparam (x0 not in it)

protected:
  void store_positions(void);                   // Store positions from vparam to unscaled_positions
  void scale_positions(unsigned int i_series);  // Replace positions from vparam with scaled versions of unscaled_positions
  void restore_positions(void);                 // Replace positions in vparam with unscaled_positions

public:
  ode_loss(Rcpp::List ode_struct, arma::mat &A, arma::mat &B, arma::mat x0, Rcpp::Nullable<Rcpp::List> param_, Rcpp::List opt_struct, std::vector<sc>* pvsc_full_, likeli* plikeli_);
  ode_loss(const ode_loss& other);
  ode_loss& operator=(const ode_loss& rhs);

  // Virtual clone and destructor (deep copying used in create_ploss, else no deep copying)
  // virtual loss* Clone() = 0;      // is supposed to make deep clone of psolver, but not pode and plikeli
  // virtual loss* deep_Clone() = 0;
  virtual ~ode_loss() {};

  // All virtuals from fun carry down, except they must can call
  void set_active_base(std::vector<param>::iterator param_it);
  /* Reason for splitting up in reduce and reduce_base is as follows:
   * 'im'-subclass furthermore has X, which needs to be reduced as well,
   * but both subclasses should call reduce_base to reduce the elements they have in common,
   * i.e, those in 'loss'
   */

  // Evaluate Steins Degrees of Freedom over a vector of lambda, calls the virtual function below
  arma::mat get_dfs(arma::vec lambdas);

  // Virtual function which returns dhat(y)/dparam (n-x-p) matrix
  virtual arma::mat get_X4dfs(std::vector<param>::iterator param_it) = 0;
};


// exact loss subclass (for exact loss optimisation, i.e., using numerical solver)
class exact : public ode_loss {
public:
  solver* psolver;    // Public so that conv_code can be extracted
  exact(Rcpp::List ode_struct, arma::mat &A, arma::mat &B, arma::mat x0, Rcpp::Nullable<Rcpp::List> param_, Rcpp::List opt_struct, std::vector<sc>* pvsc_full_, likeli* plikeli_);
  exact(const exact& other);
  exact& operator=(const exact& rhs);
  // exact* Clone();
  // exact* deep_Clone();
  ~exact() {delete psolver;};

  void evaluate(void);
  void evaluate_differential(std::vector<param>::iterator param_it);
  void set_active(std::vector<param>::iterator param_it);

  arma::mat get_X4dfs(std::vector<param>::iterator param_it);
};


// matching ode_loss subclass (does matching (integral or gradient))
class matching : public ode_loss {
protected:
  arma::mat x;    // Smoothed trajectory (without time)
  arma::vec ti;   // input time (from x)
  arma::uvec ti_jumps; // indeces where ti jumps
  arma::vec to;   // output time (from 'time')
  arma::uvec to_jumps; // indeces where to jumps

  arma::vec Y;    // matching dependent, is first filled out in subclass constructor
  arma::vec W;    // matching dependent, is first filled out in subclass constructor

  arma::vec X0;   // Initial states interpolated from x, stacked on each other
  arma::mat xout; // Interpolated values of x at 'time'

  // Functions for interpolating
  arma::vec interpolate(const double t1, const arma::vec x1, const double t2, const arma::vec x2, const double tm);
  arma::mat interpolate(const arma::vec tin, const arma::vec tout, const arma::mat xin);

  // Constructor
  matching(Rcpp::List ode_struct, arma::mat &A, arma::mat &B, arma::mat x0, Rcpp::Nullable<Rcpp::List> param_, Rcpp::List opt_struct, std::vector<sc>* pvsc_full_, likeli* plikeli_, arma::mat x_);

public:
  arma::vec match;                    // What we are supposed to match with Y (i.e., f(x(s), param) or integral of that)
  std::vector<arma::mat> dmatch;      // Has length 1 smaller than vparam, x0 is ignored
  std::vector<arma::mat> dmatch_full; // Differentials wrt param, _full only used if system is linear

  // Sets match and dmatch based on system linear or not
  void set_match(void);
  void set_dmatch(std::vector<param>::iterator param_it);

  /* The non linear methods called in case of non linear system, depends on matching method
   * set_nl_match is supposed to set match to whatever Y is matched against
   * set_nl_dmatch, however must change the matrix pointed to by dmatch_it to the differential of above
   *
   * The reason for this dmatch pointer is that the linear case will set dmatch_full in constructor,
   * while in non-linear case only dmatch will be used (and set repeatedly).
   */
  virtual void set_nl_match(void) = 0;
  virtual void set_nl_dmatch(std::vector<arma::mat>::iterator dmatch_it, std::vector<param>::iterator param_it) = 0;

  // The virtual fun functions, that basically calls set_match and plikeli
  void evaluate(void);
  void evaluate_differential(std::vector<param>::iterator param_it);
  void set_active(std::vector<param>::iterator param_it);

  arma::mat get_X4dfs(std::vector<param>::iterator param_it);

  // Get the designs and stuff
  Rcpp::List get_XYW(void);
};


// Integral matching
class integral : public matching {
  // Functions for doing integration
  void simpson(arma::vec &x1, arma::mat &y1, arma::vec x2, arma::mat &y2, std::vector<param>::iterator param_it);
  void simpson(const double t1, arma::vec &x1, arma::mat &y1, const double t2, arma::vec x2, arma::mat &y2, arma::vec tm, arma::mat xm, std::vector<param>::iterator param_it);
  arma::mat simpson(arma::vec tin, arma::vec tout, arma::mat xin, arma::mat xout, std::vector<param>::iterator param_it);

public:
  // Constructor (here you must set Y, W)
  integral(Rcpp::List ode_struct, arma::mat &A, arma::mat &B, arma::mat x0, Rcpp::Nullable<Rcpp::List> param_, Rcpp::List opt_struct, std::vector<sc>* pvsc_full_, likeli* plikeli_, arma::mat x_);

  // How to set non linear match and dmatch
  void set_nl_match(void);
  void set_nl_dmatch(std::vector<arma::mat>::iterator dmatch_it, std::vector<param>::iterator param_it);

};


# endif
