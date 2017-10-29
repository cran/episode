# ifndef ODE_SOLVER_H_
# define ODE_SOLVER_H_

#include <RcppArmadillo.h>
#include "ode_field.h"
#include "ode_butcher.h"



// State of system, (score indicates if differential ought to be updated as well and exact if it should be exact)
struct state {
  double t;
  arma::vec x;
  bool score;
  bool exact;
  arma::mat x_dparam;    // Differential wrt parameter (one at the time)

  field* pfield;
  std::vector<param>::iterator param_it;

  // Constructor, the first is useful for state_proposed, the second for starting a sensitivity solver
  state (double t_, arma::vec x_, field* pfield_ = NULL, bool score_ = false, bool exact_ = false, arma::mat x_dparam_ = arma::mat());
  state (double t_, arma::vec x_, std::vector<param>::iterator param_it_, field* pfield_ = NULL, bool exact_ = false);

  // Tells if param_it points to x0 parameter
  bool is_x0(void);

  // Reset sensitivity, which checks what iterator to use
  void reset_dparam(void);
};


/* Abstract solver base class: requires two functions
 *  - proposal (takes the current state and proposes a new state moved into the future, given param and field-pointer)
 *  - reset (resets local control parameters)
 */
class solver {
  /* Virtual proposal function which (for a given solver) takes position (t, x, x_dx0, etc) and proposes a new state of
   * the field system (moved a bit ahead). It also updates the number of steps taken.
   */
  virtual state proposal (state &state_, unsigned int &step) = 0;

protected:
  const unsigned int step_max;

public:
  unsigned int code_conv;

  solver (unsigned int step_max_ = 100);

  // Use default copy constructor, but create virtual Clone function for cloning derived classes correctly
  virtual solver* shallow_Clone() = 0;
  virtual solver* deep_Clone() = 0;
  virtual ~solver() {}

  void refresh_code (void);
  virtual void reset (void) = 0;                    // Reset local and adaptive control parameters

  void step (state &state_, double t_end);          // Move state to t_end
  arma::mat solve(field* pfield, arma::vec time);   // Solve the system at given timepoints
                                                    // Solve the system and sensitivity at given timepoints
  arma::mat solve(field* pfield, std::vector<param>::iterator param_it, arma::vec time, arma::cube &sens_dparam, bool approx);
  arma::mat sensitivity_param(field* pfield, std::vector<param>::iterator param_it, arma::vec time, std::vector<sc>::iterator sc_it, bool approx);
};


// Subclass, Runge-Kutta EMbedded Pair
class rkemp : public solver {
  const butcher b;
  const double Tol, h_init;
  double h;

  double ratio(arma::vec w, arma::vec z);
  state proposal(state &state_, unsigned int &step);

public:
  rkemp(std::string str_ = "rkf45", unsigned int step_max_ = 100, double Tol_ = 1e-6, double h_init_ = 1e-4);
  rkemp* shallow_Clone();
  rkemp* deep_Clone();
  void reset(void);
};


// Misc
solver* create_psolver(Rcpp::List ode_struct);

# endif
