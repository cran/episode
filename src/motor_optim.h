# ifndef MOTOR_OPTIM_H_
# define MOTOR_OPTIM_H_

#include <RcppArmadillo.h>
#include "motor_fun.h"


class optim {
  // If not for test_optim, everything except: pfun, constructor and optimise() should be private
protected:
  static unsigned int s_id;

public:
  unsigned int id;
  arma::vec direction;  // Current direction in which to go for parameter in focus
  fun* pfun;
  bool trace;
  double tol;

  // Type of direction (e.g., gradient descent (gd))
  const std::string dir_type;

  // Constructor
  optim (fun* pfun_, bool trace_, double tol_, std::string dir_type_);
  optim (const optim& other);
  optim& operator=(const optim& rhs);

  // Reset codes, step and tau.
  void reset_optim_param(void);

  // Get the direction for specific parameter
  void get_direction(std::vector<param>::iterator it);

  // get_and_set_active(it) uses it->preg to get active, then sets them in param and in fun
  void get_and_set_active(std::vector<param>::iterator it);

  // Backtracking that backtracks in 'direction'
  void backtrack(std::vector<param>::iterator it);

  // Evaluate convergence criteria by comparing to old loss and position
  bool conv_criteria(std::vector<param>::iterator it, double l_old, arma::vec position_old);

  // Runs at most step_cycle steps, checks param's screen and fixed for permission. Returns convergence flag.
  bool optimise(std::vector<param>::iterator it);

public:
  // Runs through parameters and calls the above on each of them
  void optimise(void);

  // Set lambda
  void set_lambda(double lambda_);
};



//=================================
// Apply to vectors of optim
//=================================

// move: moves all particles (runs optimise or screen depending on input)
void move(std::vector<optim> &voptim);

// set_lambda: sets new lambda across optims
void set_lambda(std::vector<optim> &voptim, double lambda_);

# endif
