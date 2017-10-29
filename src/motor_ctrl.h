# ifndef MOTOR_CTRL_H_
# define MOTOR_CTRL_H_

#include <RcppArmadillo.h>

struct ctrl {
  // Number of steps
  unsigned int step, step_cycle, step_max, step_screen, backtrack_max;

  // Step length
  double tau, tau_init, tau_min, tau_scale;

  // Tolerance
  double tol;

  // Convergence codes
  unsigned int backtrack_code, diff_code;

  ctrl (unsigned int step_cycle_, unsigned int step_max_, unsigned int step_screen_, unsigned int backtrack_max_, double tau_init_, double tau_min_, double tau_scale_, double tol_);
  ctrl (void);

  // Reseting stuff
  void reset_step(void);
  void reset_tau(void);
  void reset_code(void);
};



# endif
