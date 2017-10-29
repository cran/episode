//=============================================
/* motor_ctrl.cpp
*
* Content:
* - ctrl definitions
*/
//=============================================

#include "motor_ctrl.h"



//=================================
// ctrl struct
//=================================

// Constructor (no copy constructor or destructor defined, as param contains no raw pointers)
ctrl::ctrl (unsigned int step_cycle_, unsigned int step_max_, unsigned int step_screen_, unsigned int backtrack_max_, double tau_init_, double tau_min_, double tau_scale_, double tol_) :
  step(0),                        // Step counter
  step_cycle(step_cycle_),        // Maximal number of steps per cycle
  step_max(step_max_),            // Maximal number of steps in total
  step_screen(step_screen_),      // Maximal number of steps between screenings
  backtrack_max(backtrack_max_),  // Maximal number of step length reductions per step
  tau(tau_init_),                 // Step length
  tau_init(tau_init_),            // Initial step length
  tau_min(tau_min_),              // Minimal step length
  tau_scale(tau_scale_),          // Step length scale-downer
  tol(tol_),                      // Minimal distinguishable size
  backtrack_code(0),
  diff_code(0) {
    // if (step_cycle < 0 || step_max < 0 || backtrack_max < 0) Rcpp::stop("step_max, step_cycle or backtrack_max is negative.");
    if (tau_init <= tau_min) Rcpp::stop("tau_init is smaller than tau_min.");
    if (tau_min <= 0) Rcpp::stop("tau_min is non-positive.");
    if (tau_scale >= 1 || tau_scale <= 0) Rcpp::stop("tau_scale is not in (0, 1).");
    if (tol < 0) Rcpp::stop("tol is negative.");
  }

// Default constructor
ctrl::ctrl (void) :
  step(0),            // Step counter
  step_cycle(10),     // Maximal number of steps per cycle
  step_max(100),      // Maximal number of steps in total
  step_screen(10),    // Maximal number of steps between screenings
  backtrack_max(25),  // Maximal number of step length reductions per step
  tau(.1),            // Step length
  tau_init(.1),       // Initial step length
  tau_min(1e-12),     // Minimal step length
  tau_scale(.5),      // Step length scale-downer
  tol(1e-6),          // Minimal distinguishable size
  backtrack_code(0),
  diff_code(0) {
  // if (step_cycle < 0 || step_max < 0 || backtrack_max < 0) Rcpp::stop("step_max, step_cycle or backtrack_max is negative.");
  if (tau_init <= tau_min) Rcpp::stop("tau_init is smaller than tau_min.");
  if (tau_min <= 0) Rcpp::stop("tau_min is non-positive.");
  if (tau_scale >= 1 || tau_scale <= 0) Rcpp::stop("tau_scale is not in (0, 1).");
  if (tol < 0) Rcpp::stop("tol is negative.");
}

void ctrl::reset_step(void) {
  step = 0;
}

void ctrl::reset_tau(void) {
  tau = tau_init;
}

void ctrl::reset_code(void) {
  backtrack_code = 0;
  diff_code = 0;
}
