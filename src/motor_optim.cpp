//=============================================
/* motor_optim.cpp
 *
 * Content:
 * - optim and subclasses definitions
 * - functions for applying over voptim
 *
 */
//=============================================

#include "motor_optim.h"




//=================================
// optim class
//=================================

unsigned int optim::s_id;

// Constructor
optim::optim(fun* pfun_, bool trace_, double tol_ = 1e-8, std::string dir_type_ = "gd") :
  id(++s_id),
  direction(arma::zeros<arma::vec>(pfun_->vparam.begin()->differential.n_elem)),
  pfun(pfun_),
  trace(trace_),
  tol(tol_),
  dir_type(dir_type_) {

}

// Copy constructor
optim::optim(const optim& other) :
  id(other.id),
  direction(other.direction),
  pfun(other.pfun),
  trace(other.trace),
  dir_type(other.dir_type) {

}

// Assignment
optim& optim::operator=(const optim& rhs) {
  if (&rhs != this) {
    id = rhs.id;
    direction = rhs.direction;  // Current direction in which to go for parameter in focus
    pfun = rhs.pfun;
    trace = rhs.trace;

    (std::string &) dir_type = rhs.dir_type;
  }
  return *this;
}

// Reset codes, step and tau
void optim::reset_optim_param(void) {
  for (std::vector<param>::iterator param_it = pfun->vparam.begin(); param_it < pfun->vparam.end(); param_it++) {
    /* tau is clamped in between t1 = 'log-midpoint of tau_min and tau_init' and t2,
     * such log distance from t1 to tau_init is equal to log distance from t2 to tau_init
     */
    double t1 = std::exp((std::log(param_it->c.tau_init) + std::log(param_it->c.tau_min)) / 2.0);
    double t2 = std::exp(2.0 * std::log(param_it->c.tau_init) - std::log(t1));
    param_it->c.tau = std::max(t1, std::min(param_it->c.tau, t2));

    param_it->c.reset_step();
    param_it->c.reset_code();
  }
}

// Get direction in which to go
void optim::get_direction(std::vector<param>::iterator it) {
  // Evaluate fun
  pfun->evaluate(it, true);

  // Set direction
  if (dir_type == "gd") {
    direction = -it->differential;
  } else {
    Rcpp::stop("error in optim -> 'dir_type' not recognised.");
  }

  // diff_code mill
  // Prevent massive directions (if (position - prox(position, direction, tau_min)) / tol > 0.1, then no amount of bactracking can assure right tolerance)
  arma::vec minimal_change = it->position - it->preg->prox(it->position, direction, it->c.tau_min);

  if (it->c.diff_code == 2) {
    // If diff_code is 2, then it was a signal from backtracking that the action in the else-part below did not work, we remove large entries
    if (arma::norm(minimal_change, "inf") > 0.1 * it->c.tol) {
      direction.elem(arma::find(arma::abs(minimal_change) > 0.1 * it->c.tol)).zeros();
    }
  } else {
    // If we where not signaled from backtrack, but direction is large, we scale it instead
    if (arma::norm(minimal_change, "inf") > 0.1 * it->c.tol) {
      it->c.diff_code = 1;
      direction *= it->c.tol / (10 * arma::norm(minimal_change, "inf"));
    } else {
      it->c.diff_code = 0;
    }
  }

  // Note that the actions above are never mixed (which would be bad), since direction is set first.
}

// Finds and sets the active set,
void optim::get_and_set_active(std::vector<param>::iterator it) {
  // unsigned int nparam = it - pfun->vparam.begin();
  // Rcpp::Rcout << "get_and_set_active(" << nparam << ")" << std::endl;

  // Store position in position_full and update direction
  it->reset_active();
  pfun->set_active(it);
  get_direction(it);

  // Get the active
  arma::uvec active = it->preg->get_active(it->position_full, direction, it->c.tau_min);

  // Set active in param and direction
  it->set_active(active);
  direction = direction.elem(active);
  pfun->set_active(it);
}

// Backtracks in 'direction' (ASSUMES pfun->value is up to date with 'it')
void optim::backtrack(std::vector<param>::iterator it) {
  unsigned int backtrack_n = 0;           // Number of backtrack steps
  it->c.backtrack_code = 1;               // Initially 1 and set to zero if succesful

  double best_tau  = it->c.tau;
  arma::vec best_position = it->position;       // 'best's represent current best
  double best_value = pfun->value;
  const arma::vec old_position = it->position;  // 'old' represents fix current position

  bool went_up = false;   // Temporary flag to indicate if dist(best_pos, pos) < tol, to set up

  // Backtrack loop: if reaches end unsuccesfully, then it->position and l are not changed.
  while (backtrack_n < it->c.backtrack_max && it->c.tau_min < it->c.tau) {
    // Set position according to proximal operator on direction
    it->position = it->proximal(old_position, direction);

    // Evaluate loss + penalty (stored in pfun->value)
    pfun->evaluate(it, false);

    // // Rcout for seeing how it goes
    // Rcpp::Rcout.precision(12);
    // Rcpp::Rcout << "Current best:" << std::endl <<
    //   " tau \t\t" << best_tau << std::endl <<
    //     " value \t\t" << best_value <<  std::endl <<
    //       " position \t" << best_position.t() << std::endl <<
    //         "Current:" << std::endl <<
    //           " tau \t\t" << it->c.tau << std::endl <<
    //             " value \t\t" << pfun->value <<  std::endl <<
    //               " position \t" << it->position.t() << std::endl <<
    //                 " direction \t" << direction.t() << std::endl;
    // Rcpp::Rcout.precision(6);

    if (pfun->value < best_value) {
      /* Succesfull! Find new position and tau:
      *
      * First we check if this is the first backtrack step,
      *  if so
      *    We increase tau untill no sharp improvement (this will prevent a bad step from before from ruining the next steps)
      *  if not
      *    We check if a tau of in between current and previous would provide better improvement,
      *    then we check if a small tau would. In the end we choose which ever of the
      *    three choices has lowest loss.
      */

      // Store current best
      best_tau = it->c.tau;
      best_position = it->position;
      best_value = pfun->value;

      if (backtrack_n == 0) {
        // If no previous decrease of tau took place we increase tau untill we get no improvement
        it->c.tau /= it->c.tau_scale;
        it->position = it->proximal(old_position, direction);
        pfun->evaluate(it, false);

        while (backtrack_n < it->c.backtrack_max && pfun->value < best_value) {
          backtrack_n ++;

          // Save current best
          best_tau = it->c.tau;
          best_position = it->position;
          best_value = pfun->value;

          // Prepare next candidate
          it->c.tau /= it->c.tau_scale;
          it->position = it->proximal(old_position, direction);
          pfun->evaluate(it, false);
        }
        /* When it breaks out, the current best position, value and tau
         * are stored in old_position, old_value and old_tau.
         */

      } else {
        // If we did decrease tau before, then evaluate three alternatives

        // First alternative: higher value of tau (between current tau and tau / tau_scale)
        it->c.tau = best_tau * (1 + 1 / it->c.tau_scale) / 2;
        it->position = it->proximal(old_position, direction);
        pfun->evaluate(it, false);

        if (pfun->value < best_value) {
          // New current best
          best_tau = it->c.tau;
          best_position = it->position;
          best_value = pfun->value;
        }

        /* Second alternative: some lower values of tau
         * If first alternative failed the first is just next tau in backtrck, else it is halfway there
         * If it improves, we continue for some more values, else we stop
         */
        unsigned int some_lower = 0;
        it->c.tau = best_tau;
        while (some_lower < 15) {
          some_lower++;

          it->c.tau *= it->c.tau_scale;
          it->position = it->proximal(old_position, direction);
          pfun->evaluate(it, false);

          if (pfun->value < best_value) {
            // New current best
            best_tau = it->c.tau;
            best_position = it->position;
            best_value = pfun->value;

            // Continue, unless we have checked more than 15 back steps
          } else {
            // If we did not improve, then no need to continue further
            break;
          }
        }

        /* If both alternatives fail, then it falls back to first choice, as
         * old_value and old_position would never have been overwritten.
         */
      }

      // If we made it here, then backtrack code should be 0 and tau should be current best
      it->c.backtrack_code = 0;
      it->c.tau = best_tau;
      break;

    } else if (arma::norm(best_position - it->position, "inf") < it->c.tol) {
      /* The current best and the position produced by proximal are up to tolerance identical
       * All following step lengths will have the same property, so we consider it finished.
       *
       * This happens typically because the differential is numerically unstable or
       * the bounds prevents further improvements.
       */
      if (!went_up) {
        /* If we made it here, then the backtracking caused it to too small,
         * give it one more change by "resetting" tau
         */
        double proposed_tau = std::exp(std::log(best_tau) - std::log(it->c.tau_scale) * 10.0); // Go ten steps up from best_tau
        double max_tau = std::exp(std::log(it->c.tau_init) - std::log(it->c.tau_scale) * 15.0);
        it->c.tau = std::max(it->c.tau_init, std::min(proposed_tau, max_tau));
        went_up = true;
      } else {
        // We tried, so it probably is local minima
        it->c.backtrack_code = 0;
        it->c.tau = best_tau;
        break;
      }
    } else {
      // No improvement. Update counter and tau
      backtrack_n ++;
      it->c.tau *= it->c.tau_scale;
    }
  }

  /* Update position and value to current best
   *
   * tau already has been updated, since if it is small due to
   * lack of fit in backtrack, it should carry over.
   * If we said it->c.tau = old_tau, then if no improvement was ever found
   * it would not manifest itself.
   */
  it->position  = best_position;  // If no improvement encountered, then best_position == old_position
  pfun->value   = best_value;

  // Prevent tau from being < tau_min to begin with in next step
  it->c.tau = std::max(it->c.tau, it->c.tau_min / it->c.tau_scale);
}

// Convergence criterion comparing current values stored in it to old ones
bool optim::conv_criteria(std::vector<param>::iterator it, double value_old, arma::vec position_old) {
  using namespace std;
  // old one arma::norm(param - param_old, 2) + arma::norm(arma::vectorise(x0 - x0_old), 2) < tol;
  // denum counts the number of non-neglibale units in each parameter (tol_param is smallest acceptable unit)

  // Check for large directions
  double denum = abs(it->c.tau * arma::norm(direction, 2) / it->c.tol);
  arma::uvec fdirection = arma::find(direction);
  double dim = fdirection.n_elem;
  dim = std::max(dim, 1.0);

  if (denum > 1e3 * dim) {
    // Large gradient, continue search (this is signaled to get_direction through optimise)
    return false;
  } else {
    // Small gradient, then check actual improvement
    value_old = abs(value_old - pfun->value);
    position_old = arma::abs(it->position - position_old);
    denum = arma::norm(position_old, 2) / it->c.tol;
    fdirection = arma::find(position_old > it->c.tol * 0.1);
    dim = fdirection.n_elem;
    dim = std::max(dim, 1.0);

    if (denum > 1e2 * dim) {
      // For large step lengths, always continue
      return false;
    } else {
      // Else, it depends on the relative improvement in loss (note if value_old = pfun->value we return false, so that potential large gradients are taken care of)
      return (value_old <= 0 || denum <= 0) ? false : log(value_old) - log(denum) < log(tol);
    }
  }
}

// returns bool which says if converged or not
bool optim::optimise(std::vector<param>::iterator it) {
  // Flag to continue optimising
  bool cont = true;

  // If already reached maximum, set to false
  if (it->c.step >= it->c.step_max) {
    cont = false;
  }

  // Return value: has it converged,
  bool converged = false;

  if (!(it->fixed) && cont) {
    // bool bb = it->screen && (it->c.step % it->c.step_screen == 0);
    // Rcpp::Rcout << " optim::optimise -> Get direction. Screen: " << bb << std::endl;

    if (it->screen && (it->c.step % it->c.step_screen == 0)) {
      // If we screen we set new active set (and get direction + pfun->value)
      get_and_set_active(it);
    } else {
      // If not, we only acquire the direction and pfun->value
      get_direction(it);
    }
    // Rcpp::Rcout << " optim::optimise -> Got direction" << std::endl;

    // For saving last step
    arma::vec position_old = it->position;
    double value_old = pfun->value;

    // if (trace) {
    //   Rcpp::Rcout.precision(12);
    //   Rcpp::Rcout << "\t Current loss: " << std::fixed << pfun->value << std::endl;
    //   Rcpp::Rcout.precision(6);
    // }

    // Re-initialise step_n (gets added up in screen)
    unsigned int step_N = std::min(it->c.step + it->c.step_cycle, it->c.step_max);
    unsigned int step_0 = it->c.step;
    while (it->c.step < step_N) {
      value_old = pfun->value;

      // Get direction and backtrack
      if (it->c.step != step_0) {
        // For first step direction was evaluated in if(it->screen)-part

        if (it->screen && (it->c.step % it->c.step_screen == 0)) {
          // If we screen we set new active set (and get direction + pfun->value)
          get_and_set_active(it);

          // If new active set chosen, we should update old_position so converge_crit works
          position_old = it->position;
        } else {
          // If not, we only acquire the direction and pfun->value
          get_direction(it);
        }
      }
      backtrack(it);


      // Convergence criteria
      if (conv_criteria(it, value_old, position_old)) {
        it->c.step++;

        // If we made it here, it means that IF there was an issue with large gradient, it was resolved (we had succesfull backtracking)
        it->c.diff_code = 0;
        converged = true;
        break;
      } else if (pfun->value == value_old) {
        // No fruitfull backtrack
        if (it->c.diff_code == 1) {
          // First attempt to fix large gradient failed
          it->c.diff_code = 2; // Signal to get_active to do something else
        } else if (it->c.diff_code == 2) {
          // Second attempt to fix large gradient also failed
          it->c.diff_code = 3; // Signal this to jerr
          converged = true;
          break;
        } else {
          // Large gradient not detected, nothing unusual
          converged = true;
          break;
        }
      } else {
        // Neither of above --> continue

        // If we made it here, it means that there was an issue with large gradient, but it was resolved (we had succesfull backtracking)
        it->c.diff_code = 0;

        it->c.step ++;
        position_old = it->position;
      }
    }

    // if (trace) {
    //   Rcpp::Rcout.precision(12);
    //   Rcpp::Rcout << "\t Loss after " << it->c.step << " step(s): " << std::fixed << pfun->value << std::endl;
    //   Rcpp::Rcout.precision(6);
    // }
  } else {
    // If fixed it automatically has converged
    converged = true;
  }

  return converged;
}

// Runs optimise on each parameter in cycle
void optim::optimise(void) {
  // Vector telling if a parameter converged (conditioned on the others)
  arma::uvec converged = arma::zeros<arma::uvec>(pfun->vparam.size());

  // Get maximal number of cycles
  unsigned int cycle_MAX = 0;
  std::vector<param>::iterator param_it;
  for (param_it = pfun->vparam.begin(); param_it != pfun->vparam.end(); param_it++) {
    if (cycle_MAX < (unsigned int) param_it->c.step_max / param_it->c.step_cycle) {
      cycle_MAX = (unsigned int) param_it->c.step_max / param_it->c.step_cycle;
    }
  }
  cycle_MAX++;  // Interger division rounds down

  /* Flag indicating if this is potentially the last cycle (meaning previous cycle had all converged)
   * so that we may run a last screening.
   */
  bool last_cycle = false;

  // Check if any parameters are screened
  bool any_screen = false;
  for (param_it = pfun->vparam.begin(); param_it != pfun->vparam.end(); param_it++) {
    if (param_it->screen) {
      any_screen = true;
      break;
    }
  }

  // Cycle loop
  for (unsigned int i = 0; i < cycle_MAX; i++) {
    // Set all converged to false
    converged.zeros();

    // Report
    double pen = 0.0, pen_total = 0.0;
    if (trace) {
      Rcpp::Rcout.precision(12);
      Rcpp::Rcout << "Optimisation problem " << id << ". Starting cycle " << i + 1 << std::endl;
      Rcpp::Rcout.precision(6);
    }

    // Parameter loop
    for (param_it = pfun->vparam.begin(); param_it != pfun->vparam.end(); param_it++) {
      // Run optimise on that parameter and store convergence code
      converged(param_it - pfun->vparam.begin()) = optimise(param_it);

      // Report
      if (trace) {
        Rcpp::Rcout.precision(12);
        pen = param_it->preg->penalty(param_it->position);
        Rcpp::Rcout << "\t Parameter " << param_it - pfun->vparam.begin() + 1 << std::endl <<
          "\t\t Size of active set: " << param_it->active.n_elem << std::endl <<
            "\t\t Steps:\t\t" << param_it->c.step << std::endl <<
            "\t\t Loss:\t\t" << std::fixed << pfun->value - pen << std::endl <<
              "\t\t Penalty:\t" << pen << std::endl;
                // "\t\t 2Norm:\t\t" << arma::norm(param_it->position, 2) << std::endl <<
                //   "\t\t Diff:\t\t" << param_it->differential.t() <<
                //     "\t\t Dir:\t\t" << direction.t() <<
                //       "\t\t Tau:\t\t" << param_it->c.tau << std::endl;
        Rcpp::Rcout.precision(6);
        pen_total += pen;
      }
    }

    // Report total after cycle
    if (trace) {
      Rcpp::Rcout.precision(12);
      Rcpp::Rcout << "Optimisation problem " << id << ". Finished cycle " << i + 1 << std::endl <<
        "\t Loss:\t\t" << std::fixed << pfun->value - pen << std::endl <<
          "\t Penalty:\t" << std::fixed << pen_total << std::endl <<
            "\t Total:\t\t" << std::fixed << pfun->value - pen + pen_total << std::endl << "\n";
      Rcpp::Rcout.precision(6);
    }

    // Convergence check
    if (arma::all(converged) && (last_cycle || !any_screen)) {
      // If all converged and either it was marked as last cycle OR no screening then were done
      break;
    } else if (arma::all(converged) && any_screen) {
      /* Else if all converged and there are screenings (thus, last_cycle is false)
       *
       * Then it did not go through this part in the previous cycle.
       *
       * Run a last screening on all screenable parameters and mark it for next cycle.
       */
      for (param_it = pfun->vparam.begin(); param_it != pfun->vparam.end(); param_it++) {
        if (param_it->screen) {
          get_and_set_active(param_it);
        }
      }

      // Mark it for next cycle
      last_cycle = true;
    } else {
      // Otherwise not all converged, run another cycle and reset last_cycle
      last_cycle = false;
    }
  }

  // At very end we set_active to restore the last changes back to param_full
  for (param_it = pfun->vparam.begin(); param_it != pfun->vparam.end(); param_it++) {
    // Store last step in position_full
    param_it->reset_active();
    pfun->set_active(param_it);
  }

  // Save unpenalised function value in pfun->value
  pfun->evaluate();
}

// Set lambda
void optim::set_lambda(double lambda_) {
  pfun->set_lambda(lambda_);
}



//=================================
// Apply to vectors of optim
//=================================

// move: moves all particles (runs optimise or screen depending on input)
void move(std::vector<optim> &voptim) {
  // Reset optim_parameters
  std::for_each(voptim.begin(), voptim.end(), std::mem_fun_ref(&optim::reset_optim_param));

  // Move particles by optimising each of them
  std::for_each(voptim.begin(), voptim.end(), std::mem_fun_ref((void (optim::*)(void)) &optim::optimise));
}

// set_lambda: sets new lambda across optims
void set_lambda(std::vector<optim> &voptim, double lambda_) {
  std::for_each(voptim.begin(), voptim.end(), std::bind2nd(std::mem_fun_ref(&optim::set_lambda), lambda_));
}





