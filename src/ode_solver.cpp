//=============================================
/* ode_solver.cpp
 *
 * Content:
 * - state struct definition
 * - solver and subclasses definitions
 * - ode_solve
 *
 */
//=============================================


#include "ode_solver.h"

// Define cube multiplied with vector (collapses cube along 3rd dimension)
arma::mat operator*(arma::cube c, const arma::vec &v) {
  if (v.n_elem != c.n_slices) {
    Rcpp::stop("Dimension mismatch!");
  } else {
    for (unsigned int i = 0; i < v.n_elem; i ++) {
      c.slice(i) *= v(i);
    }
  }
  return arma::sum(c, 2);
}



//=================================
// state structure
//=================================

// Constructors
state::state (double t_, arma::vec x_, field* pfield_, bool score_, bool exact_, arma::mat x_dparam_) :
  t(t_), x(x_), score(score_), exact(exact_), x_dparam(x_dparam_),
  pfield(pfield_) {}

state::state (double t_, arma::vec x_, std::vector<param>::iterator param_it_, field* pfield_, bool exact_) :
  t(t_), x(x_), score(true), exact(exact_), pfield(pfield_), param_it(param_it_) {
  if (is_x0()) {
    // If the iterator points to the first parameter (x0)
    x_dparam = arma::eye(pfield->d, pfield->d);
  } else {
    x_dparam = arma::zeros(pfield->d, param_it->p);
  }
}

bool state::is_x0(void) {
  return &(*param_it) == &(*pfield->pvparam->begin());
}

// Reset sensitivity, which checks what iterator to use
void state::reset_dparam(void) {
  if (is_x0()) {
    // If the iterator points to the first parameter (x0)
    x_dparam = arma::eye(pfield->d, pfield->d);
  } else {
    x_dparam.fill(0);
  }
}

/* NOTES:
 * In all functions in solver class that uses state, e.g., solve, step, etc,
 * creates a state instance and then it survives all the way through, even in
 * step, where new proposal states are made, but they are ultimately only
 * copied over for some member variables.
 * Thus we could let pfield and iterator becomes members of state, since they are rarely copied.
 */



//=================================
// solver class
//=================================

// Constructor
solver::solver (unsigned int step_max_) :
  step_max(step_max_),
  code_conv(0) {}

// Reseting convergence code
void solver::refresh_code (void) {
  code_conv = 0;
}

/* step function which (for a given solver) takes position (t, x, x_dparam, etc) and moves the ode system to t_end.
 * The only thing (except 'reset') that is called in loss
 */
void solver::step (state &state_, double t_end) {
  // Only if positive increase we do something (else state is returned as it is)
  if (state_.t < t_end) {
    // Initialise internal time and step counter
    double t_old = state_.t;
    unsigned int step = 0;

    // Function does nothing unless current x is finite
    if (arma::is_finite(state_.x)) {
      while (step < step_max) {
        t_old = state_.t;                               // Save the old time, before it is updated by proposal

        // Rcpp::Rcout << "solver::step -> Get proposal" << std::endl;

        state state_proposed = proposal(state_, step);  // Get proposed z

        // Rcpp::Rcout << "solver::step -> Got proposal" << std::endl;

        if (arma::is_finite(state_proposed.x)) {// If proposed position x is finite proceed
          if (state_proposed.t < t_end) {         // If the new updated time point has not yet passed the t_end-mark, update and repeat
            state_.x = state_proposed.x;
            if (state_.score) {
              state_.x_dparam  = state_proposed.x_dparam;
            }
            step ++;
          } else {                                // Else linearly interpolate to t_end
            state_.t = t_end;
            state_.x += (state_proposed.x - state_.x) * (t_end - t_old) / (state_proposed.t - t_old);
            if (state_.score) {
              state_.x_dparam  += (state_proposed.x_dparam  - state_.x_dparam)  * (t_end - t_old) / (state_proposed.t - t_old);
            }
            break;
          }
        } else {                                  // Else if z proposal is not finite return nans
          state_.x.fill(arma::datum::nan);
          state_.x_dparam.fill(arma::datum::nan);
          break;
        }
      }
      if (step >= step_max) { // If failed to produce valid position value at dt
        code_conv = 1;
      }
    }
  }
}

// Using the step function we can produce any desired trajectory of the system.
arma::mat solver::solve(field* pfield, arma::vec time) {
  refresh_code();
  reset();

  unsigned int n = time.size();                     // Dimensions
  unsigned int series = 0;                          // Current number of series
  arma::mat trajectory  = arma::mat(pfield->d, n);  // Return matrix and position object
  arma::mat x0 = pfield->pvparam->begin()->get_position();
  state state_ (time(0), x0.col(series), pfield);   // By default score = false
  trajectory.col(0) = state_.x;
  // Rcpp::Rcout << "solver::solve prepared for moving" << std::endl;


  // Move-trajectory-forward-loop
  for (unsigned int i = 0; i < n - 1; i ++) {
    if (time(i + 1) < time(i)) {  // A decrease in time is convention for a new series, re-initialise
      reset();      // reset local control variables
      series ++;
      state_.t = time(i + 1);
      state_.x = x0.col(series);
      trajectory.col(i + 1) = state_.x;

      continue;
    } else {        // Otherwise, take step to next time point
      // Rcpp::Rcout << "Called step for " << i << "th time point." << std::endl;
      step(state_, time(i + 1));
      trajectory.col(i + 1) = state_.x;
    }
  }

  reset();  // local control variables are reset
  return(trajectory.t());
}

// Using the step function we can produce any desired trajectory and sensitivity of the system, not possible to reduce
arma::mat solver::solve(field* pfield, std::vector<param>::iterator param_it, arma::vec time, arma::cube &sens_dparam, bool approx) {
  // Refresh stuff
  refresh_code();
  reset();

  unsigned int n = time.size();                   // Dimensions
  unsigned int series = 0;                        // Current number of series

  arma::mat trajectory  = arma::mat(pfield->d, n);  // Return matrix and position object
  arma::mat x0 = pfield->pvparam->begin()->get_position();
  state state_ (time(0), x0.col(series), param_it, pfield, !approx);

  // Initialise trajectory and sensitivity
  trajectory.col(0) = state_.x;
  sens_dparam.clear();  sens_dparam.set_size(n, state_.x_dparam.n_rows, state_.x_dparam.n_cols);
  sens_dparam.tube(arma::span(0), arma::span::all)  = state_.x_dparam;

  // Move-trajectory-forward-loop
  for (unsigned int i = 0; i < n - 1; i ++) {
    if (time(i + 1) < time(i)) {  // A decrease in time is convention for a new series, re-initialise
      reset();      // reset local control variables
      series ++;

      // Set new state
      state_.t = time(i + 1);
      state_.x = x0.col(series);
      state_.reset_dparam();

      // Load off state
      trajectory.col(i + 1) = state_.x;
      sens_dparam.tube(arma::span(i + 1), arma::span::all)  = state_.x_dparam;

      continue;
    } else {        // Otherwise, take step to next time point
      // Move to time(i + 1)
      step(state_, time(i + 1));

      // Load off state
      trajectory.col(i + 1) = state_.x;
      sens_dparam.tube(arma::span(i + 1), arma::span::all)  = state_.x_dparam;
    }
  }

  reset();  // local control variables are reset
  return(trajectory.t());
}

// Using the step function we can produce sensitivity X, i.e., stacked version of (\partial_param x(t_j, param))_j
/*
 *
 */
arma::mat solver::sensitivity_param(field* pfield, std::vector<param>::iterator param_it, arma::vec time, std::vector<sc>::iterator sc_it, bool approx) {
  // Refresh stuff
  refresh_code();
  reset();

  unsigned int n = time.size();   // Dimensions
  unsigned int series = 0;        // Current number of series

  arma::mat x0 = pfield->pvparam->begin()->get_position();
  state state_ (time(0), x0.col(series), param_it, pfield, !approx);
  arma::mat sens_dparam = arma::zeros(pfield->d * n, param_it->p);  // Return matrix and position object

  if (state_.is_x0()) {
    sens_dparam(arma::span(0, pfield->d - 1), arma::span(0, pfield->d - 1)) = state_.x_dparam;
  } else {
    if (sc_it->exist) {
      sens_dparam(arma::span(0, pfield->d - 1), arma::span::all) = state_.x_dparam.each_row() % sc_it->col(series).t();
    } else {
      sens_dparam(arma::span(0, pfield->d - 1), arma::span::all) = state_.x_dparam;
    }
  }

  // Move-trajectory-forward-loop
  for (unsigned int i = 0; i < n - 1; i ++) {
    if (time(i + 1) < time(i)) {  // A decrease in time is convention for a new series, re-initialise
      reset();      // reset local control variables
      series ++;

      // Set new state
      state_.t = time(i + 1);
      state_.x = x0.col(series);
      state_.reset_dparam();

      // Load off state
      if (state_.is_x0()) {
        sens_dparam(arma::span((i + 1) * pfield->d, (i + 2) * pfield->d - 1), arma::span(0 + series * pfield->d, pfield->d - 1 + series * pfield->d)) = state_.x_dparam;
      } else {
        if (sc_it->exist) {
          sens_dparam(arma::span((i + 1) * pfield->d, (i + 2) * pfield->d - 1), arma::span::all) = state_.x_dparam.each_row() % sc_it->col(series).t();
        } else {
          sens_dparam(arma::span((i + 1) * pfield->d, (i + 2) * pfield->d - 1), arma::span::all) = state_.x_dparam;
        }
      }

      continue;
    } else {        // Otherwise, take step to next time point
      // Move to time(i + 1)
      step(state_, time(i + 1));

      // Load off state
      if (state_.is_x0()) {
        sens_dparam(arma::span((i + 1) * pfield->d, (i + 2) * pfield->d - 1), arma::span(0 + series * pfield->d, pfield->d - 1 + series * pfield->d)) = state_.x_dparam;
      } else {
        if (sc_it->exist) {
          sens_dparam(arma::span((i + 1) * pfield->d, (i + 2) * pfield->d - 1), arma::span::all) = state_.x_dparam.each_row() % sc_it->col(series).t();
        } else {
          sens_dparam(arma::span((i + 1) * pfield->d, (i + 2) * pfield->d - 1), arma::span::all) = state_.x_dparam;
        }
      }
    }
  }

  reset();  // local control variables are reset
  // pfield->restore();

  return(sens_dparam);
}




//=================================
// rkemp subclass
//=================================

// A bullet-proof version of ratio-calculator for choosing new h (adaptive time step)
double rkemp::ratio(arma::vec w, arma::vec z) {
  double ratio = 1;
  if (arma::is_finite(z) && arma::is_finite(w)) {
    if (arma::norm(z - w) > 1e-10 * Tol) {              // Denominator non-zero
      ratio = Tol * arma::norm(w) / arma::norm(z - w);
    } else {
      if (arma::norm(w) > 1e-10) {                      // Enumerator non-zero
        ratio = Tol * arma::norm(w) / 1e-10;
      } else {                                          // Both denominator and enumerator zero
        ratio = exp(5 * log(1 / .79));
      }
    }
  } else {                                              // If infinite/nan set to lowest possible
    ratio = .1;
  }
  ratio = std::max(.1, std::min(1e6, ratio));           // Bound between .1 and 1e6 (corresponds roughly to h = .5 * h and h = 10 * h)
  return ratio;
}

// Produces a proposed state in rkf45, also updates h, step and time in state
state rkemp::proposal(state &state_, unsigned int &step) {
  // Rcpp::Rcout << "rkemp::proposal -> Initialise xs fxs, w. pfield->d: " << state_.pfield->d << " b.t.n_elem: " << b.t.n_elem << " EXACT: " << state_.exact << std::endl;

  // Initialise
  arma::mat xs(state_.pfield->d, b.t.n_elem, arma::fill::zeros);
  arma::mat fxs(state_.pfield->d, b.t.n_elem, arma::fill::zeros);
  arma::vec w(state_.pfield->d), z(state_.pfield->d);
  double ratio_ = 1;

  // Rkf45 update for (potentially) a couple of times
  while (step < step_max) {
    Rcpp::checkUserInterrupt();     // Possible to stop calculations from R
    step ++;
    fxs.fill(0);  // Cleanse (if butcher tableaus are "upper-triangular" this step is actually not necessary, but better safe than sorry)

    // Actual Runge Kutta Embedded Pair loop
    if (state_.score & state_.exact) {
      // If we are to do exact differential we need to store xs
      xs.fill(0);
      for (unsigned int i = 0; i < b.t.n_elem; i ++) {
        xs.col(i)   = state_.x + h * fxs * b.x.col(i);
        fxs.col(i)  = state_.pfield->f(state_.t + h * b.t(i), xs.col(i));
      }
    } else {
      for (unsigned int i = 0; i < b.t.n_elem; i ++) {
        fxs.col(i)  = state_.pfield->f(state_.t + h * b.t(i), state_.x + h * fxs * b.x.col(i));
      }
    }
    w = state_.x + h * fxs * b.x.col(b.t.n_elem);
    z = state_.x + h * fxs * b.x.col(b.t.n_elem + 1);
    ratio_ = ratio(w, z);

    if (ratio_ > 1 || step >= step_max) {     // If succesfull with h, break out with the proposed update
      break;
    } else {                                  // Else try again with updated h
      h *= 0.8 * exp(log(ratio_) * .2);
    }
  }
  state_.t += h;                              // Update time (before updating h for next step), but not state_.x (as we do not know if proposal is accepted)

  state state_proposed (state_.t, z);         /* Create state object that is returned, note it has empty field pointer */
  if (state_.exact) state_proposed.exact = true;

  if (state_.score) {                         // Add score if needed
    // Rcpp::Rcout << "rkemp::proposal -> do score: " << std::endl;
    state_proposed.score = true;

    if (arma::is_finite(state_.x_dparam)) {
      if (state_.exact) {
        // If exact solution desired
        arma::cube fxs_dparam(state_.x_dparam.n_rows, state_.x_dparam.n_cols, b.t.n_elem, arma::fill::zeros);
        for (unsigned int i = 0; i < b.t.n_elem; i ++) {
          // Rcpp::Rcout << "rkemp::proposal -> exact score. xs.col(i): " << xs.col(i) << std::endl <<
          //   "param: " << state_.param_it->position.t() << std::endl;
          // Rcpp::Rcout << "rkemp::proposal -> exact score. f_dx: " << state_.pfield->f_dx(state_.t + h * b.t(i), xs.col(i)) << std::endl;
          // Rcpp::Rcout << "rkemp::proposal -> fxs_dparam*b.x" << fxs_dparam * b.x.col(i) << std::endl;
          fxs_dparam.slice(i) = state_.pfield->f_dx(state_.t + h * b.t(i), xs.col(i)) *
            (state_.x_dparam + h * (fxs_dparam * b.x.col(i)));
          if (!state_.is_x0()) {
            // Rcpp::Rcout << "rkemp::proposal -> additionals when param not x0: " << std::endl;
            fxs_dparam.slice(i) += state_.pfield->f_dparam(state_.t + h * b.t(i), xs.col(i), state_.param_it);
          }
        }
        state_proposed.x_dparam = state_.x_dparam + h * (fxs_dparam * b.x.col(b.t.n_elem + 1));  // Corresponds to z_dparam
      } else {
        // If only approximate solution desired
        // Rcpp::Rcout << "rkemp::proposal -> not exact score" << std::endl;
        state_proposed.x_dparam = state_.x_dparam + h * state_.pfield->f_dx(state_.t - h / 2., (state_.x + z) / 2.) * state_.x_dparam;
        if (!state_.is_x0()) {
          // If it is not x0 param, add the f_dparam part
          // Rcpp::Rcout << "rkemp::proposal -> additionals when param not x0: " << std::endl;
          state_proposed.x_dparam += h * state_.pfield->f_dparam(state_.t - h / 2., (state_.x + z) / 2., state_.param_it);
        }
      }
    } else {
      // If non finite, then repeat from state_ to warn future about not continuing
      state_proposed.x_dparam  = state_.x_dparam;
    }
  }
  // Rcpp::Rcout << "rkemp::proposal -> adjust step size" << std::endl;

  h *= 0.8 * exp(log(ratio_) * .2);   // Updating h for next step (before this line, h matches z, w, and k1...k6)
  return state_proposed;
}

// Constructor
rkemp::rkemp(std::string str_, unsigned int step_max_, double Tol_, double h_init_) :
  solver(step_max_),
  b(butcher(str_)),
  Tol(Tol_),
  h_init(h_init_),
  h(h_init_) {}

// Clone function
rkemp* rkemp::shallow_Clone() {
  return &(*this);  // No need to call "new", since no pointers in rkemp or solver.
}
rkemp* rkemp::deep_Clone() {
  rkemp* prkemp = new rkemp(*this);
  return prkemp;
}

void rkemp::reset(void) {
  h = h_init;
}




//=================================
// create_psolver and solve_ode
//=================================

/* ode_struct:
 *
 * mak: sparse (logical), A, B
 * plk: sparse (logical), A
 *
 * solver_ctrl: step_max, tol, h_init
 * rs: list of reg's
 */
solver* create_psolver(Rcpp::List ode_struct) {
  Rcpp::List s = Rcpp::as<Rcpp::List>(ode_struct["s"]);
  Rcpp::List ctrl = Rcpp::as<Rcpp::List>(s["ctrl"]);

  solver* psolver = NULL;

  rkemp r(Rcpp::as<std::string>(s["name"]),
    Rcpp::as<unsigned int>(ctrl["step_max"]),
    Rcpp::as<double>(ctrl["tol"]),
    Rcpp::as<double>(ctrl["h_init"]));
  solver* prkemp = &r;

  psolver = prkemp->deep_Clone();

  return psolver;
}


// [[Rcpp::export(ode_solve)]]
Rcpp::List ode_solve(Rcpp::List ode_struct, arma::mat x0, Rcpp::List param_, arma::vec time,
  bool sensitivity, bool approx) {

  // Create field pointer
  arma::mat A;
  arma::mat B;
  std::vector<param> vparam;
  field* pfield = create_pfield(ode_struct, A, B, vparam, x0.col(0));
  Rcpp::Nullable<Rcpp::List> pp = Rcpp::as< Rcpp::Nullable<Rcpp::List> >(param_[0]);
  pfield->fill_in_param(pp);

  // Prepare solver
  solver* psolver = create_psolver(ode_struct);

  // For sensitivity
  unsigned int n_series = arma::sum(arma::diff(time) < 0) + 1;

  // For trajectory
  arma::mat traj(0, pfield->d);
  arma::mat traj_block(0, pfield->d);

  // for x0 sensitiviyt
  arma::cube ret_dx0(time.n_elem, pfield->d, pfield->d);

  // for param sensitivity
  std::vector<arma::cube> ret_dparam(vparam.size() - 1);
  std::vector<arma::cube>::iterator cube_it;

  std::vector<param>::iterator param_it = vparam.begin();

  arma::cube sens_dparam; // Placeholder for differentials


  arma::uvec code_convs(n_series);

  // Jumps in time
  arma::uvec to_jumps = arma::find(arma::diff(time) < 0);
  arma::uvec vend(1); vend(0) = time.n_elem;
  to_jumps = arma::join_cols(to_jumps + 1, vend);
  vend(0) = 0;  to_jumps = arma::join_cols(vend, to_jumps);

  // Loop over series
  param_it = vparam.begin();
  Rcpp::List ret;
  // RUn solver
  // RUn solver
  if (sensitivity) {
    // Set sizes here to avoid bad alloc if sensitivity not desired
    ret_dx0 = arma::zeros<arma::cube>(time.n_elem, pfield->d, pfield->d);

    // for param sensitivity
    std::vector<param>::iterator param_it = vparam.begin();
    for (cube_it = ret_dparam.begin(); cube_it < ret_dparam.end(); cube_it++) {
      param_it++;
      *cube_it = arma::zeros<arma::cube>(time.n_elem, pfield->d, param_it->p);
    }

    // Loop over series
    for (unsigned int i_series = 0; i_series < n_series; ++i_series) {
      arma::vec to_sub = time(arma::span(to_jumps(i_series), to_jumps(i_series + 1) - 1));
      // Rcpp::Rcout << "to_sub " << to_sub.t() << std::endl;


      // fill in parameter
      vparam.begin()->position = x0.col(i_series);
      pp = Rcpp::as< Rcpp::Nullable<Rcpp::List> >(param_[i_series]);
      pfield->fill_in_param(pp);

      // x0 dparam
      param_it = vparam.begin();
      traj_block = psolver->solve(pfield, param_it, to_sub, sens_dparam, approx);



      traj = arma::join_cols(traj, traj_block);

      code_convs(i_series) = psolver->code_conv;
      psolver->code_conv = 0;

      ret_dx0.tube(arma::span(to_jumps(i_series), to_jumps(i_series + 1) - 1), arma::span::all) = sens_dparam;

      // Other parameters
      cube_it = ret_dparam.begin();
      for (unsigned int i = 0; i < vparam.size() - 1; i++) {
        param_it++;
        traj_block = psolver->solve(pfield, param_it, to_sub, sens_dparam, approx);
        cube_it->tube(arma::span(to_jumps(i_series), to_jumps(i_series + 1) - 1), arma::span::all) = sens_dparam;
        cube_it++;
      }
    }

    ret["trajectory"] = traj;
    ret["conv_code"]  = code_convs;
    ret["sens_dx0"]   = ret_dx0;

    cube_it = ret_dparam.begin();
    Rcpp::List params_(vparam.size() - 1);
    for (unsigned int i = 0; i < params_.size(); i++) {
      arma::cube c = *cube_it;
      params_[i] = c;
      cube_it++;
    }
    ret["sens_dparam"] = params_;


  } else {

    // Loop over series
    for (unsigned int i_series = 0; i_series < n_series; ++i_series) {
      arma::vec to_sub = time(arma::span(to_jumps(i_series), to_jumps(i_series + 1) - 1));

      // fill in parameter
      vparam.begin()->position = x0.col(i_series);
      pp = Rcpp::as< Rcpp::Nullable<Rcpp::List> >(param_[i_series]);
      pfield->fill_in_param(pp);

      // Reduce
      for (param_it = vparam.begin() + 1; param_it < vparam.end(); param_it++) {
        param_it->set_active(param_it->preg->get_active(param_it->position_full, arma::zeros(param_it->p_full), param_it->c.tau_min));
      }
      pfield->set_active();

      // Rcpp::Rcout << "Call solver::solve" << std::endl;

      traj_block = psolver->solve(pfield, to_sub);
      // Rcpp::Rcout << "Size of traj_block" << traj_block.n_rows << " " << traj_block.n_cols << std::endl;
      // Rcpp::Rcout << "Size of traj" << traj.n_rows << " " << traj.n_cols << std::endl;

      traj = arma::join_cols(traj, traj_block);
      code_convs(i_series) = psolver->code_conv;
      psolver->code_conv = 0;

      // Restore
      for (param_it = vparam.begin() + 1; param_it < vparam.end(); param_it++) {
        param_it->reset_active();
      }
      pfield->set_active();

    }

    ret["trajectory"] = traj;
    ret["conv_code"]  = code_convs;
  }

  return ret;
}


