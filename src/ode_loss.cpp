//=============================================
/* ode_loss.cpp
 *
 * Content:
 * - loss and subclasses definitions
 * - aim_design
 * - test_loss
 *
 *
 */
//=============================================

#include <algorithm>
#include "ode_loss.h"



//=================================
// loss class
//=================================

/* nan_guards replace weights corresponding to y = nan (missing obs) with 0. And replace those y with 0.
 * Thus if fitted values are not NaN (which could happen due to x0 or param out of bounds)
 * then that term will not affect loss. If, however, fitted value is NaN, then loss ends up NaN
 * and backtracking sends it back on track (i.e., param and x0 inside the bounds).
 */
arma::mat ode_loss::w_nan_guard (arma::mat w_, arma::mat y_) {
  w_(arma::find_nonfinite(y_)).fill(0);
  return w_;
}
arma::mat ode_loss::y_nan_guard (arma::mat y_) {
  y_.replace(arma::datum::nan, 0);
  return y_;
}

// Constructor
ode_loss::ode_loss(Rcpp::List ode_struct, arma::mat &A, arma::mat &B, arma::mat x0, Rcpp::Nullable<Rcpp::List> param_, Rcpp::List opt_struct, std::vector<sc>* pvsc_full_, likeli* plikeli_) :
  fun(std::vector<param>()),
  pvsc_full(pvsc_full_),
  plikeli(plikeli_) {
  // Scales and contexts
  vsc = *pvsc_full;

  // Field (also sets up vparam and fills it)
  pfield = create_pfield(ode_struct, A, B, vparam, x0);

  // Fill it in if it is not null
  pfield->fill_in_param(param_);

  // exact_diff (overwrites ode_struct!!)
  ode_struct = Rcpp::as<Rcpp::List>(ode_struct["rs"]);
  exact_diff = arma::zeros<arma::uvec>(ode_struct.size());
  for (unsigned int i = 0; i < ode_struct.size(); i++) {
    Rcpp::List reg_struct = Rcpp::as<Rcpp::List>(ode_struct[i]);
    exact_diff(i) = Rcpp::as<bool>(reg_struct["exact_gradient"]);
  }

  if (exact_diff.n_elem != vsc.size() + 1) {
    Rcpp::stop("error in ode_loss::ode_loss -> the number of scales must be one less than the total number of parameters.");
  }

  // Observations
  arma::mat y_ = Rcpp::as<arma::mat>(opt_struct["y"]);
  time = y_.col(0);
  y_ = y_.t();
  y = y_nan_guard(y_.rows(1, pfield->d));
  n_obs = y.n_cols;
  n_series = arma::sum(arma::diff(time) < 0) + 1;

  // Weights
  Rcpp::Nullable<arma::mat> w_ = Rcpp::as< Rcpp::Nullable<arma::mat> >(opt_struct["weights"]);
  w_exist = w_.isNotNull();
  w = w_exist ? w_nan_guard(arma::trans(Rcpp::as<arma::mat>(w_)), y) : arma::mat();

  // Unscaled positions
  unscaled_pos = std::vector<arma::vec>(vparam.size() - 1);
}

// Non-default assignment, since pfield's pointer to vec<param> must be the in fun
ode_loss::ode_loss(const ode_loss& other) :
  fun(other),
  time(other.time),
  n_obs(other.n_obs),
  n_series(other.n_series),
  pvsc_full(other.pvsc_full), // Normal copy, since it is "alive" in bull/bronc etc.
  vsc(other.vsc),
  plikeli(other.plikeli),     // Also normal, since it has no members
  y(other.y),
  w_exist(other.w_exist),
  w(other.w),
  pfield(other.pfield),
  exact_diff(other.exact_diff),
  unscaled_pos(other.unscaled_pos) {
  pfield->pvparam = &vparam;

}

// Non-default assignment, since pfield's pointer to vec<param> must be the in fun
ode_loss& ode_loss::operator=(const ode_loss& rhs) {
  // Rcpp::Rcout << "Calling ode_loss assignment" << std::endl;

  if (&rhs != this) {
    time = rhs.time;
    n_obs = rhs.n_obs;
    n_series = rhs.n_series;
    pvsc_full = rhs.pvsc_full;  // Normal copy, since it is "alive" in bull/bronc etc.
    vsc = rhs.vsc;
    plikeli = rhs.plikeli;      // Also normal, since it has no members
    y = rhs.y;
    w_exist = rhs.w_exist;
    w = rhs.w;
    pfield = rhs.pfield;
    pfield->pvparam = &vparam;  // pfield must point to the correct parameter vector
    exact_diff = rhs.exact_diff;
    unscaled_pos = rhs.unscaled_pos;
  }
  return *this;
}


// For dfs, no need to have "exclude" in, as reduce is called
// Stein DF,
arma::mat ode_loss::get_dfs(arma::vec lambdas) {
  // Save old lambda and old actives
  double lambda_old = get_lambda();
  std::vector<arma::uvec> vactive(vparam.size());
  std::vector<arma::uvec>::iterator act_it = vactive.begin();
  std::vector<param>::iterator param_it = vparam.begin();
  for (param_it = vparam.begin(); param_it < vparam.end(); param_it++) {
    // Store the old ones
    *act_it = param_it->active;

    // Set new according to non-zeros in position
    param_it->reset_active();
    param_it->set_active(arma::find(param_it->position));
    set_active(param_it); // adapt ode_loss to it (sc etc)

    act_it++;
  }

  // Return object
  arma::mat dfs = arma::zeros(vparam.size(), lambdas.n_elem);

  // Lambda loop
  for (unsigned int j = 0; j < lambdas.n_elem; j++) {
    // Set lambda on all params
    set_lambda(lambdas(j));

    param_it = vparam.begin();
    for (unsigned int i = 0; i < vparam.size(); i++) {
      if (!param_it->fixed) {
        // Rcpp::Rcout << "ode_loss::get_dfs -> calling get_X4dfs." << std::endl;
        dfs(i, j) = param_it->dfs(get_X4dfs(param_it));
      }
      param_it++;
    }
  }

  // Restore lambda and actives
  set_lambda(lambda_old);
  act_it = vactive.begin();
  for (param_it = vparam.begin(); param_it < vparam.end(); param_it++) {
    param_it->active = *act_it;
    act_it++;
  }

  return dfs;
}

// Makes sc match the active
void ode_loss::set_active_base(std::vector<param>::iterator param_it) {
  // Scales and contexts
  if (&(*param_it) != &(*vparam.begin())) {
    std::vector<sc>::iterator sc_it = vsc.begin() + (param_it - vparam.begin() - 1);
    std::vector<sc>::iterator sc_full_it = pvsc_full->begin() + (param_it - vparam.begin() - 1);
    sc_full_it->transfer(*sc_it, param_it->active);
  }

  // field
  pfield->set_active(param_it);
}

void ode_loss::store_positions(void) {
  std::vector<param>::iterator it = vparam.begin() + 1;
  for (std::vector<arma::vec>::iterator pos_it = unscaled_pos.begin(); pos_it < unscaled_pos.end(); pos_it++) {
    *pos_it = it->position;
    it++;
  }
}

void ode_loss::scale_positions(unsigned int i_series) {
  std::vector<param>::iterator it = vparam.begin() + 1;
  std::vector<sc>::iterator sc_it = vsc.begin();
  for (std::vector<arma::vec>::iterator pos_it = unscaled_pos.begin(); pos_it < unscaled_pos.end(); pos_it++) {
    it->position = *pos_it;
    if (sc_it->exist) {
      it->position %= sc_it->col(i_series);
    }
    sc_it++;
    it++;
  }
}

void ode_loss::restore_positions(void) {
  std::vector<param>::iterator it = vparam.begin() + 1;
  for (std::vector<arma::vec>::iterator pos_it = unscaled_pos.begin(); pos_it < unscaled_pos.end(); pos_it++) {
    it->position = *pos_it;
    it++;
  }
}




//=================================
// exact subclass
//=================================

exact::exact(Rcpp::List ode_struct, arma::mat &A, arma::mat &B, arma::mat x0, Rcpp::Nullable<Rcpp::List> param_, Rcpp::List opt_struct, std::vector<sc>* pvsc_full_, likeli* plikeli_) :
  ode_loss(ode_struct, A, B, x0, param_, opt_struct, pvsc_full_, plikeli_),
  psolver(create_psolver(ode_struct)) {
  // Update value
  evaluate();

  }

// Non-default assignment, since pfield's pointer to vec<param> must be the in fun
exact::exact(const exact& other) :
  ode_loss(other),
  psolver(other.psolver->deep_Clone()) {
  // Rcpp::Rcout << "Called exact copy" << std::endl;
}
exact& exact::operator=(const exact& rhs) {
  if (&rhs != this) {
    ode_loss::operator =(rhs);
    psolver  = rhs.psolver->deep_Clone();
  }
  return *this;
}


// Updates value (by applying -loglikelihood iteratively over the observations)
void exact::evaluate(void) {
  // Rcpp::Rcout << "Calling exact::evaluate" << std::endl;
  psolver->refresh_code();  // Reset convergence code in solver
  psolver->reset();         // Reset control parameters in solver (e.g., h)
  // Rcpp::Rcout << "Check that pfield points to correct parameter." << std::endl;


  if (pfield->pvparam != &vparam) {
    Rcpp::warning(" warning in exact::evaluate -> field and parameters do not agree");
  }
  // Rcpp::Rcout << "exact::evaluate -> Get x0" << std::endl;
  // Rcpp::Rcout << "Calling exact::evaluate -> resat?" << std::endl;


  // Get x0
  arma::mat x0 = vparam.begin()->get_position();
  // Rcpp::Rcout << "exact::evaluate -> Got x0" << std::endl;
  // Rcpp::Rcout << "exact::evaluate -> Store unscaled positions" << std::endl;


  // The effective parameter is eta = schur.prod(param, sc) (but param is penalised)
  store_positions();
  // Rcpp::Rcout << "exact::evaluate -> Stored unscaled positions" << std::endl;


  scale_positions(0);
  // Rcpp::Rcout << "exact::evaluate -> Scaled positions" << std::endl;


  // Prepare series
  unsigned int series = 0;                          // Current number of series
  state state_ (time(0), x0.col(series), pfield);   // By default score = false
  // Rcpp::Rcout << "exact::evaluate -> Created state. State size " << state_.x.n_elem << " Dparam size " << state_.x_dparam.size() << std::endl;

  // Loss for first observation
  value = 0.0;
  if (w_exist) {
    value += plikeli->eval(y.col(0), state_.x, w.col(0));
  } else {
    value += plikeli->eval(y.col(0), state_.x);
  }
  // Rcpp::Rcout << "exact::evaluate -> First call to plikeli" << std::endl;

  // Loop through the observations
  for (unsigned int i = 0; i < n_obs - 1; ++ i) {
    // A decrease in time is convention for a new series, re-initialise
    if (time(i + 1) < time(i)) {
      psolver->reset();
      series ++;
      state_.t = time(i + 1);
      state_.x = x0.col(series);

      // Scales positions with new scale
      scale_positions(series);
    } else {
      // Rcpp::Rcout << "exact::evaluate -> call step for " << i << "th time." << std::endl;
      // Proceed to next time point
      psolver->step(state_, time(i + 1));
    }

    // Update loss
    if (w_exist) {
      value += plikeli->eval(y.col(i + 1), state_.x, w.col(i + 1));
    } else {
      value += plikeli->eval(y.col(i + 1), state_.x);
    }
  }

  // MEAN over time points - n_series (if positive)
  if (n_obs > n_series) {
    value /= n_obs - n_series;
  }

  // Return original unscaled ones
  restore_positions();

  psolver->reset(); // VERY IMPORTANT (resets, e.g., discretisation parameter, so that the same conditions hold in next loss evaluation)
}


// Updates the differential for a given parameter (by applying -loglikelihood iteratively over the observations)
void exact::evaluate_differential(std::vector<param>::iterator param_it) {
  psolver->refresh_code();  // Reset convergence code in solver
  psolver->reset();         // Reset control parameters in solver (e.g., h)

  // Get x0
  arma::mat x0 = vparam.begin()->get_position();

  // Iterators
  std::vector<sc>::iterator sc_it;

  // The effective parameter is eta = schur.prod(param, sc) (but param is penalised)
  store_positions();
  scale_positions(0);

  // The temporary differential and gnhda (overwritten in each series)
  param_it->differential.zeros();
  arma::vec tmp_diff = param_it->differential;
  param_it->gnhda.zeros();
  arma::vec tmp_gnhda = param_it->gnhda;

  // Prepare series
  unsigned int series = 0;                          // Current number of series
  state state_ (time(0), x0.col(series), param_it, pfield, exact_diff(param_it - vparam.begin()));
  // Rcpp::Rcout << "exact::evaluate_differential -> set op state" << std::endl;

  // Differential for first observation
  if (w_exist) {
    if (state_.is_x0()) {
      tmp_diff(arma::span(0, pfield->d - 1)) += state_.x_dparam.t()  * plikeli->eval_dfit(y.col(0), state_.x, w.col(0));
      if (eval_gnhda) {
        tmp_gnhda(arma::span(0, pfield->d - 1)) += (arma::sum(state_.x_dparam % (plikeli->eval_d2fit(y.col(0), state_.x, w.col(0)) * state_.x_dparam), 0)).t();
      }
    } else {
      tmp_diff += state_.x_dparam.t()  * plikeli->eval_dfit(y.col(0), state_.x, w.col(0));
      if (eval_gnhda) {
        tmp_gnhda += (arma::sum(state_.x_dparam % (plikeli->eval_d2fit(y.col(0), state_.x, w.col(0)) * state_.x_dparam), 0)).t();
      }
    }
  } else {
    if (state_.is_x0()) {
      tmp_diff(arma::span(0, pfield->d - 1))  += state_.x_dparam.t()  * plikeli->eval_dfit(y.col(0), state_.x);
      if (eval_gnhda) {
        tmp_gnhda(arma::span(0, pfield->d - 1)) += (arma::sum(state_.x_dparam % (plikeli->eval_d2fit(y.col(0), state_.x) * state_.x_dparam), 0)).t();
      }
    } else {
      tmp_diff += state_.x_dparam.t()  * plikeli->eval_dfit(y.col(0), state_.x);
      if (eval_gnhda) {
        tmp_gnhda += (arma::sum(state_.x_dparam % (plikeli->eval_d2fit(y.col(0), state_.x) * state_.x_dparam), 0)).t();
      }
    }
  }

  // Loop through the observations
  // Rcpp::Rcout << "exact::evaluate_differential -> starting loop" << std::endl;
  for (unsigned int i = 0; i < n_obs - 1; ++ i) {
    // A decrease in time is convention for a new series, re-initialise
    if (time(i + 1) < time(i)) {
      // Scale the differentials etc with sc
      if (!state_.is_x0()) {
        sc_it = vsc.begin() + (param_it - vparam.begin() - 1);
        if (sc_it->exist) {
          tmp_diff %= sc_it->col(series);
          if (eval_gnhda) {
            tmp_gnhda %= arma::square(sc_it->col(series));
          }
        }
      }
      param_it->differential += tmp_diff;
      if (eval_gnhda) {
        param_it->gnhda += tmp_gnhda;
      }

      // Reset state etc.
      psolver->reset();
      series ++;
      state_.t = time(i + 1);
      state_.x = x0.col(series);
      state_.reset_dparam();
      tmp_diff.zeros();
      tmp_gnhda.zeros();

      // Scales positions with new scale
      scale_positions(series);
    } else {
      // Proceed to next time point
      // Rcpp::Rcout << "exact::evaluate_differential -> Calling step" << std::endl;
      // Rcpp::Rcout << "exact::evaluate_differential -> x: " << state_.x.t() << std::endl;
      // Rcpp::Rcout << "exact::evaluate_differential -> x_dparam: " << state_.x_dparam << std::endl;
      psolver->step(state_, time(i + 1));
      // Rcpp::Rcout << "exact::evaluate_differential -> Called step" << std::endl;
    }

    // Update differential
    if (w_exist) {
      if (state_.is_x0()) {
        tmp_diff(arma::span(0 + series * pfield->d, pfield->d - 1 + series * pfield->d))  += (state_.x_dparam.t() * plikeli->eval_dfit(y.col(i + 1), state_.x, w.col(i + 1)));
      } else {
        tmp_diff += (state_.x_dparam.t() * plikeli->eval_dfit(y.col(i + 1), state_.x, w.col(i + 1)));
      }
      if (eval_gnhda) {
        if (state_.is_x0()) {
          tmp_gnhda(arma::span(0 + series * pfield->d, pfield->d - 1 + series * pfield->d)) += (arma::sum(state_.x_dparam % (plikeli->eval_d2fit(y.col(i + 1), state_.x, w.col(i + 1)) * state_.x_dparam), 0)).t();
        } else {
          tmp_gnhda += (arma::sum(state_.x_dparam % (plikeli->eval_d2fit(y.col(i + 1), state_.x, w.col(i + 1)) * state_.x_dparam), 0)).t();
        }
      }
    } else {
      if (state_.is_x0()) {
        tmp_diff(arma::span(0 + series * pfield->d, pfield->d - 1 + series * pfield->d))  += (state_.x_dparam.t() * plikeli->eval_dfit(y.col(i + 1), state_.x));
      } else {
        tmp_diff  += (state_.x_dparam.t() * plikeli->eval_dfit(y.col(i + 1), state_.x));
      }
      if (eval_gnhda) {
        if (state_.is_x0()) {
          tmp_gnhda(arma::span(0 + series * pfield->d, pfield->d - 1 + series * pfield->d)) += (arma::sum(state_.x_dparam % (plikeli->eval_d2fit(y.col(i + 1), state_.x) * state_.x_dparam), 0)).t();
        } else {
          tmp_gnhda += (arma::sum(state_.x_dparam % (plikeli->eval_d2fit(y.col(i + 1), state_.x) * state_.x_dparam), 0)).t();
        }
      }
    }
  }

  // The last series was not scaled with sc or added to differential
  if (!state_.is_x0()) {
    sc_it = vsc.begin() + (param_it - vparam.begin() - 1);
    if (sc_it->exist) {
      tmp_diff %= sc_it->col(series);
      if (eval_gnhda) {
        tmp_gnhda %= arma::square(sc_it->col(series));
      }
    }
  }
  param_it->differential += tmp_diff;
  if (eval_gnhda) {
    param_it->gnhda += tmp_gnhda;
  }


  // MEAN over time points - n_series (if positive)
  if (n_obs > n_series) {
    param_it->differential  /= n_obs - n_series;
    if (eval_gnhda)  {
      param_it->gnhda /= n_obs - n_series;
    }
  }

  // Return original unscaled ones
  restore_positions();
  // unsigned int nparam = param_it - vparam.begin();
  // Rcpp::Rcout << "exact::evaluate_differential (" << nparam << ")" << param_it->differential.t() << std::endl;

  psolver->reset(); // VERY IMPORTANT (resets, e.g., discretisation parameter, so that the same conditions hold in next loss evaluation)
}

void exact::set_active(std::vector<param>::iterator param_it) {
  set_active_base(param_it);

  // exact specifics: none
}

arma::mat exact::get_X4dfs(std::vector<param>::iterator param_it) {
  // Translate to scale (sentivity_param does not use sc_it if param_it points to x0)
  std::vector<sc>::iterator sc_it = vsc.begin();
  if (&(*param_it) != &(*vparam.begin())) {
    sc_it += param_it - vparam.begin() - 1;
  }

  arma::mat X = psolver->sensitivity_param(pfield, param_it, time, sc_it, exact_diff(param_it - vparam.begin()) == 0);

  if (w_exist) {
    X.each_col() %= arma::vectorise(w);
  }

  return X;
}


//=================================
/* matching subclass
 * Assumes that x is as in R with
 * obs in rows.
 */
//=================================

matching::matching(Rcpp::List ode_struct, arma::mat &A, arma::mat &B, arma::mat x0, Rcpp::Nullable<Rcpp::List> param_, Rcpp::List opt_struct, std::vector<sc>* pvsc_full_, likeli* plikeli_, arma::mat x_) :
  ode_loss(ode_struct, A, B, x0, param_, opt_struct, pvsc_full_, plikeli_),
  x(x_.t()),
  dmatch(std::vector<arma::mat>(vparam.size() - 1)),
  dmatch_full(std::vector<arma::mat>(vparam.size() - 1)) {
  // Rcpp::Rcout << "Starting matching constructor" << std::endl;
  // Rcpp::Rcout << "in matching::matching -> param_.isNotNull() = " << param_.isNotNull() << std::endl;

  // Prepare in times (from x) and out times (y)
  ti  = arma::trans(x.row(0));
  to  = time;
  x   = x.rows(1, pfield->d);

  // Find jump time points in tout and tin
  ti_jumps  = arma::find(arma::diff(ti) < 0);
  arma::uvec vend(1);  vend(0) = ti.n_elem;       // not n_elem - 1, since in tin_sub we subtract one
  ti_jumps = arma::join_cols(ti_jumps + 1, vend); // find(diff(ti)<0) are shifted to the left and 0 and n_elem is always jump point
  vend(0) = 0;  ti_jumps = arma::join_cols(vend, ti_jumps);

  to_jumps = arma::find(arma::diff(to) < 0);
  vend(0) = to.n_elem;
  to_jumps = arma::join_cols(to_jumps + 1, vend);
  vend(0) = 0;  to_jumps = arma::join_cols(vend, to_jumps);

  // Warnings on number of series
  if (ti_jumps.n_elem != to_jumps.n_elem) {
    Rcpp::stop("Time vector in smoothed curves and time vector in y do not have the same number of decreases.");
  }
  for (std::vector<sc>::iterator sc_it = vsc.begin(); sc_it < vsc.end(); sc_it++) {
    if (sc_it->exist && sc_it->value.n_cols != n_series) {
      Rcpp::stop("Number of columns in sc does not number of contexts from data.");
    }
  }

  // Reset everything
  xout.set_size(pfield->d, to.n_elem);
  X0.clear();

  // Identify series, for each of them subset to tin and tout, get and store xout through interpolate,
  for (unsigned int i_series = 0; i_series < n_series; ++i_series) {
    arma::vec ti_sub = ti(arma::span(ti_jumps(i_series), ti_jumps(i_series + 1) - 1));  // time in
    arma::mat xi_sub = x.cols(ti_jumps(i_series), ti_jumps(i_series + 1) - 1);          // x at tin
    arma::vec to_sub = to(arma::span(to_jumps(i_series), to_jumps(i_series + 1) - 1));  // time out

    // xout
    arma::mat xout_sub = interpolate(ti_sub, to_sub, xi_sub);
    xout.cols(to_jumps(i_series), to_jumps(i_series + 1) - 1) = xout_sub;

    // X0
    arma::vec X0_block = xout_sub.col(0);
    X0 = arma::join_cols(X0, X0_block);
  }
  // xout = arma::join_rows(to, xout.t());


  // Make x0 a fixed parameter
  std::vector<param>::iterator it = vparam.begin();
  it->fixed = true;
  it->position_full = X0;
  it->position = X0.elem(it->active);

  // Rcpp::Rcout << "Finished matching constructor." << std::endl;
}

arma::vec matching::interpolate(const double t1, const arma::vec x1, const double t2, const arma::vec x2, const double tm) {
  if (tm <= std::min(t1, t2)) {
    // Left endpoint constant interpolation
    if (t1 <= t2) {
      return x1;
    } else {
      return x2;
    }
  } else if (tm >= std::max(t1, t2)) {
    // Right endpoint constant interpolation
    if (t1 <= t2) {
      return x2;
    } else {
      return x1;
    }
  } else {
    // Linear interpolation
    return (t1 == t2) ? x1 : x1 + (x2 - x1) * (tm - t1) / (t2 - t1);
  }
}
arma::mat matching::interpolate(const arma::vec tin, const arma::vec tout, const arma::mat xin) {
  // assumes tin and tout are non-decreasing!! and length(tin) and ncol(xin) matches
  // returns interpolated value of x at tout, stored in columns (uses linear interpolate, unless oob)
  if (any(arma::diff(tin) < 0)) Rcpp::stop("Internal error: in interpolate, t_in is not non-decreasing!");
  if (any(arma::diff(tout) < 0)) Rcpp::stop("Internal error: in interpolate, t_out is not non-decreasing!");

  arma::mat xout(xin.n_rows, tout.n_elem);
  unsigned int ind = 0, lind = 0, rind = 0;
  for (unsigned int i = 0; i < tout.n_elem; ++i) {
    while (tin(ind) < tout(i) && ind < tin.n_elem - 1) {
      lind = rind;
      ind++;
      rind = ind;
    }
    xout.col(i) = interpolate(tin(lind), xin.col(lind), tin(rind), xin.col(rind), tout(i));
  }
  return xout;
} // returns interpolated value xout at tout using tin and xin (linear interpolation, constant if oob)


/* Set match and dmatch
 * Note: only if linear we use the (permanent) dmatch to get match
 */
void matching::set_match(void) {
  if (arma::all(pfield->linear)) {
    // Rcpp::Rcout << "matching::set_match -> doing linear case" << std::endl;

    // Linear case
    std::vector<arma::mat>::iterator dmatch_it = dmatch.begin();
    std::vector<param>::iterator param_it = vparam.begin() + 1;
    match = (*dmatch_it) * param_it->position;
    for (dmatch_it = dmatch.begin() + 1; dmatch_it < dmatch.end(); dmatch_it++) {
      param_it++;
      match += *dmatch_it * param_it->position;
    }
  } else {
    // Non linear case
    set_nl_match();
  }
}
void matching::set_dmatch(std::vector<param>::iterator param_it) {
  if (&(*param_it) == &(*vparam.begin())) {
    Rcpp::warning("in matching::set_dmatch -> param points to x0, but x0 is fixed.");
  } else {
    if (!arma::all(pfield->linear)) {
      // Convert iterator
      std::vector<arma::mat>::iterator dmatch_it = dmatch.begin();
      dmatch_it += param_it - vparam.begin() - 1;

      // Non linear case
      set_nl_dmatch(dmatch_it, param_it);
    }
    // In linear case we do nothing (it will be set (permanently) in integral/gradient constructor by calling set_nl_dmatch)
  }
}

// With the above we can now define the evaluates from 'fun'
void matching::evaluate(void) {
  // Rcpp::Rcout << "matching::evaluate -> call set_match()" << std::endl;
  set_match();
  // Rcpp::Rcout << "matching::evaluate -> set value" << std::endl;

  if (w_exist) {
    value = plikeli->eval(Y, match, W);
  } else {
    value = plikeli->eval(Y, match);
  }

  if (n_obs > n_series) {
    value /= n_obs - n_series;
  }
}
void matching::evaluate_differential(std::vector<param>::iterator param_it) {
  if (&(*param_it) == &(*vparam.begin())) {
    // x0 is fixed in matching
    param_it->differential.zeros();
  } else {
    // Convert iterator
    std::vector<arma::mat>::iterator dmatch_it = dmatch.begin();
    dmatch_it += param_it - vparam.begin() - 1;

    // We need to get the right match
    set_match();
    set_dmatch(param_it);
    if (w_exist) {
      param_it->differential = dmatch_it->t() * plikeli->eval_dfit(Y, match, W);
      if (eval_gnhda) {
        param_it->gnhda = arma::sum(dmatch_it->t() % (plikeli->eval_d2fit(Y, match, W) * dmatch_it->t()), 0).t();
      }
    } else {
      param_it->differential = dmatch_it->t() * plikeli->eval_dfit(Y, match);
      if (eval_gnhda) {
        param_it->gnhda = arma::sum(dmatch_it->t() % (plikeli->eval_d2fit(Y, match) * dmatch_it->t()), 0).t();
      }
    }

    // Divide by effective no of obs
    if (n_obs > n_series) {
      param_it->differential /= n_obs - n_series;
      if (eval_gnhda)  {
        param_it->gnhda /= n_obs - n_series;
      }
    }
  }
}
void matching::set_active(std::vector<param>::iterator param_it) {
  set_active_base(param_it);

  // x0 does not reduce
  if (&(*param_it) != &(*vparam.begin())) {
    // Matching specifics, only if linear we actually use dmatch_full
    if (arma::all(pfield->linear)) {
      // Convert iterator
      std::vector<arma::mat>::iterator dmatch_it = dmatch.begin();
      dmatch_it += param_it - vparam.begin() - 1;
      std::vector<arma::mat>::iterator dmatch_full_it = dmatch_full.begin();
      dmatch_full_it += param_it - vparam.begin() - 1;

      *dmatch_it = dmatch_full_it->cols(param_it->active);
    }
  }
}

arma::mat matching::get_X4dfs(std::vector<param>::iterator param_it) {
  arma::mat ret;
  if (&(*param_it) == &(*vparam.begin())) {
    Rcpp::warning("matching::get_X for dfs was called on x0 parameter.");
  } else {
    // Convert to dmatch iterator
    std::vector<arma::mat>::iterator dmatch_it = dmatch.begin();
    dmatch_it += param_it - vparam.begin() - 1;

    // Set dmatch
    set_dmatch(param_it);
    ret = *dmatch_it;

    // Scale with W
    if (w_exist) {
      ret.each_col() %= W;
    }
  }
  return ret;
}

Rcpp::List matching::get_XYW(void) {
  Rcpp::List ret;

  Rcpp::List X(dmatch.size());
  std::vector<arma::mat>::iterator dmatch_it = dmatch.begin();
  for (unsigned int i = 0; i < dmatch.size(); i++) {
    // Rcpp::Rcout << "matching::get_XYW -> i = " << i << " dmatch:" << dmatch_it->t() << std::endl;
    arma::mat dm = *dmatch_it;
    X[i] = dm;
    dmatch_it++;
  }
  // Rcpp::Rcout << "matching::get_XYW -> finished X loop" << std::endl;

  ret["X"] = X;
  ret["Y"] = Y;
  ret["W"] = W;
  ret["X0"] = X0;
  ret["xout"] = arma::join_rows(to, xout.t());
  return ret;
}



//=================================
// Integral (matching) subclass
//=================================

integral::integral(Rcpp::List ode_struct, arma::mat &A, arma::mat &B, arma::mat x0, Rcpp::Nullable<Rcpp::List> param_, Rcpp::List opt_struct, std::vector<sc>* pvsc_full_, likeli* plikeli_, arma::mat x_) :
  matching(ode_struct, A, B, x0, param_, opt_struct, pvsc_full_, plikeli_, x_) {
  // Rcpp::Rcout << "Starting integral constructor" << std::endl;
  // Rcpp::Rcout << "in integral::integral -> param_.isNotNull() = " << param_.isNotNull() << std::endl;

  // Set Y and W
  Y.clear();
  W.clear();

  // Identify series, for each of them subset to tin and tout, get and store xout through interpolate,
  for (unsigned int i_series = 0; i_series < n_series; ++i_series) {
    // Rcpp::Rcout << "  Extract xout" << std::endl;

    // Extract xout part for that series
    arma::mat xout_sub = xout.cols(to_jumps(i_series), to_jumps(i_series + 1) - 1);
    // Rcpp::Rcout << "  Fill in Y" << std::endl;

    // Y
    // arma::mat Y_block = xout_sub.tail_cols(xout_sub.n_cols - 1) - xout_sub.head_cols(xout_sub.n_cols - 1);
    Y = arma::join_cols(Y, arma::vectorise(xout_sub.tail_cols(xout_sub.n_cols - 1) - xout_sub.head_cols(xout_sub.n_cols - 1)));
    // Rcpp::Rcout << "  Fill in W" << std::endl;

    // W
    if (w_exist) {
      arma::mat W_block = w(arma::span::all, arma::span(to_jumps(i_series), to_jumps(i_series + 1) - 1));
      W_block = (W_block.tail_cols(W_block.n_cols - 1) + W_block.head_cols(W_block.n_cols - 1)) / 2.0;
      W = arma::join_cols(W, arma::vectorise(W_block));
    }
  }
  // Rcpp::Rcout << "  Initialise dmatch_full if it is linear" << std::endl;

  // Initialise dmatch_full if it is linear
  if (arma::all(pfield->linear)) {
    std::vector<param>::iterator param_it = vparam.begin() + 1;
    for (std::vector<arma::mat>::iterator dmatch_full_it = dmatch_full.begin(); dmatch_full_it < dmatch_full.end(); dmatch_full_it++) {
      set_nl_dmatch(dmatch_full_it, param_it);

      // Load it onto dmatch
      set_active(param_it);
      // Rcpp::Rcout << "  dmatch_full_it:" << dmatch_full_it->t() << std::endl;
      // Rcpp::Rcout << "  dmatch_full:" << dmatch_full.begin()->t() << std::endl;
      // Rcpp::Rcout << "  dmatch:" << dmatch.begin()->t() << std::endl;

      param_it++;
    }
  }

  // Update value
  evaluate();
  // Rcpp::Rcout << "Finished integral constructor" << std::endl;

}

// Simpson rules
void integral::simpson(arma::vec &x1, arma::mat &y1, arma::vec x2, arma::mat &y2, std::vector<param>::iterator param_it) {
  // Does simpson approx of int(t1, t2, f_dparam(x(s), param)) (unless param_it->x0 then f(x(s), param))
  y1  = y2; // y2 is the new y1 (the new y2 comes from x2), x1 is old x2
  // Rcpp::Rcout << "integral::simpson -> it number " << param_it - vparam.begin() << std::endl;

  if (&(*param_it) == &(*vparam.begin())) {
    // First parameter (x0), means we integrate f(x(s))
    y2  = pfield->f(0.0, x2);
    y1 += y2 + 4 * pfield->f(0.0, (x1 + x2) / 2);
  } else {
    y2  = pfield->f_dparam(0.0, x2, param_it);
    y1 += y2 + 4 * pfield->f_dparam(0.0, (x1 + x2) / 2, param_it);
  }
  y1 /= 6;
  x1  = x2;
} // int(t1, t2, f_dparam(x(s))) claimed by in y1 * (t2 - t1) (and fx1 * (t2 - t1))

void integral::simpson(const double t1, arma::vec &x1, arma::mat &y1,
  const double t2, arma::vec x2, arma::mat &y2, arma::vec tm, arma::mat xm, std::vector<param>::iterator param_it) {
  // follows simpson convention: y2 becomes the new y1, x1 is the old x2
  arma::mat integral_ = y1;    // When passed from last simpson_mid call

  if (!tm.is_empty()) {
    simpson(x1, y1, xm.col(0), y2, param_it);  // First step: from x1 to first xm
    integral_ += (tm(0) - t1) * y1;
    for (unsigned int i = 1; i < tm.n_elem; ++i) {
      simpson(x1, y1, xm.col(i), y2, param_it);
      integral_ += (tm(i) - tm(i - 1)) * y1;
    }
  }
  simpson(x1, y1, x2, y2, param_it);  // Last step from last xm (stored in x1) to x2, called no matter if tm is empty
  if (!tm.is_empty()) {
    integral_ += (t2 - tm(tm.n_elem - 1)) * y1;
  } else {
    integral_ += (t2 - t1) * y1;
  }
  y1 = integral_; // For return
} // int(t1, t2, f(x(s))) stored in y1 based on the intermediates xm at tm. If no intermediates, it just calls simpson.

arma::mat integral::simpson(arma::vec tin, arma::vec tout, arma::mat xin, arma::mat xout, std::vector<param>::iterator param_it) {
  // Assumes tin and tout are non-decreasing and represents the whole series, xin matches tin and xout tout

  // For filling in stacked simpson integrals as dmatch
  arma::mat X;
  if (&(*param_it) == &(*vparam.begin())) {
    X.zeros(pfield->d * (tout.n_elem - 1), 1);
  } else {
    X.zeros(pfield->d * (tout.n_elem - 1), param_it->p);
  }

  // Initialise
  arma::mat y1(pfield->d, X.n_cols, arma::fill::zeros);
  arma::vec x1 = xout.col(0);
  arma::mat y2;
  if (&(*param_it) == &(*vparam.begin())) {
    y2 = pfield->f(0.0, x1);
  } else {
    y2 = pfield->f_dparam(0.0, x1, param_it);
  }

  for (unsigned int j = 0; j < tout.n_elem - 1; ++j) {
    // Find tin and xin intermediates
    arma::uvec rhs  = (tout(j) < tin);
    arma::uvec lhs  = (tin < tout(j + 1));
    arma::uvec mids = arma::find(rhs && lhs);

    // Integrate (result sent to y1)
    simpson(tout(j), x1, y1, tout(j + 1), xout.col(j + 1), y2, tin.elem(mids), xin.cols(mids), param_it);

    // Store y1 in X
    X.rows(pfield->d * j, pfield->d * (j + 1) - 1) = y1;

    // Reset y1 (if not this, then {int(tout(0), tout(i), f(x(s)))}_i for i = 0,...,end, stacked on top of each other is returned)
    y1.zeros();
  }
  return X;
} // {int(tout(i), tout(i + 1), f(x(s)))}_i for i = 0,...,end-1, stacked on top of each other if param_it->x0, else its f_dparam()


void integral::set_nl_match(void) {
  // The effective parameter is eta = schur.prod(param, sc) (but param is penalised)
  store_positions();
  scale_positions(0);
  // Rcpp::Rcout << "integral::set_nl_match -> scaled" << std::endl;

  // We use simpson integration to get int^t2_t1(f(x(s)) ds), but point to first parameter so that it use f, not f_dparam
  match.clear();

  // Identify series, for each of them subset to tin and tout, get and store xout, call simpson on it
  for (unsigned int i_series = 0; i_series < n_series; ++i_series) {
    scale_positions(i_series);

    // Sub times and trajectories
    arma::vec tin_sub = ti(arma::span(ti_jumps(i_series), ti_jumps(i_series + 1) - 1));
    arma::mat xin_sub = x.cols(ti_jumps(i_series), ti_jumps(i_series + 1) - 1);
    arma::vec tout_sub = to(arma::span(to_jumps(i_series), to_jumps(i_series + 1) - 1));
    // Rcpp::Rcout << "integral::set_nl_match -> call simpson" << std::endl;

    // fX (series part of match)
    arma::vec fX_block  = simpson(tin_sub, tout_sub, xin_sub, xout.cols(to_jumps(i_series), to_jumps(i_series + 1) - 1), vparam.begin());
    // Rcpp::Rcout << "integral::set_nl_match -> called simpson" << std::endl;

    //fX_block.each_row() %= sc_sub.col(i_series).t();
    match = arma::join_cols(match, fX_block);
  }

  // Return original unscaled ones
  restore_positions();
}

void integral::set_nl_dmatch(std::vector<arma::mat>::iterator dmatch_it, std::vector<param>::iterator param_it) {
  if (&(*param_it) == &(*vparam.begin())) {
    Rcpp::warning("error in integral::set_nl_dmatch -> x0 is not optimised in matching, param_it points to x0.");
  } else {
    // Iterators
    std::vector<sc>::iterator sc_it = vsc.begin();
    sc_it += param_it - vparam.begin() - 1;

    // The effective parameter is eta = schur.prod(param, sc) (but param is penalised)
    store_positions();
    scale_positions(0);

    // We use simpson integration to get int^t2_t1(f_dparam(x(s)) ds),
    dmatch_it->clear();
    arma::mat X(0, param_it->p);

    // Identify series, for each of them subset to tin and tout, get and store xout, call simpson on it
    for (unsigned int i_series = 0; i_series < n_series; ++i_series) {
      // Scale params
      scale_positions(i_series);

      // Sub times and trajectories
      arma::vec tin_sub = ti(arma::span(ti_jumps(i_series), ti_jumps(i_series + 1) - 1));
      arma::mat xin_sub = x.cols(ti_jumps(i_series), ti_jumps(i_series + 1) - 1);
      arma::vec tout_sub = to(arma::span(to_jumps(i_series), to_jumps(i_series + 1) - 1));

      // X (series part of dmatch)
      arma::mat X_block  = simpson(tin_sub, tout_sub, xin_sub, xout.cols(to_jumps(i_series), to_jumps(i_series + 1) - 1), param_it);
      if (sc_it->exist) {
        X_block.each_row() %= sc_it->col(i_series).t();
      }
      X = arma::join_cols(X, X_block);
    }

    // Return original unscaled ones
    restore_positions();

    *dmatch_it = X;
    // Rcpp::Rcout << "integral::set_nl_dmatch --> X " << X << std::endl;

  }
}


//=================================
// Design matrix et al for aim
//=================================

// Due to Romain
Rcpp::List resize(const Rcpp::List& x, int n){
  int oldsize = x.size();
  Rcpp::List y(n);
  for(int i = 0; i < oldsize; i++) y[i] = x[i];
  return y;
}

// [[Rcpp::export(aim_design)]]
Rcpp::List aim_design(Rcpp::List ode_struct, Rcpp::List opt_struct, Rcpp::List sc_,
  arma::mat x, Rcpp::Nullable<Rcpp::List> param_) {

  // Prepare likeli
  gaussian gg;
  likeli* plikeli_ = &gg;

  // Prepare sc
  std::vector<sc> vsc_full_;
  for (unsigned int i = 0; i < sc_.size(); i++) {
    vsc_full_.push_back(sc(Rcpp::as< Rcpp::Nullable<arma::mat> >(sc_[i])));
  }

  // Create ode_loss
  arma::mat A;
  arma::mat B;

  // Integral matching object
  arma::mat y_ = Rcpp::as<arma::mat>(opt_struct["y"]);
  unsigned int n_series = arma::sum(arma::diff(y_.col(0)) < 0) + 1;
  arma::mat x0 = x.head_rows(n_series);
  x0 = x0.cols(1, x.n_cols - 1);
  // Note x.head_rows() should represent x0, but in matching constructor it is overloaded with X0 (made from xout)
  integral integral_(ode_struct, A, B, x0.t(), param_, opt_struct, &vsc_full_, plikeli_, x);

  // If not linear we fill in dmatch
  if (!arma::all(integral_.pfield->linear)) {
    std::vector<param>::iterator param_it;
    for (param_it = integral_.vparam.begin() + 1; param_it < integral_.vparam.end(); param_it++) {
      integral_.set_dmatch(param_it);
    }
  }

  return integral_.get_XYW();
}



// //' Test loss stuff
// //'
// //' This function returns a logical vector identifying if
// //' there are leading NA, marking the leadings NA as TRUE and
// //' everything else as FALSE.
// //'
// //' @param x An integer vector
// //' @export
// // [[Rcpp::export]]
// Rcpp::List test_loss(Rcpp::List ode_struct, Rcpp::List opt_struct,
//   arma::mat x0, Rcpp::Nullable<Rcpp::List> param_, Rcpp::List sc_, arma::mat x) {
//   Rcpp::List ret;
//
//   // Prepare likeli
//   gaussian gg;
//   likeli* plikeli_ = &gg;
//   // Rcpp::Rcout << "Likeli test: " << plikeli_->eval(x.col(0), arma::ones(x.n_rows)) << ". Prepare sc" << std::endl;
//
//   // Prepare sc
//   std::vector<sc> vsc_full_;
//   for (unsigned int i = 0; i < sc_.size(); i++) {
//     vsc_full_.push_back(sc(Rcpp::as< Rcpp::Nullable<arma::mat> >(sc_[i])));
//   }
//   // Rcpp::Rcout << "Prepared sc" << std::endl;
//
//   // Create exact ode_loss (also creates psolver (via exact) and pfield (via super class ode_loss))
//   arma::mat A;
//   arma::mat B;
//
//   // integral integral_ (ode_struct, A, B, x0, param_, opt_struct, &vsc_full_, plikeli_, x);
//   exact exact_ (ode_struct, A, B, x0, param_, opt_struct, &vsc_full_, plikeli_);
//
//   // Integral matching test
//
//   // Try it out
//   fun* pfun = &exact_;
//   // Rcpp::Rcout << "Prepared exact" << std::endl;
//
//   // Try it out
//   ret["l"] = pfun->value;
//   std::vector<param>::iterator param_it = pfun->vparam.begin();
//   ret["l_dx0"] = param_it->differential;
//   param_it++;
//   ret["l_dparam"] = param_it->differential;
//
//
//   pfun->evaluate();
//
//   ret["l2"] = pfun->value;
//   param_it = pfun->vparam.begin();
//   ret["l2_dx0"] = param_it->differential;
//   param_it++;
//   ret["l2_dparam"] = param_it->differential;
//
//
//   pfun->evaluate_differential(param_it);
//
//   ret["l3"] = pfun->value;
//   param_it = pfun->vparam.begin();
//   ret["l3_dx0"] = param_it->differential;
//   param_it++;
//   ret["l3_dparam"] = param_it->differential;
//
//
//   param_it = pfun->vparam.begin();
//   pfun->evaluate_differential(param_it);
//
//   ret["l4"] = pfun->value;
//   param_it = pfun->vparam.begin();
//   ret["l4_dx0"] = param_it->differential;
//   param_it++;
//   ret["l4_dparam"] = param_it->differential;
//
//
//   param_it->position += 0.5;
//   pfun->evaluate_differential(param_it);
//
//   ret["l5"] = pfun->value;
//   param_it = pfun->vparam.begin();
//   ret["l5_dx0"] = param_it->differential;
//   param_it++;
//   ret["l5_dparam"] = param_it->differential;
//   // Rcpp::Rcout << "Finished exact loss" << std::endl;
//
//
//
//   integral integral_ (ode_struct, A, B, x0, param_, opt_struct, &vsc_full_, plikeli_, x);
//   integral_.vparam.begin()->fixed = true;
//   pfun = &integral_;
//   // Rcpp::Rcout << "Created integral" << std::endl;
//
//
//   ret["i"] = pfun->value;
//   param_it = pfun->vparam.begin();
//   ret["i_dx0"] = param_it->differential;
//   param_it++;
//   ret["i_dparam"] = param_it->differential;
//   // Rcpp::Rcout << "Evaluate fun" << std::endl;
//
//
//   pfun->evaluate();
//   // Rcpp::Rcout << "Evaluated fun" << std::endl;
//
//
//
//   ret["i2"] = pfun->value;
//   param_it = pfun->vparam.begin();
//   ret["i2_dx0"] = param_it->differential;
//   param_it++;
//   ret["i2_dparam"] = param_it->differential;
//   // Rcpp::Rcout << "Evaluate differential 2" << std::endl;
//
//
//   pfun->evaluate_differential(param_it);
//   // Rcpp::Rcout << "Evaluated differential 2" << std::endl;
//
//
//   ret["i3"] = pfun->value;
//   param_it = pfun->vparam.begin();
//   ret["i3_dx0"] = param_it->differential;
//   param_it++;
//   ret["i3_dparam"] = param_it->differential;
//   // Rcpp::Rcout << "Evaluate differential 1" << std::endl;
//
//
//   param_it = pfun->vparam.begin();
//   pfun->evaluate_differential(param_it);
//   // Rcpp::Rcout << "Evaluated differential 1" << std::endl;
//
//
//   ret["i4"] = pfun->value;
//   param_it = pfun->vparam.begin();
//   ret["i4_dx0"] = param_it->differential;
//   param_it++;
//   ret["i4_dparam"] = param_it->differential;
//   // Rcpp::Rcout << "Evaluate differential 2 at new param" << std::endl;
//
//   param_it->position += 0.5;
//   pfun->evaluate_differential(param_it);
//   // Rcpp::Rcout << "Evaluated differential 2 at new param" << std::endl;
//
//   ret["i5"] = pfun->value;
//   param_it = pfun->vparam.begin();
//   ret["i5_dx0"] = param_it->differential;
//   param_it++;
//   ret["i5_dparam"] = param_it->differential;
//   // Rcpp::Rcout << "Get XYW" << std::endl;
//
//   Rcpp::List XYW = integral_.get_XYW();
//   // Rcpp::Rcout << "Got XYW" << std::endl;
//
//   ret["XYW"] = XYW;
//
//
//   // // Scales
//   // arma::mat sc_ = sc.isNotNull() ? Rcpp::as<arma::mat>(sc) : arma::ones<arma::mat>(pode->p, x0.n_cols);
//   // arma::mat* psc = &sc_;
//   //
//   // exact exact_ (y, w, psc, plikeli, pode, psolver);
//   // loss* ploss = &exact_;
//   //
//   // ploss->update(param, x0, true);
//   //
//   // ret["l"] = ploss->l;
//   // ret["l_dparam"] = ploss->l_dparam;
//   // ret["l_dx0"] = ploss->l_dx0;
//
//   // x = x.t();
//   // im im_ (y, w, psc, plikeli, pode, psolver, x, param);
//   // ploss = &im_;
//   //
//   // ploss->update(param, x0, true);
//   //
//   // ret["l2"] = ploss->l;
//   // ret["l_dparam2"] = ploss->l_dparam;
//   // ret["l_dx02"] = ploss->l_dx0;
//   //
//   // ret["design"] = im_.get_XYW();
//   //
//   // ret["has_mak"] = is_in("mak", ode_struct.attr("class"));
//   // ret["has_plk"] = is_in("plk", ode_struct.attr("class"));
//   // ret["has_ode"] = is_in("ode", ode_struct.attr("class"));
//
//
//   return ret;
// }

