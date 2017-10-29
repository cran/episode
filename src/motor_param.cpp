//=============================================
/* motor_param.cpp
*
* Content:
* - get_span
* - param definitions
*/
//=============================================

#include "motor_param.h"

// Returns uvector holding all indices of v
arma::uvec get_span(const arma::vec &v) {
  return arma::regspace<arma::uvec>(0, v.n_elem - 1);
}


//=================================
// param class
//=================================

void param::set_coord(void) {
  if (nrow >= 1) {
    arma::uvec div = active / nrow; // Division --> column in theta, note: it remains unchanged!
    arma::uvec rem = active - arma::floor(active / nrow) * nrow;   // Remainder --> row in theta
    // The next three steps removes the gaps (corresponding to zero columns)
    div = div.n_elem > 0 ? join_cols(arma::zeros<arma::uvec>(1), arma::diff(div)) : div; // This is the change in column-number
    div = arma::clamp(div, 0, 1);
    div = arma::cumsum(div);
    coord = div * nrow + rem;
    div = arma::unique(div);
    r_sub = div.n_elem;
  } else {
    coord = arma::regspace<arma::uvec>(0, p - 1);
    r_sub = p;
  }
}

arma::mat param::get_dparam(const arma::vec v) {
  /* Assuming that the expression is f(position) = get_position() * v,
   * this will return df/dposition as matrix, but only the positions from param->coord (same length as active)
   */
  arma::mat ret = arma::zeros<arma::mat>(nrow, coord.n_elem);
  arma::uvec mat_coor = arma::regspace<arma::uvec>(0, coord.n_elem - 1) * nrow + (coord - (coord / nrow) * nrow);
  ret.elem(mat_coor) = v.elem(coord / nrow);
  return ret;
}

// Constructor (no copy constructor or destructor defined, as param contains no raw pointers)
// nrow_ = 0 means param is thought of as vector
param::param (arma::vec position_, ctrl c_, reg* preg_, bool fixed_, bool screen_, double lambda_factor_ = 1, unsigned int nrow_ = 0) :
  position_full(position_),
  active(get_span(position_full)),
  position(position_),
  differential(arma::zeros(arma::size(position))),
  gnhda(arma::zeros(arma::size(position))),
  p_full(position_full.n_elem),
  p(position.n_elem),
  fixed(fixed_),
  screen(screen_),
  c(c_),
  preg(preg_->shallow_Clone()),
  lambda_factor(lambda_factor_),
  nrow(nrow_) {
  // Coord and r_sub
  set_coord();

  // Check reg has right dimensions
  if (preg->lower.n_elem != position_.n_elem) Rcpp::stop("Length of lower bound on parameter does not match length of parameter.");
  if (preg->upper.n_elem != position_.n_elem) Rcpp::stop("Length of upper bound on parameter does not match length of parameter.");
  if (preg->v.n_elem != position_.n_elem) Rcpp::stop("Length of penalty weights does not match length of parameter.");

  // Box in the position
  set_lambda(preg->lambda);
  preg->box(position_full);
  preg->box(position);
  pen = preg->penalty(position_full);
}

param::param (Rcpp::List param_) :
  position_full(Rcpp::as<arma::vec>(param_["x"])),
  active(get_span(position_full)),
  position(position_full),
  differential(arma::zeros(arma::size(position))),
  gnhda(arma::zeros(arma::size(position))),
  p_full(position_full.n_elem),
  p(position.n_elem),
  fixed(Rcpp::as<bool>(param_["fixed"])),
  screen(Rcpp::as<bool>(param_["screen"])),
  c(ctrl(Rcpp::as<unsigned int>(param_["step_cycle"]),
    Rcpp::as<unsigned int>(param_["step_max"]),
    Rcpp::as<unsigned int>(param_["step_screen"]),
    Rcpp::as<unsigned int>(param_["backtrack_max"]),
    Rcpp::as<double>(param_["tau_init"]),
    Rcpp::as<double>(param_["tau_min"]),
    Rcpp::as<double>(param_["tau_scale"]),
    Rcpp::as<double>(param_["tol"]))),
  lambda_factor(1),
  nrow(Rcpp::as<unsigned int>(param_["nrow"])) {
  // coord and r_sub
  set_coord();

  // Regularisation pointer
  Rcpp::Nullable<arma::vec> v_ = Rcpp::as< Rcpp::Nullable<arma::vec> >(param_["v"]);
  arma::vec v = v_.isNotNull() ? Rcpp::as<arma::vec>(v_) : arma::ones(position_full.n_elem);

  Rcpp::Nullable<arma::vec> lower_ = Rcpp::as< Rcpp::Nullable<arma::vec> >(param_["lower"]);
  arma::vec lower(position_full.n_elem);
  if (lower_.isNotNull()) {
    lower = stretch_lu(Rcpp::as<arma::vec>(lower_), position_full.n_elem);
  } else {
    lower.fill(-arma::datum::inf);
  }

  Rcpp::Nullable<arma::vec> upper_ = Rcpp::as< Rcpp::Nullable<arma::vec> >(param_["upper"]);
  arma::vec upper(position_full.n_elem);
  if (upper_.isNotNull()) {
    upper = stretch_lu(Rcpp::as<arma::vec>(upper_), position_full.n_elem);
  } else {
    upper.fill(arma::datum::inf);
  }

  std::string reg_type = Rcpp::as<std::string>(param_["reg_type"]);

  // Check reg has right dimensions
  if (lower.n_elem != position_full.n_elem) Rcpp::stop("Length of lower bound on parameter does not match length of parameter.");
  if (upper.n_elem != position_full.n_elem) Rcpp::stop("Length of upper bound on parameter does not match length of parameter.");
  if (v.n_elem != position_full.n_elem) Rcpp::stop("Length of penalty weights does not match length of parameter.");

  arma::vec lambda = Rcpp::as<arma::vec>(param_["lambda"]);
  preg = create_preg(reg_type, lower, upper, v, lambda(0) * lambda_factor, Rcpp::as<double>(param_["a"]));

  // Box in the position
  preg->box(position_full);
  preg->box(position);
  pen = preg->penalty(position_full);
}

param::param (Rcpp::List reg_struct, arma::vec x, unsigned int nrow_) :
  position_full(x),
  active(get_span(position_full)),
  position(position_full),
  differential(arma::zeros(arma::size(position))),
  gnhda(arma::zeros(arma::size(position))),
  p_full(position_full.n_elem),
  p(position.n_elem),
  fixed(Rcpp::as<bool>(reg_struct["fixed"])),
  screen(Rcpp::as<bool>(reg_struct["screen"])),
  lambda_factor(Rcpp::as<double>(reg_struct["lambda_factor"])),
  nrow(nrow_) {

  Rcpp::List ctrl_struct = Rcpp::as<Rcpp::List>(reg_struct["ctrl"]);
  c = ctrl(Rcpp::as<unsigned int>(ctrl_struct["step_cycle"]),
    Rcpp::as<unsigned int>(ctrl_struct["step_max"]),
    Rcpp::as<unsigned int>(ctrl_struct["step_screen"]),
    Rcpp::as<unsigned int>(ctrl_struct["backtrack_max"]),
    Rcpp::as<double>(ctrl_struct["tau_init"]),
    Rcpp::as<double>(ctrl_struct["tau_min"]),
    Rcpp::as<double>(ctrl_struct["tau_scale"]),
    Rcpp::as<double>(ctrl_struct["tol"]));

  // coord and r_sub
  set_coord();

  // Regularisation pointer
  preg = create_preg(reg_struct, p_full);

  // Box in the position
  preg->box(position_full);
  preg->box(position);
  pen = preg->penalty(position_full);
}

param::~param() {
}

void param::set_lambda(double lambda_) {
  preg->lambda = lambda_ * lambda_factor;
}

// Apply prox operator in certain direction, used in optim::backtrack
arma::vec param::proximal(arma::vec position_, arma::vec direction) {
  return preg->prox(position_, direction, c.tau);
}

/* set_active which sets the new active set.
 * Since 'direction' in optim is what matters, differential, gnhda is not preserved.
 */
void param::set_active(arma::uvec active_) {
  // Store current value in full
  position_full.elem(active) = position;

  // Transfer value at new active
  active = active_;
  position = position_full.elem(active);
  differential.zeros(position.n_elem);
  gnhda.zeros(position.n_elem),

  // Sizes
  p = position.n_elem;

  // Coord and r_sub
  set_coord();

  // Handle reg pointer
  preg->reduce(active);
}

// if void argument, then it means the full restore. Differential, gnhda not preserved.
void param::reset_active(void) {
  // Store current value in full
  position_full.elem(active) = position;

  // Transfer value at new active
  active = get_span(position_full);
  position = position_full;
  differential.zeros(position.n_elem);
  gnhda.zeros(position.n_elem),

  // Sizes
  p = position.n_elem;

  // Coord
  set_coord();

  // Handle reg pointer
  preg->reduce();
}

/* Stein DF, assumes that ncol(X) (linearisation) and length of param->position matches (i.e., call reduce before calling this).
 * NOTE:  if weights exists, make sure to
 *            'X.each_col() %= arma::vectorise(w);'
 *        before calling this!
 */
double param::dfs(arma::mat X) {
  // Only the positive parameters
  arma::uvec upos  = arma::find(position);
  X = X.cols(upos);
  X.elem( arma::find_nonfinite(X) ).zeros();
  // Rcpp::Rcout << "param::dfs --> X.cols(upos) " << X << std::endl;

  // d²pen (multiplied by v)
  arma::uvec old_active = active;
  preg->reduce(upos);
  arma::vec va    = preg->d2pen(position.elem(upos));
  preg->reduce(old_active);
  // Rcpp::Rcout << "param::dfs --> va " << va.t() << std::endl;

  // Squared singular values
  arma::mat X_1   = X.cols(arma::find(va <= 0));
  arma::mat X_2   = X.cols(arma::find(va > 0));
  X_2.each_row() /= arma::sqrt(va.elem(arma::find(va > 0)).t());
  // Rcpp::Rcout << "param::dfs --> X_1 " << X_1 << std::endl;
  // Rcpp::Rcout << "param::dfs --> X_2 " << X_2 << std::endl;

  arma::vec d     = arma::svd(X_2);
  d %= d;

  // Add rank difference from full matrix to only those with positive d^2pen
  double ret = d.n_elem > 0 ? arma::sum(d / (d + preg->lambda)) : 0.0;
  return ret + (arma::rank(X) - arma::rank(X_2));  // Add additional parameters corresponding to negatively penalised
}
/*
 * Steins df = div(y -> hat(y)) and if
 * l(y, b) ~= .5 |y - xb|^2_2+lambda*pen(b)
 * then dy(b) = -X(X'X + V)^-1, where V = diag(lambda d^2pen), where the X with b=0 has been discarded
 * So X = (X_1, X_2) where X_1 has d²pen = 0 and X_2 has d2pen != 0
 * then dfs  = rank(X_2) + tr{(X_2'(I - Pi_1)X_2 + V_1)^-1(X_2'(I - Pi_1)X_2)}
 * set X_3 = (I - Pi_1)X_2 then
 * dfs = rank(X_2) + tr{(X_3'X_3 + V_1)^-1(X_3'X_3)}
 * X_4 = X_3 * sqrt(V_1)
 * then
 * dfs = rank(X_2) + sum_{d singular of X_4}(d_j^2 / (d_j^2 + 1))
 */



// Returns position as column-wise matrix with nrow rows, unless nrow == 0, then as a col-vector
arma::mat param::get_position(void) {
  if (nrow >= 1) {
    arma::mat theta = arma::zeros(nrow, r_sub);
    theta.elem(coord) = position;
    return theta;
  } else {
    return position;
  }
}
