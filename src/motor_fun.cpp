//=============================================
/* motor_fun.cpp
 *
 * Content:
 * - fun and subclasses definitions
 *
 */
//=============================================

#include "motor_fun.h"



//=================================
// fun class
//=================================

// Constructor
fun::fun(std::vector<param> vparam_) :
  lambda(0.0),
  vparam(vparam_),
  eval_gnhda(false),
  value(0.0){
}

//
fun::fun(const fun& other) :
  lambda(other.lambda),
  vparam(other.vparam),
  eval_gnhda(other.eval_gnhda),
  value(other.value) {

}

double fun::get_lambda(void) {
  return lambda;
}
void fun::set_lambda(double lambda_) {
  lambda = lambda_;
  for (std::vector<param>::iterator param_it = vparam.begin(); param_it < vparam.end(); param_it++) {
    param_it->set_lambda(lambda);
  }
}

void fun::evaluate(std::vector<param>::iterator it, bool differential) {
  // Set fun value
  evaluate();

  // Add penalty and optionally find differential
  if (it != vparam.end()) {
    value += it->preg->penalty(it->position);

    if (differential) {
      evaluate_differential(it);
    }
  } else {
    Rcpp::stop("error in fun::evaluate -> iterator out of bounds.");
  }
}

std::vector<arma::vec> fun::diff_at_zero(void) {
  // For storing current values and return stuff
  std::vector<arma::vec> new_diff(vparam.size());
  std::vector<arma::uvec> old_active(vparam.size());
  std::vector<arma::vec> old_diff(vparam.size());
  std::vector<arma::vec> old_gnhda(vparam.size());
  std::vector<arma::vec> old_pos_full(vparam.size());
  std::vector<arma::vec>::iterator new_diff_it = new_diff.begin(), old_diff_it = old_diff.begin(), old_gnhda_it = old_gnhda.begin(), old_pos_full_it = old_pos_full.begin();
  std::vector<arma::uvec>::iterator old_active_it = old_active.begin();

  // Set size of the above vectors
  std::vector<param>::iterator param_it = vparam.begin();
  for (param_it = vparam.begin(); param_it < vparam.end(); param_it++) {
    *new_diff_it      = arma::zeros(param_it->p_full);
    *old_active_it    = arma::zeros<arma::uvec>(param_it->active.n_elem);
    *old_diff_it      = arma::zeros(param_it->differential.n_elem);
    *old_gnhda_it     = arma::zeros(param_it->gnhda.n_elem);
    *old_pos_full_it  = arma::zeros(param_it->position_full.n_elem);

    new_diff_it++;  old_active_it++;  old_diff_it++;  old_gnhda_it++; old_pos_full_it++;
  }


  // Split into three loops (necessary, since all positions must be 0 first)
  // First loop: store old stuff and set current position to 0
  // Rcpp::Rcout << "fun::diff_at_zero -> First loop" << std::endl;
  old_diff_it = old_diff.begin(), old_gnhda_it = old_gnhda.begin(), old_pos_full_it = old_pos_full.begin();
  old_active_it = old_active.begin();
  for (param_it = vparam.begin(); param_it < vparam.end(); param_it++) {
    if (!param_it->fixed) {
      // Save old active differential etc
      *old_active_it  = param_it->active;
      *old_diff_it    = param_it->differential;
      *old_gnhda_it   = param_it->gnhda;

      // Reset active, so that all info in position is stored in position full
      param_it->reset_active();
      *old_pos_full_it  = param_it->position_full;

      // Set to zero
      param_it->position.zeros();

      // Adapt fun to new (full) active set
      set_active(param_it);
    }
    old_active_it++;  old_diff_it++;  old_gnhda_it++; old_pos_full_it++;
  }


  // Second loop: gather differentials at 0,
  new_diff_it = new_diff.begin();
  for (param_it = vparam.begin(); param_it < vparam.end(); param_it++) {
    if (!param_it->fixed) {
      // Evaluate and save differential at 0
      evaluate_differential(param_it);
      *new_diff_it  = param_it->differential;
    }
    new_diff_it++;
  }


  // Third loop: restore everything
  old_diff_it = old_diff.begin(), old_gnhda_it = old_gnhda.begin(), old_pos_full_it = old_pos_full.begin();
  old_active_it = old_active.begin();
  for (param_it = vparam.begin(); param_it < vparam.end(); param_it++) {
    if (!param_it->fixed) {
      // Restore everything
      param_it->position_full = *old_pos_full_it;
      param_it->position      = old_pos_full_it->elem(*old_active_it);
      param_it->set_active(*old_active_it);
      param_it->differential  = *old_diff_it;
      param_it->gnhda         = *old_gnhda_it;

      // Adapt fun to newly restored old active set
      set_active(param_it);
    }
    old_active_it++;  old_diff_it++;  old_gnhda_it++; old_pos_full_it++;
  }


  return new_diff;
}

arma::vec fun::lambda_max(void) {
  // Return object
  arma::vec lam_max = arma::zeros(vparam.size());

  // Get diff at 0
  std::vector<arma::vec> diffs0 = diff_at_zero();

  // Get the maximal entries
  std::vector<param>::iterator param_it = vparam.begin();
  std::vector<arma::vec>::iterator diff_it = diffs0.begin();
  for (param_it = vparam.begin(); param_it < vparam.end(); param_it++) {
    if (!param_it->fixed) {
      // Use only finites
      arma::uvec finites = arma::find_finite(*diff_it);
      arma::vec l = param_it->preg->lower_full.elem(finites);
      arma::vec u = param_it->preg->upper_full.elem(finites);
      arma::vec di = diff_it->elem(finites);

      // Useful only if bounds allow for differential to move
      arma::uvec lhs = (di > 0); arma::uvec rhs = (l < 0);
      arma::uvec useful = lhs & rhs;
      lhs = (di < 0); rhs = (u > 0);
      useful = arma::find(useful | (lhs & rhs));
      // useful = arma::find(((di > 0) & (l < 0)) | ((di < 0) & (u > 0))); //

      di = di.elem(useful);
      // Rcpp::Rcout << "fun::lambda_max -> di" << di.t() << std::endl <<
      //   "param_it->preg->v_exist = " << param_it->preg->v_exist << std::endl <<
      //     "param_it->preg->alpha = " << param_it->preg->alpha << std::endl;

      if (useful.n_elem > 0) {
        if (param_it->preg->v_exist) {
          arma::vec vv = param_it->preg->v_full.elem(finites);
          vv = vv.elem(useful);
          // Rcpp::Rcout << "fun::lambda_max -> vv " << vv.t() << std::endl;

          di = di.elem(arma::find(vv > 0));
          vv = vv.elem(arma::find(vv > 0));
          // Rcpp::Rcout << "fun::lambda_max -> di" << di.t() << std::endl <<
            // "fun::lambda_max -> vv " << vv.t() << std::endl;

          lam_max(param_it - vparam.begin()) = arma::max(arma::abs(di) / vv);
        } else {
          lam_max(param_it - vparam.begin()) = arma::max(arma::abs(di));
        }
      }

      lam_max(param_it - vparam.begin()) /= param_it->preg->alpha;
    }
    diff_it++;
  }

  return lam_max;
}

// get maximal lambda and scale lambda_factors
double fun::get_and_set_lambda_max(arma::vec &lamfac) {
  arma::vec lam_max = lambda_max();
  double lm = arma::max(lam_max);

  // lamfac
  lamfac.zeros(vparam.size());

  // Set the lambda factors and scale lambdas
  std::vector<param>::iterator param_it = vparam.begin();
  for (param_it = vparam.begin(); param_it < vparam.end(); param_it++) {
    if (lam_max(param_it - vparam.begin()) > 0) {
      param_it->lambda_factor *= lam_max(param_it - vparam.begin()) / lm;
      param_it->preg->lambda  *= lam_max(param_it - vparam.begin()) / lm;
    }

    lamfac(param_it - vparam.begin()) = param_it->lambda_factor;
  }

  return lm;
}



//=================================
// am class
//=================================

// W is always scaled so it means to 1.
am::am(arma::vec Y_, arma::vec W_, likeli* plikeli_, std::vector<param> vparam_, std::vector<fit*> vpfit_, bool centerY_) :
  fun(vparam_),
  W(W_ / arma::mean(W_)),
  W_exist(arma::any(W_ < 1.0) || arma::any(W_ > 1.0)),
  Y(centerY_ ? Y_ - arma::dot(Y_, W) / W.n_elem : Y_),
  centerY(centerY_),
  plikeli(plikeli_),
  vpfit(vpfit_) {

  if (vpfit.size() != vparam.size()) {
    Rcpp::stop("error in am -> Length of 'vparam' does not equal length of 'vfit'.");
  }

  // Fill it over, so that they match
  std::vector<fit*>::iterator pfit_it = vpfit.begin();
  std::vector<param>::iterator param_it = vparam.begin();
  for (pfit_it = vpfit.begin(); pfit_it < vpfit.end(); pfit_it++) {
    (*pfit_it)->pparam    = &(*param_it);
    (*pfit_it)->pW        = &W;
    (*pfit_it)->pW_exist  = &W_exist;

    /* Bring them up to date
     * Note: if the old pparam->position matches param_it->position,
     * then no calculations are necessary and none a carried out
     * (remember that 'fit' does a check for that)
     */
    (*pfit_it)->set_active();
    (*pfit_it)->evaluate();
    (*pfit_it)->evaluate_differential();

    if ((*pfit_it)->pparam != &(*param_it)) {
      Rcpp::stop("error in am -> param of function and param of fit do not match.");
    }

    param_it ++;
  }

  // Update value
  evaluate();
}

void am::evaluate(void) {
  // Good idea to include this here, so we can terminate from R and reset value
  Rcpp::checkUserInterrupt();
  value = 0.0;

  // Set iterator and sum of fits
  std::vector<fit*>::iterator pfit_it = vpfit.begin();
  arma::vec sum_of_fits = arma::zeros(Y.n_elem);

  // Iterate through to fill up sum of fits
  for (pfit_it = vpfit.begin(); pfit_it < vpfit.end(); pfit_it++) {
    // Call evaluate, which updates its internal fit_value
    (*pfit_it)->evaluate();

    // Add it to the rest
    sum_of_fits += (*pfit_it)->value;
  }

  // Evaluate the likelihood between fit and obs
  if (W_exist) {
    value = plikeli->eval(Y, sum_of_fits, W);
  } else {
    value = plikeli->eval(Y, sum_of_fits);
  }
}

void am::evaluate_differential(std::vector<param>::iterator it) {
  // Set iterator and sum of fits
  std::vector<fit*>::iterator pfit_it = vpfit.begin();
  arma::vec sum_of_fits = arma::zeros(Y.n_elem);

  // Iterate through to fill up sum of fits
  for (pfit_it = vpfit.begin(); pfit_it < vpfit.end(); pfit_it++) {
    // Call evaluate, which updates its internal fit_value
    (*pfit_it)->evaluate();

    // Add it to the rest
    sum_of_fits += (*pfit_it)->value;
  }

  // Translate to fit iterator
  pfit_it = vpfit.begin() + (it - vparam.begin());
  if ((*pfit_it)->pparam != &(*it)) {
    Rcpp::stop("error in am::evaluate_diff -> fit's pparam and param do not match.");
  }

  // Evaluate the differential
  (*pfit_it)->evaluate_differential();

  if (W_exist) {
    it->differential = (*pfit_it)->value_dposition.t() * plikeli->eval_dfit(Y, sum_of_fits, W);
  } else {
    it->differential = (*pfit_it)->value_dposition.t() * plikeli->eval_dfit(Y, sum_of_fits);
  }
}

void am::set_active(std::vector<param>::iterator it) {
  // Switch to fit iterator
  std::vector<fit*>::iterator pfit_it = vpfit.begin();
  pfit_it += it - vparam.begin();

  if ((*pfit_it)->pparam != &(*it)) {
    Rcpp::stop("error in am::set_active -> fit's pparam and param do not match.");
  }

  (*pfit_it)->set_active();
}


