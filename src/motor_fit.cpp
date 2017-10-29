//=============================================
/* motor_fit.cpp
*
* Content:
* - fit and subclasses definitions
*
*/
//=============================================

#include "motor_fit.h"


//=================================
// fit class
//=================================

/* fit turns parameter into fitted values
*  First order diff included
*/

// Constructor
fit::fit(param* pparam_, bool center_, arma::vec* pW_, bool* pW_exist_, bool is_linear_) :
  value(arma::zeros(pW_->n_elem)),
  value_dposition(arma::zeros(pW_->n_elem, pparam_->position.n_elem)),
  pparam(pparam_),
  pW(pW_),
  pW_exist(pW_exist_),
  pos_value(pparam_->position),
  pos_diff(pparam_->position),
  center(center_),
  is_linear(is_linear_) {
}

// Copy constructor
fit::fit(const fit& other) :
  value(other.value),
  value_dposition(other.value_dposition),
  pparam(other.pparam),
  pW(other.pW),
  pW_exist(other.pW_exist),
  pos_value(other.pos_value),
  pos_diff(other.pos_diff),
  center(other.center),
  is_linear(other.is_linear) {
}

// Assignment operator (necessary before C++11, since constant member variables cannot be assigned using default assignment)
fit& fit::operator=(const fit& rhs) {
  return *this;
}

bool fit::change(const arma::vec& refpos) {
  if (pparam->position.n_elem != refpos.n_elem) {
    return true;
  } else if (arma::norm(pparam->position - refpos, "inf") * 10.0 > pparam->c.tol) {
    return true;
  } else {
    return false;
  }
}

void fit::evaluate(void) {

  // Rcpp::Rcout << "Evaluate fit::evaluate. pos_value " << pos_value.t() << std::endl <<
    // "position: " << pparam->position.t() << std::endl;

  if (change(pos_value)) {
    pos_value = pparam->position;
    set_value();

    // Rcpp::Rcout << "New value" << value.t() << std::endl;

  }
}

void fit::evaluate_differential(void) {
  if (!is_linear) {
    if (change(pos_diff)) {
      pos_diff = pparam->position;
      set_value_dposition();
    }
  }
}



//=================================
// lfit subclass
//=================================

lfit::lfit(param* pparam_, bool center_, arma::vec* pW_, bool* pW_exist_, arma::mat X_full_) :
  fit(pparam_, center_, pW_, pW_exist_, true),
  X_full(center_ ? X_full_.each_row() - (pW_->t() / arma::sum(pW_->t())) * X_full_ : X_full_),
  X(X_full) {

  if (pparam->position_full.n_elem != X_full.n_cols) {
    Rcpp::stop("error in lfit -> number of columns in X does not match size of param vector.");
  }

  if (pparam->position_full.n_elem != X_full.n_cols) {
    Rcpp::stop("error in lfit -> number of columns in X does not match size of param vector.");
  }

  // Bring up to date
  set_value();
  set_active();   // Note this is called in stead of set_value_dposition
}

void lfit::set_value(void) {
  value = X * pos_value;
}

void lfit::set_value_dposition(void) {
  /* It is never really evaluated, since fit is flagged with is_linear = true,
   * so in evaluate_differential (the only thing that calls this function) it
   * does not do anything, but value_dposition is updated if new actives are used.
   */
  value_dposition = X;
}

// Makes X combatible with param->position
void lfit::set_active(void) {
  X = X_full.cols(pparam->active);
  value_dposition = X;
}



//=================================
// nlfit subclass
//=================================

nlfit::nlfit(param* pparam_, bool center_, arma::vec* pW_, bool* pW_exist_) :
  fit(pparam_, center_, pW_, pW_exist_, false) {

  // Dont bring up to date since nlfit is abstract
  // // Bring up to date
  // set_value();
  // set_value_dposition();
}

void nlfit::set_value(void) {
  nlfit_set_value();
  if (center) {
    value -= *pW_exist ? arma::dot(value, *pW) / pW->n_elem : arma::mean(value);
  }
}

void nlfit::set_value_dposition(void) {
  nlfit_set_value_dposition();
  if (center) {
    if (*pW_exist) {
      value_dposition.each_row() -= pW->t() * value_dposition / pW->n_elem;
    } else {
      value_dposition.each_row() -= arma::mean(value_dposition, 0);
    }
  }
}


//=================================
// loglin subclass
//=================================

// In this constructor we may modify lower/upper limits and box position
loglin::loglin(param* pparam_, bool center_, arma::vec* pW_, bool* pW_exist_, arma::vec ydot_, arma::vec y_) :
  nlfit(pparam_, center_, pW_, pW_exist_),
  ydot(ydot_),
  y(y_) {
    if (ydot.n_elem != y.n_elem) {
      Rcpp::Rcout << "error in loglin -> ydot and y does not have same length." << std::endl;
    }
    if (y.n_elem != pW->n_elem) {
      Rcpp::Rcout << "error in loglin -> y does not have same length as weights (W)." << std::endl;
    }
    if (pparam->position_full.n_elem != 1) {
      Rcpp::Rcout << "error in loglin -> number of parameters is not 1." << std::endl;
    }

    // Modify lower limits
    double lower_ = pparam->c.tol / 100.0 - arma::min(ydot / y);  // Get new lower
    if (arma::any(pparam->preg->lower_full < lower_)) {
      Rcpp::Rcout << "Warning in loglin -> lower limits changed." << std::endl;

      // Fill in
      pparam->preg->lower_full.elem( arma::find(pparam->preg->lower_full < lower_) ).fill(lower_);
      pparam->preg->lower.elem( arma::find(pparam->preg->lower < lower_) ).fill(lower_);

      // Check if lower violates upper
      if (arma::any(pparam->preg->lower_full > pparam->preg->upper_full)) {
        Rcpp::Rcout << "error in loglin -> lower limits now exceeds upper limits." << std::endl;
      }

      // Box it
      pparam->preg->box(pparam->position);
    }

    // Modify screen
    pparam->screen = false;
    if (pparam->active.n_elem != 1) {
      Rcpp::Rcout << "error in loglin -> number of active is not 1." << std::endl;
    }

    // Bring up to date
    set_active();
    nlfit_set_value();
    nlfit_set_value_dposition();
  }

// uses pos_value to set value
void loglin::nlfit_set_value(void) {
  value = arma::log(ydot + y * pos_value);
}

// uses pos_diff to set value_dposition
void loglin::nlfit_set_value_dposition(void) {
  value_dposition = y / (ydot + y * pos_diff);
}

// uses pparam->active to reduce what ever
void loglin::set_active(void) {
  // Does nothing since this parameter must not be screened.
}


//=================================
// hillA subclass
//=================================

// In this constructor we may modify lower/upper limits and box position
hillA::hillA(param* pparam_, bool center_, arma::vec* pW_, bool* pW_exist_, arma::mat x, arma::mat A) :
  nlfit(pparam_, center_, pW_, pW_exist_),
  xA_full(arma::exp(arma::log(x) * A.t())),
  xA(xA_full) {
  if (pparam->position_full.n_elem != A.n_rows) {
    Rcpp::Rcout << "error in hillA -> number of parameters is not equal to number of rows in A." << std::endl;
  }

  // Modify lower limits
  double lower_ = 0.0;
  if (arma::any(pparam->preg->lower_full < lower_)) {
    Rcpp::Rcout << "warning in hillA -> lower limits changed." << std::endl;

    // Fill in
    pparam->preg->lower_full.elem( arma::find(pparam->preg->lower_full < lower_) ).fill(lower_);
    pparam->preg->lower.elem( arma::find(pparam->preg->lower < lower_) ).fill(lower_);

    // Check if lower violates upper
    if (arma::any(pparam->preg->lower_full > pparam->preg->upper_full)) {
      Rcpp::Rcout << "error in hillA -> lower limits now exceeds upper limits." << std::endl;
    }

    // Box it
    pparam->preg->box(pparam->position);
  }

  // Bring up to date
  set_active();
  nlfit_set_value();
  nlfit_set_value_dposition();
}

// uses pos_value to set value
void hillA::nlfit_set_value(void) {
  arma::mat cxA = xA.each_row() % pos_value.t();
  value = arma::sum(arma::log(cxA + 1.0) - arma::log(cxA + pparam->c.tol / 100.0), 1);
}

// uses pos_diff to set value_dposition
void hillA::nlfit_set_value_dposition(void) {
  arma::mat cxA = xA.each_row() % pos_diff.t();
  value_dposition = xA / (cxA + 1.0) - xA / (cxA + pparam->c.tol / 100.0);
}

// uses pparam->active to reduce what ever
void hillA::set_active(void) {
  xA = xA_full.cols(pparam->active);
}


//=================================
// hillB subclass
//=================================

// In this constructor we may modify lower/upper limits and box position
hillB::hillB(param* pparam_, bool center_, arma::vec* pW_, bool* pW_exist_, arma::mat x, arma::mat B) :
  nlfit(pparam_, center_, pW_, pW_exist_),
  xB_full(arma::exp(arma::log(x) * B.t())),
  xB(xB_full) {
  if (pparam->position_full.n_elem != B.n_rows) {
    Rcpp::Rcout << "error in hillB -> number of parameters is not equal to number of rows in B." << std::endl;
  }

  // Modify lower limits
  double lower_ = 0.0;
  if (arma::any(pparam->preg->lower_full < lower_)) {
    Rcpp::Rcout << "warning in hillB -> lower limits changed." << std::endl;

    // Fill in
    pparam->preg->lower_full.elem( arma::find(pparam->preg->lower_full < lower_) ).fill(lower_);
    pparam->preg->lower.elem( arma::find(pparam->preg->lower < lower_) ).fill(lower_);

    // Check if lower violates upper
    if (arma::any(pparam->preg->lower_full > pparam->preg->upper_full)) {
      Rcpp::Rcout << "error in hillB -> lower limits now exceeds upper limits." << std::endl;
    }

    // Box it
    pparam->preg->box(pparam->position);
  }

  // Bring up to date
  set_active();
  nlfit_set_value();
  nlfit_set_value_dposition();
}

// uses pos_value to set value
void hillB::nlfit_set_value(void) {
  value = arma::sum(arma::log(xB.each_row() % pos_value.t() + 1.0), 1);
}

// uses pos_diff to set value_dposition
void hillB::nlfit_set_value_dposition(void) {
  value_dposition = xB / (xB.each_row() % pos_diff.t() + 1.0);
}

// uses pparam->active to reduce what ever
void hillB::set_active(void) {
  xB = xB_full.cols(pparam->active);
}
