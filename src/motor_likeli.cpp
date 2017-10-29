//=============================================
/* motor_likeli.cpp
*
* Content:
* - likeli and subclasses definitions
*
*/
//=============================================

#include "motor_likeli.h"


//=================================
// likeli class
//=================================

/* Likelihood functions for matching a fitted value with observation
* Score and second order differentials included
*/

// Constructor
likeli::likeli(void) {}

// Copy constructor
likeli::likeli(const likeli& other) {}

// Assignment operator (necessary before C++11, since constant member variables cannot be assigned using default assignment)
likeli& likeli::operator=(const likeli& rhs) {
  return *this;
}


//=================================
// gaussian subclass
//=================================

gaussian::gaussian(void) {}

gaussian* gaussian::shallow_Clone() {
  return &(*this);  // No need to call "new", since no pointers in likeli or gaussian.
}

double gaussian::eval(arma::vec obs, arma::vec fit, arma::vec w_i) {
  return arma::sum(arma::square(obs - fit) % w_i) / 2;
}
double gaussian::eval(arma::vec obs, arma::vec fit) {
  return arma::sum(arma::square(obs - fit)) / 2;
}

arma::vec gaussian::eval_dfit(arma::vec obs, arma::vec fit, arma::vec w_i) {
  return - (obs - fit) % w_i;
}
arma::vec gaussian::eval_dfit(arma::vec obs, arma::vec fit) {
  return - (obs - fit);
}

arma::mat gaussian::eval_d2fit(arma::vec obs, arma::vec fit, arma::vec w_i) {
  return arma::diagmat(w_i);
}
arma::mat gaussian::eval_d2fit(arma::vec obs, arma::vec fit) {
  return arma::eye(fit.n_elem, fit.n_elem);
}
