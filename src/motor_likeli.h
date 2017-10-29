# ifndef MOTOR_LIKELI_H_
# define MOTOR_LIKELI_H_

#include <RcppArmadillo.h>

// likelihood class, only holds three virtual functions, which match fit and obs
class likeli {
protected:
  likeli(void);

public:
  likeli(const likeli& other);
  likeli& operator=(const likeli& rhs);

  // Virtual destructor and clone
  virtual ~likeli() {};
  virtual likeli* shallow_Clone() = 0;

  // fit is supplied by subclass, evals fit for given obs and fitted value (hat)
  double virtual eval(arma::vec obs, arma::vec fit, arma::vec w_i) = 0;
  double virtual eval(arma::vec obs, arma::vec fit) = 0;

  arma::vec virtual eval_dfit(arma::vec obs, arma::vec fit, arma::vec w_i) = 0;
  arma::vec virtual eval_dfit(arma::vec obs, arma::vec fit) = 0;

  arma::mat virtual eval_d2fit(arma::vec obs, arma::vec fit, arma::vec w_i) = 0;
  arma::mat virtual eval_d2fit(arma::vec obs, arma::vec fit) = 0;
};

// Gaussian likelihood subclass
class gaussian : public likeli {
public:
  gaussian(void);
  gaussian* shallow_Clone();

  double eval(arma::vec obs, arma::vec fit, arma::vec w_i);
  double eval(arma::vec obs, arma::vec fit);

  arma::vec eval_dfit(arma::vec obs, arma::vec fit, arma::vec w_i);
  arma::vec eval_dfit(arma::vec obs, arma::vec fit);

  arma::mat eval_d2fit(arma::vec obs, arma::vec fit, arma::vec w_i);
  arma::mat eval_d2fit(arma::vec obs, arma::vec fit);
};

# endif
