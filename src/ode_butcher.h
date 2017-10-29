# ifndef ODE_BUTCHER_H_
# define ODE_BUTCHER_H_

#include <RcppArmadillo.h>

// Purpose: define Butcher tableaus for embedded runge kutta pairs

struct butcher {
  arma::mat x;
  arma::vec t;

  // Constructor
  butcher (std::string str_);
};

# endif
