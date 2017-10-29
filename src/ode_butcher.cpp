//=============================================
/* ode_butcher.cpp
 *
 * Content:
 * - butcher tableaus for rk embedded pairs
 */
//=============================================


#include "ode_butcher.h"

butcher::butcher(std::string str_) {
  if (str_ == "rk23") {
    // Runge-Kutta order 2/order 3 embedded pair
    x << 0 << 1 << 1./4. << 1./2. << 1./6.  << arma::endr
      << 0 << 0 << 1./4. << 1./2. << 1./6.  << arma::endr
      << 0 << 0 << 0     << 0     << 2./3.  << arma::endr;
    t << 0 << 1 << 1./2. << arma::endr;

  } else if (str_ == "bs23") {
    // Bogacki-Shampine order 2/order 3 embedded pair
    x << 0 << 1./2. << 0      << 2./9. << 7./24.  << 2./9.  << arma::endr
      << 0 << 0     << 3./4.  << 1./3. << 1./4.   << 1./3.  << arma::endr
      << 0 << 0     << 0      << 4./9. << 1./3.   << 4./9.  << arma::endr
      << 0 << 0     << 0      << 0     << 1./8.   << 0      << arma::endr;
    t << 0 << 1./2. << 3./4.  << 1     << arma::endr;

  } else if (str_ == "dp45") {
    // Dormand-Prince order 4/order 5 embedded pair
    x << 0 << 1./5. << 3./40. << 44./45.  << 19372./6561.   << 9017./3168.    << 35./384.     << 5179./57600.     << 35./384.     << arma::endr
      << 0 << 0     << 9./40. << -56./15. << -25360./2187.  << -355./33.      << 0            << 0                << 0            << arma::endr
      << 0 << 0     << 0      << 32./9.   << 64448./6561.   << 46732./5247.   << 500./1113.   << 7571./16695.     << 500./1113.   << arma::endr
      << 0 << 0     << 0      << 0        << -212./729.     << 49./176.       << 125./192.    << 393./640.        << 125./192.    << arma::endr
      << 0 << 0     << 0      << 0        << 0              << -5103./18656.  << -2187./6784. << -92097./339200.  << -2187./6784. << arma::endr
      << 0 << 0     << 0      << 0        << 0              << 0              << 11./84.      << 187./2100.       << 11./84.      << arma::endr
      << 0 << 0     << 0      << 0        << 0              << 0              << 0            << 1./40.           << 0            << arma::endr;
    t << 0 << 1./5. << 3./10. << 4./5.    << 8./9.          << 1              << 1            << arma::endr;

  } else if (str_ == "rkf45") {
    // Runge-Kutta-Fehlberg order 4/order 5 embedded pair
    x << 0 << 1./4. << 3./32. << 1932./2197.  << 439./216.    << -8./27.      << 25./216.     << 16./135.       << arma::endr
      << 0 << 0     << 9./32. << -7200./2197. << -8           << 2            << 0            << 0              << arma::endr
      << 0 << 0     << 0      << 7296./2197.  << 3680./513.   << -3544./2565. << 1408./2565.  << 6656./12825.   << arma::endr
      << 0 << 0     << 0      << 0            << -845./4104.  << 1859./4104.  << 2197./4104.  << 28561./56430.  << arma::endr
      << 0 << 0     << 0      << 0            << 0            << -11./40.     << -1./5.       << -9./50.        << arma::endr
      << 0 << 0     << 0      << 0            << 0            << 0            << 0            << 2./55.         << arma::endr;
    t << 0 << 1./4. << 3./8.  << 12./13.      << 1            << 1./2.        << arma::endr;

  } else {
    Rcpp::stop("Unrecognized Butcher tableau specification for solver!");
  }
}

/* Note there's room for improvement:
 * 1) Since first column is 0, we can avoid that.
 * 2) Some schemes follow the FSAL-principle (First Same As Last), i.e., last fxs is first fxs in next step
 * 3) Some schemes reuses the linear combination of fxs in the last steps for z
 *
 * Under 2) include FSAL-flag in struct which is sent to rkemp and in rkemp include xs in addition to fxs
 */
