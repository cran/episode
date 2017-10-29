//=============================================
/* ode_mapping.cpp
 *
 * Content:
 * - mapping and subclasses definitions
 * - numint
 *
 */
//=============================================

#include "ode_mapping.h"


// Entrywise and
inline arma::uvec operator&(const arma::uvec& lhs, arma::uvec& rhs){
  // Does bitwise AND for two uvecs, used to find in range. Returns vector of length of lhs, recycles rhs
  arma::uvec ret(arma::size(lhs), arma::fill::zeros);
  if (!lhs.is_empty()) {
    for (unsigned int i = 0; i < lhs.n_elem; ++i) {
      ret(i) = lhs(i) && rhs(i % rhs.n_elem);
    }
  }
  return ret;
}



//=================================
// mapping class
//=================================

mapping::mapping(unsigned int d_, unsigned int nrow_, unsigned int ncol_) :
  d(d_),
  nrow(nrow_),
  ncol(ncol_) {

}

arma::vec mapping::interpolate(const double t1, const arma::vec x1, const double t2, const arma::vec x2, const double tm) {
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
arma::mat mapping::interpolate(const arma::vec tin, const arma::vec tout, const arma::mat xin) {
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




// Simpson rules
void mapping::simpson(arma::vec &x1, arma::mat &y1, arma::vec x2, arma::mat &y2) {
  // Does simpson approx of int(t1, t2, f_dparam(x(s), param)) (unless param_it->x0 then f(x(s), param))
  y1  = y2; // y2 is the new y1 (the new y2 comes from x2), x1 is old x2
  // Rcpp::Rcout << "\t Calling map" << std::endl << "x1: " << x1.t() << std::endl << "x2: " << x2.t() << std::endl;
  y2  = map(x2);
  y1 += y2 + 4 * map((x1 + x2) / 2);
  y1 /= 6;
  x1  = x2;
} // int(t1, t2, map(x(s))) claimed by in y1 * (t2 - t1)

void mapping::simpson(const double t1, arma::vec &x1, arma::mat &y1,
  const double t2, arma::vec x2, arma::mat &y2, arma::vec tm, arma::mat xm) {
  // follows simpson convention: y2 becomes the new y1, x1 is the old x2
  arma::mat integral_ = y1;    // When passed from last simpson_mid call

  // Rcpp::Rcout << "\t\t Calling simpson" << std::endl;
  if (!tm.is_empty()) {
    simpson(x1, y1, xm.col(0), y2);  // First step: from x1 to first xm
    integral_ += (tm(0) - t1) * y1;
    for (unsigned int i = 1; i < tm.n_elem; ++i) {
      simpson(x1, y1, xm.col(i), y2);
      integral_ += (tm(i) - tm(i - 1)) * y1;
    }
  }
  simpson(x1, y1, x2, y2);  // Last step from last xm (stored in x1) to x2, called no matter if tm is empty
  if (!tm.is_empty()) {
    integral_ += (t2 - tm(tm.n_elem - 1)) * y1;
  } else {
    integral_ += (t2 - t1) * y1;
  }
  y1 = integral_; // For return
} // int(t1, t2, map(x(s))) stored in y1 based on the intermediates xm at tm. If no intermediates, it just calls simpson.

arma::mat mapping::simpson(arma::vec tin, arma::vec tout, arma::mat xin, arma::mat xout) {
  // Assumes tin and tout are non-decreasing and represents the whole series, xin matches tin and xout tout

  // For filling in stacked simpson integrals
  arma::mat X = arma::zeros(nrow * (tout.n_elem - 1), ncol);

  // Initialise
  arma::mat y1 = arma::zeros(nrow, X.n_cols);
  arma::vec x1 = xout.col(0);
  arma::mat y2 = map(x1);

  for (unsigned int j = 0; j < tout.n_elem - 1; ++j) {
    // Find tin and xin intermediates
    arma::uvec rhs  = (tout(j) < tin);
    arma::uvec lhs  = (tin < tout(j + 1));
    arma::uvec mids = arma::find(rhs && lhs);

    // Integrate (result sent to y1)
    // Rcpp::Rcout << "\t Calling simpson" << std::endl;
    simpson(tout(j), x1, y1, tout(j + 1), xout.col(j + 1), y2, tin.elem(mids), xin.cols(mids));

    // Store y1 in X
    X.rows(nrow * j, nrow * (j + 1) - 1) = y1;

    // Reset y1 (if not this, then {int(tout(0), tout(i), f(x(s)))}_i for i = 0,...,end, stacked on top of each other is returned)
    y1.zeros();
  }
  return X;
} // {int(tout(i), tout(i + 1), f(x(s)))}_i for i = 0,...,end-1, stacked on top of each other if param_it->x0, else its f_dparam()





arma::mat mpower(arma::vec x, const arma::mat &A) {
  if (x.n_elem != A.n_cols) {
    Rcpp::stop("error in mpower -> length of x does not match number of columns in power matrix.");
  }
  arma::mat out(A.n_rows, 1);

  if(any(x <= 0)){
    arma::uvec pos = arma::find(x > 0);
    arma::uvec non_pos = arma::find(x <= 0);
    out = arma::exp(A.cols(pos) * arma::log(x(pos)));

    /* For non positive x's if the corresponding rows of A contain positive entries
     * the result changed to 0 */
    arma::mat A_np = A.cols(non_pos);
    for (unsigned int i = 0; i < A.n_rows; ++i) {
      if (any(A_np.row(i) > 0)) {
        out(i) = 0;
      }
    }
  } else {
    out = arma::exp(A * arma::log(x));
  }

  return out.t();
}


//=================================
// power class
//=================================

power::power(arma::mat A_) :
  mapping(A_.n_cols, 1, A_.n_rows),
  A(A_) {
  if (arma::any(arma::vectorise(A) < 0)) {
    Rcpp::stop("error in power -> negative entries in A.");
  }
}

arma::mat power::map(arma::vec x) {
  return mpower(x, A);
}



//=================================
// fracpower class
//=================================

fracpower::fracpower(arma::mat A_, arma::mat B_) :
  mapping(A_.n_cols, 1, A_.n_rows),
  A(A_),
  B(B_) {
  if (arma::any(arma::vectorise(A) < 0)) {
    Rcpp::stop("error in fracpower -> negative entries in A.");
  }
  if (arma::any(arma::vectorise(B) < 0)) {
    Rcpp::stop("error in fracpower -> negative entries in B.");
  }
  if (A.n_cols != B.n_cols) {
    Rcpp::stop("error in fracpower -> A and B does not have same number of columns.");
  }
  if (A.n_rows != B.n_rows && 1 != B.n_rows) {
    Rcpp::stop("error in fracpower -> Number of rows in B not equal to 1 or that of A.");
  }
}

arma::mat fracpower::map(arma::vec x) {
  if (B.n_rows == 1) {
    arma::mat Bp = mpower(x, B);
    return mpower(x, A) / (1.0 + Bp(0));
  } else {
    return mpower(x, A) / (1.0 + mpower(x, B));
  }
}



// Creates a mapping-pointer
mapping* create_pmap(std::string type, arma::mat A, arma::mat B) {
  // Create empty pointer and fill it
  mapping* pmap = NULL;
  if (type == "power") {
    power power_(A);
    pmap = new power(power_);
  } else if (type == "fracpower") {
    fracpower frac_(A, B);
    pmap = new fracpower(frac_);
  } else {
    Rcpp::stop("type not recognised!");
  }
  return pmap;
}



//=================================
// exports
//=================================

//' Numerical integration of powers and fractions of powers via simpson rule
//'
//' @description Evaluates numerical integrals of powers or fractions of powers of a d-dimensional function x.
//'
//' @param time Vector (n) holding time points in between which the integrals are evaluated. Must be one series only (i.e., no decrements).
//' @param x Matrix (mx(d+1)) holding discretised function. First column is time, must have no decrements. The remaining columns are function values, which will be interpolated.
//' @param type String telling type, must be "power" or "fracpower".
//' @param A Matrix (pxd) holding powers in rows.
//' @param B Matrix (pxd) holding powers in rows.
//'
//' @seealso imd, aim
//'
//' @return A matrix of dimension (m-1)xp holding the row-concatinated blocks of integrals:
//' If \code{type} is "power" it evaluates the numerical integrals
//' \deqn{(\int^{t_{i+1}}_{t_i}{x(s)^A})_i}
//' where t_i are entries in first column of \code{x} and \eqn{x(s)} is constructed via a linear interpolation of \code{x}.
//' If \code{type} is "fracpower" it evaluates the numerical integrals
//' \deqn{(\int^{t_{i+1}}_{t_i}{x(s)^A / (1 + x(s)^B)})_i}
//' where the fraction is evaluated coordinate wise.
//'
//' The numerical integration uses the simpson rule using the intermediate time points in \code{time} in between any two consecutive time points in the first column of \code{x}.
//' To get more accurate integrals include more intermediate time points in \code{time}.
//'
//' @examples
//' # Trajectory of power law kinetics system
//' A <- matrix(c(1, 0, 1,
//'               1, 1, 0), byrow = TRUE, nrow = 2)
//' p <- plk(A)
//' x0 <- c(10, 4, 1)
//' theta <- matrix(c(0, -0.25,
//'                   0.75, 0,
//'                   0, -0.1), byrow = TRUE, nrow = 3)
//' Time <- seq(0, 1, by = .025)
//' traj <- numsolve(p, Time, x0, theta)
//'
//' # Example: Integrate traj(s)^A between the time points in 'ti'
//' ti <- seq(0, 1, by = .1)
//' ss <- numint(time = ti, x = traj, type = "power", A = A, B = A)
//'
//' # Example: Integrate traj(s)^A / (1 + traj(s)^B) between the time points in 'ti'
//' B <- matrix(c(0, 2, 1,
//'               2, 1, 0), byrow = TRUE, nrow = 2)
//' ss <- numint(time = ti, x = traj, type = "fracpower", A = A, B = B)
//'
//' @export
// [[Rcpp::export]]
arma::mat numint(arma::vec time, arma::mat x, std::string type, arma::mat A, arma::mat B) {
  if (A.n_cols + 1 != x.n_cols) {
    Rcpp::Rcout << "Number of columns in A: " << A.n_cols << " and number of columns in x: " << x.n_cols << std::endl;
    Rcpp::stop("error in power_integral -> number of columns in A is not one less than that of x.");
  }

  // Create mapping
  mapping* pmap = create_pmap(type, A, B);

  // Data
  arma::vec tin = x.col(0);
  arma::vec tout = time;
  arma::mat xin = arma::trans(x.cols(1, x.n_cols - 1));

  // Get xout, which is interpolated xin to time point in tout
  arma::mat xout = pmap->interpolate(tin, tout, xin);

  // Rcpp::Rcout << "Calling simpson" << std::endl;
  return pmap->simpson(tin, tout, xin, xout);
}

