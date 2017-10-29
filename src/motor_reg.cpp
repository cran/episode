//=============================================
/* motor_reg.cpp
*
* Content:
* - reg and subclasses definitions
* - create reg pointer
* - stretch_lu
* - evaluate penalty
*/
//=============================================


#include "motor_reg.h"

//=================================
// reg class
//=================================

// Constructor
reg::reg(arma::vec lower_full_, arma::vec upper_full_, arma::vec v_full_, double lambda_) :
  lower_full(lower_full_),
  upper_full(upper_full_),
  v_full(v_full_),
  lower(lower_full_),
  upper(upper_full_),
  v(v_full_),
  v_exist(arma::any(v_full_ < 1.0) || arma::any(v_full_ > 1.0)),
  lambda(lambda_) {

}

// Destructor
reg::~reg() {}

// Finds those with increment at param or already non-zero param
arma::uvec reg::get_active(arma::vec position, arma::vec direction, double tau) {
  arma::uvec position_nonzero = (position != 0);
  arma::uvec position_changes = (prox(position, direction, tau) != position);

  // for(unsigned int i = 0; i < position_nonzero.n_elem; ++i) {
  //   position_nonzero(i) = position_nonzero(i) || position_changes(i);
  // }
  return arma::find(position_nonzero || position_changes);
}

// Clamps any vector of same length as u and l between these limits
void reg::box(arma::vec &x) {
  for (unsigned int i = 0; i < x.n_elem; ++ i) {
    if (x(i) > upper(i)) {
      x(i) = upper(i);
    } else if (x(i) < lower(i)) {
      x(i) = lower(i);
    }
  }
}

// Reduce reg to only a subset of the parameters
void reg::reduce(arma::uvec active_) {
  if (v_exist) {
    v = v_full.elem(active_);
  }
  lower = lower_full.elem(active_);
  upper = upper_full.elem(active_);
}

// Restore reg to original problem
void reg::reduce(void) {
  v = v_full;
  lower = lower_full;
  upper = upper_full;
}





//=================================
// l1 subclass
//=================================

// Constructor
l1::l1(arma::vec lower_full_, arma::vec upper_full_, arma::vec v_full_, double lambda_) :
  reg(lower_full_, upper_full_, v_full_, lambda_) {
  alpha = 1.0;
}

// Destructor is just implicit (no pointer fields here)

// Clone functions
l1* l1::deep_Clone() {
  l1* pl1 = new l1(*this);
  return pl1;
}
l1* l1::shallow_Clone() {
  return &(*this);
}

double l1::penalty(arma::vec position) {
  return v_exist ? lambda * arma::dot(arma::abs(position), v) : lambda * arma::norm(position, 1);
}

arma::vec l1::prox(arma::vec position, arma::vec direction, double tau) {
  // OR, whatever is faster: x = (x - tau * (p + lambda * arma::sign(x - tau * p) % v)) % arma::conv_to<arma::vec>::from(arma::abs(x - tau * p) > lambda * tau * v);
  if (v_exist) {
    position = arma::sign(position + tau * direction) % arma::clamp(arma::abs(position + tau * direction) - lambda * tau * v, 0, arma::datum::inf);
  } else {
    position = arma::sign(position + tau * direction) % arma::clamp(arma::abs(position + tau * direction) - lambda * tau, 0, arma::datum::inf);
  }
  box(position);
  return position;
}

arma::vec l1::d2pen(arma::vec position) {
  return arma::zeros(arma::size(position));
}



//=================================
// l2 subclass
//=================================

// Constructor
l2::l2(arma::vec lower_full_, arma::vec upper_full_, arma::vec v_full_, double lambda_) :
  reg(lower_full_, upper_full_, v_full_, lambda_) {
  alpha = 0.001;
}

// Clone functions
l2* l2::deep_Clone() {
  l2* pl2 = new l2(*this);
  return pl2;
}
l2* l2::shallow_Clone() {
  return &(*this);
}

double l2::penalty(arma::vec position) {
  return v_exist ? lambda * arma::dot(arma::square(position), v) / 2 : lambda * arma::sum(arma::square(position)) / 2;
}

arma::vec l2::prox(arma::vec position, arma::vec direction, double tau) {
  if (v_exist) {
    position = position - tau * (lambda * position % v - direction);
  } else {
    position = position - tau * (lambda * position - direction);
  }
  box(position);
  return position;
}

arma::vec l2::d2pen(arma::vec position) {
  if (v_exist) {
    return v;
  } else {
    return arma::ones(arma::size(position));
  }
}



//=================================
// elnet subclass
//=================================

// Constructor
elnet::elnet(arma::vec lower_full_, arma::vec upper_full_, arma::vec v_full_, double lambda_, double a_) :
  reg(lower_full_, upper_full_, v_full_, lambda_),
  a(a_) {
  if (a > 1 || a < 0) {
    Rcpp::stop("error in elnet::elnet -> a not in [0, 1]");
  }
  alpha = std::max(0.001, a);
}

// Clone functions
elnet* elnet::deep_Clone() {
  elnet* pelnet = new elnet(*this);
  return pelnet;
}
elnet* elnet::shallow_Clone() {
  return &(*this);
}

double elnet::penalty(arma::vec position) {
  return v_exist ? lambda * (arma::dot(a * arma::abs(position) + (1 - a) * arma::square(position) / 2, v)) : lambda * (a * arma::norm(position, 1) + (1 - a) * arma::sum(arma::square(position)) / 2);
}

arma::vec elnet::prox(arma::vec position, arma::vec direction, double tau) {
  if (v_exist) {
    direction -= lambda * (1 - a) * (position % v);
    position = arma::sign(position + tau * direction) % arma::clamp(arma::abs(position + tau * direction) - lambda * a * tau * v, 0, arma::datum::inf);
  } else {
    direction -= lambda * (1 - a) * position;
    position = arma::sign(position + tau * direction) % arma::clamp(arma::abs(position + tau * direction) - lambda * a * tau, 0, arma::datum::inf);
  }
  box(position);
  return position;
}

arma::vec elnet::d2pen(arma::vec position) {
  if (v_exist) {
    return (1 - a) * v;
  } else {
    return (1 - a) * arma::ones(arma::size(position));
  }
}



//=================================
// mcp subclass
//=================================

// Helping functions
arma::vec mcp::p(arma::vec position) {
  position = arma::abs(position);
  for (arma::vec::iterator i = position.begin(); i != position.end(); ++i) {
    (*i) = ((*i) <= a * lambda) ? (*i) - std::pow(*i, 2) / (2 * a * lambda) : a * lambda / 2;
  }
  return position;
}
arma::vec mcp::pprime(arma::vec position) {
  return arma::clamp(1 - arma::abs(position) / (a * lambda), 0, arma::datum::inf);
}

// Constructor
mcp::mcp(arma::vec lower_full_, arma::vec upper_full_, arma::vec v_full_, double lambda_, double a_) :
  reg(lower_full_, upper_full_, v_full_, lambda_),
  a(a_) {
  if (a <= 1) {
    Rcpp::stop("error in mcp::mcp -> a greater than 1.");
  }
  alpha = 1.0;
}

// Clone functions
mcp* mcp::deep_Clone() {
  mcp* pmcp = new mcp(*this);
  return pmcp;
}
mcp* mcp::shallow_Clone() {
  return &(*this);
}

double mcp::penalty(arma::vec position) {
  return v_exist ? lambda * arma::dot(p(position), v) : lambda * arma::sum(p(position));
}

arma::vec mcp::prox(arma::vec position, arma::vec direction, double tau) {
  if (v_exist) {
    position = arma::sign(position + tau * direction) % arma::clamp(arma::abs(position + tau * direction) - lambda * tau * v % pprime(position), 0, arma::datum::inf);
  } else {
    position = arma::sign(position + tau * direction) % arma::clamp(arma::abs(position + tau * direction) - lambda * tau * pprime(position), 0, arma::datum::inf);
  }
  box(position);
  return position;
}

arma::vec mcp::d2pen(arma::vec position) {
  arma::vec ret = arma::conv_to<arma::vec>::from((position <= a * lambda));
  ret /= - lambda * a;
  if (v_exist) {
    return ret % v;
  } else {
    return ret;
  }
}


//=================================
// scad subclass
//=================================

// Helping functions
arma::vec scad::p(arma::vec position) {
  position = arma::abs(position);
  for (arma::vec::iterator i = position.begin(); i != position.end(); ++i) {
    if ((*i) >= lambda * a) {
      (*i) = lambda + (a - 1) * lambda / 2;
    } else if ((*i) >= lambda) {
      (*i) = lambda + ((*i) - lambda) * (1 - (*i - lambda) / (2 * lambda * (a - 1)));
    } else {

    }
  }
  return position;
}
arma::vec scad::pprime(arma::vec position) {
  return arma::clamp((a * lambda - arma::abs(position)) / (lambda * (a - 1)), 0, 1);
}

// Constructor
scad::scad(arma::vec lower_full_, arma::vec upper_full_, arma::vec v_full_, double lambda_, double a_) :
  reg(lower_full_, upper_full_, v_full_, lambda_),
  a(a_) {
  if (a <= 2) {
    Rcpp::stop("error in scad::scad -> a greater than 2.");
  }
  alpha = 1.0;
}

// Clone functions
scad* scad::deep_Clone() {
  scad* pscad = new scad(*this);
  return pscad;
}
scad* scad::shallow_Clone() {
  return &(*this);
}

double scad::penalty(arma::vec position) {
  return v_exist ? lambda * arma::dot(p(position), v) : lambda * arma::sum(p(position));
}

arma::vec scad::prox(arma::vec position, arma::vec direction, double tau) {
  if (v_exist) {
    // Rcpp::Rcout << "Call scad::prox. Position: " << position.t() << std::endl << " v " << v.t() << std::endl;
    position = arma::sign(position + tau * direction) % arma::clamp(arma::abs(position + tau * direction) - lambda * tau * v % pprime(position), 0, arma::datum::inf);
  } else {
    position = arma::sign(position + tau * direction) % arma::clamp(arma::abs(position + tau * direction) - lambda * tau * pprime(position), 0, arma::datum::inf);
  }
  box(position);
  return position;
}

arma::vec scad::d2pen(arma::vec position) {
  arma::vec ret = arma::conv_to<arma::vec>::from((position <= a * lambda) % (position >= lambda));
  ret /= - lambda * (a - 1);
  if (v_exist) {
    return ret % v;
  } else {
    return ret;
  }
}



//=================================
// none subclass
//=================================

// Constructor
none::none(arma::vec lower_full_, arma::vec upper_full_, arma::vec v_full_, double lambda_) :
  reg(lower_full_, upper_full_, v_full_, lambda_) {
  alpha = 1.0;
}

// Clone functions
none* none::deep_Clone() {
  none* pnone = new none(*this);
  return pnone;
}
none* none::shallow_Clone() {
  return &(*this);
}

double none::penalty(arma::vec position) {
  // arma::uvec positionnon0 = arma::find(position);  // old l0 reg
  // return v_exist ? lambda * arma::sum(v_sub.elem(positionnon0)) : lambda * positionnon0.n_elem;
  return 0.0;
}

arma::vec none::prox(arma::vec position, arma::vec direction, double tau) {
  position += tau * direction;
  box(position);
  return position;
}

arma::vec none::d2pen(arma::vec position) {
  return arma::zeros(arma::size(position));
}





//=================================
// Create pointer to reg object
//=================================

reg* create_preg(std::string reg_type, arma::vec lower, arma::vec upper, arma::vec v, double lambda, double a) {
  // Create regularisation pointer, handle 'a' (NULL) when making others
  reg* preg = NULL;
  if (reg_type == "l1") {
    l1 l1_ (lower, upper, v, lambda);
    reg* pl1_ = &l1_;
    preg = pl1_->deep_Clone();
  } else if (reg_type == "l2") {
    l2 l2_ (lower, upper, v, lambda);
    reg* pl2_ = &l2_;
    preg = pl2_->deep_Clone();
  } else if (reg_type == "elnet") {
    elnet elnet_ (lower, upper, v, lambda, a);
    reg* pelnet_ = &elnet_;
    preg = pelnet_->deep_Clone();
  } else if (reg_type == "scad") {
    scad scad_ (lower, upper, v, lambda, a);
    reg* pscad_ = &scad_;
    preg = pscad_->deep_Clone();
  } else if (reg_type == "mcp") {
    mcp mcp_ (lower, upper, v, lambda, a);
    reg* pmcp_ = &mcp_;
    preg = pmcp_->deep_Clone();
  } else if (reg_type == "none") {
    none none_ (lower, upper, v, lambda);
    reg* pnone_ = &none_;
    preg = pnone_->deep_Clone();
  } else {
    Rcpp::stop("reg_type not recognised!");
  }

  return preg;
}

reg* create_preg(Rcpp::List reg_struct, unsigned int p_full) {
  // Penalty factor
  Rcpp::Nullable<arma::vec> v_ = Rcpp::as< Rcpp::Nullable<arma::vec> >(reg_struct["penalty_factor"]);
  arma::vec v = v_.isNotNull() ? Rcpp::as<arma::vec>(v_) : arma::ones(p_full);

  // Lower limits
  Rcpp::Nullable<arma::vec> lower_ = Rcpp::as< Rcpp::Nullable<arma::vec> >(reg_struct["lower"]);
  arma::vec lower(p_full);
  if (lower_.isNotNull()) {
    lower = stretch_lu(Rcpp::as<arma::vec>(lower_), p_full);
  } else {
    lower.fill(-arma::datum::inf);
  }

  // Upper limits
  Rcpp::Nullable<arma::vec> upper_ = Rcpp::as< Rcpp::Nullable<arma::vec> >(reg_struct["upper"]);
  arma::vec upper(p_full);
  if (upper_.isNotNull()) {
    upper = stretch_lu(Rcpp::as<arma::vec>(upper_), p_full);
  } else {
    upper.fill(arma::datum::inf);
  }

  // Reg type
  std::string reg_type = Rcpp::as<std::string>(reg_struct["reg_type"]);

  // Check reg has right dimensions
  if (lower.n_elem != p_full) Rcpp::stop("Length of lower bound on parameter does not match length of parameter.");
  if (upper.n_elem != p_full) Rcpp::stop("Length of upper bound on parameter does not match length of parameter.");
  if (v.n_elem != p_full) Rcpp::stop("Length of penalty_factor does not match length of parameter.");

  reg* preg = create_preg(reg_type, lower, upper, v, Rcpp::as<double>(reg_struct["lambda_factor"]), Rcpp::as<double>(reg_struct["a"]));

  return preg;
}

// Takes vector and unsigned integer, if vector has length = unsigned int, repeat it out
arma::vec stretch_lu(arma::vec v, unsigned int p) {
  if (v.n_elem == 1) {
    double vv = v(0);
    v.resize(p);  v.fill(vv);
  } else {
    if (v.n_elem != p) Rcpp::stop("Length of 'lower_param' or 'upper_param' is neither 1 nor equal to number of parameters.");
  }

  return v;
}



//=================================
// evaluate penalty
//=================================

// //' Evaluate penalty function
// //'
// //' @description This function evaluates the penalty function specified via the \code{ro}-object.
// //'
// //' @param x An integer vector
// //' @export
// // [[Rcpp::export]]
// arma::vec penalty(arma::sp_mat param, Rcpp::List o, Rcpp::List r) {
//   //arma::vec penalty(arma::sp_mat param, arma::vec lambda, std::string reg_type, Rcpp::Nullable<double> a, Rcpp::Nullable<arma::vec> v) {
//
//   // Create reg pointer
//   arma::vec l = stretch_lu(Rcpp::as<arma::vec>(o["lower_param"]), Rcpp::as<unsigned int>(o["p"]));
//   arma::vec* pl = &l;
//   arma::vec u = stretch_lu(Rcpp::as<arma::vec>(o["upper_param"]), Rcpp::as<unsigned int>(o["p"]));
//   arma::vec* pu = &u;
//   reg* preg = create_preg(o, r, pl, pu);
//
//   // Check p consistency
//   if (Rcpp::as<unsigned int>(o["p"]) != param.n_rows) {
//     Rcpp::stop("The number of parameters according to 'ode'-object does not match number of rows in 'param'.");
//   }
//
//   // Prepare param and lambda
//   arma::mat param_ = arma::conv_to<arma::mat>::from(param);
//   arma::vec lambda = Rcpp::as<arma::vec>(ro_struct["lambda"]);
//
//   // Run through lambdas
//   arma::vec ret(param_.n_cols);
//   for (unsigned int i = 0; i < param_.n_cols; i++) {
//     preg->lambda = lambda(i);
//     ret(i) = preg->penalty(param_.col(i));
//   }
//
//   return ret;
// }
