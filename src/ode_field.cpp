//=============================================
/* ode_field.cpp
 *
 * Content:
 * - sc class definition
 * - matrix powers and their differential
 * - field and subclasses definitions
 * - create_pfield
 * - ode_field
 *
 */
//=============================================

#include "ode_field.h"




//=================================
// sc class
//=================================

// Constructor
sc::sc (Rcpp::Nullable<arma::mat> value_) :
  exist(value_.isNotNull()),
  value(exist ? Rcpp::as<arma::mat>(value_) : arma::mat()) {

}

// Member function: return column context (series)
arma::vec sc::col(unsigned int context) {
  if (exist) {
    if (context >= value.n_cols) {
      Rcpp::stop("error in sc::col -> context oob.");
    }
  } else {
    Rcpp::stop("error in sc::col -> col called on scales which do not exist");
  }

  return value.col(context);
}

// Transfers rows of value to another sc
void sc::transfer(sc &other, arma::uvec rows_coor) {
  if (exist) {
    other.exist = true;
    other.value = value.rows(rows_coor);
  } else {
    other.exist = false;
    other.value = arma::mat();
  }
}



//=================================
// matrix powers
//=================================

arma::vec matrix_power(arma::vec x, const arma::mat &P) {
  arma::vec out(P.n_rows);

  if(any(x <= 0)){
    arma::uvec pos = arma::find(x > 0);
    arma::uvec non_pos = arma::find(x <= 0);
    out = exp(P.cols(pos) * log(x(pos)));

    /* For non positive x's if the corresponding rows of P contain positive entries
     * the result changed to 0 */
    arma::mat P_np = P.cols(non_pos);
    for (unsigned int i = 0; i < P.n_rows; ++i) {
      if (any(P_np.row(i) > 0)) {
        out(i) = 0;
      }
    }
  } else {
    out = exp(P * log(x));
  }

  return out;
}

arma::mat dmatrix_power(arma::vec x, const arma::mat &P) {
  arma::mat out = P;
  out.each_col() %= matrix_power(x, P);

  if(any(x <= 0)){
    arma::uvec pos = arma::find(x > 0);
    arma::uvec non_pos = arma::find(x <= 0);
    arma::vec x_ = x;
    x_.elem(non_pos).ones();
    out.each_row() %= arma::trans(1 / x_);

    /* For non positive x's if the corresponding entry of P is larger than 1 or equal to 0, then
     * dx^p/dx -> 0 as x -> 0
     * if p == 1, then it converges to 1
     * if 0 < p < 1, then it converges to inf
     * if p negative it converges to -inf
     */
    arma::mat P_np = P.cols(non_pos);
    for (unsigned int j = 0; j < non_pos.n_elem; ++j) {
      x_ = x;
      x_(non_pos(j)) = 1; // Corresponds to that coordinate not entering in the matrix power
      arma::vec xP = matrix_power(x_, P);

      for (unsigned int i = 0; i < P.n_rows; ++i) {
        double a = P_np(i, j);
        if (a > 1 || a == 0) {
          out(i, non_pos(j)) = 0;
        } else if (a == 1) {
          out(i, non_pos(j)) = xP(i); // Should be multiplied with p, but it is 1.
        } else if (a > 0) {
          out(i, non_pos(j)) = arma::datum::inf;
        } else {
          out(i, non_pos(j)) = -arma::datum::inf;
        }
      }
    }
  } else {
    out.each_row() %= arma::trans(1 / x);
  }

  return out;
}



//=================================
// field class
//=================================

// Constructor (no copy constructor or destructor defined, as ode contains no raw pointers)
field::field (std::vector<param>* pvparam_, unsigned int d_, arma::uvec linear_) :
  pvparam(pvparam_), d(d_), linear(linear_) {}

void field::fill_in_param(Rcpp::Nullable<Rcpp::List> param_) {
  /* If it is null it does nothing, i.e., they remain unchanged,
   * which in most cases means they remain 0 (as this is often called
   * after create_pfield.
   */
  if (param_.isNotNull()) {
    Rcpp::List params_ = Rcpp::as<Rcpp::List>(param_);
    if ((unsigned int) params_.size() + 1 != pvparam->size()) {
      Rcpp::Rcout << "params_.size() = " << params_.size() << " and pvparam->size() - 1 = " << pvparam->size() - 1 << std::endl;
      Rcpp::stop("error in field::fill_in_param -> list of param does not match in length with pvparam.");
    }

    std::vector<param>::iterator param_it = pvparam->begin();
    for (unsigned int i = 0; i < params_.size(); i++) {
      param_it++; // Move forward (first is x0)
      arma::vec fill_in = Rcpp::as<arma::vec>(params_[i]);

      if (param_it->position_full.n_elem != fill_in.n_elem) {
        Rcpp::Rcout << "placeholder size = " << param_it->position_full.n_elem  <<
          " and provided parameter size = " << fill_in.n_elem << std::endl;
        Rcpp::stop("error -> field::fill_in_param: length of param does not match p in 'reg' object in 'ode'-object.");
      } else {
        param_it->position_full = fill_in;
        param_it->position = fill_in.elem(param_it->active);
      }
    }
  }
}

// Calls set active on each parameter
void field::set_active(void) {
  for (std::vector<param>::iterator param_it = pvparam->begin(); param_it < pvparam->end(); param_it++) {
    set_active(param_it);
  }
}



//=================================
// mak subclass
//=================================

// Constructor (note all instances of mak needs to point to the same matrix A and C, so copy constructor and destructor are default)
mak::mak(std::vector<param>* pvparam_, arma::mat* pA, arma::mat* pC) :
  field(pvparam_, pA->n_cols, arma::ones<arma::uvec>(1)),
  A(pA),
  C(pC) {
  A_sub = *A;
  C_sub = *C;
}

// Clone function
mak* mak::shallow_Clone() {
  return &(*this);
  /* Note: the only pointers that mak contains are pA and pC, and
  * since the effective versions used are A_sub and C_sub, we can keep pA and pC the same across multiple instances.
  * 1) We dont use static variables for two reasons: we may want to define two different sets of instances each with their own A and C
  * 2) Static member variables cannot be instanced in abstract classes.
  *
  * Therefore, we do the cloning as above:
  * "this" is a pointer to this object,
  * "*this" is a (shallow) copy of the value in "this" (hence only the pointer A and C are copied, not their value)
  * finally &(*this) returns the address of said copy. All in all A and C points to the same matrices.
  * DO NOT call delete on pointers to ode through instances of mak.
  *
  * Alternative is "new mak(*this)", but that would call deep copy on all the pointers in object, that is not desired.
  */
}
mak* mak::deep_Clone() {
  mak* pmak = new mak(*this);
  return pmak;
}
/* Derived desctructor (since base class has virtual destructor,
  * any delete call to instances of derived class (mak) defined as pointers to base class (ode) calls this one)
  * We use the default constructor (the one below), which implicitly calls destructors on its member variables after the call to "~mak()"
  */
mak::~mak() {

}

// Drift
arma::vec mak::f(double t, arma::vec x) {
  return C_sub * (matrix_power(x, A_sub) % (pvparam->begin() + 1)->position);
}

// State derivate
arma::mat mak::f_dx(double t, arma::vec x) {
  return (C_sub.each_row() % trans((pvparam->begin() + 1)->position)) * dmatrix_power(x, A_sub);
}

// Rate derivate
arma::mat mak::f_dparam(double t, arma::vec x, std::vector<param>::iterator it) {
  if (&(*(pvparam->begin() + 1)) == &(*it)) {
    return C_sub.each_row() % trans(matrix_power(x, A_sub));
  } else {
    Rcpp::stop("error in mak::f_dparam -> pvparam and iterator do not match.");
  }
}

// Reduce and restore system
void mak::set_active(std::vector<param>::iterator it) {
  if (&(*it) == &(*(pvparam->begin()))) {

  } else if (&(*it) == &(*(pvparam->begin() + 1))) {
    A_sub = A->rows((pvparam->begin() + 1)->active);
    C_sub = C->cols((pvparam->begin() + 1)->active);
  } else {
    Rcpp::stop("error in mak::set_active -> iterator out of bounds.");
  }
}



//=================================
// plk subclass
//=================================

// Constructor
plk::plk(std::vector<param>* pvparam_, arma::mat *pA) :
  field(pvparam_, pA->n_cols, arma::ones<arma::uvec>(1)),
  A(pA) {
  A_sub = *A;
}

// Clone function
plk* plk::shallow_Clone() {
  return &(*this);
}
plk* plk::deep_Clone() {
  plk* pplk = new plk(*this);
  return pplk;
}

// Destructor
plk::~plk() {

}

// Drift
arma::vec plk::f(double t, arma::vec x) {
  return (pvparam->begin() + 1)->get_position() * matrix_power(x, A_sub);
}

// State derivate
arma::mat plk::f_dx(double t, arma::vec x) {
  return (pvparam->begin() + 1)->get_position() * dmatrix_power(x, A_sub);
}

// Rate derivate
arma::mat plk::f_dparam(double t, arma::vec x, std::vector<param>::iterator it) {
  // Translate to fit iterator
  if (&(*(pvparam->begin() + 1)) == &(*it)) {
    return it->get_dparam(matrix_power(x, A_sub));
  } else {
    Rcpp::stop("error in plk::f_dparam -> pvparam and iterator do not match.");
  }
}

// Reduce and restore system
void plk::set_active(std::vector<param>::iterator it) {
  if (&(*it) == &(*(pvparam->begin()))) {

  } else if (&(*it) == &(*(pvparam->begin() + 1))) {
    arma::uvec coor = (pvparam->begin() + 1)->active / (pvparam->begin() + 1)->nrow;
    A_sub = A->rows(arma::unique(coor));
    if (A_sub.n_rows != (pvparam->begin() + 1)->r_sub) {
      Rcpp::stop("error in plk::set_active -> number of rows in reduced A does not match r_sub from param.");
    }
  } else {
    Rcpp::stop("error in plk::set_active -> iterator out of bounds.");
  }
}



//=================================
// ratmak subclass
//=================================

// Constructor
ratmak::ratmak(std::vector<param>* pvparam_, arma::mat *pA, arma::mat *pC) :
  field(pvparam_, pA->n_cols, arma::linspace<arma::uvec>(1, 0, 2)),
  A(pA),
  C(pC) {
  A1 = *A;
  A2 = *A;
  C_sub = *C;
}

// Clone function
ratmak* ratmak::shallow_Clone() {
  return &(*this);
}
ratmak* ratmak::deep_Clone() {
  ratmak* pratmak = new ratmak(*this);
  return pratmak;
}

// Destructor
ratmak::~ratmak() {

}

// Drift
arma::vec ratmak::f(double t, arma::vec x) {
  return C_sub * (((pvparam->begin() + 1)->get_position() * matrix_power(x, A1)) / (1 + (pvparam->begin() + 2)->get_position() * matrix_power(x, A2)));
}

// State derivate
arma::mat ratmak::f_dx(double t, arma::vec x) {
  arma::mat ret = ((pvparam->begin() + 1)->get_position().each_col() / (1 + (pvparam->begin() + 2)->get_position() * matrix_power(x, A2))) * dmatrix_power(x, A1);

  // Intermediate necessary, because dgemm in BLAS/LAPACK
  arma::mat tmp = ((pvparam->begin() + 2)->get_position().each_col() % ((pvparam->begin() + 1)->get_position() * matrix_power(x, A1) / arma::square(1 + (pvparam->begin() + 2)->get_position() * matrix_power(x, A2)))) * dmatrix_power(x, A2);
  ret -= tmp;

  return C_sub * ret;
}

// Rate derivate
arma::mat ratmak::f_dparam(double t, arma::vec x, std::vector<param>::iterator it) {
  // Translate to fit iterator
  if (&(*(pvparam->begin() + 1)) == &(*it)) {
    // Wrt K1
    arma::mat ret = C_sub.each_row() / (1 + (pvparam->begin() + 2)->get_position() * matrix_power(x, A2)).t();
    return ret * (it->get_dparam(matrix_power(x, A1)));
  } else if (&(*(pvparam->begin() + 2)) == &(*it)) {
    // Wrt K2
    arma::mat ret = C_sub.each_row() % (((pvparam->begin() + 1)->get_position() * matrix_power(x, A1)) / arma::square(1 + (pvparam->begin() + 2)->get_position() * matrix_power(x, A2))).t();
    return -ret * (it->get_dparam(matrix_power(x, A2)));
  } else {
    Rcpp::stop("error in ratmak::f_dparam -> pvparam and iterator do not match.");
  }
}

// Reduce and restore system
void ratmak::set_active(std::vector<param>::iterator it) {
  if (&(*it) == &(*(pvparam->begin()))) {

  } else if (&(*it) == &(*(pvparam->begin() + 1))) {
    arma::uvec coor = (pvparam->begin() + 1)->active / (pvparam->begin() + 1)->nrow;
    A1 = A->rows(arma::unique(coor));
    if (A1.n_rows != (pvparam->begin() + 1)->r_sub) {
      Rcpp::stop("error in ratmak::set_active -> number of rows in reduced A1 does not match r_sub from param.");
    }
  } else if (&(*it) == &(*(pvparam->begin() + 2))) {
    arma::uvec coor = (pvparam->begin() + 2)->active / (pvparam->begin() + 2)->nrow;
    A2 = A->rows(arma::unique(coor));
    if (A2.n_rows != (pvparam->begin() + 2)->r_sub) {
      Rcpp::stop("error in ratmak::set_active -> number of rows in reduced A2 does not match r_sub from param.");
    }
  } else {
    Rcpp::stop("error in ratmak::set_active -> iterator out of bounds.");
  }
}



//=================================
// rlk subclass
//=================================

// Constructor
rlk::rlk(std::vector<param>* pvparam_, arma::mat* pA1, arma::mat* pA2) :
  field(pvparam_, pA1->n_cols, arma::ones<arma::uvec>(1)),
  A1(pA1),
  A2(pA2) {
  A1_sub = *A1;
  A2_sub = *A2;
}

// Clone function
rlk* rlk::shallow_Clone() {
  return &(*this);
}
rlk* rlk::deep_Clone() {
  rlk* prlk = new rlk(*this);
  return prlk;
}

// Destructor
rlk::~rlk() {

}

// Drift
arma::vec rlk::f(double t, arma::vec x) {
  return (pvparam->begin() + 1)->get_position() * (matrix_power(x, A1_sub) / (1.0 + matrix_power(x, A2_sub)));
}

// State derivate
arma::mat rlk::f_dx(double t, arma::vec x) {
  arma::mat dd1 = dmatrix_power(x, A1_sub);
  arma::mat dd2 = dmatrix_power(x, A2_sub);
  dd1 = dd1.each_col() % (1.0 + matrix_power(x, A2_sub)) - dd2.each_col() % matrix_power(x, A1_sub);
  dd1.each_col() /= arma::square(1.0 + matrix_power(x, A2_sub));

  return (pvparam->begin() + 1)->get_position() * dd1;
}

// Rate derivate
arma::mat rlk::f_dparam(double t, arma::vec x, std::vector<param>::iterator it) {
  // Translate to fit iterator
  if (&(*(pvparam->begin() + 1)) == &(*it)) {
    return it->get_dparam(matrix_power(x, A1_sub) / (1.0 + matrix_power(x, A2_sub)));
  } else {
    Rcpp::stop("error in rlk::f_dparam -> pvparam and iterator do not match.");
  }
}

// Reduce and restore system
void rlk::set_active(std::vector<param>::iterator it) {
  if (&(*it) == &(*(pvparam->begin()))) {

  } else if (&(*it) == &(*(pvparam->begin() + 1))) {
    arma::uvec coor = (pvparam->begin() + 1)->active / (pvparam->begin() + 1)->nrow;
    A1_sub = A1->rows(arma::unique(coor));
    A2_sub = A2->rows(arma::unique(coor));
    if (A1_sub.n_rows != (pvparam->begin() + 1)->r_sub) {
      Rcpp::stop("error in rlk::set_active -> number of rows in reduced A does not match r_sub from param.");
    }
  } else {
    Rcpp::stop("error in rlk::set_active -> iterator out of bounds.");
  }
}








//=================================
// Create ode pointer from R class
//=================================

// Creates a field-pointer (and references to objects the ode-pointer points to)
field* create_pfield(Rcpp::List ode_struct, arma::mat &A, arma::mat &B, std::vector<param> &vparam, arma::mat x0) {
  // Get list of regs
  Rcpp::List rs = Rcpp::as<Rcpp::List>(ode_struct["rs"]);

  // Clear it at fill it up
  vparam.clear();

  // x0 param is pushed back
  Rcpp::List reg_struct = Rcpp::as<Rcpp::List>(rs[0]);

  unsigned int d = Rcpp::as<unsigned int>(ode_struct["d"]);
  if (d != x0.n_rows) {
    Rcpp::Rcout << "d in ode-object = " << d << ", and number of rows in x0 = " << x0.n_rows << std::endl;
    Rcpp::stop("error -> create_pfield: d in ode-object does not match number of rows in x0.");
  }

  vparam.push_back(param(reg_struct, arma::vectorise(x0), d));

  // Create empty pointer and fill it
  field* pfield = NULL;
  if (is_in("mak", ode_struct.attr("class"))) {
    if (ode_struct["sparse"]) {
      A = arma::conv_to<arma::mat>::from(Rcpp::as<arma::sp_mat>(ode_struct["A"]));
      B = arma::conv_to<arma::mat>::from(Rcpp::as<arma::sp_mat>(ode_struct["B"]));
    } else {
      A = Rcpp::as<arma::mat>(ode_struct["A"]);
      B = Rcpp::as<arma::mat>(ode_struct["B"]);
    }
    arma::mat* pA = &A;
    B = arma::trans(B - A);
    arma::mat* pC = &B;

    // Create parameter objects
    reg_struct = Rcpp::as<Rcpp::List>(rs[1]);
    vparam.push_back(param(reg_struct, arma::zeros(A.n_rows), 0));

    mak mak_ (&vparam, pA, pC);
    field* pmak_ = &mak_;
    pfield = pmak_->deep_Clone();
  } else if (is_in("plk", ode_struct.attr("class"))) {
    if (Rcpp::as<bool>(ode_struct["sparse"])) {
      A = arma::conv_to<arma::mat>::from(Rcpp::as<arma::sp_mat>(ode_struct["A"]));
    } else {
      A = Rcpp::as<arma::mat>(ode_struct["A"]);
    }
    arma::mat* pA = &A;

    // Create parameter objects
    reg_struct = Rcpp::as<Rcpp::List>(rs[1]);
    vparam.push_back(param(reg_struct, arma::zeros(A.n_elem), d));

    plk plk_ (&vparam, pA);
    field* pplk_ = &plk_;
    pfield = pplk_->deep_Clone();
  } else if (is_in("ratmak", ode_struct.attr("class"))) {
    if (ode_struct["sparseA"]) {
      A = arma::conv_to<arma::mat>::from(Rcpp::as<arma::sp_mat>(ode_struct["A"]));
    } else {
      A = Rcpp::as<arma::mat>(ode_struct["A"]);
    }
    if (ode_struct["sparseC"]) {
      B = arma::conv_to<arma::mat>::from(Rcpp::as<arma::sp_mat>(ode_struct["C"]));
    } else {
      B = Rcpp::as<arma::mat>(ode_struct["C"]);
    }
    arma::mat* pA = &A;
    B = B.t();
    arma::mat* pC = &B;

    // Create parameter objects
    reg_struct = Rcpp::as<Rcpp::List>(rs[1]);
    vparam.push_back(param(reg_struct, arma::zeros(pA->n_rows * pC->n_cols), pC->n_cols));
    reg_struct = Rcpp::as<Rcpp::List>(rs[2]);
    vparam.push_back(param(reg_struct, arma::zeros(pA->n_rows * pC->n_cols), pC->n_cols));

    ratmak ratmak_ (&vparam, pA, pC);
    field* pratmak_ = &ratmak_;
    pfield = pratmak_->deep_Clone();
  } else if (is_in("rlk", ode_struct.attr("class"))) {
    if (ode_struct["sparseA"]) {
      A = arma::conv_to<arma::mat>::from(Rcpp::as<arma::sp_mat>(ode_struct["A"]));
    } else {
      A = Rcpp::as<arma::mat>(ode_struct["A"]);
    }
    if (ode_struct["sparseB"]) {
      B = arma::conv_to<arma::mat>::from(Rcpp::as<arma::sp_mat>(ode_struct["B"]));
    } else {
      B = Rcpp::as<arma::mat>(ode_struct["B"]);
    }
    arma::mat* pA1 = &A;
    arma::mat* pA2 = &B;

    // Create parameter objects
    reg_struct = Rcpp::as<Rcpp::List>(rs[1]);
    vparam.push_back(param(reg_struct, arma::zeros(A.n_elem), d));

    rlk rlk_ (&vparam, pA1, pA2);
    field* prlk_ = &rlk_;
    pfield = prlk_->deep_Clone();
  } else {
    Rcpp::stop("ode not recognised!");
  }

  return pfield;
}




//=================================
// Evaluated field
//=================================

// [[Rcpp::export(ode_field)]]
Rcpp::List ode_field(Rcpp::List ode_struct, arma::vec x, Rcpp::Nullable<Rcpp::List> param_, bool differentials) {
  // Create ode pointer
  arma::mat A;
  arma::mat B;
  std::vector<param> vparam;

  // Rcpp::Rcout << "Calling create_pfield" << std::endl;
  field* pfield = create_pfield(ode_struct, A, B, vparam, x);

  // Fill in param values
  pfield->fill_in_param(param_);
  std::vector<param>::iterator param_it = vparam.begin();

  // Return list
  Rcpp::List ret;
  ret["f"] = pfield->f(0.0, x);
  if (differentials) {
    ret["f_dx"] = arma::sp_mat(pfield->f_dx(0.0, x));

    Rcpp::List params_(vparam.size() - 1);
    param_it = vparam.begin();
    for (unsigned int i = 0; i < params_.size(); i++) {
      param_it++;
      params_[i] = arma::sp_mat(pfield->f_dparam(0.0, x, param_it));
    }
    ret["f_dparam"] = params_;
  }

  return ret;
}


