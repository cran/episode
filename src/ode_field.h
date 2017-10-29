# ifndef ODE_FIELD_H_
# define ODE_FIELD_H_

#include <RcppArmadillo.h>

#include "motor_param.h"



/* sc (scale and context)
 * Designed for holding scales * context (column-wise product) for a given parameter
 * Has exist flag and a scale member function
 */
struct sc {
  bool exist;
  arma::mat value;

  // Constructor
  sc (Rcpp::Nullable<arma::mat> value_);

  // Returns a certain column context (series)
  arma::vec col(unsigned int context);

  // Transfer certain rows of sc to another sc instance
  void transfer(sc &other, arma::uvec rows_coor);
};


/* Calculates the matrix power x^P = (\Pi_{j}x_j^P_ij)_i
 * and the differential wrt x
 */
arma::vec matrix_power(arma::vec x, const arma::mat &P);
arma::mat dmatrix_power(arma::vec x, const arma::mat &P);



/* ode class is abstract and holds placeholder functions
 * f, f_dx, f_dparam are drift and differentials wrt parameters (x and param)
 * Reduce and restore are used to reduce the dimension of param.
 */
class field {
public:
  // Pointer to param vector (in ode_loss they must agree)
  std::vector<param>* pvparam;

  unsigned int d;     // Dimension of state
  const arma::uvec linear;  // Indicate if linear in param, has length one less than pvparam

  field (std::vector<param>* pvparam_, unsigned int d_, arma::uvec linear_);

  // Use default copy constructor, but create virtual Clone function for cloning derived classes correctly
  virtual field* shallow_Clone() = 0;
  virtual field* deep_Clone() = 0;      // deep clone for specifying correct pode in the exported functions
  virtual ~field() {}

  // Drift (pure virtual function)
  virtual arma::vec f (double t, arma::vec x) = 0;

  // State derivate (pure virtual function)
  virtual arma::mat f_dx (double t, arma::vec x) = 0;

  // Rate derivate (pure virtual function)
  virtual arma::mat f_dparam (double t, arma::vec x, std::vector<param>::iterator it) = 0;

  // For all three functions above, the vector param->position is of length p (the (potentially) reduced parameter vector)

  // Virtual active setting for the system (uses pvparam to extract active)
  virtual void set_active(std::vector<param>::iterator it) = 0;

  // Calls set_active on all params
  void set_active(void);

  // Function that fills in value from param_ list to pvparam
  void fill_in_param(Rcpp::Nullable<Rcpp::List> param_);
};


/* Subclass of field: mak (mass action kinetics) */
class mak : public field {
  const arma::mat *A, *C;
  arma::mat A_sub, C_sub;

public:
  mak(std::vector<param>* pvparam_, arma::mat* pA, arma::mat* pC);
  mak* shallow_Clone();
  mak* deep_Clone();
  ~mak();

  arma::vec f(double t, arma::vec x);
  arma::mat f_dx(double t, arma::vec x);
  arma::mat f_dparam(double t, arma::vec x, std::vector<param>::iterator it);

  void set_active(std::vector<param>::iterator it);
};



/* Subclass of field: plk (power law kinetics)
* It is like mak, but only A is specified, the remaining parameterisation arise as theta = (B-A)'diag(param)
*
* Currently all parameters are concatinated vector
*/
class plk : public field {
  const arma::mat *A;
  arma::mat A_sub;

public:
  plk(std::vector<param>* pvparam_, arma::mat *pA);
  plk* shallow_Clone();
  plk* deep_Clone();
  ~plk();

  arma::vec f(double t, arma::vec x);
  arma::mat f_dx(double t, arma::vec x);
  arma::mat f_dparam(double t, arma::vec x, std::vector<param>::iterator it);

  void set_active(std::vector<param>::iterator it);
};



/* Subclass of field: ratmak (rational mass action kinetics)
 * Instead of x^A as rate function we have rational functions:
 * f(x, K1, K2) = C (K1 x^A1 / (1 + K2 * x^A2))
 * where C dxr, A1/A2 b1xd/b2xd and K1 and K2 rxb1/rxb2 non-negative
 */
class ratmak : public field {
  const arma::mat *A, *C;
  arma::mat A1, A2, C_sub;

public:
  ratmak(std::vector<param>* pvparam_, arma::mat *pA, arma::mat *pC);
  ratmak* shallow_Clone();
  ratmak* deep_Clone();
  ~ratmak();

  arma::vec f(double t, arma::vec x);
  arma::mat f_dx(double t, arma::vec x);
  arma::mat f_dparam(double t, arma::vec x, std::vector<param>::iterator it);

  void set_active(std::vector<param>::iterator it);
};




/* Subclass of field: rlk (rational law kinetics)
 * f(x, theta) = theta (x^A1 / (1 + x^A2))
 * where theta dxr (estimable), A1 and A2 rxd
 */
class rlk : public field {
  const arma::mat *A1, *A2;
  arma::mat A1_sub, A2_sub;

public:
  rlk(std::vector<param>* pvparam_, arma::mat* pA1, arma::mat* pA2);
  rlk* shallow_Clone();
  rlk* deep_Clone();
  ~rlk();

  arma::vec f(double t, arma::vec x);
  arma::mat f_dx(double t, arma::vec x);
  arma::mat f_dparam(double t, arma::vec x, std::vector<param>::iterator it);

  void set_active(std::vector<param>::iterator it);
};




// template <typename T>
// inline bool is_in(T x, std::vector<T> v)  {
//   return std::find(v.begin(), v.end(), x) != v.end();
// }
inline bool is_in(std::string name, std::vector<std::string> names) {
  return std::find(names.begin(), names.end(), name) != names.end();
}

field* create_pfield(Rcpp::List ode_struct, arma::mat &A, arma::mat &B, std::vector<param> &vparam, arma::mat x0);

/*
 * Assumption about field subclasses: if params are all 0, then f = 0, else lambda_max not true
 */


# endif
