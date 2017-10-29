# ifndef MOTOR_FUN_H_
# define MOTOR_FUN_H_

#include "motor_likeli.h"
#include "motor_fit.h"

/* Base (and abstract) class fun
 *
 * Is used for specifying how the vector of parameters (arguments)
 * are related to the value of the function
 */
class fun {
  double lambda;

public:
  // Current value of parameters
  std::vector<param> vparam;

  // Logical used below to check if gnhda should be evaluated in evaluate_differential
  bool eval_gnhda;

  // Current value of function
  double value;

protected:
  // Constructor and destructor
  fun (std::vector<param> vparam_);
  fun (const fun& other);
  virtual ~fun() {};

public:
  double get_lambda(void);
  void set_lambda(double lambda_);

  /* Virtual function parts, i.e. evaluate function
   *
   * evaluate(void) takes current parameters (vparam), evaluate fun in those and store result in value:
   *    value = fun(vparam)
   *
   * evaluate_differential(std::vector<param>::iterator it) updates 'differential' in said param
   *    it->differential = dfun(it)/dparam
   */
  virtual void evaluate(void) = 0;
  virtual void evaluate_differential(std::vector<param>::iterator it) = 0;

  /* Non-virtual function parts
   *
   * evaluate(std::vector<param>::iterator it, bool differential) evaluates fun and
   * adds penalty and optionally updates 'differential' in said param
   *    value = fun(vparam) + penalty(it)
   *    if (differential) it->differential = dfun(it)/dparam
   */
  void evaluate(std::vector<param>::iterator it, bool differential);


  /* Virtual set_active part
   *
   * set_active(it) reduces the function procedure for parameter it by it->active
   */
  virtual void set_active(std::vector<param>::iterator it) = 0;


  /* diff_at_zero
   * Evaluates differentials at position = 0 (unless they are fixed, then they stay the same)
   * return has same length as vparam, each entry has length p_full
   */
  std::vector<arma::vec> diff_at_zero(void);

  /* lambda_max
   * Uses the above to return a vector (length of vparam) with individual maximal lambda
   * if fixed, it returns 0
   */
  arma::vec lambda_max(void);

  /* get_and_set_lambda() uses fun to get lambda_max,
   * returns max(lambda_max) and scales lambda_factors by ratio
   * the stored lambda_factors are returned by reference through lamfac
   */
  double get_and_set_lambda_max(arma::vec &lamfac);
};


/* Derived class: additive model
 *
 * E(Y|X) = mu(eta) = mu(fit_1(param_1) + fit_2(param_2) + ...) (each of them penalised differently, think group lasso)
 * V(Y|X) = W * sigma^2
 *
 * 'likeli' defines the likelihood that connects Y, fit = sum fit_i(param_i) and W
 */
class am : public fun {
public:
  arma::vec W;
  bool W_exist;

  const arma::vec Y;
  const bool centerY;

private:
  likeli* plikeli;
  std::vector<fit*> vpfit;

public:
  am(arma::vec Y_, arma::vec W_, likeli* plikeli_, std::vector<param> vparam_, std::vector<fit*> vpfit_, bool centerY_);

  void evaluate(void);
  void evaluate_differential(std::vector<param>::iterator it);
  void set_active(std::vector<param>::iterator it);
};


// Misc
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
inline arma::urowvec operator&(const arma::urowvec& lhs, arma::urowvec& rhs){
  // Does bitwise AND for two urowvecs, used to find in range. Returns vector of length of lhs, recycles rhs
  arma::urowvec ret(arma::size(lhs), arma::fill::zeros);
  if (!lhs.is_empty()) {
    for (unsigned int i = 0; i < lhs.n_elem; ++i) {
      ret(i) = lhs(i) && rhs(i % rhs.n_elem);
    }
  }
  return ret;
}

inline arma::uvec operator|(arma::uvec lhs, arma::uvec rhs){
  // Does bitwise OR for two uvecs, used to find in range. Returns vector of length of lhs, recycles rhs
  arma::uvec ret(arma::size(lhs), arma::fill::zeros);
  if (!lhs.is_empty()) {
    for (unsigned int i = 0; i < lhs.n_elem; ++i) {
      ret(i) = lhs(i) || rhs(i % rhs.n_elem);
    }
  }
  return ret;
}
inline arma::urowvec operator|(const arma::urowvec& lhs, arma::urowvec& rhs){
  // Does bitwise OR for two urowvecs, used to find in range. Returns vector of length of lhs, recycles rhs
  arma::urowvec ret(arma::size(lhs), arma::fill::zeros);
  if (!lhs.is_empty()) {
    for (unsigned int i = 0; i < lhs.n_elem; ++i) {
      ret(i) = lhs(i) || rhs(i % rhs.n_elem);
    }
  }
  return ret;
}


# endif
