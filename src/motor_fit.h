# ifndef MOTOR_FIT_H_
# define MOTOR_FIT_H_

#include <RcppArmadillo.h>
#include "motor_param.h"

// fit class, defines the param->fit mapping only holds two virtual functions
class fit {
public:
  arma::vec value;            // n-dim vector (n is dim of Y and W)
  arma::mat value_dposition;  // n-x-p matrix of differential
  param* pparam;              // Points to its corresponding parameter in am-object
  arma::vec* pW;              // Pointer to W in am
  bool* pW_exist;             // Pointer to W_exist

protected:
  arma::vec pos_value;        // Position at which value is from (used to check if value should be updated)
  arma::vec pos_diff;         // Position at which value_dposition is from
  const bool center;          // Should fit be centered?

protected:
  fit(param* pparam_, bool center_, arma::vec* pW_, bool* pW_exist_, bool is_linear_);

private:
  const bool is_linear;

public:
  fit(const fit& other);
  fit& operator=(const fit& rhs);

  // Virtual destructor and clone
  virtual ~fit() {};
  // virtual fit* shallow_Clone() = 0;

  /* Evaluates fit for its param value (position)
   * They do so by checking if pos_value (or pos_diff) matches pparam->position,
   * if so they do nothing, else they transfer pparam->position to pos_value and
   * calls set_value
   */
  void evaluate(void);              // Has length n (length Y and W)
  void evaluate_differential(void); // Dim n-x-p

  /* Virtual functions
   * which specify how to set 'value' (they must use pos_value for this purpose)
   * set_active uses pparam->active for the job
   */
  virtual void set_value(void) = 0;
  virtual void set_value_dposition(void) = 0;
  virtual void set_active(void) = 0;

private:
  // Checks if a change is needed (i.e., if refpos is different from pparam->position)
  bool change(const arma::vec& refpos);
};



// linear fit subclass
class lfit : public fit {
  const arma::mat X_full;

public:
  arma::mat X;

  lfit(param* pparam_, bool center_, arma::vec* pW_, bool* pW_exist_, arma::mat X_full_);
  // lfit* shallow_Clone();

  void set_value(void);
  void set_value_dposition(void);
  void set_active(void);
};



/* nonlinear fit subclass, this is also abstract,
 * its purpose is to group all non linear fits and
 * do the generic centering operation on the value and differentials
 * returned from its subclasses
 */
class nlfit : public fit {
protected:
  nlfit(param* pparam_, bool center_, arma::vec* pW_, bool* pW_exist_);

public:
  void set_value(void);
  void set_value_dposition(void);

  virtual void set_active(void) = 0;
  virtual void nlfit_set_value(void) = 0;
  virtual void nlfit_set_value_dposition(void) = 0;
};


/* log linear nlfit subclass
 *      fit = log(ydot + y * position)
 * with ydot, y = n-dim, position 1-dim and its lower bound will be adjusted!
 */
class loglin : public nlfit {
  const arma::vec ydot;
  const arma::vec y;

public:
  // In this constructor we may modify lower/upper limits and box position
  loglin(param* pparam_, bool center_, arma::vec* pW_, bool* pW_exist_, arma::vec ydot_, arma::vec y_);

  void set_active(void);
  void nlfit_set_value(void);
  void nlfit_set_value_dposition(void);
};

/* hillA (activating hill function) nlfit subclass
 *      fit = sum_a(log(1 + c_a * x^a) - log(tol + c_a * x^a))
 * x = n-x-d-dim, position r-dim, A = (a)_a r-x-d dim
 */
class hillA : public nlfit {
  const arma::mat xA_full;    // x^A in n-x-r dim matrix
  arma::mat xA;

public:
  // In this constructor we may modify lower/upper limits and box position
  hillA(param* pparam_, bool center_, arma::vec* pW_, bool* pW_exist_, arma::mat x, arma::mat A);

  void set_active(void);
  void nlfit_set_value(void);
  void nlfit_set_value_dposition(void);
};

/* hillB (inhibiting hill function) nlfit subclass
 *      fit = sum_b(log(1 + c_b * x^b))
 * x = n-x-d-dim, position r-dim, B = (b)_b r-x-d dim
 */
class hillB : public nlfit {
  const arma::mat xB_full;    // x^B in n-x-r dim matrix
  arma::mat xB;

public:
  // In this constructor we may modify lower/upper limits and box position
  hillB(param* pparam_, bool center_, arma::vec* pW_, bool* pW_exist_, arma::mat x, arma::mat B);

  void set_active(void);
  void nlfit_set_value(void);
  void nlfit_set_value_dposition(void);
};


# endif
