# ifndef MOTOR_REG_H_
# define MOTOR_REG_H_

#include <RcppArmadillo.h>

/* Regularisation class. Holds virtual functions:
 *  - penalty
 *  - prox operator, which given x (current position), p (negative step direction)
 *                    and tau (backtrack scale) returns the next step
 *  - d2pen (double differential of penalty)
 *
 *  Future room for improvement:
 *  - if lower, upper or v are vectors of identicals, replace with double
 */
class reg {
public:
  arma::vec lower_full, upper_full, v_full; // These are not constant or private, since nlfit needs to modify
  arma::vec lower, upper, v;
  bool v_exist;
  double lambda;    // Will be changed from outside
  double alpha;     // Represents the alpha value used in fun::lambda_max

protected:
  reg(arma::vec lower_full_, arma::vec upper_full_, arma::vec v_full_, double lambda_);

public:
  virtual reg* deep_Clone() = 0;      // deep clone for specifying correct preg in raw_maker etc
  virtual reg* shallow_Clone() = 0;   // shallow copy for optim-constructors
  virtual ~reg();

  // Virtual penalty and prox-operator (penalty also includes multiplying with lambda)
  virtual double penalty(arma::vec position) = 0;
  virtual arma::vec prox(arma::vec position, arma::vec direction, double tau) = 0;
  virtual arma::vec d2pen(arma::vec position) = 0;

  // Finds those with increment at position or already non-zero position
  arma::uvec get_active(arma::vec position, arma::vec direction, double tau);

  // Clamps any vector of same length as u and l between these limits
  void box(arma::vec &x);

  // Reduce reg to only a subset of the parameters
  void reduce(arma::uvec active_);
  void reduce(void);
};



class l1 : public reg {
public:
  l1(arma::vec lower_full_, arma::vec upper_full_, arma::vec v_full_, double lambda_);
  l1* deep_Clone();
  l1* shallow_Clone();

  double penalty(arma::vec position);
  arma::vec prox(arma::vec position, arma::vec direction, double tau);
  arma::vec d2pen(arma::vec position);
};



class l2 : public reg {
public:
  l2(arma::vec lower_full_, arma::vec upper_full_, arma::vec v_full_, double lambda_);
  l2* deep_Clone();
  l2* shallow_Clone();

  double penalty(arma::vec position);
  arma::vec prox(arma::vec position, arma::vec direction, double tau);
  arma::vec d2pen(arma::vec position);
};



class elnet : public reg {
public:
  double a;

  elnet(arma::vec lower_full_, arma::vec upper_full_, arma::vec v_full_, double lambda_, double a_);
  elnet* deep_Clone();
  elnet* shallow_Clone();

  double penalty(arma::vec position);
  arma::vec prox(arma::vec position, arma::vec direction, double tau);
  arma::vec d2pen(arma::vec position);
};



class mcp : public reg {
  arma::vec p(arma::vec position);
  arma::vec pprime(arma::vec position);

public:
  double a;   // a > 1

  mcp(arma::vec lower_full_, arma::vec upper_full_, arma::vec v_full_, double lambda_, double a_);
  mcp* deep_Clone();
  mcp* shallow_Clone();

  double penalty(arma::vec position);
  arma::vec prox(arma::vec position, arma::vec direction, double tau);
  arma::vec d2pen(arma::vec position);
};



class scad : public reg {
  arma::vec p(arma::vec position);
  arma::vec pprime(arma::vec position);

public:
  double a; // a > 2

  scad(arma::vec lower_full_, arma::vec upper_full_, arma::vec v_full_, double lambda_, double a_);
  scad* deep_Clone();
  scad* shallow_Clone();

  double penalty(arma::vec position);
  arma::vec prox(arma::vec position, arma::vec direction, double tau);
  arma::vec d2pen(arma::vec position);
};



class none : public reg {
public:
  none(arma::vec lower_full_, arma::vec upper_full_, arma::vec v_full_, double lambda_);
  none* deep_Clone();
  none* shallow_Clone();

  double penalty(arma::vec position);
  arma::vec prox(arma::vec position, arma::vec direction, double tau);
  arma::vec d2pen(arma::vec position);
};

arma::vec stretch_lu(arma::vec v, unsigned int p);

reg* create_preg(std::string reg_type, arma::vec lower, arma::vec upper, arma::vec v, double lambda, double a);
reg* create_preg(Rcpp::List reg_struct, unsigned int p_full);

# endif
