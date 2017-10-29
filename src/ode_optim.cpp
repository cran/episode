//=============================================
/* ode_optim.cpp
 *
 * Content:
 * - bull for rodeo.ode
 * - cattle for aim
 * - bronc for rodeo.aim
 *
 *
 */
//=============================================

#include "ode_loss.h"
#include "motor_optim.h"



//=====================================
/* prepare_lambda
 * checks if lambda exist in opt,
 * if so it returns that, else
 * it produces one from lambda_ratio
 * also fills up lamfac with lambda_factors
 */
//=====================================

arma::vec prepare_lambda(fun* pfun, Rcpp::List opt_struct, arma::vec &lamfac) {
  arma::vec ret;

  // Set lambda_factors according to maximal lambdas and sends lambda_factors back by reference
  double lam_max = pfun->get_and_set_lambda_max(lamfac);

  Rcpp::Nullable<arma::vec> maybe_lambda = Rcpp::as< Rcpp::Nullable<arma::vec> >(opt_struct["lambda"]);
  if (maybe_lambda.isNotNull()) {
    ret = Rcpp::as<arma::vec>(maybe_lambda);
  } else {
    if (lam_max <= 0) {
      lam_max = 1e-7;
      Rcpp::warning("Automatically generated lambda sequence may not be appropriate.");
    }

    ret = lam_max * arma::exp(arma::linspace<arma::vec>(std::log(1),
      std::log(Rcpp::as<double>(opt_struct["lambda_min_ratio"])),
      Rcpp::as<unsigned int>(opt_struct["nlambda"])));
  }

  if (!Rcpp::as<bool>(opt_struct["lambda_decrease"])) {
    ret = arma::flipud(ret);
  }

  return ret;
}



//=====================================
// bull for rodeo.ode
//=====================================

// [[Rcpp::export(bull)]]
Rcpp::List bull(Rcpp::List ode_struct, Rcpp::List opt_struct, Rcpp::List sc_,
  Rcpp::List params_, arma::cube x0s, bool trace = false) {
  Rcpp::List ret;

  /*
   * Size check
   */
  if (params_.size() != x0s.n_slices) {
    Rcpp::stop("error in bull -> number of param list and slices in x0s does not match.");
  }
  unsigned int N = params_.size();
  // Rcpp::Rcout << "N = " << N << std::endl;
  Rcpp::List ctrl_ = Rcpp::as<Rcpp::List>(opt_struct["ctrl"]);
  double tol = Rcpp::as<double>(ctrl_["tol_l"]);


  /*
   * Prepare likeli (only holds functions)
   */
  gaussian gg;
  likeli* plikeli_ = &gg;


  /*
   * Prepare sc (same for all optims)
   */
  std::vector<sc> vsc_full_;
  for (unsigned int i = 0; i < sc_.size(); i++) {
    vsc_full_.push_back(sc(Rcpp::as< Rcpp::Nullable<arma::mat> >(sc_[i])));
  }
  // Rcpp::Rcout << "Prepared sc" << std::endl;


  if (trace) Rcpp::Rcout << "Defining ODE system." << std::endl;
  arma::mat A;  // The same for all optims
  arma::mat B;


  /*
   * Functions to optimise
   */
  std::vector<exact> vodeloss;
  for (unsigned int i = 0; i < N; i++) {
    if (trace) Rcpp::Rcout << "Defining function to optimise." << std::endl;
    exact exact_ = exact(ode_struct, A, B, x0s.slice(i), Rcpp::as< Rcpp::Nullable<Rcpp::List> >(params_[i]), opt_struct, &vsc_full_, plikeli_);
    // Rcpp::List par_list = Rcpp::as<Rcpp::List>(params_[i]);
    // Rcpp::CharacterVector cc = par_list.names();
    // Rcpp::Rcout << "par_list.names() = " << cc << std::endl;
    // Rcpp::Rcout << "par_list[k] = " << Rcpp::as<arma::vec>(par_list["k"]) << std::endl;

    vodeloss.push_back(exact_);
    // bool bb = &(*exact_.pfield->pvparam->begin()) == &(*exact_.vparam.begin());
    // Rcpp::Rcout << "Are pfield and vparam agreeing: " << bb << std::endl;
    // Rcpp::Rcout << "push back" << std::endl;
  }
  // Rcpp::Rcout << "exact param size" << vodeloss[0].vparam.size() << std::endl;

  unsigned int n_param = vodeloss.begin()->vparam.size();
  std::vector<exact>::iterator loss_it = vodeloss.begin();
  // Rcpp::Rcout << "Fun value:" << loss_it->value << std::endl;
  // std::vector<param>::iterator p_it = loss_it->vparam.begin();
  // Rcpp::Rcout << " x0:" << p_it->position.t() << std::endl;
  // p_it++;
  // Rcpp::Rcout << " k:" << p_it->position.t() << std::endl;
  // loss_it->evaluate();
  // Rcpp::Rcout << "Fun value:" << loss_it->value << std::endl;
  // Rcpp::Rcout << "position: " << loss_it->vparam.begin()->position.t() << std::endl <<
  //   "position_full: " << loss_it->vparam.begin()->position_full.t() << std::endl <<
  //     "active: " << loss_it->vparam.begin()->active.t() << std::endl <<
  //       "get_position" << loss_it->vparam.begin()->get_position() << std::endl;


  /*
   * Optimisation problems
   */
  std::vector<optim> voptim;
  for (unsigned int i = 0; i < N; i++) {
    if (trace) Rcpp::Rcout << "Defining optimisation problem." << std::endl;
    voptim.push_back(optim(&(*loss_it), trace, tol, "gd"));
    voptim.back().id = i + 1;

    loss_it++;
  }
  // std::vector<optim>::iterator optim_it = voptim.begin();


  /*
   * Lambdas
   */
  arma::vec lamfac;
  arma::vec lambda = prepare_lambda(&(*vodeloss.begin()), opt_struct, lamfac);
  unsigned int n_lambda = lambda.n_elem;
  // Rcpp::Rcout << "Setting up return object." << std::endl;


  /*
   * Return objects
   */
  std::vector<arma::cube> pars;
  std::vector<param>::iterator param_it;
  for (param_it = vodeloss.begin()->vparam.begin(); param_it < vodeloss.begin()->vparam.end(); param_it++) {
    // row = param coordinate, col = lambda, slice = particle
    pars.push_back(arma::cube(param_it->p_full, n_lambda, N, arma::fill::zeros));
  }
  std::vector<arma::cube>::iterator pars_it = pars.begin();
  // Rcpp::Rcout << "Sat up return object." << std::endl;

  arma::cube dfs(n_param, n_lambda, N);   dfs.fill(arma::datum::nan);
  arma::cube codes(n_param, n_lambda, N); codes.fill(arma::datum::nan);
  arma::cube steps(n_param, n_lambda, N); steps.fill(arma::datum::nan);
  arma::mat losses(N, n_lambda);          losses.fill(arma::datum::nan);
  arma::cube pens(n_param, n_lambda, N);  pens.fill(arma::datum::nan);
  arma::cube jerr(n_param, n_lambda, N);  jerr.fill(arma::datum::nan);
  // Rcpp::Rcout << "Sat up cubes." << std::endl;


  /*
   * Fixed (used to check if dfs should be calculated for that parameter, since param->fixed may change during optim)
   */
  arma::uvec fixed_param(n_param);
  for (param_it = vodeloss.begin()->vparam.begin(); param_it < vodeloss.begin()->vparam.end(); param_it++) {
    fixed_param(param_it - vodeloss.begin()->vparam.begin()) = param_it->fixed;
  }
  // Rcpp::Rcout << "fixed_param " << fixed_param.t() << std::endl;


  /*
   * Algorithm
   */
  if (trace) Rcpp::Rcout << "Commencing optimisation." << std::endl;
  for (unsigned int i_lambda = 0; i_lambda < n_lambda; ++ i_lambda) {
    /* Resets backtracking parameters and
     * current optimal loss (thus after set_lambda),
     * then moves the particles ahead
     */
    if (trace) Rcpp::Rcout << "Setting lambda number: " << i_lambda + 1 << " (" << lambda(i_lambda) << ")" << std::endl;
    set_lambda(voptim, lambda(i_lambda));

    if (trace) Rcpp::Rcout << "Optimising." << std::endl;
    move(voptim);


    /*
     * Retrieve results
     */
    if (trace) Rcpp::Rcout << "Retrieving results." << std::endl;
    for (loss_it = vodeloss.begin(); loss_it != vodeloss.end(); loss_it++) {

      /*
       * Loop over its parameters
       */
      pars_it = pars.begin();
      for (param_it = loss_it->vparam.begin(); param_it < loss_it->vparam.end(); param_it++) {
        /*
         * Fill in estimates
         */
        pars_it->subcube(0, i_lambda, loss_it - vodeloss.begin(),
          pars_it->n_rows - 1, i_lambda, loss_it - vodeloss.begin()) = param_it->position_full;

        /*
         * Degrees of freedom
         */
        // Rcpp::Rcout << "dfs --> is it fixed? " << fixed_param(param_it - loss_it->vparam.begin()) << std::endl;
        if (fixed_param(param_it - loss_it->vparam.begin()) == 1) {
          dfs(param_it - loss_it->vparam.begin(), i_lambda, loss_it - vodeloss.begin()) = 0;
        } else {
          // Rcpp::Rcout << "Size of get_X4dfs. n_cols " << loss_it->get_X4dfs(param_it).n_cols << std::endl;
          dfs(param_it - loss_it->vparam.begin(), i_lambda, loss_it - vodeloss.begin()) = param_it->dfs(loss_it->get_X4dfs(param_it));
        }

        /*
         * Convergence codes
         */
        if (loss_it->psolver->code_conv > 0) {
          codes(param_it - loss_it->vparam.begin(), i_lambda, loss_it - vodeloss.begin()) = 4;
        } else if (param_it->c.diff_code > 0) {
          codes(param_it - loss_it->vparam.begin(), i_lambda, loss_it - vodeloss.begin()) = 3;
        } else if (param_it->c.step >= param_it->c.step_max) {
          codes(param_it - loss_it->vparam.begin(), i_lambda, loss_it - vodeloss.begin()) = 2;
        } else if (param_it->c.backtrack_code > 0) {
          codes(param_it - loss_it->vparam.begin(), i_lambda, loss_it - vodeloss.begin()) = 1;
        } else {
          codes(param_it - loss_it->vparam.begin(), i_lambda, loss_it - vodeloss.begin()) = 0;
        }

        /*
         * Other diagnostics
         */
        steps(param_it - loss_it->vparam.begin(), i_lambda, loss_it - vodeloss.begin()) = param_it->c.step;
        jerr(param_it - loss_it->vparam.begin(), i_lambda, loss_it - vodeloss.begin())  = param_it->c.backtrack_code + 10 * (param_it->c.step >= param_it->c.step_max) + 100 * param_it->c.diff_code + 1000 * loss_it->psolver->code_conv;
        pens(param_it - loss_it->vparam.begin(), i_lambda, loss_it - vodeloss.begin())  = param_it->preg->penalty(param_it->position_full);

        pars_it++;
      }

      loss_it->evaluate();
      losses(loss_it - vodeloss.begin(), i_lambda) = loss_it->value;
    }
  }

  /*
   * Fill in return list
   */
  if (trace) Rcpp::Rcout << "Finished optimisation. Preparing results." << std::endl;
  Rcpp::List particle_list(N);
  for (unsigned int i = 0; i < N; ++i) {
    Rcpp::List param_list(pars.size());
    for (pars_it = pars.begin(); pars_it < pars.end(); pars_it++) {
      param_list[pars_it - pars.begin()] = arma::conv_to<arma::sp_mat>::from(pars_it->slice(i));
    }
    particle_list[i] = param_list;
  }
  ret["params"]     = particle_list;
  ret["dfs"]        = dfs;
  ret["codes"]      = codes;
  ret["steps"]      = steps;
  ret["losses"]     = losses;
  ret["penalties"]  = pens;
  ret["jerr"]       = jerr;
  ret["lambda"]     = lambda;
  ret["lambda_fac"] = lamfac;

  // std::vector<arma::vec> diffs0 = vodeloss.begin()->diff_at_zero();
  // Rcpp::Rcout << "diffs0 x0 = " << diffs0.begin()->t() << std::endl;
  // Rcpp::Rcout << "diffs0 rate = " << (diffs0.begin() + 1)->t() << std::endl;
  // arma::vec lam_max = vodeloss.begin()->lambda_max();
  // Rcpp::Rcout << "lambda_max = " << lam_max.t() << std::endl;

  return ret;
}





//=====================================
/* cattle for aim
 *
 * difference from above:
 * params is nullable, and
 * we only have one initialisation
 * so it is only a list of vectors
 *
 * x0s -> x
 */
//=====================================

// [[Rcpp::export(cattle)]]
Rcpp::List cattle(Rcpp::List ode_struct, Rcpp::List opt_struct, Rcpp::List sc_,
  Rcpp::Nullable<Rcpp::List> param_, arma::mat x, bool trace = false) {
  Rcpp::List ret;

  /*
  * Get control list
  */
  Rcpp::List ctrl_ = Rcpp::as<Rcpp::List>(opt_struct["ctrl"]);
  double tol = Rcpp::as<double>(ctrl_["tol_l"]);


  /*
  * Prepare likeli (only holds functions)
  */
  gaussian gg;
  likeli* plikeli_ = &gg;
  // Rcpp::Rcout << "Likeli test: " << plikeli_->eval(x.col(0), arma::ones(x.n_rows)) << ". Prepare sc" << std::endl;


  /*
  * Prepare sc (same for all optims)
  */
  std::vector<sc> vsc_full_;
  for (unsigned int i = 0; i < sc_.size(); i++) {
    vsc_full_.push_back(sc(Rcpp::as< Rcpp::Nullable<arma::mat> >(sc_[i])));
  }
  // Rcpp::Rcout << "Prepared sc" << std::endl;


  if (trace) Rcpp::Rcout << "Defining ODE system." << std::endl;
  arma::mat A;  // The same for all optims
  arma::mat B;
  // Rcpp::Rcout << "in cattle -> param_.isNotNull() = " << param_.isNotNull() << std::endl;


  /*
  * Function to optimise
  */
  if (trace) Rcpp::Rcout << "Defining function to optimise." << std::endl;
  unsigned int n_series = Rcpp::as<unsigned int>(opt_struct["s"]);
  arma::mat x0 = x.head_rows(n_series);
  x0 = x0.cols(1, x.n_cols - 1);
  integral integral_(ode_struct, A, B, x0.t(), param_, opt_struct, &vsc_full_, plikeli_, x);
  // Note x.head_rows() should represent x0, but in matching constructor it is overloaded with X0 (made from xout)
  // Rcpp::Rcout << "exact param size" << vodeloss[0].vparam.size() << std::endl;
  // Rcpp::Rcout << "in cattle -> created integral. isNotNull() = " << param_.isNotNull() << std::endl;


  unsigned int n_param = integral_.vparam.size();
  // Rcpp::Rcout << "Fun value:" << loss_it->value << std::endl;
  // std::vector<param>::iterator p_it = loss_it->vparam.begin();
  // Rcpp::Rcout << " x0:" << p_it->position.t() << std::endl;
  // p_it++;
  // Rcpp::Rcout << " k:" << p_it->position.t() << std::endl;
  // loss_it->evaluate();
  // Rcpp::Rcout << "Fun value:" << loss_it->value << std::endl;
  // Rcpp::Rcout << "position: " << loss_it->vparam.begin()->position.t() << std::endl <<
  //   "position_full: " << loss_it->vparam.begin()->position_full.t() << std::endl <<
  //     "active: " << loss_it->vparam.begin()->active.t() << std::endl <<
  //       "get_position" << loss_it->vparam.begin()->get_position() << std::endl;


  /*
  * Optimisation problem
  */
  if (trace) Rcpp::Rcout << "Defining optimisation problem." << std::endl;
  optim optim_(&integral_, trace, tol, "gd");
  optim_.id = 1;


  /*
  * Lambdas
  */
  // arma::vec lambda = Rcpp::as<arma::vec>(opt_struct["lambda"]);
  arma::vec lamfac;
  arma::vec lambda = prepare_lambda(&integral_, opt_struct, lamfac);
  unsigned int n_lambda = lambda.n_elem;
  // Rcpp::Rcout << "Setting up return object." << std::endl;


  /*
  * Return objects
  */
  std::vector<arma::mat> pars;
  std::vector<param>::iterator param_it;
  for (param_it = integral_.vparam.begin(); param_it < integral_.vparam.end(); param_it++) {
    // row = param coordinate, col = lambda
    pars.push_back(arma::mat(param_it->p_full, n_lambda, arma::fill::zeros));
  }
  std::vector<arma::mat>::iterator pars_it = pars.begin();
  // Rcpp::Rcout << "Sat up return object." << std::endl;

  arma::mat dfs(n_param, n_lambda);   dfs.fill(arma::datum::nan);
  arma::mat codes(n_param, n_lambda); codes.fill(arma::datum::nan);
  arma::mat steps(n_param, n_lambda); steps.fill(arma::datum::nan);
  arma::vec losses(n_lambda);         losses.fill(arma::datum::nan);
  arma::mat pens(n_param, n_lambda);  pens.fill(arma::datum::nan);
  arma::mat jerr(n_param, n_lambda);  jerr.fill(arma::datum::nan);
  // Rcpp::Rcout << "Sat up cubes." << std::endl;


  /*
  * Fixed (used to check if dfs should be calculated for that parameter, since param->fixed may change during optim)
  */
  arma::uvec fixed_param(n_param);
  for (param_it = integral_.vparam.begin(); param_it < integral_.vparam.end(); param_it++) {
    fixed_param(param_it - integral_.vparam.begin()) = param_it->fixed;
  }
  // Rcpp::Rcout << "fixed_param " << fixed_param.t() << std::endl;


  /*
  * Algorithm
  */
  if (trace) Rcpp::Rcout << "Commencing optimisation." << std::endl;
  for (unsigned int i_lambda = 0; i_lambda < n_lambda; ++ i_lambda) {
    /* Resets backtracking parameters and
    * current optimal loss (thus after set_lambda),
    * then moves the particles ahead
    */
    if (trace) Rcpp::Rcout << "Setting lambda number: " << i_lambda + 1 << " (" << lambda(i_lambda) << ")" << std::endl;
    optim_.set_lambda(lambda(i_lambda));

    if (trace) Rcpp::Rcout << "Optimising." << std::endl;
    optim_.reset_optim_param();
    optim_.optimise();


    if (trace) Rcpp::Rcout << "Retrieving results." << std::endl;
    /*
    * Loop over parameters
    */
    pars_it = pars.begin();
    for (param_it = integral_.vparam.begin(); param_it < integral_.vparam.end(); param_it++) {
      /*
      * Fill in estimates
      */
      pars_it->col(i_lambda) = param_it->position_full;

      /*
      * Degrees of freedom
      */
      // Rcpp::Rcout << "dfs --> is it fixed? " << fixed_param(param_it - loss_it->vparam.begin()) << std::endl;
      if (fixed_param(param_it - integral_.vparam.begin()) == 1) {
        dfs(param_it - integral_.vparam.begin(), i_lambda) = 0;
      } else {
        // Rcpp::Rcout << "Size of get_X4dfs. n_cols " << loss_it->get_X4dfs(param_it).n_cols << std::endl;
        dfs(param_it - integral_.vparam.begin(), i_lambda) = param_it->dfs(integral_.get_X4dfs(param_it));
      }

      /*
      * Convergence codes
      */
      if (param_it->c.diff_code > 0) {
        codes(param_it - integral_.vparam.begin(), i_lambda) = 3;
      } else if (param_it->c.step >= param_it->c.step_max) {
        codes(param_it - integral_.vparam.begin(), i_lambda) = 2;
      } else if (param_it->c.backtrack_code > 0) {
        codes(param_it - integral_.vparam.begin(), i_lambda) = 1;
      } else {
        codes(param_it - integral_.vparam.begin(), i_lambda) = 0;
      }

      /*
      * Other diagnostics
      */
      steps(param_it - integral_.vparam.begin(), i_lambda) = param_it->c.step;
      jerr(param_it - integral_.vparam.begin(), i_lambda)  = param_it->c.backtrack_code + 10 * (param_it->c.step >= param_it->c.step_max) + 100 * param_it->c.diff_code;
      pens(param_it - integral_.vparam.begin(), i_lambda)  = param_it->preg->penalty(param_it->position_full);

      pars_it++;
    }

    integral_.evaluate();
    losses(i_lambda) = integral_.value;
  }

  /*
   * Fill in return list
   */
  if (trace) Rcpp::Rcout << "Finished optimisation. Preparing results." << std::endl;
  Rcpp::List param_list(pars.size());
  for (pars_it = pars.begin(); pars_it < pars.end(); pars_it++) {
    param_list[pars_it - pars.begin()] = arma::conv_to<arma::sp_mat>::from(*pars_it);
  }
  ret["params"]     = param_list;
  ret["dfs"]        = dfs;
  ret["codes"]      = codes;
  ret["steps"]      = steps;
  ret["losses"]     = losses;
  ret["penalties"]  = pens;
  ret["jerr"]       = jerr;
  ret["lambda"]     = lambda;
  ret["lambda_fac"] = lamfac;

  return ret;
}






//=====================================
// bronc for rodeo.aim
//=====================================

// [[Rcpp::export(bronc)]]
Rcpp::List bronc(Rcpp::List ode_struct, Rcpp::List opt_struct, Rcpp::List sc_,
  Rcpp::List params_, arma::cube x0s, arma::uvec indices,
  bool adjust_lambda, bool adjust_scales, bool adjust_weights, bool trace = false) {
  Rcpp::List ret;

  /* Format of sp_params:
   * list of lists (containing parameter vectors), so each of these lists are
   * compatible with field::fill_in_param
   */

  /*
  * Size check and prepare indices
  */
  Rcpp::List ctrl_ = Rcpp::as<Rcpp::List>(opt_struct["ctrl"]);
  double tol = Rcpp::as<double>(ctrl_["tol_l"]);
  if (params_.size() != x0s.n_slices) {
    Rcpp::stop("error in bronc -> number of param list and slices in x0s does not match.");
  }
  unsigned int N = params_.size();
  indices -= 1;  // So that it is compatible with C++-notation
  if (indices.n_elem != N) {
    Rcpp::stop("error in bronc -> number of param list and indices do not match.");
  }
  // Rcpp::Rcout << "N = " << N << std::endl;


  /*
  * Prepare likeli (only holds functions)
  */
  gaussian gg;
  likeli* plikeli_ = &gg;
  // Rcpp::Rcout << "Likeli test: " << plikeli_->eval(x.col(0), arma::ones(x.n_rows)) << ". Prepare sc" << std::endl;


  /*
  * Prepare sc (same for all optims)
  */
  std::vector<sc> vsc_full_;
  for (unsigned int i = 0; i < sc_.size(); i++) {
    vsc_full_.push_back(sc(Rcpp::as< Rcpp::Nullable<arma::mat> >(sc_[i])));
  }
  unsigned int n_param = sc_.size() + 1;
  // Rcpp::Rcout << "Prepared sc" << std::endl;


  /*
  * Function used for setting lambda
  */
  if (trace) Rcpp::Rcout << "Defining ODE system." << std::endl;
  arma::mat A;  // The same for all optims
  arma::mat B;
  exact exact_ = exact(ode_struct, A, B, x0s.slice(0), Rcpp::as< Rcpp::Nullable<Rcpp::List> >(params_[0]), opt_struct, &vsc_full_, plikeli_);


  /*
   * Lambdas
   */
  arma::vec lamfac;
  arma::vec lambda;
  if (adjust_lambda) {
    // If we are to adjust lambda, we set the one in opt_struct to null
    opt_struct["lambda"] = R_NilValue;
    if (trace) Rcpp::Rcout << "Preparing lambda sequence." << std::endl;
  }
  lambda = prepare_lambda(&exact_, opt_struct, lamfac);
  unsigned int n_lambda = lambda.n_elem;


  /*
   * Functions to optimise and forward sweep
   */
  std::vector<param>::iterator param_it;
  std::vector<exact> vodeloss;
  std::vector<exact>::iterator loss_it;
  std::vector<optim> voptim;
  std::vector<optim>::iterator optim_it;
  if (trace) Rcpp::Rcout << "Starting forward sweep." << std::endl;
  unsigned int i_indices = 0;
  for (unsigned int i_lambda = 0; i_lambda < n_lambda; i_lambda++) {
    // If i_lambda is in the indices include and move index, else move lambda,
    while(i_lambda == indices(i_indices)) {
      if (trace) Rcpp::Rcout << "Setting lambda number: " << i_lambda + 1 << " (" << lambda(i_lambda) << ")" << std::endl;
      exact_ = exact(ode_struct, A, B, x0s.slice(i_indices), Rcpp::as< Rcpp::Nullable<Rcpp::List> >(params_[i_indices]), opt_struct, &vsc_full_, plikeli_);
      vodeloss.push_back(exact_);
      voptim.push_back(optim(&(vodeloss.back()), trace, tol, "gd"));
      voptim.back().id = i_indices + 1;
      // voptim.back().trace = false;

      if (i_indices < N - 1) {
        ++i_indices;
      } else {
        break;
      }
    }

    set_lambda(voptim, lambda(i_lambda));

    if (trace) Rcpp::Rcout << "Optimising." << std::endl;
    move(voptim);
  }
  if (trace) Rcpp::Rcout << "Finished forward sweep." << std::endl;
  if (n_param != exact_.vparam.size()) {
    Rcpp::Rcout << "number of scales: " << n_param - 1 << " and number of parameters: " << exact_.vparam.size() << std::endl;
    Rcpp::stop("error bronc -> Number of scales must be one less than number of parameters.");
  }
  // Rcpp::Rcout << "vodeloss.begin()->vparam.begin()->position " << vodeloss.begin()->vparam.begin()->position.t() <<
  //   std::endl << "(vodeloss.begin()->vparam.begin() + 1)->position " << (vodeloss.begin()->vparam.begin() + 1)->position.t() <<
  //     std::endl;


  /*
   * Adjust weights
   */
  if (adjust_weights) {
    if (trace) Rcpp::Rcout << "Adjusting weights." << std::endl;

    /*
     * Get residual sum of squares and degrees of freedom
     */
    arma::mat rss(exact_.pfield->d, N);
    arma::vec dfs = arma::zeros(N);

    for (loss_it = vodeloss.begin(); loss_it < vodeloss.end(); loss_it++) {
      // squared residuals stored colum-wise
      arma::mat sqres = arma::square(loss_it->y - arma::trans(loss_it->psolver->solve(loss_it->pfield, loss_it->time)));
      if (loss_it->w_exist) {
        rss.col(loss_it - vodeloss.begin()) = arma::sum(loss_it->w % sqres, 1);
      } else {
        rss.col(loss_it - vodeloss.begin()) = arma::sum(sqres, 1);
      }

      // sum of df for each particle
      for (param_it = loss_it->vparam.begin(); param_it < loss_it->vparam.end(); param_it++) {
        // Rcpp::Rcout << "X = " << loss_it->get_X4dfs(param_it) << std::endl;
        double add_dfs = param_it->dfs(loss_it->get_X4dfs(param_it));
        if (!std::isnan(add_dfs)) dfs(loss_it - vodeloss.begin()) += add_dfs;
      }
    }


    /*
     * Reweighting part
     */
    // Rcpp::Rcout << "rss = " << rss.t() << std::endl;
    // Rcpp::Rcout << "dfs = " << dfs.t() << std::endl;
    arma::vec rec_gcv = arma::square(1 - arma::clamp(dfs / (exact_.pfield->d * exact_.n_obs), 0, 1)) / arma::sum(rss, 0).t();
    arma::vec w_const = 1 / arma::clamp(rss * (rec_gcv / arma::clamp(exact_.pfield->d * (exact_.n_obs - exact_.n_series) - dfs, 1, arma::datum::inf)), tol, arma::datum::inf);
    w_const /= arma::mean(w_const);
    // Rcpp::Rcout << "rec_gcv = " << rec_gcv.t() << std::endl;
    // Rcpp::Rcout << "w_const = " << w_const.t() << std::endl;

    for (loss_it = vodeloss.begin(); loss_it < vodeloss.end(); loss_it++) {
      // If it does not exist, it does now
      if (!loss_it->w_exist) {
        loss_it->w_exist = true;
        loss_it->w = arma::ones(arma::size(loss_it->y));
      }

      // Mean weight for each column (observation)
      arma::rowvec obsMean   = arma::mean(loss_it->w, 0);
      loss_it->w.each_col() %= w_const;                             // Rescale
      loss_it->w.each_row() %= obsMean / arma::clamp(arma::mean(loss_it->w, 0), tol, arma::datum::inf); // Obsmean remain unchanged
    }

    ret["new_weights"] = vodeloss.begin()->w.t();
  }


  /*
   * Adjust scales
   */
  std::vector<arma::vec> sc_mults(n_param - 1);
  std::vector<arma::vec>::iterator sc_mult_it = sc_mults.begin();
  if (adjust_scales) {
    if (trace) Rcpp::Rcout << "Adjusting scales." << std::endl;

    // Loop over each parameter (but not x0)
    for (unsigned int i_param = 1; i_param < n_param; i_param++) {
      // Extract gnhda for new scale (note that p == p_full, since each optim::optimise call returns to full system)
      arma::mat gnhdas = arma::zeros((exact_.vparam.begin() + i_param)->p, N);
      for (loss_it = vodeloss.begin(); loss_it < vodeloss.end(); loss_it++) {
        if (loss_it->eval_gnhda) {
          // Then already available
        } else {
          loss_it->eval_gnhda = true;
          loss_it->evaluate_differential(loss_it->vparam.begin() + i_param);
          gnhdas.col(loss_it - vodeloss.begin()) = (loss_it->vparam.begin() + i_param)->gnhda;
          loss_it->eval_gnhda = false;
        }
      }
      // Rcpp::Rcout << "Got gnhdas" << std::endl;


      // Standardise and take care of gnhdas
      for (arma::mat::iterator it = gnhdas.begin(); it < gnhdas.end(); it++) {
        if (*it <= 0) {
          *it = 0;
        } else {
          *it = 1 / (*it);
        }
      }
      double p_tol = (exact_.vparam.begin() + i_param)->c.tol / 10.0;  // Tolerance of parameter in question
      gnhdas = arma::sqrt(gnhdas);
      gnhdas.each_row() /= arma::clamp(arma::mean(gnhdas, 0), p_tol, arma::datum::inf);   // Standardise each particle
      arma::vec sc_mult = arma::mean(gnhdas, 1);    // What to scale with (average over particles)
      sc_mult /= std::max(arma::mean(sc_mult), p_tol);
      // Rcpp::Rcout << "Got sc_mult" << std::endl;

      // Adjust sc
      if ((vsc_full_.begin() + i_param - 1)->exist) {
        // Will preserve colmeans
        arma::rowvec colMean = arma::mean((vsc_full_.begin() + i_param - 1)->value, 0);
        (vsc_full_.begin() + i_param - 1)->value.each_col() %= sc_mult;
        (vsc_full_.begin() + i_param - 1)->value.each_row() %= colMean / arma::clamp(arma::mean((vsc_full_.begin() + i_param - 1)->value, 0), p_tol, arma::datum::inf);
        // Now they are preserved
      } else {
        // If not exist, create it
        (vsc_full_.begin() + i_param - 1)->exist = true;
        (vsc_full_.begin() + i_param - 1)->value = sc_mult * arma::ones<arma::rowvec>(exact_.n_series);
      }

      // Update sc in vodeloss (transfers vsc_full (changed through vsc_full_) to vsc)
      for (loss_it = vodeloss.begin(); loss_it < vodeloss.end(); loss_it++) {
        loss_it->set_active_base(loss_it->vparam.begin() + i_param);
      }

      *sc_mult_it = sc_mult;
      sc_mult_it++;
    }
  }


  /*
   * Adjust lambda
   */
  if (adjust_lambda) {
    if (trace) Rcpp::Rcout << "Adjusting lambda." << std::endl;
    lambda = prepare_lambda(&(*vodeloss.begin()), opt_struct, lamfac);
  }


  /*
  * Return objects
  */
  std::vector<arma::cube> pars;
  for (param_it = vodeloss.begin()->vparam.begin(); param_it < vodeloss.begin()->vparam.end(); param_it++) {
    // row = param coordinate, col = lambda, slice = particle
    pars.push_back(arma::cube(param_it->p_full, n_lambda, N, arma::fill::zeros));
  }
  std::vector<arma::cube>::iterator pars_it = pars.begin();
  // Rcpp::Rcout << "Sat up return object." << std::endl;

  arma::cube dfs(n_param, n_lambda, N);   dfs.fill(arma::datum::nan);
  arma::cube codes(n_param, n_lambda, N); codes.fill(arma::datum::nan);
  arma::cube steps(n_param, n_lambda, N); steps.fill(arma::datum::nan);
  arma::mat losses(N, n_lambda);          losses.fill(arma::datum::nan);
  arma::cube pens(n_param, n_lambda, N);  pens.fill(arma::datum::nan);
  arma::cube jerr(n_param, n_lambda, N);  jerr.fill(arma::datum::nan);
  // Rcpp::Rcout << "Sat up cubes." << std::endl;


  /*
  * Fixed (used to check if dfs should be calculated for that parameter, since param->fixed may change during optim)
  */
  arma::uvec fixed_param(n_param);
  for (param_it = vodeloss.begin()->vparam.begin(); param_it < vodeloss.begin()->vparam.end(); param_it++) {
    // unsigned int nn = param_it - vodeloss.begin()->vparam.begin();
    // Rcpp::Rcout << "n_param = " << n_param << ", fixed_param n. " << nn << std::endl;
    fixed_param(param_it - vodeloss.begin()->vparam.begin()) = param_it->fixed;
  }
  // Rcpp::Rcout << "fixed_param " << fixed_param.t() << std::endl;


  /*
  * Algorithm
  */
  if (trace) Rcpp::Rcout << "Commencing backward sweep." << std::endl;
  for (unsigned int i_lambda = n_lambda; i_lambda--> 0; ) {
    /* Resets backtracking parameters and
    * current optimal loss (thus after set_lambda),
    * then moves the particles ahead
    */
    if (trace) Rcpp::Rcout << "Setting lambda number: " << i_lambda + 1 << " (" << lambda(i_lambda) << ")" << std::endl;
    set_lambda(voptim, lambda(i_lambda));

    if (trace) Rcpp::Rcout << "Optimising." << std::endl;
    move(voptim);


    /*
    * Retrieve results
    */
    if (trace) Rcpp::Rcout << "Retrieving results." << std::endl;
    for (loss_it = vodeloss.begin(); loss_it != vodeloss.end(); loss_it++) {

      /*
      * Loop over its parameters
      */
      pars_it = pars.begin();
      for (param_it = loss_it->vparam.begin(); param_it < loss_it->vparam.end(); param_it++) {
        /*
        * Fill in estimates
        */
        pars_it->subcube(0, i_lambda, loss_it - vodeloss.begin(),
          pars_it->n_rows - 1, i_lambda, loss_it - vodeloss.begin()) = param_it->position_full;

        /*
        * Degrees of freedom
        */
        // Rcpp::Rcout << "dfs --> is it fixed? " << fixed_param(param_it - loss_it->vparam.begin()) << std::endl;
        if (fixed_param(param_it - loss_it->vparam.begin()) == 1) {
          dfs(param_it - loss_it->vparam.begin(), i_lambda, loss_it - vodeloss.begin()) = 0;
        } else {
          // Rcpp::Rcout << "Size of get_X4dfs. n_cols " << loss_it->get_X4dfs(param_it).n_cols << std::endl;
          dfs(param_it - loss_it->vparam.begin(), i_lambda, loss_it - vodeloss.begin()) = param_it->dfs(loss_it->get_X4dfs(param_it));
        }

        /*
        * Convergence codes
        */
        if (loss_it->psolver->code_conv > 0) {
          codes(param_it - loss_it->vparam.begin(), i_lambda, loss_it - vodeloss.begin()) = 4;
        } else if (param_it->c.diff_code > 0) {
          codes(param_it - loss_it->vparam.begin(), i_lambda, loss_it - vodeloss.begin()) = 3;
        } else if (param_it->c.step >= param_it->c.step_max) {
          codes(param_it - loss_it->vparam.begin(), i_lambda, loss_it - vodeloss.begin()) = 2;
        } else if (param_it->c.backtrack_code > 0) {
          codes(param_it - loss_it->vparam.begin(), i_lambda, loss_it - vodeloss.begin()) = 1;
        } else {
          codes(param_it - loss_it->vparam.begin(), i_lambda, loss_it - vodeloss.begin()) = 0;
        }

        /*
        * Other diagnostics
        */
        steps(param_it - loss_it->vparam.begin(), i_lambda, loss_it - vodeloss.begin()) = param_it->c.step;
        jerr(param_it - loss_it->vparam.begin(), i_lambda, loss_it - vodeloss.begin())  = param_it->c.backtrack_code + 10 * (param_it->c.step >= param_it->c.step_max) + 100 * param_it->c.diff_code + 1000 * loss_it->psolver->code_conv;
        pens(param_it - loss_it->vparam.begin(), i_lambda, loss_it - vodeloss.begin())  = param_it->preg->penalty(param_it->position_full);

        pars_it++;
      }

      loss_it->evaluate();
      losses(loss_it - vodeloss.begin(), i_lambda) = loss_it->value;
    }
  }

  /*
   * Fill in return list
   */
  if (trace) Rcpp::Rcout << "Finished backward sweep. Preparing results." << std::endl;
  Rcpp::List particle_list(N);
  Rcpp::List param_list(pars.size());
  for (unsigned int i = 0; i < N; ++i) {
    // Rcpp::Rcout << "pars.size() " << pars.size() << std::endl;
    for (pars_it = pars.begin(); pars_it < pars.end(); pars_it++) {
      // Rcpp::Rcout << "pars_it - pars.begin " << pars_it - pars.begin() << std::endl;
      param_list[pars_it - pars.begin()] = arma::conv_to<arma::sp_mat>::from(pars_it->slice(i));
    }
    // Rcpp::Rcout << "i = " << i << ", N = " <<  N << std::endl;
    particle_list[i] = param_list;
  }
  ret["params"]     = particle_list;
  ret["dfs"]        = dfs;
  ret["codes"]      = codes;
  ret["steps"]      = steps;
  ret["losses"]     = losses;
  ret["penalties"]  = pens;
  ret["jerr"]       = jerr;
  ret["lambda"]     = lambda;
  ret["lambda_fac"] = lamfac;

  if (adjust_scales) {
    Rcpp::List scales_list(n_param - 1);
    std::vector<arma::vec>::iterator sc_it = sc_mults.begin();
    for (unsigned int i = 0; i < n_param - 1; ++i) {
      scales_list[i] = *sc_it;
      sc_it++;
    }
    ret["new_scales"] = scales_list;
  }

  return ret;
}
