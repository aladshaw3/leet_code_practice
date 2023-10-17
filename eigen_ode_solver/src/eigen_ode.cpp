/*!
 *  \file eigen_ode.h eigen_ode.cpp
 *  \brief C++ Implementation of ODE solver using Eigen objects
 *  \details This file creates a set of objects to interface with Eigen solvers
 *          for creating and solving ODE systems.
 *
 *  \author Austin Ladshaw
 *  \date 12/08/2022
 *  \copyright This software was designed and built by Austin Ladshaw
 *             for ODE simulations. Copyright (c) 2022, all
 *             rights reserved.
 *
 */

#include "eigen_ode.h"
#include <time.h>
#include <cmath>

// -------------------------- EigenODE Object ---------------------------------

//Default constructor
EigenODE::EigenODE()
{
  this->functor.register_class(this);
  this->time.initializePort(0.);
  this->dt.initializePort(0.);
  this->dt_new = 0.1;
  this->dt_init = 0.1;
  this->end_time = 0.1;
  this->dt_rate_factor = 0.5;
  this->dt_min_converged = 1e-10;
  this->is_initalized = false;
  this->dt_min = 1.0e-14;
  this->dt_max = 1.0;
  this->max_eval = 100;
  this->u_tol = 1e-8;
  this->f_tol = 1e-8;
  this->df_tol = 0.;
  this->error_eps = 0.;
  this->step_bound = 100.;
  this->user_scaling = false;
  this->type = IMPLICIT;
  this->sub_type = BACKWARDS_EULER;
  this->status = NOT_STARTED;
  this->s_type = "IMPLICIT";
  this->s_sub_type = "BACKWARDS_EULER";
  this->user_dt = default_dt;
  this->derivative = default_derivative;
  this->rate = backwards_euler;
  this->coeff = future_coeff;
  this->init_console_out = default_outputs;
  this->post_console_out = default_outputs;
  this->init_file_out = default_outputs;
  this->post_file_out = default_outputs;
  this->out_file = nullptr;
  this->folder = "output";
  this->file = "EigenODE_Result";
}

//Default Destructor
EigenODE::~EigenODE()
{
  if (this->out_file != nullptr)
    fclose(this->out_file);
}

//Return functor by reference
DenseFunctorODE& EigenODE::get_functor()
{
  return this->functor;
}

//Return pointer to user data
const void* EigenODE::get_user_data()
{
	return this->user_data;
}

// Return reference to the new state
Eigen::VectorXd& EigenODE::get_u_new()
{
  return this->u_new;
}

// Return reference to the old state
Eigen::VectorXd& EigenODE::get_u_old()
{
  return this->u.getStateReference();
}

// Return copy of the old state
Eigen::VectorXd EigenODE::get_u_old_copy() const
{
  return this->u.getState();
}

// Return reference to the older ith state
Eigen::VectorXd& EigenODE::get_u_old_state(int i)
{
  return this->u.getStateReference(i);
}

// Return copy of the older ith state
Eigen::VectorXd EigenODE::get_u_old_state_copy(int i) const
{
  return this->u.getState(i);
}

// Return reference to the difference state
Eigen::VectorXd& EigenODE::get_du()
{
  return this->du;
}

// Return reference to the combined rates
Eigen::VectorXd& EigenODE::get_g()
{
  return this->g;
}

// Return reference to the coefficients
Eigen::VectorXd& EigenODE::get_co()
{
  return this->co;
}

// Return the old u value for the ith variable
double EigenODE::get_last_u(unsigned int i) const
{
  return this->u_new[i];
}

// Return the old rate evaulation for the ith variable
double EigenODE::get_last_rate(unsigned int i) const
{
  return this->g[i];
}

// Return the new state time
double EigenODE::get_time_new() const
{
  return this->time_new;
}

// Return the old state time
double EigenODE::get_time_old() const
{
  return this->time.getState();
}

// Return the older state time
double EigenODE::get_time_old_state(int i) const
{
  return this->time.getState(i);
}

// Return the new step
double EigenODE::get_dt_new() const
{
  return this->dt_new;
}

// Return the old step
double EigenODE::get_dt_old() const
{
  return this->dt.getState();
}

// Return the older step
double EigenODE::get_dt_old_state(int i) const
{
  return this->dt.getState(i);
}

// Return the intitial step
double EigenODE::get_dt_initial() const
{
  return this->dt_init;
}

// Return the max step
double EigenODE::get_dt_max() const
{
  return this->dt_max;
}

// Return the min step
double EigenODE::get_dt_min() const
{
  return this->dt_min;
}

// Return the rate factor
double EigenODE::get_dt_factor() const
{
  return this->dt_rate_factor;
}

// Return the upper limit of the rate functions for dt adjustment
double EigenODE::get_rate_upper_limit() const
{
  return this->slope_upper_limit;
}

// Return the lower limit of the rate functions for dt adjustment
double EigenODE::get_rate_lower_limit() const
{
  return this->slope_lower_limit;
}

// Return the lower limit of the dt adjustment
double EigenODE::get_dt_min_converged() const
{
  return this->dt_min_converged;
}

/// Function to return 'true' if the ith function is supposed to be steady-state
bool EigenODE::is_rate_function_steady_state(unsigned int i) const
{
  return this->is_steady_state[i];
}

// Return the ith time derivative
double EigenODE::dudt(int i) const
{
  return this->du[i]/this->dt_new;
}

// Get sub-type name
std::string EigenODE::integration_method() const
{
  return this->s_sub_type;
}

// Get type name
std::string EigenODE::integration_type() const
{
  return this->s_type;
}

// lm status
int EigenODE::get_last_lm_status() const
{
  return this->lm_status;
}

// iterations
int EigenODE::get_last_iter() const
{
  return this->iter;
}

// feval
int EigenODE::get_last_feval() const
{
  return this->f_eval;
}

// jeval
int EigenODE::get_last_jeval() const
{
  return this->df_eval;
}

// fnorm
double EigenODE::get_last_fnorm() const
{
  return this->fnorm;
}

// gnorm
double EigenODE::get_last_gnorm() const
{
  return this->gnorm;
}

// lm_param
double EigenODE::get_last_lm_param() const
{
  return this->lm_param;
}

//Return pointer for output file
FILE * EigenODE::get_file_ptr() const
{
  return this->out_file;
}

//Assign the file ptr
void EigenODE::assign_file_ptr(FILE *ptr)
{
  this->out_file = ptr;
}

//Return string of folder name
std::string EigenODE::output_folder() const
{
  return this->folder;
}

//Return string of file name
std::string EigenODE::output_file() const
{
  return this->file;
}

//Set the size of the system
void EigenODE::set_size(unsigned int size)
{
  this->max_eval = size*100;
  this->functor.set_size(size);
  this->u_new.resize(size,1);
  this->du.resize(size,1);
  this->g.resize(size,1);
  this->co.resize(size,1);
  this->scale.resize(size,1);
  this->user_func.resize(size);
  this->user_coeff.resize(size);
  this->is_steady_state.resize(size);
  this->user_scale.resize(size);
  this->user_ics.resize(size);
  for (unsigned int i=0; i<this->u.getStateSize(); i++)
  {
    this->u.getStateReference(i).resize(size,1);
    this->user_coeff[i] = default_derivative_coeff;
    this->is_steady_state[i] = false;
    this->user_scale[i] = default_derivative_coeff;
    this->user_ics[i] = default_variable_ic;
  }
}

// Set the initial time value
void EigenODE::set_initial_time(double t0)
{
  this->time.setState(t0);
  this->time_new = t0;
}

//Set the current step size
void EigenODE::set_dt(double dt)
{
  if (dt < this->dt_min)
    dt = this->dt_min*0.5;
  if (dt > this->dt_max)
    dt = this->dt_max;
  this->dt_new = dt;
  this->dt_init = dt;
  this->dt_min_converged = dt/1000.0;
}

//Set the current step size (internal version)
void EigenODE::set_dt_internal(double dt)
{
  if (dt < this->dt_min)
    dt = this->dt_min*0.5;
  if (dt > this->dt_max)
    dt = this->dt_max;
  this->dt_new = dt;
}

//Set the end time
void EigenODE::set_end_time(double t)
{
  if (t < 0.)
    t = this->dt_min;
  this->end_time = t;
}

//Set the minimum step size
void EigenODE::set_dt_min(double dt)
{
  if (dt <= DBL_EPSILON*10.0)
    dt = DBL_EPSILON*10.0;
  this->dt_min = dt;
}

//Set the minimum step size
void EigenODE::set_dt_max(double dt)
{
  if (dt <= DBL_EPSILON*100.0)
    dt = DBL_EPSILON*100.0;
	this->dt_max = dt;
}

//Set the maximum function evaluations
void EigenODE::set_max_eval(unsigned int val)
{
  if (val <= this->get_size())
    val = this->get_size()+10;
	this->max_eval = val;
}

//Set the minimum change in u states
void EigenODE::set_u_tol(double val)
{
  if (val <= DBL_EPSILON*10.0)
    val = DBL_EPSILON*10.0;
	this->u_tol = val;
}

//Set the minimum change in residuals
void EigenODE::set_f_tol(double val)
{
  if (val <= DBL_EPSILON*10.0)
    val = DBL_EPSILON*10.0;
  if (val >= 1e-4)
    val = 1e-4;
	this->f_tol = val;
}

//Set the minimum change in derivatives
void EigenODE::set_df_tol(double val)
{
  if (val <= 0.)
    val = 0.;
  if (val >= DBL_EPSILON*10.0)
    val = DBL_EPSILON*10.0;
	this->df_tol = val;
}

//Set the LevenbergMarquardt error precision
void EigenODE::set_error_precision(double val)
{
  if (val <= 0.)
    val = 0.;
  if (val >= DBL_EPSILON*10.0)
    val = DBL_EPSILON*10.0;
	this->error_eps = val;
}

//Set the user scaling option
void EigenODE::set_diag_shift_bound(double val)
{
  if (val <= 1.)
    val = 1.;
	this->step_bound = val;
}

//Set the LevenbergMarquardt diagonal shift bound
void EigenODE::set_user_scaling_option(bool opt)
{
	this->user_scaling = opt;
}

// Set the user integration option
void EigenODE::select_integration_method(integrate_subtype choice)
{
  this->sub_type = choice;
  switch (choice)
  {
    case BACKWARDS_EULER:
      this->type = IMPLICIT;
      this->derivative = default_derivative;
      this->rate = backwards_euler;
      this->coeff = future_coeff;
      this->s_type = "IMPLICIT";
      this->s_sub_type = "BACKWARDS_EULER";
      break;

    case FORWARDS_EULER:
      this->type = EXPLICIT;
      this->derivative = default_derivative;
      this->rate = forwards_euler;
      this->coeff = future_coeff;
      this->s_type = "EXPLICIT";
      this->s_sub_type = "FORWARDS_EULER";
      break;

    case TRAPEZOID_RULE:
      this->type = IMPLICIT;
      this->derivative = default_derivative;
      this->rate = trapezoid_rule;
      this->coeff = future_coeff;
      this->s_type = "IMPLICIT";
      this->s_sub_type = "TRAPEZOID_RULE";
      break;

    case BDF2:
      this->type = IMPLICIT;
      this->derivative = bdf2_derivative;
      this->rate = backwards_euler;
      this->coeff = future_coeff;
      this->s_type = "IMPLICIT";
      this->s_sub_type = "BDF2";
      break;

    case RK45:
      this->type = EXPLICIT;
      this->derivative = default_derivative;
      this->rate = runge_kutta_45;
      this->coeff = future_coeff;
      this->s_type = "EXPLICIT";
      this->s_sub_type = "RK45";
      break;

    default:
      this->sub_type = BACKWARDS_EULER;
      this->type = IMPLICIT;
      this->derivative = default_derivative;
      this->rate = backwards_euler;
      this->coeff = future_coeff;
      this->s_type = "IMPLICIT";
      this->s_sub_type = "BACKWARDS_EULER";
      break;
  }

  // Before ending, check to see is user choice is invalid
  if (this->type == EXPLICIT)
  {
    for (unsigned int i=0; i<this->get_size(); i++)
    {
      if (this->is_steady_state[i])
      {
        SIMError(ode_invalid_type_selection);
        this->sub_type = BACKWARDS_EULER;
        this->type = IMPLICIT;
        this->derivative = default_derivative;
        this->rate = backwards_euler;
        this->coeff = future_coeff;
        this->s_type = "IMPLICIT";
        this->s_sub_type = "BACKWARDS_EULER";
        return;
      }
    }
  }

}

/// Set up default rate based dt selection
void EigenODE::use_default_rate_based_dt(bool yes, double factor, double upper_limit, double lower_limit, double dt_min_con)
{
  if (yes)
  {
    this->user_dt = default_rate_based_dt;
    if (factor < 1.)
      factor = 1.;
    if (factor > 10.)
      factor = 10.;
    if (upper_limit > 1e6)
      upper_limit = 1e6;
    if (upper_limit < 2)
      upper_limit = 2;
    if (lower_limit > 1e-6)
      lower_limit = 1e-6;
    if (lower_limit < 0.5)
      lower_limit = 0.5;
    if (dt_min_con < this->dt_min)
      dt_min_con = this->dt_min*100.0;
    if (dt_min_con > this->dt_max)
      dt_min_con = this->dt_max/100.0;
    this->dt_rate_factor = factor;
    this->slope_upper_limit = upper_limit;
    this->slope_lower_limit = lower_limit;
    this->dt_min_converged = dt_min_con;
  }
}

/// Set up standard console output
void EigenODE::use_standard_console_output(bool yes)
{
  if (yes)
  {
    this->init_console_out = init_console_table_out;
    this->post_console_out = post_console_table_out;
  }
}

/// Set up standard file output
void EigenODE::use_standard_file_output(bool yes, std::string folder, std::string file)
{
  if (yes)
  {
    this->init_file_out = init_file_csv_out;
    this->post_file_out = post_file_csv_out;
    this->folder = folder;
    this->file = file;
  }
}

//Get the size of the system
unsigned int EigenODE::get_size() const
{
  return this->user_func.size();
}

//Register user data
void EigenODE::register_user_data(const void *data)
{
  this->user_data = data;
}

//Register user rate function
void EigenODE::register_rate_function(int i, double (*func) (int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj), bool is_steady )
{
  this->user_func[i] = func;
  if (is_steady)
  {
    this->is_steady_state[i] = true;
    this->user_coeff[i] = default_variable_ic;
  }
}

//Register user ceoff function
void EigenODE::register_coeff_function(int i, double (*func) (int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj) )
{
  this->user_coeff[i] = func;
}

//Register user scaling function
void EigenODE::register_scaling_function(int i, double (*func) (int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj) )
{
  this->user_scale[i] = func;
}

//Register user initial condition function
void EigenODE::register_initial_condition_function(int i, double (*func) (int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj) )
{
  this->user_ics[i] = func;
}

//Register user initial condition function
void EigenODE::register_time_stepper_function(double (*func) (const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj) )
{
  this->user_dt = func;
}

//Eval user rate function i
double EigenODE::eval_user_rate_function(int i, const Eigen::VectorXd &u, double t) const
{
	return this->user_func[i](i, u, t, this->user_data, *this);
}

//Eval user coeff function i
double EigenODE::eval_user_coeff_function(int i, const Eigen::VectorXd &u, double t) const
{
	return this->user_coeff[i](i, u, t, this->user_data, *this);
}

//Eval user scaling function i
double EigenODE::eval_user_scaling_function(int i, const Eigen::VectorXd &u, double t) const
{
	return this->user_scale[i](i, u, t, this->user_data, *this);
}

//Eval user initial condition function i
double EigenODE::eval_user_initial_condition_function(int i, const Eigen::VectorXd &u, double t) const
{
	return this->user_ics[i](i, u, t, this->user_data, *this);
}

//Evaluate current derivative
int EigenODE::eval_derivative(const Eigen::VectorXd &u)
{
	return this->derivative(u, *this);
}

//Evaluate current rate
int EigenODE::eval_rate(const Eigen::VectorXd &u)
{
	return this->rate(u, *this);
}

//Evaluate coefficients
int EigenODE::eval_coeff(const Eigen::VectorXd &u)
{
	return this->coeff(u, *this);
}

//Evaluate scaling
int EigenODE::eval_scaling(const Eigen::VectorXd &u)
{
  for (unsigned int i=0; i<this->get_size(); i++)
  {
    this->scale[i] = this->eval_user_scaling_function(i, u, this->get_time_new());
  }
  return 0;
}

// Function to initialize the system
int EigenODE::initialize_system(Eigen::LevenbergMarquardt< Eigen::NumericalDiff< DenseFunctorODE > > &lm)
{
  int success = 0;

  this->set_dt_internal(this->dt_init);
  this->dt_init = this->dt_new;
  this->dt.initializePort(DBL_MAX);
  this->time.getStateReference() = this->time_new;

  for (unsigned int i=0; i<this->get_size(); i++)
  {
    this->u_new[i] = this->eval_user_initial_condition_function(i, this->u_new, this->time_new);
  }
  this->u.initializePort(this->u_new);

  lm.setMaxfev(this->max_eval);  // Maximum calls to f()
  lm.setXtol(this->u_tol);     // Cutoff tolerance for change in x
  lm.setFtol(this->f_tol);     // Cutoff tolerance for norm in f
  lm.setGtol(this->df_tol);       // Cutoff tolerance for change in df()
  lm.setEpsilon(this->error_eps);    // Error percision
  lm.setFactor(this->step_bound);    // Step bound for diagonal shift
  lm.setExternalScaling(this->user_scaling);

  //Prepare any outputs
  success = this->init_console_out(*this);
  success = this->init_file_out(*this);
  if (success == FAILED)
  {
    SIMError(ode_failure);
    return FAILED;
  }
  if (success == USER_REQUEST)
  {
    return USER_REQUEST;
  }

  this->is_initalized = true;
  this->status = RUNNING;

  return RUNNING;
}

// Function to perform pre-step actions (compute dt, modify scaling, etc.)
int EigenODE::preprocess_step(Eigen::LevenbergMarquardt< Eigen::NumericalDiff< DenseFunctorODE > > &lm)
{
  int success = 0;

  this->set_dt_internal(this->user_dt(this->u_new, this->time_new, this->user_data, *this));
  this->time_new = this->dt_new + this->time_new;
  success = this->eval_scaling(this->u_new);
  success = RUNNING;
  this->status = RUNNING;
  lm.diag() = this->scale;

  return success;
}

// Function to evaluate a step (implicit)
int EigenODE::solve_step_implicit(Eigen::LevenbergMarquardt< Eigen::NumericalDiff< DenseFunctorODE > > &lm)
{
  int info = lm.minimize(this->u_new);
  this->lm_status = info;
  switch (info)
  {
    case -2:
      this->status = FAILED;
      break;

    case -1:
      this->status = FAILED;
      break;

    case 0:
      this->status = FAILED;
      break;

    case 5:
      this->status = FAILED;
      break;

    case 9:
      this->status = USER_REQUEST;
      break;

    default:
      this->status = RUNNING;
      break;
  }
  this->iter = lm.iterations();
  this->f_eval = lm.nfev();
  this->df_eval = lm.njev();
  this->fnorm = lm.fnorm();
  this->gnorm = lm.gnorm();
  this->lm_param = lm.lm_param();
  return this->status;
}

// Function to evaluate a step (explicitly)
int EigenODE::solve_step_explicit(Eigen::LevenbergMarquardt< Eigen::NumericalDiff< DenseFunctorODE > > &lm)
{
  //int info = lm.minimize(this->u_new);
  int success = 0;

  unsigned int N = (unsigned int)(50.0*this->dt_new/this->dt_max);
  double dt_hold = this->dt_new;
  double sub_dt = this->dt_new/(double)(N+1);
  double time_hold = this->time_new;
  this->time_new = this->time.getState();
  this->dt_new = sub_dt;

  // Perform a maximum of N+1 evenly spaced inner steps
  for (unsigned int n=0; n<N+1; n++)
  {
    this->time_new = this->time_new + sub_dt;
    success = this->eval_derivative(this->u_new);
    success = this->eval_rate(this->u_new);
    success = this->eval_coeff(this->u_new);
    for (unsigned int i=0; i<this->get_size(); i++)
    {
      this->u_new[i] = this->u_new[i] + (this->g[i] * sub_dt)/this->co[i];
    }
  }
  this->time_new = time_hold;
  this->dt_new = dt_hold;
  this->lm_status = 10;
  this->iter = N+1;
  this->f_eval = 1;
  this->df_eval = 0;
  this->fnorm = 0;
  this->gnorm = 0;
  this->lm_param = 0;
  success = RUNNING;
  this->status = RUNNING;
  return success;
}

// Function to perform post-step actions (update ports, prepare outputs, etc.)
int EigenODE::postprocess_step(Eigen::LevenbergMarquardt< Eigen::NumericalDiff< DenseFunctorODE > > &lm)
{
  int success = 0;

  this->time.setState(this->time_new);
  this->dt.setState(this->dt_new);
  this->u.setState(this->u_new);

  //Prepare any outputs
  success = this->post_console_out(*this);
  success = this->post_file_out(*this);
  success = RUNNING;
  this->status = RUNNING;

  return success;
}

// Function to solve to end of time
int EigenODE::solve_to_end(Eigen::LevenbergMarquardt< Eigen::NumericalDiff< DenseFunctorODE > > &lm)
{
  int success = 0;
  this->status = RUNNING;
  if (this->type == IMPLICIT)
  {
    success = this->solve_to_end_implicit(lm);
  }
  else
  {
    success = this->solve_to_end_explicit(lm);
  }
  this->status = (eigen_ode_status)success;

  return COMPLETE;
}

// Function to solve to end of time (implicit only)
int EigenODE::solve_to_end_implicit(Eigen::LevenbergMarquardt< Eigen::NumericalDiff< DenseFunctorODE > > &lm)
{
  int success = 0;
  this->status = RUNNING;
  // Initial setup (if not already initialized)
  if (!this->is_initalized)
  {
    success = this->initialize_system(lm);
    if (success == FAILED)
    {
      SIMError(ode_failure);
      return FAILED;
    }
    if (success == USER_REQUEST)
    {
      return USER_REQUEST;
    }
  }

  // Initial preprocess step
  success = preprocess_step(lm);
  do
  {
    success = solve_step_implicit(lm);
    if (success == USER_REQUEST)
    {
      std::cout << "WARNING: Simulation stopped via user request...\n\n";
      return USER_REQUEST;
    }
    if (success == FAILED)
    {
      // Upon a solver failure, we want to 'retry' the current
      // failed time step at a smaller dt value.
      this->set_dt_internal(this->dt_new/2.0);

      // If we cannot reduce the step anymore, then it is a true failure
      if (this->dt_new < this->dt_min)
      {
        this->status = FAILED;
        SIMError(ode_dt_min_fail);
        return this->status;
      }

      this->time_new = this->dt_new + this->time.getState();
      success = this->eval_scaling(this->u_new);
      lm.diag() = this->scale;
      this->status = RUNNING;
    }
    else
    {
      success = postprocess_step(lm);
      success = preprocess_step(lm);
    }
  } while(this->time_new < (this->end_time-this->dt_min));

  // Final step (to exactly hit 'end_time')
  if (this->time_new > this->end_time)
  {
    this->dt_new = this->end_time - (this->time_new - this->dt_new);
    this->time_new = this->end_time;
  }
  success = solve_step_implicit(lm);
  success = postprocess_step(lm);
  this->status = COMPLETE;

  return COMPLETE;
}

// Function to solve to end of time (explicit only)
int EigenODE::solve_to_end_explicit(Eigen::LevenbergMarquardt< Eigen::NumericalDiff< DenseFunctorODE > > &lm)
{
  int success = 0;
  this->status = RUNNING;
  // Initial setup (if not already initialized)
  if (!this->is_initalized)
  {
    success = this->initialize_system(lm);
    if (success == FAILED)
    {
      SIMError(ode_failure);
      return FAILED;
    }
    if (success == USER_REQUEST)
    {
      return USER_REQUEST;
    }
  }

  // Initial preprocess step
  success = preprocess_step(lm);
  do
  {
    success = solve_step_explicit(lm);
    if (success == USER_REQUEST)
    {
      std::cout << "WARNING: Simulation stopped via user request...\n\n";
      return USER_REQUEST;
    }
    if (success == FAILED)
    {
      // Upon a solver failure, we want to 'retry' the current
      // failed time step at a smaller dt value.
      this->set_dt_internal(this->dt_new/2.0);

      // If we cannot reduce the step anymore, then it is a true failure
      if (this->dt_new < this->dt_min)
      {
        this->status = FAILED;
        SIMError(ode_dt_min_fail);
        return this->status;
      }

      this->time_new = this->dt_new + this->time.getState();
      success = this->eval_scaling(this->u_new);
      lm.diag() = this->scale;
      this->status = RUNNING;
    }
    else
    {
      success = postprocess_step(lm);
      success = preprocess_step(lm);
    }
  } while(this->time_new < (this->end_time-this->dt_min));

  // Final step (to exactly hit 'end_time')
  if (this->time_new > this->end_time)
  {
    this->dt_new = this->end_time - (this->time_new - this->dt_new);
    this->time_new = this->end_time;
  }
  success = solve_step_explicit(lm);
  success = postprocess_step(lm);
  this->status = COMPLETE;

  return COMPLETE;
}

// ------------------------- End EigenODE Object ------------------------------


// -------------------------- EigenODE Helper Functions ---------------------------------

// Default derivative function (computes: du = u_new - u.getStateReference())
int default_derivative(const Eigen::VectorXd &u, EigenODE &obj)
{
  obj.get_du() = u - obj.get_u_old();
  return 0;
}

// BDF2 derivative function (computes: du = a*u_new - b*u + c*u_old)
int bdf2_derivative(const Eigen::VectorXd &u, EigenODE &obj)
{
  double t_ratio = obj.get_dt_new()/obj.get_dt_old();
  double oneplus_t_ratio = 1.+t_ratio;
  double inv_oneplus_t_ratio = 1./oneplus_t_ratio;
  double a = (1.+2.*t_ratio)*inv_oneplus_t_ratio;
  double b = oneplus_t_ratio;
  double c = (t_ratio*t_ratio)*inv_oneplus_t_ratio;
  obj.get_du() = u*a - obj.get_u_old()*b + obj.get_u_old_state(1)*c;
  return 0;
}

// 1st Order Backwards Euler rate function evaluation
int backwards_euler(const Eigen::VectorXd &u, EigenODE &obj)
{
  for (unsigned int i=0; i<obj.get_size(); i++)
  {
    obj.get_g()[i] = obj.eval_user_rate_function(i, u, obj.get_time_new());
  }
  return 0;
}

// 1st Order forwards Euler rate function evaluation
int forwards_euler(const Eigen::VectorXd &u, EigenODE &obj)
{
  for (unsigned int i=0; i<obj.get_size(); i++)
  {
    obj.get_g()[i] = obj.eval_user_rate_function(i, u, obj.get_time_new());
  }
  return 0;
}

// Trapezoid rule rate function evaluation
int trapezoid_rule(const Eigen::VectorXd &u, EigenODE &obj)
{
  for (unsigned int i=0; i<obj.get_size(); i++)
  {
    obj.get_g()[i] = (obj.eval_user_rate_function(i, u, obj.get_time_new()) +
                    obj.eval_user_rate_function(i, obj.get_u_old(), obj.get_time_old()))*0.5;
  }
  return 0;
}

// Runge-Kutta 45 rate function evaluation
int runge_kutta_45(const Eigen::VectorXd &u, EigenODE &obj)
{
  double k1, k2, k3, k4;
  double h_2 = obj.get_dt_new()/2.0;
  double h = obj.get_dt_new();
  for (unsigned int i=0; i<obj.get_size(); i++)
  {
    k1 = obj.eval_user_rate_function(i, u, obj.get_time_new());
    k2 = obj.eval_user_rate_function(i, u.array()+h_2*k1, obj.get_time_new()+h_2);
    k3 = obj.eval_user_rate_function(i, u.array()+h_2*k2, obj.get_time_new()+h_2);
    k4 = obj.eval_user_rate_function(i, u.array()+h*k3, obj.get_time_new()+h);
    obj.get_g()[i] = (k1+2.*k2+2.*k3+k4)/6.0;
  }
  return 0;
}

// Future Coeffs
int future_coeff(const Eigen::VectorXd &u, EigenODE &obj)
{
  for (unsigned int i=0; i<obj.get_size(); i++)
  {
    obj.get_co()[i] = obj.eval_user_coeff_function(i, u, obj.get_time_new());
  }
  return 0;
}

// past Coeffs
int past_coeff(const Eigen::VectorXd &u, EigenODE &obj)
{
  for (unsigned int i=0; i<obj.get_size(); i++)
  {
    obj.get_co()[i] = obj.eval_user_coeff_function(i, obj.get_u_old(), obj.get_time_old());
  }
  return 0;
}

// Default coefficient in front of derivatives
double default_derivative_coeff(int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj)
{
  return 1.0;
}

// Default variable initial condtions
double default_variable_ic(int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj)
{
  return 0.0;
}

// Default function to get next dt
double default_dt(const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj)
{
  return std::min(obj.get_dt_new()*1.5, obj.get_dt_initial());
}

// Simple function to get next dt based on approximate rates
double default_rate_based_dt(const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj)
{
  // Iterate through all rate functions to find the fastest rate
  // Adjust the dt value from the previous dt based on the
  //  fastest object rate in the system
  double dt_max_static = std::min(obj.get_dt_new()*1.5, obj.get_dt_initial());
  double dt_max = -1;
  for (unsigned int i=0; i<obj.get_size(); i++)
  {
    // Only the non-steady-state rate evaluates matter
    if (!obj.is_rate_function_steady_state(i))
    {
      //compute future rate (Absolute value)
      double rate = fabs( obj.eval_user_rate_function(i, u, obj.get_time_new()) / obj.eval_user_coeff_function(i, u, obj.get_time_new()) );

      //compare rate to upper and lower slope limits limits and compute a dt
      double dt_calc = dt_max_static;
      if (rate > obj.get_rate_upper_limit())
      {
        dt_calc = dt_max_static - (1.0/obj.get_dt_factor())*dt_max_static;
      }
      if (rate < obj.get_rate_lower_limit())
      {
        dt_calc = dt_max_static+obj.get_dt_factor()*dt_max_static;
      }
      if (dt_calc < obj.get_dt_min_converged())
      {
        dt_calc = obj.get_dt_min_converged();
      }

      // Comparse calculated dt to dt_max and adjust
      if (dt_calc < dt_max || dt_max < 0)
      {
        dt_max = dt_calc;
      }
    }
  }
  return dt_max;
}


// Default output for console (do nothing)
int default_outputs(EigenODE &obj)
{
  return 0;
}

// Return status message based on int status
std::string lm_status_message(int lm_status)
{
  std::string message;
  switch (lm_status)
  {
    case -2:
      message = "Not Started";
      break;

    case -1:
      message = "Still Running";
      break;

    case 0:
      message = "FAILED: Improper Inputs";
      break;

    case 1:
      message = "CONVERGED: Relative reduction below tolerance";
      break;

    case 2:
      message = "CONVERGED: Relative error below tolerance";
      break;

    case 3:
      message = "CONVERGED: Relative reduction and error below tolerance";
      break;

    case 4:
      message = "CONVERGED: Cosine of angles between jacobian columns have minimized";
      break;

    case 5:
      message = "FAILED: Too many function evaluations. May need to restart...";
      break;

    case 6:
      message = "CONVERGED: Absolute function residual below tolerance";
      break;

    case 7:
      message = "CONVERGED: Absolute difference between 'u' states below tolerance";
      break;

    case 8:
      message = "CONVERGED: Absolute jacobian residual below tolerance";
      break;

    case 9:
      message = "Stopped by user request";
      break;

    default:
      message = "Not Used";
      break;
  }
  return message;
}


// initial stage for console outputs
int init_console_table_out(EigenODE &obj)
{
  std::cout << " ---------------------- Results Table ------------------------ \n";
  std::cout << "\tIntegration Type:\t" << obj.integration_type() << std::endl;
  std::cout << "\tIntegration Method:\t" << obj.integration_method() << std::endl;
  std::cout << std::endl;

  std::cout << "Time";
  for (unsigned int i=0; i<obj.get_size(); i++)
  {
    std::cout << "\tu[" << i << "]";
  }
  std::cout << "\tfnorm\tgnorm\tlm_param\tf_eval\titer\tLM_status\n";

  std::cout << obj.get_time_old();
  for (unsigned int i=0; i<obj.get_size(); i++)
  {
    std::cout << "\t" << obj.get_u_new()[i];
  }
  std::cout << "\t" << 0. << "\t" << 0. << "\t" << 0. << "\t" << 0 << "\t" << 0;
  std::cout << "\t" << lm_status_message(-2);

  std::cout << std::endl;

  return 0;
}

// post-process stage for console outputs
int post_console_table_out(EigenODE &obj)
{
  std::cout << obj.get_time_new();
  for (unsigned int i=0; i<obj.get_size(); i++)
  {
    std::cout << "\t" << obj.get_u_new()[i];
  }
  std::cout << "\t" << obj.get_last_fnorm()
            << "\t" << obj.get_last_gnorm()
            << "\t" << obj.get_last_lm_param()
            << "\t" << obj.get_last_feval()
            << "\t" << obj.get_last_iter();
  std::cout << "\t" << lm_status_message(obj.get_last_lm_status());

  std::cout << std::endl;
  return 0;
}

// Standard file output (initial)
int init_file_csv_out(EigenODE &obj)
{
  int success = 0;
  // If a file pointer is not provided, then create one
  if (obj.get_file_ptr() == nullptr)
  {
    obj.assign_file_ptr( fopen( (obj.output_folder()+"/"+obj.output_file()+".csv").c_str(), "w+" ) );
  }
  // If file pointer still does not exist, this means the folder didn't exist
  if (obj.get_file_ptr() == nullptr)
  {
    success = system( ("mkdir "+obj.output_folder()).c_str() );
    obj.assign_file_ptr( fopen( (obj.output_folder()+"/"+obj.output_file()+".csv").c_str(), "w+" ) );
  }

  // Print out the standard header
  fprintf(obj.get_file_ptr(), "Time");
  for (unsigned int i=0; i<obj.get_size(); i++)
  {
    fprintf(obj.get_file_ptr(), ",u[%i]", i);
  }
  fprintf(obj.get_file_ptr(), ",fnorm,gnorm,lm_param,f_eval,iter,LM_status\n");

  // Print out the initial conditions
  fprintf(obj.get_file_ptr(), "%.6g", obj.get_time_old());
  for (unsigned int i=0; i<obj.get_size(); i++)
  {
    fprintf(obj.get_file_ptr(), ",%.6g", obj.get_u_new()[i]);
  }
  fprintf(obj.get_file_ptr(), ",%.6g,%.6g,%.6g,%i,%i,%s\n",
    0., 0., 0., 0, 0, lm_status_message(-2).c_str());
  return success;
}

// Standard file output (post)
int post_file_csv_out(EigenODE &obj)
{
  // Print out the results
  fprintf(obj.get_file_ptr(), "%.6g", obj.get_time_new());
  for (unsigned int i=0; i<obj.get_size(); i++)
  {
    fprintf(obj.get_file_ptr(), ",%.6g", obj.get_u_new()[i]);
  }
  fprintf(obj.get_file_ptr(), ",%.6g,%.6g,%.6g,%i,%i,%s\n",
    obj.get_last_fnorm(), obj.get_last_gnorm(), obj.get_last_lm_param(),
      obj.get_last_feval(), obj.get_last_iter(), lm_status_message(obj.get_last_lm_status()).c_str());
  return 0;
}

// ------------------------ END EigenODE Helper Functions -------------------------------


// ------------------- Test 1 set ---------------------------------------------
/**
            u0 --> u1
      R*d(u0)/dt = -k*u0
      d(u1)/dt = k*u0
*/
struct test1_data
{
  double k;
  double R;
};

double test1_rate1(int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj)
{
  // User function must properly typecast their own data structure (if applicable)
  test1_data *dat = (test1_data *) user_data;

  return -dat->k*u[0];
}

double test1_rate2(int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj)
{
  // User function must properly typecast their own data structure (if applicable)
  test1_data *dat = (test1_data *) user_data;

  return dat->k*u[0];
}

double test1_coeff1(int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj)
{
  // User function must properly typecast their own data structure (if applicable)
  test1_data *dat = (test1_data *) user_data;

  return dat->R;
}

double test1_ic1(int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj)
{
  return 1;
}

double test1_ic2(int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj)
{
  return 0;
}

int test1_result()
{
  int success = 0;

  EigenODE test;
  test1_data data;
  data.k = 1.0;
  data.R = 2.0;

  // Step 1 - setup size of problem and pass any data needed
  test.set_size(2);
  test.register_user_data((void*)&data);

  // Step 2 - register user functions and setup initial conditions
  test.register_rate_function(0,test1_rate1);
  test.register_rate_function(1,test1_rate2);
  test.register_coeff_function(0,test1_coeff1);
  test.register_initial_condition_function(0,test1_ic1);
  test.register_initial_condition_function(1,test1_ic2);
  test.set_initial_time(0);
  test.set_dt(0.1);
  test.set_end_time(1.0);

  // Step 3 - selection of integration method and time-steppers
  test.select_integration_method(BDF2);
  test.use_standard_console_output(true);
  test.use_standard_file_output(true);
  test.use_default_rate_based_dt(false, 2.0, 10.0, 0.1, 0.01);

  // Step 4 - Create Eigen objects to pass to solver
  Eigen::NumericalDiff< DenseFunctorODE > numDiff( test.get_functor() );
  Eigen::LevenbergMarquardt< Eigen::NumericalDiff< DenseFunctorODE > > lm(numDiff);

  // Step 5 - Call the solver
  success = test.solve_to_end(lm);

  // Test a 'Restart' of the system
  test.set_end_time(2.0);
  success = test.solve_to_end(lm);


  // Test a 2nd 'Restart' of the system
  test.set_end_time(3.25);
  test.set_dt(0.34);
  success = test.solve_to_end(lm);

  return success;
}
// -----------------------END Test 1 ------------------------------------------



// ------------------- Test 2 set ---------------------------------------------
/**
      Synchronous Buckboost
      ---------------------
      L*d(iL)/dt = vin - rL*iL - (1-d)*vc
      C*vc*d(vc)/dt = po + (1-d)*iL*vc

        where po = Vac*Iac*(1-cos(2*w*t))
              w = 2*pi*f0
*/
struct test2_data
{
  double L;
  double vin;
  double rL;
  double d;
  double C;
  double Vac;
  double Iac;
  double f0;
};

// Var iL == u[0]
double test2_rate_for_iL(int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj)
{
  // User function must properly typecast their own data structure (if applicable)
  test2_data *dat = (test2_data *) user_data;

  return dat->vin - dat->rL*u[0] - (1.-dat->d)*u[1];
}

// Var vc == u[1]
double test2_rate_for_vc(int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj)
{
  // User function must properly typecast their own data structure (if applicable)
  test2_data *dat = (test2_data *) user_data;

  // calc power
  double w = 2.*3.14159*dat->f0;
  double po = dat->Vac*dat->Iac*(1.-cos(2.*w*t));

  return po + (1.-dat->d)*u[0]*u[1];
}

//Coefficient for iL
double test2_coeff_iL(int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj)
{
  // User function must properly typecast their own data structure (if applicable)
  test2_data *dat = (test2_data *) user_data;

  return dat->L;
}

//Coefficient for vc
double test2_coeff_vc(int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj)
{
  // User function must properly typecast their own data structure (if applicable)
  test2_data *dat = (test2_data *) user_data;

  return dat->C*u[1];
}

// Initial condition for iL
double test2_ic_iL(int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj)
{
  return 0.;
}

// Initial condition for vc
double test2_ic_vc(int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj)
{
  return 400.;
}

// Run test simulation
int test2_result()
{
  int success = 0;

  EigenODE test;
  test2_data data;

  // Setup params
  data.f0 = 50.0;
  data.L = 0.0006;
  data.C = 0.0007;
  data.vin = 100.;
  data.rL = 0.01;
  data.d = 0.7;
  data.Iac = 10.;
  data.Vac = 200.;

  // Step 1 - setup size of problem and pass any data needed
  test.set_size(2);
  test.register_user_data((void*)&data);

  // Step 2 - register user functions and setup initial conditions
  test.register_rate_function(0,test2_rate_for_iL); // iL
  test.register_rate_function(1,test2_rate_for_vc); // vc

  test.register_coeff_function(0,test2_coeff_iL);
  test.register_coeff_function(1,test2_coeff_vc);

  test.register_initial_condition_function(0,test2_ic_iL);
  test.register_initial_condition_function(1,test2_ic_vc);

  test.set_initial_time(0);

  // Calc and set dt_max and dt
  double dt_max = 0.1/2./3.14159/data.f0;
  double dt = 0.1/1000.;
  test.set_dt_max(dt_max);
  test.set_dt(dt);
  test.set_end_time(2.5);

  // Step 3 - selection of integration method and time-steppers
  test.select_integration_method(FORWARDS_EULER);
  test.use_standard_file_output(true, "output", "EigenODE_SyncBuckBoost");
  test.use_default_rate_based_dt(false);

  // Step 4 - Create Eigen objects to pass to solver
  Eigen::NumericalDiff< DenseFunctorODE > numDiff( test.get_functor() );
  Eigen::LevenbergMarquardt< Eigen::NumericalDiff< DenseFunctorODE > > lm(numDiff);

  double clock_time = clock();
  // Step 5 - Call the solver
  success = test.solve_to_end(lm);

  // Update conditions and restart
  data.d = 0.6;
  test.set_end_time(5.0);
  success = test.solve_to_end(lm);

  // Update conditions and restart
  data.d = 0.7;
  data.Iac = -1.;
  test.set_end_time(7.5);
  success = test.solve_to_end(lm);

  // Update conditions and restart
  data.d = 0.6;
  test.set_end_time(10.0);
  success = test.solve_to_end(lm);

  // Update conditions and restart
  data.d = 0.7;
  data.Iac = 10.;
  test.set_end_time(12.5);
  success = test.solve_to_end(lm);

  // Update conditions and restart
  data.d = 0.6;
  test.set_end_time(15.0);
  success = test.solve_to_end(lm);

  clock_time = clock() - clock_time;
  std::cout << "Test2: Simulation time (s):\t" << (clock_time / CLOCKS_PER_SEC) << std::endl;
  return success;
}
// -----------------------END Test 2 ------------------------------------------


int Test_Eigen_ODE()
{
 int success = 0;

 success = test1_result();
 if (success != 0)
 {
   SIMError(test_failed);
   return -1;
 }
 std::cout << std::endl;

 success = test2_result();
 if (success != 0)
 {
   SIMError(test_failed);
   return -1;
 }
 std::cout << std::endl;

 return success;
}
