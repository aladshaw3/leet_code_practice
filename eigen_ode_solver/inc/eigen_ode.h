/*!
 *  \file eigen_ode.h eigen_ode.cpp
 *  \brief C++ Implementation of ODE solver using Eigen objects
 *  \details This file creates a set of objects to interface with Eigen solvers
 *          for creating and solving ODE systems.
 *
 *  \author Austin Ladshaw
 *  \date 12/08/2022
 *  \copyright This software was designed and built  by Austin Ladshaw
 *             for ODE simulations. Copyright (c) 2022, all
 *             rights reserved.
 *
 */

#pragma once
#ifndef _EIGEN_ODE_
#define _EIGEN_ODE_

#include <iostream>
#include <stdio.h>
#include <vector>
#include <cfloat>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/LevenbergMarquardt>

#include "port_template.h"
#include "sim_error.h"

/// Enumeration for the list of valid time integration types
/** The only types that have been defined are for Implicit and Explicit methods.
	Sub-type enumeration is used to denote the specific methods.*/
typedef enum {IMPLICIT, EXPLICIT} integrate_type;

/// Enumeration for the list of valid time integration subtypes
/** Theses subtypes define the specific scheme to be used. The table below gives
	a brief description of each.

	\param BACKWARDS_EULER Standard implicit method.
	\param FORWARDS_EULER Standard explicit method.
	\param TRAPEZOID_RULE Time averaged, 2nd order implicit scheme.
	\param BDF2 Backwards-Differentiation-Formula-2: 2nd order implicit method.
	\param RK45 Runge-Kutta-45: 4th order explicit method. */
typedef enum {BACKWARDS_EULER, FORWARDS_EULER, TRAPEZOID_RULE,
              BDF2, RK45} integrate_subtype;


/// Enumeration for the list of EigenODE status states
/** Theses enums will denote the current running state of the EigenODE object

	\param NOT_STARTED Solver has not started running.
	\param RUNNING Solver is currently running.
	\param COMPLETE Solver has finished at requested end_time.
	\param FAILED Solver has failed to converge.
  \param USER_REQUEST User has requested to stop in one of their functions*/
typedef enum {NOT_STARTED = -2,
              RUNNING = -1,
              COMPLETE = 0,
              FAILED = 1,
              USER_REQUEST = 2} eigen_ode_status;

/// Forward struct declaration
struct DenseFunctorODE;

/// EigenODE Base class (needed for reference in functor)
class EigenODE_Base
{
public:
  ///< Function to evaluate the current derivative
  virtual int eval_derivative(const Eigen::VectorXd &u) {return 0;};

  ///< Function to evaluate the current rate
  virtual int eval_rate(const Eigen::VectorXd &u) {return 0;};

  ///< Function to evaluate the coefficients
  virtual int eval_coeff(const Eigen::VectorXd &u) {return 0;};

  ///< Return reference to the difference state
  virtual Eigen::VectorXd &get_du() {return this->dummy;};

  ///< Return reference to the combined rates
  virtual Eigen::VectorXd &get_g() {return this->dummy;};

  ///< Return reference to the coefficients
  virtual Eigen::VectorXd &get_co() {return this->dummy;};

  ///< Return the new state time
  virtual double get_dt_new() const {return 0;}

private:
  Eigen::VectorXd dummy;
};


/// Data structure to interact with Eigen
/** This data structure is the primary structure that interacts
    with the Eigen library LevenbergMarquardt algorithm to solve
    a generalized non-linear system. We utilize that solver to
    step through a system of ODEs implicitly.
*/
struct DenseFunctorODE : Eigen::DenseFunctor<double>
{
    /* Data and Parameters */

    unsigned int vars;    /// Number of variables
    unsigned int funs;    /// Number of functions
    EigenODE_Base *ode_class;  /// Pointer to ODE class object


    /* Methods and Constructors */

    ///< Default constructor
    DenseFunctorODE(void) : Eigen::DenseFunctor<double>() {}

    ///< Assign the class pointer
    void register_class(EigenODE_Base *p) {ode_class=p;}

    ///< Return pointer to class
    EigenODE_Base * get_class() {return ode_class;}

    ///< Set the number of variables and functions
    void set_size(unsigned int size) {vars=size; funs=size;}

    ///< Return the problem size
    unsigned int get_size() {return vars;};


    /* REQUIRED Overrides for LevenbergMarquardt */

    ///< Override the inputs() function
    int inputs() const { return vars; }

    ///< Override the values() function
    int values() const { return funs; }

    ///< Non-linear residual function
    int operator()(const Eigen::VectorXd &u, Eigen::VectorXd &fvec) const
    {
      int success = 0;
      success = ode_class->eval_derivative(u);
      success = ode_class->eval_rate(u);
      success = ode_class->eval_coeff(u);
      fvec = (ode_class->get_du().array()*ode_class->get_co().array());
      fvec = fvec - ode_class->get_g() * ode_class->get_dt_new();
      return success;
    }
};


/// EigenODE Class
/** This class structure creates a C++ object that can be used to solve
    a coupled system of ODEs.
*/
class EigenODE : EigenODE_Base
{
public:
  EigenODE();												///< Default constructor
  ~EigenODE();											///< Default destructor

  DenseFunctorODE &get_functor();   ///< Function to return functor by reference

  const void *get_user_data();      ///< Function to return pointer to user data

  Eigen::VectorXd &get_u_new();             ///< Return reference to the new state
  Eigen::VectorXd &get_u_old();             ///< Return reference to the old state
  Eigen::VectorXd get_u_old_copy() const;   ///< Return copy of the old state
  Eigen::VectorXd &get_u_old_state(int i);  ///< Return reference to the ith old state
  Eigen::VectorXd get_u_old_state_copy(int i) const;  ///< Return copy of the ith old state
  Eigen::VectorXd &get_du();                ///< Return reference to the difference state
  Eigen::VectorXd &get_g();                 ///< Return reference to the combined rates
  Eigen::VectorXd &get_co();                ///< Return reference to the coefficients

  double get_last_u(unsigned int i) const;     ///< Return the last evaluated value of state u_i
  double get_last_rate(unsigned int i) const;  ///< Return the last evaluted value of rate g_i

  double get_time_new() const;             ///< Return the new state time
  double get_time_old() const;             ///< Return old state time
  double get_time_old_state(int i) const;  ///< Return ith old state time

  double get_dt_new() const;             ///< Return the new state time
  double get_dt_old() const;             ///< Return old state time
  double get_dt_old_state(int i) const;  ///< Return ith old state time
  double get_dt_initial() const;         ///< Return the initial time step
  double get_dt_max() const;             ///< Return the max time step
  double get_dt_min() const;             ///< Return the min time step

  double get_dt_factor() const;          ///< Return the dt rate change factor
  double get_rate_upper_limit() const;   ///< Return the upper limit of the rate functions for dt adjustment
  double get_rate_lower_limit() const;   ///< Return the lower limit of the rate functions for dt adjustment
  double get_dt_min_converged() const;   ///< Return the lower limit of dt when solution converges

  /// Function to return 'true' if the ith function is supposed to be steady-state
  bool is_rate_function_steady_state(unsigned int i) const;

  double dudt(int i) const;              ///< Return the ith time derivative

  std::string integration_method() const; ///< Return the string of the integration method

  std::string integration_type() const;   ///< Return the string of the integration type

  int get_last_lm_status() const;         ///< Return last valid lm_status
  int get_last_iter() const;              ///< Return last valid iteration count
  int get_last_feval() const;             ///< Return last valid func evals
  int get_last_jeval() const;             ///< Return last valid jacobian evals
  double get_last_fnorm() const;          ///< Return last valid fnorm
  double get_last_gnorm() const;          ///< Return last valid gnorm
  double get_last_lm_param() const;       ///< Return last valid lm_param

  FILE *get_file_ptr() const;             ///< Return pointer to the output file
  void assign_file_ptr(FILE *ptr);        ///< Assign the file pointer to a given pointer 
  std::string output_folder() const;      ///< Return the string for the output folder
  std::string output_file() const;        ///< Return the string for the output folder


  /// Set the size of the system of equations
  void set_size(unsigned int size);

  /// Set the initial time value
  void set_initial_time(double t0);

  /// Set the current step size
  void set_dt(double dt);

  /// Set the end time
  void set_end_time(double t);

  /// Set the minimum step size
  void set_dt_min(double dt);

  /// Set the maximum step size
  void set_dt_max(double dt);

  /// Set the maximum function evaluations
  void set_max_eval(unsigned int val);

  /// Set the minimum change in u states
  void set_u_tol(double val);

  /// Set the minimum change in residuals
  void set_f_tol(double val);

  /// Set the minimum change in derivatives
  void set_df_tol(double val);

  /// Set the LevenbergMarquardt error precision
  void set_error_precision(double val);

  /// Set the LevenbergMarquardt diagonal shift bound
  void set_diag_shift_bound(double val);

  /// Set the user scaling option
  void set_user_scaling_option(bool opt);

  /// Set the user integration option
  void select_integration_method(integrate_subtype choice);

  /// Set up default rate based dt selection
  /** This function is used to setup the simple rate-based dt
      stepper that will automatically adjust the selection of the
      time step dt based on the given information.

      Time steps are determined as follows:
          (1) Calculate the magnitudes all individual rates
          (2) If a rate (du/dt) is above the 'upper_limit', then reduce next step (dt_next = dt/factor)
          (3) If rate (du/dt) is below the 'lower_limit', then enlarge next step (dt_next = dt*factor)
          (4) If next dt is lower than dt_min_con, then set to dt_min_con
          (5) Take the smallest dt from the above calculations as the next step

      \param yes If true, this will use the rate-based time-stepper
      \param factor A factor to adjust dt by (dt*factor or dt/factor) if rates are slow or fast
      \param upper_limit A value used to determine if a rate is considered above the limit for a fast rate
      \param lower_limit A value used to determine if a rate is considered below the limit for a slow rate
      \param dt_min_con The minimum allowable time step to take (if previous solution converged)
  */
  void use_default_rate_based_dt(bool yes, double factor=2.0,
                              double upper_limit = 10., double lower_limit = 0.1,
                              double dt_min_con = 1e-4);

  /// Set up standard console output
  void use_standard_console_output(bool yes);

  /// Set up standard file output
  void use_standard_file_output(bool yes, std::string folder = "output", std::string file = "EigenODE_Result");

  /// Get the size of the system of equations
  unsigned int get_size() const;

  /// Function to register user data
  /** We cannot possibly know or anticipate what form user data will
      come in, but for this object that doesn't matter. The user should
      know the data they have and the form it is in. Within any other
      user function provided, the user is responsible for properly
      dereferencing and casting the void pointer to the appropriate
      data structure.
  */
  void register_user_data(const void *data);

  /// Function to register a user's rate function
  /** To create an ODE system, the user must register a series
      of rate functions into this object. Those rate functions
      are then used within the context of the solver to integrate
      the system over time.
  */
  void register_rate_function(int i, double (*func) (int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj), bool is_steady=false );

  /// Function to register a user's coeff function
  /** This is a function to compute a term that gets
      multiplied by an associated derivative. The reason
      this is separated out is to make the interface between
      implicit and explicit ODEs more consistent.
  */
  void register_coeff_function(int i, double (*func) (int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj) );

  /// Function to register a user's scaling function
  /** This is a function to compute a term that gets
      placed into the LevenbergMarquardt scaling vector.
      User's should properly scale their system of equations
      if it is expected that their variables will have very
      different values in magnitudes.
  */
  void register_scaling_function(int i, double (*func) (int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj) );

  /// Function to register a user's initial condition function
  /** This is a function to compute an initial condition for
      the ith variable in the system. This can be as simple
      as a constant value, or can be a function of other information
      the user provides.
  */
  void register_initial_condition_function(int i, double (*func) (int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj) );

  /// Function to register a time-stepper function
  /** This is a function to compute a time step size for each
      cycle of the integrator.
  */
  void register_time_stepper_function(double (*func) (const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj) );

  ///< Function to evaluate user rate function
  double eval_user_rate_function(int i, const Eigen::VectorXd &u, double t) const;

  ///< Function to evaluate user coeff function
  double eval_user_coeff_function(int i, const Eigen::VectorXd &u, double t) const;

  ///< Function to evaluate user scaling function
  double eval_user_scaling_function(int i, const Eigen::VectorXd &u, double t) const;

  ///< Function to evaluate user initial condition function
  double eval_user_initial_condition_function(int i, const Eigen::VectorXd &u, double t) const;

  ///< Function to evaluate the current derivative
  int eval_derivative(const Eigen::VectorXd &u);

  ///< Function to evaluate the current rate
  int eval_rate(const Eigen::VectorXd &u);

  ///< Function to evaluate the coefficients
  int eval_coeff(const Eigen::VectorXd &u);

  ///< Function to evaluate the scaling coefficients
  int eval_scaling(const Eigen::VectorXd &u);

  ///< Function to initialize the system
  int initialize_system(Eigen::LevenbergMarquardt< Eigen::NumericalDiff< DenseFunctorODE > > &lm);

  ///< Function to perform pre-step actions (compute dt, modify scaling, etc.)
  int preprocess_step(Eigen::LevenbergMarquardt< Eigen::NumericalDiff< DenseFunctorODE > > &lm);

  ///< Function to evaluate a step (implicit)
  int solve_step_implicit(Eigen::LevenbergMarquardt< Eigen::NumericalDiff< DenseFunctorODE > > &lm);

  ///< Function to evaluate a step (explicitly)
  int solve_step_explicit(Eigen::LevenbergMarquardt< Eigen::NumericalDiff< DenseFunctorODE > > &lm);

  ///< Function to perform post-step actions (update ports, prepare outputs, etc.)
  int postprocess_step(Eigen::LevenbergMarquardt< Eigen::NumericalDiff< DenseFunctorODE > > &lm);

  ///< Loop function to solve on a specific time interval
  int solve_to_end(Eigen::LevenbergMarquardt< Eigen::NumericalDiff< DenseFunctorODE > > &lm);

  ///< Loop function to solve on a specific time interval (implicit only)
  int solve_to_end_implicit(Eigen::LevenbergMarquardt< Eigen::NumericalDiff< DenseFunctorODE > > &lm);

  ///< Loop function to solve on a specific time interval (explicit only)
  int solve_to_end_explicit(Eigen::LevenbergMarquardt< Eigen::NumericalDiff< DenseFunctorODE > > &lm);

protected:
  /// Set the current step size
  void set_dt_internal(double dt);

private:
  Port<double> time;              ///< Time history (older states)
  double time_new;                ///< Current point in time
  Port<double> dt;                ///< Step size history (older states)
  double dt_new;                  ///< Current step size
  double dt_init;                 ///< Initial time step

  Port<Eigen::VectorXd> u;        ///< Solution vector history (older states)
  Eigen::VectorXd u_new;          ///< Current solution/search vector
  Eigen::VectorXd du;            ///< Current evaluated differences of u states
  Eigen::VectorXd g;             ///< Current evaluated and combined rates
  Eigen::VectorXd co;            ///< Current evaluated and coefficient function
  Eigen::VectorXd scale;         ///< Current evaluated and coefficient function


  /// Types of options
  double dt_min;                ///< Minimum allowable step size
  double dt_max;                ///< Maximum allowable step size
  unsigned int max_eval;        ///< Maximum allowable function evaluations
  double u_tol;                 ///< Allowable tolerance in change in u states
  double f_tol;                 ///< Allowable tolerance in function norm
  double df_tol;                ///< Allowable tolerance in derivative norms
  double error_eps;             ///< Error precision for LevenbergMarquardt
  double step_bound;            ///< Diagonal shift step bound
  bool user_scaling;            ///< Boolean to determine whether to use user scaling
  double end_time;              ///< Time at which to end the simulation

  double dt_rate_factor;        ///< Simple multiplicative factor for changing dt
  double slope_upper_limit;     ///< Upper limit of the slope of rate functions
  double slope_lower_limit;     ///< Lower limit of the slope of rate functions
  double dt_min_converged;      ///< Minimum dt to adjust down to if converged

  /// LevenbergMarquardt Internal Status Tracking
  int lm_status;                ///< Last valid status from LevenbergMarquardt algorithm
  int iter;                     ///< Last valid number of iterations from LevenbergMarquardt
  int f_eval;                   ///< Last valid number of function calls from LevenbergMarquardt
  int df_eval;                  ///< Last valid number of jacobian calls from LevenbergMarquardt
  double fnorm;                 ///< Last valid fnorm from LevenbergMarquardt algorithm
  double gnorm;                 ///< Last valid gnorm from LevenbergMarquardt algorithm
  double lm_param;              ///< Last valid lm_param from LevenbergMarquardt algorithm
  bool is_initalized;           ///< Check to see if model has been initialized

  integrate_type type;          ///< Type of integration scheme
  integrate_subtype sub_type;   ///< Sub-type of integration scheme
  eigen_ode_status status;      ///< Current status of the EigenODE solver object

  std::string s_type;           ///< String version of scheme
  std::string s_sub_type;       ///< String version of sub-type

  DenseFunctorODE functor;        ///< ODE functor object

  const void *user_data;          ///< Any other necessary user data
  FILE *out_file;                 ///< File for writing output to
  std::string folder;             ///< Name of folder for output file
  std::string file;               ///< Name of the file for output

  /// The set of user defined rate functions
  std::vector<double (*) (int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj)> user_func;

  /// The set of user defined coeff functions
  std::vector<double (*) (int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj)> user_coeff;

  /// List of user provided info on steady-state characteristic of the rate function
  std::vector<bool> is_steady_state;

  /// The set of user defined scaling functions
  std::vector<double (*) (int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj)> user_scale;

  /// The set of user defined initial condition functions
  std::vector<double (*) (int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj)> user_ics;

  /// User provided function to calculate next time step
  double (*user_dt) (const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj);

  int (*derivative) (const Eigen::VectorXd &u, EigenODE &obj);  ///< A function to compute a derivative

  int (*rate) (const Eigen::VectorXd &u, EigenODE &obj);        ///< A function to compute a rate

  int (*coeff) (const Eigen::VectorXd &u, EigenODE &obj);        ///< A function to coefficients

  /// Function to call at post-process stage to setup output to console
  int (*init_console_out) (EigenODE &obj);
  int (*post_console_out) (EigenODE &obj);

  /// Function to call at post-process stage to setup output to file
  int (*init_file_out) (EigenODE &obj);
  int (*post_file_out) (EigenODE &obj);

};

/// -------------------- Derivative Function Helpers -----------------------------

// Default derivative function (computes: du = u_new - u)
int default_derivative(const Eigen::VectorXd &u, EigenODE &obj);

// BDF2 derivative function (computes: du = a*u_new - b*u + c*u_old)
int bdf2_derivative(const Eigen::VectorXd &u, EigenODE &obj);


/// -------------------- Rate Function Helpers -----------------------------

// 1st Order Backwards Euler rate function evaluation
int backwards_euler(const Eigen::VectorXd &u, EigenODE &obj);

// 1st Order Forwards Euler rate function evaluation
int forwards_euler(const Eigen::VectorXd &u, EigenODE &obj);

// Trapezoid rule rate function evaluation
int trapezoid_rule(const Eigen::VectorXd &u, EigenODE &obj);

// Runge-Kutta 45 rate function evaluation
int runge_kutta_45(const Eigen::VectorXd &u, EigenODE &obj);


/// -------------------- Coefficient Function Helpers -----------------------------

// Future Coeffs
int future_coeff(const Eigen::VectorXd &u, EigenODE &obj);

// Past Coeffs
int past_coeff(const Eigen::VectorXd &u, EigenODE &obj);


/// -------------------- Default Function Helpers -----------------------------
double default_derivative_coeff(int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj);

double default_variable_ic(int i, const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj);


/// -------------------- Default Time Stepper Helpers -----------------------------
double default_dt(const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj);

double default_rate_based_dt(const Eigen::VectorXd &u, double t, const void *user_data, const EigenODE &obj);


/// -------------------- Default Output Helpers -----------------------------
int default_outputs(EigenODE &obj);

std::string lm_status_message(int lm_status);

int init_console_table_out(EigenODE &obj);

int post_console_table_out(EigenODE &obj);

int init_file_csv_out(EigenODE &obj);

int post_file_csv_out(EigenODE &obj);

/// Test function
int Test_Eigen_ODE();


#endif
