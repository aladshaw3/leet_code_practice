/*!
 *  \file sim_error.h sim_error.cpp
 *  \brief All error types are defined here
 *  \details This file defines all the different errors that may occur in any simulation
 *      in any file. 
 *  \author Austin Ladshaw
 *  \date 11/02/2022
 *  \copyright This software was designed and built by Austin Ladshaw
 *             for ODE simulations. Copyright (c) 2022, all
 *             rights reserved.
 */

#pragma once
#ifndef _SIM_ERROR_
#define _SIM_ERROR_

#include <iostream>          //Line to allow for read/write to the console using cpp functions

#ifndef SIMError
#define SIMError(i) \
{sim_error(i);         \
std::cout << "Source: " << __FILE__ << "\nLine: " << __LINE__ << std::endl;}
#endif

/// List of names for error type
typedef enum
{
generic_error,
invalid_type,
opt_no_support,
invalid_console_input,
test_failed,
system_failure,
key_not_found,
ode_failure,
ode_dt_min_fail,
ode_invalid_type_selection,
invalid_user_value,
index_out_of_bounds,
generic_fault
} error_type;

/// Error function customizes output message based on flag
/** This error function is reference in the error.cpp file, but is not called by any other
  file. Instead, all other files call the mError(i) macro that expands into this error
  function call plus prints out the file name and line number where the error occured. */
void sim_error(int flag);


/// Function to run tests on the SILError methods
int Test_SIMError();

#endif
