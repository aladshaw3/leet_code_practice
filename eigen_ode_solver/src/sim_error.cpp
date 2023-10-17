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


#include "sim_error.h"

//Error message to be displayed when exceptions are caught or problems occur
void sim_error(int flag)
{
  switch (flag)
  {
    case generic_error:
      std::cout << "\nUndefined or Generic Error!!!" << std::endl;
      break;

    case invalid_type:
      std::cout << "\nInvalid Variable Type!" << std::endl;
      break;

    case opt_no_support:
      std::cout << "\nOption given currently not supported!" << std::endl;
      break;

    case invalid_console_input:
      std::cout << "\nConsole input provided is invalid!" << std::endl;
      break;

    case test_failed:
      std::cout << "\nA Unit Test has Failed!" << std::endl;
      break;

    case system_failure:
      std::cout << "\nAn OS system failure has occured!" << std::endl;
      break;

    case key_not_found:
      std::cout << "\nGiven key not found in map object!" << std::endl;
      break;

    case ode_failure:
      std::cout << "\nEigenODE object has experienced a failure!" << std::endl;
      break;

    case ode_dt_min_fail:
      std::cout << "\nEigenODE cannot solve the system. Timestep cannot be reduced further!" << std::endl;
      break;

    case ode_invalid_type_selection:
      std::cout << "\nEigenODE cannot solve with EXPLICIT method when a rate function is designated 'steady-state'!" << std::endl;
      std::cout << "\tPlease select an IMPLICIT method instead...\n";
      break;
      
    case index_out_of_bounds:
      std::cout << "\nGiven index is out of bounds for list / vector / set!" << std::endl;
      break;
      
    case generic_fault:
      std::cout << "\nSystem is in a state that it should never have happened..." << std::endl;
      break;

    case invalid_user_value:
      std::cout << "\nUser provided input parameter is invalid!" << std::endl;
      break;

    default:
      std::cout << "\nUndefined or Generic Error!!!" << std::endl;
      break;
  }

}


/// Function to run tests on the SILError methods
int Test_SIMError()
{
  int success = 0;

  SIMError(invalid_type);

  return success;
}
