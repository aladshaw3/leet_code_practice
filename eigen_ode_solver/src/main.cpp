/*!
 *  \file main.cpp
 *  \brief Main Function
 *  \details User input provided at time of execution is used to call the ui functions
 *  \author Austin Ladshaw
 *  \date 11/02/2022
 *  \copyright This software was designed and built by Austin Ladshaw
 *             for ODE simulations. Copyright (c) 2022, all
 *             rights reserved.
 */

#include "eigen_ode.h"

int main(int argc, const char * argv[])
{
  int success = Test_Eigen_ODE();

  std::cout << "Exit Code:\t" << success << std::endl;
  return success;
}
