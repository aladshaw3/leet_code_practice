/*!
 *  \file port_template.h port_template.cpp
 *  \brief C++ Port Object
 *  \details This file creates a Port object that can be of any type. The Ports
 *          can be any inlet or outlet signal. Also used to track signal history.
 *
 *  \author Austin Ladshaw
 *  \date 11/07/2022
 *  \copyright This software was designed and built by Austin Ladshaw
 *             for ODE simulations. Copyright (c) 2022, all
 *             rights reserved.
 *
 */

#include "port_template.h"

int Test_Port()
{
 int success = 0;

 Port<int> int_port;
 int_port.initializePort(0);

 int_port.setState(1);
 int_port.setState(2);
 int_port.setState(3);

 // Make sure that the current state is correctly tracked
 if (int_port.getState() != 3)
 {
   SIMError(test_failed);
   return -1;
 }

 // Make sure that previous states are tracked
 if (int_port.getState(1) != 2)
 {
   SIMError(test_failed);
   return -1;
 }
 if (int_port.getState(2) != 1)
 {
   SIMError(test_failed);
   return -1;
 }
 // Make sure the last valid state is returned if we did not resize
 if (int_port.getState(3) != 1)
 {
   SIMError(test_failed);
   return -1;
 }

 Port<double> d_port(5);

 d_port.initializePort(2.1);

 d_port.setState(3.1);

 // Make sure that the current state is correctly tracked
 if (d_port.getState() != 3.1)
 {
   SIMError(test_failed);
   return -1;
 }

 if (d_port.getState(4) != 2.1)
 {
   SIMError(test_failed);
   return -1;
 }

 if (d_port.getStateSize() != 5)
 {
   SIMError(test_failed);
   return -1;
 }

 if (d_port.getAvgState() > 2.31 || d_port.getAvgState() < 2.29)
 {
   SIMError(test_failed);
   return -1;
 }



 return success;
}
