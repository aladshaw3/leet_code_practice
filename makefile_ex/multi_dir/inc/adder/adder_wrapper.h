/* Copyright 2022 The MathWorks, Inc. */

#ifndef _ADDER_WRAPPER_
#define _ADDER_WRAPPER_

/*
 * C Wrapper functions for interfacing "adder" class from Simulink model.
 */
extern void *createAdder();
extern void deleteAdder(void *obj);
extern int adderOutput(void *obj, int increment);

#endif /* _ADDER_WRAPPER_ */
