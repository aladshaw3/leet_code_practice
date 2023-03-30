/* Copyright 2022 The MathWorks, Inc. */

#include "adder_wrapper.h"
#include "adder_template.h"

void *createAdder()
{
    return new adder<int>();
}

void deleteAdder(void *obj)
{
    delete static_cast<adder<int> *>(obj);
}

int adderOutput(void *obj, int increment)
{
    adder<int> *adderObj = static_cast<adder<int> *>(obj);
    adderObj->add_one(increment);
    return adderObj->get_val();
}
