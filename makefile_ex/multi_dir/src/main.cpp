#include <stdio.h>
#include "adder_wrapper.h"
#include "myvec_wrapper.h"
#include <iostream>

int main()
{
    void *obj;
    obj = createAdder();

    for (int i=0; i<5; i++)
    {
        std::cout << adderOutput(obj, i) << std::endl;
    }

    deleteAdder(obj);


    void *vecobj;
    vecobj = createMyvec();

    setSize(vecobj, 5);

    for (int i=0; i<5; i++)
    {
        setValue(vecobj, i, i+i);
        std::cout << getValue(vecobj, i) << std::endl;
    }

    deleteMyvec(vecobj);

    return 0;
}
