#ifndef _MYVEC_OBJ_
#define _MYVEC_OBJ_

#include <vector>
#include "adder/adder_wrapper.h"

// Below is also valid because of makefile struct
//#include "adder_wrapper.h"

class myvec {
private:
    std::vector<int> Data;
    void *p;
public:
    myvec();
    ~myvec();
    void set_data_size(int size);
    void set_data_val(int i, int val);
    int get_data_val(int i);
};

#endif /* _MYVEC_OBJ_ */
