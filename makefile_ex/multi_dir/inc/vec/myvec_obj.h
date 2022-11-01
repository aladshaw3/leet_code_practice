#ifndef _MYVEC_OBJ_
#define _MYVEC_OBJ_

#include <vector>  

class myvec {
private:
    std::vector<int> Data; 
public:
    myvec();
    ~myvec();
    void set_data_size(int size);
    void set_data_val(int i, int val);
    int get_data_val(int i);
};

#endif /* _MYVEC_OBJ_ */