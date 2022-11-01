#include "myvec_wrapper.h"
#include "myvec_obj.h"

void *createMyvec()
{
    return new myvec();
}

void deleteMyvec(void *obj)
{
    delete static_cast<myvec *>(obj);
}

void setSize(void *obj, int size)
{
    myvec *dat = static_cast<myvec *>(obj);
    dat->set_data_size(size);
}

void setValue(void *obj, int i, int val)
{
    myvec *dat = static_cast<myvec *>(obj);
    dat->set_data_val(i, val);
}

int getValue(void *obj, int i)
{
    myvec *dat = static_cast<myvec *>(obj);
    return dat->get_data_val(i);
}