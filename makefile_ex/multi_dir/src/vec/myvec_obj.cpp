#include "myvec_obj.h"

/**********************************/
/**** Class method definitions ****/
/**********************************/

myvec::myvec()
{
    this->Data.resize(1);
}

myvec::~myvec()
{
    this->Data.clear();
}

void myvec::set_data_size(int size)
{
    this->Data.resize(size);
}

void myvec::set_data_val(int i, int val)
{
    this->Data[i]=val;
}

int myvec::get_data_val(int i)
{
    return this->Data[i];
}