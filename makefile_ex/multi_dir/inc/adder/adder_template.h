/* Copyright 2022 The MathWorks, Inc. */

#ifndef _ADDER_TEMPLATE_
#define _ADDER_TEMPLATE_

template <typename T>
class adder {
private:
    T state;
public:
    adder();
    T add_one(T increment);
    T get_val();
};

template <typename T>
adder<T>::adder() {
    state = 0;
}

template <typename T>
T adder<T>::add_one(T increment) {
    state += increment;
    return state;
}

template <typename T>
T adder<T>::get_val() {
    return state;
}

#endif /* _ADDER_TEMPLATE_ */
