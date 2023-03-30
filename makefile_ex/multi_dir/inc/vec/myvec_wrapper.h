#ifndef _MYVEC_WRAPPER_
#define _MYVEC_WRAPPER_

/*
 * C Wrapper functions for interfacing "myvec" class from Simulink model.
 */
extern void *createMyvec();
extern void deleteMyvec(void *obj);
extern void setSize(void *obj, int size);
extern void setValue(void *obj, int i, int val);
extern int getValue(void *obj, int i);

#endif /* _MYVEC_WRAPPER_ */