/*
 * nourdine.bah@crick.ac.uk
 */

#include <CL/cl.h>

int create_char_ro_buf(cl_mem* buf, size_t n, char *name, cl_context context);
int create_char_rw_buf(cl_mem* buf, size_t n, char *name, cl_context context);
int write_char_buf(cl_mem* buf, char* array, cl_command_queue queue, size_t n, char *name);
int read_char_buf(cl_mem* buf, char* array, cl_command_queue queue, size_t n, char *name);
int create_llu_rw_buf(cl_mem* buf, size_t n, char* name, cl_context context);
int write_llu_buf(cl_mem* buf, unsigned long long int* array, cl_command_queue queue, size_t n, char* name);
