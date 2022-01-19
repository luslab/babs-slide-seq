/*
 * nourdine.bah@crick.ac.uk
 */

#include <CL/cl.h>

char* get_kernel_source(char* kernel_path);
int create_krn(cl_kernel *kernel, cl_program *program, char *name);
int set_arg(
	cl_kernel *kernel, cl_mem *buf, cl_int pos, char *krn_name, char *arg_name
);

