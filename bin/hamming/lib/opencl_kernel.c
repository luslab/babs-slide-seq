/*
 * nourdine.bah@crick.ac.uk
 */

#include <stdio.h>
#include <stdlib.h>

#include <CL/cl.h>

#include "error_codes.h"

////////////////////////
char* get_kernel_source(
	char* kernel_path
)
{

	FILE *fp = fopen(kernel_path, "r");

	if (!fp) {
		printf("Error: cannot open the kernel source file\n");
		return NULL;
	}

	// get source file size
	fseek(fp, 0L, SEEK_END);
	size_t src_size = ftell(fp);
	fseek(fp, 0L, SEEK_SET);
	char *kernel_src = (char*)malloc(src_size);

	fread(kernel_src, src_size, 1, fp);

	return kernel_src;
}

///////////////
int create_krn(
	cl_kernel *kernel, cl_program *program, char *name
)
{

	cl_int error;

	*kernel = clCreateKernel(*program, name, &error);

	if ( error != CL_SUCCESS ) {
		printf("Error: cannot create %s kernel %d\n", name, error);
		return FUNCTION_FAILURE;
	} else {
		return FUNCTION_SUCCESS;
	}
}

////////////
int set_arg(
	cl_kernel *kernel, cl_mem *buf, cl_int pos, char *krn_name, char *arg_name
)
{
	cl_int error;

	error = clSetKernelArg(*kernel, pos, sizeof(cl_mem), (void*)buf);

	if ( error != CL_SUCCESS ) {
		printf(
				"Error: cannot set %s %s kernel argument %d\n",
				arg_name, krn_name, error);
		return FUNCTION_FAILURE;
	} else {
		return FUNCTION_SUCCESS;
	}
}
