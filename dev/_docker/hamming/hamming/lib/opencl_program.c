/*
 * nourdine.bah@crick.ac.uk
 */

#include <stdio.h>
#include <stdlib.h>

#include <CL/cl.h>

#include "error_codes.h"
#include "utils.h"
#include "opencl_kernel.h"

////////////////
int create_prog(
	cl_program *program, cl_context context, char *path, char *name
)
{
	cl_int error;

	// the source code and the its size
	char *src = get_kernel_source(path);
	size_t size = get_file_size(path);

	*program = clCreateProgramWithSource(
			context,
			1,
			(const char**)&src,
			(const size_t*)&size,
			&error
			);

	if ( error != CL_SUCCESS ) {
		printf("Error: cannot create %s program %d\n", name, error);
		return FUNCTION_FAILURE;
	} else {
		return FUNCTION_SUCCESS;
	}
}

///////////////
int build_prog(
	cl_program *program, cl_device_id *devices, char *name
)
{

	cl_int error;

	error = clBuildProgram(*program, 1, devices, NULL, NULL, NULL);

	if ( error != CL_SUCCESS ) {
		printf("Error: cannot build %s program %d\n", name, error);
		return FUNCTION_FAILURE;
	} else {
		return FUNCTION_SUCCESS;
	}
}
