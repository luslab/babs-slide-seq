/*
 * nourdine.bah@crick.ac.uk
 */

#include <stdio.h>
#include <CL/cl.h>
#include "error_codes.h"
#include "opencl_program.h"
#include "opencl_kernel.h"

int value(int v) { return v; }

//////////////////////
cl_kernel* new_kernel(
	cl_context context, cl_device_id* devices, char* path, char* name
)
{
	int error;

	cl_program *program = (cl_program*)malloc(sizeof(cl_program));
	cl_kernel *kernel = (cl_kernel*)malloc(sizeof(cl_kernel));

	error = create_prog(program, context, path, name);

	if ( error != 1 ) {
		printf("Error: cannot create the %s program\n", name);
		return FUNCTION_FAILURE;
	}

	error = build_prog(program, devices, name);

	if ( error != FUNCTION_SUCCESS ) {
		printf("Error: cannot build the %s program\n", name);
		return FUNCTION_FAILURE;
	}

	error = create_krn(kernel, program, name);

	if ( error != FUNCTION_SUCCESS ) {
		printf("Error: cannot create the %s kernel\n", name);
		return FUNCTION_FAILURE;
	}

	return kernel;
}
