/*
 * nourdine.bah@crick.ac.uk
 */

#include <stdio.h>
#include <CL/cl.h>
#include "error_codes.h"

///////////////////////
int create_char_ro_buf(
	cl_mem* buf, size_t n, char* name, cl_context context
)
{
	cl_int error;

	*buf = clCreateBuffer(
			context,
			CL_MEM_READ_ONLY,
			n * sizeof(char),
			NULL,
			&error
			);

	if ( error != CL_SUCCESS ) {
		printf("Error: cannot create %s buffer %d\n", name, error);
		return FUNCTION_FAILURE;
	} else {
		return FUNCTION_SUCCESS;
	}
}

///////////////////////
int create_char_rw_buf(
	cl_mem* buf, size_t n, char* name, cl_context context
)
{
	cl_int error;

	*buf = clCreateBuffer(
			context,
			CL_MEM_READ_WRITE,
			n * sizeof(char),
			NULL,
			&error
			);

	if ( error != CL_SUCCESS ) {
		printf("Error: cannot create %s buffer %d\n", name, error);
		return FUNCTION_FAILURE;
	} else {
		return FUNCTION_SUCCESS;
	}
}

//////////////////
int write_char_buf
(
	cl_mem* buf, char* array, cl_command_queue queue, size_t n, char* name
)
{
	cl_int error;
	
	error = clEnqueueWriteBuffer(
			queue,
			*buf, CL_TRUE,
			0,
			n * sizeof(char), array,
			0, NULL, NULL
			);

	if ( error != CL_SUCCESS ) {
		printf("Error: cannot write %s buffer %d\n", name, error);
		return FUNCTION_FAILURE;
	} else {
		return FUNCTION_SUCCESS;
	}
}

/////////////////
int read_char_buf
(
	cl_mem* buf, char* array, cl_command_queue queue, size_t n, char* name
)
{
	cl_int error;

	error = clEnqueueReadBuffer(
			queue, *buf, CL_TRUE,
			0, n * sizeof(char), array,
			0, NULL, NULL);

	if ( error != CL_SUCCESS ) {
		printf("Error: cannot read %s buffer %d\n", name, error);
		return FUNCTION_FAILURE;
	} else {
		return FUNCTION_SUCCESS;
	}
}

///////////////////////
int create_llu_rw_buf(
	cl_mem* buf, size_t n, char* name, cl_context context
)
{
	cl_int error;

	*buf = clCreateBuffer(
			context,
			CL_MEM_READ_WRITE,
			n * sizeof(unsigned long long int),
			NULL,
			&error
			);

	if ( error != CL_SUCCESS ) {
		printf("Error: cannot create %s buffer %d\n", name, error);
		return FUNCTION_FAILURE;
	} else {
		return FUNCTION_SUCCESS;
	}
}

///////////////////
int write_llu_buf
(
	cl_mem* buf, unsigned long long int* array, cl_command_queue queue,
	size_t n, char* name
)
{
	cl_int error;
	
	error = clEnqueueWriteBuffer(
			queue,
			*buf, CL_TRUE,
			0,
			n * sizeof(unsigned long long int), array,
			0, NULL, NULL
			);

	if ( error != CL_SUCCESS ) {
		printf("Error: cannot write %s buffer %d\n", name, error);
		return FUNCTION_FAILURE;
	} else {
		return FUNCTION_SUCCESS;
	}
}

////////////////
int read_llu_buf
(
	cl_mem* buf, unsigned long long int* array, cl_command_queue queue,
	size_t n, char* name
)
{
	cl_int error;

	error = clEnqueueReadBuffer(
			queue, *buf, CL_TRUE,
			0, n * sizeof(unsigned long long int), array,
			0, NULL, NULL);

	if ( error != CL_SUCCESS ) {
		printf("Error: cannot read %s buffer %d\n", name, error);
		return FUNCTION_FAILURE;
	} else {
		return FUNCTION_SUCCESS;
	}
}

