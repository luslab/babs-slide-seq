/*
 * nourdine.bah@crick.ac.uk
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <CL/cl.h>

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <tuple>
#include <numeric>
#include <algorithm>

#include "hamming.hpp"

const int MAX_BARCODES = 10;

struct Match
{
	Sequence read;
	char dist;
	int n;
	Sequence barcodes[MAX_BARCODES];
};

int main(int argc, char **argv)
{
	int error;

	///////////////////////////////////////////////////////////////////////////
	// READS AND BARCODES

	std::vector<std::string> reads = get_sequences(argv[1]);
	std::vector<std::string> barcodes = get_sequences(argv[2]);

	std::cout << reads[0] << ", " << barcodes[0] << std::endl;

	int N_READS = reads.size();
	int N_BARCODES = barcodes.size();

	///////////////////////////////////////////////////////////////////////////
	// LENGTH

	char LENGTH = 0;
	error = get_length(reads, barcodes, &LENGTH);
	if ( -1 == error ) {
		printf("Error: reads barcodes don't all have the same length");
		return FUNCTION_FAILURE;
	}
	if ( -2 == error ) {
		printf("Error: puck barcodes don't all have the same length");
		return FUNCTION_FAILURE;
	}
	if ( -3 == error ) {
		printf("Error: read and puck barcodes have different length");
		return FUNCTION_FAILURE;
	}

	std::cout << "Pass 1" << std::endl;

	///////////////////////////////////////////////////////////////////////////
	// MATCH RECORDS
	
	Match* matches = (Match*)malloc(N_READS*sizeof(Match));
	std::cout << "Pass 2" << std::endl;
	
	///////////////////////////////////////////////////////////////////////////
	// CONVERSION TO DIGITS
	
	Sequence* numbers = (Sequence*)malloc(sizeof(Sequence)*N_BARCODES);
	for (int i=0; i<N_BARCODES; i++)
	{
		numbers[i] = seq_to_num(barcodes[i]);
	}
	std::cout << "Pass 3" << std::endl;

	///////////////////////////////////////////////////////////////////////////
	// GPU

	cl_int err;
	cl_platform_id *platforms;
	cl_uint *num_platforms;
	cl_device_id *devices;
	cl_uint *num_devices;
	cl_context context;
	cl_command_queue queue;
	cl_uint maxD = 1; // max number of platforms to use
	cl_uint maxP = 1; // max number of devices to use

	// alloc
	platforms = (cl_platform_id*)malloc( maxP * sizeof(cl_platform_id) );
	num_platforms = (cl_uint*)malloc( maxP * sizeof(cl_uint) );
	devices = (cl_device_id*)malloc( maxD * sizeof(cl_device_id) );
	num_devices = (cl_uint*)malloc( maxD * sizeof(cl_uint) );

	// platform
	err = clGetPlatformIDs(maxP, platforms, num_platforms);
	if ( err != CL_SUCCESS ) {
		printf("Error: cannot get platform ID %d\n", err);
		return -1;
	}

	// device
	err = clGetDeviceIDs(
			platforms[0], CL_DEVICE_TYPE_GPU, maxD, devices, num_devices);
	if ( err != CL_SUCCESS ) {
		printf("Error: cannot get device ID %d\n", err);
		return -1;
	}

	// context
	context = clCreateContext(NULL, 1, devices, NULL, NULL, &err);
	if ( err != CL_SUCCESS ) {
		printf("Error: cannot create context %d\n", err);
		return -1;
	}

	// command queue
	queue = clCreateCommandQueueWithProperties(context, devices[0], NULL, &err);
	if ( err != CL_SUCCESS ) {
		printf("Error: cannot create command queue %d\n", err);
		return -1;
	}

	///////////////////////////////////////////////////////////////////////////
	// PROGRAMS
	
	cl_kernel* hamming_krn =
		new_kernel(context, devices, "cl/hamming.cl", "hamming");

	cl_kernel* min_krn =
		new_kernel(context, devices, "cl/min_dist.cl", "min_dist");

	///////////////////////////////////////////////////////////////////////////
	// SOME BUFFERS FOR HAMMING KERNEL

	// sequences numbers
	cl_mem* nbuf;
	nbuf = (cl_mem*)malloc(sizeof(cl_mem));
	error = create_llu_rw_buf(nbuf, (size_t)N_BARCODES, "nbuf", context);
	if ( error != FUNCTION_SUCCESS ) {
		printf("Error: cannot create %s buffer\n", "nbuf");
		return FUNCTION_FAILURE;
	}
	error = write_llu_buf(nbuf, numbers, queue, N_BARCODES, "nbuf");
	if ( error != FUNCTION_SUCCESS ) {
		printf("Error: cannot write %s buffer\n", "nbuf");
		return FUNCTION_FAILURE;
	}
	err = set_arg(hamming_krn, nbuf, 0, "hamming_krn", "nbuf");
	if ( err != FUNCTION_SUCCESS ) {
		return FUNCTION_FAILURE;
	}

	// the current read
	cl_mem* rbuf;
	rbuf = (cl_mem*)malloc(sizeof(cl_mem*));
	error = create_llu_rw_buf(rbuf, 1, "rbuf", context);
	if ( error != FUNCTION_SUCCESS ) {
		printf("Error: cannot create %s buffer\n", "rbuf");
		return FUNCTION_FAILURE;
	}
	err = set_arg(hamming_krn, rbuf, 1, "hamming_krn", "rbuf");
	if ( err != FUNCTION_SUCCESS ) {
		return FUNCTION_FAILURE;
	}

	// distances for the current read
	cl_mem* dbuf;
	dbuf = (cl_mem*)malloc(sizeof(cl_mem));
	error = create_char_rw_buf(dbuf, (size_t)N_BARCODES, "dbuf", context);
	if ( error != FUNCTION_SUCCESS ) {
		printf("Error: cannot create %s buffer\n", "dbuf");
		return FUNCTION_FAILURE;
	}
	err = set_arg(hamming_krn, dbuf, 2, "hamming_krn", "dbuf");
	if ( err != FUNCTION_SUCCESS ) {
		return FUNCTION_FAILURE;
	}

	// sequence length
	cl_mem* lbuf;
	lbuf = (cl_mem*)malloc(sizeof(cl_mem*));
	error = create_char_rw_buf(lbuf, 1, "lbuf", context);
	if ( error != FUNCTION_SUCCESS ) {
		printf("Error: cannot create %s buffer\n", "lbuf");
		return FUNCTION_FAILURE;
	}
	error = write_char_buf(lbuf, &LENGTH, queue, 1, "lbuf");
	if ( error != FUNCTION_SUCCESS ) {
		printf("Error: cannot write %s buffer\n", "lbuf");
		return FUNCTION_FAILURE;
	}
	err = set_arg(hamming_krn, lbuf, 3, "hamming_krn", "lbuf");
	if ( err != FUNCTION_SUCCESS ) {
		return FUNCTION_FAILURE;
	}

	///////////////////////////////////////////////////////////////////////////
	// SOME BUFFERS FOR MIN KERNEL
	
	// distances
	err = set_arg(min_krn, dbuf, 0, "min_krn", "dbuf");
	if ( err != FUNCTION_SUCCESS ) {
		return FUNCTION_FAILURE;
	}

	// current distance
	cl_mem* cbuf;
	cbuf = (cl_mem*)malloc(sizeof(cl_mem*));
	error = create_char_rw_buf(cbuf, 1, "cbuf", context);
	if ( error != FUNCTION_SUCCESS ) {
		printf("Error: cannot create %s buffer\n", "cbuf");
		return FUNCTION_FAILURE;
	}
	err = set_arg(min_krn, cbuf, 1, "min_krn", "cbuf");
	if ( err != FUNCTION_SUCCESS ) {
		return FUNCTION_FAILURE;
	}

	// min distance
	cl_mem* mbuf;
	mbuf = (cl_mem*)malloc(sizeof(cl_mem*));
	error = create_char_rw_buf(mbuf, 1, "mbuf", context);
	if ( error != FUNCTION_SUCCESS ) {
		printf("Error: cannot create %s buffer\n", "mbuf");
		return FUNCTION_FAILURE;
	}
	err = set_arg(min_krn, mbuf, 2, "min_krn", "mbuf");
	if ( err != FUNCTION_SUCCESS ) {
		return FUNCTION_FAILURE;
	}

	///////////////////////////////////////////////////////////////////////////
	
	size_t global_item_size = (size_t)N_BARCODES;
	size_t local_item_size = (size_t)1;

	char min_dist;
	char* flags;
	flags = (char*)malloc(N_BARCODES*sizeof(char));

	for (int i=0; i<reads.size(); i++)
	{
		// read
		Sequence read = seq_to_num(reads[i]);
		error = write_llu_buf(rbuf, &read, queue, 1, "rbuf");
		if ( error != FUNCTION_SUCCESS ) {
			printf("Error: cannot write %s buffer\n", "rbuf");
			return FUNCTION_FAILURE;
		}

		error = clEnqueueNDRangeKernel(
			queue, *hamming_krn,
			1, // dimension
			NULL, // global work offset
			&global_item_size, // global work size
			&local_item_size, // local work size
			0, NULL, NULL);

		// error handling
		if ( error != CL_SUCCESS ) {
			printf(
				"Error: cannot execute %s kernel %d\n",
				"hamming_krn", error
			);
			return FUNCTION_FAILURE;
		}

		///////////////////////////////////////////////////////////////////////

		min_dist = LENGTH + 1;
		error = write_char_buf(mbuf, &min_dist, queue, 1, "mbuf");
		if ( error != FUNCTION_SUCCESS ) {
			printf("Error: cannot write %s buffer\n", "mbuf");
			return FUNCTION_FAILURE;
		}
		
		for (char d=0; d<LENGTH+1; d++)
		{
			error = write_char_buf(cbuf, &d, queue, 1, "cbuf");
			if ( error != FUNCTION_SUCCESS ) {
				printf("Error: cannot write %s buffer\n", "cbuf");
				return FUNCTION_FAILURE;
			}

			error = clEnqueueNDRangeKernel(
				queue, *min_krn,
				1, // dimension
				NULL, // global work offset
				&global_item_size, // global work size
				&local_item_size, // local work size
				0, NULL, NULL);

			// error handling
			if ( error != CL_SUCCESS ) {
				printf(
					"Error: cannot execute %s kernel %d\n",
					"min_krn", error
				);
				return FUNCTION_FAILURE;
			}
		}

		error = clEnqueueReadBuffer(
				queue, *mbuf, CL_TRUE,
				0, 1*sizeof(char), &min_dist,
				0, NULL, NULL
		);
		if ( error != CL_SUCCESS ) {
			printf("Error: cannot read %s buffer\n", "mbuf");
			return FUNCTION_FAILURE;
		}

		///////////////////////////////////////////////////////////////////////
		char* dists;
		dists = (char*)malloc(sizeof(char)*N_BARCODES);
		error = read_char_buf(dbuf, dists, queue, (size_t)N_BARCODES, "dbuf");
		if ( error != FUNCTION_SUCCESS ) {
			printf("Error: cannot read %s buffer %d\n", "dbuf", error);
			return FUNCTION_FAILURE;
		}

		matches[i].read = read;
		matches[i].dist = min_dist;
		matches[i].n = 0;

		for (int j=0; j<N_BARCODES; j++)
		{
			if ( dists[j] == min_dist )
			{
				matches[i].barcodes[ matches[i].n ] = numbers[j];
				matches[i].n++;
			}
			if ( matches[i].n == MAX_BARCODES )
			{
				break;
			}
		}

		free(dists);
		///////////////////////////////////////////////////////////////////////

		if ( i % 1000 == 0 ) { std::cerr << i << "/" << N_READS << std::endl; }
		//if ( i == 10000 ) { break; } // test
	}

	///////////////////////////////////////////////////////////////////////////
	// OUTPUT

	std::ofstream csv(argv[3]);

	for (int i=0; i<N_READS; i++)
	{
		std::string bcd_str = "";
		for (int j=0; j<matches[i].n; j++)
		{
			if ( j > 0 )
			{
				bcd_str += ":";
			}
			bcd_str += num_to_seq( matches[i].barcodes[j] , LENGTH );
		}

		csv 
			<< num_to_seq( matches[i].read , LENGTH )
			<< ","
			<< +matches[i].dist
			<< ","
			<< matches[i].n
			<< ","
			<< bcd_str
			<< std::endl;
	}

	csv.close();

	///////////////////////////////////////////////////////////////////////////
	return 0;
}
