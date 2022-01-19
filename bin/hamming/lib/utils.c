/*
 * nourdine.bah@crick.ac.uk
 */

#include <stdio.h>
#include <stdlib.h>
#include "error_codes.h"

/////////////////////
size_t get_file_size(
	char* path
)
{

	FILE *fp = fopen(path, "r");

	if (!fp) {
		printf("Error: cannot open the kernel source file\n");
		return FUNCTION_FAILURE;
	}

	fseek(fp, 0L, SEEK_END);
	size_t sz = ftell(fp);

	return sz;
}
