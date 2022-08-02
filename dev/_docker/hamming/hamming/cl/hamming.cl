/*
 * nourdine.bah@crick.ac.uk
 */

 __kernel void hamming(
		 __global const unsigned long * barcodes,
		 __global unsigned long * sequence,
		 __global char* distances,
		 __global const char* seq_length
		 )
{
	int gid = get_global_id(0);

	unsigned long barcode = barcodes[gid];

	char dist = 0;
	unsigned long base1;
	unsigned long base2;
	unsigned long mask = 7;

	for (int i=0; i<seq_length[0]; i++)
	{
		base1 = *sequence & ( mask << ( i * 3 ) );
		base2 = barcode & ( mask << ( i * 3 ) );

		if ( base1 != base2 )
		{
			dist++;
		}
	}

	distances[gid] = dist;
}
