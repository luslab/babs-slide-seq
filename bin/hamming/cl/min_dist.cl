/*
 * nourdine.bah@crick.ac.uk
 */

 __kernel void min_dist(
		 __global char* distances,
		 __global char* cur_dist,
		 __global char* min_dist
		 )
{
	int gid = get_global_id(0);

	if ( cur_dist[0] <= min_dist[0] )
	{
		if ( distances[gid] == cur_dist[0] )
		{
			min_dist[0] = cur_dist[0];
		}
	}
}
