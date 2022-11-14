// Sokolova Ekaterina, 627 group
// Task 2, variant 13

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

#define MASTER 0
#define INTEGRAL_I (1.0 / 24)
#define POINTS_ON_ITER 1000

double func_F(double x, double y, double z)
{
	// we need count function only in G else F = 0
	if(x < -1 || x > 0 || y < -1 || y > 0 || z < -1 || z > 0)
		return 0;

		return (x * x * x) * (y * y) * z;
}


int main(int argc, char **argv)
{
	// single parameter is required accuracy (epsilon)
	double accur = atof(argv[1]);

	// standart variables
	int rank;
	int size;
	double ltime;
	double part_sum;
	int time_to_stop = 0;

	// initialization
	double I_0 = 0.0;
	double inaccur = 1000;
	int iter = 0;
	double t;

	// parallelization
	if(MPI_Init(&argc, &argv) != MPI_SUCCESS)
	{
		printf("MPI_Init failed\n");
		return 1;
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	unsigned int seed = rank;
	ltime = MPI_Wtime();

	if(rank == MASTER)
	{
		double summ = 0.0;
		double points[3 * POINTS_ON_ITER];

		int sendcounts[size];
		int displs[size];

		sendcounts[MASTER] = 0;
		displs[MASTER] = 0;
		// we use information about value of MASTER
		for(int i = 1; i < size; i++)
		{
			sendcounts[i] = POINTS_ON_ITER / (size - 1) + (i > POINTS_ON_ITER % (size - 1) ? 0 : 1);
			sendcounts[i] *= 3;
			displs[i] = displs[i - 1] + sendcounts[i - 1];
		}

		while(inaccur > accur)
		{
			// 1) generate
			for(int i = 0; i < 3 * POINTS_ON_ITER; i += 3)
			{
				// generate single point (x, y, z) in G
				double x = (double) rand_r(&seed) / (double) RAND_MAX - 1;
				double y = (double) rand_r(&seed) / (double) RAND_MAX - 1;
				double z = (double) rand_r(&seed) / (double) RAND_MAX - 1;

				// fill array with all points
				points[i] = x;
				points[i + 1] = y;
				points[i + 2] = z;
			}

			// 2) send to others
			MPI_Scatterv(points, sendcounts, displs, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

			// 3) receive result
			part_sum = 0.0;
			MPI_Reduce(&part_sum, &summ, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			I_0 = (I_0 * iter + summ / POINTS_ON_ITER) / (iter + 1);
			iter++;
			inaccur = fabs(INTEGRAL_I - I_0);
		}
		time_to_stop = 1;
		for(int i = 0; i < 3 * POINTS_ON_ITER; i++)
			points[i] = 1;
		MPI_Scatterv(points, sendcounts, displs, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
		ltime = MPI_Wtime() - ltime;
	}
	else
	{
		int recvcount = (POINTS_ON_ITER / (size - 1) + (rank > POINTS_ON_ITER % (size - 1) ? 0 : 1)) * 3;
		double local_points[recvcount];

		while(1)
		{
			// 1) receive
			MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, local_points, recvcount, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

			// 2) check is it time to stop
			if(local_points[0] == 1)
				break;

			// 3) count
			part_sum = 0.0;
			for(int i = 0; i < recvcount / 3; i++)
				part_sum += func_F(local_points[i], local_points[i + 1], local_points[i + 2]);

			// 4) send local part of result
			MPI_Reduce(&part_sum, NULL, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
		}

		ltime = MPI_Wtime() - ltime;
	}

	if(rank == MASTER)
	{
		MPI_Reduce(MPI_IN_PLACE, &t, 1, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
		printf("I_0: %.10f\n", I_0);
		printf("|I - I_0|: %.10f\n", inaccur);
		printf("n: %d\n", POINTS_ON_ITER * iter);
		printf("t: %.10fs\n", t);
	}
	else
		MPI_Reduce(&ltime, NULL, 1, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);

	if(MPI_Finalize() != MPI_SUCCESS)
	{
		printf("MPI_Finalize failed\n");
		return 1;
	}

	return 0;
}