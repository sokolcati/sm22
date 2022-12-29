// Sokolova Ekaterina, 627 group
// Task 3, variant 7

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <vector>
#include <iostream>

using namespace std;

#define CONST_K 1000
#define CONST_T 1.0
#define STEPS 20

#define AXIS_X 0
#define AXIS_Y 1
#define AXIS_Z 2

struct Block
{
	int min[3];
	int max[3];
	int len[3];
	int size;
};

struct idBlock
{
	int id;
	Block block;
};

typedef struct Block Block;
typedef struct idBlock idBlock;
typedef vector<double> Data;

struct Grid
{
// we may use just L and h but create 3 variables for visibility
	double L_x; // dimensional size
	double L_y;
	double L_z;

	double h_x; // dimensional step
	double h_y;
	double h_z;
	int N; // amount of dimensional points

	double T; // time range
	double ts; // time steps
	int K; // amount of time points

	Data u_vec[3];
	vector<idBlock> send;
    vector<idBlock> receive;
};

typedef struct Grid Grid;

Block newBlock(int x_min, int x_max, int y_min, int y_max, int z_min, int z_max)
{
	Block newb;

	newb.min[0] = x_min; newb.max[0] = x_max;
	newb.min[1] = y_min; newb.max[1] = y_max;
	newb.min[2] = z_min; newb.max[2] = z_max;

	newb.len[0] = x_max - x_min + 1;
	newb.len[1] = y_max - y_min + 1;
	newb.len[2] = z_max - z_min + 1;

	newb.size = newb.len[0] * newb.len[1] * newb.len[2];

	return newb;
}

idBlock newIdBlock(int i, int x_min, int x_max, int y_min, int y_max, int z_min, int z_max)
{
	idBlock newb;

	newb.block = newBlock(x_min, x_max, y_min, y_max, z_min, z_max);
	newb.id = i;

	return newb;
}

int block_counter;

void save_block(Block * blocks, int x_min, int x_max, int y_min, int y_max, int z_min, int z_max)
{
	blocks[block_counter] = newBlock(x_min, x_max, y_min, y_max, z_min, z_max);
	block_counter++;
}

double u_analytical(double x, double y, double z, double t, Grid omega)
{
	double a = M_PI * sqrt(4.0 / (omega.L_x * omega.L_x) + 4.0 / (omega.L_y * omega.L_y) + 1.0 / (omega.L_z * omega.L_z));
	return sin(2 * M_PI * x / omega.L_x + 3 * M_PI) * sin(2 * M_PI * y / omega.L_y + 2 * M_PI) * sin(M_PI * z / omega.L_z) * cos(a * t + M_PI);
}

int nested(Block block1, Block block2, int i, int j)
{
	return block2.min[i] <= block1.min[i] && block1.max[i] <= block2.max[i] &&
		   block2.min[j] <= block1.min[j] && block1.max[j] <= block2.max[j];
}

int index_converter(int x, int y, int z, Block block)
{
	return (x - block.min[0]) * block.len[1] * block.len[2] + (y - block.min[1]) * block.len[2] + (z - block.min[2]);
}

void split(int x_min, int x_max, int y_min, int y_max, int z_min, int z_max, int axis, int size, Block * blocks)
{
	if (size == 1)
	{
		save_block(blocks, x_min, x_max, y_min, y_max, z_min, z_max);
		return;
	}

	int x = x_min + (x_max - x_min) / size;
	int y = y_min + (y_max - y_min) / size;
	int z = z_min + (z_max - z_min) / size;

	if (size % 2 == 1)
	{
		switch(axis)
		{
			case AXIS_X:
				save_block(blocks, x_min, x, y_min, y_max, z_min, z_max);
				x_min = x + 1;
				axis = AXIS_Y;
				break;
			case AXIS_Y:
				save_block(blocks, x_min, x_max, y_min, y, z_min, z_max);
				y_min = y + 1;
				axis = AXIS_Z;
				break;
			case AXIS_Z:
				save_block(blocks, x_min, x_max, y_min, y_max, z_min, z);
				z_min = z + 1;
				axis = AXIS_X;
				break;
		}
		size--;
	}

	x = (x_min + x_max) / 2;
	y = (y_min + y_max) / 2;
	z = (z_min + z_max) / 2;

	switch(axis)
	{
		case AXIS_X:
			split(x_min, x, y_min, y_max, z_min, z_max, AXIS_Y, size / 2, blocks);
			split(x + 1, x_max, y_min, y_max, z_min, z_max, AXIS_Y, size / 2, blocks);
			break;
		case AXIS_Y:
			split(x_min, x_max, y_min, y, z_min, z_max, AXIS_Z, size / 2, blocks);
			split(x_min, x_max, y + 1, y_max, z_min, z_max, AXIS_Z, size / 2, blocks);
			break;
		case AXIS_Z:
			split(x_min, x_max, y_min, y_max, z_min, z, AXIS_X, size / 2, blocks);
			split(x_min, x_max, y_min, y_max, z + 1, z_max, AXIS_X, size / 2, blocks);
			break;
	}
}

vector<Data> messenger(int id, Block &block, Grid grid)
{
	Block bl;
	Data sent(bl.size);
	vector<Data> rec(grid.receive.size());
	vector<MPI_Request> request(2);
	vector<MPI_Status> status(2);

	for (int i = 0; i < grid.receive.size(); i++)
	{
cout << "prepare recv\n";
cout << grid.receive[i].block.size << "---\n";
		Data tmp = vector<double>(grid.receive[i].block.size);
cout << "not tmp\n";
		rec[i] = tmp;;
cout << "recv\n";
		bl = grid.send[i].block;
		#pragma omp parallel for
		for (int x = bl.min[0]; x <= bl.max[0]; x++)
			for (int y = bl.min[1]; y <= bl.max[1]; y++)
				for (int z = bl.min[2]; z <= bl.max[2]; z++)
					sent[index_converter(x, y, z, bl)] = grid.u_vec[id][index_converter(x, y, z, block)];
		MPI_Isend(sent.data(), bl.size, MPI_DOUBLE, grid.send[i].id, 0, MPI_COMM_WORLD, &request[0]);
		MPI_Irecv(rec[i].data(), grid.receive[i].block.size, MPI_DOUBLE, grid.receive[i].id, 0, MPI_COMM_WORLD, &request[1]);
		MPI_Waitall(2, request.data(), status.data());
	}
	return rec;
}

double value(int id, int x, int y, int z, vector<Data> &recieved, Block &block, Grid grid)
{

	if (block.min[0] <= x and x <= block.max[0] and block.min[1] <= y and y <= block.max[1] and block.min[2] <= z and z <= block.max[2])
		return grid.u_vec[id][index_converter(x, y, z, block)];

	for (int r_i = 0; r_i < grid.receive.size(); r_i++)
	{
		Block block2 = grid.receive[r_i].block;
		if (x < block2.min[0] or x > block2.max[0] or
			y < block2.min[1] or y > block2.max[1] or
			z < block2.min[2] or z > block2.max[2])
			continue;

		return recieved[r_i][index_converter(x, y, z, block2)];
	}

	return 0;
}

double analogue_laplace(int id, int x, int y, int z, vector<Data> &rec, Block &block, Grid grid)
{
	double result = 0;
	result += (value(id, x,   y-1, z,   rec, block, grid) - 2 * grid.u_vec[id][index_converter(x, y, z, block)] +
			   value(id, x,   y+1, z,   rec, block, grid)) / (grid.h_x * grid.h_x);
	result += (value(id, x-1, y,   z,   rec, block, grid) - 2 * grid.u_vec[id][index_converter(x, y, z, block)] +
			   value(id, x+1, y,   z,   rec, block, grid)) / (grid.h_y * grid.h_y);
	result += (value(id, x,   y,   z-1, rec, block, grid) - 2 * grid.u_vec[id][index_converter(x, y, z, block)] +
			   value(id, x,   y,   z+1, rec, block, grid)) / (grid.h_z * grid.h_z);
	return result;
}

void border(int id, Block block, Grid grid)
{
	// Axis X. Periodic condition

	if (block.min[0] == 0)
		#pragma omp parallel for
		for (int y = block.min[1]; y <= block.max[1]; y++)
			for (int z = block.min[2]; z <= block.max[2]; z++)
				grid.u_vec[id][index_converter(block.min[0], y, z, block)] =
					u_analytical(block.min[0] * grid.h_x, y * grid.h_y, z * grid.h_z, grid.ts, grid);

	if (block.max[0] == grid.N)
		#pragma omp parallel for
		for (int y = block.min[1]; y <= block.max[1]; y++)
			for (int z = block.min[2]; z <= block.max[2]; z++)
				grid.u_vec[id][index_converter(block.max[0], y, z, block)] = 
					u_analytical(block.max[0] * grid.h_x, y * grid.h_y, z * grid.h_z, grid.ts, grid);

	// Axis Y. Periodic condition

	if (block.min[1] == 0)
		#pragma omp parallel for
		for (int x = block.min[0]; x <= block.max[0]; x++)
			for (int z = block.min[2]; z <= block.max[2]; z++)
				grid.u_vec[id][index_converter(x, block.min[1], z, block)] = 
					u_analytical(x * grid.h_x, block.min[1] * grid.h_y, z * grid.h_z, grid.ts, grid);

	if (block.max[1] == grid.N)
		#pragma omp parallel for
		for (int x = block.min[0]; x <= block.max[0]; x++)
			for (int z = block.min[2]; z <= block.max[2]; z++)
				grid.u_vec[id][index_converter(x, block.max[1], z, block)] =
					u_analytical(x * grid.h_x, block.max[1] * grid.h_y, z * grid.h_z, grid.ts, grid);

	// Axis Z. Condition of the first kind

	if (block.min[2] == 0)
		#pragma omp parallel for
		for (int x = block.min[0]; x <= block.max[0]; x++)
			for (int y = block.min[1]; y <= block.max[1]; y++)
				grid.u_vec[id][index_converter(x, y, block.min[2], block)] = 0;

	if (block.max[2] == grid.N)
		#pragma omp parallel for
		for (int x = block.min[0]; x <= block.max[0]; x++)
			for (int y = block.min[1]; y <= block.max[1]; y++)
				grid.u_vec[id][index_converter(x, y, block.max[2], block)] = 0;
}

int main(int argc, char **argv)
{
	// parse two parameters
	int N = atoi(argv[1]);
	double L = strcmp(argv[2], "1") ? 1 : M_PI;

	Grid omega;

	omega.L_x = L;
	omega.L_y = L;
	omega.L_z = L;
	omega.N = N;
	omega.T = CONST_T;
	omega.K = CONST_K;
	omega.h_x = omega.L_x / omega.N;
	omega.h_y = omega.L_y / omega.N;
	omega.h_z = omega.L_z / omega.N;
	omega.ts = omega.T / omega.K;

	// standart variables
	int rank;
	int size;
	double time;
	double t; // result time

	double inaccur = 0;

	// parallelization
	if(MPI_Init(&argc, &argv) != MPI_SUCCESS)
	{
		printf("MPI_Init failed\n");
		return 1;
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	time = MPI_Wtime();

	// 1 - split
	block_counter = 0;
	Block blocks[size];
	split(0, omega.N, 0, omega.N, 0, omega.N, AXIS_X, size, blocks);
	Block block = blocks[rank];

	for (int i = 0; i < 3; i++)
        omega.u_vec[i].resize(block.size);

	for (int i = 0; i < size; i++)
	{
		if (i == rank)
			continue;

		// 2 - search for neighbors and forwarding
		Block local_block = blocks[i];
		if ((block.min[0] == local_block.max[0] + 1) || (local_block.min[0] == block.max[0] + 1))
		{
			int x_1 = block.min[0];
			if (x_1 != local_block.max[0] + 1)
				x_1 = block.max[0];
			int x_2 = local_block.min[0];
			if (x_2 != block.max[0] + 1)
				x_2 = local_block.max[0];

			if (nested(block, local_block, 1, 2))
			{
				omega.send.push_back(newIdBlock(i, x_1, x_1, block.min[1], block.max[1], block.min[2], block.max[2]));
				omega.receive.push_back(newIdBlock(i, x_2, x_2, block.min[1], block.max[1], block.min[2], block.max[2]));
			}
			else if (nested(local_block, block, 1, 2))
			{
				omega.send.push_back(newIdBlock(i, x_1, x_1, local_block.min[1], local_block.max[1], local_block.min[2], local_block.max[2]));
				omega.receive.push_back(newIdBlock(i, x_2, x_2, local_block.min[1], local_block.max[1], local_block.min[2], local_block.max[2]));
			}
		}

		if ((block.min[1] == local_block.max[1] + 1) || (local_block.min[1] == block.max[1] + 1))
		{
			int y_1 = block.min[1];
			if (y_1 == local_block.max[1] + 1)
				y_1 = block.max[1];
			int y_2 = local_block.min[1];
			if (y_2 == block.max[1] + 1)
				y_2 = local_block.max[1];

			if (nested(block, local_block, 0, 2))
			{
				omega.send.push_back(newIdBlock(i, block.min[0], block.max[0], y_1, y_1, block.min[2], block.max[2]));
				omega.receive.push_back(newIdBlock(i, block.min[0], block.max[0], y_2, y_2, block.min[2], block.max[2]));
			}
			else if (nested(local_block, block, 0, 2))
			{
				omega.send.push_back(newIdBlock(i, local_block.min[0], local_block.max[0], y_1, y_1, local_block.min[2], local_block.max[2]));
				omega.receive.push_back(newIdBlock(i, local_block.min[0], local_block.max[0], y_2, y_2, local_block.min[2], local_block.max[2]));
			}
		}

		if (block.min[2] == local_block.max[2] + 1 or local_block.min[2] == block.max[2] + 1)
		{
			int z_1 = block.min[2] == local_block.max[2] + 1 ? block.min[2] : block.max[2];
			int z_2 = local_block.min[2] == block.max[2] + 1 ? local_block.min[2] : local_block.max[2];

			if (nested(block, local_block, 0, 1))
			{
				omega.send.push_back(newIdBlock(i, block.min[0], block.max[0], block.min[1], block.max[1], z_1, z_1));
				omega.receive.push_back(newIdBlock(i, block.min[0], block.max[0], block.min[1], block.max[1], z_2, z_2));
			}
			else if (nested(local_block, block, 0, 1))
			{
				(omega.send).push_back(newIdBlock(i, local_block.min[0], local_block.max[0], local_block.min[1], local_block.max[1], z_1, z_1));
				(omega.receive).push_back(newIdBlock(i, local_block.min[0], local_block.max[0], local_block.min[1], local_block.max[1], z_2, z_2));
			}
		}
	}

	border(0, block, omega);
	border(1, block, omega);
cout << "Start pragma\n";
	// 3 - search solutions 
	#pragma omp parallel for
	for (int x = fmax(1, block.min[0]); x <= fmin(N - 1, block.max[0]); x++)
		for (int y = fmax(1, block.min[1]); y <= fmin(N - 1, block.max[1]); y++)
			for (int z = fmax(1, block.min[2]); z <= fmin(N - 1, block.max[2]); z++)
				omega.u_vec[0][index_converter(x, y, z, block)] = u_analytical(x * omega.h_x, y * omega.h_y, z * omega.h_z, 0, omega);

	vector<Data> rec = messenger(0, block, omega);

	#pragma omp parallel for
	for (int x = fmax(1, block.min[0]); x <= fmin(N - 1, block.max[0]); x++)
		for (int y = fmax(1, block.min[1]); y <= fmin(N - 1, block.max[1]); y++)
			for (int z = fmax(1, block.min[2]); z <= fmin(N - 1, block.max[2]); z++)
				omega.u_vec[1][index_converter(x, y, z, block)] = omega.u_vec[0][index_converter(x, y, z, block)] +
					omega.ts * omega.ts / 2 * analogue_laplace(0, x, y, z, rec, block, omega);

	for (int i = 2; i <= STEPS; i++)
	{
		rec = messenger((i + 2) % 3, block, omega);

		#pragma omp parallel for
		for (int x = fmax(1, block.min[0]); x <= fmin(N - 1, block.max[0]); x++)
			for (int y = fmax(1, block.min[1]); y <= fmin(N - 1, block.max[1]); y++)
				for (int z = fmax(1, block.min[2]); z <= fmin(N - 1, block.max[2]); z++)
					omega.u_vec[i % 3][index_converter(x, y, z, block)] = 2 * omega.u_vec[(i + 2) % 3][index_converter(x, y, z, block)] -
						omega.u_vec[(i + 1) % 3][index_converter(x, y, z, block)] +
						omega.ts * omega.ts * analogue_laplace((i + 2) % 3, x, y, z, rec, block, omega);

		border(i % 3, block, omega);
	}
cout << "4 part\n";
	// 4 - count inaccur
	double local = 0;

	#pragma omp parallel for reduction(max: local)
	for (int x = block.min[0]; x <= block.max[0]; x++)
		for (int y = block.min[1]; y <= block.max[1]; y++)
			for (int z = block.min[2]; z <= block.max[2]; z++)
					local = fmax(fabs(omega.u_vec[STEPS % 3][index_converter(x, y, z, block)] -
							u_analytical(x * omega.h_x, y * omega.h_y, z * omega.h_z, STEPS * omega.ts, omega)),
							local);
cout << "Finalize\n";
	MPI_Reduce(&local, &inaccur, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	time = MPI_Wtime() - time;
	MPI_Reduce(&time, &t, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if(rank == 0)
	{
		printf("N: %d\n", N);
		printf("inaccur: %.10f\n", inaccur);
		printf("size: %d\n", size);
		printf("t: %.10fs\n", t);
	}

	if(MPI_Finalize() != MPI_SUCCESS)
	{
		printf("MPI_Finalize failed\n");
		return 1;
	}

	return 0;
}
