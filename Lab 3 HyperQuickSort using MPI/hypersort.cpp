// mpic++ sort.cpp -o hypersort -std=c++11
// mpirun -np 8 hypersort inputxx.txt

#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <utility>
#include <stdio.h>
#include <mpi.h>
#include <sstream>
#include <cmath>


#define fr(i,n) for(int i=0; i<n; i++)
#define frr(i,a,b) for(int i=a; i<=b; i++)
#define INF (1<<30)

using namespace std;

int fun(const void *a, const void *b){
	return *(int*)a >= *(int*)b ? 1 : -1;
}


int main(int argc, char *argv[]){

	int *a; 	// Input Array
	int N; 		// Size of input array
	int P;		// #Processes
	int rank;		// My Rank

	int *ar; 	// my array
	int *l;		// left array
	int *r;		// right array
	int *rcv; 	// receive array

	int arx; 	// my array size
	int lx;		// left size
	int rx;		// right size
	int rcvx;	// receive size

	int pivot;	// median	-> to receive
	int pi;		// pivot index

	int ng; 	// #groups
	int gs; 	// group size


	////////////////////////////////////// Taking input from input file //////////////////////////////////////
	string ifile = argv[1];
	ifstream infile(ifile.c_str());
	infile >> N;
	int powersize = log2(N);
	stringstream ss;
	ss << powersize;
	string outname = "output_2power" + ss.str() + ".txt";
	
	a = (int *)malloc(N*sizeof(int));
	
	fr(i,N){
		infile >> a[i];
	}

	
    
	// utilities
	int i;
	bool islower;
	double starttime, endtime;

	MPI_Status status;
	//////////////////////////////////////////////  MPI  /////////////////////////////////////////////

	MPI_Init(&argc, &argv); 

	MPI_Comm_size(MPI_COMM_WORLD, &P); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	

	// Preparing the local arrays for each node
	arx = N/P;
	ar = (int *)malloc(arx*sizeof(int));

	fr(k,arx){
		ar[k] = a[ rank*arx + k ];
	}
	
	ng = 1;
	gs = P/ng;

	starttime = MPI_Wtime();

	// The log(n) iterations begin !!
	while(true){

		qsort(ar,arx,sizeof(int),fun);

		if(rank%gs < gs/2) islower = true; // check this ..
		else islower = false;
		i = rank/gs;	// The group number (0,1,2 ... ,ng-1)
		// The nodes in ith group are from i*gs to (i+1)*gs-1


		// break if group size is 1
		if(gs == 1) break;

		
		// Sending the pivot to its other group memebers
		// All others receiving the pivot and calculating the index of split
		if(rank%gs == 0){
			pi = (arx-1)/2;
			pivot = ar[pi];
			for(int j=i*gs+1; j<(i+1)*gs; j++){
				MPI_Send(&pivot, 1, MPI_INT, j, 1 /*Tag*/, MPI_COMM_WORLD);
			}
		}
		else{
			MPI_Recv(&pivot, 1, MPI_INT, i*gs, 1, MPI_COMM_WORLD, &status);
			
			pi = -1;
			for(int j=0; j<arx; j++){
				if(ar[j]>pivot){
					pi = j-1;
					if(j==0) pi = 0;
					break;
				}
			}
			if(pi == -1){
				pi = arx-1;
			}
		}

		// Calculate the sizes of LOWER and UPPER arrays
		// And prepare those arrays : r and s, to send and receive

		// for the lower half nodes in a group the right half is going to be exchanged with the 
		// left half of the upper half of the nodes of this group.
		// the lower node will send the size of its right half to upper node, and get upper's left size

		lx = pi+1;
		rx = arx-lx;

		if(arx==0){
			lx = rx = 0;
		}


		l = (int *)malloc(lx*sizeof(int));
		r = (int *)malloc(rx*sizeof(int));


		// Preparing the left and right arrays
		for(int j=0;j<pi+1;j++){
			l[j] = ar[j];
		}
		for(int j=pi+1;j<arx;j++){
			r[j-pi-1] = ar[j]; 
		}

		
		// First send & recv the size of the array its going to receive
		if(islower){
			MPI_Sendrecv(&rx, 1, MPI_INT, rank+gs/2, 1, &rcvx, 1, MPI_INT, rank+gs/2, 1, MPI_COMM_WORLD, &status);
			arx = lx + rcvx;
		}
		else{
			MPI_Sendrecv(&lx, 1, MPI_INT, rank-gs/2, 1, &rcvx, 1, MPI_INT, rank-gs/2, 1, MPI_COMM_WORLD, &status);
			arx = rx + rcvx;
		}

		// allocating size for the receive array
		rcv = (int *)malloc(rcvx*sizeof(int));

		
		// Time to exchange the split parts
		if(islower){
			MPI_Sendrecv(r, rx, MPI_INT, rank+gs/2, 1, rcv, rcvx, MPI_INT, rank+gs/2, 1, MPI_COMM_WORLD, &status);
		}
		else{
			MPI_Sendrecv(l, lx, MPI_INT, rank-gs/2, 1, rcv, rcvx, MPI_INT, rank-gs/2, 1, MPI_COMM_WORLD, &status);
		}

		
		free(ar);
		ar = (int *)malloc(arx*sizeof(int));

		if(islower){
			for(int j=0;j<lx;j++){
				ar[j] = l[j]; 
			}
			for(int j=lx;j<arx;j++){
				ar[j] = rcv[j-lx];
			}
		}
		else{
			for(int j=0;j<rcvx;j++){
				ar[j] = rcv[j]; 
			}
			for(int j=rcvx;j<arx;j++){
				ar[j] = r[j-rcvx];
			}
		}

		free(l);
		free(r);
		free(rcv);

		ng = ng*2;
		gs = P/ng;
	}
	////////////////////////////////////// Printing to the output file //////////////////////////////////////
	
	// To output each process first prints its share of array in the output file
	// adn then sends a signal to the next process, which once receives that signal 
	// prints its own share. This goes on till the last process
	// So effectively they all print sequentially in the output file

	// for the case of just one process
	if(P==1){
		ofstream outfile;
  		outfile.open(outname.c_str(), std::ofstream::app);
		
		fr(j,arx){
			outfile << ar[j] << " ";
		}
	}else{
		if(rank == 0){
			ofstream outfile;
	  		outfile.open(outname.c_str(), std::ofstream::app);
			
			fr(j,arx){
				outfile << ar[j] << " ";
			}
			int x = 0;
			MPI_Send(&x, 1, MPI_INT, rank+1, 1 /*Tag*/, MPI_COMM_WORLD);
		}
		else if(rank == P-1){
			int x = 0;
			MPI_Recv(&x, 1, MPI_INT, rank-1, 1, MPI_COMM_WORLD, &status);
			ofstream outfile;
	  		outfile.open(outname.c_str(), std::ofstream::app);
			
			fr(j,arx){
				outfile << ar[j] << " ";
			}
		}
		else{
			int x = 0;
			MPI_Recv(&x, 1, MPI_INT, rank-1, 1, MPI_COMM_WORLD, &status);
			ofstream outfile;
	  		outfile.open(outname.c_str(), std::ofstream::app);
			
			fr(j,arx){
				outfile << ar[j] << " ";
			}	
			MPI_Send(&x, 1, MPI_INT, rank+1, 1 /*Tag*/, MPI_COMM_WORLD);
		}
	}


	MPI_Finalize();
	


	endtime   = MPI_Wtime();
    printf("That took %f seconds\n",endtime-starttime);
}