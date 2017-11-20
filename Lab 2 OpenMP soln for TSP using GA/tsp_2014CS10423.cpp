#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <utility>
#include <omp.h>
#include <stdio.h>


#define fr(i,n) for(int i=0; i<n; i++)
#define frr(i,a,b) for(int i=a; i<=b; i++)
#define INF (1<<30)

using namespace std;


int NUMT = 0;  			// #threads
int SIZE = 0;			// #cities
int N_FITTEST = 400;	// #fittest to choose
int N_ITER = 100;		// #iterations
int POP_SIZE = 2000;	// Poplulation Size
double lengthx;			// Best path length

// Vectors
vector<int> ans; 					// Best Path
vector<vector<double>> coords;		// Coordinates of cities
vector<vector<double>> weights;		// Distance bw cities
vector<vector<int>> pop;			// Population
vector<pair<double,int>> fval;		// Fitness vals of all population, 2b sorted

bool fun(pair<double,int> p1, pair<double,int> p2){
	return p1.first < p2.first;
}


void populate(){
	vector<int> a;
	a.resize(SIZE);
	fr(i,SIZE){
		a[i] = i;
	}
	fr(i,POP_SIZE+N_FITTEST){
		random_shuffle(a.begin(), a.end());
		fr(j,SIZE){
			pop[i][j] = a[j];
		}
	}
	
}

void pmx(vector<int> &ch, int p1, int p2){

	vector<int> mapp;
	mapp.resize(SIZE,-1);
	ch.resize(SIZE,-1);

	int i1=0, i2=0, x;
	// generating the two random indexes
	while(i1==i2){
		i1 = rand()%SIZE;
		i2 = rand()%SIZE;
	}
	if(i1 > i2) swap(i1,i2);

	// making the mapping
	frr(i,i1,i2){
		mapp[pop[p2][i]] = pop[p1][i]; 
		ch[i] = pop[p2][i];
	}
	
	bool b = false;
	fr(i,SIZE){
		// transitivity
		x = pop[p1][i];
		if(i>=i1 && !b){
			b = true;
			i = i2;
			continue;
		}
		while(true){
			if(mapp[x] == -1) break;
			x = mapp[x];
		}
		ch[i] = x;

	}

	// Check child ch[] here 
}

void mutate(int p){
	int i1=0, i2=0, x;
	while(i1==i2){
		i1 = rand()%SIZE;
		i2 = rand()%SIZE;
	}
	swap(pop[p][i1], pop[p][i2]);
}


void gx(vector<int> &ch, int p1, int p2){
	vector<bool> b;
	b.resize(SIZE,false);
	ch.resize(SIZE);
	ch[0] = pop[p1][0];
	b[ch[0]] = true;
	int a,l1,l2;

	frr(i,1,SIZE-1){
		a = ch[i-1]; 						// a is the current city
		fr(j,SIZE){
			if(pop[p1][j] == a){
				l1 = pop[p1][(j+1)%SIZE]; 	// City next to a in parent1 
				break;
			}
		}
		fr(j,SIZE){
			if(pop[p2][j] == a){
				l2 = pop[p2][(j+1)%SIZE]; 	// City next to a in parent2 
				break;
			}
		}
		// Now we will find if l1 and l2 are picked already or not
		// and based on that we will select our next city to be added in child

		if(!b[l1] && b[l2]){		// l2 already taken
			ch[i] = l1;
			b[l1] = true;
		}
		else if(!b[l2] && b[l1]){	// l1 already taken
			ch[i] = l2;
			b[l2] = true;
		}
		else if(!b[l2] && !b[l1]){ 	// none taken
			if(weights[a][l1] < weights[a][l2]){
				ch[i] = l1;
				b[l1] = true;
			}
			else{
				ch[i] = l2;
				b[l2] = true;
			}
		}
		else if(b[l2] && b[l1]){	// both taken
			l1 = l2 = -1;
			fr(j,SIZE){
				if(!b[j]){
					if((rand()%20)<8){
						l1 = j;
					}
					l2 = j;
				}
			}
			if(l1 == -1) l1=l2;
			ch[i] = l1;
			b[l1] = true;
		}
	}
}


void fittest(int a[]){
	double f;

	fr(i,POP_SIZE+N_FITTEST){
		// f = 0;
		f = weights[pop[i][SIZE-1]][pop[i][0]];
		fr(j,SIZE-1){
			f += weights[pop[i][j]][pop[i][j+1]];
		}
		fval[i] = make_pair(f,i);
	}



	sort(fval.begin(),fval.end(),fun);

	fr(i,N_FITTEST){
		a[i] = fval[i].second;
	}

	if(fval[0].first < lengthx){
		lengthx = fval[0].first;
		fr(i,SIZE){
			ans[i] = pop[fval[0].second][i];
		}
	}
	
}



int main(int argc, char const *argv[]){

	double start_time = omp_get_wtime();

	// For random number generation
	time_t seconds;
	time(&seconds);
	srand((unsigned int) seconds);

	// Taking input from input file
	string x, ofile, ifile = argv[1];
	NUMT = atoi(argv[2]);
	ifstream infile(ifile);
	infile >> x >> ofile >> x >> SIZE >> x;
	cout<< "Size is " <<SIZE <<endl;
	ans.reserve(SIZE);
	coords.reserve(SIZE);
	weights.reserve(SIZE);
	fr(i,SIZE){coords[i].reserve(2);}
	fr(i,SIZE){weights[i].reserve(SIZE);}

	lengthx = INF;

	// cout<<"1"<<endl;
	fr(i,SIZE){
		infile >> x >> coords[i][0] >> coords[i][1];
	}
	infile.close();
	

	// Can optimize it. Its symmetric
	fr(i,SIZE){
		fr(j,SIZE){
			weights[i][j] = sqrt(pow(coords[i][0]-coords[j][0],2) + pow(coords[i][1]-coords[j][1],2));
			// cout << weights[i][j] << " ";
		}
		// cout << endl;
	}



	pop.resize(POP_SIZE+N_FITTEST);
	fr(i,POP_SIZE+N_FITTEST){pop[i].reserve(SIZE);}
	fval.resize(POP_SIZE+N_FITTEST);

	populate();

	int fitones[N_FITTEST];
	
	vector<vector<int>> popx;
	popx.resize(N_FITTEST);
	fr(i,N_FITTEST){
		popx[i].resize(SIZE);
	}

	
	

// FOR LOOP >>>>>>>>>>>>>>>>>>

	for(int k = 0; k < N_ITER; k++){

		fittest(fitones);


		#pragma omp parallel for num_threads(NUMT)

			for(int i=0; i<N_FITTEST; i++){
				int i1, i2;
				i1 = i2 = -1;
				while(i1==i2){
					i1 = rand()%N_FITTEST;
					i2 = rand()%N_FITTEST;
				}

				vector<int> child;
				// No shared memory till here


				if(rand()%10<8){
					gx(child, fitones[i1], fitones[i2]); // Reading pop
				}
				else{
					pmx(child, fitones[i1], fitones[i2]); // Reading pop
				}

				if(rand()%100 < 10){ 	// Reading and writing pop
					mutate(i1);
				}
				if(rand()%100 < 10){	// Reading and writing pop
					mutate(i2);
				}

				fr(j,SIZE){
					popx[i][j] = child[j]; 
					// pop[POP_SIZE+i][j] = child[j];		// Writing in last 100 elements of pop
				}
				
			}

			fr(i,N_FITTEST){
				fr(j,SIZE){
					pop[POP_SIZE+i][j] = popx[i][j];
				}
			}
		 // parallel region ends
	}

	cout << "The best path is "<< endl;
	fr(i,SIZE){
		cout<< ans[i] + 1 <<" ";
	}
	cout<< endl << endl << " The best path length is " << endl << lengthx << endl;



	double time = omp_get_wtime() - start_time;

	cout << "Time taken "<<time <<endl;


	ofstream outfile("output_"+ofile+".txt");
    outfile<< "DIMENSION : " << SIZE << endl << "TOUR_LENGTH : " << lengthx << endl << "TOUR_SECTION" << endl ;
   	
   	fr(i,SIZE){
   		outfile << ans[i] << " " ;
   	}	

   	outfile << endl << -1 << endl << "EOF" ;





	////////////////////////////////////////////////////////////////////
	// cout << "---------------------------------------------" << endl;

	// int i1=-1,i2=-1;
	// while(i1==i2){
	// 	i1 = rand()%N_FITTEST;
	// 	i2 = rand()%N_FITTEST;
	// }

	// vector<int> child;
	// pmx(child,i1,i2);

	// fr(i,SIZE){
	// 	cout<<pop[i1][i]<<" ";
	// }
	// cout<<endl;
	// fr(i,SIZE){
	// 	cout<<pop[i2][i]<<" ";
	// }
	// cout<<endl;
	// fr(i,SIZE){
	// 	cout<<child[i]<<" ";
	// }
	// cout<<endl;
	////////////////////////////////////////////////////////////////////








	return 0;
}