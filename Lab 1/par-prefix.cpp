#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <pthread.h>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <openssl/md5.h>

using namespace std;

#define MAX_SIZE 1024*1024*16
// #define NUMT 4

unsigned long *ar;
int n = 0;
int NUMT = 4;
// long long counter[NUMT];

// int rj[NUMT];


void *presum(void * t){
    int *id = (int *)t; 
    // cout<<" my id is "<<*id<<endl;
    int left, right;
    if(*id==NUMT-1){
        left = (*id)*(n/NUMT); 
        right = n;
    }
    else{
        left = (*id)*(n/NUMT); 
        right = left + (n/NUMT);    
    }
     
    for(int i=left+1; i < right; ++i){
        ar[i] += ar[i-1];  
    }
    pthread_exit(0);
}

void *finish(void *t){
    int *id = (int *)t; 
    // cout<<" my id is "<<*id<<endl;
    int left, right;
    if(*id==NUMT-1){
        left = (*id)*(n/NUMT); 
        right = n;
    }
    else{
        left = (*id)*(n/NUMT); 
        right = left + (n/NUMT);    
    }
    unsigned long val = ar[left-1];
    for(int i=left; i < right-1; ++i){
        ar[i] += val;  
    }
    pthread_exit(0);
}


int main(int argc, char *argv[]) {
	
	FILE *fptr, *outfile;

    NUMT = atoi(argv[1]);
    ifstream infile("input.txt");
//    ofstream outfile("output.txt");

	outfile = fopen("output.txt", "w");

    // ar = (int*)malloc(sizeof(int)*MAX_SIZE);    // Try by vector if it becomes fast
                                                
                                                // realloc ()  ....

                                                // long long ...

                                                // OpenMP time

    unsigned long k;
    string s;
    
    getline(infile,s);
    
    n = stoi(s);

    ar = (unsigned long*)malloc(sizeof(unsigned long)*n);

    getline(infile,s);

    int tmp = 0;

    for (int i=0; i<s.size(); ++i) {
        k=0; 
        while(i<s.size() && s[i]!=' '){
            k = k*10 + stoi(s.substr(i,1));
            i++;
        }
        ar[tmp++] = k;
    }


    // for(int i=0; i<n; ++i){
    //     cout<<ar[i]<<" ";
    // }
    // cout<<endl<<endl;
    // cout<<"n is "<<n<<endl<<endl;

    pthread_t threads[NUMT];
    int rj[NUMT];

    double start = omp_get_wtime();

    for(int i=0; i<NUMT; ++i){
        rj[i] = i ;
        auto xx = pthread_create(&threads[i], NULL, presum, &rj[i]);
        if (xx){
            printf("Something is Wrong !! %d\n", xx);
            exit(-1);
        }
    }
    for (int i=0; i<NUMT; i++){
        pthread_join(threads[i], NULL);
    }

    // cout<<"joinrd"<<endl;

    // for(int i=0; i<n; ++i){
    //     cout<<ar[i]<<" ";
    // }
    // cout<<endl<<endl;

    int l,r;

    for(int i=1;i<NUMT;i++){

        l = i*(n/NUMT);
        r = (i+1)*(n/NUMT);
        if(i==NUMT-1){
            r = n;   
        }
        
        // cout<<l<<" "<<r<<endl;
        ar[r-1] += ar[l-1];
        // cout<<"Hello"<<endl;
    }

    // for(int i=0; i<n; ++i){
    //     cout<<ar[i]<<" ";
    // }
    // cout<<endl<<endl;

    // cout<<"loop ended"<<endl;

    for(int i=1; i<NUMT; ++i){
        rj[i] = i ;
        
        auto xx = pthread_create(&threads[i], NULL, finish, rj + i);
        if(xx){
            printf("Something is Wrong !! %d\n", xx);
            exit(-1);
        }
    }

    for (int i=1; i<NUMT; i++){
        pthread_join(threads[i], NULL);
    }


    // for(int i=0; i<n; ++i){
    //     cout<<ar[i]<<" ";
    // }
    // cout<<endl<<endl;

    double time = omp_get_wtime() - start;



   // FILE *fptr, *outfile;
    int i;
    fptr = fopen("temp.txt", "w");
    // int i;
    for(i=0; i<n-1; i++)
        fprintf(fptr,"%lu ",ar[i]);
    fprintf(fptr,"%lu",ar[i]);
    fclose(fptr);
    fptr = fopen("temp.txt", "rb");
    
    MD5_CTX mdContext;
    unsigned char c[MD5_DIGEST_LENGTH];
    int bytes;
    
    unsigned char dd[1024];

    MD5_Init (&mdContext);
    while((bytes = fread(dd, 1, 1024, fptr)) != 0)
        MD5_Update (&mdContext, dd, bytes);
    MD5_Final (c,&mdContext);
    // int i;
    // for(i = 0; i < MD5_DIGEST_LENGTH; i++)
    //     printf("%02x", c[i]);
    fclose(fptr);



	fprintf(outfile, "THREADS : %d\n", NUMT );
 	fprintf(outfile, "Time : %f\n", time);
	fprintf(outfile, "Md5-sum : ");








//    cout<<" OMP Time : "<< time <<endl;

  //  outfile<<"Threads : "<<NUMT<<endl;
    //outfile<<"Time : "<<time<<endl;
    //outfile<<"Md5-sum : ";//<<endl;

    for(i = 0; i < MD5_DIGEST_LENGTH; i++)
        fprintf(outfile, "%02x", c[i]);

    //outfile<<endl;

    //outfile.close();
	fclose(outfile);
    return 0;
}
