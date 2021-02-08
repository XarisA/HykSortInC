#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi/mpi.h>
#include <omp.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "functions.h"

#define SIZE 1000000

void ParallelSelect(int *arr,int * spl,int n, int kway,int p,int rank, MPI_Comm comm);
int Hyksort(int *Ar,int *Spl,int *GlB,int n, int k,int p,int rank, MPI_Comm comm);

int SendRecieve_kway(int *Buckets,int *B,int n,int p,int rank,int k,MPI_Comm comm);

int * SampleSplitters(int *Ar,int n, int p, int rank) ;
void RankSampleSplitters(int *Ar,int *Spl,int n, int p, int rank,MPI_Comm comm);

void scan(int *A, int *B,int cnt, int p);
void bitonic_sort(int *Ar,int n,int p, int rank,MPI_Comm comm);

int main (int argc, char *argv[])
{
    int n,p,rank,i;
    int tag1,tag2,tag3;
    long N=SIZE;
    int A[SIZE];
  	int GlB[SIZE];
    int *Ar=(int *)malloc(sizeof(int)*N);

    double start,end;

    tag2=200;
    tag3=300;


    MPI_Comm comm;
    MPI_Status status;

    comm = MPI_COMM_WORLD;

    MPI_Init(NULL,NULL);
    MPI_Comm_size(comm,&p);
    MPI_Comm_rank(comm,&rank);

    int k=p;
    int *Spl=(int *)malloc(sizeof(int)*(k-1)); //Sorted Global Splitters

    n=N/p;
        MPI_Barrier(comm);
        // Table Creation
        if(rank==0) {

            //Create my Array
            init_array(A, N, 1, 100000);
            //Comment this ^ and uncomment the following for validating excel values in testCase1
            //int A[64]={61,25,63,59,7,11,21,3,36,18,50,66,93,86,22,40,85,94,14,8,21,45,7,42,57,72,28,56,9,87,72,69,63,87,28,69,49,100,23,36,17,73,1,9,58,74,1,42,68,14,50,40,11,8,82,67,79,61,75,40,100,98,60,62};
            printf("HykSort Implementation\n\n");
            // Uncomment to write to stdout
//            printf("Generating array with random int A: ");
//            print_array(A, N);
        }
        MPI_Scatter(A, n, MPI_INT, Ar, n, MPI_INT, 0, MPI_COMM_WORLD);

    if (p>1){
        // Uncomment to write to stdout
//        print_array_in_process(Ar,n,p,rank,"Split A to Ar");
        MPI_Barrier(comm);
        // Uncomment to write to stdout
            start = MPI_Wtime();
            Hyksort(Ar,Spl,GlB,n, k,p,rank,comm);
            end = MPI_Wtime();
        if (rank==0){
//            print_array_in_process(GlB, N, p, rank, "Buckets gathered into sorted input");
            printf("Elapsed time is %f\n", end - start);
        }
    }
    else if(p==1){
        printf("Sorted input array using quicksort algorithm: ");
        start = MPI_Wtime();
        local_sort(Ar,N);
        end = MPI_Wtime();
        // Uncomment to write to stdout
//        print_array(Ar, N);
        printf( "Elapsed time is %f\n", end - start );
    }

    // Uncomment the following to use BitonicSort . Also comment out Hyksort
//    if(p>1){
//        start = MPI_Wtime();
//        bitonic_sort(Ar,n,p,rank, MPI_COMM_WORLD);
//        MPI_Allgather(&Ar[0],n,MPI_INT,&GlB[0],n,MPI_INT,MPI_COMM_WORLD);
//        end = MPI_Wtime();
//        if (rank == 0){
//            printf("Sorted input array using bitonic sort algorithm: ");
//            print_array(GlB, N);
//            printf( "Elapsed time is %f\n", end - start );
//        }
//    }

    MPI_Finalize();
    return 0;
}


int Hyksort(int *Ar,int *Spl,int *GlB,int n, int k,int p,int rank, MPI_Comm comm)
{
    int N,i,j,v,w,y;
    N=n*p;
    int *B=(int *)malloc(sizeof(int)*N);

    // Local Sort & Sample Distribution
    // Each processor sorts localy
    local_sort(Ar,n);
    MPI_Barrier(comm);
    // Uncomment to write to stdout
    //print_array_in_process(Ar,n,p,rank,"Sort local Ar");

    //Choose between ParallelSelect splitters and SampleSplitters
    ParallelSelect(Ar,Spl,n,k,p,rank,comm);
    //RankSampleSplitters(Ar,Spl,n,p,rank,comm);

    // Uncomment to write to stdout
//        if(rank==0){
//            print_array_in_process(Spl, k-1, p, rank, "Sorted splitters");
//        }


    //**** Creating Buckets locally 
    //Rank(Ar,Spl)

    int *Buckets = (int *) malloc (sizeof (int) * (N + p));
    v = 0;
    w = 1;

    for (y=0; y<n; y++){
        if(v < (p-1)){
            if (Ar[y] < Spl[v])
                Buckets[((n + 1) * v) + w++] = Ar[y];
            else{
                Buckets[(n + 1) * v] = w-1;
                w=1;
                v++;
                y--;
            }
        }
        else
            Buckets[((n + 1) * v) + w++] = Ar[y];
    }
    Buckets[(n + 1) * v] = w - 1;
    MPI_Barrier(comm);
    // Uncomment to write to stdout
//    print_array_in_process(Buckets, N+p, p, rank, "Buckets");

// End of Rank Function

    int *BucketBuffer = (int *) malloc (sizeof (int) * (N + p));

    // All to All broadcast
    // Starting with MPI_Alltoall.
    MPI_Alltoall (Buckets, n + 1, MPI_INT, BucketBuffer,
                 n + 1, MPI_INT, MPI_COMM_WORLD);
   
    //Choose between ^^ MPI_Alltoall and  SendRecieve_kway
    /*
     * My implementation of AllToAll_Kway (Algorithm 3.5) found in
     * paper "Hyksort: A new Variant of Hypercube Quicksort on Distributed Memory Architectures
     */
    // while (p>1){
    //     SendRecieve_kway(Buckets,BucketBuffer, n,p,rank,4,comm);
    // }

    MPI_Barrier(comm);

    //Remove displacements
    //
    j=0;
    int count=0;
    int st;
    while(j<=k){
        st=j*(n+1);
        for(i=st+1;i<=(BucketBuffer[st]+st);i++){
            B[count]=BucketBuffer[i];
            count++;
        }
        j++;
    }

    int Bsize =0;
    for(i=0;i<N;i++){
        if(B[i]!=0){
            Bsize++;
        }
    }
    MPI_Barrier(comm);
    local_sort(B,Bsize);
    // Uncomment to write to stdout
//    print_array_in_process(BucketBuffer, N, p, rank, "Buckets redistributed");
//    print_array_in_process(B, N, p, rank, "Buckets sorted");


    int * counts = (int *)malloc(p*sizeof(int));
    int * displs = (int *)malloc(p*sizeof(int));

    MPI_Gather(&Bsize, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(comm);
    //printf("Bsize %d ",Bsize);
    //print_array_in_process(counts,p,p,rank,"recvcount ");
    if(rank==0){
        for ( i = 0; i < p; i++){
            if (i>0){
                displs[i]= displs[i-1] + counts[i-1];
            }
            else {displs[i]=0;}
        }
        // Uncomment to write to stdout
//        print_array_in_process(displs,p,p,rank,"Displacements ");
//        print_array_in_process(counts,p,p,rank,"Recievecounts ");
    }
    MPI_Barrier(comm);

    // Gather all elements in rank 0.
    MPI_Gatherv(B,Bsize,MPI_INT,GlB, counts,displs,MPI_INT,0,comm);
    MPI_Barrier(comm);

    if (rank == 0){
        //sleep(2);
        local_sort(GlB,N);
    }

    free(B);
    free(counts);
    free(displs);
    free(Buckets);
   // free(BucketBuffer);

    MPI_Barrier(comm);
}

void ParallelSelect(int *arr,int * spl, int n ,int kway, int p,int rank,MPI_Comm comm) {
    /*
     * My implementation of ParallelSelect (Algorithm 3.6) found in
     * paper "Hyksort: A new Variant of Hypercube Quicksort on Distributed Memory Architectures
     *
     */

    int npes=p;
    int i;

//    MPI_Comm_size(comm, &npes);
//    MPI_Comm_rank(comm, &rank);

    int totSize, nelem = n;
    MPI_Allreduce(&nelem, &totSize, 1,MPI_INT, MPI_SUM, comm);

    //Determine splitters. O( log(N/p) + log(p) )
    int splt_count = (1000*kway*nelem)/totSize;
    if (npes>1000*kway) splt_count = (((float)rand()/(float)RAND_MAX)*totSize<(1000*kway*nelem)?1:0);
    if (splt_count>nelem) splt_count=nelem;
    int splitters[splt_count];
    for(size_t i=0;i<splt_count;i++)
        splitters[i] = arr[rand()%nelem];

    // Gather all splitters. O( log(p) )
    int glb_splt_count;
    int glb_splt_cnts[npes];
    int glb_splt_disp[npes];
    for (i=0;i<npes;i++){
        glb_splt_disp[i]=0;
    }
    MPI_Allgather(&splt_count,1,MPI_INT, &glb_splt_cnts[0],1,MPI_INT, comm);
    scan(&glb_splt_cnts[0],&glb_splt_disp[0],npes,p);
    glb_splt_count = glb_splt_cnts[npes-1] + glb_splt_disp[npes-1];
    int glb_splitters[glb_splt_count];
    MPI_Allgatherv(&    splitters[0], splt_count, MPI_INT,
                   &glb_splitters[0], &glb_splt_cnts[0], &glb_splt_disp[0],
                   MPI_INT, comm);

    // rank splitters. O( log(N/p) + log(p) )
    int disp[glb_splt_count];
    for (i=0;i<glb_splt_count;i++){
        disp[i]=0;
    }

    if(nelem>0){
#pragma omp parallel for
        for(size_t i=0; i<glb_splt_count; i++){
            disp[i] = lower_bound(arr, nelem, glb_splitters[i]) - arr[0];
        }
    }

    int glb_disp[glb_splt_count];
    for (i=0;i<glb_splt_count;i++){
        glb_disp[i]=0;
    }
    MPI_Allreduce(&disp[0], &glb_disp[0], glb_splt_count, MPI_INT, MPI_SUM, comm);

#pragma omp parallel for
    for (int qq=0; qq<kway; qq++) {
        int* _disp = &glb_disp[0];
        int optSplitter = ((qq+1)*totSize)/(kway+1);
        for(size_t i=0; i<glb_splt_count; i++) {
            if(labs(glb_disp[i] - optSplitter) < labs(*_disp - optSplitter)) {
                _disp = &glb_disp[i];
            }
        }
        spl[qq] = glb_splitters[_disp - &glb_disp[0]];
    }

    // Uncomment to check function
//    if (rank==0){
//        print_array_in_process(spl,kway,p,rank,"Parallel splitters");
//    }
}

void RankSampleSplitters(int *Ar,int *Spl,int n, int p, int rank,MPI_Comm comm){
    /*
     * My implementation of Rank algorithm used in SampleSort (Algorithm 3.2)
     * found in paper "Hyksort: A new Variant of Hypercube Quicksort on Distributed Memory Architectures
     *
     */

    int N,j,k,l,t,x;
    N=n*p;

    int *r;                                       // Sample Splitters
    int *rAr=(int *)malloc(sizeof(int)*(N+p));   //Local histogram of sorted splitters
    int *GlrAr=(int *)malloc(sizeof(int)*(N+p)); //Global histogram of sorted splitters

    r=SampleSplitters(Ar,n,p,rank);

    // Local Histogram Creation
    // Each process counts how many numbers of Array Ar are
    // a| smaller than or equal to splitter[0],
    // b| between [splitter[i] and splitter[i+1]) , i>0 and i<(n-1)
    // c| greater than or equal to splitter[n-1]
    for (l=0;l<n;l++)
    {   rAr[l]=0;
        for (t=0;t<n+1;t++){
            if (l==0){
                if (Ar[t]<=r[l]){rAr[l]++;}
            }
            else if (l>0 && l<n){
                if (Ar[t]<=r[l] && Ar[t]>r[l-1]){
                    rAr[l]++;
                }
            }
            else if (l==n){
                if (Ar[t]>r[l-1]){
                    rAr[l]++;
                }
            }
        }
    }
    MPI_Barrier(comm);
    // Uncomment to check local histogram
//    print_array_in_process(rAr, n, p, rank, "Local Histogram of sorted_splitters");

    // Reduction
    MPI_Allreduce( rAr, GlrAr, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

    // Uncomment to check global histogram
//    if (rank==0){
//        print_array_in_process(GlrAr,n,p,rank,"Global Histogram of sorted splitters");
//    }

    MPI_Barrier(comm);

    // Find k-1 splitters from global histogram
    int jump=0;
    int cnt=0;
    for (j = 0; j < k - 1; j++) {
        for (x = jump; x < n; x++) {
            if ((cnt + GlrAr[x]) <= n) {
                cnt = cnt + GlrAr[x];
            } else {
                cnt = 0;
                Spl[j] = r[x];
                jump=x;
                break;
            }
        }
    }
}

int * SampleSplitters(int *Ar,int n, int p, int rank) {
    /*
     * My implementation of SampleSplitters (Algorithm 3.3)
     * found in paper "Hyksort: A new Variant of Hypercube Quicksort on Distributed Memory Architectures
     * The Bitonic sort used in this algorithm is forked from https://github.com/steremma/Bitonic
     */
    int offset=n/p; // How many elements between each splitter
    int i,j;
    j=0;

    int *sample_splitters=(int *)malloc(sizeof(int)*(offset));
    int *sorted_splitters;

    srand(time(NULL) + rank);

    for (i = 0; i < n; i=i+offset) {
        sample_splitters[j]=Ar[i];
        j++;
        }
    // Uncomment to check functionality of Sample splitters
//    print_array_in_process(sample_splitters,offset,p,rank,"Sample splitters");

    bitonic_sort(sample_splitters,offset,p,rank, MPI_COMM_WORLD);
    // Uncomment to check functionality of Bitonic sorted splitters
//    print_array_in_process(sample_splitters,offset,p,rank,"Bitonic sorted splitters");

    MPI_Allgather(&sample_splitters[0],offset,MPI_INT,&sorted_splitters[0],offset,MPI_INT,MPI_COMM_WORLD);

    // Uncomment to check functionality of Gathered sorted_splitters
//    if(rank==0) {
//        print_array_in_process(sorted_splitters, n, p, rank, "Gathered sorted_splitters");
//    }
    return sorted_splitters;
}
void scan(int *A, int *B,int cnt, int p){
    if(cnt<100*p){
        for(int i=1;i<cnt;i++)
            B[i]=B[i-1]+A[i-1];
        return;
    }

    int step_size=cnt/p;

#pragma omp parallel for
    for(int i=0; i<p; i++){
        int start=i*step_size;
        int end=start+step_size;
        if(i==p-1) end=cnt;
        if(i!=0)B[start]=0;
        for(int j=start+1; j<end; j++)
            B[j]=B[j-1]+A[j-1];
    }

    int *sum = (int *)malloc(p*sizeof(int));
    sum[0]=0;
    for(int i=1;i<p;i++)
        sum[i]=sum[i-1]+B[i*step_size-1]+A[i*step_size-1];

#pragma omp parallel for
    for(int i=1; i<p; i++){
        int start=i*step_size;
        int end=start+step_size;
        if(i==p-1) end=cnt;
        int sum_=sum[i];
        for(int j=start; j<end; j++)
            B[j]+=sum_;
    }
    free(sum);
}

/** Bitonic Sort **/
/*
 * Bitonic sort algorithm and the folloing helper functions
 * are forked from https://github.com/steremma/Bitonic
 */

int cmpfuncASC (const void * a, const void * b);
int cmpfuncDESC (const void * a, const void * b);

void swap(int *a, int *b);
void compute_top(int partner, int dir,int *a,int rank, int k,int n);
void compute_bottom(int partner,int dir,int *a,int rank , int k, int n);

void bitonic_sort(int *Ar,int n,int p, int rank,MPI_Comm comm) {

    int  j, k, dir;

    //start counting time
    if ((rank%2 == 0) ) qsort(Ar, n, sizeof(int), cmpfuncASC); //ASCENDING
    else qsort(Ar, n, sizeof(int), cmpfuncDESC); //DESCENDING
    MPI_Barrier(MPI_COMM_WORLD);

    for (k=1;k<p;k<<=1){
        dir = ((k*2 & rank) == 0);

        for (j=k;j>=1;j>>=1){
            int partner = rank^j;
            if(rank<partner) compute_bottom(partner,dir,Ar,rank,j,n);
            else compute_top(partner,dir,Ar,rank,j,n);
            MPI_Barrier(MPI_COMM_WORLD);
        }
        if (dir) qsort(Ar, n, sizeof(int), cmpfuncASC); //ASCENDING
        else qsort(Ar, n, sizeof(int), cmpfuncDESC); //DESCENDING
    }
}

void swap(int *a, int *b){

    int t;
    t = *a;
    *a = *b;
    *b = t;
}
/** functions used by bitonic sort**/
void compute_bottom(int partner, int dir,int *a,int rank, int k,int n){
    MPI_Status mpistat;
    int rc,i;
    int *swap_temp = (int *)malloc(n/2*sizeof(int));
    //tag = 1 means you need to proccess the message.
    //since my job is to compute the top part, im sending the bottom for the other guy.
    rc = MPI_Send(a+n/2,n/2,MPI_INT,partner,1,MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS) {
        printf ("Error idle. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    //i need to compare this to my bottom half.
    MPI_Recv(swap_temp,n/2,MPI_INT,partner,1,MPI_COMM_WORLD,&mpistat);
    if(dir){
        for(i=0;i<n/2;i++){
            if(a[i] > swap_temp[i]) swap(&a[i],&swap_temp[i]);
        }

    }
    else{
        for(i=0;i<n/2;i++){
            if(a[i] < swap_temp[i]) swap(&a[i],&swap_temp[i]);
        }
    }
    //this guy sends the results before receiving
    rc = MPI_Send(swap_temp,n/2,MPI_INT,partner,2,MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS) {
        printf ("Error idle. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    //receiving the computed top part from the partner.
    MPI_Recv(a+n/2,n/2,MPI_INT,partner,2,MPI_COMM_WORLD,&mpistat);
}
void compute_top(int partner, int dir,int *a,int rank, int k,int n) {
    MPI_Status mpistat;
    int rc,i;
    int *swap_temp = (int *)malloc(n/2*sizeof(int));
    //this guy will receive before sending. If they both attemp to send first it may lead to a deadlock.
    MPI_Recv(swap_temp,n/2,MPI_INT,partner,1,MPI_COMM_WORLD,&mpistat);
    rc = MPI_Send(a,n/2,MPI_INT,partner,1,MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS) {
        printf ("Error idle. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    if(dir){
        for(i=0;i<n/2;i++){
            if(a[i+n/2] < swap_temp[i]) swap(&a[i+n/2],&swap_temp[i]);
        }
    }
    else{
        for(i=0;i<n/2;i++){
            if(a[i+n/2] > swap_temp[i]) swap(&swap_temp[i],&a[i+n/2]);
        }
    }
    //Proccessing has ended. This guy will first receive and then send to avoid deadlock.
    MPI_Recv(a,n/2,MPI_INT,partner,2,MPI_COMM_WORLD,&mpistat);
    rc = MPI_Send(swap_temp,n/2,MPI_INT,partner,2,MPI_COMM_WORLD);
}
int cmpfuncASC (const void * a, const void * b)
{
    return ( *(int*)a - *(int*)b );
}
int cmpfuncDESC (const void * a, const void * b)
{
    return ( *(int*)b - *(int*)a );
}





int SendRecieve_kway(int *Buckets,int *BucketBuffer,int n,int p,int rank,int k,MPI_Comm comm)
{
    /*
    * My implementation of AllToAll_Kway (Algorithm 3.4) found in
    * paper "Hyksort: A new Variant of Hypercube Quicksort on Distributed Memory Architectures
    */

    int i;
    int tag1=100;
    k=4;

    MPI_Status status;
    MPI_Request req;

    int N=n*p;
    int m=p/k;
    //Οσες διεργασιες εχουν το ιδιο color θεωρουμε το ανοικουν στο ιδιο υποσυνολο
    int pr=rank;
    int precv,psend;

    if(p>1)
    {
        int color = k*pr/p;
        MPI_Request R[k];

    #pragma omp parallel for
        for (i = 0; i < k; ++i)
        {   //MPI_Request
            precv=m*((color-i)%k)+(pr%m);
            MPI_Irecv(&BucketBuffer[0], 1, MPI_INT,precv, tag1, comm, &R[i]);
        }

        for (i = 0; i < k; ++i)
        {
            precv=m*((color-i)%k)+(pr%m);
            psend=m*((color+i)%k)+(pr%m);
            MPI_Isend(&Buckets[0],1,MPI_INT,psend,tag1,comm,&R[i]);

            int j=2;
            while(i>0 && i%j==0)
            {
                merge(&R[i-j],&R[i-(j/2)],&R[i-j], i-j, i-j/2);
                j=2*j; //Do something more here
            }
            MPI_Wait(&R[i], MPI_STATUS_IGNORE);
            //MPI_WaitRecv();
        }
        MPI_Waitall(k,R, MPI_STATUSES_IGNORE);

        merge(&R[0],&R[k/2],BucketBuffer,0,k/2);

        MPI_Comm scomm;
        MPI_Comm_split(comm,color,rank,&scomm);
        pr=MPI_Comm_rank(comm,&rank);
        printf("My pr is %d",pr);
    }
}