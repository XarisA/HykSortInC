#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi/mpi.h>
#include <omp.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "functions.h"

#define SIZE 64


int * ParallelSelect(int *Ar,int n, int p, int rank) ;
int HykSort(int *Ar, int kway,int rank, MPI_Comm comm);
void bitonic_sort(int *Ar,int n,int p, int rank,MPI_Comm comm);

int main (int argc, char *argv[])
{
    int n,m,p,pr,rank,i,j,l,t,x,v,w,y;
    int tag1,tag2,tag3;
    int N=SIZE;
    //int A[SIZE];
  	int GlB[SIZE];
    int *Ar=(int *)malloc(sizeof(int)*N); 
    int *B=(int *)malloc(sizeof(int)*N);
    int *Bsorted=(int *)malloc(sizeof(int)*N);
    int *GlBBuffer = (int *) malloc (sizeof(int) * 2 * N);
    int *r=(int *)malloc(sizeof(int)*N); 
    int *rAr=(int *)malloc(sizeof(int)*N); //Local histogram of sorted splitters
    int *GlrAr=(int *)malloc(sizeof(int)*N); //Global histogram of sorted splitters
    int *Buckets = (int *) malloc (sizeof (int) * (N + p));

    int k=4;
    int *Spl=(int *)malloc(sizeof(int)*(k-1)); //Global histogram of sorted splitters
    int cnt=0;
    tag1=100;
    tag2=200;
    tag3=300;

    MPI_Comm comm;
    MPI_Status status;

    comm = MPI_COMM_WORLD;

    MPI_Init(NULL,NULL);
    MPI_Comm_size(comm,&p);
    MPI_Comm_rank(comm,&rank);

    MPI_Barrier(comm);
    // Table Creation
    if(rank==0)
    {
        n=N/p; //n=16 elements
        //int Ar[n];

        //Create my Array
        //init_array(A, N,1,100);
        int A[64]={61,25,63,59,7,11,21,3,36,18,50,66,93,86,22,40,85,94,14,8,21,45,7,42,57,72,28,56,9,87,72,69,63,87,28,69,49,100,23,36,17,73,1,9,58,74,1,42,68,14,50,40,11,8,82,67,79,61,75,40,100,98,60,62};
        printf("\t\t\t HykSort Implementation \n\n");
        printf("Generating array with random int A: ");
        print_array(A,N);
        Ar=A;

        //TODO Change the following with Broadcast
        //MPI_Bcast(&N,1,MPI_INT,0,comm);
        //Αποστολή του Ν σε όλες τις διεργασίες
        for (int target = 1; target < p; target++) {
            MPI_Send(&N, 1, MPI_INT, target, tag1, MPI_COMM_WORLD);
        }
        i= n;
        // Initial Element Distribution
        //TODO Change the following with Scater
        //Αποστολή τον πινακα σε όλες τις διεργασίες(-επεξεργαστες) (της διεύθυνσής του για την ακρίβεια).
        //Ο καθε ενας παιρνει τον ιδιο αριθμο στοιχειων
        for (int target = 1; target < p; target++) {
            MPI_Send(&A[i], n, MPI_INT, target, tag2, MPI_COMM_WORLD); 
            i += n;
        }
    }else {
        n=N/p;
        //int Ar[n];
        MPI_Recv(&N, 1, MPI_INT, 0, tag1, MPI_COMM_WORLD, &status);
        //Η κάθε διεργασία θα παραλάβει n στοιχεία και θα τα σώσει στο δικό της διάνυσμα Ar
        MPI_Recv(&Ar[0], n, MPI_INT, 0, tag2, MPI_COMM_WORLD, &status);

    }
    MPI_Barrier(comm);
    //printf("I am rank %d of %d, and I know N= %d,n = %d .",rank,p,N,n);
    print_array_in_process(Ar,n,p,rank,"Split A to Ar");
    MPI_Barrier(comm);
    // Local Sort & Sample Distribution
    // Ο καθε επεξεργαστης ταξινομει το δικο του κομματι.
    local_sort(Ar,n);
    MPI_Barrier(comm);
    print_array_in_process(Ar,n,p,rank,"Sort local Ar");
    MPI_Barrier(comm);
    r=ParallelSelect(Ar,n,p,rank);
    //print_array_in_process(r, n, p, rank, "Sorted splitters");
    MPI_Barrier(comm);

    //Βημα 2. Rank
    //Καθε διεργασια θα μετρησει ποσοι αριθμοι ειναι μικροτεροι του πρωτου Splitter (δηλ του index0)
    //και θα το σωσει στο δικο της index
    for (l=0;l<n;l++)
    {   rAr[l]=0;
        for (t=0;t<n;t++){
            if (l==0){
                if (Ar[t]<=r[l]){rAr[l]++;}
            }
            else if (l>0 && l<n-1){
                if (Ar[t]<=r[l] && Ar[t]>r[l-1] && r[l]>r[l-1]){ 
                    rAr[l]++;
                }
            }
            else if (l==n-1){
                if (Ar[t]>r[l]){ 
                    rAr[l]++;
                }
            }
        }
    }
    MPI_Barrier(comm);
    print_array_in_process(rAr, n, p, rank, "Local Histogram of sorted_splitters");

    MPI_Allreduce( rAr, GlrAr, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
    if (rank==0){
        print_array_in_process(GlrAr,n,p,rank,"Global Histogram of sorted splitters");
    }

    MPI_Barrier(comm);
    //Find k-1 splitters from global histogram

    int jump=0;
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
    
    print_array_in_process(Spl, k - 1, p, rank, "Global Splitters");

    MPI_Barrier(comm);
    //Πρεπει να φτιαξω buckets μοιραζοντας τα στοιχεια με βαση τους διαχωριστες
    //Βημα 4. M είναι το μήνυμα που πρέπει να στείλει η διεργασία με το rank pr σε κάθε άλλη διεργασία i
    //Το message του bucket που πρέπει να φτιάξει
    //Κάθε διεργασία πρέπει να χωρίσει τους αριθμούς που έχει (Ar) με βάση τους Splitters Spl
    /**** Creating Buckets locally ****/
    
    
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

    //TODO Remove The following statements

    //print_array_in_process(Buckets, N+p, p, rank, "Buckets");

    //scaterv send parts of vectors to all other processes
    //sendcounts Πόσα στοιχεία θα στείλω σε κάθε διεργασία (με την σειρά οι διεργασίες)
    //sendcounts μία θέση για κάθε διεργασία
    //displacements από ποα θέση θα ξεκινήσουμε να παίρνουμε τα στοιχεία
    //MPI_Scatterv(A,sendcounts,displs,MPI_INT,local_A,sendcounts[rank],MPI_INT,0,comm)


    int * BucketBuffer = (int *) malloc (sizeof (int) * (N + p));

    // all to all broadcast => Χρειαζομαι τον αναστροφο πινακα
    //Για αρχη Alltoallv στην συνεχεια να την αλλαξω με AlltoAllkway

    MPI_Alltoall (Buckets, n + 1, MPI_INT, BucketBuffer, 
					 n + 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Barrier(comm);
    

    //Remove displacements

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

    //sendcount πόσα στοιχεία πρέπει να πάνε σε κάθε διεργασία σε διάνυσμα
    //displacements από ποια θέση θέλουμε να στέλνουμε δεδομένα

// sendcount[i] = i + 1;
// displs1[i] = i * LENGTH;

    int Bsize =0;
    for(l=0;l<N;l++){
        if(B[l]!=0){
            Bsize++;
        }
    }
    MPI_Barrier(comm);
    //local_sort(B,Bsize);
    //print_array_in_process(B, Bsize, p, rank, "Buckets gathered");
    MPI_Barrier(comm);
    int * recvcount = (int *)malloc(N*sizeof(int));
    int * displs = (int *)malloc(N*sizeof(int));


    for (i = 0; i < N; i++) {
        recvcount[i] = i + 1;
        displs[i] = i * N;
    }

    MPI_Barrier(comm);
    MPI_Gatherv(B,recvcount[rank],MPI_INT,GlB, recvcount, displs,MPI_INT,0,comm);
    MPI_Barrier(comm);
    local_sort(GlB,N);
    if (rank == 0){
		 print_array_in_process(GlB, N, p, rank, "Buckets gathered");
    	}


    //print_array_in_process(B, N, p, rank, "True Buckets");



    free(Buckets);
    // free(Spl);
    // free(Ar);
    // free(GlrAr);
    // free(rAr);
    free(B);

    MPI_Finalize();
    return 0;
}


int * ParallelSelect(int *Ar,int n, int p, int rank) {
    int offset=n/p; //offset=4 ανα 4 στοιχεια παω να παρω splitter
    int i,j;
    j=0;

    int *sample_splitters=(int *)malloc(sizeof(int)*(offset));
    int *sorted_splitters;
    // int r[nlocal];
    // memset(r,0,nlocal*sizeof(int));

    srand(time(NULL) + rank);

    for (i = 0; i < n; i=i+offset) {
        sample_splitters[j]=Ar[i];
        j++;
        }
    print_array_in_process(sample_splitters,offset,p,rank,"Sample splitters");

    bitonic_sort(sample_splitters,offset,p,rank, MPI_COMM_WORLD);
    print_array_in_process(sample_splitters,offset,p,rank,"Bitonic sorted splitters");

    MPI_Allgather(&sample_splitters[0],offset,MPI_INT,&sorted_splitters[0],offset,MPI_INT,MPI_COMM_WORLD);

    if(rank==0) {
        print_array_in_process(sorted_splitters, n, p, rank, "Gathered sorted_splitters");
    }
    return sorted_splitters;
}

int Hyksort(int *Ar, int k, MPI_Comm comm)
{

}

    // Χωριζουμε τον πινακα του καθε επεξεργαστη με splitters.
    // αν θελω στους 3εις επεξεργαστες, βαζω δυο splitters
    //Πρεπει να επιλεξω τυχαια 4 αριθμους , οσους και τα process που εχω! αρα 2 Splitters
    //Εστω SampleSplitters και την αλλαζω στην συνεχεια με ParallelSelect


/** Bitonic Sort **/
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

