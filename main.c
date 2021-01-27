#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi/mpi.h>
#include <math.h>
#include <omp.h>
#include "functions.h"

#define SIZE 64


int ParallelSelect(int *Ar,int n, int p, int rank) ;
int HykSort(int *Ar, int kway,int rank, MPI_Comm comm);
void bitonic_sort(int *Ar,int n,int p, int rank,MPI_Comm comm);

int main (int argc, char *argv[])
{
    int n,p,rank,i,j,l,offset;
    int tag1,tag2;
    int N=SIZE;
    int A[SIZE];
    int *Ar=(int *)malloc(sizeof(int)*10);
    int *B=(int *)malloc(sizeof(int)*10);


    int k=4;

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
        init_array(A, N,1,100);
        printf("\t\t\t HykSort Implementation \n\n");
        printf("Generating array with random int A: ");
        print_array(A,N);
        Ar=A;

        //Αποστολή του Ν σε όλες τις διεργασίες
        for (int target = 1; target < p; target++) {
            MPI_Send(&N, 1, MPI_INT, target, tag1, MPI_COMM_WORLD);
        }
        i= n;
        // Initial Element Distribution
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
    ParallelSelect(Ar,n,p,rank);
    MPI_Barrier(comm);

    // all to all broadcast => Χρειαζομαι τον αναστροφο πινακα
    //Για αρχη Alltoallv στην συνεχεια να την αλλαξω με AlltoAllkway

//    int MPI_Alltoallv(const void *sendbuf, const int *sendcounts,
//                      const int *sdispls, MPI_Datatype sendtype, void *recvbuf,
//                      const int *recvcounts, const int *rdispls, MPI_Datatype recvtype, MPI_Comm comm)


    MPI_Finalize();
    return 0;
}
int ParallelSelect(int *Ar,int n, int p, int rank) {
    int offset=n/p; //offset=4 ανα 4 στοιχεια παω να παρω splitter
    int i,j,l,t;
    j=0;

    int *sample_splitters=(int *)malloc(sizeof(int)*(offset));
    int *sorted_splitters;
    int *r=(int *)malloc(sizeof(int)*n);

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
//    return sorted_splitters;
    //Αυτα να τα παω στη main
    //Καθε διεργασια θα μετρησει ποσοι αριθμοι ειναι μικροτεροι του πρωτου Splitter (δηλ του index0)
    //και θα το σωσει στο δικο της index0
    //TODO Create a local Histogram of sorted_splitters and Ar.
    for (l=1;l=n;l++)
    {
        for (t=0;t<n;t++){
            if (Ar[t]<sorted_splitters[l] && Ar[t]>sorted_splitters[l-1]) {
                r[l]++;
            }
        }
    }
    print_array_in_process(r, n, p, rank, "Local Histogram of sorted_splitters");
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

