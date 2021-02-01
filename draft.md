## To check Buckets use 

```C
if(rank==0){
MPI_Barrier(comm);
sleep(5);
print_array_in_process(Buckets, N+p, p, rank, "Buckets");
}

if(rank==1){
MPI_Barrier(comm);
sleep(10);
print_array_in_process(Buckets, N+p, p, rank, "Buckets");
}


if(rank==2){
MPI_Barrier(comm);
sleep(15);
print_array_in_process(Buckets, N+p, p, rank, "Buckets");
}


if(rank==3){
MPI_Barrier(comm);
sleep(20);
print_array_in_process(Buckets, N+p, p, rank, "Buckets");
}
```

## Hyksort ideas

```C
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi/mpi.h>
#include <omp.h>



#define SIZE 12

int A[SIZE];    //input array
int B[SIZE];    //output array
MPI_Comm comm;  //MPI communicator
int p;          //number of MPI tasks in comm
int N=SIZE;     //global number of keys in A
int n;          //local number of keys in A
int k;          //number of splitters
int s;          //the splitters
int tw;         //interconnect slowness (1/bandwidth)
int ts;         //interconnect latency
int tc;         //intranode memory slowness (1/RAM bandwidth)
int pr;         //the task id (its MPI rank)
int Ar[SIZE];   //array to be sorted (local block)
int Mrs;        //MPI message from s (send) task to r (recv) task
int partner;    //indicates task id in point-to-point exchanges

int main (int argc, char *argv[])
{

    k =4; //Change this
    int m=p/k;
    int global_sum,local_sum; //double check?
    int tag1=50;
    MPI_Status status; //Change this ??
    MPI_Request R[k];
    comm = MPI_COMM_WORLD;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &partner);
    srand(partner);

    if(partner==0)
    {
        //Create my Array
        init_array(Ar, N,1,100);
        printf("\t\t\t HykSort Implementation \n\n");
        printf("Start    : ");
        print_array(Ar,N);
        printf ( "\n\n ");
    }

    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_INT, MPI_SUM, comm);
    LocalSort(*A,N);

    //?n=N/p;
    while(p>1)
    {
        //if(k>p) k=p;
        int m=p/k;
        int pr=partner/m;
        int new_pid=partner%m;
        //s = ParallelSelect(Ar,n,N,partner);
        int color = k*pr/p;

        #pragma parallel for
        for (int i = 0; i < k; ++i) //1
        {
            //MPI_Request
            int precv=m*((color-i)%k)+(pr%m);
            MPI_Irecv(&B, N, MPI_INT,precv, tag1,comm, &R[i]);
        }
        for (int i = 0; i < k; ++i) //2
            {
                int precv=m*((color-i)%k)+(pr%m);
                int psend=m*((color+i)%k)+(pr%m);
                MPI_Isend(&B,N,MPI_INT,psend,tag1,comm,&R[i]);

                int j=2;
                while(i>0 && i%j==0) //3
                {
                    j=2*j;
                }
            MPI_Wait(&R[i],&status);
            }
        MPI_Waitall(1,&R,&status); //4

        MPI_Comm_split(comm,color,partner,comm);
        pr=MPI_Comm_rank(comm,partner);
    }
    print_array(B,N);
    printf("\n\n\nResult: %d\n", B);
    MPI_Finalize();
    return 0;
}
```

## ParallelSelect ideas

```C
#include <mpi/mpi.h>
#include <stdlib.h>

int ParallelSelect(int *A,int kway, MPI_Comm comm)
{
    int rank, npes;
    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &rank);

    int totSize,nelem = sizeof A / sizeof *A;

    MPI_Allreduce(&nelem, &totSize, 1, MPI_INT, MPI_SUM, comm);

    // Determine splitters
    int splt_count = (1000*kway*nelem)/totSize;
    if (npes>1000*kway) splt_count = (((float)rand()/(float)RAND_MAX)*totSize<(1000*kway*nelem)?1:0);
    if (splt_count>nelem) splt_count=nelem;

    int splitters[splt_count];

    for(size_t i=0;i<splt_count;i++) {
        splitters[i] = A[rand() % nelem];
    }


    // Gather all splitters
    int glb_splt_count;
    int glb_splt_cnts[npes];
    int glb_splt_disp[npes];
    for(int i = npes; i <= 0; i--){
        glb_splt_disp[i] = i;
    }
    MPI_Allgather(&splt_count,1,MPI_INT, &glb_splt_cnts[0], 1,MPI_INT, comm);
}
```

