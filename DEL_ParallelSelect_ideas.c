/*
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

*/
