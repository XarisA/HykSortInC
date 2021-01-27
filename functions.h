//
// Created by xaris on 17/01/2021.
//

#ifndef HYKSORT_INC_FUNCTIONS_H
#define HYKSORT_INC_FUNCTIONS_H

#endif //HYKSORT_INC_FUNCTIONS_H

#include <stdio.h>

int comparator (const void * p1, const void * p2)
{
    return (*(int*)p1 - *(int*)p2);
}
void local_sort(int *A, int n)
{
    qsort(A, n, sizeof(int), comparator);
}
void init_array(int *A,int n,int a,int b)
{
    //random number on every run between a and b for array of n size
    int i;

    srand(time(NULL));

    for(i=0; i<n; i++)
    {
        A[i]=a+rand()%(b-a+1);
    }
}
void print_array(int *A,int n)
{
    //This program prints the n values of an array
    int i;

    printf("[");
    for(i=0; i<n-1; i++)
        printf("%2d, ",A[i]);
    printf("%2d]",A[n-1]);
    printf ( "\n\n ");
}
#define VariableName(name) #name
void print_array_in_process(int *Ar, int n, int p, int rank , char msg[])
{
    printf("%s array in process %d of %d with %d elements : ",msg,rank,p,n);
    //This program prints the n values of an array
    int i;
    //
    printf("[");
    for(int i=0; i<n-1; i++){printf("%2d, ",Ar[i]);}
    printf("%2d]",Ar[n-1]);
    printf ( "\n\n");
}
int bs_lower_bound(int a[], int n, int x) {
    int l = 0;
    int h = n; // Not n - 1
    while (l < h) {
        int mid =  l + (h - l) / 2;
        if (x <= a[mid]) {
            h = mid;
        } else {
            l = mid + 1;
        }
    }
    return l;
}

