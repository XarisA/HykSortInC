//
// Created by xaris on 17/01/2021.
//

#ifndef HYKSORT_INC_FUNCTIONS_H
#define HYKSORT_INC_FUNCTIONS_H

#endif //HYKSORT_INC_FUNCTIONS_H

#include <stdio.h>

static int intcompare(const void *i, const void *j)
{
  if ((*(int *)i) > (*(int *)j))
    return (1);
  if ((*(int *)i) < (*(int *)j))
    return (-1);
  return (0);
}

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
int powerOf2(int n)
{
    // base cases
    // '1' is the only odd number 
    // which is a power of 2(2^0) 
    if (n == 1) 
      return 1; 
     
    // all other odd numbers are not powers of 2
    else if (n % 2 != 0 || n ==0) 
      return 0; 
     
    // recursive function call
    return powerOf2(n / 2); 
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
int lower_bound(int *a, int n, int x) {
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

int merge(int *arr1, int *arr2, int *arr3, int m, int n)
{
    int i,j,k;
    i = j = k = 0;
    for(i=0;i < m && j < n;)
    {
        if(arr1[i] < arr2[j])
        {
            arr3[k] = arr1[i];
            k++;
            i++;
        }
        else
        {
            arr3[k] = arr2[j];
            k++;
            j++;
        }
    }
    while(i < m)
    {
        arr3[k] = arr1[i];
        k++;
        i++;
    }
    while(j < n)
    {
        arr3[k] = arr2[j];
        k++;
        j++;
    }
}