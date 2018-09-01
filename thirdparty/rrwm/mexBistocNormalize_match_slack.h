#include "stdlib.h"

void bistocNorm(double* pX, int N,
        int* pIdx1, int* pID1, int* pIdx2, int* pID2,
        double* pTol, double* pDumDim, double* pDumVal, int* pMaxIter,
        double* pY)
{
    double tol = (*pTol);
    tol = tol*tol;
    double delta = tol;
    
    int i, nG1, nG2, iter;
    double *pX2 = new double[N];
    double *pTemp = new double[N];
    
    double tempSum = 0;
    
    // Num of Groups
    nG1 = pID1[N-1];
    nG2 = pID2[N-1];
   
    // MATLAB idx --> C++ idx
//     for(i = 0; i < N; i++) pIdx1[i]--;
//     for(i = 0; i < N; i++) pIdx2[i]--;
    
    iter = 0;
    
//     printf("%f\t%f\n", delta, tol);
//     printf("%d\t%d\n", iter, (*pMaxIter));
//     
    
    while(delta >= tol && iter < (*pMaxIter))
    {
        iter++;
                
        // copy the current state
        for(i = 0; i < N; i++)
            pX2[i] = pX[i];
        
        // update domain 1
        pTemp[0] = pX[pIdx1[0]];
        for(i = 1; i < N; i++)
        {
            if(pID1[i] == pID1[i-1])
                pTemp[i] = pTemp[i-1]+pX[pIdx1[i]];
            else
                pTemp[i] = pX[pIdx1[i]];
        }
        for(i = N-2; i >= 0; i--)
        {
            if(pID1[i] == pID1[i+1])
                pTemp[i] = pTemp[i+1];
        }
        
        for(i = 0; i < N; i++)
        {
            pX[pIdx1[i]] /= pTemp[i];
        }
        
        if((*pDumDim) == 1)
        {
            for(i = N-nG2; i < N; i++)
                pX[i] *= (*pDumVal);
        }
        
        /*for(i = 0; i < N; i++)
            printf("%d %d %f\n", pID1[i], pIdx1[i]+1, pX[pIdx1[i]]);*/
        
        // update domain 2
        pTemp[0] = pX[pIdx2[0]];
        for(i = 1; i < N; i++)
        {
            if(pID2[i] == pID2[i-1])
                pTemp[i] = pTemp[i-1]+pX[pIdx2[i]];
            else
                pTemp[i] = pX[pIdx2[i]];
        }
        for(i = N-2; i >= 0; i--)
        {
            if(pID2[i] == pID2[i+1])
                pTemp[i] = pTemp[i+1];
        }
        for(i = 0; i < N; i++)
        {
            pX[pIdx2[i]] /= pTemp[i];
        }
        
        if((*pDumDim) == 2)
        {
            for(i = N-nG1; i < N; i++)
                pX[i] *= (*pDumVal);
        }
        
        /*printf("\n");
        for(i = 0; i < N; i++)
            printf("%d %d %f\n", pID2[i], pIdx2[i]+1, pX[pIdx2[i]]);
        printf("\n");*/

        // check the difference for termination criterion
        delta = 0;
        for(i = 0; i < N; i++)
        {
            delta += (pX[i]-pX2[i])*(pX[i]-pX2[i]);
        }
    }
    
    // return solution
    for(i = 0; i < N; i++)
        pY[i] = pX[i];
    
//     printf("Iter: %d\n", iter);
//     
//     for(i = 0; i < N; i++)
//         printf("%f\t%f\n", pX[i], pY[i]);
    
    delete [] pX2;
    delete [] pTemp;
}