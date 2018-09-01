#include "mex.h"
#include "mexOliUtil.h"
#include "math.h"
#include "affineTransferDistanceMEX.h"
#include "stdio.h"
//  [ dist flip ] = affineTransferDistanceMEX(matchlist, affMatrix1, affMatrix2, bSingle, bReflective)
//   This generates a matrix with the 'affine symmetric transfer error' 
//   between all the possible combinations of two matches
//   -- * a version assuming the left (image 1) is a reference
//
//   * input     
//
//           matchlist   : (int32)( nMatch by 2 ) matrix having the pairs (feat1, feat2) for each match 
//           affMatrix1  : (double)( nMatch by 9 ) matrix representing the affine
//                          transformation from feat1 to a unit circle
//                          * ( 3 by 3 ) is row-vectorized into (1 by 9) 
//           affMatrix2  : (double)( nMatch by 9 ) matrix representing the affine
//                          transformation from feat2 to a unit circle 
//           bSingle     : flag for matching within a single image
//           bReflective : flag for reflective matching 
//
//   * output
//
//           distance    : ( nMatch by nMatch ) matrix of error distances
//           flip        : binary matrix notifying flips for matching 
//                         ( not needed in matching between two views )
//
//
//  The way to obtain the affMatrix1 (or 2) from the ellipse of [ x y a b c o(rientation)] 
//  -----                                        
//  A = [ a b; b d ]^(-0.5); 
//  R = [ cos(o)) -sin(o); sin(o) cos(o) ];
//  AR = A*R;
//  affMatrix(iFeat,:) = [AR(1,1), AR(1,2), feat(iFeat,1), AR(2,1), AR(2,2), feat(iFeat,2), 0, 0, 1]
//
//  Using affMatrix1 and affMatrix2, the aff homography from feat1 to feat2 can be obtained!
//

void affTransDist(int* pMatch, int nMatch, double* pAff1, int nFeat1, double* pAff2, int nFeat2, double* pBoolSing, double* pBoolRefl, double* pDist, double* pFlip)
{
    
    int i,j,idx1i,idx2i,idx1j,idx2j;
    double T1i[6], T2i[6], T12i[6], T21i[6], T1i_inv[6], T2i_inv[6], T_Refl[6], Temp[6];
    double P1j[2], P2j[2], P1j_tran[2], P2j_tran[2];
    double dist1, dist2;
    bool bSing, bRefl;
    bSing = (bool)(*pBoolSing);
    bRefl = (bool)(*pBoolRefl);
    
    double* pDistFlip;
    if(bSing) pDistFlip = new double [nMatch*nMatch];
    
    for(i = 0; i < nMatch*nMatch; i++) pDist[i] = 0;
    if(bSing)
        for(i = 0; i < nMatch*nMatch; i++) pDistFlip[i] = 0;
    T_Refl[0] = 1;
    T_Refl[1] = 0;
    T_Refl[2] = 0;
    T_Refl[3] = 0;
    T_Refl[4] = -1;
    T_Refl[5] = 0;
    
    for(i = 0; i < nMatch; i++)
    {
        idx1i = pMatch[i]-1;
        idx2i = pMatch[i+nMatch]-1;
        T1i[0] = pAff1[idx1i];
        T1i[1] = pAff1[idx1i+1*nFeat1];
        T1i[2] = pAff1[idx1i+2*nFeat1];
        T1i[3] = pAff1[idx1i+3*nFeat1];
        T1i[4] = pAff1[idx1i+4*nFeat1];
        T1i[5] = pAff1[idx1i+5*nFeat1];
        
        T2i[0] = pAff2[idx2i];
        T2i[1] = pAff2[idx2i+1*nFeat2];
        T2i[2] = pAff2[idx2i+2*nFeat2];
        T2i[3] = pAff2[idx2i+3*nFeat2];
        T2i[4] = pAff2[idx2i+4*nFeat2];
        T2i[5] = pAff2[idx2i+5*nFeat2];
        
        for(j = 0; j < nMatch; j++)
        {
            if(i == j) continue;
            idx1j = pMatch[j]-1;
            idx2j = pMatch[j+nMatch]-1;
            
            P1j[0] = pAff1[idx1j+2*nFeat1]; P1j[1] = pAff1[idx1j+5*nFeat1];
            P2j[0] = pAff2[idx2j+2*nFeat2]; P2j[1] = pAff2[idx2j+5*nFeat2];
            MatInv(T1i,T1i_inv); MatInv(T2i,T2i_inv);
            if(bRefl)
            {
                MatMul(T2i,T_Refl,Temp,1);
                MatMul(Temp,T1i_inv,T12i,1);
                MatMul(T1i,T_Refl,Temp,1);
                MatMul(Temp,T2i_inv,T21i,1);
            }
            else
            {
                //MatMul(T2i,T1i_inv,T12i,1);
                MatMul(T1i,T2i_inv,T21i,1);
            }
            //MatMul(T12i,P1j,P1j_tran,2);
            MatMul(T21i,P2j,P2j_tran,2);
            dist1 = (P1j[0]-P2j_tran[0])*(P1j[0]-P2j_tran[0])+(P1j[1]-P2j_tran[1])*(P1j[1]-P2j_tran[1]);
            //dist2 = (P2j[0]-P1j_tran[0])*(P2j[0]-P1j_tran[0])+(P2j[1]-P1j_tran[1])*(P2j[1]-P1j_tran[1]);
            //pDist[i+nMatch*j] = (sqrt(dist1)+sqrt(dist2))/2;
            pDist[i+nMatch*j] = sqrt(dist1);
            //if(bSing)
            //{
            //    MatMul(T21i,P1j,P1j_tran,2);
            //    MatMul(T12i,P2j,P2j_tran,2);
            //    dist1 = (P1j[0]-P2j_tran[0])*(P1j[0]-P2j_tran[0])+(P1j[1]-P2j_tran[1])*(P1j[1]-P2j_tran[1]);
            //    dist2 = (P2j[0]-P1j_tran[0])*(P2j[0]-P1j_tran[0])+(P2j[1]-P1j_tran[1])*(P2j[1]-P1j_tran[1]);
            //    pDistFlip[i+nMatch*j] = (sqrt(dist1)+sqrt(dist2))/2;
            //}
        }
    }

    for(j = 0; j < nMatch-1; j++)
        for(i = j+1; i < nMatch; i++)
        {
            pDist[j+nMatch*i] = (pDist[j+nMatch*i]+pDist[i+nMatch*j])/2;
            pDist[i+nMatch*j] = pDist[j+nMatch*i];
        }
    
    //if(bSing)
    //{
    //    for(j = 0; j < nMatch-1; j++)
    //        for(i = j+1; i < nMatch; i++)
    //        {
    //            pDistFlip[j+nMatch*i] = (pDistFlip[j+nMatch*i]+pDistFlip[i+nMatch*j])/2;
    //            pDistFlip[i+nMatch*j] = pDistFlip[j+nMatch*i];
    //        }
    //    for(i = 0; i < nMatch*nMatch; i++)
    //        pFlip[i] = 0;
    //    for(i = 0; i < nMatch*nMatch; i++)
    //        if(pDist[i] > pDistFlip[i])
    //            pFlip[i] = 1;
    //    for(i = 0; i < nMatch*nMatch; i++)
    //        if(pFlip[i])
    //            pDist[i] = pDistFlip[i];
    //}
    
    //if(bSing) delete [] pDistFlip;
    return;
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    enum{match, aff1, aff2, bSing, bRefl};
    enum{dist, flip};
    
    oliCheckArgNumber(nrhs,5,nlhs,2);
    
    int nMatch, nFeat1, nFeat2;
    
    int* pMatch = (int*)oliCheckArg(prhs,match,&nMatch,2,oliInt);
    double* pAff1 = oliCheckArg(prhs,aff1,&nFeat1,9,oliDouble);
    double* pAff2 = oliCheckArg(prhs,aff2,&nFeat2,9,oliDouble);
    double* pBoolSing = oliCheckArg(prhs,bSing,1,1,oliDouble);
    double* pBoolRefl = oliCheckArg(prhs,bRefl,1,1,oliDouble);
        
    plhs[dist] = mxCreateDoubleMatrix(nMatch, nMatch, mxREAL);
    double* pDist = mxGetPr(plhs[dist]);
    plhs[flip] = mxCreateDoubleMatrix(nMatch, nMatch, mxREAL);
    double* pFlip = mxGetPr(plhs[flip]);
            
    affTransDist(pMatch, nMatch, pAff1, nFeat1, pAff2, nFeat2, pBoolSing, pBoolRefl, pDist, pFlip);
}