#include "mex.h"
#include "math.h"
#include "affineTransferDistanceMEX.h"
#include "stdio.h"
//  [ edgeAttr ] = computeEdgeAttr( affMatrix )
//   This generates the matrix of features, based on normailzed polar coordinate.
//
//   * input     
//
//           frame   : (double)( 6 by nFeat ) matrix representing the affine
//                          transformation from a unit circle to a region feature 
//                          ( Translation 2 + Affine 4 )
//   * output
//
//           edgeAttr   : ( 3 by nFeat * nFeat ) matrix of edge attributes
//                        rho, theta, scale

void getGraphAttr_NormPolar(double* pFrame, int nFeat, double* pOut)
{
    int i, offset_i, j, offset_j, nDim;
	double Ti[6], Ti_inv[6];
    double Pj[2], Pj_unit[2];
    double theta, rho, scale, det_i, det_j;
    double *dst;
    
    nDim = 3;
    for(i = 0; i < nFeat; i++)
    {
        offset_i = i * 6;
        Ti[0] = pFrame[offset_i+2]; // a11
        Ti[1] = pFrame[offset_i+4]; // a12  
        Ti[2] = pFrame[offset_i];   // x
        Ti[3] = pFrame[offset_i+3]; // a21
        Ti[4] = pFrame[offset_i+5]; // a22
        Ti[5] = pFrame[offset_i+1]; // y
        det_i = abs(Ti[0]*Ti[4] - Ti[1]*Ti[3]);
        MatInv(Ti,Ti_inv);
        
        for(j = 0; j < nFeat; j++)
        {
            if(i == j) 
            {
                theta= 0;
                rho = 0;
                scale= 0;
            }
            else
            {   
                offset_j = j * 6;
                Pj[0] = pFrame[ offset_j ]; Pj[1] = pFrame[ offset_j  + 1];

                MatMul(Ti_inv,Pj,Pj_unit,2);
                theta = atan2( Pj_unit[1], Pj_unit[0]);
                rho = sqrt(Pj_unit[0]*Pj_unit[0] + Pj_unit[1]*Pj_unit[1]);

                det_j = abs(pFrame[offset_j+2]*pFrame[offset_j+5] - pFrame[offset_j+4]*pFrame[offset_j+3]);
                scale = log(det_j/det_i)/log(2.0);
                
            }
            
            //dst = pOut + j*nFeat + i;
            //*dst = theta;
            //dst = dst + nFeat*nFeat;
            //*dst = rho;
            //dst = dst + nFeat*nFeat;
            //*dst = scale;
            pOut[0 + i*nDim + j*nDim*nFeat] = rho;
            pOut[1 + i*nDim + j*nDim*nFeat] = theta;
            pOut[2 + i*nDim + j*nDim*nFeat] = scale;        
        }
    }
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double *pFrame = (double *)mxGetPr(prhs[0]);
    const mwSize *dims = mxGetDimensions(prhs[0]);
  
    mwSize out[2];
    out[0] = 3;
    out[1] = dims[1]*dims[1];
            
    plhs[0] = mxCreateNumericArray(2, out, mxDOUBLE_CLASS, mxREAL);
    double* pOut = (double *)mxGetPr(plhs[0]);
            
    getGraphAttr_NormPolar(pFrame, dims[1], pOut);
}
