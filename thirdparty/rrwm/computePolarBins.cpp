#include "mex.h"
#include "math.h"
#include "stdio.h"
//  [ dist weight ] = computeDotProdSimilarity_sym(matchlist, edgeAttr1, edgeAttr2 )
//
//   * input     
//
//           matchlist   : (int32)( 2 by nMatch ) matrix having the pairs (feat1, feat2) for each match 
//           edgeAttr1  : (double)( k by nFeat1 * nFeat1 )  
//           edgeAttr2  : (double)( k by nFeat2 * nFeat2 )  
//
//   * output
//
//           unweighted similarity : un-weight values ( nMatch by nMatch by nAttr )

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

void computePolarBins(double* pAttrE, const int* dims, double* rho_bins, int nRho_bins, double* theta_bins, int nTheta_bins, double* rho_win, int lRho_win, double* theta_win, int lTheta_win, float* pOut)
{
    
    int i,k,l_rb,l_tb,offset,i_rb,i_tb,irb_win,itb_win;
    double rho, theta;
    
    l_rb = int(floor(double(lRho_win)/2.0));
    l_tb = int(floor(double(lTheta_win)/2.0));
    
    for( i = 0; i < dims[1]; i++ )
    {
        offset = i*(nRho_bins+nTheta_bins);
        
        rho = pAttrE[0 + i*dims[0]];
        theta = pAttrE[1 + i*dims[0]];
        //scale = pAttrE[2 + i*dims[0]];
        if ( rho == 0 && theta == 0 ) 
            continue;

        // coding
        //i_rb = find(rho < rho_bins, 1, 'first');
        i_rb = 0;
        for( k = 0; k < nRho_bins; k++ )
            if( rho_bins[k] > rho )
            {
                i_rb = k;
                break;
            }
        
        //i_tb = find(theta < [ theta_bins Inf ], 1, 'first');        
        i_tb = nTheta_bins; // the final segment is combined with the first!
        for( k = 0; k < nTheta_bins; k++ )
            if( theta_bins[k] > theta )
            {
                i_tb = k;
                break;
            }
        
        // truncated binning for distance
        //irb_win = i_rb-l_rb:i_rb+l_rb;
        //ivalid = irb_win > 0 & irb_win <= nRho_bins;
        //featBin(irb_win(ivalid),i) = rho_win(ivalid);
        irb_win = i_rb-l_rb;
        for( k = 0; k < lRho_win; k++ )
            if( irb_win + k >= 0 && irb_win + k < nRho_bins )
                pOut[irb_win+k+offset] = float(rho_win[k]);
        
        // circular binning for angle
        //itb_win = i_tb-l_tb:i_tb+l_tb;
        //itb_win(itb_win<1) = itb_win(itb_win<1) + nTheta_bins;
        //itb_win(itb_win>nTheta_bins) = itb_win(itb_win>nTheta_bins) - nTheta_bins;
        //featBin(nRho_bins + itb_win, i) = theta_win(:); 
        itb_win = i_tb-l_tb;
        for( k = 0; k < lTheta_win; k++ )
        {
            if( itb_win + k < 0 )
                pOut[itb_win+k+nTheta_bins+nRho_bins+offset] = float(theta_win[k]);
            else if( itb_win + k >= nTheta_bins )
                pOut[itb_win+k-nTheta_bins+nRho_bins+offset] = float(theta_win[k]);
            else
                pOut[itb_win+k+nRho_bins+offset] = float(theta_win[k]);
        }
        
    }
    
    return;
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double *pAttrE = (double *)mxGetPr(prhs[0]);
    const int *dims0 = mxGetDimensions(prhs[0]);
            
    double *rho_bins = (double *)mxGetPr(prhs[1]);
    const int *dims1 = mxGetDimensions(prhs[1]);
    double *theta_bins = (double *)mxGetPr(prhs[2]);
    const int *dims2 = mxGetDimensions(prhs[2]);
    
    double *rho_win = (double *)mxGetPr(prhs[3]);
    const int *dims3 = mxGetDimensions(prhs[3]);
    double *theta_win = (double *)mxGetPr(prhs[4]);
    const int *dims4 = mxGetDimensions(prhs[4]);
    
    int nRho_bins, nTheta_bins, lRho_win, lTheta_win;
    int out[2];
    
    nRho_bins = dims1[0];
    nTheta_bins = dims2[0];
    lRho_win = dims3[0];
    lTheta_win = dims4[0];
    
    out[0] = nRho_bins + nTheta_bins;
    out[1] = dims0[1];
    
    plhs[0] = mxCreateNumericArray(2, out, mxSINGLE_CLASS, mxREAL);
    //plhs[0] = mxCreateNumericArray(3, out, mxDOUBLE_CLASS, mxREAL);;
    float* pOut = (float *)mxGetPr(plhs[0]);
    
    computePolarBins(pAttrE, dims0, rho_bins, nRho_bins, theta_bins, nTheta_bins, rho_win, lRho_win, theta_win, lTheta_win, pOut);
}