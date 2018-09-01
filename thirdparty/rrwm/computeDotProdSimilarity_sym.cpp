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

void computeSim(int* pMatch, int nMatch, float* pAttr1, const int* dims1, float* pAttr2, const int* dims2, float* pSimVal)
{
    
    int i,j,k,fdim,idx1i,idx2i,idx1j,idx2j, pos1, pos2, pos1_s, pos2_s, pos1_w, nFeat1, nFeat2;
    double sim, sim_s;
    
    //for(i = 0; i < nMatch; i++)
	//{
    //  offset_i = 2*i;
	//	idx1j = pMatch[offset_i];
    //    idx2j = pMatch[offset_i+1];
    //    printf("i [%d,%d] \n",idx1j,idx2j);
    //}
    //printf("wTheta: %5.3f wRho:%5.3f \n ",wTheta,wRho);
    // compute normalized polar attributes for graph 2
    //printf("[%d,%d,%d]  \n",dims1[0],dims1[1],dims1[2]);
    //printf("[%d,%d,%d]  \n",dims2[0],dims2[1],dims2[2]);
    //printf("nMatch: %d \n",nMatch);
    
    fdim = dims1[0];
    nFeat1 = int(sqrt((double)dims1[1]));
    nFeat2 = int(sqrt((double)dims2[1]));
    //nMatch = 5;
	for(i = 0; i < nMatch; i++)
	{
        idx1i = pMatch[2*i]-1;
        idx2i = pMatch[2*i+1]-1;
        //printf("match [%d,%d]\n",idx1i,idx2i);
        for(j = i+1; j < nMatch; j++)
        {
            //printf("[%d,%d] %d  \n",i,j,bFullWeight);
		
            // skip diagonals
			//if(i == j){
			//	pSimVal[nMatch*i+i] = 0;
            //    continue;
			//}
            
            idx1j = pMatch[2*j]-1;
            idx2j = pMatch[2*j+1]-1;
		
            // angle
            pos1 = idx1i+idx1j*nFeat1; // edge att pos in G1
            pos2 = idx2i+idx2j*nFeat2; // edge att pos in G2

            pos1_s = idx1j+idx1i*nFeat1; // edge att pos in G1
            pos2_s = idx2j+idx2i*nFeat2; // edge att pos in G2
            //printf("[%d,%d] %d  %f\n",idx1i,idx1j,pos1,pAttr1[pos1]);        
            //printf("[%d,%d] %d  %f\n",idx2i,idx2j,pos2,pAttr1[pos2]);        

            sim = 0;
            sim_s = 0;
            
            for(k = 0; k < fdim; k++)
            {
                sim = sim + pAttr1[k+fdim*pos1] * pAttr2[k+fdim*pos2]; 
                sim_s = sim_s + pAttr1[k+fdim*pos1_s] * pAttr2[k+fdim*pos2_s]; 
            }
            pSimVal[i+j*nMatch] = (sim + sim_s)/2.0;
            pSimVal[j+i*nMatch] = pSimVal[i+j*nMatch];
            //printf("[%d,%d] %f %f ,  %f\n",i,j,sim, sim_s, pSimVal[i+j*nMatch]);
        }
        //printf("\n");
	}
    
    return;
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    int *pMatch = (int *)mxGetPr(prhs[0]);
    const int *dims0 = mxGetDimensions(prhs[0]);
            
    float *pAttr1 = (float *)mxGetPr(prhs[1]);
    const int *dims1 = mxGetDimensions(prhs[1]);
    float *pAttr2 = (float *)mxGetPr(prhs[2]);
    const int *dims2 = mxGetDimensions(prhs[2]);
    
    int out[2];
    out[0] = dims0[1];
    out[1] = dims0[1];
    //out[2] = 3;
    
    plhs[0] = mxCreateNumericArray(2, out, mxSINGLE_CLASS, mxREAL);
    //plhs[0] = mxCreateNumericArray(3, out, mxDOUBLE_CLASS, mxREAL);;
    float* pOut = (float *)mxGetPr(plhs[0]);
    
    computeSim(pMatch, dims0[1], pAttr1, dims1, pAttr2, dims2, pOut);
}