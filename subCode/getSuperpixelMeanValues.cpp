#include <matrix.h>
#include <mex.h>   
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{    
    const int  *dim_array;
    int number_of_dims;
    int i,x,y,k;
    

    number_of_dims = mxGetNumberOfDimensions(prhs[0]);
    dim_array = mxGetDimensions(prhs[0]);

    int numARows = *dim_array; // h 
    int numACols = *(dim_array+1);// w 
    int numALays;
    if (number_of_dims ==2)
        numALays = 1;
    else
        numALays = *(dim_array+2);
            
    double *featureMap = mxGetPr(prhs[0]);
    int *superpixelMap = ( int * )mxGetPr(prhs[1]);
    
    int superpixelNum = ( unsigned int )( ( double * )mxGetData( prhs[ 2 ] ) )[ 0 ];    

    const int dims[2]={superpixelNum,numALays};
    mxArray *M = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
    
//     void *I1 = mxCalloc(superpixelNum*numALays,sizeof(float));
//     mxSetData(M,I1);
    
    const int dims2[1]={superpixelNum};
    mxArray *N = mxCreateNumericArray(1,dims2,mxDOUBLE_CLASS,mxREAL); 
    
//     void *I2 = mxCalloc(superpixelNum,sizeof(float));
//     mxSetData(N,I2);
 
    double *superpixelFeature = ( double * )mxGetData( M );
    double *superpixelN = ( double * )mxGetData( N );
    
    for (x=0;x<superpixelNum;x++){
        superpixelN[x] = 0;
        for ( k=0; k<numALays;k++){
                superpixelFeature[ superpixelNum*k + x ]   =  0;    
        }
    }
   
    for( x=0; x<numACols;x++ ){//lie
        for( y=0; y<numARows;y++ ){ //hang 
            superpixelN[  superpixelMap[numARows*x + y]-1  ]++;
            for ( k=0; k<numALays;k++){
                superpixelFeature[ superpixelNum*k + superpixelMap[numARows*x + y]-1 ] += featureMap[numARows*numACols*k + numARows*x + y];             
            }
        }
    }
    
    for (x=0;x<superpixelNum;x++){
        for ( k=0; k<numALays;k++){
            if (superpixelN[x]!=0)
                superpixelFeature[ superpixelNum*k + x ]   =  superpixelFeature[ superpixelNum*k + x ]/ superpixelN[x];    
        }
    }
    
    plhs[ 0 ] = M;
}

