/* This is largely a copy of fastranc.c, created by Rami K. Niazy.
 
  BK added a few options to monitor convergence, to allow the specification of a
  starting guess, and flipped the estimated filter so that it can be used 
  directly in a convolution operation in Matlab.

  See anc.m for an example how to call this function.

*/

#include "math.h"
#include "mex.h"

void ancEstimate(int N, double mu, int nrSamples, double *refs, double *signal, 
    double *filtered, double *noiseEstimate, double *filter,double *delta,double *initial) {

    int i, j;
    double tmp;

    // initialize the filter.
    for (i = 0; i <= N; i++) {
        filter[i]=initial[i];
        if (delta) delta[i]=0;
    }
    //loop over samples
    for (i = 0; i < nrSamples; i++) {

        if (i < N) { 
            // first N samples are untouched.
            filtered[i]=signal[i]; noiseEstimate[i]=0;
            if (delta)  delta[i] =0;               
        } else {
            // Estimate current noise by filtering the reference
            noiseEstimate[i]=0;            
            for(j = 0; j <= N; j++) {
                noiseEstimate[i] += filter[N-j]*refs[i-N+j]; 
            }
            // Subtract noise from signal
            filtered[i]=signal[i]-noiseEstimate[i];
            // Determine how much to change the filter
            for (j = 0; j <= N; j++) {
                tmp = 2*mu*filtered[i]*refs[i-N+j];
                if (delta){
                    // Keep track of the mse change
                    delta[i]+= sqrt(tmp*tmp/((double) N+1.0));
                }
                filter[N-j] += tmp;
            }            
        }
    }
}

//[filtered,noiseEstimate,filter,delta] = ancEstimate(reference,signal,N,mu,initial)
void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {

    double *refs, *signal, *filtered, *noiseEstimate, mu,*filter,*delta,*initial;
    int    N, nrSamples;

    /*  Check for proper number of arguments. */
    if (nrhs != 5) 
        mexErrMsgTxt("Five inputs required: (reference,signal,N,mu,initial) ");
    if (nlhs > 4 ) 
        mexErrMsgTxt("Four or fewer outputs required: [filtered,noiseEstimate,filter,delta]");

    /* Check to make sure the N input argument is a scalar. */
    if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
        mxGetN(prhs[2])*mxGetM(prhs[2]) != 1) 
        mexErrMsgTxt("Input N must be a scalar of type double.");
    
    /* Check to make sure the mu input argument is a scalar. */
    if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
        mxGetN(prhs[3])*mxGetM(prhs[3]) != 1) 
        mexErrMsgTxt("Input mu must be a scalar of type double.");

    /* Check that refs and signal are doulbe */
    if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]))
        mexErrMsgTxt("Inputs must be of type double.");

    /* Check to make sure length refs and d is same and are col vectors. */
    if (mxGetN(prhs[0]) != 1 ||  mxGetN(prhs[1]) != 1)
        mexErrMsgTxt("Reference and Input data must be column vectors.");
    if (mxGetM(prhs[0]) != mxGetM(prhs[1])) 
        mexErrMsgTxt("Reference and Input must be of the same length.");


    /* Get inputs */
    refs = mxGetPr(prhs[0]);
    nrSamples = mxGetM(prhs[0]);
    
    signal    = mxGetPr(prhs[1]);
    N = (int) mxGetScalar(prhs[2]);
    mu = mxGetScalar(prhs[3]);
    initial = mxGetPr(prhs[4]);      
    

    if (mxGetM(prhs[4]) != (N+1) ||  mxGetN(prhs[4]) != 1){
        mexPrintf("[%d %d]", mxGetM(prhs[4]),mxGetN(prhs[4]));
        mexErrMsgTxt("Initial filter must be column vector of size [N+1 1]");
    }

    /* Set the output pointer to the output matrix. */
    filtered = (double *) mxMalloc(nrSamples * sizeof(double));
    noiseEstimate = (double *) mxMalloc(nrSamples * sizeof(double));
    filter = (double *) mxMalloc((N+1) * sizeof(double));
    if (nlhs >3){
        delta =  (double *) mxMalloc(nrSamples * sizeof(double));
    }else{
        delta = NULL;
    }

    
    
    /* Call the C subroutine. */
    ancEstimate(N,mu, nrSamples, refs, signal, filtered, noiseEstimate, filter,delta,initial);

    /* Assign to the Matlab output variables*/
    if (nlhs>0){
        plhs[0] = mxCreateDoubleMatrix(nrSamples,1, mxREAL);
        mxSetPr(plhs[0], filtered);
    }   
    if (nlhs>1){
        plhs[1] = mxCreateDoubleMatrix(nrSamples,1, mxREAL);
        mxSetPr(plhs[1], noiseEstimate);
    }
    if (nlhs>2){
        plhs[2] = mxCreateDoubleMatrix(N+1,1, mxREAL);    
        mxSetPr(plhs[2], filter);
    }
    if (nlhs>3){
        plhs[3] = mxCreateDoubleMatrix(nrSamples, 1,mxREAL);
        mxSetPr(plhs[3], delta);      
    }
return;
}
