// File : dynam.cpp

#define S_FUNCTION_NAME InvPend_nonlinear
#define S_FUNCTION_LEVEL 2

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "simstruc.h"
#include "InvPend_nonlinear.h"

#define F (U(0)) // Force

/*====================*
 * S-function methods *
 *====================*/

/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *    The sizes information is used by Simulink to determine the S-function
 *    block's characteristics (number of inputs, outputs, states, etc.).
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, 1);  /* Number of expected parameters */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return; /* Parameter mismatch will be reported by Simulink */
    }

    ssSetNumContStates(S, 4);
    ssSetNumDiscStates(S, 0);

    if (!ssSetNumInputPorts(S, 1)) return;
     ssSetInputPortWidth(S, 0, 1);
     
    ssSetInputPortDirectFeedThrough(S, 0, 1);

    if (!ssSetNumOutputPorts(S, 1)) return;
   	ssSetOutputPortWidth(S, 0, 4);
    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    /* Take care when specifying exception free code - see sfuntmpl.doc */
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}

/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    S-function is comprised of only continuous sample time elements
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_INITIALIZE_CONDITIONS
/* Function: mdlInitializeConditions ========================================
 * Abstract:
 *    Initialize both continuous states to zero
 */
static void mdlInitializeConditions(SimStruct *S)
{
    real_T *x0 = ssGetContStates(S);
	InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S, 0);

	x0[0] = mxGetPr(ssGetSFcnParam(S, 0))[0];
	x0[1] = mxGetPr(ssGetSFcnParam(S, 0))[1];

	x0[2] = mxGetPr(ssGetSFcnParam(S, 0))[2]; 
	x0[3] = mxGetPr(ssGetSFcnParam(S, 0))[3];

}

/* Function: mdlOutputs =======================================================
 * Abstract:
 *      y = Cx + Du 
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
    real_T *x    = ssGetContStates(S);
    real_T  *y = ssGetOutputPortRealSignal(S,0);
	InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S, 0);    

	for(int i = 0; i<4; i++) y[i] = x[i];
}

#define MDL_DERIVATIVES
/* Function: mdlDerivatives =================================================
 * Abstract:
 *      xdot(x0) = x0*(1-x1^2) -x1
 *      xdot(x1) = x0
 */
static void mdlDerivatives(SimStruct *S)
{
	real_T  t = ssGetT(S);
    real_T *dx= ssGetdX(S);
    real_T *x = ssGetContStates(S);
	InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0);
    
	double M[2][2], Minv[2][2], c[2][1], tmp[2];
    double yy,dy,theta,dtheta;
    
    yy     = x[0];
    dy     = x[1];
    theta  = x[2];
    dtheta = x[3];

	*(*(M+0)+0) = Mc + Mr;
	*(*(M+0)+1) = 0.5*Mr*l*cos(theta);
	*(*(M+1)+0) = 0.5*Mr*l*cos(theta);
	*(*(M+1)+1) = Mr*SQR(l)/4.0;

	*(*(c+0)+0) = 0.5*Mr*l*SQR(dtheta) + F - b*dy;
	*(*(c+1)+0) = 0.5*Mr*l*g0*sin(theta);

	inverse(2, M, Minv);
	MatMult(2, 1, 2, *Minv, *c, tmp);

    dx[0] = dy;
	dx[1] = tmp[0];

	dx[2] = dtheta;
	dx[3] = tmp[1];   
   
}

/* Function: mdlTerminate =====================================================
 * Abstract:
 *    No termination needed, but we are required to have this routine.
 */
static void mdlTerminate(SimStruct *S)
{
}

#ifdef  MATLAB_MEX_FILE    // Is this file being compiled as a MEX-file?
#include "simulink.c"      // MEX-file interface mechanism
#else
#include "cg_sfun.h"       // Code generation registration function
#endif