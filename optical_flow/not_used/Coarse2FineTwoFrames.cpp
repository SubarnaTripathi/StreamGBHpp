
#include "project.h"
#include "OpticalFlow.h"
#include <iostream>

using namespace std;

//void Coarse2FineTwoFrames(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
void Coarse2FineTwoFrames(DImage Im1, DImage Im2, DImage warpI2, double alpha, double ratio, int minWidth, int nOuterFPIterations, int nInnerFPIterations, int nSORIterations)
{		
	//subarna
	//DImage Im1,Im2;

	//printf("width %d   height %d   nchannels %d\n",Im1.width(),Im1.height(),Im1.nchannels());
	//printf("width %d   height %d   nchannels %d\n",Im2.width(),Im2.height(),Im2.nchannels());
	
	// get the parameters: default parameter values
	/*alpha= 1;
	ratio=0.5;
	minWidth= 40;
	nOuterFPIterations = 3;
	nInnerFPIterations = 1;
	nSORIterations= 20;*/
	
	DImage vx,vy; //,warpI2;
	OpticalFlow::Coarse2FineFlow(vx,vy,warpI2,Im1,Im2,alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations);

	// output the parameters
	//subarna
	//vx.OutputToMatlab(plhs[0]);
	//vy.OutputToMatlab(plhs[1]);

	//subarna
	//if(nlhs>2)
		//warpI2.OutputToMatlab(plhs[2]);
}