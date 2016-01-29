/*
Original Code From:
Copyright (C) 2006 Pedro Felzenszwalb
Modifications (may have been made) Copyright (C) 2011,2012 
  Chenliang Xu, Jason Corso.

Modifications (may have been made) Copyright (C) 2013 
  Subarna Tripathi, UCSD.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
*/

#include <cstdio>
#include <cstdlib>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>

# ifdef WIN32
#include <time.h>
//#include "msdirent.h"
#include <direct.h>
# else
#include <dirent.h>
# endif

#include "image.h"
#include "pnmfile.h"
#include "segment-image.h"
#include "disjoint-set.h"

# include <math.h>

//subarna
#include "optical_flow\OpticalFlow.h"
#include "optical_flow\Image.h"

int main(int argc, char** argv) {
	if (argc != 11) {
		printf("%s c c_reg min sigma hie_num start_fr end_fr input output\n", argv[0]);
		printf("       c --> value for the threshold function in over-segmentation\n");
		printf("   c_reg --> value for the threshold function in hierarchical region segmentation\n");
		printf("     min --> enforced minimum supervoxel size\n");
		printf("   sigma --> variance of the Gaussian smoothing.\n");
		printf("   range --> number of frames as one subsequence (k in the paper)\n");
		printf(" hie_num --> desired number of hierarchy levels\n");
		printf(" start_fr --> starting frame index\n");
		printf(" end_fr --> end frame index\n");
		printf("   input --> input path of ppm video frames\n");
		printf("  output --> output path of segmentation results\n");
		return 1;
	}

	// Read Parameters
	float c = (float)atof(argv[1]);
	float c_reg = (float)atof(argv[2]);
	int min_size = atoi(argv[3]);
	float sigma = (float)atof(argv[4]);
	int range = atoi(argv[5]);
	int hie_num = atoi(argv[6]);
	int start_frame_num = atoi(argv[7]);
	int end_frame_num = atoi(argv[8]);
	char* input_path = argv[9];
	char* output_path = argv[10];
	if (c <= 0 || c_reg < 0 || min_size < 0 || sigma < 0 || hie_num < 0) {
		fprintf(stderr, "Unable to use the input parameters.");
		return 1;
	}

	//subarna: optical flow
	double alpha = 0.012;
	double ratio = 0.75;
	int minWidth = 20;
	int nOuterFPIterations = 7;
	int nInnerFPIterations = 1;
	int nSORIterations = 30;

	// count files in the input directory
	/*int frame_num = 0;
	struct dirent* pDirent;
	DIR* pDir;
	pDir = opendir(input_path);
	if (pDir != NULL) {
		while ((pDirent = readdir(pDir)) != NULL) {
			int len = strlen(pDirent->d_name);
			if (len >= 4) {
				if (strcmp(".ppm", &(pDirent->d_name[len - 4])) == 0)
					frame_num++;
			}
		}
	}
	//http://msdn.microsoft.com/en-us/library/windows/desktop/aa365200%28v=vs.85%29.aspx
	if (frame_num == 0) {
		fprintf(stderr, "Unable to find video frames at %s", input_path);
		return 1;
	}
	printf("Total number of frames in fold is %d\n", frame_num);

	// check if the range is right
	if (range > frame_num)
		range = frame_num;
	if (range < 1)
		range = 1;
	*/
	
	// check if the range is right
	if (range > (end_frame_num - start_frame_num))
		range = (end_frame_num - start_frame_num);
	if (range < 1)
		range = 1;

# ifdef USE_OPTICAL_FLOW
	if (range < 2)
	{
		range = 2;
		printf("range value should at least be 2 if motion values are to be used\n");
		printf("making the range value 2");
		if ( (end_frame_num - start_frame_num) < 2 )
		{
			printf("\n Not enough frames ... exiting\n\n\n\n");
			exit(1);
		}
	}
# endif

	// make the output directory
	struct stat st;
	int status = 0;
	char savepath[1024];
  	sprintf(savepath, "%s",output_path);
	if (stat(savepath, &st) != 0) {
		/* Directory does not exist */
		if (mkdir(savepath/*, S_IRWXU*/) != 0) {
			status = -1;
		}
	}
	for (int i = 0; i <= hie_num; i++) {
  		sprintf(savepath, "%s/%02d",output_path,i);
		if (stat(savepath, &st) != 0) {
			/* Directory does not exist */
			if (_mkdir(savepath/*, S_IRWXU*/) != 0) { //subarna: _mkdir
				status = -1;
			}
		}
	}
	if (status == -1) {
		fprintf(stderr,"Unable to create the output directories at %s",output_path);
		return 1;
	}
 
	// Initialize Parameters
	int last_clip = (end_frame_num - start_frame_num + 1) % range;
	int num_clip = (end_frame_num - start_frame_num + 1) / range;
	char filepath[1024];
	universe** u = new universe*[num_clip + 1];

	image<rgb>** input_first = new image<rgb>*[range];
	image<rgb>** input_middle = new image<rgb>*[range + 1];
	image<rgb>** input_last = new image<rgb>*[last_clip + 1];

	//subarna
# ifdef LUV
	DImage** input_first_LUV = new DImage*[range];
	DImage** input_middle_LUV = new DImage*[range+1];
	DImage** input_last_LUV = new DImage*[last_clip+1];
# endif

	// Time Recorder
	time_t Start_t, End_t;
	int time_task;
	Start_t = time(NULL);

	int height, width;
	
	// clip 1
	//calculate pair-wise optical flow
	//initialize memory and first frame
	printf("processing subsequence -- 0\n");
	for (int j = 0; j < range; j++) {
		sprintf(filepath, "%s/%05d.ppm", input_path, j + 1+ start_frame_num);
		input_first[j] = loadPPM(filepath);
		
		printf("load --> %s\n", filepath);

		if (j == 0 )
		{
			height = input_first[0]->height();
			width = input_first[0]->width();
		}

# ifdef LUV
		input_first_LUV[j] = new DImage (width, height,3);
		int m = 0;
		//convert to LUV and then apply bilateral filter
		for (int s = 0; s < height*width; s++)
		{     
			double x1,y1,z1;
			double x2,y2,z2;
			ccRGBtoXYZ(input_first[j]->data[s].r, input_first[j]->data[s].g, input_first[j]->data[s].b, &x1, &y1, &z1);
			ccXYZtoCIE_Luv(x1, y1, z1, &x2, &y2, &z2);
			input_first_LUV[j]->pData[m] = x2, 
			input_first_LUV[j]->pData[m+1] = y2, 
			input_first_LUV[j]->pData[m+2] = z2;
			m = m+3;
		}
		// pass input_first_LUV
# endif

	}

#ifdef USE_OPTICAL_FLOW
	DImage** D_vx_motion_first = new DImage*[range];
	DImage** D_vx_motion_middle = new DImage*[range + 1];
	DImage** D_vx_motion_last = new DImage*[last_clip + 1];

	DImage** D_vy_motion_first = new DImage*[range];
	DImage** D_vy_motion_middle = new DImage*[range + 1];
	DImage** D_vy_motion_last = new DImage*[last_clip + 1];
# endif

# ifdef USE_OPTICAL_FLOW
	DImage** D_input_first = new DImage*[range];
	DImage** D_input_middle = new DImage*[range + 1];
	DImage** D_input_last = new DImage*[last_clip + 1];

	for (int i = 0; i < range; i++ )
	{
		D_input_first[i] = new DImage(width,height,3);
		D_vx_motion_first[i] = new DImage (width,height);
		D_vy_motion_first[i] = new DImage (width,height);
	}
	image<rgb>* motion_buffer = NULL;
	DImage* D_motion_buffer = new DImage(width,height,3);

	//initialize
	int m = 0;
	int m1 = 0;
	for(int k1= 0; k1 < width; k1++)
	{
		for (int k2 = 0; k2 < height; k2++)
		{
			D_input_first[0]->pData[m1] = input_first[0]->data[m].r;
			D_input_first[0]->pData[m1+1] = input_first[0]->data[m].g;
			D_input_first[0]->pData[m1+2] = input_first[0]->data[m].b;
			m++;
			m1=m1+3;
		}
	}
	D_input_first[0]->im2double();

	for (int j = 0; j < range-1; j++) {
		// initialization for consecutive frames
		m = 0;	
		m1 = 0;
		for(int k1= 0; k1 < width; k1++)
		{
			for (int k2 = 0; k2 < height; k2++)
			{
				D_input_first[j+1]->pData[m1] = input_first[j+1]->data[m].r;
				D_input_first[j+1]->pData[m1+1] = input_first[j+1]->data[m].g;
				D_input_first[j+1]->pData[m1+2] = input_first[j+1]->data[m].b;
				m++;
				m1=m1+3;
			}
		}
		D_input_first[j+1]->im2double();

		OpticalFlow::Coarse2FineTwoFrames(*D_input_first[j], *D_input_first[j+1], &D_vx_motion_first[j], &D_vy_motion_first[j], alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations);
		
# ifdef DUMP_FLOW
		//debug: to be deleted
		image<rgb>* test_out = new image<rgb>(width, height);
		int m = 0;	
		for(int k1= 0; k1 < width; k1++)
		{
			for (int k2 = 0; k2 < height; k2++)
			{
				//test_out->data[m].r = ((int)(fabs(D_vx_motion_first[j]->pData[m])*30)>>4)<<4;
				//test_out->data[m].g = ((int)(fabs(D_vy_motion_first[j]->pData[m])*30)>>4)<<4;
				//test_out->data[m].b = 0 ; //test_out->data[m].r;

				test_out->data[m].r = 50 + 30*((int)(fabs(D_vx_motion_first[j]->pData[m])));
				test_out->data[m].g = 50 + 30*((int)(fabs(D_vy_motion_first[j]->pData[m])));
				test_out->data[m].b = 0 ; //test_out->data[m].r;
				
				m++;
			}
		}
		sprintf(filepath,"%s/motion/%05d.ppm",output_path, j + 1+ start_frame_num);
		savePPM(test_out, filepath); //"out_test/test1.ppm"
		delete test_out;
# endif

	}
	delete input_first[0];
	sprintf(filepath, "%s/%05d.ppm", input_path, range +1 + start_frame_num);
	input_first[0] = loadPPM(filepath);
	m = 0;	
	m1 = 0;
	for(int k1= 0; k1 < width; k1++)
	{
		for (int k2 = 0; k2 < height; k2++)
		{
			D_motion_buffer->pData[m1] = input_first[0]->data[m].r;
			D_motion_buffer->pData[m1+1] = input_first[0]->data[m].g;
			D_motion_buffer->pData[m1+2] = input_first[0]->data[m].b;
			m++;
			m1=m1+3;
		}
	}
	// dummy place holder
	D_motion_buffer->im2double();
	OpticalFlow::Coarse2FineTwoFrames(*D_input_first[range-2], *D_motion_buffer, &D_vx_motion_first[range-1], &D_vy_motion_first[range-1], alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations);
	
# if 0
	//copy content from D_motion_first[j-1] to D_motion_first[range-1]
	D_vx_motion_first[range-1]->copyData(*D_vx_motion_first[range-2]);
	D_vy_motion_first[range-1]->copyData(*D_vy_motion_first[range-2]);
# endif
	///////////////////////////////////////
# endif

# ifdef USE_OPTICAL_FLOW
	// frame index starts from 0
	# ifndef LUV
		u[0] = segment_image(output_path, input_first, D_vx_motion_first, D_vy_motion_first, 0, range - 1, c, c_reg, min_size,
				sigma, hie_num, NULL,start_frame_num);
	# else
		u[0] = segment_image_LUV(output_path, input_first_LUV, D_vx_motion_first, D_vy_motion_first, 0, range - 1, c, c_reg, min_size,
				sigma, hie_num, NULL,start_frame_num);
	# endif
# else
	// frame index starts from 0
	# ifndef LUV
		u[0] = segment_image(output_path, input_first,0, range - 1, c, c_reg, min_size,
				sigma, hie_num, NULL,start_frame_num);
	# else
		u[0] = segment_image_LUV(output_path, input_first_LUV, 0, range - 1, c, c_reg, min_size,
				sigma, hie_num, NULL,start_frame_num);
	# endif
# endif

	for (int j = 0; j < range; j++) {
		delete input_first[j];
# ifdef LUV
		delete input_first_LUV[j];
# endif
	}

# ifdef USE_OPTICAL_FLOW  //subarna
	for (int j = 0; j < range; j++) {	

		delete D_vx_motion_first[j];
		delete D_vy_motion_first[j];

		delete D_input_first[j];
	}
# endif

	// clip 2 -- last
	int ii;
	for (ii = 1; ii < num_clip; ii++) 
	{
		printf("processing subsequence -- %d\n", ii);

		for (int j = 0; j < range + 1; j++) 
		{
			sprintf(filepath, "%s/%05d.ppm", input_path, ii * range + j + start_frame_num);
			input_middle[j] = loadPPM(filepath);
			
			printf("load --> %s\n", filepath);

# ifdef LUV
			input_middle_LUV[j] = new DImage (width, height,3);
			int m = 0;
			//convert to LUV and then apply bilateral filter
			for (int s = 0; s < height*width; s++)
			{     
				double x1,y1,z1;
				double x2,y2,z2;
				ccRGBtoXYZ(input_middle[j]->data[s].r, input_middle[j]->data[s].g, input_middle[j]->data[s].b, &x1, &y1, &z1);
				ccXYZtoCIE_Luv(x1, y1, z1, &x2, &y2, &z2);
				input_middle_LUV[j]->pData[m] = x2, 
				input_middle_LUV[j]->pData[m+1] = y2, 
				input_middle_LUV[j]->pData[m+2] = z2;
				m = m+3;
			}
			//pass input_middle_LUV
# endif
		}

# ifdef USE_OPTICAL_FLOW
		for (int i = 0; i < range+1; i++ )
		{
			D_input_middle[i] = new DImage(width,height,3);
			D_vx_motion_middle[i] = new DImage(width,height);
			D_vy_motion_middle[i] = new DImage(width,height);
		}
		//initialize 
		m = 0;
		m1 = 0;
		for(int k1= 0; k1 < width; k1++)
		{
			for (int k2 = 0; k2 < height; k2++)
			{
				D_input_middle[0]->pData[m1] = input_middle[0]->data[m].r;
				D_input_middle[0]->pData[m1+1] = input_middle[0]->data[m].g;
				D_input_middle[0]->pData[m1+2] = input_middle[0]->data[m].b;
				m++;
				m1=m1+3;
			}
		}
		D_input_middle[0]->im2double();
		//calculate pair-wise optical flow 
		for (int j = 0; j < range; j++) {
			// initialization for consecutive frames
			m = 0;	
			m1 = 0;
			for(int k1= 0; k1 < width; k1++)
			{
				for (int k2 = 0; k2 < height; k2++)
				{
					D_input_middle[j+1]->pData[m1] = input_middle[j+1]->data[m].r;
					D_input_middle[j+1]->pData[m1+1] = input_middle[j+1]->data[m].g;
					D_input_middle[j+1]->pData[m1+2] = input_middle[j+1]->data[m].b;
					m++;
					m1 =m1+3;
				}
			}
			D_input_middle[j+1]->im2double();
			OpticalFlow::Coarse2FineTwoFrames(*D_input_middle[j], *D_input_middle[j+1], &D_vx_motion_middle[j], &D_vy_motion_middle[j], alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations);
		
# ifdef DUMP_FLOW
			//debug: to be deleted
			image<rgb>* test_out = new image<rgb>(width, height);
			int m = 0;	
			for(int k1= 0; k1 < width; k1++)
			{
				for (int k2 = 0; k2 < height; k2++)
				{
					//test_out->data[m].r = ((int)(fabs(D_vx_motion_first[j]->pData[m])*30)>>4)<<4;
					//test_out->data[m].g = ((int)(fabs(D_vy_motion_first[j]->pData[m])*30)>>4)<<4;
					//test_out->data[m].b = 0 ; //test_out->data[m].r;

					test_out->data[m].r = 50 + 40*((int)(fabs(D_vx_motion_middle[j]->pData[m])));
					test_out->data[m].g = 50 + 40*((int)(fabs(D_vx_motion_middle[j]->pData[m])));
					test_out->data[m].b = 0 ; //test_out->data[m].r;
				
					m++;
				}
			}
			sprintf(filepath,"%s/motion/%05d.ppm",output_path, i * range + j + start_frame_num);
			savePPM(test_out, filepath); //"out_test/test1.ppm"
			delete test_out;
# endif

		}
				
		// subarna: fix crash for specific frame number case
		if ( (ii == num_clip - 1 ) && (last_clip == 0) )
		{
			// for the last frame motion feature -- copy feature value from last frame, but with opposite sign
			for (int i = 0; i < height; i++)
			{
				for (int j = 0; j < width; j++)
				{
					D_vx_motion_middle[range]->pData[i*width + j] = -D_vx_motion_middle[range-1]->pData[i*width + j];
					D_vy_motion_middle[range]->pData[i*width + j] = -D_vy_motion_middle[range-1]->pData[i*width + j];
				}
			}
		}
		else
		{
			delete input_middle[0];
			sprintf(filepath, "%s/%05d.ppm", input_path, (ii+1) * range +1 + start_frame_num);
			input_middle[0] = loadPPM(filepath);		

			m = 0;	
			m1 = 0;
			for(int k1= 0; k1 < width; k1++)
			{
				for (int k2 = 0; k2 < height; k2++)
				{
					D_motion_buffer->pData[m1] = input_middle[0]->data[m].r;
					D_motion_buffer->pData[m1+1] = input_middle[0]->data[m].g;
					D_motion_buffer->pData[m1+2] = input_middle[0]->data[m].b;
					m++;
					m1=m1+3;
				}
			}
			
			// dummy place holder
			D_motion_buffer->im2double();
			OpticalFlow::Coarse2FineTwoFrames(*D_input_middle[range], *D_motion_buffer, &D_vx_motion_middle[range], &D_vy_motion_middle[range], alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations);
			/////////////////////////////////////
		}
# endif


# ifdef USE_OPTICAL_FLOW
	# ifndef LUV
		u[ii] = segment_image(output_path, input_middle, D_vx_motion_middle, D_vy_motion_middle, ii * range - 1,
				ii * range + range - 1, c, c_reg, min_size, sigma, hie_num,
				u[ii - 1], start_frame_num);
	# else
			u[ii] = segment_image_LUV(output_path, input_middle_LUV, D_vx_motion_middle, D_vy_motion_middle, ii * range - 1,
					ii * range + range - 1, c, c_reg, min_size, sigma, hie_num,
					u[ii - 1],start_frame_num);
	# endif
# else
	# ifndef LUV
		u[ii] = segment_image(output_path, input_middle, ii * range - 1,
				ii * range + range - 1, c, c_reg, min_size, sigma, hie_num,
				u[ii - 1], start_frame_num);
	# else
			u[ii] = segment_image_LUV(output_path, input_middle_LUV, ii * range - 1,
					ii* range + range - 1, c, c_reg, min_size, sigma, hie_num,
					u[ii - 1],start_frame_num);
	# endif
# endif

		if ( u[ii-1] ) delete u[ii - 1];

		////////////////////
		for (int j = 0; j < range + 1; j++) 
		{
			delete input_middle[j];			
# ifdef LUV
			delete input_middle_LUV[j];
# endif
		}

		# ifdef USE_OPTICAL_FLOW
		for (int j = 0; j < range + 1; j++) 
		{

			delete D_vx_motion_middle[j];
			delete D_vy_motion_middle[j];
			delete D_input_middle[j];
		}
		# endif
	}

	// clip last
	if (last_clip > 0) {
		printf("processing subsequence -- %d\n", num_clip);

		for (int j = 0; j < last_clip + 1; j++) 
		{
			sprintf(filepath, "%s/%05d.ppm", input_path, num_clip * range + j + start_frame_num);

			input_last[j] = loadPPM(filepath);

			printf("load --> %s\n", filepath);

# ifdef LUV
			input_last_LUV[j] = new DImage (width, height,3);
			int m=0;
			for (int s = 0; s < width*height; s++)
			{     
				double x1,y1,z1;
				double x2,y2,z2;
				ccRGBtoXYZ(input_last[j]->data[s].r, input_last[j]->data[s].g, input_last[j]->data[s].b, &x1, &y1, &z1);
				ccXYZtoCIE_Luv(x1, y1, z1, &x2, &y2, &z2);
				input_last_LUV[j]->pData[m] = x2, 
				input_last_LUV[j]->pData[m+1] = y2, 
				input_last_LUV[j]->pData[m+2] = z2;
				m = m+3;
			}
# endif
		}

# ifdef USE_OPTICAL_FLOW
		//subarna
		for (int i = 0; i < last_clip+1; i++ )
		{
			D_input_last[i] = new DImage(width,height,3);
			D_vx_motion_last[i] = new DImage(width,height);
			D_vy_motion_last[i] = new DImage(width,height);
		}
		//subarna
		//initialize
		m1 = 0;
		m = 0;
		for(int k1= 0; k1 < width; k1++)
		{
			for (int k2 = 0; k2 < height; k2++)
			{
				D_input_last[0]->pData[m1] = input_last[0]->data[m].r;
				D_input_last[0]->pData[m1+1] = input_last[0]->data[m].g;
				D_input_last[0]->pData[m1+2] = input_last[0]->data[m].b;//0
				m++;
				m1=m1+3;
			}
		}
		D_input_last[0]->im2double();
		//calculate pair-wise optical flow 
		for (int j = 0; j < last_clip; j++) {
			// initialization for consecutive frames
			m = 0;
			m1 = 0;
			for(int k1= 0; k1 < width; k1++)
			{
				for (int k2 = 0; k2 < height; k2++)
				{
					D_input_last[j+1]->pData[m1] = input_last[j+1]->data[m].r;
					D_input_last[j+1]->pData[m1+1] = input_last[j+1]->data[m].g;
					D_input_last[j+1]->pData[m1+2] = input_last[j+1]->data[m].b;
					m++;
					m1=m1+3;
				}
			}
			D_input_last[j+1]->im2double();
			OpticalFlow::Coarse2FineTwoFrames(*D_input_last[j], *D_input_last[j+1], &D_vx_motion_last[j], &D_vy_motion_last[j], alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations);

# ifdef DUMP_FLOW
			//debug: to be deleted
			image<rgb>* test_out = new image<rgb>(width, height);
			int m = 0;	
			for(int k1= 0; k1 < width; k1++)
			{
				for (int k2 = 0; k2 < height; k2++)
				{
					//test_out->data[m].r = ((int)(fabs(D_vx_motion_first[j]->pData[m])*30)>>4)<<4;
					//test_out->data[m].g = ((int)(fabs(D_vy_motion_first[j]->pData[m])*30)>>4)<<4;
					//test_out->data[m].b = 0 ; //test_out->data[m].r;

					test_out->data[m].r = 50 + 30*((int)(fabs(D_vx_motion_last[j]->pData[m])));
					test_out->data[m].g = 50 + 30*((int)(fabs(D_vx_motion_last[j]->pData[m])));
					test_out->data[m].b = 0 ; //test_out->data[m].r;
				
					m++;
				}
			}
			sprintf(filepath,"%s/motion/%05d.ppm",output_path, num_clip * range + j + start_frame_num);
			savePPM(test_out, filepath); //"out_test/test1.ppm"
			delete test_out;
# endif
		}

		// for the last frame motion feature -- copy feature value from last frame, but with opposite sign
		for (int i = 0; i < height; i++)
		{
			for (int j = 0; j < width; j++)
			{
				D_vx_motion_last[last_clip]->pData[i*width + j] = -D_vx_motion_last[last_clip-1]->pData[i*width + j];
				D_vy_motion_last[last_clip]->pData[i*width + j] = -D_vy_motion_last[last_clip-1]->pData[i*width + j];
			}
		}
		//////////////////////////////////////
# endif

# ifdef USE_OPTICAL_FLOW
	# ifndef LUV
			u[num_clip] = segment_image(output_path, input_last, D_vx_motion_last,D_vy_motion_last, num_clip * range - 1,
					num_clip * range + last_clip - 1, c, c_reg, min_size, sigma,
					hie_num, u[num_clip - 1],start_frame_num);
	# else
			u[num_clip] = segment_image_LUV(output_path, input_last_LUV, D_vx_motion_last,D_vy_motion_last, num_clip * range - 1,
					num_clip * range + last_clip - 1, c, c_reg, min_size, sigma,hie_num, u[num_clip - 1],start_frame_num);
	# endif
# else
	# ifndef LUV
			u[num_clip] = segment_image(output_path, input_last,num_clip * range - 1,
					num_clip * range + last_clip - 1, c, c_reg, min_size, sigma,
					hie_num, u[num_clip - 1],start_frame_num);
	# else
			u[num_clip] = segment_image_LUV(output_path, input_last_LUV, num_clip * range - 1,
					num_clip * range + last_clip - 1, c, c_reg, min_size, sigma,hie_num, u[num_clip - 1],start_frame_num);
	# endif
# endif


		if (u[num_clip - 1]) delete u[num_clip - 1];
		delete u[num_clip];

		for (int j = 0; j < last_clip + 1; j++) 
		{
			delete input_last[j];			
# ifdef LUV
			delete input_last_LUV[j];
# endif
		}

		//subarna
# ifdef USE_OPTICAL_FLOW
		for (int j = 0; j < last_clip + 1; j++) {

			delete D_vx_motion_last[j];
			delete D_vy_motion_last[j];
			delete D_input_last[j];
		}
		delete(D_motion_buffer);
# endif
	}

	//////////////////////////////////////
# ifdef LUV
	delete input_first_LUV;
	delete input_middle_LUV;
	delete input_last_LUV;
# endif

# ifdef USE_OPTICAL_FLOW
	delete D_vx_motion_first;
	delete D_vy_motion_first;
	delete D_vx_motion_middle;
	delete D_vy_motion_middle;
	delete D_vx_motion_last;
	delete D_vy_motion_last;
# endif
	/////////////////////////////////////

	delete[] u;

	// Time Recorder
	End_t = time(NULL);
	time_task = difftime(End_t, Start_t);
	std::ofstream myfile;
	char timefile[1024];
	sprintf(timefile, "%s/%s", output_path, "time.txt");
	myfile.open(timefile);
	myfile << time_task << endl;
	myfile.close();

	printf("Congratulations! It's done!\n");
	printf("Time_total = %d seconds\n", time_task);
	return 0;
}
