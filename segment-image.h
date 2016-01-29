/*
Original Code From:
Copyright (C) 2006 Pedro Felzenszwalb
Modifications (may have been made) Copyright (C) 2011,2012 
  Chenliang Xu, Jason Corso.

 Modifications (have been made) Copyright (C) 2013 
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

/* This is the main implementation of the streaming graph-based
 * hierarchical segmentation algorithm.
*/

#ifndef SEGMENT_IMAGE_H
#define SEGMENT_IMAGE_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
//#include <unistd.h>

#include "image.h"
#include "misc.h"
#include "filter.h"
#include "pnmfile.h"

#include "edges.h"
#include "segment-graph.h"
#include "disjoint-set.h"

#include "optical_flow/Image.h"

using namespace std;

/////////////////////////      added on April 8, 2013 : LUV and bilateral filter related     ///////////////////////////////////////////////////////
void BilateralFilter(image<rgb>* frame[], int num_frame, image<float> *smooth_r[], image<float> *smooth_g[], image<float> *smooth_b[],	int window, double sigma1, double sigma2, int height, int width);
void BilateralFilter_LUV(DImage* frame[], int num_frame, image<float> *smooth_r[], image<float> *smooth_g[], image<float> *smooth_b[],	int window, double sigma1, double sigma2, int height, int width);

static image<float> *BL_frame(image<float>* frame, int window, double sigma1, double sigma2, int height, int width);

//void bf_1dim(double* frame[], int window, double sigma1, double sigma2, int height, int width);
//void bf_3dim(double *frame[], int window, double sigma1, double sigma2, int height, int width);

void bf_3dim(image<float> **smooth_r, image<float> **smooth_g, image<float> **smooth_b, int w, double sigma1, double sigma2, int height, int width);

void compute_Gaussian_distance_weights(double* G[],int window, double sigma1);

void compute_Gaussian_intensity_weights(double* H[], double* I[], double* A, double sigma2, int a, int b, int size_r, int size_c, int i, int j, int height, int width);

double Bilateral_filter_response(double* F[], double* H[], double* I[], double* G[], int w, int size_r, int size_c, int a, int b, int i, int j);

#define JOIN_TH 1.6		//threshold for merging segments togheter
#define PURE_WHITE 255
#define MIN_REGION_BLOCK 16

#define YUV_MATCHING 1
#define LUV_MATCHING 0

/*Parameters for RGB -> XYZ -> LUV conversion*/

const double XYZ[3][3] = {	{  0.4125,  0.3576,  0.1804 },
							{  0.2125,  0.7154,  0.0721 },
							{  0.0193,  0.1192,  0.9502 }	};

enum {
  ccTA_2, ccTC_2, D50_2, D55_2, D65_2, D75_2, F2_2, F7_2, F11_2,
  ccTA_10, ccTC_10, D50_10, D55_10, D65_10, D75_10, F2_10, F7_10, F11_10,
  ccNTristimulus
};

enum { ccX, ccY, ccZ};

const double _ccTristimulusValues[ccNTristimulus][3] = 
{
      /* 2 degrees */

  /*ccTA_2*/{109.850, 100.000, 35.585},
  /*ccTC_2*/{98.074,  100.000, 118.232},
  /* D50 */ {96.422,  100.000, 82.521},
  /* D55 */ {95.682,  100.000, 92.149},
  /* D65 */ {95.047,  100.000, 108.883},
  /* D75 */ {94.972,  100.000, 122.638},
  /* F2 */  {99.187,  100.000, 67.395},
  /* F7 */  {95.044,  100.000, 108.755},
  /* F11 */ {100.966, 100.000, 64.370},

      /* 10 degrees */

  /* A */   {111.144, 100.000, 35.200},
  /* C */   {97.285,  100.000, 116.145},
  /* D50 */ {96.720,  100.000, 81.427},
  /* D55 */ {95.799,  100.000, 90.926},
  /* D65 */ {94.811,  100.000, 107.304},
  /* D75 */ {94.416,  100.000, 120.641},
  /* F2 */  {103.280, 100.000, 69.026},
  /* F7 */  {95.792,  100.000, 107.687},
  /* F11 */ {103.866, 100.000, 65.627}
};

/*Functions declaration*/

void ccRGBtoXYZ(unsigned char r, unsigned char g, unsigned char b, double* x, double* y, double* z);
void ccXYZtoCIE_Luv(double x, double y, double z, double* L, double* u, double* v);
void ccXYZtoCIE_L(double y, double* L);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/* random color */
rgb random_rgb() {
	rgb c;

	c.r = (uchar) rand();
	c.g = (uchar) rand();
	c.b = (uchar) rand();
	return c;
}

void smooth_images_LUV(DImage *im[], int num_frame, image<float> *smooth_r[],
		image<float> *smooth_g[], image<float> *smooth_b[], float sigma) {

	int width = im[0]->width();
	int height = im[0]->height();

	image<float>** r = new image<float>*[num_frame];
	image<float>** g = new image<float>*[num_frame];
	image<float>** b = new image<float>*[num_frame];
	for (int i = 0; i < num_frame; i++) {
		r[i] = new image<float>(width, height);
		g[i] = new image<float>(width, height);
		b[i] = new image<float>(width, height);
	}
	for (int i = 0; i < num_frame; i++) {
		int m = 0;
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				imRef(r[i], x, y) = im[i]->pData[m];
				imRef(g[i], x, y) = im[i]->pData[m+1];
				imRef(b[i], x, y) = im[i]->pData[m+2];
				m = m+3;
			}
		}
	}

	// smooth each color channel
	for (int i = 0; i < num_frame; i++) {
		smooth_r[i] = smooth(r[i], sigma);
		smooth_g[i] = smooth(g[i], sigma);
		smooth_b[i] = smooth(b[i], sigma);
	}
	for (int i = 0; i < num_frame; i++) {
		delete r[i];
		delete g[i];
		delete b[i];
	}
	delete[] r;
	delete[] g;
	delete[] b;
}

/* Gaussian Smoothing */
void smooth_images(image<rgb> *im[], int num_frame, image<float> *smooth_r[],
		image<float> *smooth_g[], image<float> *smooth_b[], float sigma) {

	int width = im[0]->width();
	int height = im[0]->height();

	image<float>** r = new image<float>*[num_frame];
	image<float>** g = new image<float>*[num_frame];
	image<float>** b = new image<float>*[num_frame];
	for (int i = 0; i < num_frame; i++) {
		r[i] = new image<float>(width, height);
		g[i] = new image<float>(width, height);
		b[i] = new image<float>(width, height);
	}
	for (int i = 0; i < num_frame; i++) {
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				imRef(r[i], x, y) = imRef(im[i], x, y).r;
				imRef(g[i], x, y) = imRef(im[i], x, y).g;
				imRef(b[i], x, y) = imRef(im[i], x, y).b;
			}
		}
	}

	// smooth each color channel
	for (int i = 0; i < num_frame; i++) {
		smooth_r[i] = smooth(r[i], sigma);
		smooth_g[i] = smooth(g[i], sigma);
		smooth_b[i] = smooth(b[i], sigma);
	}
	for (int i = 0; i < num_frame; i++) {
		delete r[i];
		delete g[i];
		delete b[i];
	}
	delete[] r;
	delete[] g;
	delete[] b;
}

/* Save Output */
void generate_output(bool first, char *path, int frame_id_start,
		int frame_id_end, int width, int height, universe *mess,
		int num_vertices, int num_vertices_v, int level_total, int start_frame_num) {

	int num_save = frame_id_end - frame_id_start;
	int save_start = frame_id_start + 2 + start_frame_num;
	if (first == true) {
		num_save = frame_id_end - frame_id_start + 1;
		save_start = frame_id_start + 1 + start_frame_num; //subarna
	}

# ifdef UNIQUE_COLORS
	// subarna
	char ***unique_color_flag = NULL;
	unique_color_flag = (char***)malloc(256*sizeof(char**));
	for (int r = 0; r < 256; r++)
	{
		unique_color_flag[r] = (char**)malloc(256*sizeof(char*));
		for (int g = 0; g < 256; g++)
		{
			unique_color_flag[r][g] = (char*)malloc(256*sizeof(char));
			for (int b = 0; b < 256; b++)
			{
				unique_color_flag[r][g][b] = 1;
			}
		}
	}
# endif

	char savepath[1024];
	image<rgb>** output = new image<rgb>*[num_save];
	rgb* colors = new rgb[num_vertices];
	for (int i = 0; i < num_vertices; i++)
	{
		colors[i] = random_rgb();

# ifdef UNIQUE_COLORS
		//subarna: enforce unique colors
		while (1)
		{
			if(unique_color_flag[colors[i].r][colors[i].g][colors[i].b] == 1)
			{
				unique_color_flag[colors[i].r][colors[i].g][colors[i].b] = 0;
				break;
			}
			colors[i] = random_rgb();
		}
# endif
	}

# ifdef UNIQUE_COLORS
	for (int r = 0; r < 256; r++)
	{
		for (int g = 0; g < 256; g++)
		{
			if(!unique_color_flag[r][g])
				free(unique_color_flag[r][g]);
		}
		if(!unique_color_flag[r])
			free(unique_color_flag[r]);
	}
	if(!unique_color_flag)
		free(unique_color_flag);
# endif

	rgb blank;
	blank.r = blank.g = blank.b = 0;

	// write out the ppm files.
	for (int k = 0; k <= level_total; k++) {
		for (int i = 0; i < num_save; i++) {
			// add frame index 1 to save
			sprintf(savepath, "%s/%02d/%05d.ppm", path, k, save_start + i);
			output[i] = new image<rgb>(width, height);
			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {
					int this_id = num_vertices_v + y * width + x
							+ i * (width * height);
					int comp = mess->find_in_level(this_id, k);
					if (mess->get_color_output(comp, k) == blank) {
						// color the segment if the segment is not colored yet
						imRef(output[i], x, y) = colors[comp];
						mess->set_color_output(this_id, colors[comp], k);
					} else {
						// use the segment color if it is already colored
						imRef(output[i], x, y) = mess->get_color_output(comp,
								k);
					}
				}
			}
			savePPM(output[i], savepath);
		}
		for (int i = 0; i < num_save; i++)
			delete output[i];
	}
	delete[] colors;
	delete[] output;

}

# ifdef USE_OPTICAL_FLOW
/* main operation steps in one iteration */
universe *segment_image(char *path, image<rgb> *im[], DImage *vx_motion[], DImage *vy_motion[], int frame_id_start,
		int frame_id_end, float c, float c_reg, int min_size, float sigma,
		int hie_num, universe *v, int start_frame_num)
# else
universe *segment_image(char *path, image<rgb> *im[],int frame_id_start,
		int frame_id_end, float c, float c_reg, int min_size, float sigma,
		int hie_num, universe *v, int start_frame_num)
# endif 
{
	// step 1 -- Get information
	// ----- width, heigh, frame number
	int width = im[0]->width();
	int height = im[0]->height();
	int num_frame = frame_id_end - frame_id_start + 1;

# ifdef BLF
	// for bi-lateral filtering
	double sigma1 = 1.5, sigma2 = 0.1;
	int window = 3;
	sigma1 = 1.2; //test
# endif

	// ----- first or not
	bool first;
	if (frame_id_start == 0)
		first = true;
	else
		first = false;
	
	// ----- node number
	int num_vertices = num_frame * width * height;
	int num_vertices_v = width * height;
	// ----- edge number
	int num_edges_plane = (width - 1) * (height - 1) * 2 + width * (height - 1)
			+ (width - 1) * height;
	int num_edges_layer = (width - 2) * (height - 2) * 9 + (width - 2) * 2 * 6
			+ (height - 2) * 2 * 6 + 4 * 4;
	int num_edges = num_edges_plane * (num_frame - 1)
			+ num_edges_layer * (num_frame - 1);
	
	// if it's first
	if (first == true) {
		num_edges = num_edges_plane * num_frame
				+ num_edges_layer * (num_frame - 1);
		num_vertices_v = 0;
	}
	// ----- hierarchy setup
	vector<vector<edge>*> edges_region;
	edges_region.resize(hie_num + 1);
    // ------------------------------------------------------------------

	// step 2 -- smooth images
	image<float>** smooth_r = new image<float>*[num_frame];
	image<float>** smooth_g = new image<float>*[num_frame];
	image<float>** smooth_b = new image<float>*[num_frame];

# ifndef BLF
	smooth_images(im, num_frame, smooth_r, smooth_g, smooth_b, sigma);
# else
	/* Applying Bilateral Filtering for edge preservation and surfaces smoothing */
	BilateralFilter(im, num_frame, smooth_r, smooth_g, smooth_b, window, sigma1, sigma2, height, width);
# endif
	// ------------------------------------------------------------------

	// step 3 -- build edges
	printf("start build edges\n");
	edge* edges = new edge[num_edges];


# if defined(USE_OPTICAL_FLOW) && defined(FLOW_EDGE)
	initialize_edges(edges, first, num_frame, width, height, smooth_r, smooth_g,
			smooth_b,vx_motion,vy_motion);
# else
	initialize_edges(edges, first, num_frame, width, height, smooth_r, smooth_g,
			smooth_b);
# endif

//# ifndef FLOW_EDGE 
//	initialize_edges(edges, first, num_frame, width, height, smooth_r, smooth_g,
//			smooth_b);
//# else
//	initialize_edges(edges, first, num_frame, width, height, smooth_r, smooth_g,
//			smooth_b,vx_motion,vy_motion);
//# endif

	printf("end build edges\n");
	// ------------------------------------------------------------------

	// step 4 -- build nodes
	printf("start build nodes\n");

# if defined(USE_OPTICAL_FLOW) && defined(MOTION_FEATURE)
	universe* mess = new universe(frame_id_start, frame_id_end, width, height,
			*v, num_vertices_v, smooth_r, smooth_g, smooth_b, vx_motion, vy_motion, hie_num);
# else
	universe* mess = new universe(frame_id_start, frame_id_end, width, height,
			*v, num_vertices_v, smooth_r, smooth_g, smooth_b, hie_num);
# endif


//# ifndef MOTION_FEATURE 
//	universe* mess = new universe(frame_id_start, frame_id_end, width, height,
//			*v, num_vertices_v, smooth_r, smooth_g, smooth_b, hie_num);
//# else
//	universe* mess = new universe(frame_id_start, frame_id_end, width, height,
//			*v, num_vertices_v, smooth_r, smooth_g, smooth_b, vx_motion, vy_motion, hie_num);
//# endif

	printf("end build nodes\n");
	// ------------------------------------------------------------------

	// step 5 -- over-segmentation
	printf("start over-segmentation\n");
	edges_region[0] = new vector<edge>();
	segment_graph(mess, edges_region[0], edges, num_edges, c, 0);
	// optional merging small components
	for (int i = 0; i < num_edges; i++) {
		int a = mess->find_in_level(edges[i].a, 0);
		int b = mess->find_in_level(edges[i].b, 0);
		if ((a != b)
				&& ((mess->get_size(a) < min_size)
						|| (mess->get_size(b) < min_size)))
			mess->join_noise(a, b, 0);
	}
	printf("end over-segmentation\n");
	// ------------------------------------------------------------------

	// step 6 -- hierarchical segmentation
	for (int i = 0; i < hie_num; i++) {
		printf("level = %d\n", i);
		
		// subarna
		// incremental in each hierarhcy
		min_size = min_size * 1.2;


		
# if defined(USE_OPTICAL_FLOW) && defined(HIGHEST_LEVEL_MOTION_ONLY)
		/*if ( i == hie_num -1)
			min_size = min_size*6;*/
# endif

//		// subarna: debug: to be deleted
//# ifdef HIGHEST_LEVEL_MOTION_ONLY
//		/*if ( i == hie_num -1)
//			min_size = min_size*6;*/
//# endif

		printf("start update\n");
		mess->update(i);
		printf("end update\n");

		printf("start fill edge weight\n");

# if defined(USE_OPTICAL_FLOW) && defined(HIGHEST_LEVEL_MOTION_ONLY)
		fill_edge_weight(*edges_region[i], mess, i,hie_num-1);
# else
		fill_edge_weight(*edges_region[i], mess, i);
# endif

//# ifndef HIGHEST_LEVEL_MOTION_ONLY
//		fill_edge_weight(*edges_region[i], mess, i);
//# else
//		fill_edge_weight(*edges_region[i], mess, i,hie_num-1);
//# endif

		printf("end fill edge weight\n");

		printf("start segment graph region\n");
		edges_region[i + 1] = new vector<edge>();

# if defined(USE_OPTICAL_FLOW) && defined(HIGHEST_LEVEL_MOTION_ONLY)
		segment_graph_region(mess, edges_region[i + 1], edges_region[i], c_reg, i + 1, hie_num);
# else
		segment_graph_region(mess, edges_region[i + 1], edges_region[i], c_reg, i + 1);
# endif

//# ifndef HIGHEST_LEVEL_MOTION_ONLY
//		segment_graph_region(mess, edges_region[i + 1], edges_region[i], c_reg, i + 1);
//# else
//		segment_graph_region(mess, edges_region[i + 1], edges_region[i], c_reg, i + 1, hie_num);
//# endif
		printf("end segment graph region\n");

		printf("start merging min_size\n");
		for (int it = 0; it < (int) edges_region[i]->size(); it++) {
			int a = mess->find_in_level((*edges_region[i])[it].a, i + 1);
			int b = mess->find_in_level((*edges_region[i])[it].b, i + 1);
			if ((a != b) && ((mess->get_size(a) < min_size) || (mess->get_size(b) < min_size)))
				mess->join_noise(a, b, i + 1);
		}
		printf("end merging min_size\n");

		// incremental in each hierarchy
		c_reg = c_reg * 1.4;

# if defined(USE_OPTICAL_FLOW) && defined(HIGHEST_LEVEL_MOTION_ONLY)
		/*if ( i == hie_num - 2)
			c_reg = c_reg*7;*/
# endif

//	    // subarna: debug: to be deleted
//# ifdef HIGHEST_LEVEL_MOTION_ONLY
//		/*if ( i == hie_num - 2)
//			c_reg = c_reg*7;*/
//# endif

		delete edges_region[i];
	}
	delete edges_region[hie_num];
	// ------------------------------------------------------------------

	// step 8 -- generate output
	printf("start output\n");
	generate_output(first, path, frame_id_start, frame_id_end, width, height,
			mess, num_vertices, num_vertices_v, hie_num, start_frame_num);
	printf("end output\n");
	// ------------------------------------------------------------------

	// step 9 -- clear everything and return u contains last frame info
	universe *u = new universe(*mess, num_vertices, width * height);
	delete mess;
	delete[] edges;
	for (int i = 0; i < num_frame; i++) {
		delete smooth_r[i];
		delete smooth_g[i];
		delete smooth_b[i];
	}
	delete[] smooth_r;
	delete[] smooth_g;
	delete[] smooth_b;
	return u;

}

# ifdef USE_OPTICAL_FLOW
/* main operation steps in one iteration */
universe *segment_image_LUV(char *path, DImage *im[], DImage *vx_motion[], DImage *vy_motion[], int frame_id_start,
		int frame_id_end, float c, float c_reg, int min_size, float sigma,
		int hie_num, universe *v, int start_frame_num) 
# else
universe *segment_image_LUV(char *path, DImage *im[], int frame_id_start,
		int frame_id_end, float c, float c_reg, int min_size, float sigma,
		int hie_num, universe *v, int start_frame_num) 
# endif 
{
	// step 1 -- Get information
	// ----- width, heigh, frame number
	int width = im[0]->width();
	int height = im[0]->height();
	int num_frame = frame_id_end - frame_id_start + 1;

# ifdef BLF
	// for bi-lateral filtering
	double sigma1 = 1.5, sigma2 = 0.1;
	int window = 3;
	sigma1 = 0.8; //test
# endif

	// ----- first or not
	bool first;
	if (frame_id_start == 0)
		first = true;
	else
		first = false;
	// ----- node number
	int num_vertices = num_frame * width * height;
	int num_vertices_v = width * height;
	// ----- edge number
	int num_edges_plane = (width - 1) * (height - 1) * 2 + width * (height - 1)
			+ (width - 1) * height;
	int num_edges_layer = (width - 2) * (height - 2) * 9 + (width - 2) * 2 * 6
			+ (height - 2) * 2 * 6 + 4 * 4;
	int num_edges = num_edges_plane * (num_frame - 1)
			+ num_edges_layer * (num_frame - 1);
	// if it's first
	if (first == true) {
		num_edges = num_edges_plane * num_frame
				+ num_edges_layer * (num_frame - 1);
		num_vertices_v = 0;
	}
	// ----- hierarchy setup
	vector<vector<edge>*> edges_region;
	edges_region.resize(hie_num + 1);
    // ------------------------------------------------------------------

	// step 2 -- smooth images
	image<float>** smooth_r = new image<float>*[num_frame];
	image<float>** smooth_g = new image<float>*[num_frame];
	image<float>** smooth_b = new image<float>*[num_frame];

# if (!defined(BLF)) 
	smooth_images_LUV(im, num_frame, smooth_r, smooth_g, smooth_b, sigma);
# else
	/* Applying Bilateral Filtering for edge preservation and surfaces smoothing */	
	BilateralFilter_LUV(im, num_frame, smooth_r, smooth_g, smooth_b, window, sigma1, sigma2, height, width);
# endif
	// ------------------------------------------------------------------

	// step 3 -- build edges
	printf("start build edges\n");
	edge* edges = new edge[num_edges];

# if defined(USE_OPTICAL_FLOW) && defined(FLOW_EDGE)
	initialize_edges(edges, first, num_frame, width, height, smooth_r, smooth_g,
			smooth_b,vx_motion, vy_motion);
# else
	initialize_edges(edges, first, num_frame, width, height, smooth_r, smooth_g,
			smooth_b);
# endif

//# ifndef FLOW_EDGE
//	initialize_edges(edges, first, num_frame, width, height, smooth_r, smooth_g,
//			smooth_b);
//# else
//	initialize_edges(edges, first, num_frame, width, height, smooth_r, smooth_g,
//			smooth_b,vx_motion, vy_motion);
//# endif

	printf("end build edges\n");
	// ------------------------------------------------------------------

	// step 4 -- build nodes
	printf("start build nodes\n");

# if defined(USE_OPTICAL_FLOW) && defined(MOTION_FEATURE)
	universe* mess = new universe(frame_id_start, frame_id_end, width, height,
			*v, num_vertices_v, smooth_r, smooth_g, smooth_b, vx_motion, vy_motion, hie_num);
# else
	universe* mess = new universe(frame_id_start, frame_id_end, width, height,
			*v, num_vertices_v, smooth_r, smooth_g, smooth_b, hie_num);
# endif


//# ifndef MOTION_FEATURE
//	universe* mess = new universe(frame_id_start, frame_id_end, width, height,
//			*v, num_vertices_v, smooth_r, smooth_g, smooth_b, hie_num);
//# else
//	universe* mess = new universe(frame_id_start, frame_id_end, width, height,
//			*v, num_vertices_v, smooth_r, smooth_g, smooth_b, vx_motion, vy_motion, hie_num);
//# endif
	printf("end build nodes\n");
	// ------------------------------------------------------------------

	// step 5 -- over-segmentation
	printf("start over-segmentation\n");
	edges_region[0] = new vector<edge>();
	segment_graph(mess, edges_region[0], edges, num_edges, c, 0);
	// optional merging small components
	for (int i = 0; i < num_edges; i++) {
		int a = mess->find_in_level(edges[i].a, 0);
		int b = mess->find_in_level(edges[i].b, 0);
		if ((a != b)
				&& ((mess->get_size(a) < min_size)
						|| (mess->get_size(b) < min_size)))
			mess->join_noise(a, b, 0);
	}
	printf("end over-segmentation\n");
	// ------------------------------------------------------------------

	// step 6 -- hierarchical segmentation
	for (int i = 0; i < hie_num; i++) {
		printf("level = %d\n", i);
		
		// incremental in each hierarhcy
		min_size = min_size * 1.2;

# ifdef USE_OPTICAL_FLOW
# ifdef HIGHEST_LEVEL_MOTION_ONLY
		if ( i == hie_num -1)
			min_size = min_size*6;
# endif
# endif

		printf("start update\n");
		mess->update(i);
		printf("end update\n");

		printf("start fill edge weight\n");

# if defined(USE_OPTICAL_FLOW) && defined(HIGHEST_LEVEL_MOTION_ONLY)
		fill_edge_weight(*edges_region[i], mess, i,hie_num-1);
# else
		fill_edge_weight(*edges_region[i], mess, i);
# endif

//
//# ifndef HIGHEST_LEVEL_MOTION_ONLY
//		fill_edge_weight(*edges_region[i], mess, i);
//# else
//		fill_edge_weight(*edges_region[i], mess, i,hie_num-1);
//# endif
		
		printf("end fill edge weight\n");

		printf("start segment graph region\n");
		edges_region[i + 1] = new vector<edge>();

# if defined(USE_OPTICAL_FLOW) && defined(HIGHEST_LEVEL_MOTION_ONLY)
		segment_graph_region(mess, edges_region[i + 1], edges_region[i], c_reg, i + 1, hie_num);
# else
		segment_graph_region(mess, edges_region[i + 1], edges_region[i], c_reg, i + 1);
# endif

//# ifndef HIGHEST_LEVEL_MOTION_ONLY
//		segment_graph_region(mess, edges_region[i + 1], edges_region[i], c_reg, i + 1);
//# else
//		segment_graph_region(mess, edges_region[i + 1], edges_region[i], c_reg, i + 1, hie_num);
//# endif

		printf("end segment graph region\n");

		printf("start merging min_size\n");
		for (int it = 0; it < (int) edges_region[i]->size(); it++) {
			int a = mess->find_in_level((*edges_region[i])[it].a, i + 1);
			int b = mess->find_in_level((*edges_region[i])[it].b, i + 1);
			if ((a != b) && ((mess->get_size(a) < min_size) || (mess->get_size(b) < min_size)))
				mess->join_noise(a, b, i + 1);
		}
		printf("end merging min_size\n");

		// incremental in each hierarchy
		c_reg = c_reg * 1.4;
		
		// subarna: debug: to be deleted
# if defined(USE_OPTICAL_FLOW) && defined(HIGHEST_LEVEL_MOTION_ONLY)
		if ( i == hie_num - 2)
			c_reg = c_reg*7;
# endif

//# ifdef HIGHEST_LEVEL_MOTION_ONLY
//		if ( i == hie_num - 2)
//			c_reg = c_reg*7;
//# endif

		delete edges_region[i];
	}
	delete edges_region[hie_num];
	// ------------------------------------------------------------------

	// step 8 -- generate output
	printf("start output\n");
	generate_output(first, path, frame_id_start, frame_id_end, width, height,
			mess, num_vertices, num_vertices_v, hie_num, start_frame_num);
	printf("end output\n");
	// ------------------------------------------------------------------

	// step 9 -- clear everything and return u contains last frame info
	universe *u = new universe(*mess, num_vertices, width * height);
	delete mess;
	delete[] edges;
	for (int i = 0; i < num_frame; i++) {
		delete smooth_r[i];
		delete smooth_g[i];
		delete smooth_b[i];
	}
	delete[] smooth_r;
	delete[] smooth_g;
	delete[] smooth_b;
	return u;

}


/////////////////////////      added on April 8, 2013 : LUV and bilateral filter related     ///////////////////////////////////////////////////////
void BilateralFilter(image<rgb> *im[], int num_frame, image<float> *smooth_r[],
		image<float> *smooth_g[], image<float> *smooth_b[],
		int window, double sigma1, double sigma2, int height, int width)
{
    //int i;    

    //Luv components normalization to 0:1
    double add_u, add_v, div_L, div_u, div_v;
    
    /*Parameters found previously, converting RGB(0:255,0:255,0:255) into L*u*v*/

    div_L = 100.0;                 //(max L) value = L range
    add_u = 83.079752;             //abs(min u) value (translation parameter)
    div_u = 175.053036 + add_u;    //(max u) + abs(min u) value = u range
    add_v = 134.116076;            //abs(min v) value  (translation parameter)
    div_v = 107.401365 + add_v;    //(max v) + abs(min v) value = v range

	for (int i = 0; i < num_frame; i++) {
		smooth_r[i] = new image<float>(width, height);
		smooth_g[i] = new image<float>(width, height);
		smooth_b[i] = new image<float>(width, height);
	}

	// smooth each color channel
	/*Parameters found previously, converting RGB(0:255,0:255,0:255) into L*u*v*/
	for (int i = 0; i < num_frame; i++) 
	{
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				imRef(smooth_r[i], x, y) = imRef(im[i], x, y).r/div_L;
				imRef(smooth_b[i], x, y) = (imRef(im[i], x, y).g + add_u)/div_u;
				imRef(smooth_b[i], x, y) = (imRef(im[i], x, y).b + add_v)/div_v;
			}
		}
		
		//filtering L*u*v components
		// subarna
		bf_3dim(&smooth_r[i], &smooth_g[i], &smooth_b[i], window, sigma1, sigma2, height, width);       

		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				imRef(smooth_r[i], x, y) = imRef(smooth_r[i], x, y)*100;
				imRef(smooth_g[i], x, y) = imRef(smooth_g[i], x, y)* div_u - add_u;
				imRef(smooth_b[i], x, y) = imRef(smooth_b[i], x, y)*div_v - add_v;
			}
		}
	}	
			
    return; 
}

void BilateralFilter_LUV(DImage *im[], int num_frame, image<float> *smooth_r[],
		image<float> *smooth_g[], image<float> *smooth_b[],
		int window, double sigma1, double sigma2, int height, int width)
{   
    //Luv components normalization to 0:1
    double add_u, add_v, div_L, div_u, div_v;
	int m;
	div_L = 100.0;                 //(max L) value = L range
    add_u = 83.079752;             //abs(min u) value (translation parameter)
    div_u = 175.053036 + add_u;    //(max u) + abs(min u) value = u range
    add_v = 134.116076;            //abs(min v) value  (translation parameter)
    div_v = 107.401365 + add_v;    //(max v) + abs(min v) value = v range
    
	for (int i = 0; i < num_frame; i++) {
		smooth_r[i] = new image<float>(width, height);
		smooth_g[i] = new image<float>(width, height);
		smooth_b[i] = new image<float>(width, height);
	}
	
	// smooth each color channel
	/*Parameters found previously, converting RGB(0:255,0:255,0:255) into L*u*v*/
	for (int i = 0; i < num_frame; i++) 
	{
		m = 0;
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				imRef(smooth_r[i], x, y) = (im[i]->pData[m])/div_L;
				imRef(smooth_b[i], x, y) = (im[i]->pData[m+1] + add_u)/div_u;
				imRef(smooth_b[i], x, y) = (im[i]->pData[m+2] + add_v)/div_v;
			
				m = m+3;
			}
		}
		
		//filtering L*u*v components
		bf_3dim(&smooth_r[i], &smooth_g[i], &smooth_b[i], window, sigma1, sigma2, height, width);       

		m = 0;
		for (int y = 0; y < height; y++) 
		{
			for (int x = 0; x < width; x++) 
			{
				imRef(smooth_r[i], x, y) = imRef(smooth_r[i], x, y)*100;
				imRef(smooth_g[i], x, y) = imRef(smooth_g[i], x, y)* div_u - add_u;
				imRef(smooth_b[i], x, y) = imRef(smooth_b[i], x, y)*div_v - add_v;

				/*im[i]->pData[m] = imRef(smooth_r[i], x, y);
				im[i]->pData[m+1] = imRef(smooth_g[i], x, y);
				im[i]->pData[m+2] = imRef(smooth_b[i], x, y);*/

				m=m+3;
			}
		}
	}

     return; 
}

void bf_3dim(image<float> **smooth_r, image<float> **smooth_g, image<float> **smooth_b, int w, double sigma1, double sigma2, int height, int width)
{

     double **G,**I0,**I1,**I2,**H0,**H1,**H2,**H,**F;
     double *tmp[3];
     int i,j,size_r,size_c,a,b,d,e,t,r,c;
     tmp[0] = (double*) calloc (height*width, sizeof(double));
     tmp[1] = (double*) calloc (height*width, sizeof(double));
     tmp[2] = (double*) calloc (height*width, sizeof(double));
     G = (double**) calloc (((w*2)+1), sizeof(double*));

     for (i = 0; i < ((w*2)+1); i++)
     {
         G[i] = (double*) calloc (((w*2)+1), sizeof(double));
     }

	 double *frame[3];
	 frame[0] = (double*) calloc (height*width, sizeof(double));
     frame[1] = (double*) calloc (height*width, sizeof(double));
     frame[2] = (double*) calloc (height*width, sizeof(double));

	 for (i = 0; i < height; i++)
	 {
         for (j = 0; j< width; j++)
         {    
			 frame[0][i*width+j] = (*smooth_r)->data[i*width + j]; //imRef(smooth_r,j,i).r;
			 frame[1][i*width+j] = (*smooth_g)->data[i*width + j]; //imRef(smooth_g,j,i).g;
			 frame[2][i*width+j] = (*smooth_b)->data[i*width + j]; //imRef(smooth_b,j,i).b;
		 }
	 }

     compute_Gaussian_distance_weights(G, w, sigma1);

     for (i = 0; i < height; i++)
         for (j = 0; j< width; j++)
         {           
             if ((i - w) < 0)
                a = 0;
             else a = i - w;

             if ((i + w) > (height - 1))
                e = (height - 1);
             else e = i + w;
			 
             if ((j - w) < 0)
                b = 0;
             else b = j - w;

             if ((j + w) > (width - 1))
                d = (width - 1);
             else d = j + w;
			 

             size_r = (e - a) + 1;
             size_c = (d - b) + 1;


             I0 = (double**) calloc((size_r),sizeof(double*));
             I1 = (double**) calloc((size_r),sizeof(double*));
             I2 = (double**) calloc((size_r),sizeof(double*));
             H0 = (double**) calloc((size_r),sizeof(double*));
             H1 = (double**) calloc((size_r),sizeof(double*));
             H2 = (double**) calloc((size_r),sizeof(double*));
             H = (double**) calloc((size_r),sizeof(double*));
             F = (double**) calloc((size_r),sizeof(double*));

             for (t = 0; t < size_r; t++)
             {
                 I0[t] = (double*) calloc((size_c),sizeof(double));
                 I1[t] = (double*) calloc((size_c),sizeof(double));
                 I2[t] = (double*) calloc((size_c),sizeof(double));
                 H0[t] = (double*) calloc((size_c),sizeof(double));
                 H1[t] = (double*) calloc((size_c),sizeof(double));
                 H2[t] = (double*) calloc((size_c),sizeof(double));
                 H[t] = (double*) calloc((size_c),sizeof(double));
                 F[t] = (double*) calloc((size_c),sizeof(double));
             }
			 
             compute_Gaussian_intensity_weights(H0, I0, frame[0], sigma2, a, b, size_r, size_c, i, j, height, width);
             compute_Gaussian_intensity_weights(H1, I1, frame[1], sigma2, a, b, size_r, size_c, i, j, height, width);
             compute_Gaussian_intensity_weights(H2, I2, frame[2], sigma2, a, b, size_r, size_c, i, j, height, width);

             //compute complessive H
             for (r = 0; r < size_r; r++)
                 for (c = 0; c < size_c; c++)
                     H[r][c] = exp(-((H0[r][c])*(H0[r][c])+(H1[r][c])*(H1[r][c])+(H2[r][c])*(H2[r][c]))/(2*sigma2*sigma2));

             for (t = 0; t < size_r; t++)
             {        
                 free(H0[t]); H0[t] = NULL;
                 free(H1[t]); H1[t] = NULL;
                 free(H2[t]); H2[t] = NULL;
             }

             free(H0); H0 = NULL;
             free(H1); H1 = NULL;
             free(H2); H2 = NULL;
			            
             tmp[0][((i*width)+j)] = Bilateral_filter_response(F, H, I0, G, w, size_r, size_c, a, b, i, j);
             tmp[1][((i*width)+j)] = Bilateral_filter_response(F, H, I1, G, w, size_r, size_c, a, b, i, j);
             tmp[2][((i*width)+j)] = Bilateral_filter_response(F, H, I2, G, w, size_r, size_c, a, b, i, j);             

             for (t = 0; t < size_r; t++)
             {
                 free(I0[t]); I0[t] = NULL;
                 free(I1[t]); I1[t] = NULL;
                 free(I2[t]); I2[t] = NULL;
                 free(H[t]); H[t] = NULL;
                 free(F[t]); F[t] = NULL;
             }

             free(I0); I0 = NULL;
             free(I1); I1 = NULL;
             free(I2); I2 = NULL;
             free(H);  H = NULL;
             free(F);  F = NULL;
         }

     for (i = 0; i < (height); i++)
     {
         for (j = 0; j < (width); j++)
         {
             frame[0][i*width+j] = tmp[0][i*width+j];
             frame[1][i*width+j] = tmp[1][i*width+j];
             frame[2][i*width+j] = tmp[2][i*width+j];

			 (*smooth_r)->data[i*width+j] = frame[0][i*width+j];
			 (*smooth_g)->data[i*width+j] = frame[1][i*width+j];
			 (*smooth_b)->data[i*width+j] = frame[2][i*width+j];
         }
     }    

     for (i = 0; i < ((w*2)+1); i++)
     {
         free(G[i]); G[i] = NULL;
     }

     free(G); G = NULL;
     free(tmp[0]); tmp[0] = NULL;
     free(tmp[1]); tmp[1] = NULL;
     free(tmp[2]); tmp[2] = NULL;    

     return;
}



void compute_Gaussian_distance_weights(double* G[],int window, double sigma1)
{
	 double **X,**Y;
     int i,j;
     double a,b,c;
     
     X = (double**) calloc (((window*2)+1), sizeof(double*));
     Y = (double**) calloc (((window*2)+1), sizeof(double*));

     for (i = 0; i < ((window*2)+1); i++)
     {
         X[i] = (double*) calloc (((window*2)+1), sizeof(double));
         Y[i] = (double*) calloc (((window*2)+1), sizeof(double));
     }
	 
     for (i = 0; i < ((window*2)+1); i++)
         for (j = 0; j < ((window*2)+1); j++)
         {
            X[i][j] =  (double)(j - window);
            Y[i][j] =  (double)(i - window);
         }
		         

     for (i = 0; i < ((window*2)+1); i++)
         for (j = 0; j < ((window*2)+1); j++)
         {
             a = ((X[i][j])*(X[i][j]));
             b = ((Y[i][j])*(Y[i][j]));
             c = (sigma1*sigma1);
             G[i][j] = exp(-(a+b)/(2*c));
         }


     for (i = 0; i < ((window*2)+1); i++)
     {
         free(X[i]); X[i] = NULL;
         free(Y[i]); Y[i] = NULL;
     }

     free(X); X = NULL;
     free(Y); Y = NULL;            

     return;

}

void compute_Gaussian_intensity_weights(double* H[], double* I[], double* A, double sigma2, int a, int b, int size_r, int size_c, int i, int j, int height, int width)
{

     int r,c,s,t;

     for (r = 0; r < size_r; r++)
         for (c = 0; c < size_c; c++)
         {
             t = r + a;
             s = c + b;
             I[r][c] = A[(t*width+s)];
         }

         for (r = 0; r < size_r; r++)
            for (c = 0; c < size_c; c++)
                H[r][c] = I[r][c]- A[(i*width+j)];

     return;

}

double Bilateral_filter_response(double* F[], double* H[], double* I[], double* G[], int w, int size_r, int size_c, int a, int b, int i, int j)
{

       double sum1 = 0.0, sum2 = 0.0;
       int r,c,t,s;

         
       for (r = 0; r < size_r; r++)
           for (c = 0; c < size_c; c++)
           {
               t = r + a - i + w;
               s = c + b - j + w;
               F[r][c] = H[r][c]*G[t][s];
               sum2 = sum2 + F[r][c];
               sum1 = sum1 + (F[r][c]*I[r][c]);
           }

       return (sum1/sum2);

}

/*Functions*/

void ccRGBtoXYZ(unsigned char r, unsigned char g, unsigned char b, double* x, double* y, double* z)
{
  double var_R = ((double)(r))/255;        //R = From 0 to 1
  double var_G = ((double)(g))/255;        //G = From 0 to 1
  double var_B = ((double)(b))/255;        //B = From 0 to 1

  if ( var_R > 0.04045 ) var_R = pow(( ( var_R + 0.055 ) / 1.055 ),2.4);
  else                   var_R = var_R / 12.92;

  if ( var_G > 0.04045 ) var_G = pow(( ( var_G + 0.055 ) / 1.055 ),2.4);
  else                   var_G = var_G / 12.92;

  if ( var_B > 0.04045 ) var_B = pow(( ( var_B + 0.055 ) / 1.055 ),2.4);
  else                   var_B = var_B / 12.92;

  var_R = var_R * 100;
  var_G = var_G * 100;
  var_B = var_B * 100;

  *x = var_R * 0.4124 + var_G * 0.3576 + var_B * 0.1805;
  *y = var_R * 0.2126 + var_G * 0.7152 + var_B * 0.0722;
  *z = var_R * 0.0193 + var_G * 0.1192 + var_B * 0.9505;

  return;

}



void ccXYZtoCIE_Luv(double x, double y, double z, double* L, double* u, double* v)
{
  if ((x == 0.0) && (y == 0.0) && (z == 0.0))
  {
		*L = 0.0;
		*u = 0.0;
		*v = 0.0;
		return;
  }

  double ref_X = _ccTristimulusValues[D65_2][ccX];
  double ref_Y = _ccTristimulusValues[D65_2][ccY];

  double ref_Z = _ccTristimulusValues[D65_2][ccZ];
  double ref_U = ( 4 * ref_X ) / ( ref_X + ( 15 * ref_Y ) + ( 3 * ref_Z ) );
  double ref_V = ( 9 * ref_Y ) / ( ref_X + ( 15 * ref_Y ) + ( 3 * ref_Z ) );

  double var_U = ( 4 * x ) / ( x + ( 15 * y ) + ( 3 * z ) );
  double var_V = ( 9 * y ) / ( x + ( 15 * y ) + ( 3 * z ) );
  double var_Y = y / ref_Y;

  if ( var_Y > 0.008856 ) 
	  *L = (116 * pow(var_Y,( 1./3 ))) - 16;
  else
      *L = 903.3 * var_Y;

  *u = 13 * (*L) * ( var_U - ref_U );
  *v = 13 * (*L) * ( var_V - ref_V );
  
  return;
}

void ccXYZtoCIE_L(double y, double* L)
{
  double ref_Y = _ccTristimulusValues[D65_2][ccY];
  double var_Y = y / ref_Y;

  if ( var_Y > 0.008856 ) *L = (116 * pow(var_Y,( 1./3 ))) - 16;
  else                    *L = 903.3 * var_Y;

  return;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






#endif /* SEGMENT_IMAGE_H */
