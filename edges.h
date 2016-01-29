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

/* Implements the edge data structure. */

#ifndef EDGES_H
#define EDGES_H

#include <vector>

#include "image.h"
#include "disjoint-set.h"
#include "histogram.h"

#include "optical_flow/Image.h"


using namespace std;

/* define a edge */
typedef struct {
	float w;
	int a, b;
} edge;

int validate_range(int val, int min, int max)
{
	int new_val = val;

	if (new_val < min) new_val = min;
	if (new_val > max) new_val = max;

	return new_val;
}

/* fill pixel level edges */
void generate_edge(edge *e, image<float> *r_v, image<float> *g_v,
		image<float> *b_v, image<float> *r_u, image<float> *g_u,
		image<float> *b_u, int x_v, int y_v, int z_v, int x_u, int y_u,
		int z_u) 
{
	int width = r_v->width();
	int height = r_v->height();

	e->a = y_v * width + x_v + z_v * (width * height);
	e->b = y_u * width + x_u + z_u * (width * height);
	e->w = sqrt(
          square(imRef(r_v, x_v, y_v) - imRef(r_u, x_u, y_u))
					+ square(imRef(g_v, x_v, y_v) - imRef(g_u, x_u, y_u))
					+ square(imRef(b_v, x_v, y_v) - imRef(b_u, x_u, y_u)));
}

/* fill region graph edges */
# ifndef HIGHEST_LEVEL_MOTION_ONLY
		void fill_edge_weight(vector<edge> edges_region, universe *mess, int level) {
# else
	void fill_edge_weight(vector<edge> edges_region, universe *mess, int level, int max_level) {
# endif
	float d_c = 0;
	float d_f = 0;
# ifdef TWO_FEATURES
	float d_f2 = 0;
# endif

	for (int i = 0; i < ((int)edges_region.size()); i++) {
		int a = edges_region[i].a;
		int b = edges_region[i].b;
		int a_p = mess->find_in_level(a, level);
		int b_p = mess->find_in_level(b, level);
		if (a_p != b_p) {
			float d_c = 
				mess->get_His_L(a_p)->chiSquared(*mess->get_His_L(b_p))
					+ mess->get_His_a(a_p)->chiSquared(*mess->get_His_a(b_p))
					+ mess->get_His_b(a_p)->chiSquared(*mess->get_His_b(b_p));
			edges_region[i].w = d_c;

# ifdef	HIGHEST_LEVEL_MOTION_ONLY
			if (level == max_level)
			{
				d_c = 0;
				edges_region[i].w = 0;
			}
# endif				
				
# ifdef MOTION_FEATURE	
		    d_f = mess->get_His_flow(a_p)->chiSquared(*mess->get_His_flow(b_p));
# ifdef TWO_FEATURES
			d_f2 = mess->get_His_flow_mag(a_p)->chiSquared(*mess->get_His_flow_mag(b_p));
# endif

# ifdef TWO_FEATURES
			d_c = d_c/3;	
			edges_region[i].w = (1 - (1-d_c)*(1-d_f2) /**(1-d_f)*/);
			edges_region[i].w *= edges_region[i].w;
# else
			d_c = d_c/3;	
			edges_region[i].w = (1 - (1-d_c)*(1-d_f));
			edges_region[i].w *= edges_region[i].w;
# endif
			//edges_region[i].w = d_c + d_f;
# endif

# ifdef DIST_EXP
		edges_region[i].w = (1-exp(-d_c))*(1-exp(-d_f));
# ifdef TWO_FEATURES
		edges_region[i].w *= (1-exp(-d_f2));
# endif
		edges_region[i].w = 1 - edges_region[i].w;
# endif

		} 
		else 
		{
			edges_region[i].w = 0;
		}
	}
}

/* initialize pixel level edges */
# ifndef FLOW_EDGE
void initialize_edges(edge *edges, bool first, int num_frame, int width,
		int height, image<float> *smooth_r[], image<float> *smooth_g[],
		image<float> *smooth_b[]) 
# else
void initialize_edges(edge *edges, bool first, int num_frame, int width,
		int height, image<float> *smooth_r[], image<float> *smooth_g[],
		image<float> *smooth_b[], DImage *vx_motion[], DImage *vy_motion[])
# endif 

{
	int num_frame_start;
	if (first == true) {
		num_frame_start = 0;
	} else
		num_frame_start = 1;

	int num_edges = 0;

	int xx, yy;

# ifdef FLOW_EDGE
	DImage *vx_bck_motion = new DImage (width, height);
	DImage *vy_bck_motion = new DImage (width, height);
# endif


	//subarna: connectivity on motion edge to be added
	for (int z = num_frame_start; z < num_frame; z++) {
# ifdef FLOW_EDGE
		if ( z > 0 )
		{
			for (int h = 0; h < height; h++)
			{
				for (int w = 0; w < width; w++)
				{
					int x = (int)(vx_motion[z-1]->pData[h*width+w]);
					int y = (int)(vy_motion[z-1]->pData[h*width+w]);

					int new_x = w+x;
					int new_y = h+y;

					new_x = validate_range(new_x, 0, width-1);
					new_y = validate_range(new_y, 0, height-1);

					vx_bck_motion->pData[new_y*width+new_x]= -x;
					vy_bck_motion->pData[new_y*width+new_x]= -y;
				}
			}
		}
# endif

		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				// in the same frame
				if (x < width - 1) {
					generate_edge(&edges[num_edges], smooth_r[z], smooth_g[z],
							smooth_b[z], smooth_r[z], smooth_g[z], smooth_b[z],
							x + 1, y, z, x, y, z);
					num_edges++;
				}
				if (y < height - 1) {
					generate_edge(&edges[num_edges], smooth_r[z], smooth_g[z],
							smooth_b[z], smooth_r[z], smooth_g[z], smooth_b[z],
							x, y + 1, z, x, y, z);
					num_edges++;
				}
				if ((x < width - 1) && (y < height - 1)) {
					generate_edge(&edges[num_edges], smooth_r[z], smooth_g[z],
							smooth_b[z], smooth_r[z], smooth_g[z], smooth_b[z],
							x + 1, y + 1, z, x, y, z);

					num_edges++;
				}
				if ((x < width - 1) && (y > 0)) {
					generate_edge(&edges[num_edges], smooth_r[z], smooth_g[z],
							smooth_b[z], smooth_r[z], smooth_g[z], smooth_b[z],
							x + 1, y - 1, z, x, y, z);
					num_edges++;
				}

				// to the previous frame
				if (z > 0) {
					xx = x;
					yy = y;
# ifdef FLOW_EDGE
					//subarna
					// use (z-1) as per the optical_flow image indices
					xx += vx_bck_motion->pData[y*width+x]; //use the opposite dircetion
					yy += vx_bck_motion->pData[y*width+x]; //use the opposite dircetion

					xx = validate_range(xx, 0, width-1);
					yy = validate_range(yy, 0, height-1);
# endif

					if ( xx >= 0 && xx < width && yy >= 0 && yy < height)
						generate_edge(&edges[num_edges], smooth_r[z - 1],
								smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
								smooth_g[z], smooth_b[z], xx, yy, z - 1, x, y, z); //xx may not be x && yy may not be y

					else
						generate_edge(&edges[num_edges], smooth_r[z - 1],
								smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
								smooth_g[z], smooth_b[z], x, y, z - 1, x, y, z);
					num_edges++;


					if (x > 0 && x < width - 1 && y > 0 && y < height - 1) 
					{
						// additional 8 edges
						// x - 1, y - 1
# ifdef FLOW_EDGE
						if (xx > 0 && xx < width - 1 && yy > 0 && yy < height - 1)
						{
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx - 1, yy - 1, z - 1,
									x, y, z);
							num_edges++;
							// x, y - 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx, yy - 1, z - 1, x, y,
									z);
							num_edges++;
							// x + 1, y - 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx + 1, yy - 1, z - 1,
									x, y, z);
							num_edges++;
							// x - 1, y
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx - 1, yy, z - 1, x, y,
									z);
							num_edges++;
							// x + 1, y
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx + 1, yy, z - 1, x, y,
									z);
							num_edges++;
							// x - 1, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx - 1, yy + 1, z - 1,
									x, y, z);
							num_edges++;
							// x, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx, yy + 1, z - 1, x, y,
									z);
							num_edges++;
							// x + 1, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx + 1, yy + 1, z - 1,
									x, y, z);
							num_edges++;
						}
						else
# endif
						{
							// additional 8 edges
							// x - 1, y - 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									x - 1, y - 1, z - 1,
									x, y, z);
							num_edges++;
							// x, y - 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									x, y - 1, z - 1, x, y,
									z);
							num_edges++;
							// x + 1, y - 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									x + 1, y - 1, z - 1,
									x, y, z);
							num_edges++;
							// x - 1, y
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									x - 1, y, z - 1, x, y,
									z);
							num_edges++;
							// x + 1, y
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									x + 1, y, z - 1, x, y,
									z);
							num_edges++;
							// x - 1, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									x - 1, y + 1, z - 1,
									x, y, z);
							num_edges++;
							// x, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									x, y + 1, z - 1, x, y,
									z);
							num_edges++;
							// x + 1, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									x + 1, y + 1, z - 1,
									x, y, z);
							num_edges++;
						}
					} 
					else if (x == 0 && y > 0 && y < height - 1) 
					{
#ifdef FLOW_EDGE
						// additional 5 edges
						// x, y - 1
						if (xx == 0 && yy > 0 && yy < height - 1) 
						{
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx, yy - 1, z - 1, x, y,
									z);
							num_edges++;
							// x + 1, y - 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx + 1, yy - 1, z - 1,
									x, y, z);
							num_edges++;
							// x + 1, y
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx + 1, yy, z - 1, x, y,
									z);
							num_edges++;
							// x, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx, yy + 1, z - 1, x, y,
									z);
							num_edges++;
							// x + 1, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx + 1, yy + 1, z - 1,
									x, y, z);
							num_edges++;
						}
						else
# endif
						{
							// additional 5 edges
							// x, y - 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x, y - 1, z - 1, x, y,
									z);
							num_edges++;
							// x + 1, y - 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x + 1, y - 1, z - 1,
									x, y, z);
							num_edges++;
							// x + 1, y
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x + 1, y, z - 1, x, y,
									z);
							num_edges++;
							// x, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x, y + 1, z - 1, x, y,
									z);
							num_edges++;
							// x + 1, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x + 1, y + 1, z - 1,
									x, y, z);
							num_edges++;
						}
					} 
					else if (x == width - 1 && y > 0 && y < height - 1) 
					{
						// additional 5 edges
						// x - 1, y - 1
# ifdef FLOW_EDGE
						if (xx == width - 1 && yy > 0 && yy < height - 1) 
						{
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx - 1, yy - 1, z - 1,
									x, y, z);
							num_edges++;
							// x, y - 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx, yy - 1, z - 1, x, y,
									z);
							num_edges++;
							// x - 1, y
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx - 1, yy, z - 1, x, y,
									z);
							num_edges++;
							// x - 1, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z],
									xx - 1, yy + 1, z - 1,
									x, y, z);
							num_edges++;
							// x, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx, yy + 1, z - 1, x, y,
									z);
							num_edges++;
						}
						else
# endif
						{
							// additional 5 edges
							// x - 1, y - 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
								smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
								smooth_g[z], smooth_b[z], x - 1, y - 1, z - 1,
									x, y, z);
							num_edges++;
							// x, y - 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x, y - 1, z - 1, x, y,
									z);
							num_edges++;
							// x - 1, y
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x - 1, y, z - 1, x, y,
									z);
							num_edges++;
							// x - 1, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x - 1, y + 1, z - 1,
									x, y, z);
							num_edges++;
							// x, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x, y + 1, z - 1, x, y,
									z);
							num_edges++;
						}
					} 
					else if (y == 0 && x > 0 && x < width - 1) 
					{
						// additional 5 edges
						// x - 1, y
# ifdef FLOW_EDGE
						if (yy == 0 && xx > 0 && xx < width - 1) 
						{
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx - 1, yy, z - 1, x, y,
									z);
							num_edges++;
							// x + 1, y
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx + 1, yy, z - 1, x, y,
									z);
							num_edges++;
							// x - 1, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z],
									xx - 1, yy + 1, z - 1,
									x, y, z);
							num_edges++;
							// x, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx, yy + 1, z - 1, x, y,
									z);
							num_edges++;
							// x + 1, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx + 1, yy + 1, z - 1,
									x, y, z);
							num_edges++;
						}
						else
# endif
						{
							// additional 5 edges
							// x - 1, y
							generate_edge(&edges[num_edges], smooth_r[z - 1],
								smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
								smooth_g[z], smooth_b[z], x - 1, y, z - 1, x, y,
								z);
							num_edges++;
							// x + 1, y
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x + 1, y, z - 1, x, y,
									z);
							num_edges++;
							// x - 1, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x - 1, y + 1, z - 1,
									x, y, z);
							num_edges++;
							// x, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x, y + 1, z - 1, x, y,
									z);
							num_edges++;
							// x + 1, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x + 1, y + 1, z - 1,
									x, y, z);
							num_edges++;
						}
					} 
					else if (y == height - 1 && x > 0 && x < width - 1) 
					{
						// additional 5 edges
						// x - 1, y - 1
# ifdef FLOW_EDGE
						if (yy == height - 1 && xx > 0 && xx < width - 1) 
						{
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx - 1, yy - 1, z - 1,
									x, y, z);
							num_edges++;
							// x, y - 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx, yy - 1, z - 1, x, y,
									z);
							num_edges++;
							// x + 1, y - 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx + 1, yy - 1, z - 1,
									x, y, z);
							num_edges++;
							// x - 1, y
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z],
									xx - 1, yy, z - 1, x, y,
									z);
							num_edges++;
							// x + 1, y
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx + 1, yy, z - 1, x, y,
									z);
							num_edges++;
						}
						else
# endif
						{
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x - 1, y - 1, z - 1,
									x, y, z);
							num_edges++;
							// x, y - 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x, y - 1, z - 1, x, y,
									z);
							num_edges++;
							// x + 1, y - 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x + 1, y - 1, z - 1,
									x, y, z);
							num_edges++;
							// x - 1, y
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x - 1, y, z - 1, x, y,
									z);
							num_edges++;
							// x + 1, y
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x + 1, y, z - 1, x, y,
									z);
							num_edges++;
						}
					} 
					else if (x == 0 && y == 0) 
					{
						// additional 3 edges
						// x + 1, y
# ifdef FLOW_EDGE
						if (xx == 0 && yy == 0) 
						{
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx + 1, yy, z - 1, x, y,
									z);
							num_edges++;
							// x, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z],
									xx, yy + 1, z - 1, x, y,
									z);
							num_edges++;
							// x + 1, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx + 1, yy + 1, z - 1,
									x, y, z);
							num_edges++;
						}
						else
# endif
						{
							generate_edge(&edges[num_edges], smooth_r[z - 1],
								smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
								smooth_g[z], smooth_b[z], x + 1, y, z - 1, x, y,
								z);
							num_edges++;
							// x, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x, y + 1, z - 1, x, y,
									z);
							num_edges++;
							// x + 1, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x + 1, y + 1, z - 1,
									x, y, z);
							num_edges++;
						}
					} 
					else if (x == 0 && y == height - 1) 
					{
						// additional 3 edges
						// x, y - 1
# ifdef FLOW_EDGE
						if(xx == 0 && yy == height - 1) 
						{
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx, yy - 1, z - 1, x, y,
									z);
							num_edges++;
							// x + 1, y - 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx + 1, yy - 1, z - 1,
									x, y, z);
							num_edges++;
							// x + 1, y
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx + 1, yy, z - 1, x, y,
									z);
							num_edges++;
						}
						else
# endif
						{
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x, y - 1, z - 1, x, y,
									z);
							num_edges++;
							// x + 1, y - 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x + 1, y - 1, z - 1,
									x, y, z);
							num_edges++;
							// x + 1, y
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x + 1, y, z - 1, x, y,
									z);
							num_edges++;

						}
					} 
					else if (x == width - 1 && y == 0) 
					{
						// additional 3 edges
						// x - 1, y
# ifdef FLOW_EDGE
						if (xx == width - 1 && yy == 0) 
						{
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx - 1, yy, z - 1, x, y,
									z);
							num_edges++;
							// x - 1, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx - 1, yy + 1, z - 1,
									x, y, z);
							num_edges++;
							// x, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx, yy + 1, z - 1, x, y,
									z);
							num_edges++;
						}
						else
# endif
						{
							generate_edge(&edges[num_edges], smooth_r[z - 1],
								smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
								smooth_g[z], smooth_b[z], x - 1, y, z - 1, x, y,
								z);
							num_edges++;
							// x - 1, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x - 1, y + 1, z - 1,
									x, y, z);
							num_edges++;
							// x, y + 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x, y + 1, z - 1, x, y,
									z);
							num_edges++;
						}
					} 
					else if (x == width - 1 && y == height - 1) 
					{
						// additional 3 edges
						// x - 1, y - 1
# ifdef FLOW_EDGE
						if (xx == width - 1 && yy == height - 1) 
						{
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx - 1, yy - 1, z - 1,
									x, y, z);
							num_edges++;
							// x, y - 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx, yy - 1, z - 1, x, y,
									z);
							num_edges++;
							// x - 1, y
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], 
									xx - 1, yy, z - 1, x, y,
									z);
							num_edges++;
						}
						else
# endif
						{
							generate_edge(&edges[num_edges], smooth_r[z - 1],
								smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
								smooth_g[z], smooth_b[z], x - 1, y - 1, z - 1,
								x, y, z);
							num_edges++;
							// x, y - 1
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x, y - 1, z - 1, x, y,
									z);
							num_edges++;
							// x - 1, y
							generate_edge(&edges[num_edges], smooth_r[z - 1],
									smooth_g[z - 1], smooth_b[z - 1], smooth_r[z],
									smooth_g[z], smooth_b[z], x - 1, y, z - 1, x, y,
									z);
							num_edges++;
						}
					}
				}
			}
		}
	}

# ifdef FLOW_EDGE
	delete(vx_bck_motion);
	delete(vy_bck_motion);
# endif

//	printf("num_edges = %d\n", num_edges);
}

#endif /* EDGES_H */
