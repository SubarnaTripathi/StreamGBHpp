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

/* Implements the node data structure. */

#ifndef DISJOINT_SET_H
#define DISJOINT_SET_H

#include <vector>

#include "image.h"
#include "misc.h"
#include "histogram.h"

#include "optical_flow/Image.h"

using namespace std;

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

// for color space transfer
const float ref_X = (float)95.047;
const float ref_Y = (float)100.000;
const float ref_Z = (float)108.883;

/* Each pixel has such a node structure. 
 * Once go up a level, changes cannot be made.
 * Foruther effort can be done by building a 
 * balanced searching tree.
 */
typedef struct {

	// contain its representative node in each level
	vector<int> p;

	// contain the output color in each level
	vector<rgb> color_output;

	// representative node
	bool root;

	// use to update MST mergeing condition
	float mst;
	int size;

	// record pixel level oversegmentation MST
	float mst_first;
	int size_first;

	// record the pixel color
	rgb color;
	float XYZ_X;
	float XYZ_Y;
	float XYZ_Z;
	float Lab_L;
	float Lab_a;
	float Lab_b;

	// for representative nodes
	Histogram<float> *His_L;
	Histogram<float> *His_a;
	Histogram<float> *His_b;

	// position in the 3D volume
	int x;
	int y;
	int z;

# ifdef MOTION_FEATURE
	// subarna: motion feature histograms
	float motion_u;
	float motion_v;
	float flow_orientation;
	float flow_mag;
	// for representative nodes
	Histogram<float> *His_flow;
# ifdef TWO_FEATURES
	Histogram<float> *His_flow_mag;
# endif
# endif


}node;

/* Node Forest for Subsequence/Clip (contain part info from previous one) */
class universe {
public:
	// all nodes in the forest
	vector<node> elts;

	// initialize 1st level forest
# ifndef MOTION_FEATURE
	universe(int frame_id_start, int frame_id_end, int width, int height,
			universe& v, int num_verticese_v, image<float> *r[],
			image<float> *g[], image<float> *b[], int hie_num);
# else
	universe(int frame_id_start, int frame_id_end, int width, int height,
			universe& v, int num_verticese_v, image<float> *r[],
			image<float> *g[], image<float> *b[], DImage *vx_motion[], DImage *vy_motion[], int hie_num);
# endif

	// color space transfer
	void RGBtoXYZ(int a);
	void XYZtoLab(int a);

# ifdef MOTION_FEATURE
	void UVtoOrientation(int a);
# endif

	// update an initial status of higher level
	void update(int level);

	// extract last frame of the subsequence
	universe(universe& mess, int num_vertices, int num_vertices_u);

	// destroy
	~universe();

	// find parent in one level
	int find_in_level(int a, int level);

	// semi-supervised merging
	int join(int a, int b, float mst, int level);

	// merging small segments
	void join_noise(int a, int b, int level);

	// set node parent in one level
	void set_p(int a, int p, int level) {
		elts[a].p[level] = p;
	}
	// get node parent in one level
	int get_p(int a, int level) {
		return elts[a].p[level];
	}

	// set node color in one level
	void set_color_output(int a, rgb c, int level) {
		elts[a].color_output[level] = c;
	}
	// get node color in one level
	rgb get_color_output(int a, int level) {
		return elts[a].color_output[level];
	}

	// get other information of one node
	node get_element(int a) { return elts[a]; }
	bool get_root(int a) { return elts[a].root; }
	float get_mst(int a) { return elts[a].mst; }
	int get_size(int a) { return elts[a].size; }
	Histogram<float> *get_His_L(int a) { return elts[a].His_L; }
	Histogram<float> *get_His_a(int a) { return elts[a].His_a; }
	Histogram<float> *get_His_b(int a) { return elts[a].His_b; }

# ifdef MOTION_FEATURE
	//subarna
	Histogram<float> *get_His_flow(int a) { return elts[a].His_flow; }
# ifdef TWO_FEATURES
	Histogram<float> *get_His_flow_mag(int a) { return elts[a].His_flow_mag; }
# endif
# endif

private:
	int num; // total nodes in the forest
	int board; // a board between previous clip's nodes and this clip's nodes
	int frame_id_start;
	int frame_id_end;
	int level_total;
};


# ifndef MOTION_FEATURE
universe::universe(int frame_id_start, int frame_id_end, int width, int height,
		universe& v, int num_vertices_v, image<float> *r[], image<float> *g[],
		image<float> *b[], int hie_num)
# else
universe::universe(int frame_id_start, int frame_id_end, int width, int height,
		universe& v, int num_vertices_v, image<float> *r[], image<float> *g[],
		image<float> *b[], DImage *vx_motion[], DImage *vy_motion[], int hie_num)
# endif 
{
	num = (frame_id_end - frame_id_start + 1) * width * height;

	board = num_vertices_v; // set board
	level_total = hie_num + 1;
	rgb blank;
	blank.r = blank.g = blank.b = 0;

	for (int i = 0; i < num_vertices_v; i++) {
		// fill out previous num_vertices_v nodes (previous frames)
		//elts.push_back(node());

		//subarna
		node a = node();
		elts.push_back(a);

		elts.back() = v.get_element(i);
		elts.back().mst = elts.back().mst_first;
		elts.back().size = elts.back().size_first;
# ifndef LUV
		// reset histogram
		elts.back().His_L = new Histogram<float>(20, 0, 100);
		elts.back().His_a = new Histogram<float>(20, -50, 50);
		elts.back().His_b = new Histogram<float>(20, -50, 50);
# else
		elts.back().His_L = new Histogram<float>(20, 0, 255);
		elts.back().His_a = new Histogram<float>(20, 0, 255);
		elts.back().His_b = new Histogram<float>(20, 0, 255);
# endif

# ifdef MOTION_FEATURE
# ifdef TWO_FEATURES
		//subarna: reset motion histogram
		elts.back().His_flow = new Histogram<float>(10, -30, 30);
		elts.back().His_flow_mag = new Histogram<float>(10, -30, 30); //subarna
# else
		//subarna: reset motion histogram
		elts.back().His_flow = new Histogram<float>(16, 0, 360);
# endif
# endif
	}

	for (int i = num_vertices_v; i < num; i++) {
		// initialize information of rest nodes (current frames)
		//elts.push_back(node());

		//subarna
		node a = node();
		elts.push_back(a);

		// p contains local id
		elts.back().p.resize(level_total);
		elts.back().p.front() = i;

		// forest root
		elts.back().root = true;

		// use in each iteration
		elts.back().mst = 0; // local
		elts.back().size = 1; // global

		// local position
		elts.back().z = i / (width * height);
		elts.back().y = (i % (width * height)) / width;
		elts.back().x = (i % (width * height)) % width;

		elts.back().color.r =
				imRef(r[elts.back().z], elts.back().x, elts.back().y);
		elts.back().color.g =
				imRef(g[elts.back().z], elts.back().x, elts.back().y);
		elts.back().color.b =
				imRef(b[elts.back().z], elts.back().x, elts.back().y);

# ifndef LUV
		RGBtoXYZ(i);
		XYZtoLab(i);
# else
		// already converted to RGB; no need to convert again
		elts[i].Lab_L = elts[i].color.r;
		elts[i].Lab_a = elts[i].color.g;
		elts[i].Lab_b = elts[i].color.b;
# endif

		elts.back().color_output.resize(level_total);
		for (int c_i = 0; c_i < level_total; c_i++)
			elts.back().color_output[c_i] = blank;

# ifndef LUV
		elts.back().His_L = new Histogram<float>(20, 0, 100);
		elts.back().His_a = new Histogram<float>(20, -50, 50);
		elts.back().His_b = new Histogram<float>(20, -50, 50);
# else
		elts.back().His_L = new Histogram<float>(20, 0, 255);
		elts.back().His_a = new Histogram<float>(20, 0, 255);
		elts.back().His_b = new Histogram<float>(20, 0, 255);
# endif

		//subarna
# ifdef MOTION_FEATURE
		elts.back().motion_u = 
			(float) vx_motion[elts.back().z]->pData[elts.back().x + elts.back().y*width];
		elts.back().motion_v = 
			(float) vy_motion[elts.back().z]->pData[elts.back().x + elts.back().y*width];

		UVtoOrientation(i);
		
# ifdef TWO_FEATURES
		elts.back().His_flow = new Histogram<float>(10, -30, 30); 
		elts.back().His_flow_mag = new Histogram<float>(10, -30, 30); //subarna
# else
		elts.back().His_flow = new Histogram<float>(16, 0, 360); 
# endif

# endif
	}
}

/* RGB to XYZ color space */
void universe::RGBtoXYZ(int a) {
	float var_R = elts[a].color.r / 255;
	float var_G = elts[a].color.g / 255;
	float var_B = elts[a].color.b / 255;

	if (var_R > 0.04045)
		var_R = pow(((var_R + 0.055) / 1.055), 2.4);
	else
		var_R = var_R / 12.92;

	if (var_G > 0.04045)
		var_G = pow(((var_G + 0.055) / 1.055), 2.4);
	else
		var_G = var_G / 12.92;

	if (var_B > 0.04045)
		var_B = pow(((var_B + 0.055) / 1.055), 2.4);
	else
		var_B = var_B / 12.92;

	var_R = var_R * 100;
	var_G = var_G * 100;
	var_B = var_B * 100;

	//Observer. = 2°, Illuminant = D65
	elts[a].XYZ_X = var_R * 0.412453f + var_G * 0.357580f + var_B * 0.180423f;
	elts[a].XYZ_Y = var_R * 0.212671f + var_G * 0.715160f + var_B * 0.072169f;
	elts[a].XYZ_Z = var_R * 0.019334f + var_G * 0.119193f + var_B * 0.950227f;
}

/* XYZ to Lab color space */
void universe::XYZtoLab(int a) {
	float var_X;
	float var_Y;
	float var_Z;

	// Observer= 2°, Illuminant= D65
	var_X = elts[a].XYZ_X / ref_X;
	var_Y = elts[a].XYZ_Y / ref_Y;
	var_Z = elts[a].XYZ_Z / ref_Z;
	if (var_X > 0.008856)
		var_X = pow((double) var_X, 0.333333);
	else
		var_X = (7.787 * var_X) + (16 / 116);
	if (var_Y > 0.008856)
		var_Y = pow((double) var_Y, 0.333333);
	else
		var_Y = (7.787 * var_Y) + (16 / 116);
	if (var_Z > 0.008856)
		var_Z = pow((double) var_Z, 0.333333);
	else
		var_Z = (7.787 * var_Z) + (16 / 116);

	elts[a].Lab_L = (116 * var_Y) - 16;
	elts[a].Lab_a = 500 * (var_X - var_Y);
	elts[a].Lab_b = 200 * (var_Y - var_Z);
}

// level start from 0
void universe::update(int level) {
	if (level == 0) {
		for (int i = 0; i < num; i++) {
			int p = find_in_level(i, level);
			elts[i].p[level + 1] = elts[i].p[level];
			elts[i].mst = 0;
			elts[p].His_L->addSample(elts[i].Lab_L);
			elts[p].His_a->addSample(elts[i].Lab_a);
			elts[p].His_b->addSample(elts[i].Lab_b);

			//subarna
# ifdef MOTION_FEATURE
# ifndef TWO_FEATURES
			elts[p].His_flow->addOriSample(elts[i].flow_orientation, elts[i].flow_mag);
# else
			elts[p].His_flow->addSample(elts[i].flow_orientation);
			elts[p].His_flow_mag->addSample(elts[i].flow_mag);
# endif
# endif
		}
	} else {
		for (int i = 0; i < num; i++) {
			int p = find_in_level(i, level);
			elts[i].p[level + 1] = elts[i].p[level];
			elts[i].mst = 0;
			if (i != p) {
				elts[p].His_L->mergeHistogram(*elts[i].His_L);
				elts[p].His_a->mergeHistogram(*elts[i].His_a);
				elts[p].His_b->mergeHistogram(*elts[i].His_b);

				//subarna
# ifdef MOTION_FEATURE
				elts[p].His_flow->mergeHistogram(*elts[i].His_flow);
				//elts[p].His_flow->mergeOriHistogram(*elts[i].His_flow);
# ifdef TWO_FEATURES
				elts[p].His_flow->mergeHistogram(*elts[i].His_flow_mag);
# endif
# endif			
			}
		}
	}
}

/* extract last frame of this subsequence */
universe::universe(universe& mess, int num_vertices, int num_vertices_u) {
	num = num_vertices_u;
	int offset = num_vertices - num_vertices_u;
	for (int i = 0; i < num_vertices_u; i++) {
		/*elts.push_back(mess.get_element(i + offset));
		elts.back().His_L = new Histogram<float>(*elts.back().His_L);
		elts.back().His_a = new Histogram<float>(*elts.back().His_a);
		elts.back().His_b = new Histogram<float>(*elts.back().His_b);*/

		//subarna
		node n = mess.get_element(i + offset);
		elts.push_back(n);

		Histogram<float> *t = elts.back().His_L;
		elts.back().His_L = new Histogram<float>(*t);
		t = elts.back().His_a;
		elts.back().His_a = new Histogram<float>(*t);
		t = elts.back().His_b;
		elts.back().His_b = new Histogram<float>(*t);

		//subarna
# ifdef MOTION_FEATURE
		elts.back().His_flow = new Histogram<float>(*elts.back().His_flow);
# ifdef TWO_FEATURES
		elts.back().His_flow_mag = new Histogram<float>(*elts.back().His_flow_mag);
# endif
# endif

		vector<int>::iterator it;
		for (it = elts.back().p.begin(); it != elts.back().p.end(); it++)
			*it -= offset;
	}
}

# ifdef MOTION_FEATURE
void universe::UVtoOrientation(int a)
{
	// subarna
	float motion_u = elts[a].motion_u;
	float motion_v = elts[a].motion_v;

	float param = motion_v/motion_u;
	float result = atan (param) * 180 / PI;

	/*if ((motion_v < 0) && (motion_u < 0) )
		result = 180 + result;
	else if (motion_v < 0)
		result = 360 + result;
	else if (motion_u < 0)
		result = 180 + result;*/

	if (motion_u < 0) 
		result = 180 + result;
	else if (motion_v < 0)
		result = 360 + result;
	
	elts[a].flow_orientation = result;
	elts[a].flow_mag = sqrt(motion_u*motion_u + motion_v*motion_v);
	// bin index = flow_orientation (result)
	// val in bin = flow_mag
}
# endif

universe::~universe() {
	for (int i = 0; i < num; i++) {
		delete elts[i].His_L;
		delete elts[i].His_a;
		delete elts[i].His_b;

		//subarna
# ifdef MOTION_FEATURE
		delete elts[i].His_flow;
# ifdef TWO_FEATURES
		delete elts[i].His_flow_mag;
# endif
# endif

	}
}

int universe::find_in_level(int a, int level) {
	int element = a;
	while (element != elts[element].p[level])
		element = elts[element].p[level];
	elts[a].p[level] = element;
	return element;
}

/* A semi-supervised merging.
 * Representative nodes are always assigned to the larger id node.
 * By doing this, last frame contains all alive node trees.
 */
int universe::join(int a, int b, float mst, int level) {
	rgb blank;
	blank.r = blank.g = blank.b = 0;

	if (elts[a].color_output[level] == blank
			&& elts[b].color_output[level] == blank) {
		// both node a and node b are new segment
		if (a >= b) {
			elts[b].p[level] = a;
			elts[a].size += elts[b].size;
			elts[b].root = false;
		} else {
			elts[a].p[level] = b;
			elts[b].size += elts[a].size;
			elts[a].root = false;
		}

		// syncing mst
		if (elts[a].mst >= elts[b].mst) {
			if (elts[a].mst < mst)
				elts[a].mst = mst;
			elts[b].mst = elts[a].mst;
		} else {
			if (elts[b].mst < mst)
				elts[b].mst = mst;
			elts[a].mst = elts[b].mst;
		}
		if (level == 0) {
			elts[a].mst_first = elts[a].mst;
			elts[b].mst_first = elts[a].mst;
			elts[a].size_first = elts[a].size;
			elts[b].size_first = elts[b].size;
		}

		return 0;

	} else if (elts[a].color_output[level] == blank) {
		// node a is new segment, node b is old segment
		if (a >= b) {
			elts[b].p[level] = a;
			elts[a].size += elts[b].size;
			elts[b].root = false;
		} else {
			elts[a].p[level] = b;
			elts[b].size += elts[a].size;
			elts[a].root = false;
		}

		for (int c_i = level; c_i < level_total; c_i++)
			elts[a].color_output[c_i] = elts[b].color_output[c_i];

		// syncing mst
		if (elts[a].mst >= elts[b].mst) {
			if (elts[a].mst < mst)
				elts[a].mst = mst;
			elts[b].mst = elts[a].mst;
		} else {
			if (elts[b].mst < mst)
				elts[b].mst = mst;
			elts[a].mst = elts[b].mst;
		}
		if (level == 0) {
			elts[a].mst_first = elts[a].mst;
			elts[b].mst_first = elts[a].mst;
			elts[a].size_first = elts[a].size;
			elts[b].size_first = elts[b].size;
		}

		return 0;

	} else if (elts[b].color_output[level] == blank) {
		// node a is old segment, node b is new segment
		if (a >= b) {
			elts[b].p[level] = a;
			elts[a].size += elts[b].size;
			elts[b].root = false;
		} else {
			elts[a].p[level] = b;
			elts[b].size += elts[a].size;
			elts[a].root = false;
		}

		for (int c_i = level; c_i < level_total; c_i++)
			elts[b].color_output[c_i] = elts[a].color_output[c_i];

		// syncing mst
		if (elts[a].mst >= elts[b].mst) {
			if (elts[a].mst < mst)
				elts[a].mst = mst;
			elts[b].mst = elts[a].mst;
		} else {
			if (elts[b].mst < mst)
				elts[b].mst = mst;
			elts[a].mst = elts[b].mst;
		}
		if (level == 0) {
			elts[a].mst_first = elts[a].mst;
			elts[b].mst_first = elts[a].mst;
			elts[a].size_first = elts[a].size;
			elts[b].size_first = elts[b].size;
		}
		return 0;
	} 
	else 
	{
# ifdef HIGHEST_LEVEL_MOTION_ONLY
		if (level == level_total - 1 )
		{
			if (a >= b) {
			elts[b].p[level] = a;
			elts[a].size += elts[b].size;
			elts[b].root = false;
			} else {
				elts[a].p[level] = b;
				elts[b].size += elts[a].size;
				elts[a].root = false;
			}

			// syncing mst
			if (elts[a].mst >= elts[b].mst) {
				if (elts[a].mst < mst)
					elts[a].mst = mst;
				elts[b].mst = elts[a].mst;
			} else {
				if (elts[b].mst < mst)
					elts[b].mst = mst;
				elts[a].mst = elts[b].mst;
			}		
			return 0;
		}
		else
# endif
		return 1;
	}

}

/* merging enforced minimum segments */
void universe::join_noise(int a, int b, int level) {
	if (elts[a].size >= elts[b].size) {
		if (a >= b) {
			elts[b].p[level] = a;
			elts[a].size += elts[b].size;
			elts[b].root = false;
		} else {
			elts[a].p[level] = b;
			elts[b].size += elts[a].size;
			elts[a].root = false;
		}
		for (int c_i = level; c_i < level_total; c_i++)
			elts[b].color_output[c_i] = elts[a].color_output[c_i];
	} else {
		if (a >= b) {
			elts[b].p[level] = a;
			elts[a].size += elts[b].size;
			elts[b].root = false;
		} else {
			elts[a].p[level] = b;
			elts[b].size += elts[a].size;
			elts[a].root = false;
		}
		for (int c_i = level; c_i < level_total; c_i++)
			elts[a].color_output[c_i] = elts[b].color_output[c_i];
	}

	// syncring mst
	if (elts[a].mst >= elts[b].mst) {
		elts[b].mst = elts[a].mst;
	} else {
		elts[a].mst = elts[b].mst;
	}
	if (level == 0) {
		elts[a].mst_first = elts[a].mst;
		elts[b].mst_first = elts[b].mst;
		elts[a].size_first = elts[a].size;
		elts[b].size_first = elts[b].size;
	}
}

#endif /* DISJOINT_SET_H */
