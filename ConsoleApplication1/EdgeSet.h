#pragma once
#include <vector>
#include <LibIV\libiv.h>
using namespace std;

struct point
{
	double x;
	double y;
};

struct edge
{
	int end_point1_index;
	int end_point2_index;
};

class little_polygon
{
public:
	little_polygon();
	little_polygon(int p1, int p2, int p3);
	~little_polygon();
	void initedge();
	void initshortcut(vector<v2d> short_cut_vector);
	int returnVertexsbyIndex(int index) { return this->vertexs[index]; };
	int getVertexsNumber() { return this->numberofvertexs; }
	int getShortCutNumber() { return this->numberofshort_cut; };
	int getEdgeNumber() { return this->numberofpolygon_edge; };
	edge getEdge(int index) { return this->polygon_edge[index]; };
	edge getShortCut(int index) { return this->short_cut[index]; };
	void addShortCut(vector<v2d> vertexs_vector);
	int intersect2(v2d u1, v2d u2, v2d v1, v2d v2);
	double multiply(v2d p1, v2d p2, v2d  p0) { return((p1[0] - p0[0])*(p2[1] - p0[1]) - (p2[0] - p0[0])*(p1[1] - p0[1])); }
	double determinant(double v1, double v2, double v3, double v4) { return (v1*v3 - v2*v4); }
	bool intersect(v2d aa, v2d bb, v2d cc, v2d dd);
	bool isShortCut(int c1, int c2, int c1_index, int c2_index, vector<v2d> vertexs_vector);
	bool isBoundary(int c1, int c2);
	void addEdgeByNeighbor(int index, vector<v3i>mesh_vector, vector<v3i> neighbor_vector);
	bool isPointonTriangle(int vertex_index, int mesh_index, vector<v3i>mesh_vector);
	int getAnotherPointOnTriangle(int vertex1_index, int vertex2_index, int mesh_index, vector<v3i> mesh_vector);
	void modifypolygon(vector<int> new_vertexs_vector);
private:
	int* vertexs;
	edge* exitedge;
	int numberofvertexs;
	edge* polygon_edge;
	int numberofpolygon_edge;
	edge* short_cut;
	int numberofshort_cut;
};