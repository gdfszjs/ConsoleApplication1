#pragma once
#include <string>
#include <vector>
#include <LibIV/libiv.h>
#include <GL/glut.h>
#include <fstream>
#include "GraphandNode.h"
#include "MRG.h"
using namespace std;
struct vertex
{
	double x;
	double y;
	double value;
	int range;
};

struct boundry
{
	double value;
	int end_point1_index;
	int end_point2_index;
};

struct Mesh
{
	int range;
	int* vertexs;
	int* boundrys_vector;
	int* adjance_meshs;
	int* range_node;
	double area;
};

class Element
{
public:
	Element();
	Element(vector<v2d> points_vector, vector<v3i> meshs_vector, vector<v3i> neighbors_vector, vector<v2i> edge_vector);
	~Element();
	double Heron(int a, int b, int c);
	void computeArea();
	void elementpartition();
	void elementinit(int K, Graph* e);
	void getVertexsValue(Graph* e);
	void calculateBoundryValue();
	void getHeight();
	void remeshelement();
	int isMeshAllVisited(vector<int> mesh_vector,int* visit_flag);
	void buildFinestMRG();
	void rebuildAdjanceMesh();
	void initK_rangeImage();
	void setK_range(int K) { this->K_range = K; };

	void drawElement();

	void checkVertexRange();
	void checkMeshRange();
	void checkHeight();
	void checkEdge();
	void checkMesh();

	void reorganizeEdge(int e1, int e2, vector<int> point_on_line);
	void reorganizeMesh(vector<int> temp_mesh_vector);
	void insertNewVertex(double x,double y,double value,int range);
	void insertNewMesh(int p1,int p2,int p3);
	void insertNewEdge(int e1, int e2);
	void checkVertex(int index);
	int isLineinEdgeVector(int e1,int e2);
	int isLineinFixedEdgevector(vector<vector<int>>f_e_v,int e1,int e2);

	double returnMaxValue(double a, double b, double c) { return (a > b) ? (a > c ? a : c) : (b > c ? b : c); };
	double returnMinValue(double a, double b, double c) { return (a < b) ? (a < c ? a : c) : (b < c ? b : c); };
	double returnMinFHeight() { return this->min_f_value; };
	double returnMaxFHeight() { return this->max_f_value; };
	double returnMaxHeight() { return this->max_height; };
	double returnMinHeight() { return this->min_height; };
	double returnHeight() { return this->height; };	
	int returnPoint_Number() { return this->point_number; };
	int returnMesh_Number() { return this->mesh_number; };
	int returnBoundry_Number() { return this->boundary_number; };
	vertex returnPoints(int index) { return this->points_vector[index]; };
	Mesh returnMeshs(int index) { return this->meshs_vector[index]; };
	boundry returnBoundrys(int index) { return this->boundarys_vector[index]; };
	MRG* returnMRG() { return this->mrg; };
private:

	ofstream out;

	double whole_area;
	double* height_scanning_line;
	int K_range;

	vector<int> fixing_mesh_index_vector;

	double max_f_value;
	double min_f_value;
	int min_value_index;
	int range_radius;
	double max_height;
	double min_height;
	double height;

	int point_number;
	vertex* points_vector;

	int boundary_number;
	boundry*  boundarys_vector;

	int mesh_number;
	Mesh* meshs_vector;
	int** MeshonEdge;

	MRG* mrg;
};