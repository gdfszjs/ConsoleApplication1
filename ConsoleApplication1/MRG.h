#pragma once
#include <string>
#include <vector>
#include <LibIV/libiv.h>
#include <GL/glut.h>
#include <fstream>
#include <math.h>
#include "GraphandNode.h"
struct Mnode
{
	int index;
	int ngraph;
	int nlevel;
	int nrange;
	int nnumber;

	int* Tset;
	int Tset_number;

	int* adjance_index;
	int adjance_number;

	int parent_index;

	int* child_index;
	int child_number;

	double p1;
	double p2;

	int** line_table;
	int edge_number;

	vector<int> MLIST;
	int visit_flag = -1;

	int fictitious_flag = 0;
};

class MRG
{
public:
	MRG();
	~MRG();
	void initMRG(int K_range);
	void createMnode(int graph, int level, int range, vector<int> Tset, int ** MeshonEdge_ForEachNode, int edge_number);
	void insertMnodeIndex(int index, int graph,int level, int range);

	void printLine_Queue();
	void printGraph_Number();


	void createMRG(int mesh_number);
	void getAdjanceMnodeinlevel(int level);
	int isNodeAllVisited(vector<int> node_vector, int* visit_flag);
	void getAdjanceMnodeinFinestlevel(int mesh_number);

	void calculateP12(vector<double> whole_area,vector<double> whole_len);
	void resetNGraph(int e);
	int returnLevel_Number() { return this->level_number; };
	int returnRangeByLevel(int level) { return this->range_in_each_level[level]; };
	Mnode* returnNodeByIndex(int index);
	Mnode* returnTopNode();
private:

	ofstream out4;
	//the range number in the finest level
	int K_range;

	//the node number(in mess)
	int node_number;

	//the property of each node in each level
	int* node_prop;

	//all node in the graph
	Mnode* Mnode_vector;

	//the complex graph(with the node index)
	int*** MGraph;
	//the node number(in the graph mode above)
	int** graph_number;
	//the tree number in the graph
	int level_number;
	//the range in each level
	int* range_in_each_level;
};