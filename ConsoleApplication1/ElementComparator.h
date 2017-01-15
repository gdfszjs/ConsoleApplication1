#pragma once
#include <string>
#include <vector>
#include <LibIV/libiv.h>
#include <GL/glut.h>
#include <fstream>	 
#include "MRG.h"
#include "Mesh.h"

class ElementComparator
{
public:
	ElementComparator();
	ElementComparator(Element *a, Element *b);
	~ElementComparator();
	Mnode * createFictitiousNode(int ngraph, vector<int> Tset,vector<int> ad_index,double p1,double p2);
	Mnode* adj(Mnode *m);
	double loss(Mnode *m, Mnode *n);
	double mat(Mnode *m, Mnode *n);
	double sim(Mnode *m,Mnode *n);
	double getCompareResult();
	void Initialization();
	void Matching();
	vector<Mnode*> completeNodeVector(vector<int> old_node_index_vector);
	void Unpacking(vector<int> pair_index, vector<Mnode*> one_pair);

	bool HasSameRange(Mnode *m, Mnode *n);
	bool HasParentPair(Mnode *m, Mnode *n);
	bool HasSameMLIST(Mnode *m, Mnode *n);
	bool BelongToDifferentPart(Mnode *m, Mnode *n);
	bool InOnePair(vector<Mnode*> one_pair,int node_graph,int node_index);
	void BoardCastLabel(int index);
	void extendMLISTinOneDirection(Mnode *e,int direcion, int index);

	bool isSameNode(Mnode *m, Mnode *n);
	void getNodeInformation(Mnode *m);
	void calculate_result();
private:
	ofstream out6;
	MRG* MRG_vector[2];
	MRG* R;
	MRG* S; 
	vector<Mnode*> NLIST;
	vector<vector<Mnode*>> MPAIR;
	double result;
};