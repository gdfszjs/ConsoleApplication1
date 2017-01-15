#pragma once
#include <vector>
#include "EdgeSet.h"
#include <LibIV\libiv.h>
using namespace std;
class Node
{
public:
	Node();
	Node(int index);
	Node(int index, double x_cor, double y_cor);
	~Node();
	int getNodeIndex() { return this->node_index; };
	double getValue() { return this->value; };
	double getFunction_Value() { return this->function_value; };
	double getX() { return this->x_cor; };
	double getY() { return this->y_cor; };
	void modifyValue(double new_value) { this->value = new_value; };
	void modifyFunctionValue(double new_value) { this->function_value = new_value; };
	int getAdjanceNodeNumber() { return this->adjance_node_number; };
	Node* getAdjanceNode(int index) { return this->adjance_node[index]; };
	void extendAdjanceNode(Node* e);
	void printAdjanceNode();

private:
	int node_index;
	double x_cor;
	double y_cor;
	double value;
	double function_value;
	Node** adjance_node;
	int adjance_node_number;
};

class Graph
{
public:
	Graph();
	Graph(int NODESIZE,double r);
	Graph(vector<v2d> NODE, double r);
	~Graph();
	void initAdjancePair(vector<v2i> edge_vector, vector<little_polygon*> polygon_vector);
	void checkAdjancePair();
	void insertintoVLISTwithAscendant(Node * e);
	void insertBIandAREA(int index, double area);
	double DijkstraDistanceForBI(int base_vertex);
	double det(int* p0, int* p1, int* p2){return (p1[0] - p0[0])*(p2[1] - p0[1]) - (p1[1] - p0[1])*(p2[0] - p0[0]);}
	void getResult(double Pack[240][2], double x1, double y1, double x2, double y2);
	double CalculateAreawithMelkman(vector<v2d> node_vector);
	void DijkstraDistance(int base_vertex);
	void calculateGD();
	void function_valueNormalized();
	void initCandidateBI();
	void removeCandidateBI(int vertex);
	int getBIfromCandidate();
	void findBI();
	void checkFunctionValue();
	void checkVLIST();
	void checkGD_Node();
	double returndistance(Node * a, Node *b) { return abs(sqrt(pow((a->getX() - b->getX()), 2) + pow((a->getY() - b->getY()), 2))); };
	int removeMininVLIST();
	Node* returnG_Graph(int index) { return this->GD_Graph[index]; };
private:
	double thread_hole;
	
	int g_result[240][2];//for the Malkman
	int candidate_number;
	int* candidate_BI;

	int* BI;
	double* BIarea;
	double **¦ÌMatrix;
	int BI_number;

	Node** GD_Graph;
	int node_number;

	Node** VLIST;
	int vlist_size;
};