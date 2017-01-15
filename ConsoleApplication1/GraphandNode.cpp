#include "GraphandNode.h"
#include <LibIV\libiv.h>
#include <iostream> 
#include <time.h> 
#define MAX_VALUE 10000

using namespace std;
Node::Node()
{
	this->node_index = -1;
	this->x_cor = -1;
	this->y_cor = -1;
	this->value = MAX_VALUE;
	this->function_value = -1;
	this->adjance_node = NULL;
	this->adjance_node_number = 0;
}

Node::Node(int index)
{
	this->node_index = index;
	this->x_cor = -1;
	this->y_cor = -1;
	this->value = MAX_VALUE;
	this->function_value = -1;
	this->adjance_node = NULL;
	this->adjance_node_number = 0;
}

Node::Node(int index, double x_cor, double y_cor)
{
	this->node_index = index;
	this->x_cor = x_cor;
	this->y_cor = y_cor;
	this->value = MAX_VALUE;
	this->function_value = -1;
	this->adjance_node = NULL;
	this->adjance_node_number = 0;
}

Node::~Node()
{

}

void Node::extendAdjanceNode(Node* e)
{
	if (this->adjance_node == NULL)
	{
		this->adjance_node = (Node**)malloc(sizeof(Node*) * 1);
		this->adjance_node[0] = e;
		this->adjance_node_number = 1;
	}
	else
	{
		this->adjance_node = (Node**)realloc(this->adjance_node, sizeof(Node) * (this->adjance_node_number + 1));
		this->adjance_node[this->adjance_node_number] = e;
		this->adjance_node_number++;
	}
}

void Node::printAdjanceNode()
{
	for (int i = 0; i < this->adjance_node_number; i++)
	{
		cout << this->adjance_node[i]->node_index << " ";
	}
}

Graph::Graph()
{
	this->GD_Graph = NULL;
	this->node_number = 0;

	this->vlist_size = 0;
	this->VLIST = NULL;
	this->thread_hole = 0;

	this->candidate_number = 0;
	this->candidate_BI = NULL;

	this->BI = NULL;
	this->BIarea = NULL;
	this->μMatrix = NULL;
	this->BI_number = 0;
}

Graph::Graph(int NODESIZE, double r)
{
	this->GD_Graph = new Node*[NODESIZE];
	for (int i = 0; i < NODESIZE; i++)
	{
		this->GD_Graph[i] = new Node(i);
	}
	this->node_number = NODESIZE;

	this->vlist_size = 0;
	this->VLIST = NULL;

	this->thread_hole = r;

	this->candidate_number = 0;
	this->candidate_BI = NULL;

	this->BI = NULL;
	this->BIarea = NULL;
	this->μMatrix = NULL;
	this->BI_number = 0;
}

Graph::Graph(vector<v2d> NODE, double r)
{
	this->GD_Graph = new Node*[NODE.size()];
	for (int i = 0; i < NODE.size(); i++)
	{
		this->GD_Graph[i] = new Node(i, NODE.at(i)[0], NODE.at(i)[1]);
	}
	this->node_number = NODE.size();

	this->vlist_size = 0;
	this->VLIST = NULL;

	this->thread_hole = r;

	this->candidate_number = 0;
	this->candidate_BI = NULL;

	this->BI = NULL;
	this->BIarea = NULL;
	this->μMatrix = NULL;
	this->BI_number = 0;
}

Graph::~Graph()
{
}

void Graph::initAdjancePair(vector<v2i> edge_vector, vector<little_polygon*> polygon_vector)
{
	for (int i = 0; i < edge_vector.size(); i++)
	{
		int e1 = edge_vector.at(i)[0];
		int e2 = edge_vector.at(i)[1];
		//init the first end point
		GD_Graph[e1]->extendAdjanceNode(GD_Graph[e2]);
		//init the second end point
		GD_Graph[e2]->extendAdjanceNode(GD_Graph[e1]);
	}
	for (int i = 0; i < polygon_vector.size(); i++)
	{
		for (int j = 0; j < polygon_vector.at(i)->getShortCutNumber(); j++)
		{
			int e1 = polygon_vector.at(i)->getShortCut(j).end_point1_index;
			int e2 = polygon_vector.at(i)->getShortCut(j).end_point2_index;
			//init the first end point
			GD_Graph[e1]->extendAdjanceNode(GD_Graph[e2]);
			//init the second end point
			GD_Graph[e2]->extendAdjanceNode(GD_Graph[e1]);
		}
	}
}

void Graph::checkAdjancePair()
{
	for (int i = 0; i < this->node_number; i++)
	{
		cout << GD_Graph[i]->getNodeIndex() << " ";
		GD_Graph[i]->printAdjanceNode();
		cout << endl;
	}

}

void Graph::insertintoVLISTwithAscendant(Node * e)
{
	if (this->vlist_size == 0)
	{
		this->VLIST = new Node*[1];
		this->VLIST[0] = e;
		this->vlist_size = 1;
	}
	else
	{
		Node** new_vlist = new Node*[this->vlist_size + 1];
		int i = 0;
		int j = 0;
		int k = 0;
		while (i < this->vlist_size && this->VLIST[i]->getValue() <= e->getValue())
		{
			i++;
		}
		for (; k < i; j++, k++)
		{
			new_vlist[j] = this->VLIST[k];
		}
		new_vlist[j++] = e;
		for (; k < this->vlist_size; j++, k++)
		{
			new_vlist[j] = this->VLIST[k];
		}
		delete[] this->VLIST;
		this->VLIST = new_vlist;
		this->vlist_size++;
	}
}

void Graph::insertBIandAREA(int index, double area)
{
	if (this->BI_number == 0)
	{
		this->BI = new int[1];
		this->BIarea = new double[1];

		BI[0] = index;
		BIarea[0] = area;
		this->BI_number = 1;
	}
	else
	{
		int* new_BI = new int[this->BI_number + 1];
		double* new_BIarea = new double[this->BI_number + 1];
		for (int i = 0; i < this->BI_number; i++)
		{
			new_BI[i] = this->BI[i];
			new_BIarea[i] = this->BIarea[i];
		}
		new_BI[this->BI_number] = index;
		new_BIarea[this->BI_number] = area;
		delete[] this->BI;
		delete[] this->BIarea;
		this->BI = new_BI;
		this->BIarea = new_BIarea;
		this->BI_number++;
	}
}

double Graph::DijkstraDistanceForBI(int base_vertex)
{
	vector<v2d> VertexcomposeBI;
	VertexcomposeBI.push_back(_v2d_(GD_Graph[base_vertex]->getX(), GD_Graph[base_vertex]->getY()));
	insertintoVLISTwithAscendant(GD_Graph[base_vertex]);
	removeCandidateBI(base_vertex);
	//checkVLIST();
	while (this->vlist_size != 0)
	{
		int v = removeMininVLIST();
		for (int i = 0; i < GD_Graph[v]->getAdjanceNodeNumber(); i++)
		{
			double new_value = GD_Graph[v]->getValue() + returndistance(GD_Graph[v], GD_Graph[v]->getAdjanceNode(i));
			if (new_value < GD_Graph[v]->getAdjanceNode(i)->getValue() && new_value <= this->thread_hole)
			{
				GD_Graph[v]->getAdjanceNode(i)->modifyValue(new_value);
				VertexcomposeBI.push_back(_v2d_(GD_Graph[v]->getAdjanceNode(i)->getX(), GD_Graph[v]->getAdjanceNode(i)->getY()));
				insertintoVLISTwithAscendant(GD_Graph[v]->getAdjanceNode(i));
				removeCandidateBI(GD_Graph[v]->getAdjanceNode(i)->getNodeIndex());
			}
		}
	}
	return CalculateAreawithMelkman(VertexcomposeBI);
}

void Graph::getResult(double Pack[240][2], double x1, double y1, double x2, double y2)
{
	int i, t, tmax;
	double x3, y3, R, Rmax;
	double ResultPack[240][2];
	ResultPack[0][0] = 0;
	if (Pack[0][0] <= 1)
		return;
	x3 = Pack[1][0];
	y3 = Pack[1][1];
	R = x1*y2 + x3*y1 + x2*y3 - x3*y2 - x2*y1 - x1*y3;
	Rmax = R;
	tmax = 1;
	for (i = 2; i <= Pack[0][0]; i++)
	{
		x3 = Pack[i][0];
		y3 = Pack[i][1];
		R = x1*y2 + x3*y1 + x2*y3 - x3*y2 - x2*y1 - x1*y3;
		if (R >= 0)
		{
			t = ++ResultPack[0][0];
			ResultPack[t][0] = x3;
			ResultPack[t][1] = y3;
		}
		if (R > Rmax)
		{
			Rmax = R;
			tmax = i;
		}
	}
	if (Rmax <= 0)
	{
		for (i = 1; i < ResultPack[0][0]; i++)
		{
			x3 = ResultPack[i][0];
			y3 = ResultPack[i][1];
			R = x1*y2 + x3*y1 + x2*y3 - x3*y2 - x2*y1 - x1*y3;
			if (R == 0 && !((x3 == x2&&y3 == y2) || (x3 == x1&&y3 == y1)))
			{
				t = ++g_result[0][0];
				g_result[t][0] = ResultPack[i][0];
				g_result[t][1] = ResultPack[i][1];
			}
		}
		return;
	}
	else
	{
		t = ++g_result[0][0];
		g_result[t][0] = Pack[tmax][0];
		g_result[t][1] = Pack[tmax][1];
		if (ResultPack[0][0] == 0)
			return;
	}
	getResult(ResultPack, x1, y1, Pack[tmax][0], Pack[tmax][1]);
	getResult(ResultPack, Pack[tmax][0], Pack[tmax][1], x2, y2);
}

double Graph::CalculateAreawithMelkman(vector<v2d> node_vector)
{
	double Point[240][2];//Point存所有点。  
	int i = 1, n;
	double x1, y1, x2, y2, x3, y3;
	g_result[0][0] = 0; Point[0][0] = 0;//Point的第一行第一列元素存放包里面有几个点。初始化为0。  
	n = node_vector.size();
	for (i = 1; i <= n; i++)
	{
		Point[i][0] = node_vector.at(i - 1)[0];
		Point[i][1] = node_vector.at(i - 1)[1];
	}
	Point[0][0] = i - 1;
	x1 = Point[1][0];
	y1 = Point[1][1];
	x2 = x1;
	y2 = y1;
	for (i = 2; i <= Point[0][0]; i++)
	{
		x3 = Point[i][0];
		y3 = Point[i][1];
		if (x3 < x1)
		{
			x1 = x3;
			y1 = y3;
		}
		else if (x3 > x2)
		{
			x2 = x3;
			y2 = y3;
		}
	}
	g_result[1][0] = x1;
	g_result[1][1] = y1;
	g_result[2][0] = x2;
	g_result[2][1] = y2;
	g_result[0][0] += 2;
	getResult(Point, x1, y1, x2, y2);
	getResult(Point, x2, y2, x1, y1);
	//for (i = 1; i <= g_result[0][0]; i++)
	//	printf("(%.2d,%.2d)\n", g_result[i][0], g_result[i][1]);

	double s = 0;
	i = 2;
	for (; i <= g_result[0][0] - 1; i++)
		s += det(g_result[1], g_result[i], g_result[i + 1]);
	return 0.5*fabs(s);
}

void Graph::DijkstraDistance(int base_vertex)
{
	insertintoVLISTwithAscendant(GD_Graph[base_vertex]);
	//checkVLIST();
	while (this->vlist_size != 0)
	{
		int v = removeMininVLIST();
		for (int i = 0; i < GD_Graph[v]->getAdjanceNodeNumber(); i++)
		{
			double new_value = GD_Graph[v]->getValue() + returndistance(GD_Graph[v], GD_Graph[v]->getAdjanceNode(i));
			if (new_value < GD_Graph[v]->getAdjanceNode(i)->getValue())
			{
				GD_Graph[v]->getAdjanceNode(i)->modifyValue(new_value);
				insertintoVLISTwithAscendant(GD_Graph[v]->getAdjanceNode(i));
			}
		}
	}
}

void Graph::calculateGD()
{
	//cout << "before BI=============" << endl;
	//checkGD_Node();
	findBI();
	//cout << "after BI=============" << endl;
	//checkGD_Node();
	this->μMatrix = new double*[this->BI_number];
	for (int i = 0; i < this->BI_number; i++)
	{
		this->μMatrix[i] = new double[this->node_number];
		for (int j = 0; j < this->node_number; j++)
		{
			this->μMatrix[i][j] = -1;
			this->GD_Graph[j]->modifyValue(MAX_VALUE);
		}
		this->GD_Graph[this->BI[i]]->modifyValue(0);
		DijkstraDistance(this->BI[i]);
		//checkGD_Node();
		for (int j = 0; j < this->node_number; j++)
		{
			this->μMatrix[i][j] = this->GD_Graph[j]->getValue();
		}
	}

	for (int i = 0; i < this->node_number; i++)
	{
		double μ = 0;
		for (int j = 0; j < this->BI_number; j++)
		{
			μ += this->μMatrix[j][i] * this->BIarea[j];
		}
		this->GD_Graph[i]->modifyFunctionValue(μ);
	}
	//checkFunctionValue();
	function_valueNormalized();
}

void Graph::function_valueNormalized()
{
	double max_value = this->GD_Graph[0]->getFunction_Value();
	double min_value = this->GD_Graph[0]->getFunction_Value();
	for (int i = 0; i < this->node_number; i++)
	{
		double selected_value = this->GD_Graph[i]->getFunction_Value();
		if (selected_value >= max_value)
		{
			max_value = selected_value;
		}
		if (selected_value <= min_value)
		{
			min_value = selected_value;
		}

	}
	//cout << "MAX_VALUE:" << max_value << endl;
	//cout << "MIN_VALUE:" << min_value << endl;
	for (int i = 0; i < this->node_number; i++)
	{
		//double new_N_function_value = 1.0 * (this->GD_Graph[i]->getFunction_Value() - min_value) / (max_value - min_value);
		double new_N_function_value = 1.0 * (this->GD_Graph[i]->getFunction_Value() - min_value) / max_value;
		this->GD_Graph[i]->modifyFunctionValue(new_N_function_value);
	}
}

void Graph::initCandidateBI()
{
	this->candidate_BI = new int[this->node_number];
	for (int i = 0; i < this->node_number; i++)
	{
		this->candidate_BI[i] = i;
	}
	this->candidate_number = this->node_number;
}

void Graph::removeCandidateBI(int vertex)
{
	int index = -1;
	for (int i = 0; i < this->candidate_number; i++)
	{
		if (this->candidate_BI[i] == vertex)
		{
			index = i;
		}
	}
	if (index >= 0)
	{
		int* new_candidate = new int[this->candidate_number];
		for (int i = 0; i < index; i++)
		{
			new_candidate[i] = this->candidate_BI[i];
		}
		for (int i = index + 1; i < this->candidate_number; i++)
		{
			new_candidate[i - 1] = this->candidate_BI[i];
		}
		delete[] this->candidate_BI;
		this->candidate_BI = new_candidate;
		this->candidate_number--;
	}
}

int Graph::getBIfromCandidate()
{

	int base_vertex = (rand() % (this->candidate_number - 0)) + 0;
	return candidate_BI[base_vertex];
}

void Graph::findBI()
{
	srand((unsigned)time(NULL));
	initCandidateBI();
	while (this->candidate_number != 0)
	{
		//cout << "==========================" << endl;
		//cout << "CBI number-B: " << this->candidate_number << endl;
		int base_vertex = getBIfromCandidate();
		//cout << "base_vertex: " << base_vertex << endl;
		GD_Graph[base_vertex]->modifyValue(0);
		double area = DijkstraDistanceForBI(base_vertex);
		//cout << "CBI number-A: " << this->candidate_number << endl;
		//cout << "BIandAREA: " << base_vertex << " " << area << endl;
		insertBIandAREA(base_vertex, area);
	}

	//cout << "CBI number:" << this->candidate_number << endl;
}

void Graph::checkFunctionValue()
{
	cout << " ====FUNCTION_VALUE===="<< endl;
	for (int i = 0; i < this->node_number; i++)
	{
		cout << i << " " << this->GD_Graph[i]->getFunction_Value() << endl;
	}
}

void Graph::checkVLIST()
{
	for (int i = 0; i < this->vlist_size; i++)
	{
		cout << this->VLIST[i]->getValue() << " ";
	}
	cout << endl;
}

void Graph::checkGD_Node()
{
	cout << " ====GD_VALUE====" << endl;
	for (int i = 0; i < this->node_number; i++)
	{
		cout << i << " " << this->GD_Graph[i]->getValue() << endl;
	}
}

int Graph::removeMininVLIST()
{
	int remove_index = this->VLIST[0]->getNodeIndex();
	if (this->vlist_size > 1)
	{
		Node** new_vlist = new Node*[this->vlist_size - 1];
		for (int i = 1, j = 0; i < this->vlist_size; i++, j++)
		{
			new_vlist[j] = this->VLIST[i];
		}
		delete[] this->VLIST;
		this->VLIST = new_vlist;
	}
	else
	{
		delete[] this->VLIST;
	}
	this->vlist_size--;
	return remove_index;
}
