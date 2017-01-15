#include <GL/glut.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <LibIV/libiv.h>
#include<windows.h>
#include "EdgeSet.h"
#include "GraphandNode.h"
#include "Mesh.h"
#include "ElementComparator.h"
#define KRANGHE 16
#define MIN_HEIGHT 250.0
#define MAX_HEIGHT 350.0
#define MIN_LENGTH 250.0
#define MAX_LENGTH 400.0
using namespace std;

typedef std::vector<v2d> conditions_vector;
typedef std::vector<v3i> mesh_points_vector;
typedef std::vector<little_polygon*> polygon_vector;
typedef std::vector<v3i> neighbor_vector;
typedef std::vector<v2i> edge_vector;

conditions_vector NODE;
mesh_points_vector MESH;
polygon_vector POLYGON;
neighbor_vector NEIGHBOR;
edge_vector EDGE;

ElementComparator* ec = NULL;
Graph* a_graph = NULL;
Element* element1 = NULL;
Element* element2 = NULL;
double max_height = 0;
double min_height = 0;
double height = 0;


int INDEX = 0;

double Heron(int a, int b, int c)
{
	v2d point1 = NODE.at(a);
	v2d point2 = NODE.at(b);
	v2d point3 = NODE.at(c);
	double len1 = abs(sqrt(pow((point1[0] - point2[0]), 2) + pow((point1[1] - point2[1]), 2)));
	double len2 = abs(sqrt(pow((point2[0] - point3[0]), 2) + pow((point2[1] - point3[1]), 2)));
	double len3 = abs(sqrt(pow((point3[0] - point1[0]), 2) + pow((point3[1] - point1[1]), 2)));
	double p = 0.5 * (len1 + len2 + len3);
	return sqrt(p*(p - len1)*(p - len2)*(p - len3));
}

double computeArea()
{
	double result = 0;
	for (int i = 0; i < MESH.size(); i++)
	{
		result += Heron(MESH.at(i)[0], MESH.at(i)[1], MESH.at(i)[2]);
	}
	return result;
}

void completeNode(double r)
{
	if (a_graph != NULL)
	{
		delete a_graph;
	}
	a_graph = new Graph(NODE, r);
	a_graph->initAdjancePair(EDGE, POLYGON);
	a_graph->calculateGD();
}

void savedatainfile(string filename)
{
	ofstream f1(filename);
	for (int i = 0; i < POLYGON.size(); i++)
	{
		f1 << i << " ";
		for (int j = 0; j < POLYGON.at(i)->getVertexsNumber(); j++)
		{
			f1 << POLYGON.at(i)->returnVertexsbyIndex(j) << " ";
		}
		f1 << endl;
	}
	f1.close();
}

void initpolygon_vector()
{
	vector <little_polygon*>().swap(POLYGON);
	little_polygon* a = NULL;
	for (int i = 0; i < MESH.size(); i++)
	{
		int index1 = MESH.at(i)[0];
		int index2 = MESH.at(i)[1];
		int index3 = MESH.at(i)[2];

		a = new little_polygon(index1, index2, index3);
		POLYGON.push_back(a);
	}
}

void expansionPolygon()
{
	for (int i = 0; i < POLYGON.size(); i++)
	{
		POLYGON.at(i)->addEdgeByNeighbor(i, MESH, NEIGHBOR);
		POLYGON.at(i)->addShortCut(NODE);
	}
}

void printPolygon()
{
	for (int i = 0; i < POLYGON.size(); i++)
	{
		for (int j = 0; j < POLYGON.at(i)->getVertexsNumber(); j++)
		{
			cout << POLYGON.at(i)->returnVertexsbyIndex(j) << " ";
		}
		cout << endl;
		cout << "========================================" << endl;
	}
}

void getEdgefromFile(string filename)
{
	vector <v2i>().swap(EDGE);
	ifstream im;
	const int LINE_LENGTH = 100;
	im.open(filename);
	if (!im)
	{
		cout << "--open Edge file failed!--" << endl;
	}
	else
	{
		string s;
		//cross the first line
		getline(im, s);
		while (getline(im, s))
		{
			int x = 0;
			int y = 0;
			char * strc = new char[strlen(s.c_str()) + 1];
			strcpy_s(strc, strlen(s.c_str()) + 1, s.c_str());
			char  *p = NULL, *pNext = NULL;

			p = strtok_s(strc, " ", &pNext);
			if (*p == '#')
			{
				continue;
			}

			int edge_segment_index = 1;
			while (p != NULL)
			{
				if (edge_segment_index == 2)
				{
					x = atoi(p);
				}
				else if (edge_segment_index == 3)
				{
					y = atoi(p);
				}
				p = strtok_s(NULL, " ", &pNext);
				edge_segment_index++;
			}
			EDGE.push_back(_v2i_(x, y));
			edge_segment_index = 0;
		}
	}
	im.close();
}

void getNeighborfromFile(string filename)
{
	vector <v3i>().swap(NEIGHBOR);
	ifstream im;
	const int LINE_LENGTH = 100;
	im.open(filename);
	if (!im)
	{
		cout << "--open Neighbour file failed!--" << endl;
	}
	else
	{
		string s;
		//cross the first line
		getline(im, s);
		while (getline(im, s))
		{
			int a1 = 0;
			int a2 = 0;
			int a3 = 0;
			char * strc = new char[strlen(s.c_str()) + 1];
			strcpy_s(strc, strlen(s.c_str()) + 1, s.c_str());
			char  *p = NULL, *pNext = NULL;

			p = strtok_s(strc, " ", &pNext);
			if (*p == '#')
			{
				continue;
			}

			int neighbor_point_index = 1;
			while (p != NULL)
			{
				if (neighbor_point_index == 2)
				{
					a1 = atoi(p);
				}
				else if (neighbor_point_index == 3)
				{
					a2 = atoi(p);
				}
				else if (neighbor_point_index == 4)
				{
					a3 = atoi(p);
				}
				p = strtok_s(NULL, " ", &pNext);
				neighbor_point_index++;
			}
			NEIGHBOR.push_back(_v3i_(a1, a2, a3));
			neighbor_point_index = 0;
		}
	}
	im.close();
}

void getElementfromFile(string filename)
{
	vector <v3i>().swap(MESH);
	ifstream im;
	const int LINE_LENGTH = 100;
	im.open(filename);
	if (!im)
	{
		cout << "--open Element file failed!--" << endl;
	}
	else
	{
		string s;
		//cross the first line
		getline(im, s);
		while (getline(im, s))
		{
			int a1 = 0;
			int a2 = 0;
			int a3 = 0;
			char * strc = new char[strlen(s.c_str()) + 1];
			strcpy_s(strc, strlen(s.c_str()) + 1, s.c_str());
			char  *p = NULL, *pNext = NULL;

			p = strtok_s(strc, " ", &pNext);
			if (*p == '#')
			{
				continue;
			}

			int mesh_point_index = 1;
			while (p != NULL)
			{
				if (mesh_point_index == 2)
				{
					a1 = atoi(p);
				}
				else if (mesh_point_index == 3)
				{
					a2 = atoi(p);
				}
				else if (mesh_point_index == 4)
				{
					a3 = atoi(p);
				}
				p = strtok_s(NULL, " ", &pNext);
				mesh_point_index++;
			}
			MESH.push_back(_v3i_(a1, a2, a3));
			mesh_point_index = 0;
		}
	}
	im.close();
}

void getNodefromFile(string filename)
{
	vector <v2d>().swap(NODE);
	ifstream im;
	const int LINE_LENGTH = 100;
	im.open(filename);
	if (!im)
	{
		cout << "--open Node file failed!--" << endl;
	}
	else
	{
		string s;
		//cross the first line
		getline(im, s);
		while (getline(im, s))
		{
			int x = 0;
			int y = 0;
			char * strc = new char[strlen(s.c_str()) + 1];
			strcpy_s(strc, strlen(s.c_str()) + 1, s.c_str());
			char  *p = NULL, *pNext = NULL;

			p = strtok_s(strc, " ", &pNext);
			if (*p == '#')
			{
				continue;
			}

			int line_segment_index = 1;
			while (p != NULL)
			{
				if (line_segment_index == 2)
				{
					x = atof(p);
				}
				else if (line_segment_index == 3)
				{
					y = atof(p);
				}
				p = strtok_s(NULL, " ", &pNext);
				line_segment_index++;
			}
			NODE.push_back(_v2d_(x, y));
			line_segment_index = 0;
		}
	}
	im.close();
}

void myDisplay()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3d(1.0, 1.0, 1.0);

	glFlush();

	if (element1 != NULL)
	{
		glPointSize(5.0f);
		for (int i = 0; i < element1->returnPoint_Number(); i++)
		{

			double w = 1.0 * (element1->returnPoints(i).value - element1->returnMinFHeight()) / (element1->returnMaxFHeight() - element1->returnMinFHeight());
			glColor3f((1 - w) * 1.0, w * 1.0, 0.0);
			if (element1->returnPoints(i).value == 0)
			{
				glColor3f(0.0, 0.0, 1.0);
			}
			glBegin(GL_POINTS);
			GLdouble x = element1->returnPoints(i).x;
			GLdouble y = element1->returnPoints(i).y;
			glVertex2d(x, y);
			glEnd();
		}
		glFlush();
		glColor3f(1.0, 1.0, 1.0);
		for (int i = 0; i < element1->returnBoundry_Number(); i++)
		{
			double w = 1.0 * (element1->returnBoundrys(i).value - element1->returnMinFHeight()) / (element1->returnMaxFHeight() - element1->returnMinFHeight());
			int index1 = element1->returnBoundrys(i).end_point1_index;
			int index2 = element1->returnBoundrys(i).end_point2_index;

			GLdouble x1 = element1->returnPoints(index1).x;
			GLdouble y1 = element1->returnPoints(index1).y;
			GLdouble x2 = element1->returnPoints(index2).x;
			GLdouble y2 = element1->returnPoints(index2).y;

			glColor3f((1 - w) * 1.0, w * 1.0, 0.0);
			glBegin(GL_LINES);
			glVertex2d(x1, y1);
			glVertex2d(x2, y2);
			glEnd();
			glColor3f(1.0, 1.0, 1.0);


		}
		glFlush();
	}
	glColor3f(1.0, 1.0, 1.0);
	glFlush();
}

int main(int argc, char *argv[])
{
	cout << "Krang: "<< KRANGHE << endl;
	//element 1
	string s = "D:\\MRGTrueFile\\1\\2\\";
	getNodefromFile(s + "node.node");
	getElementfromFile(s + "element.ele");
	getNeighborfromFile(s + "neighbor.neig");
	getEdgefromFile(s + "edge.edge");

	double area = 0;
	area = computeArea();
	cout << "ELEMENT 1" << endl;
	cout << "filename: " << s << endl;
	cout << "area before partition:"<< area << endl;
	initpolygon_vector();
	expansionPolygon();
	completeNode(sqrt(0.05 * area));

	element1 = new Element(NODE, MESH, NEIGHBOR, EDGE);
	element1->elementinit(KRANGHE, a_graph);
	element1->elementpartition();
	element1->remeshelement();
	element1->rebuildAdjanceMesh();
	element1->buildFinestMRG();

	//element 2
	string ss = "D:\\MRGTrueFile\\1\\3\\";
	getNodefromFile(ss + "node.node");
	getElementfromFile(ss + "element.ele");
	getNeighborfromFile(ss + "neighbor.neig");
	getEdgefromFile(ss + "edge.edge");

	area = computeArea();
	cout << "ELEMENT 2" << endl;
	cout << "filename: " << ss << endl;
	cout << "area before partition:" << area << endl;
	initpolygon_vector();
	expansionPolygon();
	completeNode(sqrt(0.05 * area));

	element2 = new Element(NODE, MESH, NEIGHBOR, EDGE);
	element2->elementinit(KRANGHE, a_graph);
	element2->elementpartition();
	element2->remeshelement();
	element2->rebuildAdjanceMesh();
	element2->buildFinestMRG();

	ec = new ElementComparator(element1, element2);
	double result = ec->getCompareResult();
	cout <<"compare result: "<< result << endl;

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(600, 600);
	glutCreateWindow("第一个OpenGL程序");

	glutDisplayFunc(&myDisplay);
	gluOrtho2D(MIN_LENGTH, MAX_LENGTH, MIN_HEIGHT, MAX_HEIGHT);
	glutMainLoop();
	return 0;
}
