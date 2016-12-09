#include <GL/glut.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <LibIV/libiv.h>
#include<windows.h>
#include "EdgeSet.h"
using namespace std;

typedef std::vector<v2d> conditions_vector;
typedef std::vector<v3i> mesh_points_vector;
typedef std::vector<little_polygon*> polygon_vector;
typedef std::vector<v3i> neighbor_vector;

conditions_vector NODE;
mesh_points_vector MESH;
polygon_vector POLYGON;
neighbor_vector NEIGHBOR;

int INDEX = 0;
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

void getNeighborfromFile(string filename)
{
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
}

void getElementfromFile(string filename)
{
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
			//cout << s << endl;
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
}

void getNodefromFile(string filename)
{
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
			//cout << s << endl;
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
}

void myDisplay()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glPointSize(5.0f);
	glBegin(GL_POINTS);
	for (int i = 0; i < NODE.size(); i++)
	{
		GLdouble x = NODE.at(i)[0];
		GLdouble y = NODE.at(i)[1];

		glVertex2d(x, y);

	}
	glEnd();

	for (int i = 0; i < MESH.size(); i++)
	{
		int index1 = MESH.at(i)[0];
		int index2 = MESH.at(i)[1];
		int index3 = MESH.at(i)[2];

		GLdouble x1 = NODE.at(index1)[0];
		GLdouble y1 = NODE.at(index1)[1];
		GLdouble x2 = NODE.at(index2)[0];
		GLdouble y2 = NODE.at(index2)[1];
		GLdouble x3 = NODE.at(index3)[0];
		GLdouble y3 = NODE.at(index3)[1];
		glBegin(GL_LINE_LOOP);
		glVertex2d(x1, y1);
		glVertex2d(x2, y2);
		glVertex2d(x3, y3);
		glEnd();
	}
	glFlush();
	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINE_LOOP);
	if (INDEX >= POLYGON.size())
	{
		INDEX = 0;
	}
	int vertex_number = POLYGON.at(INDEX)->getVertexsNumber();
	for (int i = 0; i < vertex_number; i++)
	{
		int i1 = i % vertex_number;
		int i2 = (i + 1) % vertex_number;
		GLdouble x1 = NODE.at(POLYGON.at(INDEX)->returnVertexsbyIndex(i1))[0];
		GLdouble y1 = NODE.at(POLYGON.at(INDEX)->returnVertexsbyIndex(i1))[1];
		GLdouble x2 = NODE.at(POLYGON.at(INDEX)->returnVertexsbyIndex(i2))[0];
		GLdouble y2 = NODE.at(POLYGON.at(INDEX)->returnVertexsbyIndex(i2))[1];
		glVertex2d(x1, y1);
		glVertex2d(x2, y2);
		
	}
	glEnd();
	glFlush();

	glColor3f(1.0, 0.0, 0.0);
	for (int i = 0; i < POLYGON.at(INDEX)->getShortCutNumber(); i++)
	{
		glBegin(GL_LINES);
		GLdouble x1 = NODE.at(POLYGON.at(INDEX)->getShortCut(i).end_point1_index)[0];
		GLdouble y1 = NODE.at(POLYGON.at(INDEX)->getShortCut(i).end_point1_index)[1];
		GLdouble x2 = NODE.at(POLYGON.at(INDEX)->getShortCut(i).end_point2_index)[0];
		GLdouble y2 = NODE.at(POLYGON.at(INDEX)->getShortCut(i).end_point2_index)[1];
		glVertex2d(x1, y1);
		glVertex2d(x2, y2);
		glEnd();
		glFlush();
	}
	INDEX++;

	glColor3d(1.0, 1.0, 1.0);
	glFlush();
	//Sleep(1000);
	glutPostRedisplay();
}

int main(int argc, char *argv[])
{
	string s = "D:\\MRGTrueFile\\1\\3\\";
	getNodefromFile(s + "node.node");
	getElementfromFile(s + "element.ele");
	getNeighborfromFile(s + "neighbor.neig");
	initpolygon_vector();
	expansionPolygon();
	//printPolygon();
	savedatainfile("D:\\polygoncheck.txt");
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(1000, 1000);
	glutCreateWindow("第一个OpenGL程序");
	glutDisplayFunc(&myDisplay);
	gluOrtho2D(0.0, 600.0, 0.0, 600.0);
	glutMainLoop();
	return 0;
}
