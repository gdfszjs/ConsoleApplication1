#include "EdgeSet.h"

little_polygon::little_polygon()
{
	this->vertexs = NULL;
	this->numberofvertexs = 0;
	this->polygon_edge = NULL;
	this->numberofpolygon_edge = 0;
	this->short_cut = NULL;
	this->numberofshort_cut = 0;

}

little_polygon::little_polygon(int p1, int p2, int p3)
{
	//init vertexs
	this->vertexs = (int*)malloc(sizeof(int) * 3);
	this->vertexs[0] = p1;
	this->vertexs[1] = p2;
	this->vertexs[2] = p3;
	this->numberofvertexs = 3;

	//init exit edge
	this->exitedge = (edge*)malloc(sizeof(edge) * 3);
	this->exitedge[0].end_point1_index = p2;
	this->exitedge[0].end_point2_index = p3;
	this->exitedge[1].end_point1_index = p1;
	this->exitedge[1].end_point2_index = p3;
	this->exitedge[2].end_point1_index = p1;
	this->exitedge[2].end_point2_index = p2;

	//init polygon edges
	initedge();

	//init short cuts
	this->short_cut = NULL;
	this->numberofshort_cut = 0;
}

little_polygon::~little_polygon()
{
	delete[] this->polygon_edge;
	delete[] this->vertexs;
	delete[] this->short_cut;
	this->vertexs = NULL;
	this->numberofvertexs = 0;
	this->polygon_edge = NULL;
	this->numberofpolygon_edge = 0;
	this->short_cut = NULL;
	this->numberofshort_cut = 0;
}

void little_polygon::initedge()
{
	this->polygon_edge = (edge*)malloc(sizeof(edge) * this->numberofvertexs);
	for (int i = 0; i < this->numberofvertexs; i++)
	{
		if (i == this->numberofvertexs - 1)
		{
			this->polygon_edge[i].end_point1_index = this->vertexs[i];
			this->polygon_edge[i].end_point2_index = this->vertexs[0];
		}
		else 
		{
			this->polygon_edge[i].end_point1_index = this->vertexs[i];
			this->polygon_edge[i].end_point2_index = this->vertexs[i + 1];
		}
	}
	this->numberofpolygon_edge = this->numberofvertexs;
}

void little_polygon::initshortcut(vector<v2d> short_cut_vector)
{
	this->short_cut = (edge*)malloc(sizeof(edge) * short_cut_vector.size());
	for (int i = 0; i < short_cut_vector.size(); i++)
	{
		this->short_cut[i].end_point1_index = short_cut_vector.at(i)[0];
		this->short_cut[i].end_point2_index = short_cut_vector.at(i)[1];
	}
	this->numberofshort_cut = short_cut_vector.size();
}

void little_polygon::addShortCut(vector<v2d> vertexs_vector)
{
	vector<v2d> short_cut_vector;
	for (int i = 0; i < this->numberofvertexs; i++)
	{
		int candidate1 = this->vertexs[i];
		for (int j = i + 1; j < this->numberofvertexs; j++)
		{
			int candidate2 = this->vertexs[j];
			if (!isBoundary(candidate1, candidate2))
			{
				if (isShortCut(candidate1, candidate2, i, j, vertexs_vector))
				{
					short_cut_vector.push_back(_v2d_(candidate1, candidate2));
				}
			}
		}
	}
	initshortcut(short_cut_vector);
}

int little_polygon::intersect2(v2d u1, v2d u2, v2d v1, v2d v2) {

	return((max(u1[0], u2[0]) >= min(v1[0], v2[0])) &&

		(max(v1[0], v2[0]) >= min(u1[0], u2[0])) &&

		(max(u1[1], u2[1]) >= min(v1[1], v2[1])) &&

		(max(v1[1], v2[1]) >= min(u1[1], u2[1])) &&		
		(multiply(v1, u2, u1)*multiply(u2, v2, u1) >= 0) &&

		(multiply(u1, v2, v1)*multiply(v2, u2, v1) >= 0)); 
}

bool little_polygon::intersect(v2d aa, v2d bb, v2d cc, v2d dd)
{
	double delta = determinant(bb[0] - aa[0], cc[0] - dd[0], bb[1] - aa[1], cc[1] - dd[1]);
	if (delta <= (1e-6) && delta >= -(1e-6))  // delta=0，表示两线段重合或平行  
	{
		return false;
	}
	double namenda = determinant(cc[0] - aa[0], cc[0] - dd[0], cc[1] - aa[1], cc[1] - dd[1]) / delta;
	if (namenda>1 || namenda<0)
	{
		return false;
	}
	double miu = determinant(bb[0] - aa[0], cc[0] - aa[0], bb[1] - aa[1], cc[1] - aa[1]) / delta;
	if (miu>1 || miu<0)
	{
		return false;
	}
	return true;
}

bool little_polygon::isShortCut(int c1, int c2,int c1_index,int c2_index, vector<v2d> vertexs_vector)
{
	double ax = vertexs_vector[c1][0];
	double ay = vertexs_vector[c1][1];
	double bx = vertexs_vector[c2][0];
	double by = vertexs_vector[c2][1];

	double c1x = vertexs_vector[this->exitedge[c1_index].end_point1_index][0];
	double c1y = vertexs_vector[this->exitedge[c1_index].end_point1_index][1];
	double d1x = vertexs_vector[this->exitedge[c1_index].end_point2_index][0];
	double d1y = vertexs_vector[this->exitedge[c1_index].end_point2_index][1];

	double c2x = vertexs_vector[this->exitedge[c2_index].end_point1_index][0];
	double c2y = vertexs_vector[this->exitedge[c2_index].end_point1_index][1];
	double d2x = vertexs_vector[this->exitedge[c2_index].end_point2_index][0];
	double d2y = vertexs_vector[this->exitedge[c2_index].end_point2_index][1];

	v2d aa = _v2d_(ax,ay);
	v2d bb = _v2d_(bx, by);

	v2d cc1 = _v2d_(c1x, c1y);
	v2d dd1 = _v2d_(d1x, d1y);

	v2d cc2 = _v2d_(c2x, c2y);
	v2d dd2 = _v2d_(d2x, d2y);

	if (intersect2(aa,bb,cc1,dd1) != 0)
	{
		if (intersect2(aa,bb,cc2,dd2) != 0)
		{
			return true;
		}
		else return false;
	}
	else return false;
	
	
}

bool little_polygon::isBoundary(int c1, int c2)
{
	for (int i = 0; i < this->numberofpolygon_edge; i++)
	{
		int e1 = this->polygon_edge[i].end_point1_index;
		int e2 = this->polygon_edge[i].end_point2_index;
		if ((c1 == e1 && c2 == e2) || (c1 == e2 && c2 == e1))
		{
			return true;
		}
	}
	return false;
}

void little_polygon::addEdgeByNeighbor(int index, vector<v3i>mesh_vector, vector<v3i> neighbor_vector)
{
	vector<int> new_vertexs_vector;
	int end_point1 = 0;
	int end_point2 = 0;
	for (int i = 0; i < this->numberofvertexs; i++)
	{
		end_point1 = this->vertexs[i];
		if (i == this->numberofvertexs - 1)
		{
			end_point2 = this->vertexs[0];
		}
		else end_point2 = this->vertexs[i + 1];

		//cout << "one edge: "<< end_point1 << " " << end_point2 << endl;

		int selected_mesh_index = 0;
		for (int j = 0; j < 3; j++)
		{
			selected_mesh_index = neighbor_vector.at(index)[j];
			if (selected_mesh_index != -1)
			{
				if (isPointonTriangle(end_point1, selected_mesh_index, mesh_vector) && isPointonTriangle(end_point2, selected_mesh_index, mesh_vector))
				{
					int new_vertex_index = getAnotherPointOnTriangle(end_point1, end_point2, selected_mesh_index, mesh_vector);
					new_vertexs_vector.push_back(end_point1);
					new_vertexs_vector.push_back(end_point2);
					new_vertexs_vector.push_back(new_vertex_index);
				}
			}
		}
	}
	modifypolygon(new_vertexs_vector);
	//initedge();
}

bool little_polygon::isPointonTriangle(int vertex_index, int mesh_index, vector<v3i> mesh_vector)
{
	for (int i = 0; i < 3; i++)
	{
		if (vertex_index == mesh_vector.at(mesh_index)[i])
		{
			return true;
		}
	}
	return false;
}

int little_polygon::getAnotherPointOnTriangle(int vertex1_index, int vertex2_index, int mesh_index, vector<v3i> mesh_vector)
{
	for (int i = 0; i < 3; i++)
	{
		if (mesh_vector.at(mesh_index)[i] != vertex1_index && mesh_vector.at(mesh_index)[i] != vertex2_index)
		{
			return mesh_vector.at(mesh_index)[i];
		}
	}
}

void little_polygon::modifypolygon(vector<int> new_vertexs_vector)
{
	
	for (int i = 0; i < new_vertexs_vector.size(); i = i + 3)
	{
		int start_vertex = new_vertexs_vector.at(i);
		int end_vertex = new_vertexs_vector.at(i + 1);
		int new_vertex = new_vertexs_vector.at(i + 2);
		int* new_vertex_vector = (int*)malloc(sizeof(int) * (this->numberofvertexs + 1));
		edge* new_exit_edge_vector = (edge*)malloc(sizeof(edge) * (this->numberofvertexs + 1));
		edge* new_polygon_edge_vector = (edge*)malloc(sizeof(edge) * (this->numberofpolygon_edge + 2));
		for (int j = 0; j < this->numberofpolygon_edge; j++)
		{
			new_polygon_edge_vector[j] = this->polygon_edge[j];
		}
		for (int j = 0,k = 0; j < this->numberofvertexs; j++,k++)
		{
			if (this->vertexs[j] == start_vertex)
			{
				new_vertex_vector[k] = this->vertexs[j];
				new_exit_edge_vector[k] = this->exitedge[j];

				new_vertex_vector[k + 1] = new_vertex;
				new_exit_edge_vector[k + 1].end_point1_index = start_vertex;
				new_exit_edge_vector[k + 1].end_point2_index = end_vertex;

				new_polygon_edge_vector[this->numberofpolygon_edge].end_point1_index = start_vertex;
				new_polygon_edge_vector[this->numberofpolygon_edge].end_point2_index = new_vertex;
				new_polygon_edge_vector[this->numberofpolygon_edge + 1].end_point1_index = end_vertex;
				new_polygon_edge_vector[this->numberofpolygon_edge + 1].end_point2_index = new_vertex;
				k++;
			}
			else
			{
				new_vertex_vector[k] = this->vertexs[j];
				new_exit_edge_vector[k] = this->exitedge[j];
			}
		}
		this->numberofpolygon_edge += 2;
		this->numberofvertexs++;
		delete this->vertexs;
		delete this->exitedge;
		delete this->polygon_edge;
		this->polygon_edge = new_polygon_edge_vector;
		this->vertexs = new_vertex_vector;
		this->exitedge = new_exit_edge_vector;
	}
}
