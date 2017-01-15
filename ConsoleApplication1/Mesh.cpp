#include "Mesh.h"
#define MAX_VALUE -10000
#define MIN_VALUE 10000
#define OMIGA 0.5
Element::Element()
{
	this->points_vector = 0;
	this->point_number = -1;

	this->meshs_vector = 0;
	this->mesh_number = -1;

	this->height = 0;
	this->max_height = MAX_VALUE;
	this->min_height = MIN_VALUE;
	this->max_f_value = MAX_VALUE;
	this->min_f_value = MIN_VALUE;

	this->K_range = 0;
}

Element::Element(vector<v2d> points_vector, vector<v3i> meshs_vector, vector<v3i> neighbors_vector, vector<v2i> edge_vector)
{
	this->out.open("D:\\out.txt");
	//points
	this->points_vector = new vertex[points_vector.size()];
	for (int i = 0; i < points_vector.size(); i++)
	{
		this->points_vector[i].x = points_vector.at(i)[0];
		this->points_vector[i].y = points_vector.at(i)[1];
		this->points_vector[i].value = -1;
	}
	this->point_number = points_vector.size();
	//boundarys
	this->boundarys_vector = new boundry[edge_vector.size()];
	for (int i = 0; i < edge_vector.size(); i++)
	{
		this->boundarys_vector[i].end_point1_index = edge_vector.at(i)[0];
		this->boundarys_vector[i].end_point2_index = edge_vector.at(i)[1];
	}
	this->boundary_number = edge_vector.size();
	//meshs
	this->meshs_vector = new Mesh[meshs_vector.size()];
	for (int i = 0; i < meshs_vector.size(); i++)
	{
		//vertex
		this->meshs_vector[i].vertexs = new int[3];
		for (int j = 0; j < 3; j++)
		{
			this->meshs_vector[i].vertexs[j] = meshs_vector.at(i)[j];
		}
		//neighbor
		this->meshs_vector[i].adjance_meshs = new int[3];
		for (int j = 0; j < 3; j++)
		{
			this->meshs_vector[i].adjance_meshs[j] = neighbors_vector.at(i)[j];
		}
		//edge
		this->meshs_vector[i].boundrys_vector = new int[3];
		for (int j = 0; j < 3; j++)
		{
			int e1 = this->meshs_vector[i].vertexs[j];
			int e2 = 0;
			if (j == 2) e2 = this->meshs_vector[i].vertexs[0];
			else e2 = this->meshs_vector[i].vertexs[j + 1];
			for (int k = 0; k < this->boundary_number; k++)
			{
				boundry selected_bondry = this->boundarys_vector[k];
				int e3 = selected_bondry.end_point1_index;
				int e4 = selected_bondry.end_point2_index;
				if (e1 == e3 && e2 == e4)
				{
					this->meshs_vector[i].boundrys_vector[j] = k;
				}
				else if (e2 == e3 && e1 == e4)
				{
					this->meshs_vector[i].boundrys_vector[j] = k;
				}
			}
		}
		//area
		this->meshs_vector[i].area = -1;
		this->meshs_vector[i].range = -1;
	}
	this->mesh_number = meshs_vector.size();

	this->height = 0;
	this->max_height = MAX_VALUE;
	this->min_height = MIN_VALUE;
	this->max_f_value = MAX_VALUE;
	this->min_f_value = MIN_VALUE;

	this->K_range = 0;
}

Element::~Element()
{
	if (this->meshs_vector != 0) delete[] this->meshs_vector;
	if (this->points_vector != 0) delete[] this->points_vector;
}

double Element::Heron(int a, int b, int c)
{
	v2d point1 = _v2d_(this->points_vector[a].x, this->points_vector[a].y);
	v2d point2 = _v2d_(this->points_vector[b].x, this->points_vector[b].y);
	v2d point3 = _v2d_(this->points_vector[c].x, this->points_vector[c].y);
	double len1 = abs(sqrt(pow((point1[0] - point2[0]), 2) + pow((point1[1] - point2[1]), 2)));
	double len2 = abs(sqrt(pow((point2[0] - point3[0]), 2) + pow((point2[1] - point3[1]), 2)));
	double len3 = abs(sqrt(pow((point3[0] - point1[0]), 2) + pow((point3[1] - point1[1]), 2)));
	double p = 0.5 * (len1 + len2 + len3);
	return sqrt(p*(p - len1)*(p - len2)*(p - len3));
}

void Element::computeArea()
{
	double result = 0;
	for (int i = 0; i < this->mesh_number; i++)
	{
		this->meshs_vector[i].area = Heron(this->meshs_vector[i].vertexs[0], this->meshs_vector[i].vertexs[1], this->meshs_vector[i].vertexs[2]);
		result += this->meshs_vector[i].area;
	}
	
	cout << "area after partition:"<< result << endl;
}

void Element::elementpartition()
{
	double upper_limit = 0;
	double lower_limit = 0;
	for (int i = 0; i < this->mesh_number; i++)
	{
		int range_number[3] = { -1,-1,-1 };
		double v1_h = this->points_vector[this->meshs_vector[i].vertexs[0]].value;
		double v2_h = this->points_vector[this->meshs_vector[i].vertexs[1]].value;
		double v3_h = this->points_vector[this->meshs_vector[i].vertexs[2]].value;
		for (int j = 0; j < this->K_range; j++)
		{
			double l_h = this->height_scanning_line[j];
			double h_h = this->height_scanning_line[j + 1];
			if (v1_h >= l_h && v1_h <= h_h)
			{
				if (v1_h == l_h)
				{
					this->points_vector[this->meshs_vector[i].vertexs[0]].range = j;
				}
				else if (v1_h == h_h)
				{
					this->points_vector[this->meshs_vector[i].vertexs[0]].range = j + 1;
				}
				else
				{
					this->points_vector[this->meshs_vector[i].vertexs[0]].range = j;
				}
			}
			if (v2_h >= l_h && v2_h <= h_h)
			{
				if (v2_h == l_h)
				{
					this->points_vector[this->meshs_vector[i].vertexs[1]].range = j;
				}
				else if (v2_h == h_h)
				{
					this->points_vector[this->meshs_vector[i].vertexs[1]].range = j + 1;
				}
				else
				{
					this->points_vector[this->meshs_vector[i].vertexs[1]].range = j;
				}
			}
			if (v3_h >= l_h && v3_h <= h_h)
			{
				if (v3_h == l_h)
				{
					this->points_vector[this->meshs_vector[i].vertexs[2]].range = j;
				}
				else if (v3_h == h_h)
				{
					this->points_vector[this->meshs_vector[i].vertexs[2]].range = j + 1;
				}
				else
				{
					this->points_vector[this->meshs_vector[i].vertexs[2]].range = j;
				}
			}

			if (v1_h >= l_h && v1_h <= h_h)
			{
				if (v2_h >= l_h && v2_h <= h_h)
				{
					if (v3_h >= l_h && v3_h <= h_h)
					{
						this->meshs_vector[i].range = j;
					}
				}
			}


		}
		if (this->meshs_vector[i].range == -1)
		{
			this->fixing_mesh_index_vector.push_back(i);
		}
		/*if (range_number[0] != range_number[1] != range_number[2])
		{

		}
		for (int j = 0; j < 3; j++)
		{
			this->points_vector[this->meshs_vector[i].vertexs[j]].range = range_number[j];
		}*/


	}
	checkVertexRange();
	checkMesh();
}

void Element::elementinit(int K, Graph* e)
{
	getVertexsValue(e);
	calculateBoundryValue();
	setK_range(K);
	getHeight();
	out << "element init done!" << endl;
}

void Element::getVertexsValue(Graph* e)
{
	for (int i = 0; i < this->point_number; i++)
	{
		this->points_vector[i].value = e->returnG_Graph(i)->getFunction_Value();
		if (this->points_vector[i].value >= this->max_f_value)
		{
			this->max_f_value = this->points_vector[i].value;
		}
		if (this->points_vector[i].value <= this->min_f_value)
		{
			this->min_f_value = this->points_vector[i].value;
		}
	}
}

void Element::calculateBoundryValue()
{
	double w = 0.5;
	for (int i = 0; i < this->boundary_number; i++)
	{
		int v1 = this->boundarys_vector[i].end_point1_index;
		int v2 = this->boundarys_vector[i].end_point2_index;
		double f_v1 = this->points_vector[v1].value;
		double f_v2 = this->points_vector[v2].value;
		this->boundarys_vector[i].value = 1.0 * OMIGA * f_v1 + 1.0 * (1 - OMIGA) * f_v2;
	}
}

void Element::getHeight()
{
	for (int i = 0; i < this->point_number; i++)
	{
		double new_value = this->points_vector[i].value;
		if (new_value >= this->max_height)
		{
			this->max_height = new_value;
		}
		if (new_value <= this->min_height)
		{
			this->min_height = new_value;
			this->min_value_index = i;
		}
	}
	this->height = this->max_height - this->min_height;

	this->height_scanning_line = new double[this->K_range + 1];
	for (int i = 0; i <= this->K_range; i++)
	{
		this->height_scanning_line[i] = this->min_height + (i * 1.0 * this->height / this->K_range);
	}
	checkHeight();
}

void Element::remeshelement()
{

	out << "before fixing: " << this->mesh_number << endl;
	out << "fixing number: " << this->fixing_mesh_index_vector.size() << endl;
	out << "====EDGE VECTOR BEFORE FIX====\n";
	checkEdge();
	vector<vector<int>> fixed_edge;
	for (int i = 0; i < this->fixing_mesh_index_vector.size(); i++)
	{
		//get new points
		vector<int> point_on_line1;//the point on the edge(lowest to middle)
		vector<int> point_on_line2;//the point on the edge(middle to higest)
		vector<int> point_on_line3;//the point on the edge(lowest to higest)

		//get the mesh need to fixed
		Mesh editing_mesh = this->meshs_vector[this->fixing_mesh_index_vector.at(i)];
		double max_value = MAX_VALUE;
		int max_index = -1;
		double min_value = MIN_VALUE;
		int min_index = -1;
		int new_number = 0;
		int node_edge_flag = -1;
		for (int j = 0; j < 3; j++)
		{
			if (this->points_vector[editing_mesh.vertexs[j]].value >= max_value)
			{
				max_value = this->points_vector[editing_mesh.vertexs[j]].value;
				max_index = j;
			}
			if (this->points_vector[editing_mesh.vertexs[j]].value <= min_value)
			{
				min_value = this->points_vector[editing_mesh.vertexs[j]].value;
				min_index = j;
			}
		}
		//sort the point by the value
		int point[3] = { editing_mesh.vertexs[min_index],editing_mesh.vertexs[3 - (min_index + max_index)],editing_mesh.vertexs[max_index] };
		out << "====POINT ON OLD MESH====" << endl;
		for (int j = 0; j < 3; j++)
		{
			checkVertex(point[j]);
		}
		//split the old edge into some new edge
		//
		//lowest to middle[0]~[1]
		vector<int> one_fixed_edge;
		node_edge_flag = isLineinFixedEdgevector(fixed_edge, point[0], point[1]);
		//if it is a new edge(not in the fixed_edge_vector)
		if (node_edge_flag == -1)
		{
			one_fixed_edge.push_back(point[0]);
			one_fixed_edge.push_back(point[1]);

			point_on_line1.push_back(point[0]);
			new_number = this->points_vector[point[1]].range - this->points_vector[point[0]].range;
			if (this->points_vector[point[1]].value - this->height_scanning_line[this->points_vector[point[1]].range] == 0)
			{
				//if the end point is acutly on the range,the new node number need minus 1;
				new_number--;
			}
			for (int j = 0; j < new_number; j++)
			{
				double new_value = this->height_scanning_line[this->points_vector[point[0]].range + j + 1];
				int new_range = this->points_vector[point[0]].range + j + 1;
				double new_w = 1.0 * (new_value - this->points_vector[point[1]].value) / (this->points_vector[point[0]].value - this->points_vector[point[1]].value);
				double new_x = this->points_vector[point[0]].x * new_w + this->points_vector[point[1]].x * (1 - new_w);
				double new_y = this->points_vector[point[0]].y * new_w + this->points_vector[point[1]].y * (1 - new_w);
				insertNewVertex(new_x, new_y, new_value, new_range);
				one_fixed_edge.push_back(this->point_number - 1);
				point_on_line1.push_back(this->point_number - 1);
			}
			point_on_line1.push_back(point[1]);
			//the line has been fixed,need to change the edge vector
			reorganizeEdge(point_on_line1.front(), point_on_line1.back(), point_on_line1);
			fixed_edge.push_back(one_fixed_edge);
		}
		else
		{
			point_on_line1.push_back(point[0]);
			for (int j = 2; j < fixed_edge.at(node_edge_flag).size(); j++)
			{
				point_on_line1.push_back(fixed_edge.at(node_edge_flag).at(j));
			}
			point_on_line1.push_back(point[1]);
		}
		//
		//middle to higest[1]~[2]
		vector<int> second_fixed_edge;
		node_edge_flag = isLineinFixedEdgevector(fixed_edge, point[1], point[2]);
		//if it is a new edge(not in the fixed_edge_vector)
		if (node_edge_flag == -1)
		{
			second_fixed_edge.push_back(point[1]);
			second_fixed_edge.push_back(point[2]);

			point_on_line2.push_back(point[1]);
			new_number = this->points_vector[point[2]].range - this->points_vector[point[1]].range;
			if (this->points_vector[point[2]].value - this->height_scanning_line[this->points_vector[point[2]].range] == 0)
			{
				//if the end point is acutly on the range,the new node number need minus 1;
				new_number--;
			}
			for (int j = 0; j < new_number; j++)
			{
				double new_value = this->height_scanning_line[this->points_vector[point[1]].range + j + 1];
				int new_range = this->points_vector[point[1]].range + j + 1;
				double new_w = 1.0 * (new_value - this->points_vector[point[2]].value) / (this->points_vector[point[1]].value - this->points_vector[point[2]].value);
				double new_x = this->points_vector[point[1]].x * new_w + this->points_vector[point[2]].x * (1 - new_w);
				double new_y = this->points_vector[point[1]].y * new_w + this->points_vector[point[2]].y * (1 - new_w);
				insertNewVertex(new_x, new_y, new_value, new_range);
				second_fixed_edge.push_back(this->point_number - 1);
				point_on_line2.push_back(this->point_number - 1);
			}
			point_on_line2.push_back(point[2]);
			//the line has been fixed,need to change the edge vector
			reorganizeEdge(point_on_line2.front(), point_on_line2.back(), point_on_line2);
			fixed_edge.push_back(second_fixed_edge);
		}
		else
		{
			point_on_line2.push_back(point[1]);
			for (int j = 2; j < fixed_edge.at(node_edge_flag).size(); j++)
			{
				point_on_line2.push_back(fixed_edge.at(node_edge_flag).at(j));
			}
			point_on_line2.push_back(point[2]);
		}
		//
		//lowest to higest[0]~[2]
		vector<int> third_fixed_edge;
		node_edge_flag = isLineinFixedEdgevector(fixed_edge, point[0], point[2]);
		//if it is a new edge(not in the fixed_edge_vector)
		if (node_edge_flag == -1)
		{
			third_fixed_edge.push_back(point[0]);
			third_fixed_edge.push_back(point[2]);

			point_on_line3.push_back(point[0]);
			new_number = this->points_vector[point[2]].range - this->points_vector[point[0]].range;
			if (this->points_vector[point[2]].value - this->height_scanning_line[this->points_vector[point[2]].range] == 0)
			{
				//if the end point is acutly on the range,the new node number need minus 1;
				new_number--;
			}
			for (int j = 0; j < new_number; j++)
			{
				double new_value = this->height_scanning_line[this->points_vector[point[0]].range + j + 1];
				int new_range = this->points_vector[point[0]].range + j + 1;
				double new_w = 1.0 * (new_value - this->points_vector[point[2]].value) / (this->points_vector[point[0]].value - this->points_vector[point[2]].value);
				double new_x = this->points_vector[point[0]].x * new_w + this->points_vector[point[2]].x * (1 - new_w);
				double new_y = this->points_vector[point[0]].y * new_w + this->points_vector[point[2]].y * (1 - new_w);
				insertNewVertex(new_x, new_y, new_value, new_range);
				third_fixed_edge.push_back(this->point_number - 1);
				point_on_line3.push_back(this->point_number - 1);
			}
			point_on_line3.push_back(point[2]);
			//the line has been fixed,need to change the edge vector
			reorganizeEdge(point_on_line3.front(), point_on_line3.back(), point_on_line3);
			fixed_edge.push_back(third_fixed_edge);
		}
		else
		{
			point_on_line3.push_back(point[0]);
			for (int j = 2; j < fixed_edge.at(node_edge_flag).size(); j++)
			{
				point_on_line3.push_back(fixed_edge.at(node_edge_flag).at(j));
			}
			point_on_line3.push_back(point[2]);
		}
		
		//get the new point on old mesh
		vector<int> left_edge_vector;//[0]~[1]~[2]
		vector<int> right_edge_vector;//[0]~[2]
		//init left_edge_vector
		for (int j = 0; j < point_on_line1.size(); j++)
		{
			left_edge_vector.push_back(point_on_line1.at(j));
		}
		for (int j = 1; j < point_on_line2.size(); j++)
		{
			left_edge_vector.push_back(point_on_line2.at(j));
		}
		//init right_edge_vector
		out << "======NEW POINT ON OLD MESH======" << endl;
		for (int j = 0; j < point_on_line3.size(); j++)
		{
			right_edge_vector.push_back(point_on_line3.at(j));
		}
		//output point on mesh
		out << "left_edge_vector: ";
		for (int j = 0; j < left_edge_vector.size(); j++)
		{
			out << left_edge_vector.at(j) << " ";
		}
		out << endl;
		out << "right_edge_vector:";
		for (int j = 0; j < right_edge_vector.size(); j++)
		{
			out << right_edge_vector.at(j) << " ";
		}
		out << endl;
		//start to remesh the old mesh
		int end_flag = 0;
		int si1 = 0;
		int si2 = 0;
		int ei1 = si1 + 1;
		int ei2 = si2 + 1;
		vector<int> temp_mesh_vector;
		temp_mesh_vector.push_back(this->fixing_mesh_index_vector.at(i));
		while (!end_flag)
		{
			int start_index1 = left_edge_vector.at(si1);
			int start_index2 = right_edge_vector.at(si2);
			int end_index1 = left_edge_vector.at(ei1);
			int end_index2 = right_edge_vector.at(ei2);
			out << "======ONE SCAN AREA======" << endl;
			out << "start_index1:" << start_index1 << " " << this->points_vector[start_index1].value << endl;
			out << "start_index2:" << start_index2 << " " << this->points_vector[start_index2].value << endl;
			out << "end_index1:" << end_index1 << " " << this->points_vector[end_index1].value << endl;
			out << "end_index2:" << end_index2 << " " << this->points_vector[end_index2].value << endl;
			//the first new mesh
			if (start_index1 == start_index2)
			{
				//push back a new mesh
				out << "======ONE EXTEND======" << endl;
				out << "A1:" << start_index1 << " " << this->points_vector[start_index1].value << endl;
				out << "A2:" << end_index1 << " " << this->points_vector[end_index1].value << endl;
				out << "A3:" << end_index2 << " " << this->points_vector[end_index2].value << endl;
				temp_mesh_vector.push_back(start_index1);
				temp_mesh_vector.push_back(end_index1);
				temp_mesh_vector.push_back(end_index2);
				//insertNewMesh(start_index1, end_index1, end_index2);
				if (this->points_vector[end_index1].value != this->points_vector[end_index2].value)
				{
					int another_index1 = end_index1;
					int another_index2 = left_edge_vector.at(ei1 + 1);
					int another_index3 = end_index2;
					//push back a new mesh
					out << "======ONE EXTEND======" << endl;
					out << "A1:" << another_index1 << " " << this->points_vector[another_index1].value << endl;
					out << "A2:" << another_index2 << " " << this->points_vector[another_index2].value << endl;
					out << "A3:" << another_index3 << " " << this->points_vector[another_index3].value << endl;
					temp_mesh_vector.push_back(another_index1);
					temp_mesh_vector.push_back(another_index2);
					temp_mesh_vector.push_back(another_index3);
					//insertNewMesh(another_index1, another_index2, another_index3);
					si1 = ei1;
					ei1 = ei1 + 1;
				}
			}
			//the last new mesh
			else if (end_index1 == end_index2)
			{
				//push back a new mesh
				out << "======ONE EXTEND======" << endl;
				out << "A1:" << start_index1 << " " << this->points_vector[start_index1].value << endl;
				out << "A2:" << start_index2 << " " << this->points_vector[start_index2].value << endl;
				out << "A3:" << end_index1 << " " << this->points_vector[end_index1].value << endl;
				temp_mesh_vector.push_back(start_index1);
				temp_mesh_vector.push_back(start_index2);
				temp_mesh_vector.push_back(end_index1);
				//insertNewMesh(start_index1, start_index2, end_index1);
				end_flag = 1;
				//end the remesh
			}
			else
			{
				//push back a new mesh
				out << "======ONE EXTEND======" << endl;
				out << "A1:" << start_index1 << " " << this->points_vector[start_index1].value << endl;
				out << "A2:" << end_index1 << " " << this->points_vector[end_index1].value << endl;
				out << "A3:" << end_index2 << " " << this->points_vector[end_index2].value << endl;
				temp_mesh_vector.push_back(start_index1);
				temp_mesh_vector.push_back(end_index1);
				temp_mesh_vector.push_back(end_index2);
				//insertNewMesh(start_index1, end_index1, end_index2);
				//push back a new mesh
				out << "======ONE EXTEND======" << endl;
				out << "A1:" << start_index1 << " " << this->points_vector[start_index1].value << endl;
				out << "A2:" << start_index2 << " " << this->points_vector[start_index2].value << endl;
				out << "A3:" << end_index2 << " " << this->points_vector[end_index2].value << endl;
				temp_mesh_vector.push_back(start_index1);
				temp_mesh_vector.push_back(start_index2);
				temp_mesh_vector.push_back(end_index2);
				//insertNewMesh(start_index1, start_index2, end_index2);
				if (ei2 == right_edge_vector.size() - 1)
				{
					end_flag = 1;
				}
				else if (this->points_vector[end_index1].value != this->points_vector[end_index2].value)
				{
					int another_index1 = end_index1;
					int another_index2 = left_edge_vector.at(ei1 + 1);
					int another_index3 = end_index2;
					//push back a new mesh
					out << "======ONE EXTEND======" << endl;
					out << "A1:" << another_index1 << " " << this->points_vector[another_index1].value << endl;
					out << "A2:" << another_index2 << " " << this->points_vector[another_index2].value << endl;
					out << "A3:" << another_index3 << " " << this->points_vector[another_index3].value << endl;
					temp_mesh_vector.push_back(another_index1);
					temp_mesh_vector.push_back(another_index2);
					temp_mesh_vector.push_back(another_index3);
					//insertNewMesh(another_index1, another_index2, another_index3);
					si1 = ei1;
					ei1 = ei1 + 1;
				}
			}
			//push the line point
			si1++;
			ei1++;
			si2++;
			ei2++;
		}
		reorganizeMesh(temp_mesh_vector);
	}
	out << "====EDGE VECTOR AFTRE FIX====\n";
	checkEdge();
	checkMesh();
	checkMeshRange();
	out.close();
}

int Element::isMeshAllVisited(vector<int> mesh_vector, int * visit_flag)
{
	for (int i = 0; i < mesh_vector.size(); i++)
	{
		if (visit_flag[i] == -1)
		{
			return i;
		}
	}
	return -1;
}

void Element::buildFinestMRG()
{
	ofstream out3("D:\\out3.txt");

	computeArea();
	vector<double> whole_area;
	vector<double> whole_len;
	mrg = new MRG();
	mrg->initMRG(this->K_range);
	//init visit flag
	int* visit_flag = new int[this->mesh_number];
	vector<int> temp_mesh_queue;
	for (int i = 0; i < this->mesh_number; i++)
	{
		visit_flag[i] = -1;
		temp_mesh_queue.push_back(i);
	}
	//start to create node
	int unvisited_index = isMeshAllVisited(temp_mesh_queue, visit_flag);
	while (unvisited_index != -1)
	{
		out3 << "====ONE NODE====" << endl;
		int node_first_index = temp_mesh_queue.at(unvisited_index);
		int selected_range = this->meshs_vector[node_first_index].range;
		out3 << "range: " << selected_range << endl;

		vector<int> candidate_mesh;
		vector<int> one_node;

		visit_flag[node_first_index] = 1;
		candidate_mesh.push_back(node_first_index);

		while (candidate_mesh.size() != 0)
		{
			int selected_index = candidate_mesh.at(0);
			one_node.push_back(selected_index);
			std::vector<int>::iterator it = candidate_mesh.begin();
			candidate_mesh.erase(it);

			out3 << selected_index << " " << this->meshs_vector[selected_index].range << endl;

			for (int i = 0; i < 3; i++)
			{
				int adjance_index = this->meshs_vector[selected_index].adjance_meshs[i];
				if (adjance_index != -1)
				{
					if (visit_flag[adjance_index] == -1)
					{
						if (this->meshs_vector[adjance_index].range == selected_range)
						{
							visit_flag[adjance_index] = 1;
							candidate_mesh.push_back(adjance_index);
						}
					}
				}
			}
		}
		//get edge queue
		vector<int> edge_queue;
		for (int i = 0; i < one_node.size(); i++)
		{
			for (int j = 0; j < 3; j++)
			{
				int edge_index = this->meshs_vector[one_node.at(i)].boundrys_vector[j];
				if (edge_index != -1)
				{
					if (find(edge_queue.begin(), edge_queue.end(), edge_index) != edge_queue.end()) 
					{
						//do nothing
					}
					else 
					{
						edge_queue.push_back(edge_index);
					}
				}
			}
		}
		int ** MeshonEdge_ForEachNode = new int*[edge_queue.size()];
		for (int i = 0; i < edge_queue.size(); i++)
		{
			MeshonEdge_ForEachNode[i] = new int[2];
			MeshonEdge_ForEachNode[i][0] = this->MeshonEdge[edge_queue.at(i)][0];
			MeshonEdge_ForEachNode[i][1] = this->MeshonEdge[edge_queue.at(i)][1];
		}
		double node_area = 0;
		double max_len = MAX_VALUE;
		double min_len = MIN_VALUE;
		for (int i = 0; i < one_node.size(); i++)
		{
			int selected_mesh_index = one_node.at(i);
			node_area += this->meshs_vector[selected_mesh_index].area;
			for (int j = 0; j < 3; j++)
			{
				int selected_point_index = this->meshs_vector[selected_mesh_index].vertexs[j];
				double selected_value = this->points_vector[selected_point_index].value;
				if (selected_value >= max_len)
				{
					max_len = selected_value;
				}
				if (selected_value <= min_len)
				{
					min_len = selected_value;
				}
			}
		}
		out3 << "node_area: " << node_area << " max - min: " << max_len - min_len << endl;
		whole_area.push_back(node_area);
		whole_len.push_back(max_len - min_len);

		mrg->createMnode(0,0, selected_range, one_node, MeshonEdge_ForEachNode, edge_queue.size());
		unvisited_index = isMeshAllVisited(temp_mesh_queue, visit_flag);
	}

	mrg->createMRG(this->mesh_number);
	mrg->calculateP12(whole_area, whole_len);
	mrg->printGraph_Number();

	out3 << "end partition!"<< endl;
	out3.close();
}

void Element::rebuildAdjanceMesh()
{
	ofstream out2("D:\\out2.txt");
	MeshonEdge = new int*[this->boundary_number];
	int* adjance_number = new int[this->mesh_number];
	for (int i = 0; i < this->mesh_number; i++)
	{
		adjance_number[i] = 0;

	}
	for (int i = 0; i < this->boundary_number; i++)
	{
		MeshonEdge[i] = new int[2];
		MeshonEdge[i][0] = -1;
		MeshonEdge[i][1] = -1;
	}
	for (int i = 0; i < this->mesh_number; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (MeshonEdge[meshs_vector[i].boundrys_vector[j]][0] == -1)
			{
				MeshonEdge[meshs_vector[i].boundrys_vector[j]][0] = i;
			}
			else MeshonEdge[meshs_vector[i].boundrys_vector[j]][1] = i;
		}
	}
	out2 << "====MESH ON EDGE====" << endl;
	for (int i = 0; i < this->boundary_number; i++)
	{
		out2 << i << " "<< MeshonEdge[i][0] << " " << MeshonEdge[i][1] << endl;
	}

	for (int i = 0; i < this->boundary_number; i++)
	{
		int m1 = MeshonEdge[i][0];
		int m2 = MeshonEdge[i][1];
		out2 << "====ONE EDGE====" << endl;
		out2 << i << " " << m1 << " " << m2 << endl;
		if (m2 != -1)
		{
			//for m1
			meshs_vector[MeshonEdge[i][0]].adjance_meshs[adjance_number[MeshonEdge[i][0]]] = MeshonEdge[i][1];
			adjance_number[MeshonEdge[i][0]]++;
			//for m2
			meshs_vector[MeshonEdge[i][1]].adjance_meshs[adjance_number[MeshonEdge[i][1]]] = MeshonEdge[i][0];
			adjance_number[MeshonEdge[i][1]]++;
		}
		else
		{
			//for m1
			meshs_vector[MeshonEdge[i][0]].adjance_meshs[adjance_number[MeshonEdge[i][0]]] = MeshonEdge[i][1];
			adjance_number[MeshonEdge[i][0]]++;
		}
	}

	out2 << "====ADJANCE MESH====" << endl;
	for (int i = 0; i < this->mesh_number; i++)
	{
		out2 << i << " ";
		for (int j = 0; j < 3; j++)
		{
			out2 << meshs_vector[i].adjance_meshs[j] << " ";
		}
		out2 << endl;
	}
	out2.close();
}

void Element::initK_rangeImage()
{
	double p1 = this->points_vector[this->min_value_index].value;
	double p1_x = this->points_vector[this->min_value_index].x;
	double p1_y = this->points_vector[this->min_value_index].y;
}

void Element::drawElement()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glPointSize(5.0f);
	glBegin(GL_POINTS);
	for (int i = 0; i < this->point_number; i++)
	{
		GLdouble x = this->points_vector[i].x;
		GLdouble y = this->points_vector[i].y;
		glVertex2d(x, y);
	}
	glEnd();

	for (int i = 0; i < this->mesh_number; i++)
	{
		int index1 = this->meshs_vector[i].vertexs[0];
		int index2 = this->meshs_vector[i].vertexs[1];
		int index3 = this->meshs_vector[i].vertexs[2];

		GLdouble x1 = this->points_vector[index1].x;
		GLdouble y1 = this->points_vector[index1].y;
		GLdouble x2 = this->points_vector[index2].x;
		GLdouble y2 = this->points_vector[index2].y;
		GLdouble x3 = this->points_vector[index3].x;
		GLdouble y3 = this->points_vector[index3].y;
		glBegin(GL_LINE_LOOP);
		glVertex2d(x1, y1);
		glVertex2d(x2, y2);
		glVertex2d(x3, y3);
		glEnd();
	}
}

void Element::checkVertexRange()
{
	out << "===VERTEX_RANGE_VALUE===" << endl;
	for (int i = 0; i < this->point_number; i++)
	{
		out << i << " " << this->points_vector[i].x << " " << this->points_vector[i].y << " " << this->points_vector[i].value << " " << this->points_vector[i].range << endl;
	}
}

void Element::checkMeshRange()
{
	for (int i = 0; i < this->mesh_number; i++)
	{
		int range_number[3] = { -1,-1,-1 };
		double v1_h = this->points_vector[this->meshs_vector[i].vertexs[0]].value;
		double v2_h = this->points_vector[this->meshs_vector[i].vertexs[1]].value;
		double v3_h = this->points_vector[this->meshs_vector[i].vertexs[2]].value;
		for (int j = 0; j < this->K_range; j++)
		{
			double l_h = this->height_scanning_line[j];
			double h_h = this->height_scanning_line[j + 1];
			if (v1_h >= l_h && v1_h <= h_h)
			{
				if (v2_h >= l_h && v2_h <= h_h)
				{
					if (v3_h >= l_h && v3_h <= h_h)
					{
						this->meshs_vector[i].range = j;
					}
				}
			}
		}
	}

	out << "===MESH_RANGE_VALUE===" << endl;
	for (int i = 0; i < this->mesh_number; i++)
	{
		out << i << " " << this->meshs_vector[i].range << endl;
	}
}

void Element::checkHeight()
{
	out << "===RANGE_HEIGHT===" << endl;
	for (int i = 0; i <= this->K_range; i++)
	{
		out << i << " " << this->height_scanning_line[i] << endl;
	}
}

void Element::checkEdge()
{
	out << "===EDGE VECTOR===" << endl;
	for (int i = 0; i < this->boundary_number; i++)
	{
		out << i << " " << this->boundarys_vector[i].end_point1_index << " " << this->boundarys_vector[i].end_point2_index << " " << this->boundarys_vector[i].value << endl;
	}
}

void Element::checkMesh()
{
	out << "===MESH INFORMATION===" << endl;
	for (int i = 0; i < this->mesh_number; i++)
	{
		out << i << " " << this->meshs_vector[i].vertexs[0] << " " << this->meshs_vector[i].vertexs[1] << " " << this->meshs_vector[i].vertexs[2] << " " << " " << this->meshs_vector[i].boundrys_vector[0] << " " << this->meshs_vector[i].boundrys_vector[1] << " " << this->meshs_vector[i].boundrys_vector[2] << " " << this->meshs_vector[i].range << endl;
	}
}

void Element::reorganizeEdge(int e1, int e2, vector<int> point_on_line)
{
	//get the old edge index
	int edge_index = -1;
	edge_index = isLineinEdgeVector(e1, e2);
	boundry* new_boundarys_vector = new boundry[this->boundary_number + point_on_line.size() - 1 - 1];
	for (int i = 0; i < this->boundary_number; i++)
	{
		new_boundarys_vector[i] = this->boundarys_vector[i];

	}
	//let a new edge to replace the old edge(in order to not change the other edge's side)
	out << "====A NEW EDGE INSERT====\n";
	out << edge_index << " " << point_on_line.at(0) << " " << point_on_line.at(1) << endl;
	new_boundarys_vector[edge_index].end_point1_index = point_on_line.at(0);
	new_boundarys_vector[edge_index].end_point2_index = point_on_line.at(1);
	new_boundarys_vector[edge_index].value = this->points_vector[point_on_line.at(0)].value * OMIGA * 1.0 + this->points_vector[point_on_line.at(1)].value * (1 - OMIGA) * 1.0;
	//other new edge follow the vector end
	int new_index = this->boundary_number;
	for (int i = 1; i < point_on_line.size() - 1; i++)
	{
		out << "====A NEW EDGE INSERT====\n";
		out << this->boundary_number << " " << point_on_line.at(i) << " " << point_on_line.at(i + 1) << endl;
		new_boundarys_vector[new_index].end_point1_index = point_on_line.at(i);
		new_boundarys_vector[new_index].end_point2_index = point_on_line.at(i + 1);
		new_boundarys_vector[new_index].value = this->points_vector[point_on_line.at(i)].value * OMIGA * 1.0 + this->points_vector[point_on_line.at(i + 1)].value * (1 - OMIGA) * 1.0;
		new_index++;
	}
	delete[] this->boundarys_vector;
	this->boundarys_vector = new_boundarys_vector;
	this->boundary_number = this->boundary_number + point_on_line.size() - 1 - 1;
}

void Element::reorganizeMesh(vector<int> temp_mesh_vector)
{
	int new_mesh_number = (temp_mesh_vector.size() - 1) / 3;
	Mesh* new_mesh_vector = new Mesh[this->mesh_number + new_mesh_number - 1];
	for (int i = 0; i < this->mesh_number; i++)
	{
		new_mesh_vector[i] = this->meshs_vector[i];
	}
	//one new mesh replace the old mesh(in order not to influence other meshes)
	int old_index = temp_mesh_vector.front();
	//vertex
	new_mesh_vector[old_index].vertexs = new int[3];
	new_mesh_vector[old_index].vertexs[0] = temp_mesh_vector.at(1);
	new_mesh_vector[old_index].vertexs[1] = temp_mesh_vector.at(2);
	new_mesh_vector[old_index].vertexs[2] = temp_mesh_vector.at(3);
	//neighbor
	new_mesh_vector[old_index].adjance_meshs = new int[3];
	for (int j = 0; j < 3; j++)
	{
		new_mesh_vector[old_index].adjance_meshs[j] = -1;
	}
	//edge
	new_mesh_vector[old_index].boundrys_vector = new int[3];
	int e1 = 0;
	int e2 = 0;
	int edge_flag = -1;
	//the first edge
	e1 = new_mesh_vector[old_index].vertexs[0];
	e2 = new_mesh_vector[old_index].vertexs[1];
	edge_flag = isLineinEdgeVector(e1, e2);
	if (edge_flag == -1)
	{
		insertNewEdge(e1, e2);
		new_mesh_vector[old_index].boundrys_vector[0] = this->boundary_number - 1;
	}
	else new_mesh_vector[old_index].boundrys_vector[0] = edge_flag;
	//the second edge
	e1 = new_mesh_vector[old_index].vertexs[1];
	e2 = new_mesh_vector[old_index].vertexs[2];
	edge_flag = isLineinEdgeVector(e1, e2);
	if (edge_flag == -1)
	{
		insertNewEdge(e1, e2);
		new_mesh_vector[old_index].boundrys_vector[1] = this->boundary_number - 1;
	}
	else new_mesh_vector[old_index].boundrys_vector[1] = edge_flag;
	//the third edge
	e1 = new_mesh_vector[old_index].vertexs[0];
	e2 = new_mesh_vector[old_index].vertexs[2];
	edge_flag = isLineinEdgeVector(e1, e2);
	if (edge_flag == -1)
	{
		insertNewEdge(e1, e2);
		new_mesh_vector[old_index].boundrys_vector[2] = this->boundary_number - 1;
	}
	else new_mesh_vector[old_index].boundrys_vector[2] = edge_flag;
	//area
	new_mesh_vector[old_index].area = Heron(new_mesh_vector[old_index].vertexs[0], new_mesh_vector[old_index].vertexs[1],new_mesh_vector[old_index].vertexs[2]);
	new_mesh_vector[old_index].range = -1;

	//other new mesh follow the vector
	int new_index = this->mesh_number;
	for (int i = 4; i < temp_mesh_vector.size(); i = i + 3)
	{
		//vertex
		new_mesh_vector[new_index].vertexs = new int[3];
		new_mesh_vector[new_index].vertexs[0] = temp_mesh_vector.at(i);
		new_mesh_vector[new_index].vertexs[1] = temp_mesh_vector.at(i + 1);
		new_mesh_vector[new_index].vertexs[2] = temp_mesh_vector.at(i + 2);
		//neighbor
		new_mesh_vector[new_index].adjance_meshs = new int[3];
		for (int j = 0; j < 3; j++)
		{
			new_mesh_vector[new_index].adjance_meshs[j] = -1;
		}
		//edge
		new_mesh_vector[new_index].boundrys_vector = new int[3];
		int e1 = 0;
		int e2 = 0;
		int edge_flag = -1;
		//the first edge
		e1 = new_mesh_vector[new_index].vertexs[0];
		e2 = new_mesh_vector[new_index].vertexs[1];
		edge_flag = isLineinEdgeVector(e1, e2);
		if (edge_flag == -1)
		{
			insertNewEdge(e1, e2);
			new_mesh_vector[new_index].boundrys_vector[0] = this->boundary_number - 1;
		}
		else new_mesh_vector[new_index].boundrys_vector[0] = edge_flag;
		//the second edge
		e1 = new_mesh_vector[new_index].vertexs[1];
		e2 = new_mesh_vector[new_index].vertexs[2];
		edge_flag = isLineinEdgeVector(e1, e2);
		if (edge_flag == -1)
		{
			insertNewEdge(e1, e2);
			new_mesh_vector[new_index].boundrys_vector[1] = this->boundary_number - 1;
		}
		else new_mesh_vector[new_index].boundrys_vector[1] = edge_flag;
		//the third edge
		e1 = new_mesh_vector[new_index].vertexs[0];
		e2 = new_mesh_vector[new_index].vertexs[2];
		edge_flag = isLineinEdgeVector(e1, e2);
		if (edge_flag == -1)
		{
			insertNewEdge(e1, e2);
			new_mesh_vector[new_index].boundrys_vector[2] = this->boundary_number - 1;
		}
		else new_mesh_vector[new_index].boundrys_vector[2] = edge_flag;
		//area
		new_mesh_vector[new_index].area = Heron(new_mesh_vector[new_index].vertexs[0], new_mesh_vector[new_index].vertexs[1], new_mesh_vector[new_index].vertexs[2]);
		new_mesh_vector[new_index].range = -1;
		new_index++;
	}

	delete[] this->meshs_vector;
	this->meshs_vector = new_mesh_vector;
	this->mesh_number = this->mesh_number + new_mesh_number - 1;
}

void Element::insertNewVertex(double x, double y, double value, int range)
{
	out << "insert a new vertex: " << x << " " << y << " " << value << " " << range << endl;
	vertex* new_points_vector = NULL;
	new_points_vector = new vertex[this->point_number + 1];
	for (int i = 0; i < this->point_number; i++)
	{
		new_points_vector[i] = this->points_vector[i];
	}
	new_points_vector[this->point_number].range = range;
	new_points_vector[this->point_number].value = value;
	new_points_vector[this->point_number].x = x;
	new_points_vector[this->point_number].y = y;
	this->point_number++;
	delete[] this->points_vector;
	this->points_vector = new_points_vector;
}

void Element::insertNewMesh(int p1, int p2, int p3)
{
	Mesh* new_mesh_vector = new Mesh[this->mesh_number + 1];
	for (int i = 0; i < this->mesh_number; i++)
	{
		new_mesh_vector[i] = this->meshs_vector[i];
	}

	//vertex
	new_mesh_vector[this->mesh_number].vertexs = new int[3];
	new_mesh_vector[this->mesh_number].vertexs[0] = p1;
	new_mesh_vector[this->mesh_number].vertexs[1] = p2;
	new_mesh_vector[this->mesh_number].vertexs[2] = p3;
	//neighbor
	new_mesh_vector[this->mesh_number].adjance_meshs = new int[3];
	for (int j = 0; j < 3; j++)
	{
		new_mesh_vector[this->mesh_number].adjance_meshs[j] = -1;
	}
	//edge
	new_mesh_vector[this->mesh_number].boundrys_vector = new int[3];
	int e1 = 0;
	int e2 = 0;
	int edge_index = -1;
	//the first edge
	e1 = p1;
	e2 = p2;
	edge_index = isLineinEdgeVector(e1, e2);
	if (edge_index == -1)
	{
		insertNewEdge(e1, e2);
		new_mesh_vector[this->mesh_number].boundrys_vector[0] = this->boundary_number;
	}
	else new_mesh_vector[this->mesh_number].boundrys_vector[0] = edge_index;
	//the second edge
	e1 = p2;
	e2 = p3;
	edge_index = isLineinEdgeVector(e1, e2);
	if (edge_index == -1)
	{
		insertNewEdge(e1, e2);
		new_mesh_vector[this->mesh_number].boundrys_vector[1] = this->boundary_number;
	}
	else new_mesh_vector[this->mesh_number].boundrys_vector[1] = edge_index;
	//the third edge
	e1 = p1;
	e2 = p3;
	edge_index = isLineinEdgeVector(e1, e2);
	if (edge_index == -1)
	{
		insertNewEdge(e1, e2);
		new_mesh_vector[this->mesh_number].boundrys_vector[2] = this->boundary_number;
	}
	else new_mesh_vector[this->mesh_number].boundrys_vector[2] = edge_index;
	//area
	new_mesh_vector[this->mesh_number].area = Heron(p1,p2,p2);
	new_mesh_vector[this->mesh_number].range = -1;
	delete[] this->meshs_vector;
	this->meshs_vector = new_mesh_vector;
	this->mesh_number++;
}

void Element::insertNewEdge(int e1, int e2)
{

	boundry* new_boundarys_vector = new boundry[this->boundary_number + 1];
	for (int i = 0; i < this->boundary_number; i++)
	{
		new_boundarys_vector[i] = this->boundarys_vector[i];
	}
	new_boundarys_vector[this->boundary_number].end_point1_index = e1;
	new_boundarys_vector[this->boundary_number].end_point2_index = e2;
	new_boundarys_vector[this->boundary_number].value = this->points_vector[e1].value * OMIGA * 1.0 + this->points_vector[e2].value * (1 - OMIGA) * 1.0;
	delete[] this->boundarys_vector;
	this->boundarys_vector = new_boundarys_vector;
	this->boundary_number++;

}

void Element::checkVertex(int index)
{
	vertex e = this->points_vector[index];
	out << "a vertex: " << index << " " << e.x << " " << e.y << " " << e.value << " " << e.range << endl;
}

int Element::isLineinEdgeVector(int e1, int e2)
{
	for (int i = 0; i < this->boundary_number; i++)
	{
		int se1 = this->boundarys_vector[i].end_point1_index;
		int se2 = this->boundarys_vector[i].end_point2_index;
		if (se1 == e1 && se2 == e2 || se1 == e2 && se2 == e1)
		{
			return i;
		}
	}
	return -1;
}

int Element::isLineinFixedEdgevector(vector<vector<int>> f_e_v, int e1, int e2)
{
	for (int i = 0; i < f_e_v.size(); i++)
	{
		int se1 = f_e_v.at(i).at(0);
		int se2 = f_e_v.at(i).at(1);
		if (se1 == e1 && se2 == e2 || se1 == e2 && se2 == e1)
		{
			return i;
		}
	}
	return -1;
}
