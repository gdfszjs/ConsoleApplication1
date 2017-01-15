#include "MRG.h"

MRG::MRG()
{
}

MRG::~MRG()
{
}

void MRG::initMRG(int K_range)
{
	out4.open("D:\\out4.txt");
	this->K_range = K_range;
	this->level_number = log(K_range) / log(2) + 1;
	this->node_number = 0;
	this->node_prop = new int[this->level_number];
	for (int i = 0; i < this->level_number; i++)
	{
		this->node_prop[i] = -1;
	}
	this->Mnode_vector = NULL;

	this->MGraph = new int**[this->level_number];
	this->range_in_each_level = new int[this->level_number];
	for (int i = 0; i < this->level_number; i++)
	{
		int range_number = this->K_range / pow(2, i);
		this->range_in_each_level[i] = range_number;
		this->MGraph[i] = new int*[range_number];
		for (int j = 0; j < range_number; j++)
		{
			this->MGraph[i][j] = NULL;
		}
	}

	this->graph_number = new int*[this->level_number];
	for (int i = 0; i < this->level_number; i++)
	{
		int range_number = this->K_range / pow(2, i);
		this->graph_number[i] = new int[range_number];
		for (int j = 0; j < range_number; j++)
		{
			this->graph_number[i][j] = 0;
		}
	}
}

void MRG::createMnode(int graph, int level, int range, vector<int> Tset, int ** MeshonEdge_ForEachNode, int edge_number)
{
	if (this->node_number == 0)
	{
		this->Mnode_vector = new Mnode[1];

		this->Mnode_vector[this->node_number].index = this->node_number;
		this->Mnode_vector[this->node_number].ngraph = graph;
		this->Mnode_vector[this->node_number].nlevel = level;
		this->Mnode_vector[this->node_number].nrange = range;
		this->Mnode_vector[this->node_number].nnumber = this->graph_number[level][range];
		insertMnodeIndex(this->node_number, graph, level, range);
		this->Mnode_vector[this->node_number].Tset = new int[Tset.size()];
		this->Mnode_vector[this->node_number].Tset_number = Tset.size();
		for (int i = 0; i < Mnode_vector[this->node_number].Tset_number; i++)
		{
			this->Mnode_vector[this->node_number].Tset[i] = Tset.at(i);
		}
		this->Mnode_vector[this->node_number].adjance_index = NULL;
		this->Mnode_vector[this->node_number].adjance_number = 0;
		this->Mnode_vector[this->node_number].child_index = NULL;
		this->Mnode_vector[this->node_number].child_number = 0;
		this->Mnode_vector[this->node_number].parent_index = -1;
		this->Mnode_vector[this->node_number].p1 = 0;
		this->Mnode_vector[this->node_number].p2 = 0;

		this->Mnode_vector[this->node_number].edge_number = edge_number;
		this->Mnode_vector[this->node_number].line_table = new int*[edge_number];
		for (int i = 0; i < edge_number; i++)
		{
			this->Mnode_vector[this->node_number].line_table[i] = new int[2];
			this->Mnode_vector[this->node_number].line_table[i][0] = MeshonEdge_ForEachNode[i][0];
			this->Mnode_vector[this->node_number].line_table[i][1] = MeshonEdge_ForEachNode[i][1];
		}

		this->node_number++;
		this->graph_number[level][range]++;
	}
	else
	{
		int old_node_number = this->node_number;
		Mnode* new_Mnode = new Mnode[old_node_number + 1];
		for (int i = 0; i < old_node_number; i++)
		{
			new_Mnode[i] = Mnode_vector[i];
		}

		new_Mnode[old_node_number].index = old_node_number;
		new_Mnode[old_node_number].ngraph = graph;
		new_Mnode[old_node_number].nlevel = level;
		new_Mnode[old_node_number].nrange = range;
		new_Mnode[old_node_number].nnumber = this->graph_number[level][range];
		insertMnodeIndex(old_node_number, graph, level, range);
		new_Mnode[old_node_number].Tset = new int[Tset.size()];
		new_Mnode[old_node_number].Tset_number = Tset.size();
		for (int i = 0; i < new_Mnode[old_node_number].Tset_number; i++)
		{
			new_Mnode[old_node_number].Tset[i] = Tset.at(i);
		}
		new_Mnode[old_node_number].adjance_index = NULL;
		new_Mnode[old_node_number].adjance_number = 0;
		new_Mnode[old_node_number].child_index = NULL;
		new_Mnode[old_node_number].child_number = 0;
		new_Mnode[old_node_number].parent_index = -1;
		new_Mnode[old_node_number].p1 = 0;
		new_Mnode[old_node_number].p2 = 0;

		new_Mnode[old_node_number].edge_number = edge_number;
		new_Mnode[old_node_number].line_table = new int*[edge_number];
		for (int i = 0; i < edge_number; i++)
		{
			new_Mnode[old_node_number].line_table[i] = new int[2];
			new_Mnode[old_node_number].line_table[i][0] = MeshonEdge_ForEachNode[i][0];
			new_Mnode[old_node_number].line_table[i][1] = MeshonEdge_ForEachNode[i][1];
		}

		delete[] this->Mnode_vector;
		this->Mnode_vector = new_Mnode;

		this->node_number++;
		this->graph_number[level][range]++;
	}
}

void MRG::insertMnodeIndex(int index, int graph, int level, int range)
{
	if (this->graph_number[level][range] == 0)
	{
		this->MGraph[level][range] = new int[1];
		this->MGraph[level][range][0] = index;
	}
	else
	{
		int * new_MGraph = new int[this->graph_number[level][range] + 1];
		for (int i = 0; i < this->graph_number[level][range]; i++)
		{
			new_MGraph[i] = this->MGraph[level][range][i];
		}
		new_MGraph[this->graph_number[level][range]] = index;
		delete[] this->MGraph[level][range];
		this->MGraph[level][range] = new_MGraph;
	}
}

void MRG::printLine_Queue()
{
	for (int i = 0; i < this->level_number; i++)
	{
		for (int j = 0; j < this->range_in_each_level[i]; j++)
		{
			out4 << "====LEVEL:" << i << "  RANGE:" << j << "====" << endl;
			out4 << "NODE NUMBER:" << this->graph_number[i][j] << endl;
			for (int k = 0; k < this->graph_number[i][j]; k++)
			{
				out4 << "NODE INDEX:" << this->MGraph[i][j][k] << endl;
				for (int l = 0; l < this->Mnode_vector[this->MGraph[i][j][k]].edge_number; l++)
				{
					out4 << "LINE QUEUE:" << this->Mnode_vector[this->MGraph[i][j][k]].line_table[l][0] << " " << this->Mnode_vector[this->MGraph[i][j][k]].line_table[l][1] << endl;
				}
			}
		}
	}
}

void MRG::printGraph_Number()
{
	for (int i = 0; i < this->level_number; i++)
	{
		for (int j = 0; j < this->range_in_each_level[i]; j++)
		{
			out4 << "====LEVEL:" << i << "  RANGE:" << j << "====" << endl;
			out4 << "NODE NUMBER:" << this->graph_number[i][j] << endl;
			for (int k = 0; k < this->graph_number[i][j]; k++)
			{
				out4 << "NODE INDEX:" << this->MGraph[i][j][k] << endl;
				for (int l = 0; l < this->Mnode_vector[this->MGraph[i][j][k]].adjance_number; l++)
				{
					out4 << "ADJANCE NODE INDEX:" << this->Mnode_vector[this->MGraph[i][j][k]].adjance_index[l] << endl;
				}
			}
		}
	}
}

void MRG::createMRG(int mesh_number)
{
	this->getAdjanceMnodeinFinestlevel(mesh_number);
	for (int i = 1; i < this->returnLevel_Number(); i++)
	{
		this->getAdjanceMnodeinlevel(i);
	}
}

void MRG::getAdjanceMnodeinlevel(int level)
{
	ofstream out5("D:\\out5.txt");
	int node_number = 0;
	int last_level = level - 1;
	for (int j = 0; j < this->range_in_each_level[last_level]; j++)
	{
		node_number += this->graph_number[last_level][j];
	}
	out5 << "node_number: " << node_number << endl;

	//init visit flag
	int* visit_flag = new int[this->node_number];
	vector<int> temp_node_queue;
	for (int i = 0; i < this->node_number; i++)
	{
		visit_flag[i] = -1;
	}
	
	for (int j = 0; j < this->range_in_each_level[last_level]; j++)
	{
		for (int k = 0; k < this->graph_number[last_level][j]; k++)
		{
			int temp_node_index = this->Mnode_vector[this->MGraph[last_level][j][k]].index;
			out5 << " selected node: " << temp_node_index << endl;
			temp_node_queue.push_back(temp_node_index);
		}
	}
	//start to create node
	int unvisited_index = isNodeAllVisited(temp_node_queue, visit_flag);
	while (unvisited_index != -1)
	{
		out5 << "====ONE NODE====" << endl;
		int node_first_index = temp_node_queue.at(unvisited_index);
		int selected_range = this->Mnode_vector[node_first_index].nrange / 2;
		out5 << "old_range: " << this->Mnode_vector[node_first_index].nrange << " new_range: " << selected_range << endl;

		vector<int> candidate_node;
		vector<int> one_node;

		visit_flag[node_first_index] = 1;
		candidate_node.push_back(node_first_index);

		while (candidate_node.size() != 0)
		{
			int selected_index = candidate_node.at(0);
			one_node.push_back(selected_index);
			std::vector<int>::iterator it = candidate_node.begin();
			candidate_node.erase(it);

			out5 << selected_index << " " << this->Mnode_vector[selected_index].nrange << " " << this->Mnode_vector[selected_index].nrange / 2 << endl;

			for (int i = 0; i < this->Mnode_vector[selected_index].adjance_number; i++)
			{
				int adjance_index = this->Mnode_vector[selected_index].adjance_index[i];
				if (adjance_index != -1)
				{
					if (visit_flag[adjance_index] == -1)
					{
						int geted_range = this->Mnode_vector[adjance_index].nrange / 2;
						if (geted_range == selected_range)
						{
							visit_flag[adjance_index] = 1;
							candidate_node.push_back(adjance_index);
						}
					}
				}
			}
		}
		vector<int> edge_queue;
		for (int i = 0; i < one_node.size(); i++)
		{
			for (int j = 0; j < this->Mnode_vector[one_node.at(i)].adjance_number; j++)
			{
				int ad_index = this->Mnode_vector[one_node.at(i)].adjance_index[j];
				int se_index = this->Mnode_vector[one_node.at(i)].index;
				edge_queue.push_back(se_index);
				edge_queue.push_back(ad_index);
			}
		}
		int edge_number = edge_queue.size() / 2;
		int ** MeshonEdge_ForEachNode = new int*[edge_number];
		for (int i = 0; i < edge_number; i++)
		{
			MeshonEdge_ForEachNode[i] = new int[2];
			MeshonEdge_ForEachNode[i][0] = edge_queue.at(2 * i);
			MeshonEdge_ForEachNode[i][1] = edge_queue.at(2 * i + 1);
		}

		createMnode(0, level, selected_range, one_node, MeshonEdge_ForEachNode, edge_number);
		unvisited_index = isNodeAllVisited(temp_node_queue, visit_flag);
	}
	printGraph_Number();
	//find adjance node

	int* node_prop = new int[this->node_number];
	int i = 0;
	for (int j = 0; j < this->range_in_each_level[level]; j++)
	{
		for (int k = 0; k < this->graph_number[level][j]; k++)
		{
			int one_Tset_number = this->Mnode_vector[this->MGraph[level][j][k]].Tset_number;
			for (int l = 0; l < one_Tset_number; l++)
			{
				int node_index = this->Mnode_vector[this->MGraph[level][j][k]].Tset[l];
				node_prop[node_index] = this->Mnode_vector[this->MGraph[level][j][k]].index;
			}
		}
	}
	out5 << "====MESH PROP====" << endl;
	for (int j = 0; j < node_number; j++)
	{
		out5 << j << " " << node_prop[j] << endl;
	}
	out5 << "====START FIND AD INDEX====" << endl;
	for (int j = 0; j < this->range_in_each_level[level]; j++)
	{
		for (int k = 0; k < this->graph_number[level][j]; k++)
		{
			vector<int> ad_index;

			int edge_number = this->Mnode_vector[this->MGraph[level][j][k]].edge_number;
			int selected_index = this->Mnode_vector[this->MGraph[level][j][k]].index;
			out5 << "====FIND ADJANCE INDEX====" << endl;
			out5 << "selected_index: " << selected_index << endl;
			for (int l = 0; l < edge_number; l++)
			{
				int i1 = this->Mnode_vector[this->MGraph[level][j][k]].line_table[l][0];
				int i2 = this->Mnode_vector[this->MGraph[level][j][k]].line_table[l][1];
				int i1p = node_prop[i1];
				int i2p = node_prop[i2];
				out5 << "ONE EDGE: " << i1 << " " << i2 << " " << i1p << " " << i2p << endl;
				if (node_prop[i1] == selected_index)
				{
					if (i2 != -1 && node_prop[i2] != selected_index)
					{
						if (find(ad_index.begin(), ad_index.end(), node_prop[i2]) != ad_index.end())
						{
							//do nothing
						}
						else
						{
							ad_index.push_back(node_prop[i2]);
						}
					}
				}
				else
				{
					if (i1 != -1 && node_prop[i1] != selected_index)
					{
						if (find(ad_index.begin(), ad_index.end(), node_prop[i1]) != ad_index.end())
						{
							//do nothing
						}
						else
						{
							ad_index.push_back(node_prop[i1]);
						}
					}
				}
			}
			this->Mnode_vector[this->MGraph[level][j][k]].adjance_number = ad_index.size();
			this->Mnode_vector[this->MGraph[level][j][k]].adjance_index = new int[ad_index.size()];
			for (int l = 0; l < ad_index.size(); l++)
			{
				out5 << "SELECTED INDEX:" << ad_index.at(l) << endl;
				this->Mnode_vector[this->MGraph[level][j][k]].adjance_index[l] = ad_index.at(l);
			}

		}
	}
	out5.close();
}

int MRG::isNodeAllVisited(vector<int> node_vector, int * visit_flag)
{
	for (int i = 0; i < node_vector.size(); i++)
	{
		if (visit_flag[node_vector.at(i)] == -1)
		{
			return i;
		}
	}
	return -1;
}

void MRG::getAdjanceMnodeinFinestlevel(int mesh_number)
{
	int* mesh_prop = new int[mesh_number];
	int i = 0;
	for (int j = 0; j < this->range_in_each_level[i]; j++)
	{
		for (int k = 0; k < this->graph_number[i][j]; k++)
		{
			int node_number = this->Mnode_vector[this->MGraph[i][j][k]].Tset_number;
			for (int l = 0; l < node_number; l++)
			{
				int mesh_index = this->Mnode_vector[this->MGraph[i][j][k]].Tset[l];
				mesh_prop[mesh_index] = this->Mnode_vector[this->MGraph[i][j][k]].index;
			}
		}
	}
	out4 << "====MESH PROP====" << endl;
	for (int j = 0; j < mesh_number; j++)
	{
		out4 << j << " " << mesh_prop[j] << endl;
	}
	out4 << "====START FIND AD INDEX====" << endl;
	for (int j = 0; j < this->range_in_each_level[i]; j++)
	{
		for (int k = 0; k < this->graph_number[i][j]; k++)
		{
			vector<int> ad_index;
			
			int edge_number = this->Mnode_vector[this->MGraph[i][j][k]].edge_number;
			int selected_index = this->Mnode_vector[this->MGraph[i][j][k]].index;
			out4 << "====FIND ADJANCE INDEX====" << endl;
			out4 << "selected_index: " << selected_index << endl;
			for (int l = 0; l < edge_number; l++)
			{
				int i1 = this->Mnode_vector[this->MGraph[i][j][k]].line_table[l][0];
				int i2 = this->Mnode_vector[this->MGraph[i][j][k]].line_table[l][1];
				int i1p = mesh_prop[i1];
				int i2p = mesh_prop[i2];
				out4 << "ONE EDGE: " << i1 << " " << i2 << " " << i1p << " " << i2p << endl;
				if (mesh_prop[i1] == selected_index)
				{
					if (i2 != -1 && mesh_prop[i2] != selected_index)
					{
						if (find(ad_index.begin(), ad_index.end(), mesh_prop[i2]) != ad_index.end())
						{
							//do nothing
						}
						else
						{
							ad_index.push_back(mesh_prop[i2]);
						}
					}
				}
				else
				{
					if (i1 != -1 && mesh_prop[i1] != selected_index)
					{
						if (find(ad_index.begin(), ad_index.end(), mesh_prop[i1]) != ad_index.end())
						{
							//do nothing
						}
						else
						{
							ad_index.push_back(mesh_prop[i1]);
						}
					}
				}
			}
			this->Mnode_vector[this->MGraph[i][j][k]].adjance_number = ad_index.size();
			this->Mnode_vector[this->MGraph[i][j][k]].adjance_index = new int[ad_index.size()];
			for (int l = 0; l < ad_index.size(); l++)
			{
				out4 << "SELECTED INDEX:" << ad_index.at(l) << endl;
				this->Mnode_vector[this->MGraph[i][j][k]].adjance_index[l] = ad_index.at(l);
			}
			
		}
	}
}

void MRG::calculateP12(vector<double> whole_area, vector<double> whole_len)
{
	double areaS = 0;
	double lenS = 0;
	for (int i = 0; i < whole_area.size(); i++)
	{
		areaS += whole_area.at(i);
		lenS += whole_len.at(i);
	}

	for (int i = 0; i < whole_area.size(); i++)
	{
		double this_area = 1.0 / this->level_number * (whole_area.at(i) / areaS);
		double this_len = 1.0 / this->level_number * (whole_len.at(i) / lenS);
		this->Mnode_vector[i].p1 = this_area;
		this->Mnode_vector[i].p2 = this_len;
	}
	for (int i = 1; i < this->level_number; i++)
	{
		for (int j = 0; j < this->range_in_each_level[i]; j++)
		{
			for (int k = 0; k < this->graph_number[i][j]; k++)
			{
				double this_area = 0;
				double this_len = 0;
				for (int l = 0; l  < this->Mnode_vector[this->MGraph[i][j][k]].Tset_number; l++)
				{
					this_area += this->Mnode_vector[this->Mnode_vector[this->MGraph[i][j][k]].Tset[l]].p1;
					this_len += this->Mnode_vector[this->Mnode_vector[this->MGraph[i][j][k]].Tset[l]].p2;
				}
				this->Mnode_vector[this->MGraph[i][j][k]].p1 = this_area;
				this->Mnode_vector[this->MGraph[i][j][k]].p2 = this_len;
			}
		}
	}

	for (int i = this->level_number - 1; i > 0; i--)
	{
		for (int j = 0; j < this->range_in_each_level[i]; j++)
		{
			for (int k = 0; k < this->graph_number[i][j]; k++)
			{
				int selected_index = this->MGraph[i][j][k];
				int child_number = this->Mnode_vector[selected_index].Tset_number;
				this->Mnode_vector[selected_index].child_number = child_number;
				this->Mnode_vector[selected_index].child_index = new int[child_number];
				for (int l = 0; l < child_number; l++)
				{
					int selected_child_index = this->Mnode_vector[selected_index].Tset[l];
					this->Mnode_vector[selected_child_index].parent_index = selected_index;
					this->Mnode_vector[selected_index].child_index[l] = selected_child_index;
				}
			}
		}
	}

}

void MRG::resetNGraph(int e)
{
	for (int i = 0; i < this->level_number; i++)
	{
		for (int j = 0; j < this->range_in_each_level[i]; j++)
		{
			for (int k = 0; k < this->graph_number[i][j]; k++)
			{
				this->Mnode_vector[this->MGraph[i][j][k]].ngraph = e;
			}
		}
	}
}

Mnode * MRG::returnNodeByIndex(int index)
{
	return &(this->Mnode_vector[index]);
}

Mnode * MRG::returnTopNode()
{
	return &(this->Mnode_vector[this->MGraph[this->level_number - 1][0][0]]);
}
