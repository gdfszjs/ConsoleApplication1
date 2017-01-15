#include "ElementComparator.h"
#include <math.h>
#define w 0.5
#define MAX_VALUE -10000
#define MIN_VALUE 10000
ElementComparator::ElementComparator()
{
	this->R = NULL;
	this->S = NULL;
	this->result = -1;
}

ElementComparator::ElementComparator(Element *a, Element *b)
{
	out6.open("D:\\out6.txt");

	this->MRG_vector[0] = a->returnMRG();
	this->MRG_vector[1] = b->returnMRG();
	this->R = a->returnMRG();
	this->S = b->returnMRG();

	this->R->resetNGraph(0);
	this->S->resetNGraph(1);

	this->result = 0;
}

ElementComparator::~ElementComparator()
{

}

Mnode * ElementComparator::createFictitiousNode(int ngraph, vector<int> Tset, vector<int> ad_index, double p1, double p2)
{
	Mnode * e = new Mnode;
	e->ngraph = ngraph;
	e->Tset = new int[Tset.size()];
	for (int i = 0; i < Tset.size(); i++)
	{
		e->Tset[i] = Tset[i];
	}
	e->Tset_number = Tset.size();
	e->adjance_index = new int[ad_index.size()];
	for (int i = 0; i < ad_index.size(); i++)
	{
		e->adjance_index[i] = ad_index[i];
	}
	e->adjance_number = ad_index.size();
	e->p1 = p1;
	e->p2 = p2;
	e->fictitious_flag = 1;
	return e;
}

Mnode* ElementComparator::adj(Mnode *m)
{
	int ngraph = m->ngraph;
	vector<int> ad_index;
	vector<int> Tset;
	double p1 = 0;
	double p2 = 0;
	for (int i = 0; i < m->adjance_number; i++)
	{
		int one_ad_index = m->adjance_index[i];
		ad_index.push_back(one_ad_index);
		Mnode * ad_node = MRG_vector[ngraph]->returnNodeByIndex(one_ad_index);
		p1 += ad_node->p1;
		p2 += ad_node->p2;
	}
	return createFictitiousNode(ngraph, Tset, ad_index, p1, p2);
}

double ElementComparator::loss(Mnode * m, Mnode * n)
{
	return 0.5 * (sim(m, m) + sim(n, n)) - sim(m, n);
}

double ElementComparator::mat(Mnode * m, Mnode * n)
{
	return -loss(m, n) - loss(adj(m), adj(n));
}

double ElementComparator::sim(Mnode * m, Mnode * n)
{
	return w * min(m->p1, n->p1) + (1 - w) * min(m->p2, n->p2);
}

double ElementComparator::getCompareResult()
{
	out6 << "========START COMPARE========" << endl;
	Initialization();
	Matching();
	calculate_result();
	out6.close();
	return this->result;
}

void ElementComparator::Initialization()
{
	NLIST.push_back(this->R->returnTopNode());
	NLIST.push_back(this->S->returnTopNode());
}

void ElementComparator::Matching()
{
	Mnode *left_node = NULL;
	Mnode *right_node = NULL;

	while (NLIST.size() != 0)
	{
		//notice that pair has not only two node
		vector<Mnode*> one_pair;
		//notice that inde pair has not only two node index
		vector<int> pair_index;

		double max_sim = MAX_VALUE;
		int max_index = MAX_VALUE;

		int left_index = 0;
		int right_index = 0;

		out6 << "====LEFT NODE====" << endl;
		for (int i = 0; i < NLIST.size(); i++)
		{
			Mnode *selected_node = NLIST.at(i);
			double present_sim = sim(selected_node, selected_node);
			if (present_sim >= max_sim)
			{
				max_sim = present_sim;
				max_index = i;
			}
		}
		if (max_index == -1)
		{
			out6 << "====MAX INDEX -1 SELECTED FALIED!====" << endl;
		}

		left_node = NLIST.at(max_index);
		one_pair.push_back(left_node);
		pair_index.push_back(max_index);
		getNodeInformation(left_node);

		vector<int> right_node_candidate_index_vector;
		out6 << "====CANDIDATE NODE====" << endl;
		for (int i = 0; i < NLIST.size(); i++)
		{

			Mnode *candidate_node = NLIST.at(i);
			
			if (BelongToDifferentPart(left_node, candidate_node))
			{
				if (HasSameRange(left_node, candidate_node))
				{
					if (HasParentPair(left_node, candidate_node))
					{
						if (HasSameMLIST(left_node, candidate_node))
						{
							right_node_candidate_index_vector.push_back(i);

						}
					}
				}
			}
		}

		if (right_node_candidate_index_vector.size() == 1)
		{
			out6 << "=====ONE RIGHT NODE=====" << endl;
			int right_index = right_node_candidate_index_vector.at(0);
			right_node = NLIST.at(right_index);

			getNodeInformation(right_node);

			one_pair.push_back(right_node);
			pair_index.push_back(right_index);
		}
		else if (right_node_candidate_index_vector.size() == 0)
		{
			out6 << "=====NO RIGHT NODE=====" << endl;
		}
		else
		{
			out6 << "=====LOTS RIGHT NODE=====" << endl;

			int max_candidate_index = -1;
			double max_value = MAX_VALUE;
			vector<Mnode*> completed_Node_Vector = completeNodeVector(right_node_candidate_index_vector);
			for (int i = 0; i < completed_Node_Vector.size(); i++)
			{
				Mnode * candidate_node = completed_Node_Vector[i];
				getNodeInformation(candidate_node);
				double new_value = mat(left_node, candidate_node);
				out6 << "THIS NODE VALUE: "<< new_value << endl;
				if (new_value >= max_value)
				{
					max_value = new_value;
					max_candidate_index = i;
				}
			}

			out6 << "=====GET THE LARGEST RIGHT NODE=====" << endl;
			right_node = completed_Node_Vector.at(max_candidate_index);

			getNodeInformation(right_node);

			one_pair.push_back(right_node);
			if (right_node->fictitious_flag == 1)
			{
				for (int i = 0; i < right_node->child_number; i++)
				{
					pair_index.push_back(right_node->child_index[i]);
				}
			}
			else pair_index.push_back(right_node_candidate_index_vector[max_candidate_index]);
		}
		//only there is a node pair,we do insert, boardcast and unpack
		if (one_pair.size() > 1)
		{
			out6 << "====THERE IS A NODE PAIR====" << endl;
			//set the visit flag
			for (int i = 0; i < one_pair.size(); i++)
			{
				Mnode *e = one_pair[i];
				if (e->fictitious_flag == 1)
				{
					for (int j = 0; j < e->child_number; j++)
					{
						NLIST[e->child_index[j]]->visit_flag = 1;
					}
				}
				else e->visit_flag = 1;
			}
			MPAIR.push_back(one_pair);
			BoardCastLabel(MPAIR.size() - 1);
			Unpacking(pair_index, one_pair);
		}
		else
		{
			out6 << "====ONLY LEFT NODE====" << endl;
			//so just remove it
			std::vector<Mnode*>::iterator it;
			//remove the pair node
			for (int i = 0; i < pair_index.size(); i++)
			{
				int selected_node_index = pair_index.at(i);
				it = this->NLIST.begin() + selected_node_index;
				NLIST.erase(it);
				for (int j = i + 1; j < pair_index.size(); j++)
				{
					int repair_index = pair_index.at(j);
					if (repair_index < selected_node_index)
					{
						//do nothing
					}
					else
					{
						pair_index[j] = repair_index - 1;
					}
				}
			}
		}
	}
}
vector<Mnode*> ElementComparator::completeNodeVector(vector<int> old_node_index_vector)
{
	vector<vector<int>> a;
	for (int i = 0; i < old_node_index_vector.size(); i++)
	{
		Mnode * candidate_node = NLIST[old_node_index_vector[i]];
		int present_index = old_node_index_vector[i];
		for (int j = 0; j < candidate_node->adjance_number; j++)
		{
			int ad_index = candidate_node->adjance_index[j];
			int k = 0;
			for (; k < a.size(); k++)
			{
				if (ad_index == a.at(k)[0])
				{
					a.at(k).push_back(present_index);
					break;
				}
			}
			if (k >= a.size())
			{
				vector<int> b;
				b.push_back(ad_index);
				b.push_back(present_index);
				a.push_back(b);
			}
		}
	}
	vector<Mnode*> new_node_vector;
	for (int i = 0; i < old_node_index_vector.size(); i++)
	{
		new_node_vector.push_back(NLIST[old_node_index_vector[i]]);
	}
	for (int i = 0; i < a.size(); i++)
	{
		if (a[i].size() >= 3)
		{
			int ngraph = NLIST[old_node_index_vector[0]]->ngraph;
			//store the node index actully in it's MRG(as the Tset node)
			vector<int> Tset;
			vector<int> ad_index;
			//stor the node index in the NLIST(as the child node)
			vector<int> child_index;
			double p1 = 0;
			double p2 = 0;
			int duplicate_flag = 0;
			//calculate the usual property for the fictitious node
			for (int j = 1; j < a[i].size(); j++)
			{
				Mnode * inner_node = NLIST[a[i][j]];
				Tset.push_back(inner_node->index);
				for (int k = 0; k < inner_node->adjance_number; k++)
				{
					for (int l = 0; l < ad_index.size(); l++)
					{
						if (inner_node->adjance_index[k] == ad_index[l])
						{
							duplicate_flag = 1;
							break;
						}
					}
					if (duplicate_flag == 0)
					{
						ad_index.push_back(inner_node->adjance_index[k]);
					}
					duplicate_flag = 0;
				}
				child_index.push_back(a[i][j]);
				p1 += inner_node->p1;
				p2 += inner_node->p2;
			}
			Mnode * new_c_node = createFictitiousNode(ngraph, Tset, ad_index, p1, p2);
			//additional property for the fictitious node
			new_c_node->child_index = new int[child_index.size()];
			for (int j = 0; j < child_index.size(); j++)
			{
				new_c_node->child_index[j] = child_index[j];
			}
			new_c_node->child_number = child_index.size();
			//push back the fictitious node
			new_node_vector.push_back(new_c_node);
		}
	}
	return new_node_vector;
}

void ElementComparator::Unpacking(vector<int> pair_index, vector<Mnode*> one_pair)
{
	std::vector<Mnode*>::iterator it;
	//remove the pair node
	for (int i = 0; i < pair_index.size(); i++)
	{
		int selected_node_index = pair_index.at(i);
		it = this->NLIST.begin() + selected_node_index;
		NLIST.erase(it);
		for (int j = i + 1; j < pair_index.size(); j++)
		{
			int repair_index = pair_index.at(j);
			if (repair_index < selected_node_index)
			{
				//do nothing
			}
			else
			{
				pair_index[j] = repair_index - 1;
			}
		}
	}
	//then unpack them
	for (int i = 0; i < one_pair.size(); i++)
	{
		Mnode * e = one_pair.at(i);
		if (e->fictitious_flag == 1)
		{
			for (int j = 0; j < e->Tset_number; j++)
			{
				int Tset_index = e->Tset[j];
				Mnode *ee = MRG_vector[e->ngraph]->returnNodeByIndex(Tset_index);
				for (int k = 0; k < ee->child_number; k++)
				{
					int child_index = ee->child_index[k];
					NLIST.push_back(MRG_vector[ee->ngraph]->returnNodeByIndex(child_index));
				}
			}
		}
		else
		{
			for (int j = 0; j < e->child_number; j++)
			{
				int child_index = e->child_index[j];
				NLIST.push_back(MRG_vector[e->ngraph]->returnNodeByIndex(child_index));
			}
		}
	}
}

bool ElementComparator::HasSameRange(Mnode *m, Mnode *n)
{
	if (m->nlevel == n->nlevel)
	{
		if (m->nrange == n->nrange)
		{
			return true;
		}
		else return false;
	}
	else return false;
}

bool ElementComparator::HasParentPair(Mnode *m, Mnode *n)
{

	int left_parent_graph = m->ngraph;
	int left_parent_index = m->parent_index;
	int right_parent_graph = n->ngraph;
	int right_parent_index = n->parent_index;

	if (left_parent_index == -1 && right_parent_index == -1)
	{
		return true;
	}

	for (int i = 0; i < MPAIR.size(); i++)
	{
		vector<Mnode*> one_pair = MPAIR.at(i);
		if (InOnePair(one_pair, left_parent_graph, left_parent_index) && (one_pair, right_parent_graph, right_parent_index))
		{
			return true;
		}
	}
	return false;
}

bool ElementComparator::HasSameMLIST(Mnode *m, Mnode *n)
{
	if (m->MLIST.size() != n->MLIST.size())
	{
		return false;
	}
	else
	{
		for (int i = 0; i < m->MLIST.size(); i++)
		{
			if (m->MLIST.at(i) != n->MLIST.at(i))
			{
				return false;
			}
		}
		return true;
	}
}

bool ElementComparator::BelongToDifferentPart(Mnode *m, Mnode *n)
{
	if (m->ngraph == n->ngraph)
	{
		return false;
	}
	else return true;
}

bool ElementComparator::InOnePair(vector<Mnode*> one_pair, int node_graph, int node_index)
{
	for (int i = 0; i < one_pair.size(); i++)
	{
		Mnode *e = one_pair.at(i);
		if (e->fictitious_flag == 1)
		{
			for (int j = 0; j < e->Tset_number; j++)
			{
				Mnode *ee = MRG_vector[e->ngraph]->returnNodeByIndex(e->Tset[j]);
				if (ee->ngraph == node_graph)
				{
					if (ee->index == node_index)
					{
						return true;
					}
				}
			}
		}
		else
		{
			if (e->ngraph == node_graph)
			{
				if (e->index == node_index)
				{
					return true;
				}
			}
		}
	}
	return false;
}

void ElementComparator::BoardCastLabel(int index)
{
	out6 << "======START BOARDCAST======" << endl;
	vector<Mnode*> selected_node_pair = MPAIR.at(index);
	for (int i = 0; i < selected_node_pair.size(); i++)
	{
		Mnode *e = selected_node_pair.at(i);
		if (e->fictitious_flag == 1)
		{
			for (int j = 0; j < e->child_number; j++)
			{
				extendMLISTinOneDirection(NLIST[e->child_index[i]], 1, index);
				extendMLISTinOneDirection(NLIST[e->child_index[i]], -1, index);
			}
		}
		else
		{
			//upper
			extendMLISTinOneDirection(e, 1, index);
			//lower
			extendMLISTinOneDirection(e, -1, index);
		}
	}

}

void ElementComparator::extendMLISTinOneDirection(Mnode* e, int direcion, int index)
{
	if (e->ngraph == 0)
	{
		for (int i = 0; i < e->adjance_number; i++)
		{
			int ad_index = e->adjance_index[i];
			Mnode * adjance_node = this->R->returnNodeByIndex(ad_index);
			if (adjance_node->nrange == e->nrange + direcion)
			{
				if (adjance_node->visit_flag == -1)
				{
					adjance_node->MLIST.push_back(index);
				}
				extendMLISTinOneDirection(adjance_node, direcion, index);
			}
		}
	}
	else
	{
		for (int i = 0; i < e->adjance_number; i++)
		{
			int ad_index = e->adjance_index[i];
			Mnode * adjance_node = this->S->returnNodeByIndex(ad_index);
			if (adjance_node->nrange == e->nrange + direcion)
			{
				adjance_node->MLIST.push_back(index);
				extendMLISTinOneDirection(adjance_node, direcion, index);
			}
		}
	}
}

bool ElementComparator::isSameNode(Mnode * m, Mnode * n)
{
	if (m->ngraph == n->ngraph)
	{
		if (m->index == n->index)
		{
			return true;
		}
		else return false;
	}
	else return false;
}

void ElementComparator::getNodeInformation(Mnode * m)
{
	if (m->fictitious_flag == 0)
	{
		out6 << "======ONE NODE======" << endl;
		out6 << "GRAPH: " << m->ngraph << endl;
		out6 << "LEVEL: " << m->nlevel << endl;
		out6 << "RANGE: " << m->nrange << endl;
		out6 << "NUMBER: " << m->nnumber << endl;
		out6 << "P1: " << m->p1 << endl;
		out6 << "p2: " << m->p2 << endl;
		out6 << "======END======" << endl;
	}
	else
	{
		out6 << "======FICTITIOUS NODE======" << endl;
		out6 << "GRAPH: " << m->ngraph << endl;
		for (int i = 0; i < m->Tset_number; i++)
		{
			out6 << "====NODES UNDER FICTITIOUS NODE====" << endl;
			getNodeInformation(MRG_vector[m->ngraph]->returnNodeByIndex(m->Tset[i]));
		}
		out6 << "P1: " << m->p1 << endl;
		out6 << "p2: " << m->p2 << endl;
		out6 << "======END======" << endl;
	}
}

void ElementComparator::calculate_result()
{
	for (int i = 0; i < MPAIR.size(); i++)
	{
		double one_result = 0;
		Mnode* left_node = MPAIR[i][0];
		Mnode* right_node = MPAIR[i][1];
		one_result = sim(left_node, right_node);
		this->result += one_result;
	}
}
