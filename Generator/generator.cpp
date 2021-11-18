#include "generator.h"

#define first_child 0
#define second_child 1
#define third_child 2

//static variables (private variables to the scope of this file)
//we use static and not make it a global variable because we may need to run only (generator.cpp) or (SAT.cpp) without running the other, that's why we choose the variable to be static and private to the scope of each class)
string path;
static unsigned int hierarchy_level;
static unsigned int no_children_per_SIB;
static unsigned int initial_configuration;
static unsigned int SIBbased_MUXbased_Mix;

static unsigned int module_no;
static unsigned int total_no_SIBs;           //this number represent number of network's SRs
static unsigned int total_no_instruments;    //this number represent number of network's (In, Out, SC)
static unsigned int Expected_no_nodes;

static vector<vertex> NW_vertices;
static vector<_selection_clause> NW_clauses;
static vector<_SDG> NW_SDG;
static vector<_smv> NW_smv;
static vector<_graph> NW_graph;
static vector<string> NW_TDRs;


vertex::vertex(const string& directory, unsigned int children_per_SIB, unsigned int level_of_hierarchy, unsigned int init_conf, unsigned int SIB_MUX_Mix)
{
	//we can't use "initializer_list" instead: https://stackoverflow.com/questions/18479295/member-initializer-does-not-name-a-non-static-data-member-or-base-class/18479406
	path = directory;
	hierarchy_level = level_of_hierarchy;
	no_children_per_SIB = children_per_SIB;
	initial_configuration = init_conf;
	SIBbased_MUXbased_Mix = SIB_MUX_Mix;
}

vertex::vertex(const string& id, unsigned int index, unsigned int level, const string& type, const string& clause, const string& address_control_bit = "", unsigned int length = 1, bool isStateElement = false)
	: index(index), level(level), type(type), selection_clause(clause), address_control_bit(address_control_bit), length((!isStateElement) ? 0 : length), isStateElement(isStateElement)	//length((!isStateElement) ? 0 : length) where if the vertex represent a multiplexer or a fanout element in the network(!isStateElement) then its a virtual element in the network with (0) length.
{
	//Note: copy constructor would be called in initilizer_list{: clause(clause)}, (https://stackoverflow.com/questions/9200827/why-does-passing-by-reference-involve-a-copy-constructor)

	if (type != "SI" && type != "SO")
		this->id = type + id + to_string(index); //adding (index) to "id" will work in making unique identifiers
	else
		this->id = id; //this condition for only TDI and TDO

	if (type == "SR")
		this->predecessor.reserve(2); //(2) because until now we are dealing only with "MUX(2*1)" in NW_generation, so the MUX's inputs or SR's predecessors willn't be more than two.
	else if (type == "Br")
		this->successor.reserve(2); //(2) for the same reason where (Br) node indicate that there will be branching next, and this until "now" is just related to MUX's inputs, which are only two.
}

vertex::vertex(const vertex& x)
	: index(x.index), id(x.id), level(x.level), predecessor(x.predecessor), successor(x.successor), isStateElement(x.isStateElement), length(x.length), selection_clause(x.selection_clause), type(x.type), address_control_bit(x.address_control_bit)
{
	printf("Vertex copied!!\n");
}

vertex::~vertex()
{
	//check for any memory leaks
	std::vector<unsigned int>().swap(predecessor);
	std::vector<unsigned int>().swap(successor);

	reset_vectors();
}

void vertex::reserve_vectors_memory()
{
	//reset (module_no) in case it was used inside a loop, to stop increament on old values!!
	module_no = 0;

	//here I did only a mathematical prediction for the number of (NW_vertices, _SCT, NW_SCT_to_SAT_predicates, NW_SAT_clauses)
	//however other vectors' sizes are predicted using try and error!!
	//mathimatical predictions:
	total_no_SIBs = 0;           //this number represent number of network's SRs
	total_no_instruments = 0;    //this number represent number of network's (In, Out, SC)

	for (unsigned int l = 1; l <= hierarchy_level; l++)
		total_no_SIBs += myPow(no_children_per_SIB, l);
	total_no_instruments = myPow(no_children_per_SIB, hierarchy_level);
	Expected_no_nodes = (total_no_SIBs * 5) + total_no_instruments + 2;  //(*5)--> because "each SIB" in the network contains 5 nodes (Br, A, I0, I1, SR) except the SIBs of the last hierarchy needs one more node for (in, Out, SC) nodes, (+2)--> for (SI, SO) nodes

	NW_vertices.reserve(Expected_no_nodes);
	NW_clauses.reserve(Expected_no_nodes / 3);
	NW_SDG.reserve(Expected_no_nodes);
	NW_smv.reserve(Expected_no_nodes / 2);
	NW_graph.reserve(Expected_no_nodes * 2 / 3);
	NW_TDRs.reserve(total_no_instruments);
}

void vertex::reset_vectors()
{
	/*
	https://www.techiedelight.com/delete-vector-free-memory-cpp/
	If you want to reset you vector back to a empty state then we can use the swap trick to swap the contents of the vector into a temporary that will get destroyed and free the memory
	This will create a temporary empty vector, swap it with the one you have so now it is empty and then destroy everything in the temporary vector.
	std::vector<string>().swap(NW_SAT_clauses);
	OR
	NW_SAT_clauses.clear();
	NW_SAT_clauses.shrink_to_fit();
	*/

	//release the memory associated with each vector 
	//this method should not be called except if I was looping and I need to reset old vectors
	std::vector<vertex>().swap(NW_vertices);
	std::vector<_selection_clause>().swap(NW_clauses);
	std::vector<_SDG>().swap(NW_SDG);
	std::vector<_smv>().swap(NW_smv);
	std::vector<_graph>().swap(NW_graph);
	std::vector<string>().swap(NW_TDRs);
}

void vertex::generate_NW()
{
	if ((bal_unbal_flat == 1) && (preMux_postMux_mixMux == 1) && (SIBbased_MUXbased_Mix == 1))
		generate_NW_struct1();

	generate_outputs();
	print_vectors_size_capacity_comparison();
}

void vertex::generate_NW_struct1()
{
	reserve_vectors_memory(); //This way of vector's reserve prediction is only associated with balanced networks (NW_struct1), use 'reserve' keyword for better performance to minimize number of vectors expansions

	insert_SI_vertex();
	unsigned int SR_index = insert_pSIB(NW_vertices.size() - 1, 1, "TRUE", "");	//we start from level (1)
	insert_SO_vertex(SR_index);
}

unsigned int vertex::insert_pSIB(unsigned int predecessor_index, unsigned int level, const string& parent_clause, const string& address_control_bit)
{
	//I have three main while loops
	//1- while (no of modules < no_children_per_SIB)
	//2- while (required herarchy level is not reached yet)
	//3- while (no of children for max hierachy level < no_children_per_SIB)
	//however, applying recursion will be more effective than using the while loop

	//Forward path: check on the level of hierarchy
	//Backward path: check on the no. of children
	unsigned int SIB_children;
	unsigned int last_level_children;
	unsigned int range, length;
	string id;

	unsigned int next_level;
	string next_level_addressControlBit, next_level_clause;

	//there is essential difference between "pBr_index" and "Br_index"; the former refers to the index parent vertex of "Br"/predecessor vertex, while the last refers to the "Br" vertex itself
	unsigned int pBr_index; //pBr_index: the index of "Br" vertex of the Previous_level, "p" for Parent, where pBr_vertex appears in the graph order as the Parent for the current "Br" vertex
	unsigned int Br_index;  //Br_index:  the index of "Br" vertex of the Current_level, and it is used to set/save the index of "Br" to use it in next level vertices generation
	unsigned int SR_index;  //SR_index:  the index of "SR" vertex of the Current_level
	unsigned int pSR_index; //pSR_index:  the index of "SR" vertex of the Previous_level, I need this index after the completion of generating all vertices of a deep level and when I return back to some of the lower levels

	if (level < hierarchy_level)
	{
		SIB_children = 0;
		pBr_index = predecessor_index;

		while (SIB_children < no_children_per_SIB)
		{
			if (level == 1) //means we are in first level since the next is level(2)
				module_no++;

			//update "id" before any "Br" vertex insertion
			id = "-M" + to_string(module_no) + "_" + to_string(level) + "/";  //M refers to Module, "/" instead of "." to be suitable also with Rene's code where I was having an error in (SR-M1_1.1.SDG), where there was a conflict with other keywords like ".o, .out, .SI"

			//I have to adjust only the level of "Br_" node with every new call, because all other subsequent nodes for the same SIB are deriving the same "Br_level"
			insert_Br_vertex(pBr_index, id, level, parent_clause);
			Br_index = NW_vertices.size() - 1;

			//prepare next level clauses
			next_level = level + 1;
			next_level_addressControlBit = "SR" + id + to_string(Br_index + 3);
			next_level_clause = parent_clause == "TRUE" ? next_level_addressControlBit : next_level_addressControlBit + "^" + parent_clause; //where SR_index = Br_index + 3

			insert_AUX_vertex(Br_index, id, "!" + next_level_clause, next_level_addressControlBit); //size(NW_vertices) - 1: Br_index
			insert_I0_vertex(Br_index + 1, id); //size(NW_vertices) - 1: AUX_index
			insert_SR_vertex(Br_index + 2, id, parent_clause, address_control_bit);
			SR_index = Br_index + 3;		//(+3) because SR_vertex comes after three insertions from inserting BR_vertex [AUX, I0, SR]

			SIB_children++;

			pSR_index = insert_pSIB(Br_index, next_level, next_level_clause, next_level_addressControlBit);

			//AFTER "RETURN" from RECURSION
			insert_I1_vertex(pSR_index, SR_index, id);
			pBr_index = SR_index;
		}
		return SR_index; //I need this index to be used finally to assign last_SR_index to "SO" node
	}
	else if (level == hierarchy_level)
	{
		last_level_children = 0;
		pBr_index = predecessor_index;

		while (last_level_children < no_children_per_SIB)
		{
			if (hierarchy_level == 1)
				module_no++;

			//range = myPow(2, upper_bound_length) - myPow(2, lower_bound_length) + 1;
			//length = rand() % range + myPow(2, lower_bound_length); //this eqaution to get the random number using "specific number of Bits" (10 bits, 5 bits) or we can use: rand() % (int)((pow(2, upper_bound_length) - pow(2, lower_bound_length)) + 1) + pow(2, lower_bound_length); for specific value (max, min)
			length = 20;

			//update only "id", index to the parent of Br node "pBr_index", without any change in "next_level" since all last level vertices inherit the same final "next_level" 
			id = "-M" + to_string(module_no) + "_" + to_string(level) + "/";  //M refers to Module, "/" instead of "." to be suitable also with Rene's code where I was having an error in (SR-M1_1.1.SDG), where it may conflicts with other keywords like ".o, .out, .SI"

			insert_Br_vertex(pBr_index, id, level, parent_clause);
			Br_index = NW_vertices.size() - 1;

			next_level_addressControlBit = "SR" + id + to_string(Br_index + 3);
			next_level_clause = parent_clause == "TRUE" ? next_level_addressControlBit : next_level_addressControlBit + "^" + parent_clause;

			insert_AUX_vertex(Br_index, id, "!" + next_level_clause, next_level_addressControlBit);
			insert_I0_vertex(Br_index + 1, id);
			insert_SR_vertex(Br_index + 2, id, parent_clause, address_control_bit);
			SR_index = Br_index + 3;

			last_level_children++;

			if (last_level_children == 1) //means this is the first child
				insert_TDR_vertex(Br_index, id, length, next_level_clause, "In", next_level_addressControlBit); //(In, Out, SC) all have the same level as "Br" vertex, size(NW_vertices) - 2: since  have (Br->AUX->In) so the difference to "Br" vertex is two previous pushes

			else if (last_level_children == 2)
				insert_TDR_vertex(Br_index, id, length, next_level_clause, "Out", next_level_addressControlBit);

			else
				insert_TDR_vertex(Br_index, id, length, next_level_clause, "SC", next_level_addressControlBit);

			pSR_index = NW_vertices.size() - 1; //"size(NW_vertices) - 1" represent index of (In, Out, SC)
			insert_I1_vertex(pSR_index, SR_index, id);
			pBr_index = SR_index;
		}

		return SR_index;
	}

	return -1; //no need for this, this is only to overcome the warning message "not all control paths return a value"
}

void vertex::insert_SI_vertex()
{
	NW_vertices.emplace_back("TDI", NW_vertices.size(), 1, "SI", "TRUE");
}

void vertex::insert_SO_vertex(unsigned int predecessor_index)
{
	NW_vertices.emplace_back("TDO", NW_vertices.size(), 1, "SO", "TRUE");
	//adjust predecessor of "SO" to point to "SO" vertex
	NW_vertices[predecessor_index].successor.emplace_back(NW_vertices.size() - 1); //size(NW_vertices)-1 represent last index
	// set predecessor pointer of the "SO" vertex to point to the previous node / predecessor_vertex
	NW_vertices.back().predecessor.emplace_back(predecessor_index);
}

void vertex::insert_Br_vertex(unsigned int predecessor_index, const string& id, unsigned int Br_level, const string& selection_clause)
{
	NW_vertices.emplace_back(id, NW_vertices.size(), Br_level, "Br", selection_clause);
	NW_vertices[predecessor_index].successor.emplace_back(NW_vertices.size() - 1);
	NW_vertices.back().predecessor.emplace_back(predecessor_index);
}

void vertex::insert_AUX_vertex(unsigned int predecessor_index, const string& id, const string& selection_clause, const string& address_control_bit)
{
	//one more step is to adjust (address_control_bit) of each "SR" vertices
	NW_vertices.emplace_back(id, NW_vertices.size(), NW_vertices[predecessor_index].level, "AUX", selection_clause, address_control_bit);
	NW_vertices[predecessor_index].successor.emplace_back(NW_vertices.size() - 1);
	NW_vertices.back().predecessor.emplace_back(predecessor_index);
}

void vertex::insert_TDR_vertex(unsigned int predecessor_index, const string& id, unsigned int length, const string& selection_clause, const string& type, const string& address_control_bit)
{
	//one more step is to adjust (address_control_bit) of each "SR" vertices
	NW_vertices.emplace_back(id, NW_vertices.size(), NW_vertices[predecessor_index].level, type, selection_clause, address_control_bit, length, true);
	NW_vertices[predecessor_index].successor.emplace_back(NW_vertices.size() - 1);
	NW_vertices.back().predecessor.emplace_back(predecessor_index);

	NW_TDRs.emplace_back(type + id + to_string(NW_vertices.size() - 1));
}

void vertex::insert_I0_vertex(unsigned int predecessor_index, const string& id)
{
	NW_vertices.emplace_back("[I0]" + id, NW_vertices.size(), NW_vertices[predecessor_index].level, "MUX", NW_vertices[predecessor_index].selection_clause);
	NW_vertices[predecessor_index].successor.emplace_back(NW_vertices.size() - 1);
	NW_vertices.back().predecessor.emplace_back(predecessor_index);
}

void vertex::insert_I1_vertex(unsigned int predecessor_index, unsigned int successor_index, const string& id)
{
	NW_vertices.emplace_back("[I1]" + id, NW_vertices.size(), NW_vertices[predecessor_index].level, "MUX", NW_vertices[predecessor_index].selection_clause);

	//adjust successor of "predecessor" node
	NW_vertices[predecessor_index].successor.emplace_back(NW_vertices.size() - 1);

	//adjust predecessor and successor of "I1"
	NW_vertices.back().predecessor.emplace_back(predecessor_index);
	NW_vertices.back().successor.emplace_back(successor_index);

	//adjust Second predecessor of "SR"
	NW_vertices[successor_index].predecessor.emplace_back(NW_vertices.size() - 1);
}

void vertex::insert_SR_vertex(unsigned int predecessor_index, const string& id, const string& selection_clause, const string& address_control_bit)
{
	NW_vertices.emplace_back(id, NW_vertices.size(), NW_vertices[predecessor_index].level, "SR", selection_clause, address_control_bit, 1, true); //"" for empty address_control_bit
	NW_vertices[predecessor_index].successor.emplace_back(NW_vertices.size() - 1);
	NW_vertices.back().predecessor.emplace_back(predecessor_index);
}

void vertex::generate_NW_clauses()
{
	for (size_t index = 0, e = NW_vertices.size(); index < e; index++)  //no need to call size() every iteration, //https://en.cppreference.com/w/cpp/types/size_t ,//https://stackoverflow.com/questions/4849678/c-for-loop-size-type-vs-size-t
	{
		if (NW_vertices[index].isStateElement)
			NW_clauses.emplace_back(NW_vertices[index].id, NW_vertices[index].length, initial_configuration, NW_vertices[index].selection_clause); //clause vector will be updated with a vector of clauses through (split_selection_clause_into_vectorOfClauses) method, which is called through the structure's constructor. Initial_configuration represent (initial state and reset state)
	}
}

void vertex::generate_NW_SDG()
{
	string SDG_predecessors = "";
	for (size_t i = NW_vertices.size() - 1; i > 0; i--) //start from SO
	{
		if (NW_vertices[i].type != "AUX")
		{
			//successor in SDG is the predecessor in NW-connection
			if (NW_vertices[i].predecessor.size() > 1) //like "SR" vertex
			{
				SDG_predecessors += get_predecessor_SDG_vertex(NW_vertices[i].predecessor[first_child]);
				SDG_predecessors += ", " + get_predecessor_SDG_vertex(NW_vertices[i].predecessor[second_child]);
				NW_SDG.emplace_back(NW_vertices[i].id, SDG_predecessors, NW_vertices[i].address_control_bit, NW_vertices[i].length);
			}
			else
			{
				SDG_predecessors += get_predecessor_SDG_vertex(NW_vertices[i].predecessor[first_child]);
				NW_SDG.emplace_back(NW_vertices[i].id, SDG_predecessors, NW_vertices[i].address_control_bit, NW_vertices[i].length);
			}

			SDG_predecessors.clear();
		}
		//else: continue searching for SDG/"non AUX" vertex

	}
	//for i==0; "SI" node there is no predecessor vertices, so I need to push an "empty" value
	NW_SDG.emplace_back(NW_vertices[0].id, SDG_predecessors, "", 1);
}

void vertex::generate_NW_smv()
{
	string selectControlInput;
	for (size_t i = 0, e = NW_vertices.size(); i < e; i++)
	{
		if (NW_vertices[i].isStateElement) //(SR, In, Out, SC)
		{
			selectControlInput = (NW_vertices[i].selection_clause == "TRUE") ? "TRUE" : "sel_" + NW_vertices[i].id + ".o";

			if (NW_vertices[i].type == "SR") //CSU Register
			{
				NW_smv.emplace_back(NW_vertices[i].id, "m-" + get_vertexID(NW_vertices[i].id) + ".o", "SUC", "", selectControlInput);
				NW_smv.emplace_back("m-" + get_vertexID(NW_vertices[i].id), get_predecessor_smv_vertex(NW_vertices[i].predecessor[first_child]) + ", " + get_predecessor_smv_vertex(NW_vertices[i].predecessor[second_child]), "MUX", NW_vertices[i].address_control_bit + ".ToSel", "");
			}
			else //Shift Register
			{
				NW_smv.emplace_back(NW_vertices[i].id, get_predecessor_smv_vertex(NW_vertices[i].predecessor[first_child]), "ShiftRegister_" + to_string(NW_vertices[i].length), "", selectControlInput);
			}
		}
	}
}

string vertex::get_vertexID(const string& id)
{
	return id.substr(id.find("-") + 1); //get after "-" to the end; 
}

string vertex::get_predecessor_smv_vertex(unsigned int index)
{
	string IDwE; //ID with extension
	while ((!NW_vertices[index].isStateElement) && (NW_vertices[index].type != "SI"))
		index = NW_vertices[index].predecessor[first_child];

	//it means for this special index "NW_vertices[index]" is either a (stateElement) or (SI) vertex 
	if (NW_vertices[index].isStateElement)
		IDwE = NW_vertices[index].id + ".SDO";
	else
		IDwE = "SG.TDI";

	return IDwE;
}

string& vertex::get_predecessor_SDG_vertex(unsigned int index)
{
	while (NW_vertices[index].type == "AUX")
		index = NW_vertices[index].predecessor[first_child]; //here I just care about the first child; because it is expexcted in NW_struct1 to have only one parent for these passes IDs.

	return NW_vertices[index].id;
}

void vertex::set_SAT_vertex_predecessors(unsigned int index)
{
	unsigned int predecessor_index;
	//here I need to use for loop and not (while) as in (get_predecessor_smv_vertex) method; because I need here to apply recursion on each predecessor child and not to take only the "first_child) as in (get_predecessor_smv_vertex)
	for (size_t i = 0, e = NW_vertices[index].predecessor.size(); i < e; i++)
	{
		predecessor_index = NW_vertices[index].predecessor[i];
		if ((!NW_vertices[predecessor_index].isStateElement) && (NW_vertices[predecessor_index].type != "AUX") && (NW_vertices[predecessor_index].type != "SI"))
			set_SAT_vertex_predecessors(predecessor_index);
		else
			NW_graph.back().predecessors.emplace_back(NW_vertices[predecessor_index].id);
	}
}

void vertex::set_SAT_vertex_successors(unsigned int index)
{
	unsigned int successor_index;
	for (size_t i = 0, e = NW_vertices[index].successor.size(); i < e; i++)
	{
		successor_index = NW_vertices[index].successor[i];
		if ((!NW_vertices[successor_index].isStateElement) && (NW_vertices[successor_index].type != "AUX") && (NW_vertices[successor_index].type != "SO"))
			set_SAT_vertex_successors(successor_index);
		else
			NW_graph.back().successors.emplace_back(NW_vertices[successor_index].id);
	}
}

void vertex::generate_NW_graph()
{
	//generate_NW_graph from the general (NW_vertices) graph constructed by the NW_generator.

	for (size_t index = 0, e = NW_vertices.size(); index < e; index++)
	{
		if (NW_vertices[index].type != "MUX" && NW_vertices[index].type != "Br") //(SR, In, Out, SC), AUX, SI, SO 
		{
			NW_graph.emplace_back(NW_vertices[index].id, NW_vertices[index].selection_clause, NW_vertices[index].address_control_bit, NW_vertices[index].length);
			NW_graph.back().predecessors.reserve(2);					//(2) because for (2*1) MUXes we have only two inputs, and node's predecessor can either (two) if it drives MUX's inputs or (one) for the remaining nodes. here I want to reserve these vectors only once, and then reuse them after clearing it each loop, because of that I reserved the max size that I could ever need.
			NW_graph.back().successors.reserve(hierarchy_level + 1);	//(hierarchy_level + 1) because the "max" branching which could ever happen, occurs when the node get branched into number of AUX nodes, each associated with a MUX from different hierarchy_level and one more branch for the (In/Out/SC). That's why I shouldn't need more than (hierarchy_level + 1) by max.

			set_SAT_vertex_successors(index);
			set_SAT_vertex_predecessors(index);
		}
	}
}

void vertex::generate_outputs()
{
	generate_NW_clauses();	//used to build the SCT
	generate_NW_SDG();		//used to build the Active Scan Path
	generate_NW_smv();		//used as an input to SAT_solution implemented by "Rene"
	generate_NW_graph();	//used to build SAT_clauses for SAT_solution implemented by "Me"

	print_outputs();
}

void vertex::print_outputs()
{
	//0-(NOT NEEDED) NW clauses associated with structured retargetting
	//print_file(NW_Vertices, "NW_nodes.txt");
	//1- NW clauses associated with structured retargetting
	print_file(NW_Clauses, "NW_clauses.txt");
	//2- NW SDG connections associated with structured retargetting
	//print_file(NW_SelDepGraph, "NW_SDG.txt");
	//3- NW smv associated with Rene's executable file
	print_file(NW_SMV, "NW_smv.smv");
	//4- NW graph to be used in SAT retargeting, where predicates are constructed from NW_graph
	print_file(NW_Graph, "NW_graph.txt");
	//5- NW TDRs to be used in SAT retargeting, in case I needed to loop on all NW TDrs
	//print_file(NW_Instruments, "NW_TDRs.txt");
}

void vertex::print_file(_Type option, const string& input_file)
{
	FILE *file;
	string output_text = "";
	string file_name = path + input_file;

	if ((file = fopen(file_name.c_str(), "w")) != NULL)
	{
		switch (option)
		{
		case NW_Vertices:
		{
			string selection_clause;

			output_text += "{id, successors, level, isStateElement, length, selection_clause, type, address_control_bit}\n\n{";
			for (size_t i = 0, e1 = NW_vertices.size(); i < e1; i++)
			{
				output_text += "{" + NW_vertices[i].id + ", {";
				for (size_t j = 0, e2 = NW_vertices[i].successor.size(); j < e2; j++)
				{
					if (j != 0)
						output_text += ", ";
					output_text += NW_vertices[NW_vertices[i].successor[j]].id;
				}
				output_text += "}, " + to_string(NW_vertices[i].level) + ", " + to_string(NW_vertices[i].isStateElement) + ", " + to_string(NW_vertices[i].length) + ", \"" + NW_vertices[i].selection_clause + "\", " + NW_vertices[i].type + ", " + NW_vertices[i].address_control_bit;

				if (i != e1 - 1)
					output_text += "},\n";
			}
			output_text += "}};";
			break;
		}
		case NW_Clauses:
		{
			output_text += to_string(NW_clauses.size() - NW_TDRs.size()) + "\n";  //this number is used while re-generating the NW_clauses to set the "NW_stateElements"
			output_text += to_string(NW_TDRs.size()) + "\n{";  //this number is used while re-generating the NW_clauses to set the "NW_selectElements"

			//{string reg_id, string selectionClause, int reg_len, int reset_val}
			//{ "SCB_8", "SCB_10 ^ !SCB_9", 1, 0 },
			//PLEASE note that: for "NW_struct1"/"SIB_based" networks, we have no any clause's negation like (!SR_2/1)
			for (size_t i = 0, e1 = NW_clauses.size(); i < e1; i++)
			{
				output_text += "{ \"" + NW_clauses[i].id + "\", \"";
				for (size_t j = 0, e2 = NW_clauses[i].clause.size(); j < e2; j++)
				{
					output_text += NW_clauses[i].clause[j].clause;
					if (j != e2 - 1)
						output_text += "^";

				}
				output_text += "\", " + to_string(NW_clauses[i].length) + ", " + to_string(NW_clauses[i].reset_value);
				if (i != e1 - 1)
					output_text += " },\n";
			}
			output_text += " }};";
			break;
		}
		case NW_SelDepGraph:
		{
			//(TDI, TDO) instead of (SI, SO) and there is no AUX nodes
			//{string reg_id, vector<string> next_OrderedNodes, string selection_control, int length}
			//{ "TDO", { "SCB_10" }, "", 1 },
			//{ "I3",{ "I0_SCB_3", "I1_SCB_3" }, "SCB_3", 20 },
			//{ "TDI",{}, "", 1 },

			output_text += "{";
			for (size_t i = 0, e = NW_SDG.size(); i < e; i++)
			{
				output_text += "{ \"" + NW_SDG[i].id + "\", { \"" + NW_SDG[i].successors + "\" }, \"" + NW_SDG[i].address_control_bit + "\", " + to_string(NW_SDG[i].length);

				if (i != e - 1)
					output_text += " },\n";
			}
			output_text += " }};";
			break;
		}
		case NW_SMV:
		{
			//here I have to generate two files (smv, pdl)
			output_text += "MODULE main\n";
			output_text += "VAR\n";
			output_text += "SG : SigGen()\n\n";
			output_text += "TDO : SIGNAL (" + NW_vertices[NW_vertices[NW_vertices.size() - 1].predecessor[0]].id + ".SDO);\n\n";

			for (size_t i = 0, e = NW_smv.size(); i < e; i++)
			{
				// SR-M1_1/1: SUC(m-M1_1/1.o, SG.SCK, TRUE, SG.ShiftEn, SG.UpdateEn); 
				if (NW_smv[i].type == "SUC")
					output_text += NW_smv[i].id + ": SUC(" + NW_smv[i].predecessors + ", SG.SCK, " + NW_smv[i].selectControlInput + ", SG.ShiftEn, SG.UpdateEn);\n";

				//m-M1_1/1: MUX(SG.TDI, SR-M1_2/2.SDO, SR-M1_1/1.ToSel);
				else if (NW_smv[i].type == "MUX")
					output_text += NW_smv[i].id + ": MUX" + to_string(count(NW_smv[i].predecessors.begin(), NW_smv[i].predecessors.end(), ',') + 1) + "(" + NW_smv[i].predecessors + ", " + NW_smv[i].address_control_bit + ");\n";

				//In - M1_2 / 1: ShiftRegister_10(SG.TDI, SG.SCK, TRUE, SG.ShiftEn, SG.UpdateEn);
				else //NW_smv[i].type == "ShiftRegister"
					output_text += NW_smv[i].id + ": " + NW_smv[i].type + "(" + NW_smv[i].predecessors + ", SG.SCK, " + NW_smv[i].selectControlInput + ", SG.ShiftEn, SG.UpdateEn);\n";
			}
			//I need to write the pdl file also
			//iWrite In-M1_2/1 0b0
			FILE *second_file;
			string second_file_name = path + input_file.substr(0, input_file.find(".")) + ".pdl";
			if ((second_file = fopen(second_file_name.c_str(), "w")) != NULL)
			{
				fprintf(second_file, "iWrite %s 0b0\n", NW_TDRs[NW_TDRs.size() - 1].c_str()); //here we set initially (target_reg) in 'pdl' file to be the Last/Deepest TDR in the chain with the longest retargeting path.
				fclose(second_file);
			}
			else printf("%s \n", "Unable to open pdl file");

			break;
		}
		case NW_Graph:
		{
			output_text += to_string(NW_graph.size()) + "\n{";  //this number is used while re-generating NW_graph from the input file in [SAT.cpp/ reBuild_NW_graph()] method
			//{string id, string selection_clause, string address_control_bit, unsigned int length, vector <string> predecessor, vector <string> successor}

			for (size_t i = 0, e1 = NW_graph.size(); i < e1; i++)
			{
				output_text += "{ \"" + NW_graph[i].id + "\", \"" + NW_graph[i].selection_clause + "\", \"" + NW_graph[i].address_control_bit + "\", " + to_string(NW_graph[i].length) + ", {";

				for (size_t j = 0, e2 = NW_graph[i].predecessors.size(); j < e2; j++)
				{
					output_text += NW_graph[i].predecessors[j];
					if (j != e2 - 1)
						output_text += ",";
				}
				output_text += "}, {";
				for (size_t j = 0, e2 = NW_graph[i].successors.size(); j < e2; j++)
				{
					output_text += NW_graph[i].successors[j];
					if (j != e2 - 1)
						output_text += ",";
				}
				output_text += "}";

				if (i != e1 - 1)
					output_text += " },\n";
			}
			output_text += " }};";
			break;
		}
		case NW_Instruments:
		{
			output_text += to_string(NW_TDRs.size()) + "\n";

			for (size_t i = 0, e1 = NW_TDRs.size(); i < e1; i++)
				output_text += NW_TDRs[i] + "\n";
			break;
		}
		}
		fprintf(file, "%s", output_text.c_str());
		fclose(file);
	}
	else printf("%s \n", "Unable to open file");
}

void vertex::print_vectors_size_capacity_comparison()
{
	FILE *file;
	string file_name = path + "NW_vectors_sizes.txt"; //NW_vectors_sizes has important information like (size of NW_graph, size of NW_clauses, hierarchy_level, ..)

	if ((file = fopen(file_name.c_str(), "w")) != NULL)		//  "w" to overwrite, previous content will be deleted. https://stackoverflow.com/questions/2393345/how-to-append-text-to-a-text-file-in-c
	{
		fprintf(file, "NO of Children per SIB = %d\n", no_children_per_SIB);
		fprintf(file, "Hierarchy = %d\n", hierarchy_level);
		fprintf(file, "NO of TDRs = %d\n", total_no_instruments);

		//here for capacity I can't say "Capacity= %d, NW_vertices.capacity()" because it will return the current capacity. However, I want the capacity which I reserved initially!
		fprintf(file, "For NW_vertices: Size= %lu,\t Capacity= %d\n", NW_vertices.size(), Expected_no_nodes); //where .size() is of type: long unsigned int (https://stackoverflow.com/questions/4033018/how-to-print-an-unsigned-long-int-with-printf-in-c/4033039)
		fprintf(file, "For NW_graph: Size= %lu,\t Capacity= %d\n", NW_graph.size(), (Expected_no_nodes * 2 / 3));
		fprintf(file, "For NW_clauses: Size= %lu,\t Capacity= %d\n", NW_clauses.size(), (Expected_no_nodes / 3));
		fprintf(file, "For NW_SDG: Size= %lu,\t Capacity= %d\n", NW_SDG.size(), Expected_no_nodes);
		fprintf(file, "For NW_smv: Size= %lu, \t Capacity= %d\n", NW_smv.size(), (Expected_no_nodes / 2));
		fprintf(file, "For NW_TDRs: Size= %lu, \t Capacity= %d\n", NW_TDRs.size(), total_no_instruments);
		fprintf(file, "////////////////////////////////////////////////////////////////////////////////////\n");

		/*
		printf("For NW_vertices: Reservation = %s\n", (Expected_no_nodes >= NW_vertices.size()) ? "True" : "FALSE!!!!!!!!!");
		printf("For NW_clauses: Reservation = %s\n", (Expected_no_nodes / 3 >= NW_clauses.size()) ? "True" : "FALSE!!!!!!!!!");
		printf("For NW_SDG: Reservation = %s\n", (Expected_no_nodes >= NW_SDG.size()) ? "True" : "FALSE!!!!!!!!!");
		printf("For NW_smv: Reservation = %s\n", (Expected_no_nodes / 2 >= NW_smv.size()) ? "True" : "FALSE!!!!!!!!!");
		printf("For NW_graph: Reservation= %s\n", (Expected_no_nodes / 2)>= NW_graph.size()) ? "True" : "FALSE!!!!!!!!!");
		printf("////////////////////////////////////////////////////////////////////////////////////\n");
		*/

		fclose(file);
	}
	else printf("%s \n", "Unable to open pdl file");
}