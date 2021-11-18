#ifndef Generator           //https://www.dreamincode.net/forums/topic/138027-struct-type-redefinition/
#define Generator           //https://www.youtube.com/watch?v=RU5JUHAiR18

#ifdef _MSC_VER				//https://www.acodersjourney.com/top-10-c-header-file-mistakes-and-how-to-fix-them/
#pragma once
#endif  // _MSC_VER

#include "structs.h"

extern unsigned int max_no_TDRS;

class vertex
{
	unsigned int index; //carry index of vertex in the "NW_vertices" vector.
	string id;
	unsigned int level;
	vector <unsigned int> predecessor; //we need both predecessor and successor pointers for SAT to check (SP_SS, SP_MS, MP_SS, MP_MS), (S,M,P,S): (single, multiple, predecesssor, successor) 
	vector <unsigned int> successor; //<unsigned int> these vector carry the Indices of successor and predeccor vertices; The trouble with using vector<vertex*> is that, whenever the vector goes out of scope unexpectedly(like when an exception is thrown), the vector cleans up after yourself, but this will only free the memory it manages for holding the pointer, not the memory you allocated for what the pointers are referring to.(https://stackoverflow.com/questions/1361139/how-to-avoid-memory-leaks-when-using-a-vector-of-pointers-to-dynamically-allocat/1361227)
	bool isStateElement; //is it SCB or CSU_register and not only a shift_register (here we assume that all SRs and TDRs are state elements). 
	unsigned int length;
	string selection_clause; //I need vector in order to accomulate on parent clauses 
	string type; //SI, SO, Br, AUX, MUX, SR, In, Out, SC
	string address_control_bit;

public:
	static unsigned int bal_unbal_flat; //static variables are only initialized once through the whole class with private accessing through only this class 
	static unsigned int preMux_postMux_mixMux;
	static unsigned int lower_bound_length;
	static unsigned int upper_bound_length;

	vertex(const string &directory, unsigned int no_children_per_SIB, unsigned int hierarchy_level, unsigned int initial_configuration, unsigned int SIBbased_MUXbased_Mix); //this constructor is called in the main function to set all the required parameters to generate a NW.
	~vertex();
	vertex(const string& id, unsigned int index, unsigned int level, const string& type, const string& clause, const string& address_control_bit, unsigned int length, bool isStateElement);
	vertex(const vertex& x);

	void reserve_vectors_memory(); //Very Important method used to minimze the number of copied objects and to adjust the capacity of used vectors for better performance
	void reset_vectors();

	unsigned int insert_pSIB(unsigned int predecessor_index, unsigned int level, const string& parent_clause, const string& next_level_addressControlBit);
	unsigned int insert_SIBp(unsigned int predecessor_index, unsigned int level, const string& parent_clause, const string& next_level_addressControlBit);
	void insert_SI_vertex();
	void insert_SO_vertex(unsigned int predecessor_index);
	void insert_Br_vertex(unsigned int predecessor_index, const string& id, unsigned int level, const string& selection_clause); //pVertex_index: parent vertex index
	void insert_AUX_vertex(unsigned int predecessor_index, const string& id, const string& selection_clause, const string& address_control_bit);
	void insert_TDR_vertex(unsigned int predecessor_index, const string& id, unsigned int length, const string& selection_clause, const string& type, const string& address_control_bit); //type could be ("In", "Out", "SC").
	void insert_I0_vertex(unsigned int predecessor_index, const string& id);
	void insert_I1_vertex(unsigned int predecessor_index, unsigned int successor_index, const string& id);
	void insert_SR_vertex(unsigned int predecessor_index, const string& id, const string& selection_clause, const string& address_control_bit);

	string get_vertexID(const string& id); //I can't return (string&) because I need to concatenate the returned string.
	string get_predecessor_smv_vertex(unsigned int index); //please note that for smv, SDG: the method is called "get_predecessor_..."/without (s), while for SAT: the method is called "get_predecessor(s)_..."; and this is because in smv and in SDG I need only one predecessor while in SAT I may have multiple predecessors and I need to return them all !!
	string& get_predecessor_SDG_vertex(unsigned int index);
	void set_SAT_vertex_predecessors(unsigned int index);
	void set_SAT_vertex_successors(unsigned int index);
	
	void generate_NW();
	void generate_NW_struct1();
	void generate_NW_clauses();
	void generate_NW_SDG();
	void generate_NW_smv();
	void generate_NW_graph();					//this "graph" would be traversed for structural_SAT objective.
	void generate_outputs();
	
	void print_vectors_size_capacity_comparison();
	void print_outputs();
	void print_file(_Type option, const string& input_file);
};


#endif //Generator          //https://en.wikipedia.org/wiki/Pragma_once