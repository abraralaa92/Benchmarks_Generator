#include "generator.h"


#define LOG(x) std::cout << x << std::endl;

using namespace std;

//Static variables of class vertex()
unsigned int vertex::bal_unbal_flat;
unsigned int vertex::preMux_postMux_mixMux;
unsigned int vertex::lower_bound_length;
unsigned int vertex::upper_bound_length;

//global variables
unsigned int max_no_TDRS;

int main(int argc, char** argv) //https://www.tutorialspoint.com/cprogramming/c_command_line_arguments.htmargc 
//int main()
{
	//This part to initialize static variables
	vertex::bal_unbal_flat = 1;
	vertex::preMux_postMux_mixMux = 1;
	vertex::lower_bound_length = 5;
	vertex::upper_bound_length = 10;
	max_no_TDRS = 10000;
	
	//int argc = 1;
	//const char *inputs[] = { "", "2", "2", "0" }; //{type, path, initial_conf, write_read, UB-CSUs, Optimal_First} OR {type, path, hierarchy, children_per_SIB, initial_conf, write_read, UB-CSUs, Optimal_First}
	//const char** argv = inputs;
	
	unsigned int no_children_per_SIB, hierarchy_level, initial_configuration, SIBbased_MUXbased_Mix, expected_number_of_instruments;
	string path;

	if (argc > 1)
	{
		no_children_per_SIB = stoi(argv[1]);
		hierarchy_level = stoi(argv[2]);
		initial_configuration = stoi(argv[3]);
		SIBbased_MUXbased_Mix = 1;

		path = "./Benchmarks_" + to_string(max_no_TDRS) + "_TDRs/" + argv[1] + "_children_per_SIB/" + argv[2] + "_hierarchies/";

		vertex rootA(path, no_children_per_SIB, hierarchy_level, initial_configuration, SIBbased_MUXbased_Mix);
		rootA.generate_NW();
		printf("Generator Finish !!!!\n");

	}
	else 
	{
		initial_configuration = 0;
		SIBbased_MUXbased_Mix = 1;

		for (no_children_per_SIB = 1; no_children_per_SIB <= 10; no_children_per_SIB++)
		{
			for (hierarchy_level = 1; hierarchy_level <= 10; hierarchy_level++)
			{
				path = "./Benchmarks_" + to_string(max_no_TDRS) + "_TDRs/" + to_string(no_children_per_SIB) + "_children_per_SIB/" + to_string(hierarchy_level) + "_hierarchies/";

				//don't "Generate" this network if number of its TDRs exceed the limit
				expected_number_of_instruments = myPow(no_children_per_SIB, hierarchy_level);
				if (expected_number_of_instruments <= max_no_TDRS)
				{
					vertex rootA(path, no_children_per_SIB, hierarchy_level, initial_configuration, SIBbased_MUXbased_Mix);
					rootA.generate_NW();
				}
				else
					break; //in case I reached the max number of TDRs, then stop looping on (hierarchy_level) and move to the next (children_per_SIB)
			}
			printf("(%d) children per SIB has been generated with different level of hierarchies!!!!\n", no_children_per_SIB);
		}
	}
	
	/*
	char x;
	printf("please enter a char to exit");
	scanf("%c", &x);
	*/

	return 0;
}

