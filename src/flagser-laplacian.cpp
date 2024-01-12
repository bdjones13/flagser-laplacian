#include <fstream>
#include <iostream>

#define ASSEMBLE_REDUCTION_MATRIX
#define INDICATE_PROGRESS
#define SKIP_APPARENT_PAIRS
#define USE_ARRAY_HASHMAP
#define USE_CELLS_WITHOUT_DIMENSION
#define SORT_COLUMNS_BY_PIVOT
// #define WITH_HDF5
// #define USE_COEFFICIENTS
// #define MANY_VERTICES

#include "../include/parameters.h"
#include "../include/persistence_real.h"

//
// Compute directed flag complex homology
//

#include "../include/complex/directed_flag_complex_computer_real.h"
// #include "../include/complex/directed_flag_complex_in_memory_computer.h"

#include "../include/usage/flagser.h"

template <class T>
#ifdef RETRIEVE_PERSISTENCE
std::vector<persistence_computer_t<T>>
#else
void
#endif
compute_homology(filtered_directed_graph_t& graph, const flagser_parameters_real& params, const flagser_parameters& params_standard) {

	std::vector<filtered_directed_graph_t> subgraphs{graph};
	if (params.split_into_connected_components) { subgraphs = graph.get_connected_subgraphs(2); }

	// auto output = get_output<T>(params);
	size_t component_number = 1;
	for (auto subgraph : subgraphs) {
		T complex(subgraph, params_standard);

		// output->set_complex(&complex);
		if (params.split_into_connected_components) {
			// if (component_number > 1) output->print("\n");
			// output->print("## Path component number ");
			// output->print(std::to_string(component_number));
			// output->print("\n");

#ifdef INDICATE_PROGRESS
			std::cout << "\033[K";
#endif
			if (component_number > 1) std::cout << "\n";
			std::cout << "# Path component number " << component_number << std::endl;

			component_number++;
		}


		real_persistence_computer_t<decltype(complex)> persistence_computer(complex, 
									//output.get(),
								 	params.max_entries,
									// params.max_filtration,
									params.out_prefix);
		persistence_computer.compute_persistent_spectra(params.min_dimension, params.max_dimension);
        std::cout << "computing" << std::endl;
	}
}

int main(int argc, char** argv) {
	try {
        auto start_flagser_L = std::chrono::high_resolution_clock::now();
        
		auto arguments = parse_arguments(argc, argv);

		auto positional_arguments = get_positional_arguments(arguments);
		auto named_arguments = get_named_arguments(arguments);

		if (named_arguments.find("help") != named_arguments.end()) { print_usage_and_exit(-1); }

		auto params = flagser_parameters_real(named_arguments);
        auto params_standard = params.real_params_to_standard();
		if (positional_arguments.size() == 0) { print_usage_and_exit(-1); }
		const char* input_filename = positional_arguments[0];

		filtered_directed_graph_t graph = read_filtered_directed_graph(input_filename, params_standard); // pass custom params to graph so don't have to modify graph class

        compute_homology<real_directed_flag_complex_computer::real_directed_flag_complex_computer_t>(graph, params, params_standard);
        auto end_flagser_L = std::chrono::high_resolution_clock::now();
        auto duration_flagser_L = std::chrono::duration_cast<std::chrono::milliseconds>(end_flagser_L - start_flagser_L);
        std::cout << "duration flagser laplacian (ms):" << duration_flagser_L.count() << std::endl;
	} catch (const std::exception& e) { std::cout << e.what() << std::endl; }
}
