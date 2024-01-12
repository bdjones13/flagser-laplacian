#pragma once

#include <cassert>
#include <deque>
#include <iostream>
#include <queue>
#include <set>
// #include <eigen3/Eigen/src/Core/IO>
#include <eigen3/unsupported/Eigen/SparseExtra>
#include "definitions.h"

#include <chrono>

#include <PersistentLaplacians>

struct real_entry_t {
	index_t index;
	real_coefficient_t coefficient;
	real_entry_t(index_t _index, real_coefficient_t _coefficient) : index(_index), coefficient(_coefficient) {}
	real_entry_t(index_t _index) : index(_index), coefficient(1) {}
	real_entry_t() : index(0), coefficient(1) {}
};

real_entry_t real_make_entry(index_t _index, real_coefficient_t _coefficient) { return real_entry_t(_index, _coefficient); }
index_t get_index(real_entry_t e) { return e.index; }
real_coefficient_t get_coefficient(real_entry_t e) { return e.coefficient; } // changed type index_t to real_coefficient_t 2023-03-03 Ben Jones
void set_coefficient(real_entry_t& e, const real_coefficient_t c) { e.coefficient = c; }
bool operator==(const real_entry_t& e1, const real_entry_t& e2) {
	return get_index(e1) == get_index(e2) && get_coefficient(e1) == get_coefficient(e2);
}

std::ostream& operator<<(std::ostream& stream, const real_entry_t& e) {
	stream << get_index(e) << ":" << get_coefficient(e);
	return stream;
}

const real_entry_t& get_entry(const real_entry_t& e) { return e; }

// template <typename Entry> struct smaller_index {
// 	bool operator()(const Entry& a, const Entry& b) { return get_index(a) < get_index(b); }
// };

class real_filtration_index_t : public std::pair<value_t, index_t> {
public:
	real_filtration_index_t() : std::pair<value_t, index_t>() {}
	real_filtration_index_t(std::pair<value_t, index_t> p) : std::pair<value_t, index_t>(p) {}
};
value_t get_filtration(real_filtration_index_t i) { return i.first; }
index_t get_index(real_filtration_index_t i) { return i.second; }

class real_filtration_entry_t : public std::pair<value_t, real_entry_t> {
public:
	real_filtration_entry_t(std::pair<value_t, real_entry_t> p) : std::pair<value_t, real_entry_t>(p) {}
	real_filtration_entry_t(real_entry_t e) : std::pair<value_t, real_entry_t>(value_t(0), e) {}
	real_filtration_entry_t() : real_filtration_entry_t(0) {}
	real_filtration_entry_t(value_t _filtration, index_t _index, real_coefficient_t _coefficient)
	    : std::pair<value_t, real_entry_t>(_filtration, real_make_entry(_index, _coefficient)) {}
	real_filtration_entry_t(real_filtration_index_t _filtration_index, real_coefficient_t _coefficient)
	    : std::pair<value_t, real_entry_t>(get_filtration(_filtration_index),
	                                  real_make_entry(get_index(_filtration_index), _coefficient)) {}
	real_filtration_entry_t(real_filtration_index_t _filtration_index) : real_filtration_entry_t(_filtration_index, 1) {}
};

const real_entry_t& get_entry(const real_filtration_entry_t& p) { return p.second; }
real_entry_t& get_entry(real_filtration_entry_t& p) { return p.second; }
index_t get_index(const real_filtration_entry_t& p) { return get_index(get_entry(p)); }
real_coefficient_t get_coefficient(const real_filtration_entry_t& p) { return get_coefficient(get_entry(p)); }
const value_t& get_filtration(const real_filtration_entry_t& p) { return p.first; }
void set_coefficient(real_filtration_entry_t& p, const real_coefficient_t c) { set_coefficient(get_entry(p), c); }

template <typename Entry> struct real_smaller_filtration_or_smaller_index {
	bool operator()(const Entry& a, const Entry& b) {
		return (get_filtration(a) < get_filtration(b)) ||
		       ((get_filtration(a) == get_filtration(b)) && (get_index(a) < get_index(b)));
	}
};

template <typename ValueType> class real_compressed_sparse_matrix {
	std::deque<size_t> bounds;
	std::deque<ValueType> entries;

public:
	size_t size() const { return bounds.size(); }

	void clear() {
		bounds.clear();
		bounds.shrink_to_fit();
		entries.clear();
		entries.shrink_to_fit();
	}

	typename std::deque<ValueType>::const_iterator cbegin(size_t index) const {
		assert(index < size());
		return index == 0 ? entries.cbegin() : entries.cbegin() + bounds[index - 1];
	}

	typename std::deque<ValueType>::const_iterator cend(size_t index) const {
		assert(index < size());
		return entries.cbegin() + bounds[index];
	}

	template <typename Iterator> void append_column(Iterator begin, Iterator end) {
		for (Iterator it = begin; it != end; ++it) { entries.push_back(*it); }
		bounds.push_back(entries.size());
	}

	void append_column() { bounds.push_back(entries.size()); }

	void push_back(ValueType e) {
		assert(0 < size());
		entries.push_back(e);
		++bounds.back();
	}

	void pop_back() {
		assert(0 < size());
		entries.pop_back();
		--bounds.back();
	}

	template <typename Collection> void append_column(const Collection collection) {
		append_column(collection.cbegin(), collection.cend());
	}

    void print() {

		std::cout << "printing compressed_sparse_matrix" << std::endl;
		for(real_filtration_entry_t entry : entries)
			std::cout << "entry="<< get_entry(entry) << " index=" <<get_index(entry) << " coeff="<< get_coefficient(entry) << std::endl;
		for(size_t bound_a : bounds)
			std::cout << "bounds=" << bound_a << std::endl;
		
		std::cout << "end compressed sparse matrix" << std::endl;
		
	}
	void print_uncompressed(const int rows) {
		//compressed sparse matrix doesn't hold information about the actual size of the matrix
		int this_size = size();
		std::cout << "printing compressed_sparse_matrix as an uncompressed matrix in MATLAB format" << std::endl;
		int cols = bounds.size();
		std::vector< std::vector<real_coefficient_t> > uncompressed_matrix(rows,std::vector<real_coefficient_t>(cols,0.0));
		
		
		for(auto j = 0; j < this_size; j++){
			for(auto it = cbegin(j); it != cend(j); ++it){
				std::cout << *it << "..." << std::flush;
				uncompressed_matrix[get_index(*it)][j] = get_coefficient(*it);
			}
		}
		std::cout<< "\n[";
		for(int i = 0; i < rows; i++){
			for(int j = 0; j < this_size; j++){
				std::cout << uncompressed_matrix[i][j] << " "; 
			}
			std::cout << ";" << std::endl << std::flush;
		}
		std::cout << "]" << std::endl << std::flush;
		

	}

	SparseMatrix to_Eigen(int rows, int columns,int dim){
		SparseMatrix sparse(rows, columns);
		int this_size = size();
		sparse.reserve(Eigen::VectorXi::Constant(columns, dim + 1));//columns or rows?
		for(auto j = 0; j < this_size; j++){
			for(auto it = cbegin(j); it != cend(j); ++it){
				sparse.insert(get_index(*it),j) = get_coefficient(*it); 
			}
		}
		return sparse;
	}
};


template <typename Complex> class real_persistence_computer_t {
private:
	Complex& complex;
	size_t max_entries;
    std::string out_prefix;

    unsigned short min_dimension;
	unsigned short max_dimension;
	unsigned short top_dimension;

    std::vector<std::vector<real_filtration_index_t>> all_filtration_pairs;
	std::vector<value_t> total_filtration;
	std::vector<std::vector<int>> indices_of_filtered_boundaries;		

    std::vector<SparseMatrix> coboundaries;
	std::vector<SparseMatrix> sorted_coboundaries;
	std::vector<SparseMatrix> sorted_boundaries;
    std::vector<std::vector<std::vector<real_coefficient_t>>> spectra;

    std::vector<std::vector<int>> persistent_betti;
	std::vector<std::vector<real_coefficient_t>> persistent_lambda;

public:
	real_persistence_computer_t(Complex& _complex, size_t _max_entries = std::numeric_limits<size_t>::max(), 
		                        std::string _out_prefix = "") : complex(_complex), max_entries(_max_entries),
                                out_prefix(_out_prefix){}


    void compute_persistent_spectra(unsigned short min_dim = 0,
									unsigned short max_dim = std::numeric_limits<unsigned short>::max(),
									std::string out_prefix = ""){
		min_dimension = min_dim;
		max_dimension = max_dim;
		compute_coboundaries(min_dimension, max_dimension);
		sort_coboundaries();
		
		make_boundaries();
		
		set_total_filtration(); 
		
		bool save_boundaries = true;
		if (save_boundaries){
			write_filtrations();
			for (int dim = 0; dim < top_dimension; dim++){
				write_eigen_sparse(sorted_boundaries[dim], "./matrix_" + std::to_string(dim));
			}
		}
		setup_eigenval_storage();
        // TODO: convert to std::vector<SparseMatrixInt> boundaries, std::vector<std::vector<filtration_type>> filtrations
        std::vector<Eigen::SparseMatrix<int,Eigen::ColMajor>> boundaries_int;
        for (int i = 0; i < (int) sorted_boundaries.size(); i++){
            boundaries_int.push_back(sorted_boundaries[i].cast<int>());
        }
        std::vector<std::vector<double>> all_filtrations;
        for (int i = 0; i < (int) all_filtration_pairs.size(); i++){
            std::vector<real_filtration_index_t> filtration_pairs_in_dim = all_filtration_pairs[i];
            std::vector<double> filtrations_in_dim;
            for (int j = 0; j < (int) filtration_pairs_in_dim.size(); j++){
                filtrations_in_dim.push_back(get_filtration(filtration_pairs_in_dim[j]));
            }
            all_filtrations.push_back(filtrations_in_dim);
        }

        auto start_pl = std::chrono::high_resolution_clock::now();
        PersistentLaplacians::PersistentLaplacian pl(boundaries_int, all_filtrations);

        std::vector<std::tuple<int, double, double, std::vector<real_coefficient_t>>> computed_results = pl.spectra();
        
        for (int i = 0; i < (int) computed_results.size(); i++){
            int dim = std::get<0>(computed_results[i]);
            std::vector<real_coefficient_t> current_spectra = std::get<3>(computed_results[i]);
            spectra[dim].push_back(current_spectra);

        }
        auto end_pl = std::chrono::high_resolution_clock::now();
        auto duration_pl = std::chrono::duration_cast<std::chrono::milliseconds>(end_pl - start_pl);
        std::cout << "duration pl (ms):" << duration_pl.count() << std::endl;

		compute_summaries();
		store_spectra();	
		store_spectra_summary();
	}

    void write_filtrations(){
		for (int i = 0; i < (int) all_filtration_pairs.size(); i++){

			std::vector<real_filtration_index_t> filtration_pairs = all_filtration_pairs[i];
			std::ofstream outstream("./" + out_prefix + "_filtrations_" + std::to_string(i) + ".txt");			
			for (int j = 0; j < (int) filtration_pairs.size(); j++){
				real_filtration_index_t current_filtration = filtration_pairs[j];
				outstream << get_filtration(current_filtration) <<",";

			}
			outstream << std::endl;
			outstream.close();
		}
	}

    void store_spectra(){
		std::cout << "\nBegin writing all spectra to files...\n" << std::flush;
		for (int i = 0; i < (int) spectra.size(); i++){
			std::cout <<"Writing spectra of dimension " << i << " to file " << out_prefix << "_spectra_" << i << ".txt\n";
			std::ofstream outstream("./" + out_prefix + "_spectra_" + std::to_string(i) + ".txt");
			
			std::vector<std::vector<real_coefficient_t>> current_dim = spectra[i];
			for (int j = 0; j < (int) current_dim.size(); j++){
				std::vector<real_coefficient_t> current_filtration = current_dim[j];
				for (int k = 0; k < (int) current_filtration.size(); k++){
					outstream << current_filtration[k] << " ";
				}
				outstream << std::endl;
			}
			outstream.close();
		}
	}

    void compute_summaries(){
		for (int dim = min_dimension; dim <= top_dimension; dim++){
			persistent_lambda.push_back(std::vector<real_coefficient_t>());
			persistent_betti.push_back(std::vector<int>());
		}
		double thresshold = pow(0.1,8);
		for (int i = 0; i < (int) spectra.size(); i++){
			std::vector<std::vector<real_coefficient_t>> current_dim = spectra[i];

			for (int j = 0; j < (int) current_dim.size(); j++){
				// NB: stop loop one before last filtration!
				// 	the last filtration step is a dummy filtration
				//	for conveniently computing the final state of the complex
				std::vector<real_coefficient_t> current_filtration = current_dim[j];
				int current_betti = 0;
				bool added = false;
				for (int k = 0; k < (int) current_filtration.size(); k++){
					if (current_filtration[k] > thresshold){
						persistent_lambda[i].push_back(current_filtration[k]);
						persistent_betti[i].push_back(current_betti);
						added = true;
						break;
					}
					current_betti++;
				}
				// if (current_betti == 0){
				// 	persistent_betti[i].push_back(0);
				// 	if (current_filtration.size() > 0){
				// 		persistent_lambda[i].push_back(current_filtration[0]);
				// 	} else{
				// 		persistent_lambda[i].push_back(0);
				// 	}
				if (!added){
					if (current_betti == 0){
						persistent_betti[i].push_back(0);
						if (current_filtration.size() > 0){
							persistent_lambda[i].push_back(current_filtration[0]);
						} else {
							persistent_lambda[i].push_back(0);
						}
					} else { //never encountered a nonzero eigenvalue
						persistent_betti[i].push_back(current_betti);
						persistent_lambda[i].push_back(0);
					}
				

				} //end !added
			} //end j
		} // end i
	}
	void store_spectra_summary(){
		//write all the filtration, betti_k^{i,i+1}, and lambda_k^{i,i+1}
		std::ofstream outstream("./" + out_prefix + "_spectra_summary.txt");
		
		//column headers of tab-separated data
		outstream << "i\tfiltration";
		for (int dim = min_dimension; dim <= top_dimension; dim++){
			outstream << "\tbetti_" << dim;
		}
		for (int dim = min_dimension; dim <= top_dimension; dim++){
			outstream << "\tlambda_" << dim;
		}
		outstream << std::endl;
		//output values
		for (int filtration_index = 0; filtration_index < (int) total_filtration.size()-1; filtration_index++){
			outstream << filtration_index << "\t" << total_filtration[filtration_index];
			for (int dim = min_dimension; dim <= top_dimension; dim++){
				outstream << "\t" << persistent_betti[dim][filtration_index];
			}
			for (int dim = min_dimension; dim <= top_dimension; dim++){
				outstream << "\t" << persistent_lambda[dim][filtration_index];
			}
			outstream << std::endl;
		}
		outstream.close();
	}

    void compute_coboundaries(unsigned short min_dimension, unsigned short max_dimension){
		for (auto dimension = 0u; dimension <= max_dimension; ++dimension) {

			std::vector<real_filtration_index_t> filtration_pairs;
			
			//build the set of real_filtration_index_t to sort
			index_t num_cells = index_t(complex.number_of_cells(dimension));
			for(index_t index = 0; index < num_cells; ++index){
				value_t filtration = complex.filtration(dimension, index);
				filtration_pairs.push_back(std::make_pair(filtration, index));
			}

			std_algorithms::sort(filtration_pairs.begin(), filtration_pairs.end(),
		                     real_smaller_filtration_or_smaller_index<real_filtration_index_t>());

			all_filtration_pairs.push_back(filtration_pairs);

			complex.prepare_next_dimension(dimension);//DEBUG TODO: is this in the right place
			coboundaries.push_back(complex.get_coboundary_as_Eigen());
			if (complex.is_top_dimension() || dimension >= max_dimension) {
				top_dimension = dimension;
				// output->remaining_homology_is_trivial();
				break;
			}
		}
	}

	void sort_coboundaries() {
		// set sorted_coboundaries = coboundary matrix ordered by filtration value then index
		for(index_t dimension = min_dimension; dimension < top_dimension; dimension++){
			if (complex.number_of_cells(dimension+1) == 0)
				break;
			SparseMatrix current = coboundaries[dimension];
			SparseMatrix temp = reorder_columns(current, all_filtration_pairs[dimension]);
			temp = reorder_rows(temp, all_filtration_pairs[dimension+1]);
			sorted_coboundaries.push_back(temp);
		}

	}

	SparseMatrix reorder_columns(SparseMatrix a,std::vector<real_filtration_index_t> order){
		int num_cols = a.cols();
		int num_rows = a.rows();
		SparseMatrix b(num_rows, num_cols);
		for (int i = 0; i < num_cols; i++){
			b.col(i) = a.col(get_index(order[i]));
		}
		return b;
	}

	SparseMatrix reorder_rows(SparseMatrix a, std::vector<real_filtration_index_t> order){
		//https://stackoverflow.com/questions/57858014/permute-columns-of-matrix-in-eigen
		//Since the SparseMatrix type is column-major, the variable a.row(i) is read-only. We must use
		//  permutation matrices to reorder the rows of the column-major matrix.
		Eigen::VectorXi indices(a.rows());
		for (int i = 0; i < (int) indices.size(); i++){
			/* Eigen documentation (https://eigen.tuxfamily.org/dox/classEigen_1_1PermutationMatrix.html) as of 2023-05-08
					"The indices array has the meaning that the permutations sends each integer i to indices[i]."
					This means it sends the ith row to the indices[i]-th row
					i.e. b.row(indices[i]) = a.row(i)
				We have the vector named order sorted to mean that we want to send order[j]th row of a  to the jth row of b
					i.e. b.row(j) = a.row(order[j])
				Then "indices[order[i]] = i" will give us a Permutation matrix that for each k:
					send kth row of a to the to indices[k]th row of b , 
						here k=order[i], and indices[order[i]] = i, so we send the order[i]th row of a to the ith row of b  
						This is exactly the meaning we describe above for what we intend the vector order to mean.
				A more readable algorithm would be:
						indices[i] = order.index_of(i)
					but index_of implementations are O(n) and we need to loop over all elements, 
					which would take O(n^2). Here n=number of simplices in current dimension. 
			*/
			indices[get_index(order[i])] = i;
		}

		Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P;
		P.indices() = indices;
		SparseMatrix b = P*a;
		return b;

	}

	SparseMatrix boundary_at_filtration(int dim, int filtration_index){
		//get B_{dim}^{filtration_index} (where filtration_index is the filtration_index-th filtration step, not the filtration value )
		int next_row_index = indices_of_filtered_boundaries[dim-1][filtration_index];
		int next_col_index = indices_of_filtered_boundaries[dim][filtration_index];
		
		return sorted_boundaries[dim-1].block(0,0,next_row_index+1, next_col_index+1);
	}

    void set_indices_of_filtered_boundaries(){

		std::vector<std::vector<int>> temp_indices_of_filtered_boundaries(top_dimension+1,std::vector<int>(total_filtration.size()));
		for (auto dimension = 0u; dimension <= top_dimension; ++dimension) {
			std::vector<real_filtration_index_t> filtration_pairs = all_filtration_pairs[dimension];
			for (int j = 0; j < (int) total_filtration.size(); j++){
				int k = 0;
				while (k < (int) filtration_pairs.size() && get_filtration(filtration_pairs[k]) <= total_filtration[j]){//this is should be re-written for clarity.
					k++;
				}
				temp_indices_of_filtered_boundaries[dimension][j] = k-1;
			}
		}
		indices_of_filtered_boundaries = temp_indices_of_filtered_boundaries;
	}
	void set_total_filtration(){
		std::set<real_coefficient_t> total_filtration_set;
		for (int i = 0; i < (int) all_filtration_pairs.size(); i++){
			std::vector<real_filtration_index_t> current = all_filtration_pairs[i];
			for (int j = 0; j < (int) current.size(); j++){
				total_filtration_set.insert(get_filtration(current[j]));		
			}
		}		
		for (auto& it: total_filtration_set ){		
			total_filtration.push_back(it);
		}		
		// we also want to capture the final state, 
		// e.g. if the filtration values were 0 1 2
		// and the persistence intervals in dimension 0 would be [0, ) and [1, 2) and [2, ), 
		// we need to capture things that are born in 2 and don't die by 2.
		// The simple way to program this is to add a dummy filtration value of 3, because all cells will be in 
		// both C_2 and C_3, so we will have L_n^{2,3} = L_n^2.
		// This is not the most efficient way to achieve this.
		// total_filtration is already in sorted order, so we only need to append a value greater than the last 
		int count_filtrations = total_filtration.size();
		real_coefficient_t last_filtration = total_filtration[count_filtrations-1];		
        total_filtration.push_back(last_filtration+1.0);

		set_indices_of_filtered_boundaries();		
	}

	void make_boundaries() {
		for (int i = 0; i < (int) sorted_coboundaries.size(); i++){
			sorted_boundaries.push_back(sorted_coboundaries[i].transpose());
		}
	}

    void setup_eigenval_storage(){
		for (int i = 0; i <= top_dimension; i++){
			spectra.push_back(std::vector<std::vector<real_coefficient_t>>());
		}
	}

	void write_eigen_sparse(SparseMatrix m, std::string filename){
        Eigen::saveMarket(m,filename);
    }

    SparseMatrix read_eigen_sparse(std::string filename){
        SparseMatrix m;
        Eigen::loadMarket(m, filename);
        return m;
    }



};
