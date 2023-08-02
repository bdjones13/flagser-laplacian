#pragma once

// #define INDICATE_PROGRESS
// #define USE_COEFFICIENTS
// #define USE_GOOGLE_HASHMAP


#include <cassert>
#include <deque>
#include <iostream>
#include <queue>
#include <set>

#include "definitions.h"
#include "output/base.h"

#include "Eigen/Eigen/src/Core/IO.h"
// #include "MatlabEngine.hpp"
// #include "MatlabDataArray.hpp"
// // #include "MatlabCppSharedLib.hpp"
// #include "Eigen/Eigen/Sparse"
// #include "Eigen/Eigen/Dense"
// typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SparseMatrix;

#ifdef USE_ARRAY_HASHMAP
typedef std::deque<index_t> pivot_column_index_t;
// const index_t INVALID_INDEX = std::numeric_limits<index_t>::max();
#else
typedef fast_hash_map<index_t, index_t> pivot_column_index_t;
#endif



// std::vector<coefficient_t> multiplicative_inverse_vector(const coefficient_t m) {
// 	std::vector<coefficient_t> inverse(m);
// 	inverse[1] = 1;
// 	// m = a * (m / a) + m % a
// 	// Multipying with inverse(a) * inverse(m % a):
// 	// 0 = inverse(m % a) * (m / a)  + inverse(a)  (mod m)
// 	for (coefficient_t a = 2; a < m; ++a) inverse[a] = m - (inverse[m % a] * (m / a)) % m;
// 	return inverse;
// }

// TODO: Make the packed attribute work again
struct real_entry_t {
	// index_t index : 8 * (sizeof(index_t) - sizeof(coefficient_t));
	index_t index;
	real_coefficient_t coefficient;
	real_entry_t(index_t _index, real_coefficient_t _coefficient) : index(_index), coefficient(_coefficient) {}
	real_entry_t(index_t _index) : index(_index), coefficient(1) {}
	real_entry_t() : index(0), coefficient(1) {}
};
// } __attribute__((packed));

// static_assert(sizeof(real_entry_t) == sizeof(index_t), "size of real_entry_t is not the same as index_t");

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

//commented Ben Jones 2023-03-13, can probably be removed
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

//commented Ben Jones 2023-03-13, can probably be removed
template <typename Entry> struct real_greater_filtration_or_smaller_index {
	bool operator()(const Entry& a, const Entry& b) {
		return (get_filtration(a) > get_filtration(b)) ||
		       ((get_filtration(a) == get_filtration(b)) && (get_index(a) < get_index(b)));
	}
};

template <typename Entry> struct real_smaller_filtration_or_smaller_index {
	bool operator()(const Entry& a, const Entry& b) {
		return (get_filtration(a) < get_filtration(b)) ||
		       ((get_filtration(a) == get_filtration(b)) && (get_index(a) < get_index(b)));
	}
};


class real_filtered_union_find {
	std::vector<index_t> parent;
	std::vector<uint8_t> rank;
	const std::vector<value_t> filtration;

public:
	real_filtered_union_find(const std::vector<value_t>& _filtration)
	    : parent(_filtration.size()), rank(_filtration.size(), 0), filtration(_filtration) {
		for (size_t i = 0ul; i < _filtration.size(); ++i) parent[i] = index_t(i);
	}

	index_t find(index_t x) {
		index_t y = x, z = parent[y];
		while (z != y) {
			y = z;
			z = parent[y];
		}
		y = parent[x];
		while (z != y) {
			parent[x] = z;
			x = y;
			y = parent[x];
		}
		return z;
	}
	void link(index_t x, index_t y) {
		x = find(x);
		y = find(y);
		if (x == y) return;
		if (filtration[x] < filtration[y] || (filtration[x] == filtration[y] && rank[x] > rank[y]))
			parent[y] = x;
		else {
			parent[x] = y;
			if (rank[x] == rank[y]) ++rank[y];
		}
	}
};

template <typename Heap>
real_filtration_entry_t pop_pivot(Heap& column
#ifdef USE_COEFFICIENTS
                             ,
                             coefficient_t modulus
#endif
) {
	if (column.empty())
		return real_filtration_entry_t(-1);
	else {
		auto pivot = column.top();

#ifdef USE_COEFFICIENTS
		coefficient_t coefficient = 0;
		do {
			coefficient = (coefficient + get_coefficient(column.top()));// Ben Jones 2023-03-03// % modulus;
			column.pop();

			if (coefficient == 0) {
				if (column.empty())
					return real_filtration_entry_t(-1);
				else
					pivot = column.top();
			}
		} while (!column.empty() && get_index(column.top()) == get_index(pivot));
		if (get_index(pivot) != -1) { set_coefficient(pivot, coefficient); }
#else
		column.pop();
		while (!column.empty() && get_index(column.top()) == get_index(pivot)) {
			column.pop();
			if (column.empty())
				return real_filtration_entry_t(-1);
			else {
				pivot = column.top();
				column.pop();
			}
		}
#endif
		return pivot;
	}
}

template <typename Heap>
real_filtration_entry_t get_pivot(Heap& column
#ifdef USE_COEFFICIENTS
                             ,
                             coefficient_t modulus
#endif
) {
	real_filtration_entry_t result = pop_pivot(column
#ifdef USE_COEFFICIENTS
	                                      ,
	                                      modulus
#endif
	);
	if (get_index(result) != -1) column.push(result);
	return result;
}

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
		std::vector< std::vector<double> > uncompressed_matrix(rows,std::vector<double>(cols,0.0));
		
		
		for(auto j = 0; j < this_size; j++){
			for(auto it = cbegin(j); it != cend(j); ++it){
				std::cout << *it << "..." << std::flush;
				uncompressed_matrix[get_index(*it)][j] = get_coefficient(*it); // 1
				// uncompressed_matrix[it][j] 
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
		
		
		// for(int i = 0; i < rows; i++){
		// 	for (int j = 0; j < cols; j++){
		// 		uncompressed_matrix[i][j] = ;
		// 	}
		// }
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

template <typename Heap> void push_entry(Heap& column, index_t i, real_coefficient_t c, value_t filtration) {
	real_entry_t e = real_make_entry(i, c);
	column.push(std::make_pair(filtration, e));
}

// This class is just an ordinary priority queue, but once the
// queue gets too long (because a lot of faces are inserted multiple
// times) it starts collecting the coefficients and only inserting each
// new face once
template <class Container, class Comparator>
class real_priority_queue_t : public std::priority_queue<real_filtration_entry_t, Container, Comparator> {
	std::unordered_map<index_t, coefficient_t> coefficients;

	static const real_filtration_entry_t dummy;
// 	// TODO: Enable dynamic switching to the dense version for coefficients
// #ifdef USE_COEFFICIENTS
	bool use_dense_version = true;
// #else
// 	bool use_dense_version = false;
// #endif
	// coefficient_t modulus;
	size_t dense_threshold;

public:
	real_priority_queue_t(size_t _dense_threshold)
	    : dense_threshold(_dense_threshold) {}

	void push(const real_filtration_entry_t& value) {
		if (use_dense_version) {
			// If we already have this value: update the count and don't push it again
			auto p = coefficients.find(get_index(value));
			if (p != coefficients.end()) {
// #ifdef USE_COEFFICIENTS
				p->second = (p->second + get_coefficient(value)) ;//% modulus;
// #else
// 				p->second = !p->second;
// #endif
				return;
			}
		}

		std::priority_queue<real_filtration_entry_t, Container, Comparator>::push(value);

		if (use_dense_version) coefficients.insert(std::make_pair(get_index(value), get_coefficient(value)));

#ifndef USE_COEFFICIENTS
		if (!use_dense_version &&
		    std::priority_queue<real_filtration_entry_t, Container, Comparator>::size() >= dense_threshold)
			use_dense_version = true;
#endif
	}

	void pop() {
		// Don't use this, only allow get_pivot
		throw std::exception();
	}

	real_filtration_entry_t pop_pivot() {
		remove_trivial_coefficient_entries();
		if (std::priority_queue<real_filtration_entry_t, Container, Comparator>::empty())
			return dummy;
		else {
			auto pivot = get_top();

#ifdef USE_COEFFICIENTS
			coefficient_t coefficient = 0;
			do {
				coefficient = (coefficient + get_coefficient(get_top()));// % modulus;
				safe_pop();
				remove_trivial_coefficient_entries();

				if (coefficient == 0) {
					if (std::priority_queue<real_filtration_entry_t, Container, Comparator>::empty())
						return dummy;
					else
						pivot = get_top();
				}
			} while (!std::priority_queue<real_filtration_entry_t, Container, Comparator>::empty() &&
			         get_index(get_top()) == get_index(pivot));
			if (get_index(pivot) != -1) { set_coefficient(pivot, coefficient); }
#else
			safe_pop();
			while (!std::priority_queue<real_filtration_entry_t, Container, Comparator>::empty() &&
			       get_index(std::priority_queue<real_filtration_entry_t, Container, Comparator>::top()) ==
			           get_index(pivot)) {
				safe_pop();
				remove_trivial_coefficient_entries();

				if (std::priority_queue<real_filtration_entry_t, Container, Comparator>::empty())
					return dummy;
				else {
					pivot = get_top();
					safe_pop();
				}
			}
#endif
			return pivot;
		}
	}

	real_filtration_entry_t get_pivot() {
		real_filtration_entry_t result = pop_pivot();
		if (get_index(result) != -1) { push(result); }
		return result;
	}

private:
	inline real_filtration_entry_t get_top() {
		auto pivot = std::priority_queue<real_filtration_entry_t, Container, Comparator>::top();

#ifdef USE_COEFFICIENTS
		if (use_dense_version) {
			auto e = coefficients.find(get_index(pivot));
			if (e != coefficients.end()) set_coefficient(pivot, e->second);
		}
#endif

		return pivot;
	}

	inline void safe_pop() {
		if (use_dense_version) {
			auto e =
			    coefficients.find(get_index(std::priority_queue<real_filtration_entry_t, Container, Comparator>::top()));
			if (e != coefficients.end()) coefficients.erase(e);
		}
		std::priority_queue<real_filtration_entry_t, Container, Comparator>::pop();
	}

	inline void remove_trivial_coefficient_entries() {
		if (use_dense_version) {
			if (std::priority_queue<real_filtration_entry_t, Container, Comparator>::size() == 0) return;
			auto p =
			    coefficients.find(get_index(std::priority_queue<real_filtration_entry_t, Container, Comparator>::top()));
#ifdef USE_COEFFICIENTS
			while (p != coefficients.end() ){//&& p->second % modulus == 0) {
#else
			while (p != coefficients.end() && p->second == false) {
#endif
				coefficients.erase(p);
				std::priority_queue<real_filtration_entry_t, Container, Comparator>::pop();
				p = coefficients.find(get_index(std::priority_queue<real_filtration_entry_t, Container, Comparator>::top()));
			}
		}
	}
};
template <class Container, class Comparator>
const real_filtration_entry_t
    real_priority_queue_t<Container, Comparator>::dummy(real_filtration_entry_t(std::make_pair(value_t(0.0f), real_entry_t(-1))));

#ifdef SORT_COLUMNS_BY_PIVOT
template <typename Entry, typename Complex> struct real_greater_filtration_or_better_pivot_or_smaller_index {
	real_greater_filtration_or_better_pivot_or_smaller_index(Complex& _complex) : complex(_complex) {}
	bool operator()(Entry a, Entry b) const {
		// First order by the filtration value
		if (get_filtration(a) > get_filtration(b)) return true;
		if (get_filtration(a) < get_filtration(b)) return false;

		const auto& ta = get_coboundary_size_and_gap_and_pivot(a);
		const auto& tb = get_coboundary_size_and_gap_and_pivot(b);

		// Then the number of non-trivial coboundary entries
		if (std::get<0>(ta) < std::get<0>(tb)) return true;
		if (std::get<0>(ta) > std::get<0>(tb)) return false;

		// Then order by the better pivoting
		if (std::get<2>(ta) < std::get<2>(tb)) return true;
		if (std::get<2>(ta) > std::get<2>(tb)) return false;

		if (std::get<1>(ta) > std::get<1>(tb)) return true;
		if (std::get<1>(ta) < std::get<1>(tb)) return false;

		// Finally, order by their indices
		return get_index(a) < get_index(b);
	}

private:
	Complex& complex;

	// A column is considered to be a better pivot if the jump from pivot to the next
	// non-trivial element is as big as possible. This prevents accidentally inserting
	// non-trivial elements just below the pivot, which sometimes creates very long
	// reduction chains.
	// The second sort criterium is for it to be small because the small pivots will be
	// used the most.
	std::tuple<size_t, size_t, index_t> get_coboundary_size_and_gap_and_pivot(real_filtration_entry_t a) const {
		// Look at the first two gaps of the pivot and the next element
		index_t pivot = 0;
		size_t gap_after_pivot = 0;
		auto iterator = complex.coboundary(a);
		size_t coboundary_size = 0;
		while (iterator.has_next()) {
			coboundary_size++;
			index_t next_index = get_index(iterator.next().second);
			if (next_index > pivot) {
				gap_after_pivot = next_index - pivot;
				pivot = next_index;
			}
		}

		return std::make_tuple(coboundary_size, gap_after_pivot, pivot);
	}
};
#endif

template <typename Complex> class real_persistence_computer_t {
private:
	Complex& complex;
	output_t<Complex>* output;
	value_t max_filtration;
	size_t max_entries;
	index_t euler_characteristic = 0;
	bool print_betti_numbers_to_console = true;


	unsigned short min_dimension;
	unsigned short max_dimension;
	unsigned short top_dimension;
	// coefficient_t modulus = 0;//2;
	std::vector<coefficient_t> multiplicative_inverse;
	std::deque<real_filtration_index_t> columns_to_reduce;

	std::vector<std::vector<real_filtration_index_t>> all_filtration_pairs;
	std::vector<value_t> total_filtration;
	std::vector<std::vector<int>> indices_of_filtered_boundaries;			
	
	std::vector<SparseMatrix> coboundaries;
	std::vector<SparseMatrix> sorted_coboundaries;
	std::vector<SparseMatrix> sorted_boundaries;
	std::vector<SparseMatrix> Laplacians;

	SparseMatrix B_a;
	SparseMatrix B_ab;

	// std::unique_ptr<matlab::engine::MATLABEngine> m_matlab_engine; //TODO: re-enable matlab
	// std::vector<std::pair<std::vector<real_compressed_sparse_matrix<real_entry_t>>, std::vector<size_t> >> coboundaries;

public:
	real_persistence_computer_t(Complex& _complex, output_t<Complex>* _output,
	                       size_t _max_entries = std::numeric_limits<size_t>::max(), //int _modulus = 2,
	                       value_t _max_filtration = std::numeric_limits<value_t>::max())
	    : complex(_complex), output(_output), max_filtration(_max_filtration), max_entries(_max_entries) {
			//m_matlab_engine = matlab::engine::startMATLAB();//TODO: re-enable matlab
		}

	    //   multiplicative_inverse(multiplicative_inverse_vector(modulus)) {

	

	void set_print_betti_numbers(bool print_betti_numbers) { print_betti_numbers_to_console = print_betti_numbers; }

	void compute_persistent_spectra(unsigned short min_dim = 0,
									unsigned short max_dim = std::numeric_limits<unsigned short>::max()){
		min_dimension = min_dim;
		max_dimension = max_dim;

		compute_coboundaries(min_dimension, max_dimension);
		// std::cout << "\n\% coboundaries ... \n";
		// for(int i = 0; i < (int) coboundaries.size(); i++){
		// 	print_Eigen_Sparse(coboundaries[i]);
		// }
		// std::cout << "\n\% sorting coboundaries ... \n";
		sort_coboundaries();
		// std::cout << "\n\% sorted coboundaries ... \n";
		// for(int i = 0; i < (int) sorted_coboundaries.size(); i++){
		// 	std::cout <<"\ncbdy" << i << "=";
		// 	print_Eigen_Sparse(sorted_coboundaries[i]);
		// }
		// std::cout << "\n\% end sorted coboundaries ... \n";
		make_boundaries();
		set_total_filtration();
		print_all_filtration_pairs();

		std::cout << "\n%matlab list for spectra\n" << "betti0=[];\n"<< "betti1=[];\n"<< "betti2=[]; betti3=[]\n";  
		for (int i = 0; i < (int) total_filtration.size()-1; i++){
		
			double filtration = total_filtration[i];
			double next_filtration = total_filtration[i+1];
			
			for(int dim = min_dim; dim <= top_dimension; dim++){
				std::cout << "\n\% B_" << filtration << "[" << dim-1 << "] with cols in dim=" << dim << ", rows in dim-1=" << dim-1;
				std::cout << "\nB_a=" << std::flush;
				if (dim != min_dim){
					boundary_at_filtration(dim,filtration,i);

					print_Eigen_Sparse(B_a);
				} 
				else {

					std::cout << "zeros(0," << indices_of_filtered_boundaries[dim][i]+1 << ");";
				}			

				if (dim != top_dimension){	
					boundary_between_filtrations(dim, filtration, next_filtration,i); 
				}
				else{
					
					std::cout << "\nB_qp1_L=zeros(" << indices_of_filtered_boundaries[dim][i] << ",0);\n";
				}
				std::cout << "\n\% dim = " << dim << "filtration=" << filtration << ", next_filtration=" << next_filtration;
				std::cout << "\n evals=PL(B_qp1_L,n_qL,n_qK,B_a);\n" << std::flush;
				std::cout << "betti" << dim << "=[betti" << dim << "; nnz(~evals)];\n"; //matlab code to build a list of betti numbers, like "betti1 = [betti1; 2];"
			}
		}
		std::cout << "\n\% print matlab lists for spectra\n" << "betti0\n"<< "betti1\n"<< "betti2\n betti3\n";	
	}

	void compute_persistence(unsigned short min_dimension = 0,
	                         unsigned short max_dimension = std::numeric_limits<unsigned short>::max(),
	                         bool check_euler_characteristic = true) {
		compute_zeroth_persistence(min_dimension, max_dimension);// TODO: check this
		compute_higher_persistence(min_dimension, max_dimension);// TODO: check this
		compute_laplacians();// TODO: min dimension and max dimension
		for(int i = 0; i < (int) Laplacians.size(); i++){
			print_Eigen_Sparse(Laplacians[i]);
		}

		std::vector<std::vector<double>> spectra;
		for (int dim = min_dimension; dim < max_dimension-1; dim++){
			std::cout << "compute spectra at dimension=" << dim << "\n" << std::flush;
            spectra.push_back(compute_spectra(dim,0));
			std::cout << "finish compute spectra at dimension=" << dim << "\n" << std::flush;
            if (spectra[dim].size() == 0)
                break;
		}

		for (int dim = min_dimension; dim < max_dimension; dim++){
					std::cout << "spectra[" << dim << "]=[";
			std::vector<double> cur_spectra = spectra[dim];
			
			for(int i = 0; i < (int) cur_spectra.size(); i++){
				std::cout << cur_spectra[i] << ", ";
			}
			std::cout << "]\n" << std::flush;
            if (cur_spectra.size() == 0)
                break;
		}

	
		complex.finished();
		output->finished(check_euler_characteristic);

		// Sanity check whether there were any problems computing the homology
		bool computed_full_homology = min_dimension == 0 && max_dimension == std::numeric_limits<unsigned short>::max();
		if (check_euler_characteristic && computed_full_homology && max_entries == std::numeric_limits<size_t>::max()) {
			index_t cell_euler_characteristic = 0;
			for (size_t i = 0; i <= complex.top_dimension(); i++) {
				cell_euler_characteristic += (i % 2 == 1 ? -1 : 1) * index_t(complex.number_of_cells(index_t(i)));
			}

			if (print_betti_numbers_to_console)
				std::cout << "The Euler characteristic is given by: " << cell_euler_characteristic << std::endl;

			if (cell_euler_characteristic != euler_characteristic) {
				std::cerr << "ERROR: The homological Euler characteristic (which is " << euler_characteristic
				          << ") differs from the cellular Euler characteristic (which is " << cell_euler_characteristic
				          << "), apparently there is an error in the program." << std::endl;
			}
		}


	}

	void print_Eigen_Sparse(SparseMatrix m){
		Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "];");
		std::cout << Eigen::MatrixXd(m).format(OctaveFmt)  << std::endl;
	}
	void compute_laplacians(){
		std::cout <<"begin computing laplacians" << std::endl;
		//first and last laplacian are different
		//note we start with coboundaries, not boundaries
		Laplacians.push_back(compute_zeroth_laplacian());
		int num_laplacians = (int) coboundaries.size(); 
		for(int i = 0; i < num_laplacians-1; i++){
			Laplacians.push_back(compute_middle_laplacian(i));
		}
		Laplacians.push_back(compute_last_laplacian());
		
		std::cout <<"end computing laplacians" << std::endl;
	}

	SparseMatrix compute_zeroth_laplacian(){
		//note we start with coboundaries, not boundaries
		SparseMatrix d1 = coboundaries[0];
		return d1.transpose()*d1;
	}

	SparseMatrix compute_middle_laplacian(int i){
		//note we start with coboundaries, not boundaries
		SparseMatrix di = coboundaries[i];
		SparseMatrix dnext = coboundaries[i+1];
		return di*di.transpose() + dnext.transpose()*dnext;
	}

	SparseMatrix compute_last_laplacian(){
		//note we start with coboundaries, not boundaries
		SparseMatrix dn = coboundaries.back();
		return dn*dn.transpose();
	}

	void print_all_filtration_pairs(){
		std::cout << "\n\% all_filtration_pairs={" << std::flush;
		for (int i = 0; i < (int) all_filtration_pairs.size(); i++){
			std::vector<real_filtration_index_t> filtration_pairs = all_filtration_pairs[i];
			std::cout << "\n\%filtration_pairs_" << i << "=[" << std::flush;
			for (int j = 0; j < (int) filtration_pairs.size(); j++){
				std::cout <<"(" << get_filtration(filtration_pairs[j])<<"," << get_index(filtration_pairs[j]) << "),";
			}
			std::cout << "]" << std::flush;
		}
		std::cout << "}(end all_filtration_pairs)" << std::flush << std::endl;
	}

	std::vector<double> compute_spectra(int dim, int num_eigenvals){

		//most of this code is identical to https://github.com/wangru25/HERMES/blob/main/src/snapshot.cpp
		// it needed eval=eigenvalue return variable, evec=eigenvector return variable, int matrix_size, and vector<double>& for row, col, and val, and int es 
		// vector for row, col, and val is a sparse matrix representation probably, 
		// matrix_size is used for setting up "vms", es for "ves" "vee matrix size" and "vee es", whatever es is
		// es is the number of eigenvalues
		int matrix_size = complex.number_of_cells(dim);
        std::cout << "matrix size = " << matrix_size << std::endl;
		if (matrix_size <= 0){
            std::vector<double> dummy;
            return dummy;
        }
        std::vector<double> row, col, val;
		// row.reserve(size); //TODO: similarly for col and val - efficiency improvements
		SparseMatrix L = Laplacians[dim];

		for(int i=0; i < L.outerSize(); ++i){
			for(SparseMatrix::InnerIterator iter(L,i); iter; ++iter){
				row.push_back(static_cast<double>(iter.row()+1));
				col.push_back(static_cast<double>(iter.col()+1));
				val.push_back(static_cast<double>(iter.value()));
			}
		}
		std::vector<double> eigenvalues;
	    std::vector<ColumnVector> eigenvectors;
		// int es = num_eigenvals; //TODO: reconcile the variables
		
		if (num_eigenvals == 0)
			//compute all
			num_eigenvals = matrix_size;
		

		//TODO: re-enable matlab
		// std::vector<double> vms = {static_cast<double>(matrix_size)};
		// std::vector<double> ves = {static_cast<double>(num_eigenvals)};
		// using namespace matlab;

		// data::ArrayFactory factory;
		
		// data::TypedArray<double> mms = factory.createArray<std::vector<double>::iterator>({ 1, 1 }, vms.begin(), vms.end());
		// data::TypedArray<double> mrow = factory.createArray<std::vector<double>::iterator>({ row.size(), 1 }, row.begin(), row.end());
		// data::TypedArray<double> mcol = factory.createArray<std::vector<double>::iterator>({ col.size(), 1 }, col.begin(), col.end());
		// data::TypedArray<double> mval = factory.createArray<std::vector<double>::iterator>({ val.size(), 1 }, val.begin(), val.end());
		// data::TypedArray<double> m_num_eigenvals = factory.createArray<std::vector<double>::iterator>({ 1, 1 }, ves.begin(), ves.end());
		
		// m_matlab_engine->setVariable(u"size", std::move(mms));
		// m_matlab_engine->setVariable(u"row", std::move(mrow));
		// m_matlab_engine->setVariable(u"col", std::move(mcol));
		// m_matlab_engine->setVariable(u"val", std::move(mval));
		// m_matlab_engine->setVariable(u"es", std::move(m_num_eigenvals));
		// m_matlab_engine->eval(u"A=sparse(row, col, val, size, size);");
		// m_matlab_engine->eval(u"A=A+speye(size)*1e-12;");
		// m_matlab_engine->eval(u"if size > es \n num = es; \n else \n num = size; \n end \n");
		// m_matlab_engine->eval(u"[V,D]=eigs((A+A')/2, num, 'smallestabs');");
		// data::TypedArray<double> dd = m_matlab_engine->getVariable(u"D");
		// data::TypedArray<double> vv = m_matlab_engine->getVariable(u"V");

		// eigenvalues.reserve(dd.getDimensions()[0]);
		// for(int i=0; i<(int)dd.getDimensions()[0]; ++i){
		// 	eigenvalues.push_back(dd[i][i]);
		// }


		return eigenvalues;

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
			if (complex.is_top_dimension()) {
				top_dimension = dimension;
				output->remaining_homology_is_trivial();
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

	void boundary_at_filtration(int dim, double filtration, int filtration_index){
		//set B_a to be B_{dim}^{filtration}
		//need to get the max row and max column
		int next_row_index = indices_of_filtered_boundaries[dim-1][filtration_index];
		int next_col_index = indices_of_filtered_boundaries[dim][filtration_index];
		// int next_row_index = index_of_filtered_boundary(dim-1,filtration, filtration_index);
		// int next_col_index = index_of_filtered_boundary(dim, filtration, filtration_index);
		
		B_a = sorted_boundaries[dim-1].block(0,0,next_row_index+1, next_col_index+1);
	}

	void boundary_between_filtrations(int dim, double a, double b, int filtration_index){
		//Set the B_{p}^{a,b} representing \left(i_{p-1}^{a,b}\right)^* \circ d_{p}^b \circ \iota_{p}^{a,b}
		// This is B_p^b with columns removed (zeroed out) if the boundary of that column is not in C_{p-1}^a
		// So we go through all columns that are new in B_p^b and not in B_p^a
		int b_row_index = indices_of_filtered_boundaries[dim][filtration_index+1];
		int b_col_index = indices_of_filtered_boundaries[dim+1][filtration_index+1];
		int a_row_index = indices_of_filtered_boundaries[dim][filtration_index];	
		// // int b_row_index = index_of_filtered_boundary(dim,b, filtration_index+1); 
		// // int b_col_index = index_of_filtered_boundary(dim+1, b, filtration_index+1);
		// int a_row_index = index_of_filtered_boundary(dim,a, filtration_index);

		B_ab = sorted_boundaries[dim].block(0,0,b_row_index+1, b_col_index+1); 
		// the last two arguments of .block are the *size* of the matrix, which would occur at index + 1 (https://eigen.tuxfamily.org/dox/group__TutorialBlockOperations.html)

		std::cout <<"\n n_qL=" << b_row_index+1 << ";\n n_qK=" << a_row_index+1 << ";\n" << std::flush;//the +1 is because of matlab indexing
		std::cout << "\n B_qp1_L=" << std::flush;
		print_Eigen_Sparse(B_ab);
		std::cout << ";\n" << std::flush;
	
	}
protected:

	void set_indices_of_filtered_boundaries(){

		std::vector<std::vector<int>> temp_indices_of_filtered_boundaries(top_dimension+1,std::vector<int>(total_filtration.size()));
		for (auto dimension = 0u; dimension <= top_dimension; ++dimension) {
			std::vector<real_filtration_index_t> filtration_pairs = all_filtration_pairs[dimension];
			for (int j = 0; j < (int) total_filtration.size(); j++){
				int k = 0;
				while (get_filtration(filtration_pairs[k]) <= total_filtration[j] && k < (int) filtration_pairs.size()){//this is should be re-written for clarity.
					k++;
				}
				temp_indices_of_filtered_boundaries[dimension][j] = k-1;
			}
		}
		indices_of_filtered_boundaries = temp_indices_of_filtered_boundaries;
	}
	void set_total_filtration(){
		std::set<double> total_filtration_set;
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
		double last_filtration = total_filtration[count_filtrations-1];
		total_filtration.push_back(last_filtration+1.0);
		

		set_indices_of_filtered_boundaries();
	}

	void make_boundaries() {
		for (int i = 0; i < (int) sorted_coboundaries.size(); i++){
			sorted_boundaries.push_back(sorted_coboundaries[i].transpose());
		}
	}

	void compute_zeroth_persistence(unsigned short min_dimension, unsigned short) {
		complex.prepare_next_dimension(0);
		coboundaries.push_back(complex.get_coboundary_as_Eigen());

	
//To just get boundary matrix don't need persistence, Ben Jones 2023-03-13
// 		// Only compute this if we actually need it
// 		if (min_dimension > 1) return;

// #ifdef INDICATE_PROGRESS
// 		std::cout << "\033[K"
// 		          << "computing persistent homology in dimension 0" << std::flush << "\r";
// #endif

// 		long long betti_number = 0;

// 		size_t n = complex.number_of_cells(0);
// 		real_filtered_union_find dset(complex.vertex_filtration());
// 		std::vector<real_filtration_index_t> edges;
// 		index_t number_of_edges = index_t(complex.number_of_cells(1));
// 		for (index_t index = 0; index < number_of_edges; index++) {
// 			value_t filtration = complex.filtration(1, index);
// 			if (filtration <= max_filtration) edges.push_back(std::make_pair(filtration, index));
// 		}
// 		std_algorithms::sort(edges.rbegin(), edges.rend(), real_greater_filtration_or_smaller_index<real_filtration_index_t>());

// 		// Let the output class know that we are now computing zeroth degree barcodes
// 		output->computing_barcodes_in_dimension(0);

// 		for (auto e : edges) {
// 			const auto vertices = complex.vertices_of_edge(get_index(e));
// 			index_t u = dset.find(vertices.first), v = dset.find(vertices.second);

// 			if (u != v) {
// 				// Only output bars if we are interested in zeroth homology
// 				const auto filtration_u = complex.filtration(0, u);
// 				const auto filtration_v = complex.filtration(0, v);
// 				dset.link(u, v);
// 				if (min_dimension == 0 && get_filtration(e) > std::max(filtration_u, filtration_v)) {
// 					// Check which vertex is merged into which other vertex.
// 					const auto f = dset.find(u) == u ? filtration_v : filtration_u;
// 					output->new_barcode(f, get_filtration(e));

// 				}
// 			} else {
// 				columns_to_reduce.push_back(e);
// 			}
// 		}
// 		std::reverse(columns_to_reduce.begin(), columns_to_reduce.end());

// 		// If we don't care about zeroth homology, then we can stop here
// 		if (min_dimension == 1) return;

// 		for (index_t index = 0; index < index_t(n); ++index) {
// 			if (dset.find(index) == index) {
// 				output->new_infinite_barcode(complex.filtration(0, index));
// 				betti_number++;

// 			}
// 		}

// 		// Report the betti number back to the complex and the output
// 		complex.computation_result(0, betti_number, 0);
// 		output->betti_number(betti_number, 0);

// 		if (print_betti_numbers_to_console) {
// 			std::cout << "\033[K"
// 			          << "The dimensions of the real homology of the full complex are:" << std::endl
// 			          << std::endl
// 			          << "dim H_0 = " << betti_number << std::endl;
// 		}
// 		euler_characteristic += index_t(betti_number);
	}

	void compute_higher_persistence(unsigned short min_dimension, unsigned short max_dimension) {
		// max dimension is 65535 for some reason. Actually that's fine, but it should stop before then
		std::cout << "min dimension = " << min_dimension << " max_dimension = " << max_dimension << std::endl;
		for (auto dimension = 1u; dimension <= max_dimension; ++dimension) {
			complex.prepare_next_dimension(dimension);
			coboundaries.push_back(complex.get_coboundary_as_Eigen());

// 			if (dimension + 1 == min_dimension) {
// 				// Here we need to reduce *all* cells because we did not compute anything of smaller dimension
// 				index_t number_of_cells = index_t(complex.number_of_cells(dimension));
// 				for (index_t index = 0; index < number_of_cells; index++) {
// 					columns_to_reduce.push_back(std::make_pair(complex.filtration(dimension, index), index));
// 				}
// 			}

// 			if (dimension + 1 < min_dimension) continue;

// 			output->computing_barcodes_in_dimension(dimension);

// 			sort_columns();

// #ifdef INDICATE_PROGRESS
// 			std::cout << "\033[K"
// 			          << "computing persistent homology in dimension " << dimension << std::flush << "\r";
// #endif
// #ifdef USE_ARRAY_HASHMAP
// 			pivot_column_index_t pivot_column_index(complex.number_of_cells(dimension + 1), INVALID_INDEX);
// #else
// 			pivot_column_index_t pivot_column_index;
// 			pivot_column_index.reserve(complex.number_of_cells(dimension + 1));
// #endif

// 			auto betti = real_compute_pairs(dimension, pivot_column_index, dimension >= min_dimension);
// 			if (dimension >= min_dimension) {
// 				complex.computation_result(dimension, betti.first, betti.second);

// 				output->betti_number(betti.first, betti.second);
// 				euler_characteristic += (dimension & 1 ? -1 : 1) * betti.first;

// 				if (print_betti_numbers_to_console) {
// 					std::cout << "\033[K"
// 					          << "dim H_" << dimension << " = " << betti.first;
// 					if (betti.second > 0) { std::cout << " (skipped " << betti.second << ")"; }
// 					std::cout << std::endl;
// 				}
// 			} else if (int(dimension) == min_dimension - 1 && print_betti_numbers_to_console &&
// 			           max_entries < std::numeric_limits<size_t>::max()) {
// 				std::cout << "\033[K"
// 				          << "# Skipped columns in dimension " << dimension << ": " << betti.second << std::endl;
// 			}
// 			if (dimension < max_dimension) assemble_columns_to_reduce(dimension, pivot_column_index);

// 			// Stop early
			if (complex.is_top_dimension()) {
				output->remaining_homology_is_trivial();
				break;
			}
		}
	}

	
	void assemble_columns_to_reduce(index_t dimension, pivot_column_index_t& pivot_column_index) {
		index_t num_cells = index_t(complex.number_of_cells(dimension + 1));

		columns_to_reduce.clear();

#ifdef INDICATE_PROGRESS
		std::cout << "\033[K"
		          << "assembling " << num_cells << " columns" << std::flush << "\r";
#endif

		for (index_t index = 0; index < num_cells; ++index) {
			if (
#ifdef USE_ARRAY_HASHMAP
			    pivot_column_index[index] == INVALID_INDEX
#else
			    pivot_column_index.find(index) == pivot_column_index.end()
#endif
			) {
				value_t filtration = complex.filtration(dimension + 1, index);
				if (filtration <= max_filtration) { columns_to_reduce.push_back(std::make_pair(filtration, index)); }
#ifdef INDICATE_PROGRESS
				if ((index + 1) % 100000 == 0)
					std::cout << "\033[K"
					          << "assembled " << columns_to_reduce.size() << " out of " << (index + 1) << "/"
					          << num_cells << " columns" << std::flush << "\r";
#endif
			}
		}
#ifdef INDICATE_PROGRESS
		std::cout << "\033[K";
#endif
	}

	void sort_columns() {
#ifdef INDICATE_PROGRESS
		std::cout << "\033[K"
		          << "sorting " << columns_to_reduce.size() << " columns" << std::flush << "\r";
#endif

#ifdef SORT_COLUMNS_BY_PIVOT
		std_algorithms::sort(
		    columns_to_reduce.begin(), columns_to_reduce.end(),
		    real_greater_filtration_or_better_pivot_or_smaller_index<real_filtration_index_t, decltype(complex)>(complex));
#else
		std_algorithms::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
		                     real_greater_filtration_or_smaller_index<real_filtration_index_t>());
#endif
#ifdef INDICATE_PROGRESS
		std::cout << "\033[K";
#endif
	}


//TODO: reimplement persistence pairs
//Removed temporarily Ben Jones 2023-03-13
	std::pair<index_t, index_t> real_compute_pairs(index_t, pivot_column_index_t& pivot_column_index,
	                                          bool generate_output = true) {
		return std::make_pair(0,0); //Ben Jones, 2023-03-13 to skip this function
// 		index_t betti = 0;
// 		index_t betti_error = 0;
// #ifdef INDICATE_PROGRESS
// 		auto verbose_logging_threshold = size_t(columns_to_reduce.size() * 0.90);
// #endif

// #ifdef ASSEMBLE_REDUCTION_MATRIX
// 		compressed_sparse_matrix<real_filtration_entry_t> reduction_coefficients;
// #else
// #ifdef USE_COEFFICIENTS
// 		std::vector<real_filtration_entry_t> reduction_coefficients;
// #endif
// #endif


// 		std::vector<real_filtration_entry_t> coface_entries;

// 		for (auto i = 0ul; i < columns_to_reduce.size(); ++i) {
// 			auto column_to_reduce = columns_to_reduce[i];

// #ifdef ASSEMBLE_REDUCTION_MATRIX
// 			std::priority_queue<real_filtration_entry_t, std::deque<real_filtration_entry_t>, smaller_index<real_filtration_entry_t>>
// 			    reduction_column;
// #endif

// 			real_priority_queue_t<std::deque<real_filtration_entry_t>, real_greater_filtration_or_smaller_index<real_filtration_entry_t>>
// 			    working_coboundary(columns_to_reduce.size());

// 			value_t filtration = get_filtration(column_to_reduce);

// #ifdef INDICATE_PROGRESS
// 			if ((i + 1) % 10000 == 0 || (i >= verbose_logging_threshold && (i + 1) % 1000 == 0)) {
// 				std::cout << "\033[K"
// 				          << "reducing column " << i + 1 << "/" << columns_to_reduce.size() << " (filtration "
// 				          << filtration << ", infinite bars: " << betti;
// 				if (betti_error > 0) std::cout << " (skipped " << betti_error << ")";
// 				std::cout << ")" << std::flush << "\r";
// 			}
// #endif

// 			index_t j = i;

// 			// start with a dummy pivot entry with coefficient -1 in order to initialize
// 			// working_coboundary with the coboundary of the simplex with index column_to_reduce
// 			real_filtration_entry_t pivot(0, -1, -1); //last param was -1 + modulus

// #ifdef ASSEMBLE_REDUCTION_MATRIX
// 			// initialize reduction_coefficients as identity matrix
// 			reduction_coefficients.append_column();
// 			reduction_coefficients.push_back(real_filtration_entry_t(column_to_reduce, 1));
// #else
// #ifdef USE_COEFFICIENTS
// 			reduction_coefficients.push_back(real_filtration_entry_t(column_to_reduce, 1));
// #endif
// #endif

// #ifndef SKIP_APPARENT_PAIRS
// 			bool might_be_apparent_pair = (i == size_t(j));
// #endif

// 			size_t iterations = 0;
// 			do {
// 				const coefficient_t factor = modulus - get_coefficient(pivot);

// #ifdef ASSEMBLE_REDUCTION_MATRIX
// 				auto coeffs_begin = reduction_coefficients.cbegin(j), coeffs_end = reduction_coefficients.cend(j);
// #else
// #ifdef USE_COEFFICIENTS
// 				auto coeffs_begin = &reduction_coefficients[j], coeffs_end = &reduction_coefficients[j] + 1;
// #else
// 				auto coeffs_begin = &columns_to_reduce[j], coeffs_end = &columns_to_reduce[j] + 1;
// #endif
// #endif
// 				for (auto it = coeffs_begin; it != coeffs_end; ++it) {
// 					real_filtration_entry_t cell = *it;
// 					set_coefficient(cell, get_coefficient(cell) * factor % modulus);

// #ifdef ASSEMBLE_REDUCTION_MATRIX
// 					reduction_column.push(cell);
// #endif
// 					coface_entries.clear();

// 					auto coboundary = complex.coboundary(cell);

// 					while (coboundary.has_next()) {
// 						real_filtration_entry_t coface = coboundary.next();

// 						if (get_filtration(coface) <= max_filtration) {
// #ifndef SKIP_APPARENT_PAIRS
// 							coface_entries.push_back(coface);
// 							if (might_be_apparent_pair && (get_filtration(cell) == get_filtration(coface))) {
// #ifdef USE_ARRAY_HASHMAP
// 								if (pivot_column_index[get_index(coface)] == INVALID_INDEX)
// #else
// 								if (pivot_column_index.find(get_index(coface)) == pivot_column_index.end())
// #endif
// 								{
// 									pivot = coface;
// 									goto found_persistence_pair;
// 								}
// 								might_be_apparent_pair = false;
// 							}
// #else
// 							iterations++;
// 							working_coboundary.push(coface);
// #endif
// 						}
// 					}

// #ifndef SKIP_APPARENT_PAIRS
// 					for (auto e : coface_entries) {
// 						iterations++;
// 						working_coboundary.push(e);
// 					}
// #endif
// 				}

// 				if (iterations > max_entries) {
// 					// Abort, this is too expensive
// 					if (generate_output) output->skipped_column(filtration);

// 					betti_error++;
// 					break;
// 				}

// 				pivot = working_coboundary.get_pivot();

// 				if (get_index(pivot) != -1) {
// #ifdef USE_ARRAY_HASHMAP
// 					auto pivot_column_idx = pivot_column_index[get_index(pivot)];

// 					if (pivot_column_idx != INVALID_INDEX) {
// 						j = pivot_column_idx;
// 						continue;
// 					}
// #else
// 					auto pair = pivot_column_index.find(get_index(pivot));

// 					if (pair != pivot_column_index.end()) {
// 						j = pair->second;
// 						continue;
// 					}
// #endif
// 				} else {
// 					if (generate_output) {
// 						output->new_infinite_barcode(filtration);
// 						betti++;

// 					}
// 					break;
// 				}

// #ifndef SKIP_APPARENT_PAIRS
// 			found_persistence_pair:
// #endif
// 				value_t death = get_filtration(pivot);
// 				if (generate_output && filtration != death) {
// 					output->new_barcode(filtration, death);

// 				}

// #ifdef USE_ARRAY_HASHMAP
// 				pivot_column_index[get_index(pivot)] = i;
// #else
// 				pivot_column_index.insert(std::make_pair(get_index(pivot), i));
// #endif

// #ifdef USE_COEFFICIENTS
// 				const coefficient_t inverse = multiplicative_inverse[get_coefficient(pivot)];
// #endif

// #ifdef ASSEMBLE_REDUCTION_MATRIX
// 				// replace current column of reduction_coefficients (with a single diagonal 1 entry)
// 				// by reduction_column (possibly with a different entry on the diagonal)
// 				reduction_coefficients.pop_back();
// 				while (true) {
// #ifdef USE_COEFFICIENTS
// 					real_filtration_entry_t e = pop_pivot(reduction_column, modulus);
// #else
// 					real_filtration_entry_t e = pop_pivot(reduction_column);
// #endif
// 					if (get_index(e) == -1) break;
// #ifdef USE_COEFFICIENTS
// 					set_coefficient(e, inverse * get_coefficient(e) % modulus);
// 					assert(get_coefficient(e) > 0);
// #endif
// 					reduction_coefficients.push_back(e);
// 				}
// #else
// #ifdef USE_COEFFICIENTS
// 				reduction_coefficients.pop_back();
// 				reduction_coefficients.push_back(real_filtration_entry_t(column_to_reduce, inverse));
// #endif
// #endif
// 				break;
// 			} while (true);
// 		}

// #ifdef INDICATE_PROGRESS
// 		std::cout << "\033[K";
// #endif
		
// 		return std::make_pair(betti, betti_error);
	}
};
