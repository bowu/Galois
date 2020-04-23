#ifndef MINER_HPP_
#define MINER_HPP_
#include "embedding.h"

#ifdef LARGE_SIZE
template <typename InTy = unsigned, typename OutTy = unsigned>
inline std::vector<OutTy> prefix_sum(const std::vector<InTy> &in) {
	std::vector<OutTy> sums(in.size() + 1);
	OutTy total = 0;
	for (size_t n = 0; n < in.size(); n++) {
		sums[n] = total;
		total += (OutTy)in[n];
	}
	sums[in.size()] = total;
	return sums;
}

template <typename InTy = unsigned, typename OutTy = unsigned>
inline std::vector<OutTy> parallel_prefix_sum(const std::vector<InTy> &in) {
    const size_t block_size = 1<<20;
    const size_t num_blocks = (in.size() + block_size - 1) / block_size;
    std::vector<OutTy> local_sums(num_blocks);
	// count how many bits are set on each thread
	galois::do_all(galois::iterate((size_t)0, num_blocks), [&](const size_t& block) {
		OutTy lsum = 0;
		size_t block_end = std::min((block + 1) * block_size, in.size());
		for (size_t i=block * block_size; i < block_end; i++)
			lsum += in[i];
		local_sums[block] = lsum;
	});
	std::vector<OutTy> bulk_prefix(num_blocks+1);
	OutTy total = 0;
	for (size_t block=0; block < num_blocks; block++) {
		bulk_prefix[block] = total;
		total += local_sums[block];
	}
	bulk_prefix[num_blocks] = total;
	std::vector<OutTy> prefix(in.size() + 1);
	galois::do_all(galois::iterate((size_t)0, num_blocks), [&](const size_t& block) {
		OutTy local_total = bulk_prefix[block];
		size_t block_end = std::min((block + 1) * block_size, in.size());
		for (size_t i=block * block_size; i < block_end; i++) {
			prefix[i] = local_total;
			local_total += in[i];
		}
	});
	prefix[in.size()] = bulk_prefix[num_blocks];
	return prefix;
}

#else
template <typename InTy = unsigned, typename OutTy = unsigned>
inline galois::gstl::Vector<OutTy> parallel_prefix_sum(const galois::gstl::Vector<InTy> &in) {
	galois::gstl::Vector<OutTy> sums(in.size() + 1);
	OutTy total = 0;
	for (size_t n = 0; n < in.size(); n++) {
		sums[n] = total;
		total += (OutTy)in[n];
	}
	sums[in.size()] = total;
	return sums;
} 
#endif

class Miner {
public:
	Miner() {}
	virtual ~Miner() {}
	// insert single-edge embeddings into the embedding queue (worklist)
	inline void insert(EmbeddingQueueType &queue) {
		galois::do_all(galois::iterate(graph->begin(), graph->end()),
			[&](const GNode& src) {
				#ifdef ENABLE_LABEL
				auto& src_label = graph->getData(src);
				#endif
				for (auto e : graph->edges(src)) {
					GNode dst = graph->getEdgeDst(e);
					#ifndef USE_DAG
					if(src >= dst) continue;
					#endif
					#ifdef ENABLE_LABEL
					auto& dst_label = graph->getData(dst);
					#endif
					EmbeddingType new_emb;
					#ifdef ENABLE_LABEL
					new_emb.push_back(ElementType(src, 0, src_label));
					new_emb.push_back(ElementType(dst, 0, dst_label));
					#else
					new_emb.push_back(ElementType(src));
					new_emb.push_back(ElementType(dst));
					#endif
					queue.push_back(new_emb);
				}
			},
			galois::chunk_size<CHUNK_SIZE>(), galois::steal(),
			galois::loopname("InitEmbQueue")
		);
		if(show) queue.printout_embeddings(0);
	}
	inline unsigned intersect(unsigned a, unsigned b) {
		return intersect_merge(a, b);
	}
	inline unsigned intersect_dag(unsigned a, unsigned b) {
		return intersect_dag_merge(a, b);
	}

protected:
	Graph *graph;
	unsigned max_size;
	std::vector<unsigned> degrees;
	std::vector<BYTE> is_wedge; // indicate a 3-vertex embedding is a wedge or chain (v0-cntered or v1-centered)

	#ifdef USE_QUERY_GRAPH
	std::vector<VertexId> matching_order;
	std::vector<VertexId> matching_order_map;
	std::vector<VertexId> automorph_group_id;
	// Read the preset file to hardcode the presets
	void read_presets() {
		std::ifstream ifile;
		ifile.open(preset_filename);
		if (!ifile) printf("Error in reading file %s\n", preset_filename.c_str());
		VertexId x;
		for (size_t i = 0; i< max_size; ++i) {
			ifile >> x;
			matching_order[i] = x;
			if(debug) std::cout << "matching_order[" << i << "] = " << x << "\n";
		}
		for (size_t i = 0; i < max_size; ++i) {
			ifile >> x;
			matching_order_map[i] = x;
			if(debug) std::cout << "matching_map[" << i << "] = " << x << "\n";
		}
		for (size_t i = 0; i < max_size; ++i) {
			ifile >> x;
			automorph_group_id[i] = x;
			if(debug) std::cout << "automorph_group_id[" << i << "] = " << x << "\n";
		}
		ifile.close();
	}
	#endif
	template <typename EmbeddingTy = VertexEmbedding>
	inline bool is_automorphism_dag(unsigned n, const EmbeddingTy& emb, unsigned idx, VertexId dst) {
		//if (dst <= emb.get_vertex(0)) return true;
		for (unsigned i = 0; i < n; ++i) if (dst == emb.get_vertex(i)) return true;
		for (unsigned i = 0; i < idx; ++i) if (is_connected_dag(dst, emb.get_vertex(i))) return true;
		//for (unsigned i = idx+1; i < n; ++i) if (dst < emb.get_vertex(i)) return true;
		return false;
	}
	template <typename EmbeddingTy = VertexEmbedding>
	inline bool is_vertexInduced_automorphism(unsigned n, const EmbeddingTy& emb, unsigned idx, VertexId dst) {
		//unsigned n = emb.size();
		// the new vertex id should be larger than the first vertex id
		if (dst <= emb.get_vertex(0)) return true;
		// the new vertex should not already exist in the embedding
		for (unsigned i = 1; i < n; ++i)
			if (dst == emb.get_vertex(i)) return true;
		// the new vertex should not already be extended by any previous vertex in the embedding
		for (unsigned i = 0; i < idx; ++i)
			if (is_connected(emb.get_vertex(i), dst)) return true;
		// the new vertex id should be larger than any vertex id after its source vertex in the embedding
		for (unsigned i = idx+1; i < n; ++i)
			if (dst < emb.get_vertex(i)) return true;
		return false;
	}
	inline unsigned find_motif_pattern_id(unsigned n, unsigned idx, VertexId dst, const VertexEmbedding& emb, unsigned pos = 0) {
		unsigned pid = 0;
		if (n == 2) { // count 3-motifs
			pid = 1; // 3-chain
			if (idx == 0) {
				if (is_connected(emb.get_vertex(1), dst)) pid = 0; // triangle
				#ifdef USE_WEDGE
				else if (max_size == 4) is_wedge[pos] = 1; // wedge; used for 4-motif
				#endif
			}
		} else if (n == 3) { // count 4-motifs
			unsigned num_edges = 1;
			pid = emb.get_pid();
			if (pid == 0) { // extending a triangle
				for (unsigned j = idx+1; j < n; j ++)
					if (is_connected(emb.get_vertex(j), dst)) num_edges ++;
				pid = num_edges + 2; // p3: tailed-triangle; p4: diamond; p5: 4-clique
			} else { // extending a 3-chain
				assert(pid == 1);
				std::vector<bool> connected(3, false);
				connected[idx] = true;
				for (unsigned j = idx+1; j < n; j ++) {
					if (is_connected(emb.get_vertex(j), dst)) {
						num_edges ++;
						connected[j] = true;
					}
				}
				if (num_edges == 1) {
					pid = 0; // p0: 3-path
					unsigned center = 1;
					#ifdef USE_WEDGE
					if (is_wedge[pos]) center = 0;
					#else
					center = is_connected(emb.get_vertex(1), emb.get_vertex(2)) ? 1 : 0;
					#endif
					if (idx == center) pid = 1; // p1: 3-star
				} else if (num_edges == 2) {
					pid = 2; // p2: 4-cycle
					unsigned center = 1;
					#ifdef USE_WEDGE
					if (is_wedge[pos]) center = 0;
					#else
					center = is_connected(emb.get_vertex(1), emb.get_vertex(2)) ? 1 : 0;
					#endif
					if (connected[center]) pid = 3; // p3: tailed-triangle
				} else {
					pid = 4; // p4: diamond
				}
			}
		} else { // count 5-motif and beyond
			pid = find_motif_pattern_id_eigen(n, idx, dst, emb);
		}
		return pid;
	}
	unsigned get_degree(Graph *g, VertexId vid) {
		return std::distance(g->edge_begin(vid), g->edge_end(vid));
	}
	void degree_counting() {
		degrees.resize(graph->size());
		galois::do_all(galois::iterate(graph->begin(), graph->end()), [&] (GNode v) {
			degrees[v] = std::distance(graph->edge_begin(v), graph->edge_end(v));
		}, galois::loopname("DegreeCounting"));
	}
	inline unsigned intersect_merge(unsigned src, unsigned dst) {
		unsigned count = 0;
		for (auto e : graph->edges(dst)) {
			GNode dst_dst = graph->getEdgeDst(e);
			for (auto e1 : graph->edges(src)) {
				GNode to = graph->getEdgeDst(e1);
				if (dst_dst == to) {
					count += 1;
					break;
				}
				if (to > dst_dst) break;
			}
		}
		return count;
	}
	inline unsigned intersect_dag_merge(unsigned p, unsigned q) {
		unsigned count = 0;
		auto p_start = graph->edge_begin(p);
		auto p_end = graph->edge_end(p);
		auto q_start = graph->edge_begin(q);
		auto q_end = graph->edge_end(q);
		auto p_it = p_start;
		auto q_it = q_start;
		int a;
		int b;
		while (p_it < p_end && q_it < q_end) {
			a = graph->getEdgeDst(p_it);
			b = graph->getEdgeDst(q_it);
			int d = a - b;
			if (d <= 0) p_it ++;
			if (d >= 0) q_it ++;
			if (d == 0) count ++;
		}
		return count;
	}
	inline unsigned intersect_search(unsigned a, unsigned b) {
		if (degrees[a] == 0 || degrees[b] == 0) return 0;
		unsigned count = 0;
		unsigned lookup = a;
		unsigned search = b;
		if (degrees[a] > degrees[b]) {
			lookup = b;
			search = a;
		} 
		Graph::edge_iterator begin = graph->edge_begin(search, galois::MethodFlag::UNPROTECTED);
		Graph::edge_iterator end = graph->edge_end(search, galois::MethodFlag::UNPROTECTED);
		for (auto e : graph->edges(lookup)) {
			GNode key = graph->getEdgeDst(e);
			if(binary_search(key, begin, end)) count ++;
		}
		return count;
	}
	inline bool is_all_connected_except(unsigned dst, unsigned pos, const BaseEmbedding &emb) {
		unsigned n = emb.size();
		bool all_connected = true;
		for(unsigned i = 0; i < n; ++i) {
			if (i == pos) continue;
			unsigned from = emb.get_vertex(i);
			if (!is_connected(from, dst)) {
				all_connected = false;
				break;
			}
		}
		return all_connected;
	}
	inline bool is_all_connected_except_dag(unsigned dst, unsigned pos, const BaseEmbedding &emb) {
		unsigned n = emb.size();
		bool all_connected = true;
		for(unsigned i = 0; i < n; ++i) {
			if (i == pos) continue;
			unsigned from = emb.get_vertex(i);
			if (!is_connected_dag(dst, from)) {
				all_connected = false;
				break;
			}
		}
		return all_connected;
	}
	inline bool is_all_connected(unsigned dst, const BaseEmbedding &emb, unsigned end, unsigned start = 0) {
		assert(start >= 0 && end > 0);
		bool all_connected = true;
		for(unsigned i = start; i < end; ++i) {
			unsigned from = emb.get_vertex(i);
			if (!is_connected(from, dst)) {
				all_connected = false;
				break;
			}
		}
		return all_connected;
	}
	inline bool is_all_connected_dag(unsigned dst, const BaseEmbedding &emb, unsigned end, unsigned start = 0) {
		assert(start >= 0 && end > 0);
		bool all_connected = true;
		for(unsigned i = start; i < end; ++i) {
			unsigned from = emb.get_vertex(i);
			if (!is_connected_dag(dst, from)) {
				all_connected = false;
				break;
			}
		}
		return all_connected;
	}
	inline bool is_all_connected_dag(unsigned dst, const std::vector<VertexId> &emb, unsigned end, unsigned start = 0) {
		assert(start >= 0 && end > 0);
		bool all_connected = true;
		for(unsigned i = start; i < end; ++i) {
			unsigned from = emb[i];
			if (!is_connected_dag(dst, from)) {
				all_connected = false;
				break;
			}
		}
		return all_connected;
	}

	// check if vertex a is connected to vertex b in a undirected graph
	inline bool is_connected(unsigned a, unsigned b) {
		if (degrees[a] == 0 || degrees[b] == 0) return false;
		unsigned key = a;
		unsigned search = b;
		if (degrees[a] < degrees[b]) {
			key = b;
			search = a;
		} 
		auto begin = graph->edge_begin(search, galois::MethodFlag::UNPROTECTED);
		auto end = graph->edge_end(search, galois::MethodFlag::UNPROTECTED);
		//return serial_search(key, begin, end);
		return binary_search(key, begin, end);
	}
	inline int is_connected_dag(unsigned key, unsigned search) {
		if (degrees[search] == 0) return false;
		auto begin = graph->edge_begin(search, galois::MethodFlag::UNPROTECTED);
		auto end = graph->edge_end(search, galois::MethodFlag::UNPROTECTED);
		//return serial_search(key, begin, end);
		return binary_search(key, begin, end);
	}
	inline bool serial_search(unsigned key, Graph::edge_iterator begin, Graph::edge_iterator end) {
		for (auto offset = begin; offset != end; ++ offset) {
			unsigned d = graph->getEdgeDst(offset);
			if (d == key) return true;
			if (d > key) return false;
		}
		return false;
	}
	inline bool binary_search(unsigned key, Graph::edge_iterator begin, Graph::edge_iterator end) {
		Graph::edge_iterator l = begin;
		Graph::edge_iterator r = end-1;
		while (r >= l) { 
			Graph::edge_iterator mid = l + (r - l) / 2; 
			unsigned value = graph->getEdgeDst(mid);
			if (value == key) return true;
			if (value < key) l = mid + 1; 
			else r = mid - 1; 
		} 
		return false;
	}
	inline int binary_search(unsigned key, Graph::edge_iterator begin, int length) {
		if (length < 1) return -1;
		int l = 0;
		int r = length-1;
		while (r >= l) { 
			int mid = l + (r - l) / 2; 
			unsigned value = graph->getEdgeDst(begin+mid);
			if (value == key) return mid;
			if (value < key) l = mid + 1; 
			else r = mid - 1; 
		} 
		return -1;
	}
	inline void gen_adj_matrix(unsigned n, const std::vector<bool> &connected, Matrix &a) {
		unsigned l = 0;
		for (unsigned i = 1; i < n; i++)
			for (unsigned j = 0; j < i; j++)
				if (connected[l++]) a[i][j] = a[j][i] = 1;
	}
	// calculate the trace of a given n*n matrix
	inline MatType trace(unsigned n, Matrix matrix) {
		MatType tr = 0;
		for (unsigned i = 0; i < n; i++) {
			tr += matrix[i][i];
		}
		return tr;
	}
	// matrix mutiplication, both a and b are n*n matrices
	inline Matrix product(unsigned n, const Matrix &a, const Matrix &b) {
		Matrix c(n, std::vector<MatType>(n));
		for (unsigned i = 0; i < n; ++i) { 
			for (unsigned j = 0; j < n; ++j) { 
				c[i][j] = 0; 
				for(unsigned k = 0; k < n; ++k) {
					c[i][j] += a[i][k] * b[k][j];
				}
			} 
		} 
		return c; 
	}
	// calculate the characteristic polynomial of a n*n matrix A
	inline void char_polynomial(unsigned n, Matrix &A, std::vector<MatType> &c) {
		// n is the size (num_vertices) of a graph
		// A is the adjacency matrix (n*n) of the graph
		Matrix C;
		C = A;
		for (unsigned i = 1; i <= n; i++) {
			if (i > 1) {
				for (unsigned j = 0; j < n; j ++)
					C[j][j] += c[n-i+1];
				C = product(n, A, C);
			}
			c[n-i] -= trace(n, C) / i;
		}
	}
	inline void get_connectivity(unsigned n, unsigned idx, VertexId dst, const VertexEmbedding &emb, std::vector<bool> &connected) {
		connected.push_back(true); // 0 and 1 are connected
		for (unsigned i = 2; i < n; i ++)
			for (unsigned j = 0; j < i; j++)
				if (is_connected(emb.get_vertex(i), emb.get_vertex(j)))
					connected.push_back(true);
				else connected.push_back(false);
		for (unsigned j = 0; j < n; j ++) {
			if (j == idx) connected.push_back(true);
			else if (is_connected(emb.get_vertex(j), dst))
				connected.push_back(true);
			else connected.push_back(false);
		}
	}
	// eigenvalue based approach to find the pattern id for a given embedding
	inline unsigned find_motif_pattern_id_eigen(unsigned n, unsigned idx, VertexId dst, const VertexEmbedding& emb) {
		std::vector<bool> connected;
		get_connectivity(n, idx, dst, emb, connected);
		Matrix A(n+1, std::vector<MatType>(n+1, 0));
		gen_adj_matrix(n+1, connected, A);
		std::vector<MatType> c(n+1, 0);
		char_polynomial(n+1, A, c);
		bliss::UintSeqHash h;
		for (unsigned i = 0; i < n+1; ++i)
			h.update((unsigned)c[i]);
		return h.get_value();
	}
};

#endif // MINER_HPP_
