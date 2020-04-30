/*
 * This file belongs to the Galois project, a C++ library for exploiting parallelism.
 * The code is being released under the terms of the 3-Clause BSD License (a
 * copy is located in LICENSE.txt at the top-level directory).
 *
 * Copyright (C) 2018, The University of Texas at Austin. All rights reserved.
 * UNIVERSITY EXPRESSLY DISCLAIMS ANY AND ALL WARRANTIES CONCERNING THIS
 * SOFTWARE AND DOCUMENTATION, INCLUDING ANY WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR ANY PARTICULAR PURPOSE, NON-INFRINGEMENT AND WARRANTIES OF
 * PERFORMANCE, AND ANY WARRANTY THAT MIGHT OTHERWISE ARISE FROM COURSE OF
 * DEALING OR USAGE OF TRADE.  NO WARRANTY IS EITHER EXPRESS OR IMPLIED WITH
 * RESPECT TO THE USE OF THE SOFTWARE OR DOCUMENTATION. Under no circumstances
 * shall University be liable for incidental, special, indirect, direct or
 * consequential damages or loss of profits, interruption of business, or
 * related expenses which may arise from use of Software or Documentation,
 * including but not limited to those resulting from defects in Software and/or
 * Documentation, or loss or inaccuracy of data of any kind.
 */

#include "galois/Galois.h"
#include "galois/Reduction.h"
#include "galois/Bag.h"
#include "galois/Timer.h"
#include "galois/graphs/LCGraph.h"
#include "galois/ParallelSTL.h"
#include "llvm/Support/CommandLine.h"
#include "Lonestar/BoilerPlate.h"

#include "galois/runtime/Profile.h"
#include "def.h"
#include "AutoMiner.h"

#include <boost/iterator/transform_iterator.hpp>

#include <utility>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>

thread_local std::vector<VSPtr> AutoMiner::vsBuf;
thread_local std::vector<uint32_t*> AutoMiner::vsMemBuf;
thread_local std::vector<AutoMiner::PerLevelContext> AutoMiner::path;
thread_local unsigned int AutoMiner::curLevel= 0;
thread_local VSPtr AutoMiner::tempVS;

const char* name = "Pattern Mining";
const char* desc = "Mine an arbitrary pattern in a graph";
const char* url  = 0;

enum Algo {
    nodeiterator,
    edgeiterator
};

enum MiningApp {
  clique = 0,
  motif,
  pattern
};

namespace cll = llvm::cl;
static cll::opt<std::string>
        inputFilename(cll::Positional, cll::desc("<input file>"), cll::Required);
static cll::opt<Algo> algo(
        "algo", cll::desc("Choose an algorithm:"),
        cll::values(clEnumValN(Algo::nodeiterator, "nodeiterator", "Node Iterator"),
                    clEnumValN(Algo::edgeiterator, "edgeiterator", "Edge Iterator (default)")
        ),
        cll::init(Algo::edgeiterator));

static cll::opt<MiningApp> miningApp(
        "a", cll::desc("Choose a mining app:"),
        cll::values(clEnumValN(MiningApp::clique, "clique", "Clique counting"),
                    clEnumValN(MiningApp::motif, "motif", "Motif counting"),
                    clEnumValN(MiningApp::pattern, "pattern", "Pattern Matching")
        ),
        cll::init(MiningApp::clique));

static cll::opt<unsigned int> patternSize("k", cll::desc("Set a pattern size"), cll::init(0));

// static cll::opt<unsigned int> numThreads("t", cll::desc("Set the number of threads"), cll::init(0));

static cll::opt<std::string> patternFileName("patternFile", cll::desc("Specify the pattern filename"), cll::value_desc("pattern filename"), cll::init(""));

typedef galois::graphs::LC_CSR_Graph<uint32_t, void>::with_numa_alloc<
        true>::type ::with_no_lockable<true>::type Graph;

typedef Graph::GraphNode GNode;
/**
 * Like std::lower_bound but doesn't dereference iterators. Returns the first
 * element for which comp is not true.
 */
template <typename Iterator, typename Compare>
Iterator lowerBound(Iterator first, Iterator last, Compare comp) {
  Iterator it;
  typename std::iterator_traits<Iterator>::difference_type count, half;
  count = std::distance(first, last);
  while (count > 0) {
    it   = first;
    half = count / 2;
    std::advance(it, half);
    if (comp(it)) {
      first = ++it;
      count -= half + 1;
    } else {
      count = half;
    }
  }
  return first;
}

/**
 * std::set_intersection over edge_iterators.
 */
template <typename G>
size_t countEqual(G& g, typename G::edge_iterator aa,
                  typename G::edge_iterator ea, typename G::edge_iterator bb,
                  typename G::edge_iterator eb) {
  size_t retval = 0;
  while (aa != ea && bb != eb) {
    typename G::GraphNode a = g.getEdgeDst(aa);
    typename G::GraphNode b = g.getEdgeDst(bb);
    if (a < b) {
      ++aa;
    } else if (b < a) {
      ++bb;
    } else {
      retval += 1;
      ++aa;
      ++bb;
    }
  }
  return retval;
}

template <typename G>
struct LessThan {
  G& g;
  typename G::GraphNode n;
  LessThan(G& g, typename G::GraphNode n) : g(g), n(n) {}
  bool operator()(typename G::edge_iterator it) { return g.getEdgeDst(it) < n; }
};

template <typename G>
struct GreaterThanOrEqual {
  G& g;
  typename G::GraphNode n;
  GreaterThanOrEqual(G& g, typename G::GraphNode n) : g(g), n(n) {}
  bool operator()(typename G::edge_iterator it) {
    return !(n < g.getEdgeDst(it));
  }
};

template <typename G>
struct DegreeLess : public std::binary_function<typename G::GraphNode,
                                                typename G::GraphNode, bool> {
  typedef typename G::GraphNode N;
  G* g;
  DegreeLess(G& g) : g(&g) {}

  bool operator()(const N& n1, const N& n2) const {
    return std::distance(g->edge_begin(n1), g->edge_end(n1)) <
           std::distance(g->edge_begin(n2), g->edge_end(n2));
  }
};
template <typename G>
struct DegreeGreater : public std::binary_function<typename G::GraphNode,
                                                typename G::GraphNode, bool> {
  typedef typename G::GraphNode N;
  G* g;
  DegreeGreater(G& g) : g(&g) {}

  bool operator()(const N& n1, const N& n2) const {
    return std::distance(g->edge_begin(n1), g->edge_end(n1)) >
           std::distance(g->edge_begin(n2), g->edge_end(n2));
  }
};
template <typename G>
struct GetDegree
    : public std::unary_function<typename G::GraphNode, ptrdiff_t> {
  typedef typename G::GraphNode N;
  G* g;
  GetDegree(G& g) : g(&g) {}

  ptrdiff_t operator()(const N& n) const {
    return std::distance(g->edge_begin(n), g->edge_end(n));
  }
};

template <typename GraphNode, typename EdgeTy>
struct IdLess {
  bool
  operator()(const galois::graphs::EdgeSortValue<GraphNode, EdgeTy>& e1,
             const galois::graphs::EdgeSortValue<GraphNode, EdgeTy>& e2) const {
    return e1.dst < e2.dst;
  }
};

void nodeIteratingAlgo(Graph& graph) {

  galois::GAccumulator<size_t> numTriangles;

  galois::do_all(
      galois::iterate(graph),
      [&](const GNode& n) {
      },
      galois::chunk_size<CHUNK_SIZE>(), galois::steal(),
      galois::loopname("nodeIteratingAlgo"));

}

void edgeIteratingAlgo(Graph& graph) {
  struct WorkItem {
    GNode src;
    GNode dst;
    WorkItem(const GNode& a1, const GNode& a2) : src(a1), dst(a2) {}
  };

  galois::InsertBag<WorkItem> items;
  galois::GAccumulator<size_t> numTriangles;

  galois::do_all(galois::iterate(graph),
      [&](GNode n) {
      for (Graph::edge_iterator edge :
          graph.out_edges(n, galois::MethodFlag::UNPROTECTED)) {
      GNode dst = graph.getEdgeDst(edge);
      if (n < dst)
      items.push(WorkItem(n, dst));
      }
      },
      galois::loopname("Initialize"));

  galois::do_all(
      galois::iterate(items),
      [&](const WorkItem& w) {
      },
      galois::loopname("edgeIteratingAlgo"),
      galois::chunk_size<CHUNK_SIZE>(),
      galois::steal());
}

void makeGraph(Graph& graph, const std::string& triangleFilename) {
    typedef galois::graphs::FileGraph G;
    typedef G::GraphNode N;

    G initial, permuted;

    initial.fromFileInterleaved<void>(inputFilename);

    // Getting around lack of resize for deque
    std::deque<N> nodes;
    std::copy(initial.begin(), initial.end(), std::back_inserter(nodes));


    /* Sort by degree:
     *  DegreeLess: Sorts in the ascending order of node degrees
     *  DegreeGreater: Sorts in the descending order of the node degrees
     *
     *  The order of sorting has a huge impact on performance
     *  For this algorithm, sorting in descending order delivers the
     *  best performance due to the way ties are broken.
     */
    galois::ParallelSTL::sort(nodes.begin(), nodes.end(), DegreeGreater<G>(initial));

    std::deque<N> p;
    std::copy(nodes.begin(), nodes.end(), std::back_inserter(p));
    // Transpose
    size_t idx = 0;
    for (N n : nodes) {
        p[n] = idx++;
    }

    galois::graphs::permute<void>(initial, p, permuted);
    galois::do_all(galois::iterate(permuted),
                   [&](N x) { permuted.sortEdges<void>(x, IdLess<N, void>()); });

    std::cout << "Writing new input file: " << triangleFilename << "\n";
    permuted.toFile(triangleFilename);
    galois::graphs::readGraph(graph, permuted);
}

void readGraph(Graph& graph) {
  if (inputFilename.find(".gr.triangles") !=
      inputFilename.size() - strlen(".gr.triangles")) {
    // Not directly passed .gr.triangles file
    std::string triangleFilename = inputFilename + ".triangles";
    std::ifstream triangleFile(triangleFilename.c_str());
    if (!triangleFile.good()) {
      // triangles doesn't already exist, create it
      galois::StatTimer Trelabel("GraphRelabelTimer");
      galois::gPrint("WARNING: Sorted graph does not exist; Relabelling and Creating a sorted graph\n");
      Trelabel.start();
      makeGraph(graph, triangleFilename);
      Trelabel.stop();
    } else {
      // triangles does exist, load it
      galois::graphs::readGraph(graph, triangleFilename);
    }
  } else {
    galois::graphs::readGraph(graph, inputFilename);
  }

  size_t index = 0;
  for (GNode n : graph) {
    graph.getData(n) = index++;
  }
}

int main(int argc, char** argv) {
  galois::SharedMemSys G;
  LonestarStart(argc, argv, name, desc, url);
  galois::setActiveThreads(numThreads);

  Graph graph;

  galois::StatTimer Tinitial("GraphReadingTime");
  Tinitial.start();
  readGraph(graph);
  Tinitial.stop();

//   galois::preAlloc(numThreads + 16 * (graph.size() + graph.sizeEdges()) /
//       galois::runtime::pagePoolSize());
//   galois::reportPageAlloc("MeminfoPre");

  galois::StatTimer TGeneration("PlanGenerationTime");
  TGeneration.start();
  std::vector<Graphlet> all;
  switch(miningApp) {
    case clique:
      if(patternSize == 0) {
        std::cerr << "Pattern size must be a positive integer.\n";
        exit(1);
      }
      std::cout << "Clique counting with pattern size " << patternSize << std::endl;
      all.push_back(Graphlet::clique(patternSize));
      break;
    case motif:
      if(patternSize == 0) {
        std::cerr << "Pattern size must be a positive integer.\n";
        exit(1);
      }
      std::cout << "Motif counting with pattern size " << patternSize << std::endl;
      all = Graphlet::all_connected(patternSize);
      break;
    case pattern:
      if(patternFileName == "") {
        std::cerr << "Must specify a pattern file name.\n";
        exit(1);
      }
      std::cout << "Single pattern matching with file: " << patternFileName << std::endl;
      all.push_back(Graphlet::from_file(patternFileName));
      break;
    default:
      std::cerr << "Must specify a mining app.\n";
  }

  int total = 0;
  std::vector<std::vector<ExecutionPlan>> planss;
  for(unsigned int i=0; i< all.size(); ++i){
    std::vector<ExecutionPlan> plans;
    MultiRed mr(all.at(i),plans);
    planss.push_back(plans);
    total += plans.size();
  }

  //best plans
  MultiRestPlan mrp(0);

  for(unsigned int i=0;i<planss.size();++i){
    std::vector<ExecutionPlan> plans = planss.at(i);
    double bestcomplex = std::numeric_limits<double>::infinity();
    int bindex = 0;
    int windex = 0;
    double worstcomplex = 0;
    for(unsigned int j=0;j<plans.size();++j){
      RestPlan test(plans[j],i);
      double comple = test.time_complexity();
      //std::cout<<"plan "<<j<<" for graphlet "<<i<<" has expected complexity "<<comple<<std::endl;
      if(comple<bestcomplex){
        bindex=j;
        bestcomplex = comple;
      }
      if(comple>worstcomplex){
        windex=j;
        worstcomplex = comple;
      }
    }
    std::cout<<"Chose plan "<<bindex<<" as best for graphlet "<<i<<std::endl;
    mrp.add_ex_plan(plans[bindex],i);
  }

  TGeneration.stop();

  uint64_t maxDegree = 0;
  for(GNode node : graph) {
    Graph::edge_iterator ii = graph.edge_begin(node,galois::MethodFlag::UNPROTECTED);
    Graph::edge_iterator ie = graph.edge_end(node,galois::MethodFlag::UNPROTECTED);
    uint32_t numEdges = ie - ii;
    if(maxDegree < numEdges)
      maxDegree = numEdges;
  }

  galois::StatTimer TExe("ExecutionTime");
  TExe.start();
  std::vector<size_t> counters;
  AutoMiner am(&graph, mrp, maxDegree);
  counters = am.count();
  TExe.stop();

  std::cout << "Counters: ";
  for(auto c : counters) {
    std::cout << c << ", ";
  }
  std::cout<<"\b\b\n";

//   galois::reportPageAlloc("MeminfoPost");
  return 0;
}

