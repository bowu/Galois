#ifndef PATTERN_MINING_H
#define PATTERN_MINING_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <unistd.h>
#include <map>
#include <unordered_map>
#include <set>
#include <memory>

#include "debug.h"
#include "galois/Galois.h"
#include "galois/Reduction.h"
#include "galois/Bag.h"
#include "galois/ParallelSTL.h"

#include "galois/runtime/Profile.h"
#include "def.h"
#include "VertexSet.h"

// #define UINT32_MAX std::numeric_limits<uint32_t>::max()


typedef galois::graphs::LC_CSR_Graph<uint32_t, void>::with_numa_alloc<true>::type 
          ::with_no_lockable<true>::type Graph;
typedef Graph::GraphNode GNode;

typedef VertexSet<uint32_t,true,true> TopLevelVS;
typedef VertexSet<uint32_t,false,true> OtherLevelVS;

typedef std::map<RestSet,std::unique_ptr<TopLevelVS> > TopLevelVSMap;
typedef std::map<RestSet,std::unique_ptr<OtherLevelVS> > OtherLevelVSMap;

constexpr static const unsigned CHUNK_SIZE  = 64u;

//debug-begin
std::ostream& operator<<(std::ostream& os, TopLevelVS& vs)
{
  std::cerr << "{";
  std::copy(vs.begin(), vs.end(), std::ostream_iterator<uint32_t>(std::cerr, ", "));
  std::cerr << "}\n";
  return os;
}

std::ostream& operator<<(std::ostream& os, OtherLevelVS& vs)
{
  std::cerr << "{";
  std::copy(vs.begin(), vs.end(), std::ostream_iterator<uint32_t>(std::cerr, ", "));
  std::cerr << "}\n";
  return os;
}

std::ostream& operator<<(std::ostream& os, TopLevelVSMap& vs)
{
  for( const auto& e : vs ) {
    std::cerr << "key: " << e.first;
    std::cerr << "value: " << *(e.second);
  }
  return os;
}

std::ostream& operator<<(std::ostream& os, OtherLevelVSMap& vs)
{
  for( const auto& e : vs ) {
    std::cerr << "key: " << e.first;
    std::cerr << "value: " << *(e.second);
  }
  return os;
}
template<class T>
void print_set(std::set<T> &v)
{
  std::copy(v.begin(), v.end(), std::ostream_iterator<T>(std::cerr, " "));
  std::cerr << std::endl;
}

std::ostream & operator<<(std::ostream & os, const RestSet& rs) 
{
  os << "RestSet: ins: (";
  std::copy(rs.ins.begin(), rs.ins.end(), std::ostream_iterator<int>(std::cerr, " "));
  os << ") out: (";
  std::copy(rs.out.begin(), rs.out.end(), std::ostream_iterator<int>(std::cerr, " "));
  os << ") restriction: (";
  for(auto i : rs.res_chain)
    os << i << ",";
  os << "), depth = " << rs.depth << ", key = " << rs.key << ", op = " << rs.op;
  return os;
}

std::ostream & operator<<(std::ostream & os, RestSet& rs) 
{
  os << "RestSet: ins: (";
  std::copy(rs.ins.begin(), rs.ins.end(), std::ostream_iterator<int>(std::cerr, " "));
  os << ") out: (";
  std::copy(rs.out.begin(), rs.out.end(), std::ostream_iterator<int>(std::cerr, " "));
  os << ") restriction: (";
  for(auto i : rs.res_chain)
    os << i << ",";
  os << "), depth = " << rs.depth << ", key = " << rs.key << ", op = " << rs.op;
  return os;
}

std::ostream & operator<<(std::ostream & os, RestPlan & rp)
{
  std::cerr << "rest: (";
  //print_vector<int>(rp.rest);
  std::cerr << ") \n loopns: \n";
  for(auto rs : rp.loopons)
    std::cerr << rs << "\n";
  //std::cerr << "depends:\n";
  int i=0; 
  for(auto d : rp.depends) {
    std::cerr << "depends[" << i++ << "]: ";
    for(auto rs : d) {
      std::cerr << "[" << rs << "]";
    }
    std::cerr << std::endl;
  }
  std::cerr << std::endl;
  return os;
}

//debug-end

enum MiningAlgo {
  breadthFirst = 0,
  depthFirst,
  balanced
};

class AutoMiner {
public:
  AutoMiner(Graph* g, MultiRestPlan& p, uint32_t d) : graph(g), plan(p), maxDegree(d) { 
    accums.resize(p.numPlans()); 
    for(auto p : MultiRestPlan::allRestSets) {
      std::cerr << "rs: " << p.first << "\n";
      std::cerr << "key: " << p.second << "\n";
    }
  }

  ~AutoMiner() {
//     for(size_t i=0; i<vsMemBuf.size(); ++i)
//       free(vsMemBuf[i]);
  }

  std::vector<size_t> count() {
    std::vector<size_t> result;
    switch(algo) {
      case depthFirst:
        depthFirstMine();
        break;
      case breadthFirst:
        break;
      case balanced:
        break;
      default:
        GALOIS_DIE("Choose a valid algorithm");
    }
    for(auto& c : accums)
      result.push_back(c.reduce());
    return result;
  }

  //TODO: To implement this function, we need to modify MultiRestPlan to not use counters
  //in the last level
  void enumerateToFile(std::string fname) const {
  }

  //Enumerate readable instances to stream
  void enumerateToStream(std::ostream& stream)
  {

  }

  void setMemBudget(std::vector<uint64_t> v) {
    assert(v.size() == plan.totalDepth);
    memBudget.resize(plan.totalDepth);
    std::copy(v.begin(), v.end(), memBudget.begin()); 
  }

  void setAlg(MiningAlgo a) { algo = a; }

private:
   
  void genMemBudget() {
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    long totalMemSize = pages * page_size - reserved;
    memBudget.resize(plan.totalDepth);
    memBudget[0] = 0;
    assert(memBudget.size() > 1);
    std::fill(memBudget.begin()+1,memBudget.end(),totalMemSize/(memBudget.size()-1));
  }

  std::unique_ptr<OtherLevelVS> computeVertexSet(uint32_t node, const RestSet& rs) {
//     std::cerr << "\ncomputing vertex set for node: " << node << ", with rest set: " << rs;
    std::unique_ptr<TopLevelVS> set2 = getVSFromNode(node,std::numeric_limits<uint32_t>::max());
//     RestSet parent = rs.parent();
    int  bound     = rs.bound;
    std::unique_ptr<OtherLevelVS> result = std::make_unique<OtherLevelVS>(vsMemBuf[rs.key],0,maxDegree,-1);
    if(rs.depth == 1) {
//       std::cerr<< "\nCurrently in topMap:" << topMap;
//       std::cerr << "\nparent is: " << parent;
      //std::unique_ptr<TopLevelVS>& set1 = topMap.at(parent);
//       std::cerr << "\n set1: " << *set1 << " set 2: " << *set2;

      switch(rs.op) {
        case RestSet::intersect:
        {
          std::unique_ptr<TopLevelVS>& set1 = topVsBuf[rs.parentKey];
          if(bound == -1) {
            set1->intersectNoBound(*result, *set2); 
          } else {
            set1->intersectWithBound(*result, *set2, path[bound]); 
          }
          break;
        }
        case RestSet::noParentDifference:
        {
          std::unique_ptr<TopLevelVS> s = getVSFromNode(path[rs.out[0]], std::numeric_limits<uint32_t>::max());
          if(bound == -1)
            set2->differenceNoBound(*result, *s);
          else
            set2->differenceWithBound(*result, *s, path[bound]); 
          break;
        }
        case RestSet::difference:
        {
          std::unique_ptr<TopLevelVS>& set1t = topVsBuf[rs.parentKey];
          if(bound == -1) {
            set1t->differenceNoBound(*result, *set2); 
          } else {
            set1t->differenceWithBound(*result, *set2, path[bound]); 
          }
          break;
        }
        default:
          std::cerr << "op is not supported";
          exit(1);
      }
      return result;
    } else {
//       std::cerr<< "\nCurrently in otherMap:" << otherMap;
      //std::unique_ptr<OtherLevelVS>& set1 = otherMap.at(parent);
//       std::cerr << "\nparent is: " << parent;
//       std::cerr << "\n set1: " << *set1 << " set 2: " << *set2;
      std::unique_ptr<OtherLevelVS> result = std::make_unique<OtherLevelVS>(vsMemBuf[rs.key],0,maxDegree,-1);

      switch(rs.op) {
        case RestSet::intersect:
        {
          std::unique_ptr<OtherLevelVS>& set1 = otherVsBuf[rs.parentKey];
          if(bound == -1) {
            set1->intersectNoBound(*result, *set2); 
          } else {
            set1->intersectWithBound(*result, *set2, path[bound]); 
          }
          break;
        }
        case RestSet::noParentDifference:
        {
          std::unique_ptr<TopLevelVS> s = getVSFromNode(path[rs.out[0]], std::numeric_limits<uint32_t>::max());
          set2->differenceNoBound(*result, *s);
          for(size_t i=1; i<rs.out.size(); ++i)
          {
            s = getVSFromNode(path[rs.out[i]], std::numeric_limits<uint32_t>::max());
            result->differenceNoBound(*result, *s);
          }
          if(bound != -1)
            result->bound(path[bound]);
          break;
        }
        case RestSet::difference: 
        {
          std::unique_ptr<OtherLevelVS>& set1t = otherVsBuf[rs.parentKey];
          if(bound == -1) {
            set1t->differenceNoBound(*result, *set2); 
          } else {
            set1t->differenceWithBound(*result, *set2, path[bound]); 
          }
          break;
        }
        default:
          std::cerr << "op is not supported";
          exit(1);
      }
      return result;
    }
  }

  //void computeVertexSetSize(uint32_t node, RestSet rs, uint32_t index, TopLevelVSMap& topMap, OtherLevelVSMap& otherMap) {
  void computeVertexSetSize(uint32_t node, RestSet rs, uint32_t index) {
//     std::cerr << "\ncomputing vertex set size for node: " << node << ", with rest set: " << rs << std::endl;
    std::unique_ptr<TopLevelVS> set2 = getVSFromNode(node,std::numeric_limits<uint32_t>::max());
    int  bound     = rs.bound;
    if(rs.depth == 1) {
//       std::cerr<< "Currently in topMap:" << topMap << std::endl;
      //std::unique_ptr<TopLevelVS>& set1 = topMap.at(parent);
//       std::cerr << "\n set1: " << *set1 << " set 2: " << *set2;
//       std::cerr << "\n counter: " << accums[0].reduce();
      switch(rs.op) {
        case RestSet::intersect:
        {
          std::unique_ptr<TopLevelVS>& set1 = topVsBuf[rs.parentKey];
          if(bound == -1) {
            accums[index] += set1->intersectNoBoundSize(*set2); 
          } else {
            accums[index] += set1->intersectWithBoundSize(*set2, path[bound]); 
          }
          break;
        }
        case RestSet::noParentDifference:
        {
          std::unique_ptr<TopLevelVS> s = getVSFromNode(path[rs.out[0]], std::numeric_limits<uint32_t>::max());
          if(bound == -1)
            accums[index] += set2->differenceNoBoundSize(*s);
          else
            accums[index] += set2->differenceWithBoundSize(*s, path[bound]);
          break;
        }
        case RestSet::difference:
        {
          std::unique_ptr<TopLevelVS>& set1t = topVsBuf[rs.parentKey];
          if(bound == -1) {
            accums[index] += set1t->differenceNoBoundSize(*set2); 
          } else {
            accums[index] += set1t->differenceWithBoundSize(*set2, path[bound]); 
          }
          break;
        }
        default:
          std::cerr << "op is not supported";
          exit(1);
      }
//       std::cerr << "\n counter changes to : " << accums[0].reduce();
    } else {
//       std::cerr<< "Currently in otherMap:" << otherMap << std::endl;
      //std::unique_ptr<OtherLevelVS>& set1 = otherMap.at(parent);
//       std::cerr << "\n set1: " << *set1 << " set 2: " << *set2;
      switch(rs.op) {
        case RestSet::intersect:
        {
          std::unique_ptr<OtherLevelVS>& set1 = otherVsBuf[rs.parentKey];
          if(bound == -1) {
            accums[index] += set1->intersectNoBoundSize(*set2); 
          } else {
            accums[index] += set1->intersectWithBoundSize(*set2, path[bound]); 
          }
          break;
        }
        case RestSet::noParentDifference:
        {
          std::unique_ptr<TopLevelVS> s = getVSFromNode(path[rs.out[0]], std::numeric_limits<uint32_t>::max());
          set2->differenceNoBound(*tempVS, *s);
          for(size_t i=1; i<rs.out.size(); ++i)
          {
            s = getVSFromNode(path[rs.out[i]], std::numeric_limits<uint32_t>::max());
            tempVS->differenceNoBound(*tempVS, *s);
          }
          if(bound != -1)
            tempVS->bound(path[bound]);
          accums[index] += tempVS->setSize;
          break;
        }
        case RestSet::difference:
        {
//           std::cerr << "\nrs.depth" << rs.depth << "rs.key: " << rs.key << "rs.parentKey: " << rs.parentKey <<"\n";
          std::unique_ptr<OtherLevelVS>& set1t = otherVsBuf[rs.parentKey];
          if(bound == -1) {
            accums[index] += set1t->differenceNoBoundSize(*set2); 
          } else {
            accums[index] += set1t->differenceWithBoundSize(*set2, path[bound]); 
          }
//           std::cerr << "********";
          break;
        }
        default:
          std::cerr << "op is not supported";
          exit(1);
      }
    }
  }

  //void depthRecurse(GNode node, std::map<RestSet,MultiRestPlan*>& m, TopLevelVSMap& topMap, OtherLevelVSMap& otherMap){
  void depthRecurse(GNode node, std::map<RestSet,MultiRestPlan*>& m){
    path[++depth] = node;
    //std::cout << node << ", ";
    for(auto p : m)
    {
      RestSet loopon = p.first;
      std::set<RestSet>& atlev = p.second->atlev;
      std::map<RestSet,int>& counters = p.second->counters;

//       std::cerr << "Recurse on: << " << node << ", loopon is: " << loopon;

      //otherMap[loopon] = computeVertexSet(node, loopon, topMap, otherMap); 
      //if(otherVsBuf[loopon.key] == nullptr)
      otherVsBuf[loopon.key] = computeVertexSet(node, loopon);
//       std::cerr << "\n depth: " << depth << std::endl << "loopon: " << loopon << ", key: " << loopon.key << std::endl;
      for(const RestSet rs : atlev) {
        //otherMap[rs] = computeVertexSet(node, rs, topMap, otherMap);
//         std::cerr << "atlev: " << rs << std::endl;
        otherVsBuf[rs.key] = computeVertexSet(node, rs);
      }

      //std::cerr<< otherMap;

      for(uint32_t *ptr = otherVsBuf[loopon.key]->begin(); ptr!=otherVsBuf[loopon.key]->end(); ++ptr) {
        depthRecurse(*ptr, p.second->children);
        for(auto cp : counters) {
          computeVertexSetSize(*ptr, cp.first, cp.second);
        }
      }

      //otherMap.erase(loopon);
      otherVsBuf[loopon.key] = nullptr;
      for(RestSet rs : atlev) {
        //otherMap.erase(rs);
        otherVsBuf[rs.key] = nullptr;;
      }
    }
    --depth;
  }
  
  inline std::unique_ptr<TopLevelVS> getVSFromNode(GNode node, GNode upper)
  {
    Graph::edge_iterator ii = graph->edge_begin(node,galois::MethodFlag::UNPROTECTED);
    Graph::edge_iterator ie = graph->edge_end(node,galois::MethodFlag::UNPROTECTED);
    
    uint32_t numEdges = ie - ii;
    //assert(numEdges != 0);
    uint32_t* ptr = graph->getEdgeDstPtr(ii);
    std::unique_ptr<TopLevelVS> result;
    //std::cout << "\nnode: " << node << ", upper: " << upper << "*ptr = " << *ptr << "\n";
    if(upper == std::numeric_limits<GNode>::max())
      result = std::make_unique<TopLevelVS>(ptr,numEdges,node);
    else
      result = std::make_unique<TopLevelVS>(ptr,numEdges,node,upper);
    //std::cout << "getVSFromNode result: " << *result << "\n";
    return result;
  }

  void depthFirstMine()
  {
    galois::do_all(
      galois::iterate(*graph),
      [&](const GNode& n) {
//         TopLevelVSMap topVsMap;
        //may improve performance to use a vector of maps
//         std::vector<OtherLevelVSMap> vsMap(totalDepth-1);
//         OtherLevelVSMap vsMap; 
        //at top level, we should have 1 or 2 restsets

        threadInit();
        for(RestSet rs : plan.atlev) {
          //std::cout << "*********************\n";
          if(rs.res_chain[0] == -1) {
            //topVsMap[rs] = getVSFromNode(n, std::numeric_limits<GNode>::max());
            topVsBuf[rs.key] = getVSFromNode(n, std::numeric_limits<GNode>::max());
          } else {
            //topVsMap[rs] = getVSFromNode(n,n);
            topVsBuf[rs.key] = getVSFromNode(n,n);
          }
          //std::cout << rs;
          //std::cout << *(topVsMap[rs]);
          //std::cout << "*********************\n";
        }


        //std::cerr<< topVsMap;
//         for( const auto& e : topVsMap) {
//           std::cerr << "key: " << e.first;
//           std::cerr << "value: " << *(e.second);
//         }

        MultiRestPlan* mp;
//         std::cout << "\nStart: " << n << ", ";
        path[0] = n;

        for( auto c : plan.children) {
          RestSet loopon = c.first;
          std::map<RestSet,int>& counters = c.second->counters;
          mp = c.second;
          std::unique_ptr<TopLevelVS>& loop = topVsBuf[loopon.key];
          for(uint32_t *ptr = loop->begin(); ptr!=loop->end(); ++ptr) {
            //depthRecurse(*ptr, mp->children, topVsMap, vsMap);
            path[++depth] = *ptr;
            for(const RestSet rs : mp->atlev) {
              //otherMap[rs] = computeVertexSet(node, rs, topMap, otherMap);
//               std::cerr << "atlev: " << rs << std::endl;
              otherVsBuf[rs.key] = computeVertexSet(*ptr, rs);
            }
            depthRecurse(*ptr, mp->children);
            for(auto cp : counters) {
              //computeVertexSetSize(*ptr, cp.first, cp.second, topVsMap, vsMap);
              computeVertexSetSize(*ptr, cp.first, cp.second);
            }
            --depth;
          }
        }
      },
      galois::chunk_size<CHUNK_SIZE>(), galois::steal(),
      galois::loopname("depthFirstAlgorithm"));
  }

  void threadInit() {
    path.resize(plan.totalDepth);
    topVsBuf.resize(MultiRestPlan::allRestSets.size());
    otherVsBuf.resize(MultiRestPlan::allRestSets.size());
    vsMemBuf.resize(otherVsBuf.size());
    depth = 0;
    tempVS = std::make_unique<OtherLevelVS>(maxDegree);
    for(size_t i=0; i<vsMemBuf.size(); ++i)
      vsMemBuf[i] = (uint32_t*)malloc(maxDegree*sizeof(uint32_t));
  }

  Graph* graph;
  static thread_local std::vector<std::unique_ptr<TopLevelVS> > topVsBuf;
  static thread_local std::vector<std::unique_ptr<OtherLevelVS> > otherVsBuf;
  static thread_local std::vector<uint32_t*> vsMemBuf;
  static thread_local std::vector<uint32_t> path;
  static thread_local unsigned int depth;
  static thread_local std::unique_ptr<OtherLevelVS> tempVS;
  MultiRestPlan& plan; 
  MiningAlgo algo = depthFirst;
  std::vector<uint64_t> memBudget;
  std::vector<galois::GAccumulator<size_t>> accums;
  uint64_t reserved = 2*1024*1024*1024L; //by default reserve 2G space for data other than vertexset;
  uint32_t maxDegree;
};

#endif
