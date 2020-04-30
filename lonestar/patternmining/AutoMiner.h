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

typedef VertexSet<uint32_t> VSU32;
typedef std::unique_ptr<VSU32> VSPtr;

constexpr static const unsigned CHUNK_SIZE  = 64u;

//debug-begin
std::ostream& operator<<(std::ostream& os, VSU32& vs)
{
  std::cerr << "{";
  std::copy(vs.begin(), vs.end(), std::ostream_iterator<uint32_t>(std::cerr, ", "));
  std::cerr << "}\n";
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

  struct PerLevelContext {
    //which loopon it's in
    MultiRestPlan* plan;        
    //which loopon child (if any) it's processing
    std::map<RestSet,MultiRestPlan*>::iterator looponIt; 
    std::map<RestSet,MultiRestPlan*>::iterator looponEnd; 
    //which vertex in the loopon set it's processing
    GNode* vi; 
    //end of the set
    GNode* ve; 
  };


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

  VSPtr computeVertexSet(uint32_t node, const RestSet& rs) {
//     std::cerr << "\ncomputing vertex set for node: " << node << ", with rest set: " << rs;
    if(rs.depth == 0) {
      VSPtr result;
      if(rs.res_chain[0] == -1)
        result = getVSFromNode(node,UINT32_MAX);
      else if(rs.res_chain[0] == 0)
        result = getVSFromNode(node, node);
      else {
        std::cerr << "rs_chain[0] = " << rs.res_chain[0] << std::endl;
        exit(1);
      }
        
      return result;
    }

    VSPtr set2 = getVSFromNode(node,UINT32_MAX);
    int  bound     = rs.bound;
    VSPtr result = std::make_unique<VSU32>(vsMemBuf[rs.key],0,maxDegree,-1);

    switch(rs.op) {
      case RestSet::intersect:
      {
        VSPtr& set1 = vsBuf[rs.parentKey];
        if(bound == -1) {
          set1->intersectNoBound(*result, *set2); 
        } else {
          set1->intersectWithBound(*result, *set2, *path[bound].vi); 
        }
        break;
      }
      case RestSet::noParentDifference:
      {
        //TODO:check problem when depth = 1
        VSPtr s = getVSFromNode(*path[rs.out[0]].vi,UINT32_MAX); 
        set2->differenceNoBound(*result, *s);
        for(size_t i=1; i<rs.out.size(); ++i)
        {
          s = getVSFromNode(*path[rs.out[i]].vi, UINT32_MAX);
          result->differenceNoBound(*result, *s);
        }
        if(bound != -1)
          result->bound(*path[bound].vi);
        break;
      }
      case RestSet::difference: 
      {
        VSPtr& set1t = vsBuf[rs.parentKey];
        if(bound == -1) {
          set1t->differenceNoBound(*result, *set2); 
        } else {
          set1t->differenceWithBound(*result, *set2, *path[bound].vi); 
        }
        break;
      }
      default:
        std::cerr << "op is not supported";
        exit(1);
    }
    return result;
  }

  void computeVertexSetSize(uint32_t node, RestSet rs, uint32_t index) {
//     std::cerr << "\ncomputing vertex set size for node: " << node << ", with rest set: " << rs << std::endl;
    VSPtr set2 = getVSFromNode(node,UINT32_MAX);
//     std::cerr << *path[0].vi << ", " << node << std::endl;
    int  bound     = rs.bound;
    switch(rs.op) {
      case RestSet::intersect:
      {
        VSPtr& set1 = vsBuf[rs.parentKey];
        if(bound == -1) {
          accums[index] += set1->intersectNoBoundSize(*set2); 
        } else {
          accums[index] += set1->intersectWithBoundSize(*set2, *path[bound].vi); 
        }
        break;
      }
      case RestSet::noParentDifference:
      {
        VSPtr s = getVSFromNode(*path[rs.out[0]].vi, UINT32_MAX);
        set2->differenceNoBound(*tempVS, *s);
        for(size_t i=1; i<rs.out.size(); ++i)
        {
          s = getVSFromNode(*path[rs.out[i]].vi, UINT32_MAX);
          tempVS->differenceNoBound(*tempVS, *s);
        }
        if(bound != -1)
          tempVS->bound(*path[bound].vi);
        accums[index] += tempVS->setSize;
        break;
      }
      case RestSet::difference:
      {
        VSPtr& set1t = vsBuf[rs.parentKey];
        if(bound == -1) {
          accums[index] += set1t->differenceNoBoundSize(*set2); 
        } else {
          accums[index] += set1t->differenceWithBoundSize(*set2, *path[bound].vi); 
        }
        break;
      }
      default:
        std::cerr << "op is not supported";
        exit(1);
    }
  }

  inline VSPtr getVSFromNode(GNode node, GNode upper)
  {
    Graph::edge_iterator ii = graph->edge_begin(node,galois::MethodFlag::UNPROTECTED);
    Graph::edge_iterator ie = graph->edge_end(node,galois::MethodFlag::UNPROTECTED);

    //std::cerr << node << ": ";
//     Graph::edge_iterator it;
//     for(it = ii;it!=ie;++it)
//       std::cerr <<graph->getEdgeDst(it) << " "; 
//     std::cerr << std::endl;
    
    uint32_t numEdges = ie - ii;
    //assert(numEdges != 0);
    uint32_t* ptr = graph->getEdgeDstPtr(ii);
    VSPtr result;
    //std::cout << "\nnode: " << node << ", upper: " << upper << "*ptr = " << *ptr << "\n";
    if(upper == std::numeric_limits<GNode>::max())
      result = std::make_unique<VSU32>(ptr,numEdges,node);
    else
      result = std::make_unique<VSU32>(ptr,numEdges,node,upper);
    //std::cout << "getVSFromNode result: " << *result << "\n";
    return result;
  }

  void depthFirstMine()
  {
    galois::do_all(
      galois::iterate(*graph),
      [&](const GNode& n) {

        threadInit();

        curLevel = 0;
        VSPtr vs0 = std::make_unique<VSU32>(1);
        *(vs0->ptr) = n;
        path[0].vi = vs0->begin();
        path[0].ve = vs0->end();
        path[0].plan = &plan;
        bool terminate = false; 
        while(!terminate) {
          GNode node = *(path[curLevel].vi);
          MultiRestPlan *p = path[curLevel].plan; 
          for(const RestSet& rs : p->atlev) {
            vsBuf[rs.key] = computeVertexSet(node, rs); 
          }

          for(const auto& cp : p->counters) {
            computeVertexSetSize(node, cp.first, cp.second);
          }

          bool hasVertexInLoopon = false;
          for(std::map<RestSet,MultiRestPlan*>::iterator ii = p->children.begin(); ii!=p->children.end(); ii++) {
            const RestSet& loopon = ii->first;
            if(vsBuf[loopon.key] == nullptr) {
              vsBuf[loopon.key] = computeVertexSet(node, loopon);
            }
            if(vsBuf[loopon.key]->setSize != 0) {
              curLevel++;
              path[curLevel].looponIt = ii;
              path[curLevel].looponEnd = p->children.end();
              path[curLevel].plan = ii->second;
              path[curLevel].vi = vsBuf[loopon.key]->begin();
              path[curLevel].ve = vsBuf[loopon.key]->end();
              hasVertexInLoopon = true;
              break;
            }
          }

          //advance to the next node
          if(!hasVertexInLoopon){
            if(curLevel == 0) break;
            //no loopon at this level
            GNode *vi = path[curLevel].vi;
            GNode *ve = path[curLevel].ve;
            //handle duplicate dst
            while(vi+1 != ve) {
              if(*(vi+1) == *vi) {
                ++vi;
                continue;
              }
              else
                break;
            }
            ++vi;
            path[curLevel].vi = vi;

            while(vi == ve) {
              //move this before if to test size first
              while(++path[curLevel].looponIt != path[curLevel].looponEnd) {
                const RestSet& loopon = path[curLevel].looponIt->first;
                vsBuf[loopon.key] = computeVertexSet(*vi, loopon);
                if(vsBuf[loopon.key]->setSize != 0) {
                  path[curLevel].vi = vsBuf[loopon.key]->begin();
                  path[curLevel].ve = vsBuf[loopon.key]->end();
                  break;
                }
              }
              //all loopons are processed
              curLevel--;
              if(curLevel <= 0) {
                terminate = true;
                break;
              }
              vi = ++path[curLevel].vi;
              ve = path[curLevel].ve;
            }
          } 
        }
      },
      galois::chunk_size<CHUNK_SIZE>(), galois::steal(),
      galois::loopname("depthFirstAlgorithm"));
  }

  void threadInit() {
    path.resize(plan.totalDepth);
    vsBuf.resize(MultiRestPlan::allRestSets.size());
    vsMemBuf.resize(vsBuf.size());
    tempVS = std::make_unique<VSU32>(maxDegree);
    for(size_t i=0; i<vsMemBuf.size(); ++i)
      vsMemBuf[i] = (uint32_t*)malloc(maxDegree*sizeof(uint32_t));
  }

  Graph* graph;
  static thread_local std::vector<VSPtr> vsBuf;
  static thread_local std::vector<uint32_t*> vsMemBuf;
  static thread_local std::vector<PerLevelContext> path;
  static thread_local unsigned int curLevel;
  static thread_local VSPtr tempVS;
  MultiRestPlan& plan; 
  MiningAlgo algo = depthFirst;
  std::vector<uint64_t> memBudget;
  std::vector<galois::GAccumulator<size_t>> accums;
  uint64_t reserved = 2*1024*1024*1024L; //by default reserve 2G space for data other than vertexset;
  uint32_t maxDegree;
};

#endif
