#ifndef VERTEX_SET_H
#define VERTEX_SET_H
#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <memory>
#include <numa.h>

template<class NodeTy, bool isFromEdgeList=false, bool useNumaAlloc=false>
class VertexSet {
public:
  typedef NodeTy SizeType;

  template <bool _isFromEdgeList>
  struct with_from_edge_list {
    typedef VertexSet<NodeTy,_isFromEdgeList,useNumaAlloc> type;
  };

  template <bool _useNumaAlloc>
  struct with_numa_alloc {
    typedef VertexSet<NodeTy,isFromEdgeList,_useNumaAlloc> type;
  };

  //use an unprotected raw pointer to improve performance
  NodeTy *ptr=nullptr; 
  SizeType setSize=0;
  SizeType bufSize=0;
  NodeTy fromVid=-1; //prevent (0, 1, 0) to be recognized as a wedge
private: 
  bool   isLocalAlloc=false;

  VertexSet()=delete;
  VertexSet(const VertexSet&)=delete;
  VertexSet& operator=(const VertexSet&)=delete;

  void boundWithSize(SizeType s, NodeTy upper) {
    if(s > 64){
      NodeTy idx_l = -1;
      NodeTy idx_r = s;
      while(idx_r-idx_l > 1){
        NodeTy idx_t = (idx_l+idx_r)/2;
        if(ptr[idx_t] < upper) idx_l = idx_t;
        else idx_r = idx_t;
      }
      setSize = idx_l+1;
    }else{
      NodeTy idx_l = 0;

      while(idx_l < s && ptr[idx_l] < upper) ++idx_l;
      setSize = idx_l;
    }
  }

  //The following multiple functions need refactoring, but it should not affect
  //performance
  SizeType intersectBufNoBoundSize(NodeTy* otherPtr, SizeType otherSize) const
  {
    SizeType idx_l = 0, idx_r = 0, idx_out = 0;
    while(idx_l < setSize && idx_r < otherSize) {
      NodeTy left = ptr[idx_l];
      NodeTy right = otherPtr[idx_r];
      if(left <= right) idx_l++;
      if(right <= left) idx_r++;
      if(left == right) idx_out++;
    }
    return idx_out;
  }

  SizeType intersectBufNoBound(NodeTy* outPtr, NodeTy* otherPtr, SizeType otherSize)
  {
    SizeType idx_l = 0, idx_r = 0, idx_out = 0;
    while(idx_l < setSize && idx_r < otherSize) {
      NodeTy left = ptr[idx_l];
      NodeTy right = otherPtr[idx_r];
      if(left <= right) idx_l++;
      if(right <= left) idx_r++;
      if(left == right) outPtr[idx_out++] = left;
    }
    return idx_out;  
  }

  SizeType intersectBufWithBound(NodeTy* outPtr, NodeTy* otherPtr, SizeType otherSize, NodeTy upper)
  {
    SizeType idx_l = 0, idx_r = 0, idx_out = 0;
    while(idx_l < setSize && idx_r < otherSize) {
      NodeTy left = ptr[idx_l];
      NodeTy right = otherPtr[idx_r];
      if(left >= upper) break;
      if(right >= upper) break;
      if(left <= right) idx_l++;
      if(right <= left) idx_r++;
      if(left == right) outPtr[idx_out++] = left;
    }
    return idx_out;
  }

  NodeTy intersectBufWithBoundSize(NodeTy* otherPtr, SizeType otherSize, NodeTy upper) const {
    SizeType idx_l = 0, idx_r = 0, idx_out = 0;
    while(idx_l < setSize && idx_r < otherSize) {
      NodeTy left = ptr[idx_l];
      NodeTy right = otherPtr[idx_r];
      if(left >= upper) break;
      if(right >= upper) break;
      if(left <= right) idx_l++;
      if(right <= left) idx_r++;
      if(left == right) idx_out++;
    }
    return idx_out;
  }

  SizeType differenceBufNoBound(NodeTy* outPtr, NodeTy* otherPtr, SizeType otherSize, NodeTy otherVid)
  {
    SizeType idx_l = 0, idx_r = 0, idx_out = 0;
    while(idx_l < setSize && idx_r < otherSize) {
      NodeTy left = ptr[idx_l];
      NodeTy right = otherPtr[idx_r];
      if(left <= right) idx_l++;
      if(right <= left) idx_r++;
      if(isFromEdgeList) {
        if(left < right && left != otherVid) {
          outPtr[idx_out++] = left;
        }
      } else {
        if(left < right) {
          outPtr[idx_out++] = left;
        }
      }
    }
    while(idx_l < setSize) {
      NodeTy left = ptr[idx_l];
      idx_l++;
      if(isFromEdgeList) {
        if(left != otherVid) {
          outPtr[idx_out++] = left;
        }
      } else {
          outPtr[idx_out++] = left;
      }
    }
    return idx_out;
  }
      
  SizeType differenceBufNoBoundSize(NodeTy* otherPtr, SizeType otherSize, NodeTy otherVid) const {
    SizeType idx_l = 0, idx_r = 0, idx_out = 0;
    while(idx_l < setSize && idx_r < otherSize) {
      NodeTy left = ptr[idx_l];
      NodeTy right = otherPtr[idx_r];
      if(left <= right) idx_l++;
      if(right <= left) idx_r++;
      if(isFromEdgeList) {
        if(left < right && left != otherVid) {
          idx_out++;
        }
      } else {
        if(left < right)
          idx_out++;
      }
    }
    while(idx_l < setSize) {
      NodeTy left = ptr[idx_l];
      idx_l++;
      if(isFromEdgeList) {
        if(left != otherVid) {
          idx_out++;
        }
      } else {
          idx_out++;
      }
    }
    return idx_out;
  }

  SizeType differenceBufWithBound(NodeTy* outPtr, NodeTy* otherPtr, SizeType otherSize, NodeTy otherVid, NodeTy upper) {
    SizeType idx_l = 0, idx_r = 0, idx_out = 0;
    while(idx_l < setSize && idx_r < otherSize) {
      NodeTy left = ptr[idx_l];
      NodeTy right = otherPtr[idx_r];
      if(left >= upper) break;
      if(right >= upper) break;
      if(left <= right) idx_l++;
      if(right <= left) idx_r++;
      if(isFromEdgeList) {
        if(left < right && left != otherVid) outPtr[idx_out++] = left;
      } else {
        if(left < right) outPtr[idx_out++] = left;
      }
    }
    while(idx_l < setSize) {
      NodeTy left = ptr[idx_l];
      if(left >= upper) break;
      idx_l++;
      if(isFromEdgeList) {
        if(left != otherVid) {
          outPtr[idx_out++] = left;
        }
      } else {
          outPtr[idx_out++] = left;
      }
    }
    return idx_out;
  }

  SizeType differenceBufWithBoundSize(NodeTy *otherPtr, SizeType otherSize, NodeTy otherVid, NodeTy upper) const {
    SizeType idx_l = 0, idx_r = 0, idx_out = 0;
    while(idx_l < setSize && idx_r < otherSize) {
      NodeTy left = ptr[idx_l];
      NodeTy right = otherPtr[idx_r];
      if(left >= upper) break;
      if(right >= upper) break;
      if(left <= right) idx_l++;
      if(right <= left) idx_r++;
      if(isFromEdgeList) {
        if(left < right && left != otherVid) idx_out++;
      } else {
        if(left < right) idx_out++;
      }
    }
    while(idx_l < setSize) {
      NodeTy left = ptr[idx_l];
      if(left >= upper) break;
      idx_l++;
      if(isFromEdgeList) {
        if(left != otherVid) {
          idx_out++;
        }
      } else {
          idx_out++;
      }
    }
    return idx_out;
  }

public:
  //to store temporary resultant set
  VertexSet(SizeType s)
  {
    setSize = s;
    bufSize = s;
    if(useNumaAlloc) {
//       std::cout << "allocating numa array...\n";
      ptr = (NodeTy*)numa_alloc_local(sizeof(NodeTy) * s);
    } else {
      ptr = new NodeTy[s];
    }
  }

  VertexSet(NodeTy* edgeBegin, NodeTy edgeListSize, NodeTy vid) : ptr(edgeBegin), setSize(edgeListSize), bufSize(edgeListSize), fromVid(vid), isLocalAlloc(false) { }

  VertexSet(NodeTy* edgeBegin, NodeTy edgeListSize, NodeTy vid, NodeTy upper) : ptr(edgeBegin), bufSize(edgeListSize), fromVid(vid), isLocalAlloc(false) { 
    ptr = edgeBegin;
    boundWithSize(edgeListSize, upper);
  }

  VertexSet(VertexSet &other, NodeTy up) 
  {
    VertexSet(other.ptr, other.setSize, other.fromVid, up);
  }

  ~VertexSet() {
    if(isLocalAlloc) {
      if(useNumaAlloc)
        numa_free(ptr, sizeof(NodeTy)*bufSize);
      else
        delete[] ptr; 
    }
  }

  NodeTy* begin() { return ptr; }
  NodeTy* end() { return ptr+setSize; }

  void bound(NodeTy upper) {
    boundWithSize(setSize, upper);
  }


  VertexSet<NodeTy,false,useNumaAlloc>& intersectNoBound(VertexSet<NodeTy,false,useNumaAlloc> &dst, const VertexSet<NodeTy,true,useNumaAlloc> &other) {
    SizeType s = intersectBufNoBound(dst.ptr, other.ptr, other.setSize);
    dst.setSize = s;
    return dst;
  }

  VertexSet<NodeTy,false,useNumaAlloc>& intersectNoBound(VertexSet<NodeTy,false,useNumaAlloc> &dst, const VertexSet<NodeTy,false,useNumaAlloc> &other) {
    SizeType s = intersectBufNoBound(dst.ptr, other.ptr, other.setSize);
    dst.setSize = s;
    return dst;
  }

  SizeType intersectNoBoundSize(const VertexSet<NodeTy,true,useNumaAlloc> &other) const {
    return intersectBufNoBoundSize(other.ptr, other.setSize); 
  }

  SizeType intersectNoBoundSize(const VertexSet<NodeTy,false,useNumaAlloc> &other) const {
    return intersectBufNoBoundSize(other.ptr, other.setSize); 
  }

  VertexSet<NodeTy,false,useNumaAlloc>& intersectWithBound(VertexSet<NodeTy,false,useNumaAlloc> &dst, const VertexSet<NodeTy,true,useNumaAlloc> &other, NodeTy up) {
    dst.setSize = intersectBufWithBound(dst.ptr, other.ptr, other.setSize, up);
    return dst;
  }

  VertexSet<NodeTy,false,useNumaAlloc>& intersectWithBound(VertexSet<NodeTy,false,useNumaAlloc> &dst, const VertexSet<NodeTy,false,useNumaAlloc> &other, NodeTy up) {
    dst.setSize = intersectBufWithBound(dst.ptr, other.ptr, other.setSize, up);
    return dst;
  }

  SizeType intersectWithBoundSize(const VertexSet<NodeTy,true,useNumaAlloc> &other, NodeTy up) const {
    return intersectBufWithBoundSize(other.ptr, other.setSize, up);
  }
  
  SizeType intersectWithBoundSize(const VertexSet<NodeTy,false,useNumaAlloc> &other, NodeTy up) const {
    return intersectBufWithBoundSize(other.ptr, other.setSize, up);
  }

  VertexSet<NodeTy,false,useNumaAlloc>& differenceNoBound(VertexSet<NodeTy,false,useNumaAlloc> &dst, const VertexSet<NodeTy,true,useNumaAlloc> &other) {
    SizeType s = differenceBufNoBound(dst.ptr, other.ptr, other.setSize, other.fromVid);
    dst.setSize = s;
    return dst;
  }

  VertexSet<NodeTy,false,useNumaAlloc>& differenceNoBound(VertexSet<NodeTy,false,useNumaAlloc> &dst, const VertexSet<NodeTy,false,useNumaAlloc> &other) {
    SizeType s = differenceBufNoBound(dst.ptr, other.ptr, other.setSize, other.fromVid);
    dst.setSize = s;
    return dst;
  }

  SizeType differenceNoBoundSize(const VertexSet<NodeTy,true,useNumaAlloc> &other) const {
    return differenceBufNoBoundSize(other.ptr, other.setSize, other.fromVid); 
  }

  SizeType differenceNoBoundSize(const VertexSet<NodeTy,false,useNumaAlloc> &other) const {
    return differenceBufNoBoundSize(other.ptr, other.setSize, other.fromVid); 
  }

  VertexSet<NodeTy,false,useNumaAlloc>& differenceWithBound(VertexSet<NodeTy,false,useNumaAlloc> &dst, const VertexSet<NodeTy,true,useNumaAlloc> &other, NodeTy up) {
    dst.setSize = differenceBufWithBound(dst.ptr, other.ptr, other.setSize, other.fromVid, up);
    return dst;
  }

  VertexSet<NodeTy,false,useNumaAlloc>& differenceWithBound(VertexSet<NodeTy,false,useNumaAlloc> &dst, const VertexSet<NodeTy,false,useNumaAlloc> &other, NodeTy up) {
    dst.setSize = differenceBufWithBound(dst.ptr, other.ptr, other.setSize, other.fromVid, up);
    return dst;
  }

  SizeType differenceWithBoundSize(const VertexSet<NodeTy,true,useNumaAlloc> &other, NodeTy up) const {
    return differenceBufWithBoundSize(other.ptr, other.setSize, other.fromVid, up);
  }
  
  SizeType differenceWithBoundSize(const VertexSet<NodeTy,false,useNumaAlloc> &other, NodeTy up) const {
    return differenceBufWithBoundSize(other.ptr, other.setSize, other.fromVid, up);
  }

  void print() {
    std::copy(ptr, ptr+setSize, std::ostream_iterator<NodeTy>(std::cerr, " "));
    std::cerr << "\n fromVid: " << fromVid;
    std::cerr << std::endl;
  }
};


//TODO: Figure out why it doesn't compile
// template<typename NodeTy, bool useNumaAlloc>
// inline VertexSet<NodeTy,false,useNumaAlloc> makeVertexSet(NodeTy setSize) {
//  // typedef typename VertexSet<NodeTy>::with_from_edge_list<false>::type typename ::with_numa_alloc<useNumaAlloc>::type VS;
//   typedef typename VertexSet<NodeTy>::with_from_edge_list<false>::type VS;
//   return VS(setSize);
// }

// inline template<typename NodeTy>
// VertexSet<NodeTy,true,false> makeVertexSetFromEdgeList(NodeTy *edgeBegin, NodeTy edgeListSize, NodeTy vid) {
//   
//   return VertexSet<NodeTy>(edgeBegin, edgeListSize, vid);
// }
// 
// inline template<typename NodeTy>
// VertexSet<NodeTy,true,false> makeVertexSetFromEdgeListWithBound(NodeTy* ptr, NodeTy s, NodeTy bufSize, NodeTy vid, NodeTy upper) 
// {
//   return VertexSet<NodeTy,true>(ptr, s, bufSize, vid, upper);
// }
// 
// //This function should not be used to create vertex set from other which is created from
// //an edge list. Otherwise, (0, 1, 0) may be taken as a wedge.
// //This function can save space when multiple vertex sets are based on the same vertex set
// //but with different bounds
// inline template<typename NodeTy>
// VertexSet<NodeTy,false,false> makeVertexSetFromOtherWithBound(VertexSet &other, NodeTy upper)
// {
//   return VertexSet<NodeTy,false>(other.ptr, other.setSize, other.bufSize, -1, upper);
// }
 
#endif
