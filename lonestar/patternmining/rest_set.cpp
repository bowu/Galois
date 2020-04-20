#include "def.h"
#include<iostream>

RestSet::RestSet(const std::vector<int>& i, const std::vector<int> &o,int rest) : ins(i), out(o), restrict(rest) {
  depth =0;
  if(i.size()>0) depth = std::max(depth, i[i.size()-1]);
  if(o.size()>0) depth = std::max(depth, o[o.size()-1]);
  varname = var_name();
}

bool RestSet::operator<(const RestSet& other) const{

  //ins size, outs size, ins values, outs values, restriction.
  if(other.ins.size() < ins.size())return true;
  if(other.ins.size() > ins.size())return false;
  if(other.out.size() < out.size())return true;
  if(other.out.size() > out.size())return false;

  for(int i=0;i<ins.size();++i){
    if(other.ins.at(i) < ins.at(i)) return true;
    if(other.ins.at(i) > ins.at(i)) return false;
  }
  for(int i=0;i<out.size();++i){
    if(other.out.at(i) < out.at(i)) return true;
    if(other.out.at(i) > out.at(i)) return false;
  }
  return restrict < other.restrict;//unrestricted happen first 
  // return false;
  //  return other.varname < varname;
}

std::string RestSet::var_name(){
  
  std::ostringstream oss;
  int ini =0;
  int oui =0;
  int i=0;
  /*
  if(ins.size() == 1 && restrict == -1 && out.size() == 0){
    oss<<"g.N(v0)";
    return oss.str();
    }*/
  for(i=0;i<ins.size()+out.size();++i){
    if(oui<out.size() &&out[oui]==i){
      oss<<"n"<<i;
      ++oui;
    }
    else{
      oss<<"y"<<i;
      ++ini;
    }
    //if(i==restrict) oss<<"f"<<restrict;
  }
  //ignore -1 restrict
  if(restrict>=0) oss<<"f"<<restrict;
  return oss.str();
}

//adjusting this function can change how performance ends up going.
RestSet RestSet::parent() const{
  //std::cout<<"CALLING PARENT FOR "<<varname<<std::endl;
  //highest variable value remaining
  int lastvar = ins[ins.size()-1];//
  if(out.size()>0)lastvar = std::max(lastvar,out[out.size()-1]);
  if(restrict>=0){
    return RestSet(ins,out,-1).parent();
  }
  if(ins[ins.size()-1] == lastvar){
    std::vector<int> nins(ins);
    nins.pop_back();
    return RestSet(nins,out,restrict);
  }
  else {
    std::vector<int> nout(out);
    nout.pop_back();
    return RestSet(ins,nout,restrict);
  }
  //return nullptr;
}

bool RestSet::valid(){
  return ins.size()>0;
}


//available is what is available at that level already
void RestSet::append_calc_to_stream(std::ostream& oss,bool numb,std::set<RestSet> &available) const{
  RestSet par = parent();
  if(!numb) {
    NOGPU(oss << "VertexSet " << varname << " = ");
  }
  if(!numb && restrict!=-1){
    RestSet testpar(ins,out,-1);
    if(available.count(testpar)){
      oss<<"bounded("<<testpar.varname<<",v"<<restrict<<");\n";
      return;
    }
  }
  //std::cout<<varname<<" is a child of "<<par.varname<<std::endl;
  if(par.ins.size() == 0) { // has transient unnamed intermediate sets
    //bounded for 0
    if(par.out.size()==0 && restrict == 0){
      oss<<"bounded(y0,v0);\n";
      return;
    }
    //we can't use the parent, it doesn't exist
    //just compute ourselves
    //oss<<"VertexSet "<<varname<< " = ";
    for(int i=0;i<out.size();++i) {
      oss << "difference_" << ((numb && i==0) ? "num(" : "set(");
    }
    if(out.size()==0) {
      NOGPU(oss << "g.N(v0)");
    } else oss<< "y" << ins[0];
    for(int i=0;i<out.size();++i){
      oss<<", y"<<out[i]<<"";
      if(i==0 && restrict!=-1) oss<<", v"<<restrict;
      oss<<")";
    }
    if(numb&& out.size()==0)oss<<".size()";
    oss<<";\n";
    return;
  }
  //if this is a bounded, and parent isn't, the difference is the bounding.
  /*if(restrict!=-1 &&  ){
    oss<<"VertexSet "<< varname << " = "
       << "bounded("<<par.varname<<", v"<<restrict<<");\n";
    return;
    }*/
  
  //this set is a difference

  if(out.size()!=par.out.size()){
    oss<< "difference_"<< (numb ? "num(" : "set(")
       << par.varname << ", y" << out[out.size()-1];
    if(restrict!=-1){
      oss<<", v"<<restrict; 
    }
    oss<<");\n";
    return;
  }
  //this set is an intersection
  if(ins.size()!=par.ins.size()){
    oss<< "intersection_" << (numb?"num(":"set(")
       <<par.varname<<", y"<<ins[ins.size()-1];
    if(restrict!=-1){
      oss<<", v"<<restrict; 
    }
    oss<<");\n";
    return;
  }  
}

double RestSet::data_complexity_ignoring_restrictions() const{
  std::set<int> inset(ins.begin(),ins.end());
  std::set<int> outset(out.begin(),out.end());
  std::pair<std::set<int>,std::set<int>> dar(inset,outset);
  double xp = expected_size(dar);
  if(std::isnan(xp)){
    std::cout<<"Nan from size" <<ins.size()<<" "<<out.size()<<" "<<restrict<<std::endl;
  }
  return expected_size(dar);
}

double RestSet::time_complexity_ignoring_restrictions() const {
  std::set<RestSet> emp;
  return time_complexity_ignoring_restrictions(false,emp);
}
double RestSet::time_complexity_ignoring_restrictions(bool numb,std::set<RestSet> &available) const{
  /*std::set<int> inset(ins.begin(),ins.end());
  std::set<int> outset(out.begin(),out.end());
  std::pair<std::set<int>,std::set<int>> dar(inset,outset);
  */
  if(!numb && restrict!=-1){
    RestSet testpar(ins,out,-1);
    if(available.count(testpar)==1){
      //basically nothing, we just do it bounded
      return 0;//std::min(64,restset.data_complexity_ignoring_restrictions());
    }
  }
  RestSet par = parent();
  if(par.ins.size()==0){
    //we compute here.
    //global average degree time 1 + expected size one intersect + expected size 2+ ... expected size num outs -1
    // = (global average degree - expected size numouts intersect)/(1-(1-p))
    double p = (double)GLOBAL_AVERAGE_DEGREE/(double)GLOBAL_VERTEX_COUNT;
    return GLOBAL_AVERAGE_DEGREE*(1-pow(1-p, out.size()))/p; 
  }
  
  //  std::cout<<data_complexity_ignoring_restrictions();
  return parent().data_complexity_ignoring_restrictions();
}
