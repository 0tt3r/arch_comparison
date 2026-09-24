// Matrix-free fixed-(N_up,N_down) Hubbard evolution. No full-sector sparse matrix.
// Copy of hubbard_sector.cpp (hash-pinned by earlier results) with the site limit raised to 18.
// Build: g++ -O3 -std=c++17 -fopenmp cpp/hubbard_strip.cpp -o build/hubbard_strip
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <omp.h>
using C = std::complex<double>;
using V = std::vector<C>;
struct Pair { int a,b,sign; };
struct Link { int index,sign; };
struct Sector {
  int nx,ny,n,ku,kd; size_t nu,nd,dim;
  std::vector<std::pair<int,int>> sites,edges;
  std::vector<int> up,down,ui,di;
  std::vector<std::vector<Link>> hu,hd;
  size_t initial_index;
  Sector(int x,int y):nx(x),ny(y),n(x*y) {
    if(n>18 || n<2) throw std::runtime_error("Require 2 <= sites <= 18");
    int iu=0;
    for(int j=0;j<ny;j++) for(int k=0;k<nx;k++) {
      int i=j%2 ? nx-1-k:k; sites.push_back({i,j});
      if((i+j)%2==0) iu|=1<<(sites.size()-1);
    }
    ku=__builtin_popcount(unsigned(iu)); kd=n-ku;
    ui.assign(1<<n,-1);di=ui;
    for(int m=0;m<(1<<n);m++) {
      if(__builtin_popcount(unsigned(m))==ku) {ui[m]=up.size();up.push_back(m);}
      if(__builtin_popcount(unsigned(m))==kd) {di[m]=down.size();down.push_back(m);}
    }
    nu=up.size();nd=down.size();dim=nu*nd;
    initial_index=size_t(ui[iu])*nd+di[((1<<n)-1)^iu];
    for(int i=0;i<n;i++) for(auto xy: {std::make_pair(sites[i].first+1,sites[i].second),
                                      std::make_pair(sites[i].first,sites[i].second+1)}) {
      auto it=std::find(sites.begin(),sites.end(),xy);
      if(it!=sites.end()) edges.push_back({i,int(it-sites.begin())});
    }
    hu=links(up,ui);hd=links(down,di);
  }
  static int sign(int m,int i,int j) {
    int lo=std::min(i,j),hi=std::max(i,j);
    unsigned between=((1u<<hi)-1)^((1u<<(lo+1))-1);
    return __builtin_popcount(unsigned(m)&between)%2 ? -1:1;
  }
  std::vector<std::vector<Link>> links(const std::vector<int>& masks,const std::vector<int>& lookup) {
    std::vector<std::vector<Link>> result(masks.size());
    for(size_t a=0;a<masks.size();a++) for(auto e:edges) {
      int m=masks[a],i=e.first,j=e.second;
      if(((m>>i)^(m>>j))&1) result[a].push_back({lookup[m^(1<<i)^(1<<j)],-sign(m,i,j)});
    }
    return result;
  }
  V initial() const {V p(dim);p[initial_index]=1.;return p;}
  double potential(size_t a,size_t b) const {
    return 4.*(__builtin_popcount(unsigned(up[a]&down[b]))-.5*(ku+kd)+.25*n);
  }
  void hamiltonian(const V& p,V& out) const {
    #pragma omp parallel for schedule(static)
    for(size_t a=0;a<nu;a++) {
      C* q=out.data()+a*nd; const C* row=p.data()+a*nd;
      for(size_t b=0;b<nd;b++) q[b]=potential(a,b)*row[b];
      for(auto l:hu[a]) {
        const C* other=p.data()+size_t(l.index)*nd;
        for(size_t b=0;b<nd;b++) q[b]+=double(l.sign)*other[b];
      }
      for(size_t b=0;b<nd;b++) for(auto l:hd[b]) q[b]+=double(l.sign)*row[l.index];
    }
  }
  V exact(const std::string& coefficients,double scale) const {
    std::ifstream f(coefficients);double re,im;std::vector<C> coef;
    while(f>>re>>im) coef.emplace_back(re,im);
    if(coef.size()<2) throw std::runtime_error("Missing Chebyshev coefficients");
    V prev=initial(),cur(dim),next(dim),out(dim);
    hamiltonian(prev,cur);
    #pragma omp parallel for schedule(static)
    for(size_t i=0;i<dim;i++) {cur[i]/=scale;out[i]=coef[0]*prev[i]+coef[1]*cur[i];}
    for(size_t k=2;k<coef.size();k++) {
      hamiltonian(cur,next);
      #pragma omp parallel for schedule(static)
      for(size_t i=0;i<dim;i++) {next[i]=2./scale*next[i]-prev[i];out[i]+=coef[k]*next[i];}
      prev.swap(cur);cur.swap(next);
    }
    return out;
  }
  void onsite(V& p,double dt,int part=0) const {
    // part=0 all sites; 1 site zero only; 2 all except site zero.
    std::vector<C> phase(2*n+1);
    for(int k=-n;k<=n;k++) phase[k+n]=std::exp(C(0,-dt*k));
    #pragma omp parallel for schedule(static)
    for(size_t a=0;a<nu;a++) for(size_t b=0;b<nd;b++) {
      int v=int(potential(a,b));int first=((up[a]&1)==(down[b]&1))?1:-1;
      if(part==1) v=first; else if(part==2) v-=first;
      p[a*nd+b]*=phase[v+n];
    }
  }
  void pair_gate(V& p,int i,int j,double dt,int spin,bool hadamard=false) const {
    const auto& masks=spin ? down:up;const auto& lookup=spin ? di:ui;
    std::vector<Pair> pairs;std::vector<int> both;
    for(size_t a=0;a<masks.size();a++) {
      int m=masks[a];
      if((m&(1<<i)) && !(m&(1<<j))) pairs.push_back({int(a),lookup[m^(1<<i)^(1<<j)],sign(m,i,j)});
      if(hadamard && (m&(1<<i)) && (m&(1<<j))) both.push_back(a);
    }
    double co=std::cos(dt),si=std::sin(dt),root=1/std::sqrt(2.);
    if(spin==0) {
      #pragma omp parallel for schedule(static)
      for(size_t k=0;k<pairs.size();k++) {
        auto t=pairs[k];C* a=p.data()+size_t(t.a)*nd;C* b=p.data()+size_t(t.b)*nd;
        for(size_t z=0;z<nd;z++) {C x=a[z],y=b[z];
          if(hadamard) {a[z]=root*(x+double(t.sign)*y);b[z]=root*(double(t.sign)*x-y);}
          else {a[z]=co*x+C(0,si*t.sign)*y;b[z]=co*y+C(0,si*t.sign)*x;}
        }
      }
      if(hadamard) {
        #pragma omp parallel for schedule(static)
        for(size_t k=0;k<both.size();k++) for(size_t z=0;z<nd;z++) p[size_t(both[k])*nd+z]*=-1.;
      }
    } else {
      #pragma omp parallel for schedule(static)
      for(size_t k=0;k<nu;k++) {
        C* row=p.data()+k*nd;
        for(auto t:pairs) {C x=row[t.a],y=row[t.b];
          if(hadamard) {row[t.a]=root*(x+double(t.sign)*y);row[t.b]=root*(double(t.sign)*x-y);}
          else {row[t.a]=co*x+C(0,si*t.sign)*y;row[t.b]=co*y+C(0,si*t.sign)*x;}
        }
        if(hadamard) for(int b:both) row[b]*=-1.;
      }
    }
  }
  void group(V& p,int offset,double dt) const {
    for(int y=-offset;y<ny;y+=2) for(int x=-offset;x<nx;x+=2) {
      std::vector<int> v;
      for(auto xy:{std::make_pair(x,y),std::make_pair(x+1,y),std::make_pair(x+1,y+1),std::make_pair(x,y+1)}) {
        auto it=std::find(sites.begin(),sites.end(),xy);
        if(it!=sites.end()) v.push_back(it-sites.begin());
      }
      for(int spin=0;spin<2;spin++) {
        if(v.size()==2) pair_gate(p,v[0],v[1],dt,spin);
        if(v.size()==4) {
          pair_gate(p,v[2],v[0],0,spin,true);pair_gate(p,v[1],v[3],0,spin,true);
          pair_gate(p,v[1],v[2],2*dt,spin);
          pair_gate(p,v[2],v[0],0,spin,true);pair_gate(p,v[1],v[3],0,spin,true);
        }
      }
    }
  }
  V trotter(double tau,int r,const std::string& algorithm) const {
    V p=initial();double dt=tau/r;
    for(int step=0;step<r;step++) {
      if(algorithm=="plaquette") {
        onsite(p,dt/2);group(p,0,dt/2);group(p,1,dt);group(p,0,dt/2);onsite(p,dt/2);
      } else if(algorithm=="termwise") {
        onsite(p,dt/2,1);
        for(auto e:edges) for(int s=0;s<2;s++) pair_gate(p,e.first,e.second,dt/2,s);
        onsite(p,dt,2);
        for(auto it=edges.rbegin();it!=edges.rend();it++) for(int s=1;s>=0;s--) pair_gate(p,it->first,it->second,dt/2,s);
        onsite(p,dt/2,1);
      } else throw std::runtime_error("Unknown ordering");
    }
    return p;
  }
  void report(const V& p) const {
    double norm=0;std::vector<double> D(n,0);
    #pragma omp parallel
    {
      double local_norm=0;std::vector<double> local(n,0);
      #pragma omp for schedule(static)
      for(size_t a=0;a<nu;a++) for(size_t b=0;b<nd;b++) {
        double weight=std::norm(p[a*nd+b]);local_norm+=weight;
        unsigned mask=up[a]&down[b];
        while(mask) {int bit=__builtin_ctz(mask);local[bit]+=weight;mask&=mask-1;}
      }
      #pragma omp critical
      {norm+=local_norm;for(int i=0;i<n;i++) D[i]+=local[i];}
    }
    std::cout<<"\"norm_squared\":"<<norm<<",\"doublons\":[";
    for(int i=0;i<n;i++) std::cout<<(i?",":"")<<D[i];
    std::cout<<"]";
  }
};
int main(int argc,char** argv) {
  try {
    if(argc<7) throw std::runtime_error("usage: solver nx ny tau exact|plaquette|termwise|benchmark r coefficients [state.bin]");
    auto start=std::chrono::steady_clock::now();Sector model(std::stoi(argv[1]),std::stoi(argv[2]));
    double tau=std::stod(argv[3]);std::string mode=argv[4];int r=std::stoi(argv[5]);
    double scale=2*model.edges.size()+model.n;V p;
    if(mode=="exact") p=model.exact(argv[6],scale);
    else if(mode=="benchmark") {
      p=model.initial();V q(model.dim);
      for(int k=0;k<r;k++) {model.hamiltonian(p,q);p.swap(q);}
    } else {if(r<1) throw std::runtime_error("Require r>0");p=model.trotter(tau,r,mode);}
    auto end=std::chrono::steady_clock::now();
    std::cout<<std::setprecision(17)<<"{\"dimension\":"<<model.dim<<",\"threads\":"<<omp_get_max_threads()
             <<",\"seconds\":"<<std::chrono::duration<double>(end-start).count()<<",\"spectral_bound\":"<<scale<<",";
    model.report(p);std::cout<<"}\n";
    if(argc>7) {std::ofstream f(argv[7],std::ios::binary);f.write(reinterpret_cast<const char*>(p.data()),p.size()*sizeof(C));if(!f) throw std::runtime_error("State write failed");}
  } catch(const std::exception& e) {std::cerr<<e.what()<<"\n";return 1;}
}
