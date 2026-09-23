// Reuse the validated matrix-free kernels without changing the frozen solver.
#define main frozen_hubbard_main
#include "hubbard_sector.cpp"
#undef main

struct ProductSector : Sector {
  ProductSector(int x,int y,int iu,int id):Sector(x,y) {
    if(iu<0 || id<0 || (iu>>n) || (id>>n)) throw std::runtime_error("Invalid input masks");
    ku=__builtin_popcount(unsigned(iu));kd=__builtin_popcount(unsigned(id));
    up.clear();down.clear();ui.assign(1<<n,-1);di=ui;
    for(int m=0;m<(1<<n);m++) {
      if(__builtin_popcount(unsigned(m))==ku) {ui[m]=up.size();up.push_back(m);}
      if(__builtin_popcount(unsigned(m))==kd) {di[m]=down.size();down.push_back(m);}
    }
    nu=up.size();nd=down.size();dim=nu*nd;
    initial_index=size_t(ui[iu])*nd+di[id];hu=links(up,ui);hd=links(down,di);
  }
  std::vector<std::pair<int,int>> read_order(const std::string& path) const {
    std::ifstream f(path);int i,j;std::vector<std::pair<int,int>> order;
    while(f>>i>>j) order.push_back({i,j});
    auto sorted=order,expected=edges;
    for(auto* v:{&sorted,&expected}) {for(auto& e:*v) if(e.first>e.second) std::swap(e.first,e.second);std::sort(v->begin(),v->end());}
    if(sorted!=expected) throw std::runtime_error("Order must contain each physical edge once");
    return order;
  }
  V circuit(double tau,int r,const std::string& path) const {
    auto order=read_order(path);V p=initial();double dt=tau/r;
    for(int k=0;k<r;k++) {
      for(auto e:order) for(int spin=0;spin<2;spin++) pair_gate(p,e.first,e.second,dt/2,spin);
      onsite(p,dt);
      for(auto it=order.rbegin();it!=order.rend();it++) for(int spin=1;spin>=0;spin--) pair_gate(p,it->first,it->second,dt/2,spin);
    }
    return p;
  }
  V evolve(const V& input,const std::string& file,double scale) const {
    std::ifstream f(file);double re,im;std::vector<C> coef;
    while(f>>re>>im) coef.emplace_back(re,im);
    if(coef.size()<2) throw std::runtime_error("Missing Chebyshev coefficients");
    V prev=input,cur(dim),next(dim),out(dim);
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
  void observables(const V& p) const {
    double norm=0;std::vector<double> D(n),pu(nu),pd(nd);
    #pragma omp parallel
    {
      double ln=0;std::vector<double> ld(n),lpd(nd);
      #pragma omp for schedule(static)
      for(size_t a=0;a<nu;a++) {
        double row=0;
        for(size_t b=0;b<nd;b++) {
          double w=std::norm(p[a*nd+b]);row+=w;lpd[b]+=w;
          unsigned mask=up[a]&down[b];
          while(mask) {int i=__builtin_ctz(mask);ld[i]+=w;mask&=mask-1;}
        }
        pu[a]=row;ln+=row;
      }
      #pragma omp critical
      {norm+=ln;for(int i=0;i<n;i++) D[i]+=ld[i];for(size_t b=0;b<nd;b++) pd[b]+=lpd[b];}
    }
    std::vector<double> zu(n,1),zd(n,1);
    for(int i=0;i<n;i++) {
      for(size_t a=0;a<nu;a++) if((up[a]>>i)&1) zu[i]-=2*pu[a];
      for(size_t b=0;b<nd;b++) if((down[b]>>i)&1) zd[i]-=2*pd[b];
    }
    std::cout<<"\"norm_squared\":"<<norm<<",\"doublons\":[";
    for(int i=0;i<n;i++) std::cout<<(i?",":"")<<D[i];
    std::cout<<"],\"z_up\":[";
    for(int i=0;i<n;i++) std::cout<<(i?",":"")<<zu[i];
    std::cout<<"],\"z_down\":[";
    for(int i=0;i<n;i++) std::cout<<(i?",":"")<<zd[i];
    std::cout<<"]";
  }
};

int main(int argc,char** argv) {
  try {
    if(argc<10) throw std::runtime_error("usage: willow_sector nx ny upmask downmask tau exact|circuit|trajectory r_or_times coeff order [state.bin]");
    ProductSector model(std::stoi(argv[1]),std::stoi(argv[2]),std::stoi(argv[3]),std::stoi(argv[4]));
    double tau=std::stod(argv[5]),scale=2*model.edges.size()+model.n;
    std::string mode=argv[6];int r=std::stoi(argv[7]);
    auto start=std::chrono::steady_clock::now();V p;
    auto report=[&](double t){
      double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
      std::cout<<std::setprecision(17)<<"{\"dimension\":"<<model.dim<<",\"threads\":"<<omp_get_max_threads()
               <<",\"tau\":"<<t<<",\"seconds\":"<<elapsed<<",\"spectral_bound\":"<<scale<<",";
      model.observables(p);std::cout<<"}\n"<<std::flush;
    };
    if(mode=="trajectory") {
      if(r<1) throw std::runtime_error("Require positive sample count");
      p=model.initial();
      for(int k=1;k<=r;k++) {p=model.evolve(p,argv[8],scale);report(k*tau);}
    } else {
      if(mode=="exact") p=model.exact(argv[8],scale);
      else if(mode=="circuit" && r>0) p=model.circuit(tau,r,argv[9]);
      else throw std::runtime_error("Unknown mode or invalid steps");
      report(tau);
    }
    if(argc>10) {std::ofstream f(argv[10],std::ios::binary);f.write(reinterpret_cast<const char*>(p.data()),p.size()*sizeof(C));if(!f) throw std::runtime_error("State write failed");}
  } catch(const std::exception& e) {std::cerr<<e.what()<<"\n";return 1;}
}
