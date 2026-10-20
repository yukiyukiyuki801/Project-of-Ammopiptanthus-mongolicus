#include <Rcpp.h>
#include <geodesic.h>
#include <algorithm>
#include <array>
#include <vector>
#include <cmath>
using namespace Rcpp;
struct Box {std::array<double,7> lo,hi; int begin,end,left=-1,right=-1;};
class Search {
 public:
 NumericMatrix t,q; NumericVector knots,coeff; bool geo; std::vector<int> order; std::vector<Box> nodes; std::array<double,7> scale; struct geod_geodesic earth; double best,tol=1e-10,near; int winner,ties; long long evals=0; int qi;
 Search(NumericMatrix tt,NumericMatrix qq,NumericVector kk,NumericVector cc,bool gg):t(tt),q(qq),knots(kk),coeff(cc),geo(gg){geod_init(&earth,6378137,1/298.257223563);order.resize(t.nrow());for(int i=0;i<t.nrow();i++)order[i]=i;for(int d=0;d<7;d++){double lo=t(0,d),hi=lo;for(int i=1;i<t.nrow();i++){lo=std::min(lo,t(i,d));hi=std::max(hi,t(i,d));}scale[d]=std::max(hi-lo,1e-14);}build(0,t.nrow());}
 int build(int a,int b){Box z;z.begin=a;z.end=b;for(int d=0;d<7;d++){z.lo[d]=z.hi[d]=t(order[a],d);for(int k=a+1;k<b;k++){z.lo[d]=std::min(z.lo[d],t(order[k],d));z.hi[d]=std::max(z.hi[d],t(order[k],d));}}int id=nodes.size();nodes.push_back(z);if(b-a>24){int d=0;double spread=-1;for(int j=0;j<7;j++){if(!geo && j>=5)continue;double s=(z.hi[j]-z.lo[j])/scale[j];if(s>spread){spread=s;d=j;}}int mid=(a+b)/2;std::nth_element(order.begin()+a,order.begin()+mid,order.begin()+b,[&](int x,int y){return t(x,d)<t(y,d)||(t(x,d)==t(y,d)&&x<y);});int l=build(a,mid),r=build(mid,b);nodes[id].left=l;nodes[id].right=r;}return id;}
 double spline(double v,double a,double b,double c){if(v<=a)return 0;if(v>=c)return 1;if(v<b)return (v-a)*(v-a)/((b-a)*(c-a));return 1-(c-v)*(c-v)/((c-b)*(c-a));}
 double geographic(double d){if(!geo)return 0;return coeff[0]*spline(d,knots[0],knots[0],knots[1])+coeff[1]*spline(d,knots[0],knots[1],knots[2])+coeff[2]*spline(d,knots[1],knots[2],knots[2]);}
 double lower(int id){const Box&z=nodes[id];double v=0;for(int d=0;d<5;d++){double x=q(qi,d);v+=std::max(0.0,std::max(z.lo[d]-x,x-z.hi[d]));}if(geo){double x=std::max(0.0,std::max(z.lo[5]-q(qi,5),q(qi,5)-z.hi[5]));double y=std::max(0.0,std::max(z.lo[6]-q(qi,6),q(qi,6)-z.hi[6]));v+=geographic(std::hypot(x,y));}return -std::expm1(-v);}
 double score(int j){evals++;double v=0;for(int d=0;d<5;d++)v+=std::abs(q(qi,d)-t(j,d));if(geo)v+=geographic(std::hypot(q(qi,5)-t(j,5),q(qi,6)-t(j,6)));return -std::expm1(-v);}
 void walk(int id,bool tiepass){if(lower(id)>best+(tiepass?tol:0)+2e-14)return;const Box&z=nodes[id];if(z.left<0){for(int k=z.begin;k<z.end;k++){int j=order[k];double v=score(j);if(!tiepass){best=std::min(best,v);}else if(std::abs(v-best)<=tol){ties++;double meters;geod_inverse(&earth,q(qi,6),q(qi,5),t(j,6),t(j,5),&meters,NULL,NULL);double km=meters/1000;if(km<near-tol||(std::abs(km-near)<=tol&&(winner<0||j<winner))){near=km;winner=j;}}}return;}double a=lower(z.left),b=lower(z.right);if(a<=b){walk(z.left,tiepass);walk(z.right,tiepass);}else{walk(z.right,tiepass);walk(z.left,tiepass);}}
 NumericMatrix run(){NumericMatrix out(q.nrow(),5);for(qi=0;qi<q.nrow();qi++){if(qi%1000==0)checkUserInterrupt();best=1;near=R_PosInf;winner=-1;ties=0;long long before=evals;walk(0,false);walk(0,true);out(qi,0)=best;out(qi,1)=winner+1;out(qi,2)=near;out(qi,3)=ties;out(qi,4)=double(evals-before);}return out;}
};
// [[Rcpp::export]]
NumericMatrix exact_gdm_search(NumericMatrix targets,NumericMatrix queries,NumericVector geo_knots,NumericVector geo_coeff,bool geographic=true){if(targets.ncol()!=7||queries.ncol()!=7||geo_knots.size()!=3||geo_coeff.size()!=3)stop("Unexpected model dimensions");Search s(targets,queries,geo_knots,geo_coeff,geographic);return s.run();}
