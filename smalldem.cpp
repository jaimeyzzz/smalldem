#include <cmath>
#include <cstdio>
#include <vector>
inline double clampRad(double a) { 
    return a<-M_PI?clampRad(a+2.0*M_PI):(a>M_PI?clampRad(a-2.0*M_PI):a);
}
struct Vec {
    double x,y;
    Vec(double x=0.0,double y=0.0) : x(x),y(y) {}
    Vec operator+(const Vec&a) const {return Vec(x+a.x,y+a.y); }
    Vec operator-(const Vec&a) const {return Vec(x-a.x,y-a.y); }
    Vec operator*(double a) const { return Vec(x*a,y*a); }
    Vec operator/(double a) const { return Vec(x/a,y/a); }
    double len() { return std::sqrt(x*x+y*y); }
    double dot(const Vec&a) { return x*a.x+y*a.y; }
};
struct Bond {
    int s,j;double l;Vec d;double tn,tt;
    Bond(int j,int s,double l,Vec d,double tn,double tt)
        :j(j),s(s),l(l),d(d),tn(tn),tt(tt){};
};
struct Particle {
    double im,iI,r,q,w,wMid,T;Vec x,v,vMid,F;int c;
    std::vector<Bond> bonds;
    Particle(double im,double iI,double r,Vec x,int c)
        :im(im),iI(iI),r(r),x(x),c(c),wMid(0),w(0),q(0){}
};
int main() {
    int numOfFrames=50;double fps=60,r=0.02,kn=1e7;
    std::vector<Particle> pts;
    for(double x=0.1;x<0.9;x+=2*r)for(double y=0.4;y<0.45;y+=2*r)
        pts.push_back(Particle(1.3e-4/r/r,2.6e-4/r/r/r/r,r,Vec(x,y),1));
    for(double x=0.7;x<0.9;x+=2*r)for(double y=0.1;y<0.3;y+=2*r)
        pts.push_back(Particle(1.3e-4/r/r,2.6e-4/r/r/r/r,r,Vec(x,y),2));
    for(int i=0;i<pts.size();i++)for(int j=0;j<pts.size();j++)if(i!=j){
        Vec l=(pts[j].x-pts[i].x);
        double overlap=pts[i].r+pts[j].r-l.len();
        if(overlap>=-0.1*r)
            pts[i].bonds.push_back(Bond(j,0,2.0*r,l/l.len(),.07*kn,.07*kn));
    }
    for(double x=r;x<1.0;x+=2*r){
        pts.push_back(Particle(0.0,0.0,r,Vec(x,0.0),0));
        pts.push_back(Particle(0.0,0.0,r,Vec(x,1.0),0));
        pts.push_back(Particle(0.0,0.0,r,Vec(0.0,x),0));
        pts.push_back(Particle(0.0,0.0,r,Vec(1.0,x),0));
    }
    for (double x=0.2;x<0.8;x+=2*r)for(double y=0.5;y<0.7;y+=2*r)
        pts.push_back(Particle(1.3e-4/r/r,2.6e-4/r/r/r/r,r,Vec(x,y),3));
    int s=(int)(1.0/fps/std::sqrt(7.5e3*r*r/kn))*10;double dt=1.0/fps/(double)s;
    FILE*fp=std::fopen("output.tex","w");
    fprintf(fp,"\\documentclass{article}\\usepackage{tikz}");
    fprintf(fp,"\\usepackage{caption}\n\\begin{document}\n");
    for (int frameIdx=0;frameIdx<numOfFrames;frameIdx++){
        for (int iter=0;iter<s;iter++){
            for(int i=0;i<pts.size();i++){
                pts[i].F=Vec();pts[i].T=0.0;
                pts[i].x=pts[i].x+pts[i].vMid*dt;
                pts[i].q=clampRad(pts[i].q+pts[i].wMid*dt);
            }
            for(int i=0;i<pts.size();i++)if(pts[i].im>1e-6){
                for(int j=0;j<pts.size();j++)if(i != j){
                    Vec lij=pts[i].x-pts[j].x;
                    double o=pts[i].r+pts[j].r-lij.len();
                    if(o<=1e-12)continue;
                    Vec n=lij/lij.len();
                    double a=1.4*std::sqrt(kn/(pts[i].im+pts[j].im));
                    pts[i].F=pts[i].F+n*(kn*o+a*(pts[j].v-pts[i].v).dot(n));
                }
                for(int iter=0;iter<pts[i].bonds.size();iter++){
                    Bond& b=pts[i].bonds[iter];if(b.s!=0)continue;
                    Vec l=pts[b.j].x-pts[i].x,n=l/l.len(),t(-n.y, n.x);
                    double dl=l.len()-b.l,qb=atan2(b.d.y,b.d.x)-atan2(n.y,n.x);
                    double ti=clampRad(qb+pts[i].q),tj=clampRad(qb+pts[b.j].q);
                    Vec Fn=n*kn*dl,Ft=t*-kn/3.0*r*r/l.len()*(ti + tj);
                    double T=kn/6.0*r*r*(tj-3.0*ti); 
                    if ((dl>0.0&&(Fn.len()/2.0/r+std::abs(kn/2.0*(tj-ti)))>b.tn)
                        ||Ft.len()/2.0/r>b.tt) {b.s=1;continue;}
                    pts[i].F=pts[i].F+Fn+Ft;pts[i].T=pts[i].T+T;
                }
            }
            for(int i=0; i<pts.size();i++){
                Vec acc=pts[i].F*pts[i].im+(pts[i].im>1e-6?Vec(0,-9.8):Vec());
                pts[i].v=pts[i].vMid+acc*0.5*dt;pts[i].vMid=pts[i].vMid+acc*dt;
                pts[i].w=pts[i].wMid+pts[i].T*pts[i].iI*0.5*dt;
                pts[i].wMid=pts[i].wMid+pts[i].T*pts[i].iI*dt;
            }
        }
        fprintf(fp,"\\begin{figure}\\begin{tikzpicture}[scale=12.5]\n");
        fprintf(fp,"\\draw[darkgray] (0,0)--(0,1)--(1,1)--(1,0)--(0,0);");
        const char*cols[]={"darkgray","olive","cyan","orange"};
        for(int i=0;i<pts.size();i ++){
            fprintf(fp,"\\draw[gray,fill=%s] (%f,",cols[pts[i].c],pts[i].x.x);
            fprintf(fp,"%f)circle [radius=%f];",pts[i].x.y,pts[i].r);
        }
        fprintf(fp,"\\end{tikzpicture}\n");printf("write frame %d\n",frameIdx);
        fprintf(fp,"\\caption{frame %d}\n\\end{figure}\n",frameIdx);
    }
    fprintf(fp,"\\end{document}");fclose(fp);
}