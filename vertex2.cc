// C++ program for simulation of 3D-morphogenesis through nonuniform sheet growth regurated by morphogen intensity distribution over the sheet.
// Written by Hiroshi C. Ito (2024)
// Email: hiroshibeetle@gmail.com
//
// copyright (C) 2024 Hiroshi C. Ito
// This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License Version 2 as 
// published by the Free Software Foundation.
// http://www.r-project.org/Licenses/

//Compilation: g++ [this program file]
//Execution: ./a.out -o output

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <cstring>
#include <fstream>
using namespace std;

#define PI 3.14159265358979323846
const double Rand =  RAND_MAX;
double theta=(0.0)*PI;

/*_/_/_/_/_/_/_/_/_/ Model prameters _/_/_/_/_/_/_/_/*/

const int nn=129; /* number of vertices along the X- or Y-axis*/
const int nfid=2*3*(nn-1)*(nn-1); /* number of faces */
const int ncon=4*(nn-1)*(nn-1); /* number of springs */


double range_max=2.0; /* maximum value for X or Y corrdinates of vertices */
double range_min=-2.0; /* minimum value for X or Y corrdinates of vertices */

int tcool=600000; /* start time for cooling (end of spring extension) */
int TT  =1250000; /* total time period for simulation */
double develop_time=80; /* time for sheet growth (for development of growth metric) */
double dt=0.0008;  /* time step size */
int view_interval=1000; /* output interval */

double kk=10.0; /* spring constant during spring extension */
double kk_cool=10; /* spring constant during cooling */

double amp_cf=1.0*(17.0*17)/(25*25); /* bending elasticity */
double satu=0.0001; /* saturation constant for bending elesticity formula */


double upg0=3.5; /* water pressure during spring extension */
double upg_cool=0.03;  /* water pressure during cooling */
double shift_upg=50000; /* steepness of switching from upg0 to upg_cool */
double upg=upg0;

double upgg=-8.0*5; /* strength and direction of the force on corners of the sheet */


double d_gradx_cen=0.0; /* curvature for central horn along the X-axis  */
double d_grady_cen=-0.15*0.4; /* curvature for central horn along the Y-axis  */
double amp_tw_cen=0.0; /* degree of twist for central horn */


double d_gradx_side=-0.0; /* curvature for side horns along the X-axis  */
double d_grady_side=0.05; /* curvature for side horns along the X-axis  */
double amp_tw_side=0.0125; /* degree of twist for side horns */


double amp_curve_shape=1.0; /* adjustment parameter for horn curvatures */
double amp_horn=0.05; /* adjustment parameter for horn sizes */

double ampxy_init=1.0; /* initial sheet extension along the x- and y-axes */


int flag_cf_along=1;
int flag_edge=1;
double eps=1e-5;




double x[nn][nn],y[nn][nn],z[nn][nn],w0[nn][nn],w1[nn][nn];
double v[nn*nn][3],v0[nn*nn][3],F[nn*nn][3],vnorm[nn*nn][3],vb[nn][nn][3],vnormb[nn][nn][3];
double vl[ncon][3],vu[ncon][3],Fcon[ncon][3];
double flen[ncon],vlen[ncon];
double cf[nn][nn][3];

int edgemask[nn*nn],ancmask[nn*nn], vid[nn][nn];
double met[nn*nn][2];
double metric[nn*nn][3];
int fid0[nn-1][nn-1][3], fid1[nn-1][nn-1][3];
int fid[2*3*(nn-1)*(nn-1)];

int con[2][ncon];
double len0[ncon],len00[ncon],len1[ncon],len[ncon],dlen[ncon],nabla2[ncon];

const int nface=(nn-1)*(nn-1)*2;
int nv[nface];
double vv[nfid][3];

double dd;

double dxe[nn-1][nn][3],dxn[nn-2][nn][3];
double dye[nn][nn-1][3],dyn[nn][nn-2][3];

double dxn0[nn-2][nn][3],dxn1[nn-2][nn][3];
double dyn0[nn][nn-2][3],dyn1[nn][nn-2][3];
double dxlenx[nn-1][nn];
double dyleny[nn][nn-1];


/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
/*_/_/ Primitive functions for vertices and faces                 _/_/*/
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/


inline void cerr3(double a[3]){
  cerr<<a[0]<<" "<<a[1]<<" "<<a[2]<<endl;
}
inline void set3(double a[3],double b[3]){
  a[0]=b[0];
  a[1]=b[1];
  a[2]=b[2];
}
inline void set3(double a[3],double b){
  a[0]=b;
  a[1]=b;
  a[2]=b;
}

inline void scale3(double a[3],double b[3],double c){
  a[0]=b[0]*c;
  a[1]=b[1]*c;
  a[2]=b[2]*c;
}

inline int id(int i,int j){
  return i+nn*j;
}
inline void add3(double a[3],double b[3],double c[3]){
  a[0]=b[0]+c[0];
  a[1]=b[1]+c[1];
  a[2]=b[2]+c[2];
}

inline void add3q(double a[3],double b[3],double c[3]){
  a[0]+=b[0]+c[0];
  a[1]+=b[1]+c[1];
  a[2]+=b[2]+c[2];
}
inline void add3(double a[3],double b[3],double c){
  a[0]=b[0]+c;
  a[1]=b[1]+c;
  a[2]=b[2]+c;
}


inline void sub3(double a[3],double b[3],double c[3]){
  a[0]=b[0]-c[0];
  a[1]=b[1]-c[1];
  a[2]=b[2]-c[2];
}

inline void mult3(double a[3],double b[3],double c[3]){
  a[0]=b[0]*c[0];
  a[1]=b[1]*c[1];
  a[2]=b[2]*c[2];
}


inline double get_length(double xv[3]){
  return sqrt(xv[0]*xv[0]+xv[1]*xv[1]+xv[2]*xv[2]);
}

void set_vert_id(int vid[nn][nn]){
  int count=0;
  for(int j=0;j<nn;j++){
	for(int i=0;i<nn;i++){
	  vid[i][j]=count;
	  count++;
	}
  }
  
}

void set_face_id(void){
  int count=0;
  for(int i=0;i<nn-1;i++){
    for(int j=0;j<nn-1;j++){
      if(count%2==0){
	fid0[j][i][0]=vid[j][i];
	fid0[j][i][1]=vid[j+1][i];
	fid0[j][i][2]=vid[j][i+1];
      }else{
	fid0[j][i][0]=vid[j][i];
	fid0[j][i][1]=vid[j+1][i];
	fid0[j][i][2]=vid[j+1][i+1];
      }
      
      count++;
    }
  }

  count=0;
  for(int i=0;i<nn-1;i++){
    for(int j=0;j<nn-1;j++){
      if(count%2==0){
	fid1[j][i][0]=vid[j+1][i];
	fid1[j][i][1]=vid[j+1][i+1];
	fid1[j][i][2]=vid[j][i+1];
      }else{
	fid1[j][i][0]=vid[j][i];
	fid1[j][i][1]=vid[j+1][i+1];
	fid1[j][i][2]=vid[j][i+1];
      }
      count++;
    }
  }

   count=0;
  for(int i=0;i<nn-1;i++){
    for(int j=0;j<nn-1;j++){

	con[0][count]=i+nn*j;
	con[1][count]=(i+1)+nn*j;
	count++;
	con[0][count]=i+nn*j;
	con[1][count]=i+nn*(j+1);
	count++;
	con[0][count]=i+nn*j;
	con[1][count]=(i+1)+nn*(j+1);
	count++;
	con[0][count]=(i+1)+nn*j;
	con[1][count]=i+nn*(j+1);
	count++;

	

    }
  }
  
  
}

void set_v(double v[nn*nn][3], double x[nn][nn],double y[nn][nn], double z[nn][nn]){
  int count=0;
    for(int j=0; j<nn;j++){
      for(int i=0; i<nn;i++){
	v[count][0]=x[i][j];
	v[count][1]=y[i][j];
	v[count][2]=z[i][j];
	count++;
      }
    }
}



void set_vv(double vv[nfid][3],double v[nn*nn][3]){
  for(int i=0; i<nfid;i++){
    vv[i][0]=v[fid[i]][0];
    vv[i][1]=v[fid[i]][1];
    vv[i][2]=v[fid[i]][2];
  }
}

void out_nv_vv(std::ofstream oout){
  oout<<nface<<endl;
  for(int i=0;i<nface;i++){
    oout<<nv[i]<<" ";
    if((i+1)%nn==0)oout<<endl;
  }
  oout<<endl<<endl;
  
  oout<<nfid<<endl;

  
  for(int i=0;i<nfid;i++){
    oout<<vv[i][0]<<" "<<vv[i][1]<<" "<<vv[i][2]<<endl;
  }


}




void out_nv_vv_norm(std::ofstream oout){

  oout<<nface<<endl;
  for(int i=0;i<nface;i++){
    oout<<nv[i]<<" ";
    if((i+1)%nn==0)oout<<endl;
  }
  oout<<endl<<endl;
  
  oout<<nfid<<endl;

  
  for(int i=0;i<nfid;i++){
    oout<<vv[i][0]<<" "<<vv[i][1]<<" "<<vv[i][2]<<endl;
  }


  oout<<endl<<endl;
  
  oout<<nn*nn<<endl;

 
   for(int i=0;i<nn*nn;i++){
    oout<<v[i][0]<<" "<<v[i][1]<<" "<<v[i][2]<<endl;
  }

   for(int j=0;j<nn;j++){
   for(int i=0;i<nn;i++){
     for(int k=0;k<3;k++)oout<<cf[i][j][k]<<" ";
     oout<<endl;
   }
   }

}
  
double expf(double x,double y,double sx,double sy){
  return exp(-(x*x)/(2.0*sx*sx))*exp(-(y*y)/(2.0*sy*sy));
}

void expf_v(double z[nn][nn], double x[nn][nn],double y[nn][nn],double sx,double sy){
    for(int j=0;j<nn;j++){
      for(int i=0;i<nn;i++){
	z[i][j]=expf(x[i][j],y[i][j],sx,sy);
      }
    }
}


double exp_area(double x, double y, double dx,double dy, double rx, double ry,double sd){
  double dis;
  dis=(x-dx)*(x-dx)/(rx*rx)+(y-dy)*(y-dy)/(ry*ry);
  if(dis<=1){
    return 1;
  }else{
    return 
      exp(-0.5*(sqrt(dis)-1)*(sqrt(dis)-1)/(sd*sd));
  }
}


double expfa(double x,double y,double a0,double a1,double dx0,double dx1,double dy0,double dy1,double sx0,double sx1,double sy0,double sy1, int ndiv){
      double value=0.0;
      double aa,ddx,ddy,ssx,ssy;
      for(int k=0;k<=ndiv;k++){
	aa=k*(a1-a0)/double(ndiv) + a0;
	ddx=k*(dx1-dx0)/double(ndiv) + dx0;
	ddy=k*(dy1-dy0)/double(ndiv) + dy0;
	ssx=k*(sx1-sx0)/double(ndiv) + sx0;
	ssy=k*(sy1-sy0)/double(ndiv) + sy0;
	
	value+=aa*expf(x-ddx,y-ddy,ssx,ssy);
      }
      return value;
}


void flat2mesh(double vb[nn][nn][3],double v[nn*nn][3]){
  int count=0;
 
  for(int j=0;j<nn;j++){
    for(int i=0;i<nn;i++){
      vb[i][j][0]=v[count][0];
      vb[i][j][1]=v[count][1];
      vb[i][j][2]=v[count][2];
      }
      count++;
    }
}


void mesh2flat(double v[nn*nn][3],double vb[nn][nn][3]){
  int count=0;
  for(int j=0;j<nn;j++){
    for(int i=0;i<nn;i++){
      	v[count][0]=vb[i][j][0];
	v[count][1]=vb[i][j][1];
	v[count][2]=vb[i][j][2];
	count++;
    }
  }
  
}

double get_area(double dv2[3], double dv1[3]){
  double norm[3];
  norm[0]=dv1[1]*dv2[2]-dv1[2]*dv2[1];
  norm[1]=dv1[2]*dv2[0]-dv1[0]*dv2[2];
  norm[2]=dv1[0]*dv2[1]-dv1[1]*dv2[0];
  return get_length(norm);
  
}
void face_norm_tri(double norm[3],double v0[3],double v1[3], double v2[3]){
  double dv1[3],dv2[3];
  double leng;

  sub3(dv1,v1,v0);
  sub3(dv2,v2,v0);
  
  
    norm[0]=dv1[1]*dv2[2]-dv1[2]*dv2[1];
    norm[1]=dv1[2]*dv2[0]-dv1[0]*dv2[2];
    norm[2]=dv1[0]*dv2[1]-dv1[1]*dv2[0];
  
  
  leng=get_length(norm);

  norm[0]= -1*norm[0]/leng;
  norm[1]= -1*norm[1]/leng;
  norm[2]= -1*norm[2]/leng;
}


bool inside(int i, int j){
  if((-1<i)&&(i<nn)&&(-1<j)&&(j<nn)){
    return true;
  }else{
    return false;
  }
  
}

void vert_normal(void){
  double normals[4][3], v0[3],v1[3],v2[3],vn[3];
  double vnlen;
  int count,vcount;
  
  count=0;  
  for(int j=0;j<nn;j++){
    for(int i=0;i<nn;i++){
      vcount=0;
      if(inside(i,j-1)&&inside(i-1,j)){
	
	face_norm_tri(normals[vcount],v[i+j*nn],v[i+(j-1)*nn],v[i-1+j*nn]);

	vcount++;
	
      }

      if(inside(i+1,j)&&inside(i,j-1)){
	face_norm_tri(normals[vcount],v[i+j*nn],v[i+1+j*nn],v[i+(j-1)*nn]);
	vcount++;
      }
      if(inside(i,j+1)&&inside(i+1,j)){
	face_norm_tri(normals[vcount],v[i+j*nn],v[i+(j+1)*nn],v[i+1+j*nn]);
	vcount++;
      }
      if(inside(i-1,j)&&inside(i,j+1)){
	face_norm_tri(normals[vcount],v[i+j*nn],v[i-1+j*nn],v[i+(j+1)*nn]);
	vcount++;
      }

      vn[0]=0;vn[1]=0;vn[2]=0;
      for(int k=0;k<vcount;k++){
	for(int l=0;l<3;l++)vn[l]+=normals[k][l];
	
      }


      
      vnlen=get_length(vn);
      
      
      for(int k=0;k<3;k++){
	vn[k]=vn[k]/vnlen;
	vnorm[count][k]=vn[k];
      }
      count++;
    }
  }

}



double in_prod(double a[3], double b[3]){
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

void out_mat(double z[nn][nn], std::ofstream oout){
  for(int j=0;j<nn;j++){
    for(int i=0;i<nn;i++){
      oout<<z[i][j]<<" ";
    }
    oout<<endl;
  }
}
void out_vert(double x[nn][nn],double y[nn][nn],double z[nn][nn], std::ofstream oout){
  for(int j=0;j<nn;j++){
    for(int i=0;i<nn;i++){
      oout<<x[i][j]<<" "<<y[i][j]<<" "<<z[i][j]<<endl;
    }
    oout<<endl;
  }
}

void set_len0(){
  double dx, dy, dz;
  for(int i=0;i<ncon;i++){
    dx=v[con[1][i]][0]-v[con[0][i]][0];
    dy=v[con[1][i]][1]-v[con[0][i]][1];
    dz=v[con[1][i]][2]-v[con[0][i]][2];
    			   
    len0[i]=sqrt(dx*dx+dy*dy+dz*dz);
    len00[i]=sqrt(dx*dx+dy*dy);
  }
}



inline void get_ev(double a[3],double b[3],double c[3]){
  double le;
  a[0]=b[0]-c[0];
  a[1]=b[1]-c[1];
  a[2]=b[2]-c[2];
  le=get_length(a);
  a[0]/=le;
  a[1]/=le;
  a[2]/=le;
  
  
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
/*_/_/ Functions for calculation of forces on vertices            _/_/*/
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/



void calc_curv_force(void){
  double lenbuf;
  double vbuf[3],vbuf1[3],vbuf2[3];
  double dxe0[3],dye0[3],dxe1[3],dye1[3],dxe0_dxe1,dye0_dye1;
 
  for(int j=0;j<nn;j++){
    for(int i=0;i<nn-1;i++){
      sub3(vbuf,v[i+1+j*nn],v[i+j*nn]);
      lenbuf=get_length(vbuf);
      scale3(vbuf,vbuf,1.0/lenbuf);
      
      set3(dxe[i][j],vbuf);
      dxlenx[i][j]=lenbuf;
    }
  }
  
  for(int j=0;j<nn;j++){
    for(int i=0;i<nn-2;i++)sub3(dxn[i][j],dxe[i+1][j],dxe[i][j]);
  }
  
  for(int j=0;j<nn-1;j++){
    for(int i=0;i<nn;i++){
      sub3(vbuf,v[i+(j+1)*nn],v[i+j*nn]);
      lenbuf=get_length(vbuf);
      //scale3(vbuf,vbuf,lenbuf);
      scale3(vbuf,vbuf,1.0/lenbuf);
      
      set3(dye[i][j],vbuf);
      dyleny[i][j]=lenbuf;
    }
  }

  for(int j=0;j<nn-2;j++){
    for(int i=0;i<nn;i++)sub3(dyn[i][j],dye[i][j+1],dye[i][j]);
  }
  

  for(int j=0;j<nn;j++){
    for(int i=0;i<nn-2;i++){
      set3(dxe0,dxe[i][j]);
      set3(dxe1,dxe[i+1][j]);

      dxe0_dxe1=in_prod(dxe0,dxe1);
      for(int k=0;k<3;k++){
      dxn0[i][j][k]=(0.5/(1+dxe0_dxe1+satu))*(dxe1[k]-1.0*dxe0_dxe1*dxe0[k]);
      dxn1[i][j][k]=(0.5/(1+dxe0_dxe1+satu))*(dxe0_dxe1*dxe1[k]-dxe0[k]);

      }
    }
  }


  for(int j=0;j<nn-2;j++){
    for(int i=0;i<nn;i++){
      set3(dye0,dye[i][j]);
      set3(dye1,dye[i][j+1]);
      
      dye0_dye1=in_prod(dye0,dye1);
      for(int k=0;k<3;k++){
	dyn0[i][j][k]=(0.5/(1+dye0_dye1+satu))*(dye1[k]-1.0*dye0_dye1*dye0[k]);
	dyn1[i][j][k]=(0.5/(1+dye0_dye1+satu))*(dye0_dye1*dye1[k]-dye0[k]);
    }
  }
  }


  
  for(int j=0;j<nn;j++)for(int i=0;i<nn;i++)set3(cf[i][j],0.0);

  for(int j=0;j<nn;j++){
    for(int i=0;i<nn-2;i++){
      for(int k=0;k<3;k++){
	vbuf1[k]=0.5*dxn0[i][j][k]/dxlenx[i][j];
	vbuf2[k]=0.5*dxn1[i][j][k]/dxlenx[i+1][j];
	
	cf[i+1][j][k]+=vbuf1[k]+vbuf2[k];
	cf[i][j][k]-=vbuf1[k];
	cf[i+2][j][k]-=vbuf2[k];
      }
    }
  }

   for(int j=0;j<nn-2;j++){
    for(int i=0;i<nn;i++){
      for(int k=0;k<3;k++){
	vbuf1[k]=0.5*dyn0[i][j][k]/dyleny[i][j];
	vbuf2[k]=0.5*dyn1[i][j][k]/dyleny[i][j+1];
	
	cf[i][j+1][k]+=vbuf1[k]+vbuf2[k];
	cf[i][j][k]-=vbuf1[k];
	cf[i][j+2][k]-=vbuf2[k];
     
      }
    }
   }


  
}

void force(double len[]){
  double cflen,disx,disy,bufvx1[3],bufvx2[3],bufvy1[3],bufvy2[3];
  
  for(int i=0;i<nn*nn;i++){
    set3(F[i],0.0);
  }
  
  for(int i=0;i<ncon;i++){
    sub3(vl[i],v[con[1][i]],v[con[0][i]]);
    vlen[i]=get_length(vl[i]);
    
    scale3(vu[i],vl[i],1.0/vlen[i]);

  }
  
  for(int i=0;i<ncon;i++){
        flen[i]=-kk*(vlen[i]-len[i])/len[i];
    
    scale3(Fcon[i],vu[i],flen[i]);

   
  }

  
  for(int i=0;i<ncon;i++){
    for(int k=0;k<3;k++){
       F[con[1][i]][k]+=Fcon[i][k];    
            F[con[0][i]][k]-=Fcon[i][k];
    }
  }
  
  
  for(int i=0;i<nn*nn;i++){
    if(edgemask[i]==1)F[i][2]+=upgg*dd*dd;

  }
  
  vert_normal();

  for(int j=1;j<nn-1;j++){
    for(int i=1;i<nn-1;i++){
         sub3(bufvx1,v[(i+1)+nn*j],v[(i-1)+nn*j]);
	 sub3(bufvy1,v[i+nn*(j+1)],v[i+nn*(j-1)]);
    
	 	 for(int k=0;k<3;k++)F[i+nn*j][k]+=vnorm[i+nn*j][k]*upg*(get_area(bufvx1,bufvy1));


    }

  }
  
  calc_curv_force();
  for(int j=0;j<nn;j++){
    for(int i=0;i<nn;i++){
      if(flag_cf_along==1){
	for(int k=0;k<3;k++)F[i+j*nn][k]+=amp_cf*cf[i][j][k];
      }else{
      cflen=in_prod(cf[i][j],vnorm[i+nn*j]);
      for(int k=0;k<3;k++)F[i+j*nn][k]+=amp_cf*cflen*vnorm[i+nn*j][k];
      }
    }
  }
  
}



/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
/*_/_/ Functions for sheet growth regulated by morphogen distribution _/_/*/
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

double sx0=0.5;
double sy0=0.5;
double dx0=0.2;
double dy0=0.5;
double a11=0.1;
double dx11=0.2;
double dy11=-1.0;
double sx11=0.5;
double sy11=1.0;


inline void rotate(double &x0, double &y0){
  double x,y;
  x=x0*cos(theta)+y0*sin(theta);
  y=-1*x0*sin(theta)+y0*cos(theta);
  x0=x;
  y0=y;
 
}

double morpha_head(double x, double y){
  rotate(x,y);  
  
  double dx111=0.2*7;
  double dy111=1.2;
  double value=0;
  
  // value=2*0.7*expfa(x,y,0.3,0.2,0,0,0.5*dy0,0.5*dy0,sx0,0.07,sy0,0.07,20);
  //value=3*0.7*expfa(x,y,0.3,0.2,0,0,0.6*dy0,0.6*dy0,sx0*0.7,0.07,sy0*0.8,0.07,20);
  value=3.5*0.7*expfa(x,y,0.3,0.2,0,0,0.6*dy0,0.6*dy0,sx0*0.7,0.07,sy0*0.8,0.07,20);
    
  return value;
}
double morpha(double x, double y){  
  rotate(x,y);  
  
  double dx111=0.2*7;
  double dy111=1.2;
  double value=0;
  value+=6*0.8*exp(-x*x/(2*1.3*1.3)-(y-1.0)*(y-1.0)/(2*0.5*0.5));
  value+=6*0.8*exp(-x*x/(2*0.7*0.7)-(y-0.3)*(y-0.3)/(2*1*1));

  value+=3*1.5*exp(-x*x/(2*1*1)-(y+1.7)*(y+1.7)/(2*1*1));
  

  value+=0.5*exp(-(x-dx111)*(x-dx111)/(2*0.7*0.7)-(y+dy111)*(y+dy111)/(2*0.7*0.7));
  value+=0.5*exp(-(x+dx111)*(x+dx111)/(2*0.7*0.7)-(y+dy111)*(y+dy111)/(2*0.7*0.7));
  
  value+=1.5*exp(-(x-dx111)*(x-dx111)/(2*0.4*0.4)-(y+dy111)*(y+dy111)/(2*0.4*0.4));
  value+=1.5*exp(-(x+dx111)*(x+dx111)/(2*0.4*0.4)-(y+dy111)*(y+dy111)/(2*0.4*0.4));
 
   value+=0.3*expfa(x,y,0.3,0.2, 0.0, 0.0,-0.5,-0.5,0.5*sx0,0.07,0.5*sy0,0.07,20);
   
  return value;
}

void set_init(void){
  dd=(range_max-range_min)/double(nn-1);
  for(int j=0;j<nn;j++){
    for(int i=0;i<nn;i++){
      x[i][j]=i*dd+range_min;
      y[i][j]=j*dd+range_min;
      z[i][j]=0.0;
      v0[i+nn*j][0]=x[i][j];
      v0[i+nn*j][1]=y[i][j];
      v0[i+nn*j][2]=0;
      
     
    }
  }


  set_vert_id(vid);
  set_face_id();

  for(int i=0;i<(nn-1)*(nn-1)*2;i++){
    nv[i]=3;
  }

  
  
    int count=0;
  for(int j=0;j<nn;j++){
    for(int i=0;i<nn;i++){
      edgemask[count]=0;
      if(flag_edge==1){
	if((i==0)||(i==(nn-1))||(j==0)||(j==(nn-1)))edgemask[count]=1;
	count++;
      }
      if(flag_edge==2){
      if(((i==0)||(i==(nn-1)))&&((j==0)||(j==(nn-1))))edgemask[count]=1;
      count++;
      }
      
    }
  }
}



void develop(double nabla2[],double (*morphogen)(double,double),double d_gradx,double d_grady, double amp_tw){
  double gx,gy,gradxy[2],gradx,grady,prod,nabla,len00,gradex,gradey,gradlen,dgradlen,dgex,dgey,cv;
  double vlxy[3],vlxy_tw[3],bufv0[3],bufv1[3],bufv[3];
  
  for(int i=0;i<ncon;i++){
    set3(bufv1,v0[con[1][i]]);
    set3(bufv0,v0[con[0][i]]);
    sub3(bufv,bufv1,bufv0);
    len00=get_length(bufv);
    get_ev(vlxy,bufv1,bufv0);
  
    vlxy_tw[0]=vlxy[0]*cos(amp_tw)+vlxy[1]*sin(amp_tw);
    vlxy_tw[1]=-1*vlxy[0]*sin(amp_tw)+vlxy[1]*cos(amp_tw);
 
    gx=0.5*(bufv0[0]+bufv1[0]);
    gy=0.5*(bufv0[1]+bufv1[1]);
    gradx=amp_horn*(morphogen(gx+0.5*eps,gy)-morphogen(gx-0.5*eps,gy))/eps;
    grady=amp_horn*(morphogen(gx,gy+0.5*eps)-morphogen(gx,gy-0.5*eps))/eps;
    gradxy[0]=gradx;
    gradxy[1]=grady;
    gradlen=sqrt(gradx*gradx+grady*grady);
    gradex=gradx/gradlen;
    gradey=grady/gradlen;

    dgradlen=sqrt(d_gradx*d_gradx+d_grady*d_grady);
    dgex=d_gradx/dgradlen;
    dgey=d_grady/dgradlen;

    cv=1+dgex*gradex+dgey*gradey;
    
    prod=vlxy_tw[0]*gradxy[0]+vlxy_tw[1]*gradxy[1];


    nabla=prod*(1+dgradlen*exp(amp_curve_shape*log(cv+1e-10)));
    //nabla=prod*(1+dgradlen*cv);
    
    nabla2[i]+=nabla*nabla;
  }
}

void calc_met(double (*morphogen)(double,double),double d_gradx,double d_grady, double amp_tw){
  double gx,gy,gradxy[2],gradx,grady,prod,nabla,len00,dgx;
  double vlxy[3],vlxy_tw[3],bufv0[3],bufv1[3],bufv[3];
  double cv,cvb,gradlen,gradex,gradey,dgradlen,dgex,dgey;
  
  for(int i=0;i<nn;i++){
      for(int j=0;j<nn;j++){
    gx=x[i][j];
    gy=y[i][j];
    dgx=1.0;


    gradx=amp_horn*(morphogen(gx+0.5*eps,gy)-morphogen(gx-0.5*eps,gy))/eps;
    grady=amp_horn*(morphogen(gx,gy+0.5*eps)-morphogen(gx,gy-0.5*eps))/eps;
    /*
     gradxy[0]=gradx*exp((dgx*d_gradx*gradx+d_grady*grady));
    gradxy[1]=grady*exp((dgx*d_gradx*gradx+d_grady*grady));
    */
    gradlen=sqrt(gradx*gradx+grady*grady);
    gradex=gradx/gradlen;
    gradey=grady/gradlen;

    dgradlen=sqrt(d_gradx*d_gradx+d_grady*d_grady);
    dgex=d_gradx/dgradlen;
    dgey=d_grady/dgradlen;

    cv=1+dgex*gradex+dgey*gradey;
    //cvb=(1+dgradlen*cv);
    cvb=(1+dgradlen*exp(amp_curve_shape*log(cv+1e-10)));
    gradxy[0]=gradx*cvb;
    gradxy[1]=grady*cvb;
    
    
    vlxy_tw[0]=gradxy[0]*cos(amp_tw)+gradxy[1]*sin(amp_tw);
    vlxy_tw[1]=-1*gradxy[0]*sin(amp_tw)+gradxy[1]*cos(amp_tw);
  

    met[i+j*nn][0]+=vlxy_tw[0];
      met[i+j*nn][1]+=vlxy_tw[1];

        metric[i+j*nn][0]+=vlxy_tw[0]*vlxy_tw[0];
    metric[i+j*nn][1]+=vlxy_tw[0]*vlxy_tw[1];
    metric[i+j*nn][2]+=vlxy_tw[1]*vlxy_tw[1];
    
  }
}
}



double morpha_deko(double x,double y){
    rotate(x,y);  

    double value=0;
  
      value+=0.8*exp(-x*x/(2*0.4*0.4)-(y-1.0)*(y-1.0)/(2*0.4*0.4))+2.0*exp(-x*x/(2*1*1)-(y+1.7)*(y+1.7)/(2*1*1))+0.8*exp(-(x-7.5*dx11)*(x-7.5*dx11)/(2*1*1)-(y+1.7)*(y+1.7)/(2*1*1))+0.8*exp(-(x+7.5*dx11)*(x+7.5*dx11)/(2*1*1)-(y+1.7)*(y+1.7)/(2*1*1));

    value+=0.3*expfa(x,y,0.3,0.2, 0.7*dx11, 0.7*dx11,dy11,dy11,sx0,sx0,sy0,sy0,20);
    value+=0.3*expfa(x,y,0.3,0.2,-0.7*dx11,-0.7*dx11,dy11,dy11,sx0,sx0,sy0,sy0,20);
    
    value+=0.3*expfa(x,y,0.3,0.2, 3*dx11, 3*dx11,dy11*0.9,dy11*0.9,sx0,0.1,sy0,0.1,20);
    value+=0.3*expfa(x,y,0.3,0.2,-3*dx11,-3*dx11,dy11*0.9,dy11*0.9,sx0,0.1,sy0,0.1,20);
     
    value+=0.7*expfa(x,y,0.3,0.2, 6*dx11, 6*dx11,dy11*0.8,dy11*0.8,sx0,0.07,sy0,0.07,20);
    value+=0.7*expfa(x,y,0.3,0.2,-6*dx11,-6*dx11,dy11*0.8,dy11*0.8,sx0,0.07,sy0,0.07,20);
    
    value+=0.4*expfa(x,y,0.3,0.2, 6*dx11, 8*dx11,dy11*0.7,dy11*0.7,sx0,0.07,sy0,0.07,20);
    value+=0.4*expfa(x,y,0.3,0.2,-6*dx11,-8*dx11,dy11*0.7,dy11*0.7,sx0,0.07,sy0,0.07,20);

  return value;

}

double morpha_sideL(double x,double y){
    rotate(x,y);  
    
    double value=0.0;
    value+=3.0*expfa(x,y,0.3,0.2,-6*dx11,-6*dx11,dy11*0.8,dy11*0.8,0.3,0.07,0.3,0.07,20);
    return value;
}

double morpha_sideR(double x,double y){
    rotate(x,y);  
    double value=0.0;
    value+=3.0*expfa(x,y,0.3,0.2,6*dx11,6*dx11,dy11*0.8,dy11*0.8,0.3,0.07,0.3,0.07,20);
    
    return value;
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
/*_/_/ Functions for output                                       _/_/*/
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/


void out_v_con_len_file(void){
       ofstream ofi("out_v_len.dat");

       ofi<<nn<<endl;


  
       for(int j=0;j<nn;j++){
	 for(int i=0;i<nn;i++){
	   ofi<<v[i+nn*j][0]<<" "<<v[i+nn*j][1]<<" "<<v[i+nn*j][2]<<endl;
	 }
	 ofi<<endl;
       }
       ofi<<endl;


       ofi<<ncon<<endl;
       for(int i=0;i<ncon;i++)ofi<<con[0][i]<<endl;
       ofi<<endl;
       for(int i=0;i<ncon;i++)ofi<<con[1][i]<<endl;
       ofi<<endl;
       for(int i=0;i<ncon;i++)ofi<<vlen[i]<<endl;
       ofi<<endl;

       for(int i=0;i<ncon;i++)ofi<<len[i]<<endl;
       ofi<<endl;
	

}

void out_met_file(std::string met_filename){
       ofstream ofi(met_filename);
       ofi<<nn<<endl;
       for(int j=0;j<nn;j++){
	 for(int i=0;i<nn;i++){
	   ofi<<metric[i+nn*j][0]<<" "<<metric[i+nn*j][1]<<" "<<metric[i+nn*j][2]<<endl;
	 }
	 ofi<<endl;
       }
       ofi<<endl;

}


/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
/*_/_/    Main function     _/_/*/
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

double tim;
int outcount;

main(int argc, char *argv[]){
int opt;
double buf;
   std::string ofile_name="gomi";
  while((opt=getopt(argc, argv,
		   "a:b:c:r:pe:f:g:h:i:j:k:l:m:n:o:")) != EOF){
   switch(opt){
   case 'o':
    ofile_name=optarg;
    cerr<<ofile_name<<endl;
     break;
   }
 }
 
   std::string outcount_filename=ofile_name+"_outcount.dat";
   std::string met_filename=ofile_name+"_met.dat";
   std::string out_filename=ofile_name+".dat";

   cerr<<out_filename<<endl;
   cerr<<outcount_filename<<endl;
   cerr<<met_filename<<endl;
   
   ofstream oout(out_filename);
   
  outcount=0;

  
  tim=0;
  set_init();

  for(int i=0;i<ncon;i++)nabla2[i]=0.0;
  
  develop(nabla2,morpha,d_gradx_cen,d_grady_cen,amp_tw_cen);
  develop(nabla2,morpha_deko,-0*d_gradx_side,d_grady_side,-0*amp_tw_side);
  develop(nabla2,morpha_sideL,0.15,0,-amp_tw_side);
  develop(nabla2,morpha_sideR,-0.15,0,amp_tw_side);
  develop(nabla2,morpha_head,d_gradx_cen,1.5*d_grady_cen,amp_tw_cen);
   
 for(int i=0;i<nn*nn;i++){
   met[i][0]=0.0;
   met[i][1]=0.0;
   metric[i][0]=0.0;
   metric[i][1]=0.0;
   metric[i][2]=0.0;
 }
 
 calc_met(morpha,d_gradx_cen,d_grady_cen,amp_tw_cen);
 calc_met(morpha_deko,-0*d_gradx_side,d_grady_side,-0*amp_tw_side);
 calc_met(morpha_sideL,0.15,0,-amp_tw_side);
 calc_met(morpha_sideR,-0.15,0,amp_tw_side);  
 calc_met(morpha_head,d_gradx_cen,1.5*d_grady_cen,amp_tw_cen);

 for(int i=0;i<nn*nn;i++){
   metric[i][0]=1.0+develop_time*metric[i][0];
   metric[i][1]=develop_time*metric[i][1];
   metric[i][2]=1.0+develop_time*metric[i][2];
 }
 

  /* initial setting for sheet extension */
  set_v(v,x,y,w0);
  set_len0();

 
  for(int i=0;i<nn*nn;i++){
    v[i][0]*=ampxy_init;
    v[i][1]*=ampxy_init;
    
  }
 set_vv(vv,v);
 
 for(int i=0;i<ncon;i++){
   dlen[i]+=sqrt(nabla2[i]*develop_time+1);
   len[i]=dlen[i]*len00[i];
 }

 for(int i=0;i<ncon;i++)len1[i]=len0[i];
 

   oout<<nn<<endl;

    /* sheet extension */
   
   for(int t=0;t<=TT;t++){
     upg=(upg0-upg_cool)/(1+exp((t-tcool)/shift_upg))+upg_cool;
     if(t==tcool){
       cerr<<"cooling!"<<endl;
          kk=kk_cool;
     }
   
    if(t<=1.0*tcool){
      for(int i=0;i<ncon;i++){
		len1[i]=log(len0[i])+(log(len[i])-log(len0[i]))*(double(t)/double(1.0*tcool));
	len1[i]=exp(len1[i]);
      }
    }
    
     
    force(len1);
    
    double gF[3]={0,0,0};
    for(int i=0;i<nn*nn;i++){
      for(int k=0;k<2;k++)gF[k]+=F[i][k];
    }

    for(int k=0;k<3;k++)gF[k]/=double(nn*nn);
    
    
    for(int i=0;i<nn*nn;i++){
      for(int k=0;k<3;k++)v[i][k]+=(F[i][k]-gF[k])*dt;
      
      if((edgemask[i]==1)){
	if(v[i][2]>0.0)v[i][2]=0.0;
      }
    }
    tim+=dt;

    if(t%view_interval==0){
    set_vv(vv,v);
    cerr<<t<<" "<<int(100*t/tcool)<<"%  total:"<<int(100*t/TT)<<"% "<<endl;
    cerr<<"gF:"<<gF[0]<<" "<<gF[1]<<" "<<gF[2]<<endl;
    oout<<t<<" "<<tcool<<" "<<TT<<endl;
   
    for(int j=0;j<nn;j++){
      for(int i=0;i<nn;i++){
	oout<<v[i+nn*j][0]<<" "<<v[i+nn*j][1]<<" "<<v[i+nn*j][2]<<endl;
      }
      oout<<endl;
    }
    oout<<endl;

  
    ofstream outcount_file(outcount_filename);
     outcount_file<<outcount<<endl;
     outcount_file.close();
     outcount++;
     //out_v_con_len_file();
     out_met_file(met_filename);
     
    }
   }
   oout.close();
   }
