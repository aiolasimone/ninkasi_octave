#include <octave/oct.h>
#include <octave/Cell.h>
#include <iostream>
#include <stdio.h>
#include <stdbool.h>
#ifdef __cplusplus
extern "C"
{
#endif
#include <omp.h>
#include <string.h>
#ifdef __cplusplus
}  /* end extern "C" */
#endif

using namespace std;


/*--------------------------------------------------------------------------------*/

double get_value(octave_value val)
{
  NDArray myptr=val.array_value();
  return myptr(0,0);

}


/*--------------------------------------------------------------------------------*/
double convert_val(const char *c, int fmt)
{
  char ii[8];
  double val;
  float *f;
  double *d;
  short int *s;
  long *l;
  int *myi;
  switch(fmt) {
  case -4:
    f=(float *)ii;
    ii[0]=c[3];
    ii[1]=c[2];
    ii[2]=c[1];
    ii[3]=c[0];
    val=(double)(*f);
    break;
  case 2:
    s=(short int *)ii;
    ii[0]=c[1];
    ii[1]=c[0];
    val=(double)(*s);
    break;
    //case 1:
    //val=(double)(*c);
    //break;
  case 4:
    //val=(int)(*c);
    myi=(int *)ii;
    ii[0]=c[3];
    ii[1]=c[2];
    ii[2]=c[1];
    ii[3]=c[0];
    val=(double)(*myi);
    break;
  case -8:
    d=(double *)ii;
    ii[0]=c[7];
    ii[1]=c[6];
    ii[2]=c[5];
    ii[3]=c[4];
    ii[4]=c[3];
    ii[5]=c[2];
    ii[6]=c[1];
    ii[7]=c[0];
    val=(double)(*d);
    break;
  case 8:
    //val=(long)(*s);
    l=(long *)ii;
    ii[0]=c[7];
    ii[1]=c[6];
    ii[2]=c[5];
    ii[3]=c[4];
    ii[4]=c[3];
    ii[5]=c[2];
    ii[6]=c[1];
    ii[7]=c[0];
    val=(double)(*l);
    break;
  case 1:
    val=(int)(c[0]);
    break;
  }
  
  return val;

}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (convert_char_to_num_fits, args, nargout, "Convert a block of unsigned chars into floats/doubles/ints etc..\n")
{
  charMatrix cm=args(0).char_matrix_value();
  char *c=cm.fortran_vec();
  int myformat=(int)get_value(args(1));
  int nb=myformat;
  if (nb<0)
    nb*=-1;

  dim_vector dm=cm.dims();
  //printf("dimensions are %d %d\n",dm(0),dm(1));
  //dim_vector dd(2);
  int nrow=dm(0)/nb;
  int ncol=dm(1);
  //dd(0)=nrow;
  //dd(1)=ncol;
  dim_vector dd(nrow,ncol);
  Matrix mat(dd);
  double *mm=mat.fortran_vec();

  int nelem=nrow*ncol;
  for (int i=0;i<nelem;i++)
    mm[i]=convert_val(c+i*nb,myformat);

  return octave_value(mat);

}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (convert_fits_bintable, args, nargout, "Read TOD data into memory.\n")
{
  charMatrix cm=args(0).char_matrix_value();
  char *c=cm.fortran_vec();
  Matrix my_format=args(1).matrix_value();
  int nrow=(int)get_value(args(2));

  dim_vector dm=my_format.dims();
  int ncol=dm(0)*dm(1);

  //dim_vector dd(2);
  //dd(0)=ncol;
  //dd(1)=nrow;
  dim_vector dd(ncol,nrow);
  Matrix mat(dd);
  double *mm=mat.fortran_vec();
  int *ifmt=(int *)malloc(sizeof(int)*ncol);
  double *dfmt=my_format.fortran_vec();
  int nbyte=0;

  for (int i=0;i<ncol;i++) {
    ifmt[i]=dfmt[i];
    nbyte+=fabs(ifmt[i]);
  }

  int *offset=(int *)malloc(sizeof(int)*ncol);
  offset[0]=0;
  for (int i=1;i<ncol;i++) {
    int doff=fabs(ifmt[i-1]);
    if (doff>10)
      doff-=10;
    offset[i]=offset[i-1]+doff;
    //offset[i]=offset[i-1]+fabs(ifmt[i-1]);
  }

  //#pragma omp parallel for shared(nrow,ncol,mm,c,nbyte,offset,ifmt) default(none)
  for (int i=0;i<nrow;i++) 
    for (int j=0;j<ncol;j++) {
      mm[i*ncol+j]=convert_val(c+i*nbyte+offset[j],ifmt[j]);
      
    }
  
  free(ifmt);

  return octave_value(mat);

}
