#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

/*################################################################*/
void factorial_cache_init(double **tab,int *tab_size,int size) {
  int i;
  double val;
  if(!(*tab)) {
    (*tab)=(double*)malloc(size*sizeof(double));
    (*tab_size)=size;
    (*tab)[0]=1;
    val=1;
    for(i=1;i<size;i++) {
      val*=i;
      (*tab)[i]=val;
    }
  }    
  else if(size>(*tab_size)) {
    int tab_old_size=(*tab_size);
    (*tab)=(double*)realloc((*tab),size*sizeof(double));
    (*tab_size)=size;
    val=(*tab)[tab_old_size-1];
    for(i=tab_old_size;i<size;i++) {
      val*=i;
      (*tab)[i]=val;
    }
  }
}

/*################################################################*/
/*double factorial(int n)
{
  double ret;
  if(n>factorial_cache_size) factorial_cache_init(n);
  ret=factorial_cache[n];
  return ret; 
}*/

/*################################################################*/
void factorial_cache_free(double **tab) {
  free(*tab);
}
