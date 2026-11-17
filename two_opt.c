#include <math.h>
// 2-opt local search on tour (in place). d row-major n*n.
void twoopt(int n,const double*d,int*t){
  int improved=1;
  while(improved){ improved=0;
    for(int i=0;i<n-1;i++) for(int j=i+2;j<n;j++){
      int a=t[i],b=t[i+1],c=t[j],e=t[(j+1)%n]; if(e==a) continue;
      double delta=d[a*n+c]+d[b*n+e]-d[a*n+b]-d[c*n+e];
      if(delta< -1e-12){ for(int l=i+1,r=j;l<r;l++,r--){int x=t[l];t[l]=t[r];t[r]=x;} improved=1; }
    }}
}
