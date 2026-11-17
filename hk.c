#include <stdlib.h>
#include <math.h>
#include <float.h>
// Held-Karp; d: n*n row-major; out tour order (n ints). returns cost
double hk(int n, const double *d, int *tour){
  int full = 1<<(n-1); int m=n-1;
  double *dp = malloc(sizeof(double)*(size_t)full*m);
  int *par = malloc(sizeof(int)*(size_t)full*m);
  for(size_t i=0;i<(size_t)full*m;i++) dp[i]=DBL_MAX;
  for(int j=0;j<m;j++){ dp[(size_t)(1<<j)*m+j]=d[0*n+j+1]; par[(size_t)(1<<j)*m+j]=-1;}
  for(int S=1;S<full;S++) for(int j=0;j<m;j++){ if(!(S>>j&1)) continue; double v=dp[(size_t)S*m+j]; if(v==DBL_MAX) continue;
     for(int k=0;k<m;k++){ if(S>>k&1) continue; int T=S|(1<<k); double w=v+d[(j+1)*n+k+1];
       if(w<dp[(size_t)T*m+k]){dp[(size_t)T*m+k]=w; par[(size_t)T*m+k]=j;} } }
  double best=DBL_MAX; int bj=-1; int S=full-1;
  for(int j=0;j<m;j++){ double w=dp[(size_t)S*m+j]+d[(j+1)*n]; if(w<best){best=w;bj=j;} }
  int pos=n-1; int j=bj; tour[0]=0;
  while(j>=0){ tour[pos--]=j+1; int pj=par[(size_t)S*m+j]; S&=~(1<<j); j=pj; }
  free(dp); free(par); return best;
}
