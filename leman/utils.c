#include <stdio.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>

double mem(void)
{
    FILE                   *file;
    char                   proc[256];
    int                    mm, rss;

  sprintf(proc,"/proc/%d/statm",(int)getpid());
  if (!(file = fopen(proc,"r"))) {
    return -1;
  }
  fscanf(file,"%d %d",&mm,&rss);
  fclose(file);
//  printf("pid %d: %d, %d, pg = %d\n",(int)getpid(),mm,rss,getpagesize());
  return ((double)rss) * ((double)getpagesize());
}

double second(void)
{
    struct timeval tp;
    struct timezone tzp;
    int i;
    i = gettimeofday(&tp,&tzp);
    return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

//
//       Fortran interface with *single underscore*
//
double mem_(void)
{
    FILE                   *file;
    char                   proc[256];
    int                    mm, rss;

  sprintf(proc,"/proc/%d/statm",(int)getpid());
  if (!(file = fopen(proc,"r"))) {
    return -1;
  }
  fscanf(file,"%d %d",&mm,&rss);
  fclose(file);
  return ((double)rss) * ((double)getpagesize());
}

double second_(void)
{
    struct timeval tp;
    struct timezone tzp;
    int i;
    i = gettimeofday(&tp,&tzp);
    return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}
