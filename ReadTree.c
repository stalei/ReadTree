//#  © Shahram Talei @ 2020 The University of Alabama
//# Main functions are copied from Sage semi-analytic model. See Croton. et al. 2006 and 2016 for more information.
//#you can redistribute it and/or modify
//#it under the terms of the GNU General Public License as published by
//#the Free Software Foundation; either version 3 of the License, or
//#(at your option) any later version.
//#You should have received a copy of the GNU General Public License
//#along with this program.  If not, see <http://www.gnu.org/licenses/>.
//# This code is distributed as is and there is no warranty or technical support
//# This code reads the merger tree



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

void *mymalloc_s(size_t n);
void myfree_s(void *);

#define NOUT 1
char   SimulationDir_s[512];
int LastSnapShotNr;
int Ntrees;
//int totNHalos;
int *TreeNHalos;
//int TreeFirstHalo;
int *TreeNgals[NOUT];
int *TreeFirstHalo;
int MaxGals;
int FoF_MaxGals;

/////////////////////////////////

struct halo_data_s
{
  // merger tree pointers
  int Descendant;
  int FirstProgenitor;
  int NextProgenitor;
  int FirstHaloInFOFgroup;
  int NextHaloInFOFgroup;

  // properties of halo
  int Len;
  float M_Mean200, Mvir, M_TopHat;  // Mean 200 values (Mvir=M_Crit200)
  float Pos[3];
  float Vel[3];
  float VelDisp;
  float Vmax;
  float Spin[3];
  long long MostBoundID;

  // original position in subfind output
  int SnapNum;
  int FileNr;
  int SubhaloIndex;
  float SubHalfMass;
}
*Halo_s;

struct halo_aux_data_s   // auxiliary halo data
{
  int DoneFlag;
  int HaloFlag;
  int NGalaxies;
  int FirstGalaxy;
}
*HaloAux_s;


/////////////////////////////////


#ifndef MINIMIZE_IO

size_t myfread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  return fread(ptr, size, nmemb, stream);
}

size_t myfwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  return fwrite(ptr, size, nmemb, stream);
}

int myfseek(FILE * stream, long offset, int whence)
{
  return fseek(stream, offset, whence);
}


#else


size_t myfwrite(void *ptr, size_t size, size_t nmemb, FILE * fd)
{
  size_t n;

  n = (size_t) fd;

  if(n == 4)
  {
    if(offset_galaxydata + size * nmemb > maxstorage_galaxydata)
    {
      printf("out of space\n");
      //ABORT(1212);
    }

    memcpy(ptr_galaxydata + offset_galaxydata, ptr, size * nmemb);
    offset_galaxydata += size * nmemb;
    if(offset_galaxydata > filled_galaxydata)
      filled_galaxydata = offset_galaxydata;
  }

  if(n >= 10)
  {
    n -= 10;

    if(offset_galsnapdata[n] + size * nmemb > maxstorage_galsnapdata[n])
    {
      printf("out of space in galaxy storage for snapfile=%zu\n", n);
      //ABORT(1212);
    }

    memcpy(ptr_galsnapdata[n] + offset_galsnapdata[n], ptr, size * nmemb);
    offset_galsnapdata[n] += size * nmemb;
    if(offset_galsnapdata[n] > filled_galsnapdata[n])
      filled_galsnapdata[n] = offset_galsnapdata[n];
  }

  return size * nmemb;
}



size_t myfread(void *ptr, size_t size, size_t nmemb, FILE * fd)
{
  size_t n;

  n = (size_t) fd;

  if(n == 1)
  {
    memcpy(ptr, ptr_treedata + offset_treedata, size * nmemb);
    offset_treedata += size * nmemb;
  }

  return size * nmemb;
}



int myfseek(FILE * fd, long offset, int whence)
{
  size_t n;

  n = (size_t) fd;

  if(n == 1)
    offset_treedata += offset;

  if(n == 2)
    offset_auxdata += offset;

  if(n == 3)
    offset_dbids += offset;

  if(n == 4)
    offset_galaxydata += offset;

  if(n >= 10)
  {
    n -= 10;
    offset_galsnapdata[n] += offset;
  }

  return 0;
}



void load_all_treedata(int filenr)
{
  FILE *fd;
  char buf[1000];
  int ret;
  struct stat filestatus;
  size_t bytes;

  sprintf(buf, "%s/treedata/trees_%03d.%d", SimulationDir_s, LastSnapShotNr, filenr);

  ret = stat(buf, &filestatus);

  if(ret != 0)			/* seems not to exist */
  {
    printf("can't open file `%s'\n", buf);
    //ABORT(1);
  }

  fd = fopen(buf, "r");

  bytes = filestatus.st_size;

  ptr_treedata = mymalloc_s(bytes);
  offset_treedata = 0;

  printf("reading %s... (%lu bytes)\n", buf, bytes);
  fflush(stdout);

  fread(ptr_treedata, 1, bytes, fd);

  printf("done\n");
  fflush(stdout);

  fclose(fd);

  bytes = filestatus.st_size * ALLOCPARAMETER / 20.0;

  for(ret = 0; ret < NOUT; ret++)
  {
    ptr_galsnapdata[ret] = mymalloc_s(bytes);

    offset_galsnapdata[ret] = 0;
    filled_galsnapdata[ret] = 0;
    maxstorage_galsnapdata[ret] = bytes;
  }
}



#endif


//////////////////////

void load_tree_table(int filenr)
{
  int i, n, totNHalos;
  char buf[1000];
  FILE *fd;

  sprintf(buf, "%s/treedata/trees_%03d.%d", SimulationDir_s, LastSnapShotNr, filenr);
  printf("file_loadTable:`%s'\n",buf);
  if(!(fd = fopen(buf, "r")))
  {
    printf("can't open file `%s'\n", buf);
    //ABORT(1);

  }
  myfseek(fd, 0, SEEK_SET); //move to the begining of the file
  myfread(&Ntrees, 1, sizeof(int), fd);
  myfread(&totNHalos, 1, sizeof(int), fd);
  printf("Ntrees - NHalos in load:%d-%d\n",Ntrees,totNHalos);
  TreeNHalos = mymalloc_s(sizeof(int) * Ntrees);
  TreeFirstHalo = mymalloc_s(sizeof(int) * Ntrees);

  for(n = 0; n < NOUT; n++)
    TreeNgals[n] = mymalloc_s(sizeof(int) * Ntrees);
  myfread(TreeNHalos, Ntrees, sizeof(int), fd);

//#ifndef MINIMIZE_IO
  fclose(fd);
//#endif

  if(Ntrees)
    TreeFirstHalo[0] = 0;
  for(i = 1; i < Ntrees; i++)
    TreeFirstHalo[i] = TreeFirstHalo[i - 1] + TreeNHalos[i - 1];

  for(n = 0; n < NOUT; n++)
  {
    for(i = 0; i < Ntrees; i++)
      TreeNgals[n][i] = 0;

    //sprintf(buf, "%s/sage_out/%s_z%1.3f_%d", OutputDir_s, FileNameGalaxies, ZZ[ListOutputSnaps[n]], filenr);

    //if(!(fd = fopen(buf, "w")))
    //{
    //  printf("can't open file `%s'\n", buf);
      //ABORT(1);
    //}
    //fclose(fd);
    //TotGalaxies[n] = 0;
  }
}



void free_tree_table(void)
{
  int n;

  for(n = NOUT - 1; n >= 0; n--)
    myfree_s(TreeNgals[n]);

  myfree_s(TreeFirstHalo);
  myfree_s(TreeNHalos);

#ifdef MINIMIZE_IO
  for(n = NOUT - 1; n >= 0; n--)
    myfree_s(ptr_galsnapdata[n]);
  myfree_s(ptr_treedata);
#endif
}



void load_tree(int filenr, int nr)
{
  int i;
  FILE *fd;

#ifndef MINIMIZE_IO
  char buf[1000];
#endif


#ifdef MINIMIZE_IO
  fd = (FILE *) 1;
  offset_treedata = 0;
#else
  sprintf(buf, "%s/treedata/trees_%03d.%d", SimulationDir_s, LastSnapShotNr, filenr);
  if(!(fd = fopen(buf, "r")))
  {
    printf("can't open file `%s'\n", buf);
    //ABORT(1);
  }
#endif

  myfseek(fd, sizeof(int) * (2 + Ntrees), SEEK_CUR);
  myfseek(fd, sizeof(struct halo_data_s) * TreeFirstHalo[nr], SEEK_CUR);

  Halo_s = mymalloc_s(sizeof(struct halo_data_s) * TreeNHalos[nr]);

  myfread(Halo_s, TreeNHalos[nr], sizeof(struct halo_data_s), fd);
///////
int h=3;
printf("file:%d, tree:%d,Descendant:%d, FirstProgenitor:%d, NextProgenitor:%d\n ",filenr,nr,Halo_s[h].Descendant, Halo_s[h].FirstProgenitor,Halo_s[h].NextProgenitor);
//////
/*
int Descendant;
int FirstProgenitor;
int NextProgenitor;
int FirstHaloInFOFgroup;
int NextHaloInFOFgroup;

*/
#ifndef MINIMIZE_IO
  fclose(fd);
#endif

  //MaxGals = (int)(MAXGALFAC * TreeNHalos[nr]);
  //if(MaxGals < 10000)
  //  MaxGals = 10000;

  FoF_MaxGals = 10000;

  HaloAux_s = mymalloc_s(sizeof(struct halo_aux_data_s) * TreeNHalos[nr]);
  //HaloGal = mymalloc_s(sizeof(struct GALAXY) * MaxGals);
  //Gal = mymalloc_s(sizeof(struct GALAXY) * FoF_MaxGals);

  for(i = 0; i < TreeNHalos[nr]; i++)
  {
    HaloAux_s[i].DoneFlag = 0;
    HaloAux_s[i].HaloFlag = 0;
    HaloAux_s[i].NGalaxies = 0;
  }

}



void free_galaxies_and_tree(void)
{
  //myfree_s(Gal);
  //myfree_s(HaloGal);
  myfree_s(HaloAux_s);
  myfree_s(Halo_s);


}





int main()
{
//////////////////
int filenr,tree, FirstFile=0,LastFile=0;
//char bufz0[1000];
sprintf(SimulationDir_s,"/media/shahram/SD/Sample100Mpc/717");
LastSnapShotNr=264;
//FILE *fd;
for(filenr = FirstFile; filenr <= LastFile; filenr++)
  {
    // if(filenr = FirstFile)
    // {
    /*
         printf("Main: filenr = %d\n", filenr);
         printf("LastSnapShotNr = %d\n", LastSnapShotNr);
         sprintf(filelen, "%s/treedata/Treefile_%03d.txt", SimulationDir_s, LastSnapShotNr);
         //treep = fopen(filelen, "w+");
         //fprintf(treep, "%s/treedata/trees_%03d.%d\n", SimulationDir_s, LastSnapShotNr, filenr);
         //fclose(treep);
    // }
    */
/*
       sprintf(bufz0, "%s/treedata/trees_%03d.%d", SimulationDir_s, LastSnapShotNr, filenr);
       printf("file:`%s'\n",bufz0);
    if(!(fd = fopen(bufz0, "r")))
      continue;  // tree file does not exist, move along
    else
      fclose(fd);
*/
    //sprintf(bufz0, "%s/sage_out/%s_z%1.3f_%d", OutputDir_s, FileNameGalaxies, ZZ[ListOutputSnaps[0]], filenr);
    //if(stat(bufz0, &filestatus) == 0)	 // seems to exist, move along
    //  continue;

    //if((fd = fopen(bufz0, "w")))
    //  fclose(fd);
    load_tree_table(filenr);
printf("load tree table done! Ntrees:%d\n",Ntrees);
///////////////////////////////////////////////

    for(tree = 0; tree < Ntrees; tree++)
    {

    //  if(gotXCPU)
  //      ABORT(5);
     // printf("tree = %d\n", tree);
      if(tree % 10000 == 0)
      {
        // printf("\ttask: %d\tnode: %s\tfile: %i\ttree: %i of %i\n", ThisTask, ThisNode, filenr, tree, Ntrees);

         printf("\ttask:file: %i\ttree: %i of %i\n", filenr, tree, Ntrees);
         fflush(stdout);
      }
     printf("load tree\n");
      load_tree(filenr, tree);

      ///////////////////////////////////////////////////


     // printf("gsl_rng_set\n");
      //gsl_rng_set(random_generator, filenr * 100000 + tree);
      //NumGals = 0;
      //GalaxyCounter = 0;
      //printf("construct galaxsies\n");
      //for(halonr = 0; halonr < TreeNHalos[tree]; halonr++)
      //  if(HaloAux_s[halonr].DoneFlag == 0)
      //  construct_galaxies(halonr);

      //printf("save galaxies\n");
      //save_galaxies(filenr, tree);
     // printf("free_galaxies\n");
      //free_galaxies_and_tree();
    }
   // printf("finalize galaxies\n");
    //finalize_galaxy_file(filenr);
   // printf("free table\n");
    free_tree_table();

    printf("\ndone file %d\n\n", filenr);

  }



printf("Test done!\n");
//////////////////
  return 0;
}
