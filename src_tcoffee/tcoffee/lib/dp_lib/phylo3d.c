
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <search.h>
#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"

Alignment *phylo3d (Alignment *inA, Constraint_list *CL)
{
  //will only keep the sequences whose structure is known
  Alignment *A=aln2trim3d(inA, CL);
  p3D *D=fill_p3D(A, CL);
  int a;
  
  if (!A)printf_exit ( EXIT_FAILURE,stderr, "\nERROR: No contact information could be gathered for any sequence  [FATAL]");  
  
  
  
  //Filter columns containing too many gaps
  D->N=filter_columns_with_gap  (D->col, A, D->max_gap);
  
  //Filter columns with distances
  if (D->maxd<MY_EPSILON)D->maxd=scan_maxd(D);
  D->N=filter_columns_with_dist (A,D->pos, D->col,D->dm3d,D->maxd);
      
  //Compute the first tree or DM
  makerep(D,0);//The first replicates is always the original column list;
  aln2dm(D,A);
  addtree(D,A);
  
  for (a=1; a<=D->replicates;a++)
    {
      makerep(D,1);
      aln2dm(D,A);
      addtree(D,A);
    }
  inA->Tree=A->Tree;
  return inA;    
}

double scan_maxd (p3D *D)
{
  NT_node RT, T;
  float rf, brf;
  float bmaxd=1500;
  Sequence *S=aln2seq(D->A);
  Alignment *A=D->A;
  int scan3D_min=(getenv("min_maxd_4_TCOFFEE"))?atof(getenv("min_maxd_4_TCOFFEE")):5;
  int scan3D_max=(getenv("max_maxd_4_TCOFFEE"))?atof(getenv("max_maxd_4_TCOFFEE")):D->extremed;
  
  int a;

  if (!getenv ("REFERENCE_TREE"))
    {
      RT= compute_cw_tree(D->A);
      fprintf ( stderr, "\n!#ref_tree= [cw] -- %s",tree2string (RT));
    }
  else if ( getenv ("REFERENCE_TREE") && isfile(getenv ("REFERENCE_TREE")))
    {
      RT=main_read_tree (getenv ("REFERENCE_TREE"));
      fprintf ( stderr, "\n!#ref_tree=%s", getenv ("REFERENCE_TREE"));
    }
  else
    {
      RT= compute_cw_tree(D->A);
      vfclose (tree2file (RT, S, "newick",vfopen(getenv ("REFERENCE_TREE"), "w")));
      fprintf ( stderr, "\n!# ref_tree=%s [cw]", getenv ("REFERENCE_TREE"));
    }
  
  fprintf ( stderr, "\n!# scan3D_max=%d", scan3D_max);
  
  brf=0;
  bmaxd=D->extremed;
 
  for (a=scan3D_min; a<scan3D_max; a++)
    {
      static char *treeF=vtmpnam (NULL);
      D->maxd=(double)a*100;
      makerep(D,0);
      filter_columns_with_dist (A,D->pos, D->colrep, D->dm3d, D->maxd);
      
      if (aln2dm(D,A))
	{
	  dist2nj_tree (D->dm, A->name, A->nseq, treeF);
	  T=main_read_tree(treeF);
	  rf=simple_tree_cmp(RT,T, S, 1);
	  fprintf ( stderr, "\n!# Threshold = %3d Angstrom ==> %6.2f %% Similiarity with ref_tree -- Nsites: %d", a, rf, D->nsites) ;
	  
	  if ( rf>brf)
	    {
	      brf=rf;
	      bmaxd=(double)a;
	    }
	}
      else
	{
	  fprintf ( stderr, "\n!# Threshold = %3d Angstrom ==> Missing Values in the distance matrix", a) ;
	}
    }
  if (brf>0)
    {
      fprintf ( stderr, "\n!# Optimal Threshold: %d Angstrom  ==> %.2f %% RF similarity with ref_tree\n", (int)bmaxd, brf);
    }
  else
    {
      fprintf ( stderr, "\n!# WARNING -- Missing Values -- Could not find any suitable threshold - Use max value %d Angstrom -use a higher value for scan3D_max ", (int)bmaxd);
    }
  free_sequence (S, -1);
  
  return bmaxd*100;//in picometers
}

Alignment *phylo3d_gt (Alignment *inA, Constraint_list *CL)
{
  Alignment *A=aln2trim3d(inA, CL);
  p3D *D=fill_p3D(A, CL);
  Sequence *S;
  Alignment *B;
  
  int n, s1, s2, tot, l;
  int **pos;
  char *align_method=getenv ("align_method_4_TCOFFEE");
  if (!A)printf_exit ( EXIT_FAILURE,stderr, "\nERROR: No contact information could be gathered for any sequence  [FATAL]");
  
  //set the parameters
  
  l=2*A->len_aln;
  D->colrep=declare_int ((l*l-l)/2+1,2);

  S=aln2seq(A);

  if (D->maxd<MY_EPSILON)D->maxd=D->extremed*100+1;


  tot=((A->nseq*A->nseq)-A->nseq)/2;
  for (n=0,s1=0; s1<A->nseq-1; s1++)
    {
      int rs1, rs2;
      rs1=name_is_in_hlist (A->name[s1], S->name, S->nseq);
      for ( s2=s1+1; s2<A->nseq; s2++, n++)
	{
	  int  *buf_pos_s1, *buf_pos_s2;
	  char *buf_s1, *buf_s2;
	  
	  rs2=name_is_in_hlist (A->name[s2], S->name, S->nseq); 

	  if (!align_method)B=align_two_sequences (S->seq[s1], S->seq[s2], "blosum62mt", -8, -1, "myers_miller_pair_wise");
	  else B=align_two_structures  (S, rs1, rs2,align_method);
	

	  
	  pos=aln2pos_simple(B,2);
	  msa2column_list (B,D->colrep);
	  
	  //buffer data
	  buf_s1=A->seq_al[s1];
	  buf_s2=A->seq_al[s2];
	  buf_pos_s1=D->pos[s1];
	  buf_pos_s2=D->pos[s2];
	  
	  //swapp data
	  A->seq_al[s1]=B->seq_al[0];
	  A->seq_al[s2]=B->seq_al[1];
	  D->pos[s1]=pos[0];
	  D->pos[s2]=pos[1];
	  
	  
	  
	  
	  D->dm[s1][s2]=D->dm[s2][s1]=pair2dist(D,s1,s2);
	  output_completion (stderr,n,tot,1, "Guide Tree Computation");
	  
	  //restaure data
	  A->seq_al[s1]=buf_s1;
	  A->seq_al[s2]=buf_s2;
	  D->pos[s1]=buf_pos_s1;
	  D->pos[s2]=buf_pos_s2;
	  free_int (pos, -1);
	  free_aln (B);
	}
    }
  
  addtree (D, A);
  inA->Tree=A->Tree;
  return inA;
}

int free_3dD (p3D*D)
{
  vfree (D->tree_mode);
  free_int (D->pos, -1);
  free_arrayN((void***)D->dm3d,3);
  free_arrayN((void**)D->dm,2);
  free_arrayN ((void**)D->col,2);
  vfree (D);
  return 1;
}

p3D* fill_p3D (Alignment *A, Constraint_list *CL)
{
  int a;
  p3D *D=(p3D*)vcalloc ( sizeof (p3D), 1);
  D->A=A;
  D->CL=CL;
  D->S =CL->S;
  
  if (getenv ("max_gap_4_TCOFFEE")){D->max_gap  =atofgetenv("max_gap_4_TCOFFEE");}
  else D->max_gap=0.5;
  
  if (getenv ("TREE_MODE_4_TCOFFEE"))D->tree_mode=csprintf (D->tree_mode, "%s",getenv("TREE_MODE_4_TCOFFEE"));
  else D->tree_mode=csprintf (D->tree_mode, "nj");
  if (getenv ("REPLICATES_4_TCOFFEE"))
    {
      if ( strm (getenv ("REPLICATES_4_TCOFFEE"), "columns"))D->replicates=A->len_aln;
      else D->replicates=atoigetenv ("REPLICATES_4_TCOFFEE");

    }
  else D->replicates=0;

  

  if (getenv ("enb_4_TCOFFEE")){D->enb  =atoi("enb_4_TCOFFEE");}
  else D->enb=3;

    
  if (getenv ("columns4tree_4_TCOFFEE"))D->col=file2column_list (getenv ("columns4tree_4_TCOFFEE"), NULL);
  else D->col=msa2column_list (A, NULL);
  

  

  D->pos=aln2pos_simple (A, A->nseq);
  D->dm3d=aln2dm3d (A, CL, &D->extremed);
  if (getenv ("maxd_4_TCOFFEE"))
    {
      if (strm (getenv ("maxd_4_TCOFFEE"), "scan"))D->maxd=-1;
      else D->maxd  =(double)atof(getenv("maxd_4_TCOFFEE"))*100;
    }
  else D->maxd=D->extremed*100+1;
  


  if (D->maxd>=0 && D->maxd<MY_EPSILON)D->maxd=D->extremed*100;
  
  D->dm =declare_double(A->nseq, A->nseq);
  D->used_site=(int*)vcalloc (A->len_aln, sizeof (int));
  
  if (A->Tree) free_aln (A->Tree); A->Tree=NULL;
  A->Tree=declare_aln2(D->replicates+1, 0);
  (A->Tree)->dmF_list=(char**)vcalloc (D->replicates+1, sizeof (char*));
  for (a=0; a<=D->replicates; a++)(A->Tree)->dmF_list[a]=vtmpnam (NULL);
  (A->Tree)->nseq=0;
  return D;
}
int col2n (int **col)
{
  int i=0;
  if (!col) return 0;
  else if (!col[0])return 0;
  
  while (col[i++][0]!=-1);
  return i-1;
}

p3D * makerep (p3D *D, int mode)
{
  D->colrep=col2colrep(D->col, D->colrep, D->N, mode);
  return D;
}
int **col2colrep (int **colin,int **colout, int ni, int mode)
{
  int i;
  
  if (!colout)colout=declare_int (ni+1,2);
  for (i=0; i<ni; i++)
    {
      int ri=(mode)?rand ()%ni:i;
      colout[i][0]=colin[ri][0];
      colout[i][1]=colin[ri][1];
    }
  colout[i][0]=-1;
  return colout;
}


Alignment * addtree (p3D *D,Alignment *A)
{
  static char *treeF=vtmpnam (NULL);
  static int ml=get_longest_string (A->name, A->nseq, NULL, NULL)+1;
  FILE *fp;
  int s1, s2;
  Alignment *TREEA=A->Tree;
  

  if (A->nseq>2)dist2nj_tree (D->dm, A->name, A->nseq, treeF);
  TREEA->seq_al[TREEA->nseq]=file2string(treeF);
  sprintf (TREEA->name[TREEA->nseq], "%d",TREEA->nseq+1); 	    
  TREEA->seq_al[TREEA->nseq]=file2string(treeF);
  fp=vfopen (TREEA->dmF_list[TREEA->nseq], "w");
  
  fprintf ( fp, "%d \n", A->nseq);
  for ( s1=0; s1<A->nseq;s1++)
    {
      fprintf (fp, "%-*.*s ", ml,ml,A->name[s1]);
      for (s2=0; s2<A->nseq; s2++)
	fprintf (fp, "%6.3f ", (float)((double)D->dm[s1][s2])/(float)100);
      fprintf (fp, "\n");
    }
  fprintf (fp, "\n");
  if (getenv ("PRINT_NSITES"))fprintf (fp, "!# NSITES: %d\n", D->nsites);
  vfclose (fp);
  TREEA->nseq++;
  return TREEA;
}
int filter_columns_with_dist(Alignment *B, int **pos,int **col, int***dm, double maxd)
{
  if ( getenv ("strict_maxd_4_TCOFFEE"))return filter_columns_with_dist_strict(B,pos,col,dm,maxd);
  else return filter_columns_with_dist_relaxed(B,pos,col,dm,maxd);
} 
  
int filter_columns_with_dist_strict(Alignment *B, int **pos,int **col, int***dm, double maxd)
{
  //Keep every pair of columns with ALL distances fulfilling the condition*/
  int i,ni,s;
  i=ni=0;
  if (maxd<MY_EPSILON)return col2n(col);
 
  while (col[i][0]!=-1)
    {
      int count;
      int c1=col[i][0];
      int c2=col[i][1];
      
      for (count=0,s=0; s<B->nseq; s++)
	{
	  int r1=pos[s][c1]-1;
	  int r2=pos[s][c2]-1;
	  if (r1<r2){int rb=r1; r1=r2;r2=rb;}
	  
	  if (r1>=0 && r2>=0 && dm[s][r1-1][r2-1]<=maxd)count++;
	}
      if ( count==B->nseq)
	{
	  col[ni][0]=c1;
	  col[ni][1]=c2;
	  ni++;
	}
      i++;
    }
  col[ni][0]=-1;
  return ni;
}
int filter_columns_with_dist_relaxed(Alignment *B, int **pos,int **col, int***dm, double maxd)
{
  //Keep every pair of columns with at least one distance below the threshold
  int i,ni,s;
  i=ni=0;
  if (maxd<MY_EPSILON)return col2n(col);
  
  while (col[i][0]!=-1)
    {
      int c1=col[i][0];
      int c2=col[i][1];
      
      for (s=0; s<B->nseq; s++)
	{
	  int r1=pos[s][c1]-1;
	  int r2=pos[s][c2]-1;
	  if (r1<r2){int rb=r1; r1=r2;r2=rb;}
	  
	  if (r1>=0 && r2>=0 && dm[s][r1][r2]<=maxd)
	    {
	      col[ni][0]=c1;
	      col[ni][1]=c2;
	      ni++;
	      break;
	    }
	}
      i++;
    }
  col[ni][0]=-1;
  return ni;
}
int filter_columns_with_gap (int **col, Alignment *B, float max_gap)
{
  float *gap=(float*)vcalloc (B->len_aln, sizeof (float));
  int c,s, i, ni;

  for (c=0; c<B->len_aln; c++)
    {
      for (s=0; s<B->nseq; s++)if (B->seq_al[s][c]=='-')gap[c]+=1;
      gap[c]/=(float)B->nseq;
    }
  i=ni=0;
  while (col[i][0]!=-1)
    {
      int c1=col[i][0];
      int c2=col[i][1];
      if (gap[c1]<=max_gap && gap[c2]<=max_gap)
	{
	  col[ni][0]=c1;
	  col[ni][1]=c2;
	  ni++;
	}
      i++;
    }
  vfree (gap);
  col[ni][0]=-1;
  return ni;
}
	  
	  

 int** msa2column_list (Alignment *B, int **col)
{
  int i, j, n;
  if (!col)col=declare_int(((B->len_aln*B->len_aln)-B->len_aln)/2+1,2);
  for (n=0,i=1; i<B->len_aln-1; i++)
    for (j=i+1; j<B->len_aln; j++, n++)
      {
	col[n][0]=i;
	col[n][1]=j;
      }
  col[n][0]=-1;
  return col;
}
 int   **  file2column_list (char *file, int **list)
{
  int i, n;
  char ***l;

  if (!(l=file2list (file, " \t")))return NULL;
    
  i=0;
  while (l[i++]);
  if (!list)list=declare_int (i+1, 3);
  while (l[i])
    {
      
      int n=atoi(l[i][0]);
      if (n>=3 && l[i][1][0]!='#')
	{
	  list[n][0]=atoi (l[i][1])-1;
	  list[n][1]=atoi (l[i][2])-1;
	  n++;
	}
      i++;
    }  
  list[n][0]=-1;
  free_arrayN((void***)l, 3);
  return list;
}
   
int*** aln2dm3d (Alignment *A, Constraint_list*CL, double *max)
{
  int s, r, sA, i;
  Sequence *S=CL->S;
  int ***dm=(int***)vcalloc (A->nseq, sizeof (int**));
  max[0]=0;
  for (sA=0; sA<A->nseq; sA++)
    {
      s=name_is_in_hlist (A->name[sA], S->name, S->nseq);
      dm[sA]=(int**)vcalloc(S->len[s]+1,sizeof (int*));
      for (r=1; r<=S->len[s]; r++)dm[sA][r]=(int*)vcalloc (r, sizeof(int));
      for (r=1; r<=S->len[s]; r++)
	{
	  for (i=1; i<CL->residue_index[s][r][0]; i+=ICHUNK)
	    {
	      int nr =CL->residue_index[s][r][i+R2];
	      int we =CL->residue_index[s][r][i+WE];
	      int ns =CL->residue_index[s][r][i+SEQ2];
	      
	      if (s!=ns){HERE ("Warning: Contact library contains inter-sequence data contacts: %d %d",i,ns);}
	      if (r>nr)dm[sA][r-1][nr-1]=we;
	      else dm[sA][nr-1][r-1]=we;
	      if ((double)we>max[0])max[0]=(double)we;
	    }
	}
    }
  max[0]/=100;//put max back into Angstrom
  return dm;
}


Alignment *aln2trim3d (Alignment *A, Constraint_list *CL)
{
  //Keep only sequences with contact information
  Alignment *B=NULL;
  static char *tmpF=vtmpnam (NULL);
  int ns, s,r,rs;
  FILE*fp;
  if (!CL) return B;
  
  fp=vfopen (tmpF, "w");
  for (ns=0,s=0; s<A->nseq; s++)
    {
      if ((rs=name_is_in_hlist (A->name[s], (CL->S)->name, (CL->S)->nseq)!=-1))
	{
	  for (r=1; r<=(CL->S)->len[s];r++)
	    {
	      if (CL->residue_index[s][r][0])
		{
		  fprintf (fp, ">%s\n%s\n",A->name[s], A->seq_al[s]);
		  ns++;
		  break;
		}
	    }
	}
    }
  vfclose (fp);
  if (ns)B=quick_read_fasta_aln (B,tmpF);
  return B;
}

int aln2dm (p3D *D, Alignment *A)
{
  int rv=1;
  int s1, s2;
  int c, t;
    
  for (c=0; c<A->len_aln; c++)D->used_site[c]=0;
  
  for (s1=0; s1<A->nseq-1; s1++)
    for (s2=s1+1; s2<A->nseq; s2++)
      {
	D->dm[s1][s2]=D->dm[s2][s1]=pair2dist (D,s1, s2);
	if (D->dm[s1][s2]<-99)rv=0;
      }
  for (D->nsites=0,c=0; c<A->len_aln; c++)
    D->nsites+=D->used_site[c];
  
  return rv;
}
double pair2dist(p3D *D, int s1, int s2)
{ 
  //s1 and s2 match the D->dm3d 
  // check the offseet of the residue index - starts at 1 or 0
  double score=0;
  double max=0;
  double rscore, rmax;
  double w1, w2;
  int i;
  int ns=0;
  if (s1>s2){int sb=s1;s1=s2;s2=sb;}

  i=0;
  while (D->colrep[i][0]!=-1)
    {
      int c1=D->colrep[i][0];
      int c2=D->colrep[i][1];
      i++;
      if (c1<c2){int cb=c1;c1=c2;c2=cb;}
      

      int r11=D->pos[s1][c1]-1;
      int r12=D->pos[s1][c2]-1;
      int r21=D->pos[s2][c1]-1;
      int r22=D->pos[s2][c2]-1;
      
      if (r11<0 || r12<0)continue;
      if (r21<0 || r22<0)continue;
      
      if (abs((r12-r11))<D->enb)continue;
      if (abs((r22-r21))<D->enb)continue;
      
      
      
      w1=(double)D->dm3d[s1][r11][r12];
      w2=(double)D->dm3d[s2][r21][r22];
      
      if (w1<MY_EPSILON || w2<MY_EPSILON)continue;
      if (w1>D->maxd || w2>D->maxd)continue;

      ns++;
      phylo3d2score (w1, w2, &rscore, &rmax);
      
      score+=rscore;
      max+=rmax;
      D->used_site[c1]=D->used_site[c2]=1;

    }
  
  if (max<MY_EPSILON)score=-100;
  else score=(double)100*((double)1 - (score/max));
  return score;
}
double phylo3d2score (double w1, double w2, double *rscore, double *rmax)
{
  double we=0;
  double sc=0;
    
  static int   distance_mode=atoigetenv ("THREED_TREE_MODE");
  static double distance_modeE=atofgetenv ("THREED_TREE_MODE_EXP");
  static int no_weights=atoigetenv ("THREED_TREE_NO_WEIGHTS");
  
  if ( distance_modeE<MY_EPSILON)distance_modeE=3;
			   
			   
   //first attempt-- Major issue because non symetrical and therefore not a distance
   if (!distance_mode)
     {
       static int warn;
       we=w1;
       sc=((MIN((w1/w2),(w2/w1))));
       if (warn==0)
	 {
	   add_warning ( stderr, "distance_mode==0 should not be used");
	   warn=1;
	 }
       
     }
   //same as before but symetrical: the distance ratio weighted by the average distance
   else if (distance_mode==1)
     {
       we=(w1+w2)/2;
       sc=((MIN((w1/w2),(w2/w1))));
     }
   //the absolute difference of distance normalized by the average distance
   else if (distance_mode==2)
     {
       
       we=(w1+w2)/2;
       sc=(double) 1-(FABS((w1-w2))/we);
       
     }
   //the absolute difference of distance normalized by the average distance and weighted by the average distance
   else if (distance_mode==3)
     {
       we=(w1+w2)/2;
       sc=(double)1-(FABS((w1-w2))/we);
     }
   else if (distance_mode ==4)
     {
       we=((w1>w2)?w1:w2);
       sc=(double) 1-(FABS((w1-w2))/we);
       
     }
   else if (distance_mode ==5)
     {
       we=(w1+w2);
       sc=(double)1-(FABS((w1-w2))/we);
       
     }
   else if (distance_mode ==6)
     {
       we=1;
       sc=(double)(FABS((w1-w2)));
     }
   
   if (no_weights)
     we=1;
   
   //Compute the score and its normalization
   //If 
   rmax[0]=we;
   rscore[0]=pow(sc,distance_modeE)*we;
   return rscore[0]/rmax[0];
   
}



 
