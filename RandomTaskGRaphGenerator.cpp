# include <stdio.h>
# include <stdlib.h>
# include <iostream>
# include <fstream>
# include <math.h>
//# define NOP 5    /*NOP is the number of processors*/

using namespace std;

int pNOP[4]={2,3,4,5}; //added line 9 and 10 instead of line 5
int NOP;
float shapeparameter[3] =   {0.5,1,2};
float ccr[5]            =   {0.1,0.5,1,5,10};
int *wid,*h,v=0;
FILE *fp;

fstream output_file; //for the input in heft.c

int task();                                  /*Declaration of all the functions*/
float ratio();
float graphshape();
int avgcomp();
int width(float,int);
int out_degree(int);
int below(int);
float compcost(int);
int commcost(int ,float);
void output(int*,float*);

int main()
    {
    //srand(time(NULL));
    int n,k,wDAG,height,nog,*outdegree,i,r,j,temp,*cost;
    float alpha,ccratio,*comp;
    
    fp=fopen("RandomTaskGRaphs.txt","w");         /*Output File*/
    output_file.open("Input.txt",ios::out);
    
    printf("How many graphs you want to generate?\n");
    scanf("%d",&nog);    
    
    for(k=0;k<nog;k++)
          {
          v=0;
          NOP=pNOP[rand()%4]; //generate number of Processors
          n=task();
            printf("%d",n);
            
          ccratio=ratio();
          alpha=graphshape();
          wDAG=avgcomp();
          
          height=(int)(sqrt(n)/alpha);
          
          fprintf(fp,"\n\nGraph # %d \nn=%d, alpha=%f,CC Ratio=%f ,Average Computational Cost=%d,\theight=%d \n\n\n",k+1,n,alpha,ccratio,wDAG,height);
                    
          h=(int*)calloc(height,sizeof(int));
          for(i=0;i<height;i++)                          /*Storing no of edges in each level into an array*/
                  {
                  h[i]=width(alpha,n);
                  fprintf(fp,"width of level %d=%d\n",i+1,h[i]);
                  v=v+h[i];
                  }
          
          fprintf(fp,"\nFinal number of vertices=%d\n\n",v);
          v=v/2;//jsut for fun
          
         output_file<<v<<" "<<NOP<<endl; //First line of heft_input
          cost      =(int*)calloc(v*v,sizeof(int));     /*Communication Cost matrix*/
          outdegree =(int*)calloc(v,sizeof(int));           /*Array to store outdegree of each vertex*/
          comp      =(float*)calloc(v*NOP,sizeof(float));   /*Computational Cost matrix*/
    
          for(i=0;i<v*v;i++)
               cost[i]=-1;
    
          for(i=0;i<v;i++)                                  /*Storing outdegree of each vertex in an array*/
               {
               outdegree[i]=out_degree(i);        
               fprintf(fp,"outdegree of %d=%d\n",i+1,outdegree[i]);
               }
          
         for(i=0;i<v;i++)              
               {
               for(j=0;j<=outdegree[i];j++)                        
                 {
                 if(below(i+1))
                               temp=(rand()%below(i+1)+(v-below(i+1)));/*Randomly selecting vertices to be connected with given vertex*/
                               cost[i+temp*v]=commcost(wDAG,ccratio);     /*Storing the value of Communication cost between Vertices*/
                 }          
               }
          
          for(i=0;i<NOP;i++)
               for(j=0;j<v;j++)
                       comp[i+j*NOP]=compcost(wDAG);        /*Storing the value of Computation cost of a process for each task*/
          output(cost,comp);
          
          }
    fclose(fp);
    output_file.close();
    cout<<"\nYour Graphs have been saved in text file named RandomGraphs.txt"<<endl;
    return 123;
    }
    
int task()                    /*Function to randomly generate number of tasks*/
    {
    int r,x;
    r=rand()%10;
  
    x=(r+1)*10;
    return(x);         
    }
    
float ratio()                 /*Function to randomly determine Communication to Computation Ratio*/
      {
      int temp;
      temp=rand()%5;
      return(ccr[temp]);     
      }

float graphshape()           /*Function to randomly determine Shape Parameter*/
    {
    int temp;
    temp=rand()%3;
    return(shapeparameter[temp]);
    }

int avgcomp()                    /*Function to randomly determine average computational time*/                    
    {
    int r,cos;
    r=rand()%4;
    cos=(r+1)*10;
    return(cos);         
    }

int width(float alpha,int x)     /*Function to randomly determine number of vertices at each level*/
    {
    int i,mean,temp;
    mean=(int)(alpha*sqrt(x));
    wid=(int*)calloc(2*mean-1,sizeof(int));
    for(i=0;i<(2*mean-1);i++)
          wid[i]=i+1;
    temp=rand()%(2*mean-1);
    return(wid[temp]);
    }    

int out_degree(int i)            /*Function to randomly determine number of outgoing edges*/
    {
    int r,temp;
    temp=below(i+1);
    if(!temp)return 0;
    else    return((rand()%temp)+1);     
    }  

int below(int i)                  /*Function to Calculate number of vertices at level lower than any of the given vetex*/
    {
    int j=0,temp=0;
    while(temp<i)
         temp=temp+h[j++]; 
    return(v-temp);
    }

float compcost(int wDAG)          /*Function to randomly determine Computation cost of any Task*/
    {
    int r;
    float cos;
    cos=(rand()%(2*wDAG)+1);
    return(cos);
    }

int commcost(int wDAG,float ccratio) /*Function to randomly determine Communication cost between any two given Tasks*/
    {
    int cos,cDAG;
    
    cDAG=int(wDAG*ccratio);
    cos=(rand()%cDAG+1);
    return(cos);
    }
    
void output(int *cost,float *comp)   /*Function to Print output in a Text File*/
     {
     int i,j;
     fprintf(fp,"\nComputation matrix is\n\n");
     fprintf(fp,"---------\t");
     
     for(i=0;i<NOP;i++)
           fprintf(fp,"*Pr%d*\t",(i+1));
    
     fprintf(fp,"\n");      
     
     for(j=0;j<v;j++)
          {
            fprintf(fp,"Task # %d\t",(j+1));
            for(i=0;i<NOP;i++)
            {
              fprintf(fp,"%d\t",(int)(comp[j+i*v]+1));
              output_file<<(int)(comp[j+i*v]+1)<<"\t"; //computation matrix of tasks and processors
            }
            output_file<<endl;
            fprintf(fp,"\n");          
          }        
        // data communication rate between the processors
      for(int q=0;q<NOP;++q)
      {
        for(int w=0;w<NOP;++w)
        {
          int p=1;
          if(q==w) p=0;
          output_file<<p<<" ";
        }
        output_file<<endl;
      }



     fprintf(fp,"\nCommunictaion matrix is\n\n");
     fprintf(fp,"--From / To---\t");
     for(i=0;i<v;i++)
                fprintf(fp,"#%d\t",i+1);              
    
     fprintf(fp,"\n");  
     for(i=0;i<v;i++)
           {
            fprintf(fp,"From Task # %d\t",(i+1));
            for(j=0;j<v;j++)
            {
              fprintf(fp,"%d\t",cost[i+j*v]);
              output_file<<cost[i+j*v]<<"\t";  //communication matrix between the tasks
            }
            output_file<<endl;
            fprintf(fp,"\n");
           }
     cout<<endl;
     }
