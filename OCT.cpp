#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <map>
#define MAX 1000
using namespace std;
int parent_proc[MAX];													
int OUT[MAX][MAX];	
float tp[MAX]={0};												//processor availability
//float time=0;
float aft[MAX]={0};													// clock time for processor														
int comm_mat[MAX][MAX]={0}; 										//process time on diff. processor
int adj[MAX][MAX]={0};												//adjancecy matrix
int p_matrix[MAX][MAX],p_temp[MAX][MAX];						//processor matrix
int nodes=0;													//number of nodes
int n_proc=0;													//number of processors
float rank[MAX][MAX]={0.0};										//oct ranking
float est[MAX][MAX]={0};											// average of all rank on processor
vector<float> ready_list;
vector<int> n_ready;										
ifstream ifs;
void SORT();
void remove_duplicate();												
inline void set_aft(){for(int i=0;i<n_ready.size();i++)aft[i]=0;}
inline void set_rank_zero(int i){for(int j=0;j<nodes;j++)rank[i][j]=0.0;}
inline double round(double d){ return floor(d + 0.5);}			//round off float 
void read_adj();												//to read adj. matrix
void read_processor();											//to read processor matrix
void read_weight();												//to read processor time matrix
bool check_end_node(int);										//check for end nodes
void set_end_nodes();											//sets rank of end nodes to zero
void identity(const int I[MAX][MAX],int w,int O[MAX][MAX]);		
void OCT(int a,int b,int A[MAX][MAX]);
void addMatrix(const int a[MAX][MAX],int b[MAX][MAX]);
void row_wise_min();											// to find row wise min of OUT matrix 				
void assign_rank(int);
float EFT(int ni,int pj,int previous);
void display();													//To display result
void N_entry();													//To find entry nodes in  a graph
void add_mul_mat(float add,float mul,int a[][MAX]);
void Oeft();
void push_children(int i);
int highest_oct();
int g_level(int ni=0,int rank=0){								//find rank of a node in a graph
	if(ni<=0) return rank;
	else if(ni==1)	return 1;
	int save=0;
	for(int i=0;i<nodes;i++)
		if(adj[i][ni]!=-1)
			save=i;
		return g_level(save,rank+1);
}
int main(int argc,char *argv[])
{
	if(argv[1]==NULL){cout<<"\nError No File Supplied!\nUsage: [./a.out file_name.txt]\n";return 1;}
	ifs.open(argv[1]);									/* read from a file */
	ifs>>nodes;
	ifs>>n_proc;
	read_weight();
	read_processor();
	read_adj();ifs.close();										/* read endds */
	set_end_nodes();
	for(int i=nodes-1;i>=0;i--)
	{
		if(check_end_node(i)==true)	set_rank_zero(i);			/* If a terminal node set rank to 0.0*/
		
		else{													/* if not a terminal node*/
			for(int j=i+1;j<nodes;j++)							/*check for all successor nodes*/
			{
				if(adj[i][j]!=-1)								/*if a successor */
				{
					identity(p_matrix,adj[i][j],p_temp);		/* for Cij */
					OCT(i,j,OUT);								/* OCT(Tj,Pw)+W(Tj,Pw)*/	
					addMatrix(p_temp,OUT);						/* OCT(Tj,Pw) + W(Tj,Pw)+ Cij*/ 
					row_wise_min();								/* min(OCT(Tj,Pw) + W(Tj,Pw)+ Cij) */
					assign_rank(i); 							/* max(min(OCT(Tj,Pw) + W(Tj,Pw)+ Cij)) */
				}
			}
		}
	}
	/* calculation of RANKoct */ 
	display();
	cout<<endl;
	for(int i=0;i<nodes;i++){cout<<"Process: "<<i+1;
		int temp=0;
		for(int j=0;j<n_proc;j++){
			cout<<"\t"<<rank[i][j];temp+=rank[i][j];
		}
		ready_list.push_back((float)temp/n_proc);
		cout<<"\t";printf("%.1f \n",(float)temp/n_proc);
	}
	Oeft();										//PEFT Algorithm
	return 0;
}
void Oeft()
{
	int parent=-1,comp=10000;
	int ni=0,save=0;
	int p=0,pro=0;
	float oc;int previous=0;
	N_entry();set_aft();
	while(n_ready.size())
	{
		int test_case=10000;
		
		cout<<endl;
		oc=10000;
		//ni=highest_oct();
		ni=n_ready[0];
		for(int b=0;b<nodes;b++)	if(adj[b][ni]!=-1 && parent!=b)	parent=b;		
		for(int j=0;j<n_proc;j++)
		{
			//int c=n_ready[0];
			if((EFT(ni,j,parent_proc[parent])+rank[ni][j])<oc){
				oc=EFT(ni,j,parent_proc[parent])+rank[ni][j];
				pro=j;	
				comp=EFT(ni,j,parent_proc[parent]);	
			}
			cout<<"EFT: "<<EFT(ni,j,parent_proc[parent])<<endl;
			cout<<"P"<<ni+1<<" RANKu: "<<(EFT(ni,j,parent_proc[parent])+rank[ni][j])<<"\tProcessor: "<<j+1<<endl;
			
		}

		previous=pro;		
		cout<<endl;
		int i=n_ready[0];
		parent_proc[i]=pro;
		tp[pro]=comp;
		aft[i]=comp;
		cout<<"AFT: "<<i+1<<"   "<<aft[i]<<endl;
		save=n_ready[0];														//assign AFT
		n_ready.erase(n_ready.begin());									//pop first element 
		
		for(int l=0;l<n_proc;l++)	cout<<"P[ "<<l+1<<" ]: "<<tp[l]<<endl;
		
		push_children(i);
		SORT();
		int rank=0;
	//	for(int ml=0;ml<nodes;ml++)
	//		cout<<"Node: "<<ml<<" RANK:"<<set_level(ml,0)<<endl;
		for(int kk=0;kk<n_ready.size();kk++)	cout<<"T"<<n_ready[kk]+1<<" ";
		printf("Processor Assigned: %d\t Process: %d\t\n",pro+1,i+1);
		cout<<"-----__________________________________________-----\n";
	getchar();
	}
}

void push_children(int i)
{
	for(int j=i+1;j<nodes;j++)
		if(adj[i][j]!=-1)	n_ready.push_back(j);
}

int highest_oct()
{
	int save=-1;
	int temp=-1;
	for(int i=0;i<n_ready.size();i++)
		if(ready_list[n_ready[i]]>temp){
		 temp=ready_list[n_ready[i]];
			save=i;
	}
n_ready.erase(n_ready.begin()+save);
return save;
}

float EFT(int ni,int pj,int previous)
{
	int temp[MAX]={0};
	float mt;		
	for(int i=0;i<=ni;i++)
	{
		if(adj[i][ni]!=-1)
		{
			for(int j=0;j<n_proc;j++)
			{
				mt=aft[i]+p_matrix[parent_proc[i]][j]*adj[i][ni];
				if(temp[j]<=mt)
					temp[j]=mt;
			}
		}
	}
	for(int i=0;i<n_proc;i++){
		if(tp[i]>temp[i])	est[ni][i]=tp[i];
		else	est[ni][i]=temp[i];
	}
	/*
		for(int j=0;j<=ni;j++)
		{
			if(adj[j][ni]!=-1)								//if  a pred. //
			{	
				add_mul_mat(aft[j],adj[j][ni],p_temp);			/// AFT(Nm)+Cm,i 
			}
			else continue;
			for(int z=0;z<n_proc;z++)
			{
				if(tp[z]>=p_temp[previous][z])	est[ni][z]=tp[z];
				else			est[ni][z]=p_temp[previous][z];		// EST 
			}
		}

	*/
		return est[ni][pj]+comm_mat[ni][pj];	/* EFT */
}

void add_mul_mat(float add,float mul,int a[][MAX])
{
	float mx=0;
	for(int i=0;i<n_proc;i++)
		for(int j=0;j<n_proc;j++)
		{
			a[i][j]=(p_matrix[i][j]*mul)+add;
			//if(a[previous][j]>mx)	mx=(float)a[previous][j];
		}
	//return mx;
}

void N_entry()
{
	bool flag=true;
	for(int i=0;i<nodes;i++)
	{
		for(int j=0;j<nodes;j++)
			if(adj[j][i]!=-1)
			{	
				flag=false;
				break;
			}
		if(flag==true)	
		{
			n_ready.push_back(i);
			flag=true;
			for(int k=0;k<n_proc;k++)
				est[i][k]=0.0;
		}	
	}
}

void assign_rank(int a)
{
	for(int i=0;i<n_proc;i++)
		if(rank[a][i]<OUT[i][0])
			rank[a][i]=OUT[i][0];
}
void row_wise_min()
{
	for(int i=0;i<n_proc;i++)
	{
		int temp=1000;
		for(int j=0;j<n_proc;j++)
			if(OUT[i][j]<temp)	temp=OUT[i][j];
		OUT[i][0]=temp;
	}
}
void OCT(int a,int b,int A[MAX][MAX])
{
	for(int i=0;i<n_proc;i++)
		for(int j=0;j<n_proc;j++)
			A[i][j]=rank[b][j]+comm_mat[b][j];								/* OCT(Tj,Pw)+W(Tj,Pw)*/
		
}
void addMatrix(const int a[MAX][MAX],int b[MAX][MAX])						/* adds nXn Matrix */ 
{
	for(int i=0;i<n_proc;i++)
		for(int j=0;j<n_proc;j++)
			b[i][j]+=a[i][j];
}
void identity(const int I[MAX][MAX],int w,int O[MAX][MAX])
{
	for(int i=0;i<n_proc;i++)
		for(int j=0;j<n_proc;j++)
			O[i][j]=w*I[i][j];
}
void set_end_nodes()
{
	for(int i=0;i<nodes;i++)
	{
		if(check_end_node(i)==true)	{
			for(int j=0;j<nodes;j++) rank[i][j]=0.0;
			//cout<<i<<" "<<rank[i]<<" ";
		}
	}
}
bool check_end_node(int n)
{	for(int j=0;j<nodes;j++)
		if(adj[n][j]!=-1)	return false;
	return true;	

}
void read_adj()
{	for(int i=0;i<nodes;i++)
		for(int j=0;j<nodes;j++)
			ifs>>adj[i][j];
}
void read_processor()
{	for (int i = 0; i <n_proc; i++)
		for(int j=0;j<n_proc;j++)
			ifs>>p_matrix[i][j];
}
void read_weight()
{	for(int i=0;i<nodes;i++)
		for(int j=0;j<n_proc;j++)
			ifs>>comm_mat[i][j];
}
void display()
{
	cout<<"Nodes: "<<nodes<<endl<<"Processor: "<<n_proc<<endl<<"Communication Cost Matrix\n";
	for(int i=0;i<nodes;i++)
	{
		for(int j=0;j<n_proc;j++)
			cout<<comm_mat[i][j]<<"\t";
		cout<<endl;
	}
	cout<<"\nAdj Matrix\n";
	for(int i=0;i<nodes;i++){
		for(int j=0;j<nodes;j++)
			cout<<adj[i][j]<<" ";
		cout<<endl;
	}
	cout<<"\nProcessor Matrix\n";
	for (int i = 0; i <n_proc; i++){
		for(int j=0;j<n_proc;j++)
			cout<<p_matrix[i][j]<<"\t";
		cout<<endl;
	}
}
void SORT()
{
	for(int i=0;i<n_ready.size();i++)
		for(int j=i+1;j<n_ready.size();j++)		
			if(ready_list[n_ready[i]]<ready_list[n_ready[j]])
			{
				int temp=n_ready[i];
				n_ready[i]=n_ready[j];
				n_ready[j]=temp;
	}remove_duplicate();
}
void remove_duplicate()
{
	for(int i=0;i<n_ready.size();i++)
		for(int j=i+1;j<n_ready.size();j++)
			if(n_ready[i]==n_ready[j])	n_ready.erase(n_ready.begin()+j);
}














/*	sample input
10 3
14 16 9
13 19 18
11 13 19
13 8 17
12 13 10
13 16 9
7 15 11
5 11 14
18 12 20
21 7 16
0 1 1 
1 0 1
1 1 0
-1 18 12 9 11 14 -1 -1 -1 -1
-1 -1 -1 -1 -1 -1 -1 19 16 -1
-1 -1 -1 -1 -1 -1 23 -1 -1 -1
-1 -1 -1 -1 -1 -1 -1 27 23 -1
-1 -1 -1 -1 -1 -1 -1 -1 13 -1
-1 -1 -1 -1 -1 -1 -1 15 -1 -1
-1 -1 -1 -1 -1 -1 -1 -1 -1 17
-1 -1 -1 -1 -1 -1 -1 -1 -1 11
-1 -1 -1 -1 -1 -1 -1 -1 -1 13
-1 -1 -1 -1 -1 -1 -1 -1 -1 -1
*/
