#include <iostream>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <cstdlib>
#include <algorithm>
#define MAX 1000
using namespace std;

const float r=0.3;
int nodes,n_proc,p_matrix[MAX][MAX],adj[MAX][MAX],processor_assigned[MAX];
ifstream ifs;
float aft[MAX],weight[MAX][MAX],rank[MAX],rank_proposed[MAX],EFT[MAX][MAX]={0.0},EST[MAX][MAX]={0.0},tp[MAX]={0.0};
vector<int> ready_list;

float weight_ni(int ni);
float weight_abstract(int);
void read();
void display();
float max_nj_succ(int ni);
void algo();
void est(int i);
int Pwik(int p);
bool sort_R(int i,int j){return rank_proposed[i]>rank_proposed[j];}

int main(int argc,char *argv[])
{system("cls");cout<<"\nTask Scheduling For Heterogenous Computing Systems (CrossMark)\n----------------------------------------------------------------\n\n\n";
	if(argv[1]==NULL){cout<<"\nError No File Supplied!\nUsage: [./a.out file_name.txt]\n";return 1;}
	
	ifs.open(argv[1]);									/* read from a file */
	read();ifs.close();display();
	for(int i=0;i<nodes;i++)
		rank[i]=weight_ni(i);
		cout<<"RANK(Proposed)\n\n";								/* display Rank  */
	for(int i=nodes-1;i>=0;i--){
		rank_proposed[i]=max_nj_succ(i)+rank[i];
		printf("Node[%d]\t %0.2f\n",i+1,rank_proposed[i]);
	}
	algo();
	return 0;
}
void algo()
{
	float min=1000;int pro;
	cout<<"\nReady Queue ";
	for(int i=0;i<nodes;i++)
		ready_list.push_back(i);

	sort(ready_list.begin(),ready_list.end(),sort_R);
	for(int i=0;i<nodes;i++)	cout<<ready_list[i]+1<<"\t";cout<<endl;
		cout<<"\n________________________________________________________________________________\n\n";
	while(ready_list.size())
	{
		cout<<"\tPROCESS\t"<<ready_list[0]+1<<endl;
		min=1000;
		est(ready_list[0]);
		cout<<"\nEST\t";
		for(int z=0;z<n_proc;z++)
			cout<<EST[ready_list[0]][z]<<"\t";
		cout<<endl<<"EFT\t";
		
		for(int i=0;i<n_proc;i++)
		{
			EFT[ready_list[0]][i]=weight[ready_list[0]][i]+EST[ready_list[0]][i];
			cout<<EFT[ready_list[0]][i]<<"\t";
			if(EFT[ready_list[0]][i]<=min){
				pro=i;
				min=EFT[ready_list[0]][i];							//selection of min. EFT
			}			
		}//cout<<"\nWeight: "<<weight[ready_list[0]][0]<<"  "<<weight[ready_list[0]][1]<<endl;
		//cout<<endl<<"PWIK: "<<Pwik(ready_list[0])<<"MIN: "<<weight[ready_list[0]][pro]<<"\n";

		if(weight[ready_list[0]][pro]<=Pwik(ready_list[0])){
			processor_assigned[ready_list[0]]=pro;
			aft[ready_list[0]]=EFT[ready_list[0]][pro];
			tp[pro]=EFT[ready_list[0]][pro];
		}
		else{//cout<<weight_ni(ready_list[0])<<" "<<weight_abstract(ready_list[0])<<endl;
			if((weight_ni(ready_list[0])/weight_abstract(ready_list[0]))>=r)				//cross threshold
			{cout<<"\n\nCross-Over Ratio: "<<weight_ni(ready_list[0])/weight_abstract(ready_list[0])<<endl;
				float max=EFT[ready_list[0]][0];
				for(int i=0;i<n_proc;i++)
					if(max<=EFT[ready_list[0]][i]){
						pro=i;max=EFT[ready_list[0]][i];
					}
				processor_assigned[ready_list[0]]=pro;
				tp[pro]=aft[ready_list[0]]=EFT[ready_list[0]][pro];
				cout<<"\n\nCross-Over Detected\n";
			}
			else{
					processor_assigned[ready_list[0]]=pro;
			aft[ready_list[0]]=EFT[ready_list[0]][pro];
			tp[pro]=EFT[ready_list[0]][pro];
					
			}
		}
		printf("\nActual Finish Time:\t%0.1f\nProcessor Selected:\t%d",aft[ready_list[0]],pro+1);
	ready_list.erase(ready_list.begin());
	cout<<endl<<"Processor State: "<<tp[0]<<"\t"<<tp[1]<<endl;
	getchar();
	cout<<"\n_____________________________________________________________\n\n";
	}
}

void est(int ni)
{
	float temp[MAX]={0};
	float mt;		
	for(int i=0;i<nodes;i++)
	{
		if(adj[i][ni]!=-1)
		{
			for(int j=0;j<n_proc;j++)
			{
				mt=aft[i]+p_matrix[processor_assigned[i]][j]*adj[i][ni];
				if(temp[j]<=mt)
					temp[j]=mt;
			}
		}
	}
	for(int i=0;i<n_proc;i++){
		if(tp[i]>temp[i])	EST[ni][i]=tp[i];
		else	EST[ni][i]=temp[i];
	}
		
}

float weight_abstract(int p)
{
	float min=10000;
	float max=0;
	for(int i=0;i<n_proc;i++)
	{
		if(min>=EFT[p][i])
			min=EFT[p][i];
		if(max<=EFT[p][i])
			max=EFT[p][i];
	}
	return	(max-min)/(max/min); 
}

float max_nj_succ(int ni)
{
	float temp=0.0;
	for(int i=0;i<nodes;i++)
		if(adj[ni][i]!=-1){
			if((rank_proposed[i]+adj[ni][i])>temp)
				temp=rank_proposed[i]+adj[ni][i];
		}
	return temp;
}
void read()
{
	ifs>>nodes;
	ifs>>n_proc;
	for(int i=0;i<nodes;i++)
		for(int j=0;j<n_proc;j++)
			ifs>>weight[i][j];
	for (int i = 0; i <n_proc; i++)
		for(int j=0;j<n_proc;j++)
			ifs>>p_matrix[i][j];
	for(int i=0;i<nodes;i++)
		for(int j=0;j<nodes;j++)
			ifs>>adj[i][j];
}

float weight_ni(int ni)
{
	float min=weight[ni][0],max=weight[ni][0];
	for(int i=0;i<n_proc;i++)
	{
		if(min>weight[ni][i])
			min=weight[ni][i];
		if(max<weight[ni][i])
			max=weight[ni][i];
	}
	return float(((max-min))/(max/min));
}

void display()
{
	
	cout<<"Nodes: "<<nodes<<"\tProcessor: "<<n_proc<<endl<<"\nProcessing Cost Matrix\n";
	for(int i=0;i<nodes;i++)
	{
		for(int j=0;j<n_proc;j++)
			cout<<weight[i][j]<<"\t";
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
		cout<<endl<<endl;
	}
}
int Pwik(int p)
{
	int min=weight[p][0];
	for(int i=0;i<n_proc;i++)
		if(min>weight[p][0])
			min=weight[p][0];
	return min;
}