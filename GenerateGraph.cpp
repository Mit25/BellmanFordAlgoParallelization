//This program generates random Adjecency matrix for a Graph 
#include <bits/stdc++.h>
using namespace std;

int main(int argc,char* argv[]){
	if(argc==2);
	else{
		cout<<"Invalid Argunments"<<endl;
		return 0;
	}
	char *tmp=argv[1];
	string fname(tmp);
	fname="Input/"+fname;
	FILE *fp1=fopen(fname.c_str(),"w");
	int n,MaxEdgeWeight;
	cin>>n>>MaxEdgeWeight;
	int *adjM=(int *)malloc(n*n*sizeof(int));

	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			adjM[i*n+j]=rand()%MaxEdgeWeight;
			if(adjM[i*n+j]==(MaxEdgeWeight)-1){
				adjM[i*n+j]=INT_MAX;
			}
		}
	}

	fprintf(fp1,"%d\n",n);
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			fprintf(fp1,"%d ", adjM[i*n+j]);
		}
		fprintf(fp1,"\n");
	}

	fclose(fp1);
    return 0;
}