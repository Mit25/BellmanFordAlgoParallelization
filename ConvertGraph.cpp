//This program convert adjecency list into adjecency matrix
#include <bits/stdc++.h>
using namespace std;

int main(int argc,char* argv[]){
	if(argc==3);
	else{
		cout<<"Invalid Argunments"<<endl;
		return 0;
	}
	char *t1=argv[1];
	string in(t1);
	in="GraphInput/"+in;
	char *t2=argv[2];
	string out(t2);
	out="Input/"+out;
	FILE *fp=fopen(in.c_str(),"r");
	FILE *fp1=fopen(out.c_str(),"w");
	int n;int e;
	fscanf(fp,"%d %d ",&n,&e);
	int adjM[n][n];

	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			adjM[i][j]=INT_MAX;
		}
	}

	for(int i=0;i<e;i++){
		int u,v,w;
		u--;v--;
		fscanf(fp,"%d %d %d ",&u,&v,&w);
		if(u<n && v<n && u>=0 && v>=0){
			adjM[u][v]=w;
			adjM[v][u]=w;
			//cout<<u<<" "<<v<<" "<<w<<endl;fflush(stdout);
		}
	}

	fprintf(fp1,"%d\n",n);
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			fprintf(fp1,"%d ", adjM[i][j]);
		}
		fprintf(fp1,"\n");
	}

	fclose(fp);fclose(fp1);
    return 0;
}