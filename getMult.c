/* finding the eigenvalues of a complex matrix */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define lineBuffer 1000000

void main(int argc, char *argv[])	
{	
	int i, j, k, g, P, N, L;
	char *next, *ptr, sline[lineBuffer], input_dir[100], vec_file[100], pcs_file[100], geno_file[100];
	double r;
	FILE *fp;

	P=18288;
	L=10;
	N=64014;
	
	if(argc!=5){
		printf("Need to specify parameters\n");
		return;
	}

	sscanf(argv[1], "%d", &P);
	printf("Size of correlation matrix: %d\n", P);

	sscanf(argv[2], "%d", &L);
	printf("Looking for the top %d eigen values\n", L);

	sscanf(argv[3], "%d", &N);
	printf("Total number of samples: %d\n", N);

	sscanf(argv[4], "%s", &input_dir);

	strcpy(vec_file, input_dir);
	strcat(vec_file, "vector.txt");
	strcpy(pcs_file, input_dir);
	strcat(pcs_file, "pcs.txt");
	strcpy(geno_file, input_dir);
	strcat(geno_file, "geno.raw");

	printf("Read eigen vectors from %s\n", vec_file);
	printf("Write PCS to %s\n", pcs_file);
	
	double *B=(double*) malloc(P*N*sizeof(double));
	printf("Start to load matrix B\n");
	fp = fopen(vec_file, "r");
	i=0;
	while(fgets(sline, lineBuffer, fp)){ 
		ptr = sline;
		for (j=0; j<L; j++){
			next = strchr(ptr, '\t');
			if(next !=NULL) *next='\0';
			sscanf(ptr, "%lg", &r);
			B[j+i*L]=r;
			ptr = next+1;
		}
		i++;
		if(i==P) break;
	}
	printf("Done loading matrix B.\n");
	fclose(fp);

	double *X=(double*) malloc(N*P*sizeof(double));
	printf("Start to load matrix X\n");
	fp = fopen(geno_file, "r");
	i=0;
	fgets(sline, lineBuffer, fp);
	while(fgets(sline, lineBuffer, fp)){
		ptr = sline;
		for(j=0; j<6; j++){
			next = strchr(ptr, ' ');
			if(next !=NULL) *next='\0';
			ptr = next+1;
		}
		for (j=0; j<P; j++){
			next = strchr(ptr, ' ');
			if(next !=NULL) *next='\0';
			if(strcmp(ptr, "NA")==0)
				g=-1;
			else
				sscanf(ptr, "%d", &g);
			X[j+i*P]=g;
			ptr = next+1;
		}
		i++;
		if(i%1000==0){
			printf ("%d individuals read\n", i);
		}
		if(i==N) break;
	}
	printf("Done loading the matrix X.\n");
	fclose(fp);
	
	for (j=0; j<P; j++){
		double mean=0;
		double var=0;
		int cnt=0;
		for ( i=0; i<N; i++){
			if(X[j+P*i]>-0.5){
				mean+=X[j+P*i];
				var += X[j+P*i]*X[j+P*i];
				cnt++;
			}
		}
		mean/=cnt;
		var/=cnt;
		var-=mean*mean;
		var=sqrt(var);
		for ( i=0; i<N; i++){
			if(X[j+P*i]>-0.5){
				if(var==0){
//					printf("Error\n", j);
//					printf("%d, %d, %g, %g\n", j, cnt, mean, var);
					var=1;
				}
				X[j+P*i]= (X[j+P*i]-mean)/var;
			}else{
				X[j+P*i]=0;
			}
		}
		if(j%2000==0)
			printf ("%d SNPs scaled\n", j);
		
//		printf("%d, %d, %g, %g\n", j, cnt, mean, var);
	}
	printf("Start to perform multiplication and write to file\n");
	
	FILE *fp_vector=fopen(pcs_file, "w" );
	for ( i=0; i<N; i++){
		for (j=0; j<L; j++){
			double sum=0;
			for(k=0; k<P; k++){
				sum+=X[P*i+k]*B[L*k+j];
//				printf("%d, %d, %d, %g, %g, %g\n", i, j, k,X[P*i+k], B[L*k+j], sum);
			}
			fprintf(fp_vector, "%f", sum);
			if(j==L-1){
				fprintf(fp_vector, "\n");
			}else{
				fprintf(fp_vector, "\t");
			}
		}
	}
	fclose(fp_vector);

}
