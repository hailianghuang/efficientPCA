/* finding the eigenvalues of a complex matrix */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define lineBuffer 1000000

void main(int argc, char *argv[])	
{	
	int N, M, i, j, info, LWORK, LIWORK, IL, IU, *IWORK, v; //v is not used
	FILE *fp, *fp_value, *fp_vector;
	char *next, *ptr, sline[lineBuffer], JOBZ, UPLO, RANGE, input_dir[100], cor_file[100], val_file[100], vec_file[100];
	double r, ABSTOL, *WORK;
	double imem=0;
	
	if(argc!=4){
		printf("Need to specify parameters\n");
		return;
	}

	sscanf(argv[1], "%d", &N);
	printf("Size of correlation matrix: %d\n", N);

	sscanf(argv[2], "%d", &M);
	printf("Looking for the top %d eigen values\n", M);

	sscanf(argv[3], "%s", &input_dir);

	strcpy(cor_file, input_dir);
	strcat(cor_file, "cor.ld");

	strcpy(val_file, input_dir);
	strcat(val_file, "value.txt");

	strcpy(vec_file, input_dir);
	strcat(vec_file, "vector.txt");

	printf("Read correlation matrix from %s\n", cor_file);
	printf("Write eigen values to %s\n", val_file);
	printf("Write eigen vectors to %s\n", vec_file);

	IU=N;
	IL=N-M+1;
	JOBZ ='V';
	RANGE='I';
	UPLO ='U';
	ABSTOL=1E-6;
			
	printf("Allocating memory for matrix A (%dx%dx%d): ", N, N, sizeof(double));
	imem+=N*N*sizeof(double);
	double *A=(double*) malloc(N*N*sizeof(double));
	if(A==NULL){
		printf("Failed\n");
		return;
	}else
		printf("Success\n");

	printf("Allocating memory for matrix E (%dx1x%d): ", N, sizeof(double));
	imem+=N*sizeof(double);
	double *E=(double*) malloc(N*sizeof(double));
	if(E==NULL){
		printf("Failed\n");
		return;
	}else
		printf("Success\n");

	printf("Allocating memory for matrix Z (%dx%dx%d): ", N,M, sizeof(double));
	imem+=N*M*sizeof(double);
	double *Z = (double*) malloc(M*N*sizeof(double));
	if(Z==NULL){
		printf("Failed\n");
		return;
	}else
		printf("Success\n");

	printf("Allocating memory for matrix ISUPPZ (%dx%dx%d): ", M, 2, sizeof(int));
	imem+=2*M*sizeof(int);
	int *ISUPPZ = (int*) malloc(2*M*sizeof(int));
	if(ISUPPZ==NULL){
		printf("Failed\n");
		return;
	}else
		printf("Success\n");

	printf("Total Memory usage is %g MB\n", imem/1024/1024);
	
	printf("Start to load data\n");
	fp = fopen(cor_file, "r");
	i=0;
	while(fgets(sline, lineBuffer, fp)){
		ptr = sline;
		for (j=0; j<N; j++){
	    	next = strchr(ptr, ' ');
			if(next !=NULL) *next='\0';
			sscanf(ptr, "%lg", &r);
			A[j+i*N]=r;
//fortran is column major, C is row major
//			AT[j+size*i]=A[j][i]
			ptr = next+1;
//			if(j==N-1)
//				printf("%g\t",r);
		}
		if(i%2000==0)
			printf("%d SNPs read\n", i);
		i++;
		if(i==N) break;
	}
	fclose(fp);
	printf("Done loading the data, start to solve the eigen system.\n");

	LWORK=-1;
	WORK=(double*) malloc(sizeof(double));
	LIWORK=-1;
	IWORK=(int*) malloc(sizeof(int));

	dsyevr_(&JOBZ, &RANGE, &UPLO, &N, A, &N, &v, &v, &IL, &IU, &ABSTOL, &M, E, Z, &N,ISUPPZ, WORK, &LWORK, IWORK, &LIWORK,&info);

	printf("Optimal WORK size is %d\n", (int)WORK[0]);
	printf("Optimal IWORK size is %d\n", IWORK[0]);
	
	LWORK=WORK[0]*2;
	LIWORK=IWORK[0]*2;
	
	free(WORK);
	free(IWORK);

	printf("Allocating memory for matrix WORK (%dx1x%d): ", LWORK, sizeof(double));
	imem+=LWORK*sizeof(double);
	WORK=(double*) malloc(LWORK*sizeof(double));
	if(WORK==NULL){
		printf("Failed\n");
		return;
	}else
		printf("Success\n");	

	printf("Allocating memory for matrix IWORK (%dx1x%d): ", LIWORK, sizeof(int));
	imem+=LIWORK*sizeof(int);
	IWORK=(int*) malloc(LIWORK*sizeof(int));
	if(IWORK==NULL){
		printf("Failed\n");
		return;
	}else
		printf("Success\n");	
			
	printf("Total Memory usage is %g MB\n", imem/1024/1024);
	dsyevr_(&JOBZ, &RANGE, &UPLO, &N, A, &N, &v, &v, &IL, &IU, &ABSTOL, &M, E, Z, &N,ISUPPZ, WORK, &LWORK, IWORK, &LIWORK,&info);
	
	printf("Done calculations, start to write the result for the top %d eigen values.\n", M);

	fp_value=fopen(val_file, "w" ); 
	fp_vector=fopen(vec_file, "w" );
	if (info==0){
		for (i=M-1; i>=0; i--)
			fprintf(fp_value, "%f\n", E[i]);
		for ( j=0; j<N; j++){
			if(j%500==0)
				printf("%d eigen vectors processed\n", j);
			for (i=M-1; i>=0; i--){
				fprintf(fp_vector, "%f", Z[j+N*i]);
				if(i==0){
					fprintf(fp_vector, "\n");
				}else{
					fprintf(fp_vector, "\t");
				}
			}
		}
	}else 
		printf("An error occured: %d", info);
	fclose(fp_value);
	fclose(fp_vector);
	printf("Successfully completed\n");
	//free others as well
	free(A);	
}
