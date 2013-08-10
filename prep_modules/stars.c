//http://en.wikiversity.org/wiki/C_source_code_to_find_the_median_and_mean
// Median and mean
 
#include<stdio.h>

void main() 
{
int i,k,i1,k1,nx,ny,n;
float mean(int,int[]);
float median(int,int[]);
float R;

//1. Load data from the file into 3 massives x,y,I:
FILE * myfile;
myfile = fopen("stars.txt","r");

fscanf("%d\t%d\t%d",&nx,&ny,&n);
float hwfm[n];
for(i=0;i<n;i++){
    fscanf("%d\t%d\t%.3f",&x1[i],y1[i],hwfm[i]);}
fclose(myfile);

//2. Masking
FILE * myfile1;
FILE * myfile2;
myfile1 = fopen("badpix.txt","a");
myfile2 = fopen("stars.txt","w");
for(i1=0;i1<n;i1++){
	for(k=0;k<ny;k++){
		for(i=0;i<nx;i++){
			if (x1[i1]-i)*(x1[i1]-i)-(y1[i1]-k)*(y1[i1]-k)<=hfwm[i1]:
				fprintf(myfile1,"%i\t%i",i,k)		 
				fprintf(myfile2,"%i\t%i",i,k)		
fclose(myfile1);
fclose(myfile2);
}

float mean(int m, int a[]) {
    int sum=0, i;
    for(i=0; i<m; i++)
        sum+=a[i];
    return((float)sum/m);
}
 
 
float median(int n, float x[]) {
    float temp;
    int i, j;
    // the following two loops sort the array x in ascending order
    for(i=0; i<n-1; i++) {
        for(j=i+1; j<n; j++) {
            if(x[j] < x[i]) {
                // swap elements
                temp = x[i];
                x[i] = x[j];
                x[j] = temp;
            }
        }
    }
 
    if(n%2==0) {
        // if there is an even number of elements, return mean of the two elements in the middle
        return((x[n/2] + x[n/2 - 1]) / 2.0);
    } else {
        // else return the element in the middle
        return x[n/2];
    }
}
