//http://en.wikiversity.org/wiki/C_source_code_to_find_the_median_and_mean
// Median and mean

#include <iostream>
#include <math.h>
#include <stdio.h>
using namespace std;

int main()
{


int i,k,i1,k1,nx,ny,m;
float mean(int,int[]);
float median(int,float[]);
int R;

//1. Load data from the file into 3 massives x,y,I:
FILE *fp;
fp=fopen("image.txt", "r");

fscanf(fp,"%d\t%d\t%d",&nx,&ny,&R);

float I[nx*ny], Inew;
for(i=0;i<nx*ny;i++){
    fscanf(fp,"%.3f",&I[i]);}
fclose(fp);

//2. Median filter
FILE *fk;
fk = fopen("image_new.txt","w");
for(k=0;k<ny;k++){
	for(i=0;i<nx;i++){
		m = 0;
		float I1[nx*ny];
		for(k1=0;k1<ny;k1++){		 
			for(i1=0;i1<nx;i1++){
				if ((i-i1)*(i-i1) + (k-k1)*(k-k1)<=R*R){
					I1[m] = I[k1*nx+i1]; 
					m=m+1;}
					    }
				     }			
		Inew = median(m,I1);
		fprintf(fk,"%.3f\n",Inew);}
		}
fclose(fk);
}





/*
float mean(int m, int a[]) {
    int sum=0, i;
    for(i=0; i<m; i++)
        sum+=a[i];
    return((float)sum/m);
}
*/
 
 
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
