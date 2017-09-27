//comma separated value (CSV)
FILE *fp = new fopen("test.csv", "wt");
double mat[100][1000];

int nrows = 100;
int ncols = 10000;
int i, j;

for (i = 0; i < nrows; i++){
   for(j = 0; j < ncols; j++){
      
      fprintf(fp, "%lf", mat[i][j]);
      if(j==ncols-1){
         fprintf(fp, "\n");
      } else {
         fprintf(fp, ",");
      }
   }
}

fclose(fP);
