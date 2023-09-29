

#include "pdef_dynamics_1.1_space.h"

//  basic simulation variables:
int Z =50;                     // number of patches
double a=1;                       //shape parameter a of beta distribution for Mainland-Island Scenario
double b=3;                       //shape parameter b of beta distribution for Mainland-Island Scenario
int n=5;                            // number of clusters for the smallworld scenario
double sigma=0.08;                   // sigma for gaussian distribution in smallworld scenario (the smaller sigma, the higher the clustering)



char* directory = getenv("OUTPUTDIR"); 		// diretory for output files
char* bash_seed= getenv("SEED");            //read in seed from bash
int seed = atoi(bash_seed);              // different replicates need different random numbers!


//  ********** Main function **********
int main(int argc, char* argv[])
{
    int i;

    gsl_rng_default_seed = seed;                  // seed the rng
    gsl_rng *r=gsl_rng_alloc(gsl_rng_default);
    

    Z = 50 ; 	    // number of patches

    
    web_calc(r);					// initialise model

    gsl_rng_free(r);

    return(0);
}

//  ********** generate a network, simulate population dynamics, and evaluate the final network **********
double *web_calc(gsl_rng *r)
{

    gsl_matrix *Loc1       = gsl_matrix_calloc(Z,3);
    gsl_matrix *Loc2       = gsl_matrix_calloc(Z,3);
    gsl_matrix *Loc3       = gsl_matrix_calloc(Z,3);
    gsl_matrix *Loc4       = gsl_matrix_calloc(Z,3);
    gsl_matrix *Loc5       = gsl_matrix_calloc(Z,3);
    gsl_matrix *Loc6       = gsl_matrix_calloc(Z,3);
    

    
    printf("1 Create Landscapes: Z = %d\n",Z);

    
    Generate_RGG_structure(r,Loc1,Loc2,Loc3,Loc4,Loc5,Loc6, a,b,n,sigma);
    printf("4 Landscape generated\n");
    
   
    output(Loc1,Loc2,Loc3,Loc4,Loc5,Loc6);
    printf("8 Output generated\n");

 

//  *********** Free memories ***********
  
    gsl_matrix_free(Loc1);
    gsl_matrix_free(Loc2);
    gsl_matrix_free(Loc3);
    gsl_matrix_free(Loc4);
    gsl_matrix_free(Loc5);
    gsl_matrix_free(Loc6);


    return 0;
}



//***************************************************************************************
//  basic RGG network structure
//***************************************************************************************
static void Generate_RGG_structure(gsl_rng *r, gsl_matrix *Loc1, gsl_matrix *Loc2, gsl_matrix *Loc3, gsl_matrix *Loc4, gsl_matrix *Loc5, gsl_matrix *Loc6, double a, double b, int n, double sigma)
{

   
    int m,i,j,k,temp3;
    double pos_x,pos_y,temp1,temp2,temp4;
    pos_x=0;
    pos_y=0;
    temp1=0;
    temp2=0;
    temp3=0;
    temp4=0;
    
    gsl_vector *X1         = gsl_vector_calloc(Z);
    gsl_vector *X2         = gsl_vector_calloc(Z);
    gsl_vector *X3         = gsl_vector_calloc(Z);
    gsl_vector *Y1         = gsl_vector_calloc(Z);
    gsl_vector *Y2         = gsl_vector_calloc(Z);
    gsl_vector *Y3         = gsl_vector_calloc(Z);
    
    gsl_vector *AreaR      = gsl_vector_calloc(Z);
    gsl_vector *AreaG      = gsl_vector_calloc(Z);
    
    gsl_permutation *Sort1 = gsl_permutation_alloc(Z);
    gsl_permutation *Sort2 = gsl_permutation_alloc(Z);
    gsl_permutation *Sort3 = gsl_permutation_alloc(Z);
    
    
    for (m=0; m<Z; m++)
    {
        temp3= 1e5*(1+gsl_rng_uniform_int(r,100));       // Sets Patch Area
        gsl_vector_set(AreaR,m,temp3);
    }
    
    gsl_vector_memcpy(AreaG,AreaR);
    gsl_sort_vector(AreaG);
    gsl_vector_reverse(AreaG);
    for (j=0;j<Z;j++)

    
        for (m=0; m<Z; m++)
        {
            
            gsl_vector_set(X1,m,gsl_rng_uniform(r));       // Sets X location of Patch m to random value
            gsl_vector_set(Y1,m,gsl_rng_uniform(r));       // Sets Y location of Patch m to random value
            
        }
        
        gsl_sort_vector_index(Sort1,X1);
        
        gsl_matrix_set_col(Loc1,0,X1);
        gsl_matrix_set_col(Loc1,1,Y1);
        
        gsl_matrix_memcpy(Loc2,Loc1);
       
    for (j=0;j<Z;j++)
        {
        
            int index = gsl_permutation_get(Sort1,j);
            
            
        gsl_matrix_set(Loc1,index,2,gsl_vector_get(AreaR,j));
        gsl_matrix_set(Loc2,index,2,gsl_vector_get(AreaG,j));
        }
        
    
    

        for (i=0; i<n; i++)
        {
            pos_x = gsl_rng_uniform(r);
            pos_y = gsl_rng_uniform(r);
            for (m=0; m<Z/n; m++)
            {
                temp1=gsl_ran_gaussian(r,sigma)+pos_x;
                temp2=gsl_ran_gaussian(r,sigma)+pos_y;
                
                if(temp1>0 && temp2>0 && temp1<1 && temp2<1)
                {
                    gsl_vector_set(X2,m+i*(Z/n),temp1);    // Sets X location of Patch m to random value
                    gsl_vector_set(Y2,m+i*(Z/n),temp2);    // Sets Y location of Patch m to random value
                }
                else
                {
                    m--;
                }
                
            }
        }
    
    gsl_sort_vector_index(Sort2,X2);
    
    gsl_matrix_set_col(Loc3,0,X2);
    gsl_matrix_set_col(Loc3,1,Y2);
    
    gsl_matrix_memcpy(Loc4,Loc3);
    
    for (j=0;j<Z;j++)
    {
        
        int index = gsl_permutation_get(Sort2,j);
        
        
        gsl_matrix_set(Loc3,index,2,gsl_vector_get(AreaR,j));
        gsl_matrix_set(Loc4,index,2,gsl_vector_get(AreaG,j));
    }

       
        for (m=0; m<Z; m++)
        {
            gsl_vector_set(X3,m, gsl_ran_beta(r,a,b));       // Sets X location of Patch m to random value
            gsl_vector_set(Y3,m, gsl_ran_beta(r,a,b));       // Sets Y location of Patch m to random value
        }
        gsl_sort_vector_index(Sort3,X3);
        
        gsl_matrix_set_col(Loc5,0,X3);
        gsl_matrix_set_col(Loc5,1,Y3);
        
        gsl_matrix_memcpy(Loc6,Loc5);
    for (j=0;j<Z;j++)
    {
        
        int index = gsl_permutation_get(Sort3,j);
        
        gsl_matrix_set(Loc5,index,2,gsl_vector_get(AreaR,j));
        gsl_matrix_set(Loc6,index,2,gsl_vector_get(AreaG,j));
    }
    
    gsl_vector_free(AreaR);
    gsl_vector_free(AreaG);
    
    gsl_vector_free(X1);
    gsl_vector_free(X2);
    gsl_vector_free(X3);
    gsl_vector_free(Y1);
    gsl_vector_free(Y2);
    gsl_vector_free(Y3);
    
    gsl_permutation_free(Sort1);
    gsl_permutation_free(Sort2);
    gsl_permutation_free(Sort3);
    return;
}



//***************************************************************************************
//  write the adjacency matrix, bodymasses, and biomasses to a file (in matrix/vector-format)
//***************************************************************************************
static void output(gsl_matrix *Loc1,gsl_matrix *Loc2,gsl_matrix *Loc3,gsl_matrix *Loc4,gsl_matrix *Loc5,gsl_matrix *Loc6)
{
    int i;
    
    i=0;

    FILE *file1;
    FILE *file2;
    FILE *file3;
    FILE *file4;
    FILE *file5;
    FILE *file6;
    



    char str[99];
    char landscape_out1[99] = "landscape_";
    char landscape_out2[99] = "landscape_";
    char landscape_out3[99] = "landscape_";
    char landscape_out4[99] = "landscape_";
    char landscape_out5[99] = "landscape_";
    char landscape_out6[99] = "landscape_";
    
    
    strcpy(landscape_out1,directory);
    strcat(landscape_out1,"/");
    strcat(landscape_out1,"RGGR_");
    strcat(landscape_out1,bash_seed);
    strcat(landscape_out1,".out");
    
    strcpy(landscape_out2,directory);
    strcat(landscape_out2,"/");
    strcat(landscape_out2,"RGGC_");
    strcat(landscape_out2,bash_seed);
    strcat(landscape_out2,".out");
    
    strcpy(landscape_out3,directory);
    strcat(landscape_out3,"/");
    strcat(landscape_out3,"SWR_");
    strcat(landscape_out3,bash_seed);
    strcat(landscape_out3,".out");
    
    strcpy(landscape_out4,directory);
    strcat(landscape_out4,"/");
    strcat(landscape_out4,"SWC_");
    strcat(landscape_out4,bash_seed);
    strcat(landscape_out4,".out");
    
    strcpy(landscape_out5,directory);
    strcat(landscape_out5,"/");
    strcat(landscape_out5,"MIR_");
    strcat(landscape_out5,bash_seed);
    strcat(landscape_out5,".out");
    
    strcpy(landscape_out6,directory);
    strcat(landscape_out6,"/");
    strcat(landscape_out6,"MIC_");
    strcat(landscape_out6,bash_seed);
    strcat(landscape_out6,".out");

    
    file1 = fopen(landscape_out1,"w");
    fprintf(file1,"%s, %s, %s, %s\n","patch", "X.Koord", "Y.Koord", "Area");
    for(i=0; i<Z; i++)
    {
        fprintf(file1,"%d, %.9g, %.9g, %.9g\n",i ,gsl_matrix_get(Loc1,i,0) ,gsl_matrix_get(Loc1,i,1),gsl_matrix_get(Loc1,i,2));
    }
    fclose(file1);
    
    file2 = fopen(landscape_out2,"w");
    fprintf(file2,"%s, %s, %s, %s\n","patch", "X.Koord", "Y.Koord", "Area");
    for(i=0; i<Z; i++)
    {
        fprintf(file1,"%d, %.9g, %.9g, %.9g\n",i ,gsl_matrix_get(Loc2,i,0) ,gsl_matrix_get(Loc2,i,1),gsl_matrix_get(Loc2,i,2));
    }
    fclose(file2);
    
    file3 = fopen(landscape_out3,"w");
    fprintf(file3,"%s, %s, %s, %s\n","patch", "X.Koord", "Y.Koord", "Area");
    for(i=0; i<Z; i++)
    {
        fprintf(file3,"%d, %.9g, %.9g, %.9g\n",i ,gsl_matrix_get(Loc3,i,0) ,gsl_matrix_get(Loc3,i,1),gsl_matrix_get(Loc3,i,2));
    }
    fclose(file3);
    
    file4 = fopen(landscape_out4,"w");
    fprintf(file4,"%s, %s, %s, %s\n","patch", "X.Koord", "Y.Koord", "Area");
    for(i=0; i<Z; i++)
    {
        fprintf(file4,"%d, %.9g, %.9g, %.9g\n",i ,gsl_matrix_get(Loc4,i,0) ,gsl_matrix_get(Loc4,i,1),gsl_matrix_get(Loc4,i,2));
    }
    fclose(file4);
    
    file5 = fopen(landscape_out5,"w");
    fprintf(file5,"%s, %s, %s, %s\n","patch", "X.Koord", "Y.Koord", "Area");
    for(i=0; i<Z; i++)
    {
        fprintf(file5,"%d, %.9g, %.9g, %.9g\n",i ,gsl_matrix_get(Loc5,i,0) ,gsl_matrix_get(Loc5,i,1),gsl_matrix_get(Loc5,i,2));
    }
    fclose(file5);
    
    file6 = fopen(landscape_out6,"w");
    fprintf(file6,"%s, %s, %s, %s\n","patch", "X.Koord", "Y.Koord", "Area");
    for(i=0; i<Z; i++)
    {
        fprintf(file6,"%d, %.9g, %.9g, %.9g\n",i ,gsl_matrix_get(Loc6,i,0) ,gsl_matrix_get(Loc6,i,1),gsl_matrix_get(Loc6,i,2));
    }
    fclose(file6);
    
    return;

}
