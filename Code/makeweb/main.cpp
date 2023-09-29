



#include "pdef_dynamics_1.1_space.h"

//  basic simulation variables:
int Web_Id=0;
int S;                        // total number of species
int S_c;                      // number consumers (animals, predators)
int S_b;                      // number basal species (plants)
int f_basal;
int N=2;                      // number of nutrients

//  parameters that determine the network topology
double zeta_c=6;              // log_10 of the mean of the consumer (animal) body masses
double sigma_c=3;             // width of the distribution of consumer (animal) body masses
double cutoff_c=1e5;          // half relative cutoff for distribution of consumer body masses
double zeta_b=5;              // log_10 of the mean of the basal (plant) body masses
double sigma_b=3;             // width of the distribution of basal (plant) body masses
double cutoff_b=1e5;          // half relative cutoff for distribution of basal body masses

double m_p_min = 0;           // minimla and maximal log10 body masses in case of uniform distributions
double m_a_min = 2;
double m_p_max = 3;
double m_a_max = 6;

double cutoff=0.01;           // cutoff of the Ricker curve for setting a link between predator and prey
double R_opt=100;             // optimal predator-prey body-mass ratio
double g;                  // width of the Ricker curve (higher g-value -> narrower curve)  DRAW UNIFORM (1-10)


//  parameters of the functional response
double a_0=50;                // scaling factor for the attack rate
double a_c=0.46;                   // exponent for predator body-mass scaling of the attack rate (value is drawn from distribution)
//double mean_a_c = 0.47;       // (mean for the above distribution)
//double sigma_a_c = 0.04;      // (standard deviation for the above distribution)
double a_p=0.15;                   // exponent for prey body-mass scaling of the attack rate (value is drawn from distribution)
//double mean_a_p = 0.15;       // (mean for the above distribution)
//double sigma_a_p = 0.03;      // (standard deviation for the above distribution)
double a_plant = 20;          // constant attack rate coefficient for plants

double h_0=0.4;               // scaling factor for the handling time
double h_c=-0.48;                   // exponent for predator body-mass scaling of the handling time (value is drawn from distribution)
//double mean_h_c = -0.48;      // (mean for the above distribution)
//double sigma_h_c = 0.03;      // (standard deviation for the above distribution)
double h_p=-0.66;                   // exponent for prey body-mass scaling of the handling time (value is drawn from distribution)
//double mean_h_p = -0.66;      // (mean for the above distribution)
//double sigma_h_p = 0.02;      // (standard deviation for the above distribution)

double hill=1.5;                  // Hill coefficient (value is drawn from distribution)
//double mu_hill=1.5;           // (mean for the above distribution)
//double sigma_hill=0.2;        // (standard deviation for the above distribution)
double C_interference=0.8;        // interference competition (value is drawn from distribution)
//double mu_Cint=0.8;           // (mean for the above distribution)
//double sigma_Cint=0.2;        // (standard deviation for the above distribution)

double x_resp_a = 0.314;      // intercept of animal respiration rates
double x_resp_p = 0.138;      // intercept of producer respiration rates
double e_p=0.45;              // assimilation efficiency for plant resources
double e_a=0.85;              // assimilation efficiency for animal resources

//  parameters of the nutrient model
double C1=1;                  // content of first nutrient in plants
double C2=0.5;                // content of second nutrient in plants
double D_N=0.25;              // nutrient turnover rate
double K_min=0.1;             // minimal nutrient uptake half saturation density
double K_rel=0.1;             // width of interval for nutrient uptake half saturation densities
double Nc1=10;               // Nutrient 1 supply concentration
double Nc2=10;             // Nutrient 2 supply concentration

int UNIFORM_MASSES = 1;
int CANNIBALISM = 0;


// parameters for plant migration dynamics
double emigr_a_S_b=0.1;                         // max. emigration rate
double emigr_b_S_b=10;                         // shape of curve for emigration
//double mean_emigr_b_S_b=0;                 // (mean for the above distribution)
//double sigma_emigr_b_S_b=10;                 // (standard deviation for the above distribution)

// parameters for animal migration dynamics
double emigr_a_S_c=0.1;                         // max. emigration rate
double emigr_b_S_c=10;                         // shape of curve for emigration
//double mean_emigr_b_S_c=15;                 // (mean for the above distribution)
//double sigma_emigr_b_S_c=5;                 // (standard deviation for the above distribution)

const double MIN_MIGRATION = 1e-10;         // minimum migration threshold

// disperal parameters for RGG
double D_0;                     // minimal dispersal distance -- defined in Generate_RGG_structure
double eps = 0.05;              // scaling factor of body mass for maximal dispersal distance d_max_i=D_0*m_i^eps
double theta = 1;               // shape of curve for dispersal success



char* directory = getenv("OUTPUTDIR");         // diretory for output files
char* bash_seed= getenv("SEED");            //read in seed from bash
int seed = atoi(bash_seed);              // different replicates need different random numbers!






int main(int argc, char* argv[])
{

    
    gsl_rng_default_seed = seed;                  // seed the rng
    gsl_rng *r=gsl_rng_alloc(gsl_rng_default);
    
        
    
    S_b = gsl_rng_uniform_int(r,3)+10;                      // number of basal (plant) species (5-15)
        f_basal= gsl_rng_uniform_int(r,2)+3;                    // fraction of basal spp. (1:2 ... 1:5)
        S_c = S_b * f_basal;                                   // number of consumer (animal) species
        S= S_b + S_c;                                           // total number of species
    
   
        g=gsl_rng_uniform_int(r,5)+9;                          // With of rickers curve
    
    
    web_calc(r);                    // initialise model
    
    gsl_rng_free(r);
    
    return(0);
}

//  ********** generate a network, simulate population dynamics, and evaluate the final network **********
double *web_calc(gsl_rng *r)
{
    
    //gsl_rng_default_seed = constseed;                  // seed the rng with a constant
    //gsl_rng *f=gsl_rng_alloc(gsl_rng_default);
    
    gsl_matrix *Ap=gsl_matrix_calloc(S,S);                     // adjacency matrix
    gsl_vector *mass=gsl_vector_calloc(S);                     // mean body masses of the species
    gsl_vector *Basalvec=gsl_vector_calloc(S);                  // Vector for basal species
    gsl_matrix *A=gsl_matrix_calloc(S,S);                                 // attack rates
    gsl_matrix *H=gsl_matrix_calloc(S,S);                                 // handling times
    gsl_matrix *Up=gsl_matrix_calloc(S,N);                              // nutrient uptake half saturation densities
    gsl_vector *Disp = gsl_vector_calloc(S);                             // max dispersal distance for each species
    
    printf("1 Start simulation: S_b = %d; S_c = %d\n",S_b,S_c);
    
    pdef_structure(r,Ap,mass,Basalvec);                                   // generate random network structure
    printf("2 Foodweb generated\n");
    
    set_parameters(r,Ap,mass,A,H,Up,Disp);
    printf("3 Parameters set\n");
    
    output(Ap,mass,Basalvec,Up,Disp);
    printf("4 Output generated\n");
   
    //  *********** Free memories ***********
    gsl_matrix_free(Ap);
    gsl_vector_free(mass);
    gsl_vector_free(Basalvec);
    gsl_matrix_free(A);
    gsl_matrix_free(H);
    gsl_matrix_free(Up);
    gsl_vector_free(Disp);
    
    return 0;
}






static void pdef_structure(gsl_rng *r, gsl_matrix *Ap, gsl_vector *mass, gsl_vector *Basalvec)
{
    
    int i,j,flag=0;
    double zeta_act;
    double temp1,temp2,a_max,R;
    
    
    
    while(flag==0)
    {
        flag=1;
        
        gsl_matrix_set_zero(Ap);
        

        

            gsl_vector_set_all(mass,0);
            gsl_vector_set_all(Basalvec,0);
            
            //  ********* determine body masses *********
        
        
                for(i = 0; i< S_c; i++)
                gsl_vector_set(mass,S_b+i,pow(10.,gsl_ran_flat(r,m_a_min,m_a_max)));
                gsl_sort_vector(mass);
                
                gsl_vector_view mass_b_vec=gsl_vector_subvector(mass,0,S_b);
                gsl_vector *mass_b=&mass_b_vec.vector;
                
                for(i = 0; i< S_b; i++)
                {
                gsl_vector_set(mass_b,i,pow(10.,gsl_ran_flat(r,m_p_min,m_p_max)));
                gsl_vector_set(Basalvec,i,1);
                }
                gsl_sort_vector(mass_b);
               
        
        
        
        //  ********** fill the adjacency matrix with Ricker attack rates *************
        for(i=0; i<S_c; i++)
        {
            temp1=gsl_vector_get(mass,S_b+i);
            a_max=pow(gsl_vector_get(mass,S_b+i)/R_opt,0.25);
            
            for(j=0; j<S; j++)
            {
                temp2=gsl_vector_get(mass,j);
                R = temp1/temp2;
                if(pow((R / R_opt) * exp(1 - (R / R_opt)), g) >= cutoff)            //#
                gsl_matrix_set(Ap,S_b+i,j,1);
            }
            
            
            if(!CANNIBALISM)
            gsl_matrix_set(Ap,S_b+i,S_b+i,0);                               // remove canibalistic links if necessary
            
            gsl_vector_view tempp=gsl_matrix_row(Ap,S_b+i);                   // reject networks with consumers or predators without prey
            flag=flag*(1-gsl_vector_isnull(&tempp.vector));
        }
        
        for(i=0; i<S_b; i++)
        {
            gsl_vector_view tempp = gsl_matrix_column(Ap,i);                  // reject networks with uncontrolled basal species
            flag = flag*(1-gsl_vector_isnull(&tempp.vector));
        }
    }
    
    return;
}



//  ********** write all parameters required for the dynamics to the array 'params' **********
static void set_parameters(gsl_rng *r, gsl_matrix *Ap, gsl_vector *mass, gsl_matrix *A, gsl_matrix *H, gsl_matrix *Up, gsl_vector *Disp)
{
    int i,j,k;
    double temp1,temp2,temp3,temp4,temp5,temp6,temp7,R,D_max_scale, D_max_max;
    
    gsl_matrix_memcpy(A,Ap);
    
    
    for(i=0; i<S; i++)
    {
        if(i >= S_b)                                                        // the following lines are only for non-basal species
        {
            gsl_vector_view tempp=gsl_matrix_row(A,i);
            temp1=gsl_blas_dasum(&tempp.vector);                              // temp1 stores the number of prey species of predator i
            gsl_vector_scale(&tempp.vector,1/temp1);                          // reduce attack rates for generalists
            
            temp1 = pow(gsl_vector_get(mass,i),a_c);
            
            for(j=0; j<S; j++)
            {
                if(j < S_b)
                    temp2 = a_plant;
                else
                    temp2 = pow(gsl_vector_get(mass,j),a_p);
                R = gsl_vector_get(mass,i) / gsl_vector_get(mass,j);            // predator-prey body-mass ratio for Ricker curve
                temp3 = temp2 * temp1 * pow((R / R_opt) * exp(1 - (R / R_opt)),g);
                gsl_matrix_set(A,i,j,a_0*temp3*gsl_matrix_get(A,i,j));          // attack rates
                
                temp4 = h_0 * pow(gsl_vector_get(mass,i),h_c) * pow(gsl_vector_get(mass,j),h_p);
                gsl_matrix_set(H,i,j,temp4);                                    // handling times
            }
        }
        
        //printf("Test A-Matrix:\n");
        //show_matrix(A,S,S);
        
        
        //printf("Test H-Matrix:\n");
        //show_matrix(H,S,S);

    }
    
    for(i=0; i<S_b; i++)                                                  // the following lines are only for basal species
    {
        for(j=0; j<N; j++)
        {
            temp7=K_min+K_rel*gsl_rng_uniform(r);         // nutrient uptake half saturation densities
            gsl_matrix_set(Up,i,j,temp7);
        }
    }

    
    
    //printf("Test Up-Matrix:\n");
    //show_matrix(Up,S_b,N);
    
// Set species specific maximal dispersal distance.
   
    
    
    D_max_max=0.4;      // maximal distance for the largest species
    D_max_scale=D_max_max/pow(1e6,eps);  // scaling factor acording to D_max_max for D_max for maximal body size of 1e12 for consumers
    
    for (i=0; i<S_b; i++)
    {
        
        temp6= gsl_rng_uniform(r)*D_max_max;
        
        gsl_vector_set(Disp, i, temp6);        // Calculate max dispersal distance for basal species
        
    }
    
    for (i=S_b; i<S; i++)
        gsl_vector_set(Disp, i, D_max_scale*pow(gsl_vector_get(mass,i), eps));  // Calculate max dispersal distance for consumer species
    

    
    return;
}


static void show_matrix(gsl_matrix *M, int rows, int cols)
{
    
    int i, j;
    for(i=0; i<(rows); i++)
    {
        for(j=0; j<(cols); j++)
        {
            
            if(j==cols-1)
                printf("%.3g \n",
                       (double)gsl_matrix_get(M, i, j));
            else
                printf("%.3g \t",
                       (double)gsl_matrix_get(M, i, j));
            
        }
    }
    return;
}

//***************************************************************************************
//  write the adjacency matrix, bodymasses, and biomasses to a file (in matrix/vector-format)
//***************************************************************************************
static void output(gsl_matrix *Ap, gsl_vector *mass, gsl_vector *Basalvec, gsl_matrix *Up, gsl_vector * Disp)
{
    int i,j;
    
    FILE *file1,*file2, *file3;
    
    char str[99];
    char web_out[99] = "web_";
    char mass_out[99] = "masses_";
    char params_out[99] = "params_";
    
    strcpy(web_out,directory);
    strcpy(mass_out,directory);
    strcpy(params_out,directory);
    
  //  strcat(web_out,"/output/");
  //  strcat(mass_out,"/output/");
  //  strcat(params_out,"/output/");
    
    strcat(web_out,"web_");
    strcat(mass_out,"BodyMass_");
    strcat(params_out,"params_");
    
    strcat(web_out,bash_seed);
    strcat(mass_out,bash_seed);
    strcat(params_out,bash_seed);
    
    strcat(web_out,".out");
    strcat(mass_out,".out");
    strcat(params_out,".out");
    
    file1 = fopen(web_out,"w");
    file2 = fopen(mass_out,"w");
    file3 = fopen(params_out,"w");

    
    for(i=0; i<S; i++)
    {
        for(j=0; j<S; j++)
            fprintf(file1,"%g ",gsl_matrix_get(Ap,i,j));
        fprintf(file1,"\n");
    }
    
   
    
    
    
    
    
    // mass.out
    fprintf(file2, "%s, %s, %s, %s, %s\n","body.mass", "if.basal.spp", "D_Max_i","Nutr.halfsat.dens1","Nurt.halfsat.dens2");
    

        
        for(i=0; i<S; i++)
        {
            fprintf(file2,"%.9g, %.9g, %.9g, %.9g, %.9g\n", gsl_vector_get(mass,i), gsl_vector_get(Basalvec,i), gsl_vector_get(Disp,i),gsl_matrix_get(Up,i,0),gsl_matrix_get(Up,i,1));
        
    
        }


    fprintf(file3, "%s, %s, %s\n", "Number.of.Spp", "Fraction.Basal", "Rickers.shape");
    fprintf(file3, "%d, %d, %.9g\n", S, f_basal, g);

            
            
            
    fclose(file1);
    fclose(file2);
    fclose(file3);
    
    return;
}
