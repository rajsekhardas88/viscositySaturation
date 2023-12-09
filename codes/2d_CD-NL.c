/******************************************************************************************************
  CELL MODEL 
  THIS MC SIMULATOR USES: NEIGHBOR LISTS + CELL SUBDIVISION
 *******************************************************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <mpi.h>
#include <time.h>
#define PI                          3.141592653589793
/*---------------------------- CONSTANTS IN THE SYSTEM ------------------------------------------------*/
#define	N	  	        			        500
#define	DOUBLE_N			        	    500.0
#define NA					                500
#define	TIME_STEP           				10.0
#define	HALF_TIME_STEP		        	0.5*TIME_STEP
#define DENSITY 				            0.90
#define KB                          1.38e-11

#define	CUTOFF_SQRD				         (1.0*1.00)
#define	DR					               0.40
#define	MAX_NEBZ				           1000
#define	MAX_CELLS     	 			     (N/2) 
#define	MAX_BUDS_PER_CELL    			 1000
#define TAU_T       					     5.0
#define	TWO_PI			      		     6.28318530717958
#define LARGE                      1
#define SMALL                      0
#define DIVIDE_TO                  2
#define D0					               0.010
/*-------------------- MINIMIZER SHIT ---------------------------------------------------------------------*/
#define	MAX_ITERATIONS			    	 500000
#define	TOL					               1e-18 		
#define	LINMIN_MAX_ITERATIONS	     40
#define	LINMIN_G1				           2.0 // factors for growing and shrinking the interval- don't change
#define	LINMIN_G2				           1.25
#define	LINMIN_G3				           0.5
#define	LAST_X_DEFAULT				     0.0001
/*------------------------------------ BORING STUFF FROM NUMERICAL RECIPES -----------------------------------*/
#define	IA					               16807
#define	IM					               2147483647
#define	AM					               (1.0/IM)
#define	IQ					               127773
#define	IR					               2836
#define	NTAB					             32
#define	NDIV					             (1+(IM-1)/NTAB)
#define	RNMX					             (1.0-DBL_EPSILON)
#define	ITMAX					             200
#define	EPS					               3.0e-7
#define	FPMIN					             1.0e-30
#define	UNI					               ((double)rand()/((double)RAND_MAX + 1.0)) 
#define SWAP(a,b) 				         temp=(a);(a)=(b);(b)=temp;
#define M 					               7
#define NSTACK 					           50
/*-------------------------------------------------------------------------------------------------------*/
//oid calculateForces();
void calculateForcesCD();
void initializeSystem();
void updateNebzLists();
int min();

/*-------------------------- GLOBALS VARIABLES -----------------------------------------------------------*/
double wx;
double wz;
double wy;
double LENGTH;
double HALF_LENGTH;
double L; 			/*	length of <square> box 		*/
double invL;			/*	inverse of length		*/
double V;			/*	volume				*/
double T;			/*	temperature			*/
double inst_T;			/* 	instantaneous temperature	*/
double t;			/*	time				*/
double u;			/*	potential energy		*/
double particleU;		/* 	energy of single particle	*/
double kinetic;			/*	kinetic energy			*/
double rx[N];			/*	x component of position		*/
double ry[N];			/*	y component of position		*/
double rxUnFolded[N];			/*	x component of position		*/
double ryUnFolded[N];			/*	y component of position		*/
double px[N];			/*	x component of momentum		*/
double py[N];			/*	y component of momentum		*/
double fx[N];			/*	x component of force		*/
double fy[N];			/*	y component of force		*/
int    type[N]; 			/*	type large or small		*/
double R[N];
double gammaMt[N];

double expFactor;
double sqrtVariance;	
double eta[N];
double sqrtoneByR;
double deltaR;
double hijt;
double Aijt;
double nui = 0.5;
double Ei[N];
double fad[N];
double GAMMAM[N];
double NU;
double MAXDIA;
/*-------------- restart job options -------------------------------------------------------------------*/
int restartEq;
int restartRun;
int iterFinish;
/*---------------------------- PARAMETER SFOR THE SYSTEM -----------------------------------------------*/
/*------------------------------------- FOR NEBZ LIST -----------------------------------------------------*/
int nebz[N][MAX_NEBZ]; 		//maximum we give here MAX_NEBZ nebz. 
int numOfNebz[N]; 		//number of elements in nebz[N][MAX_NEBZ]
double maxD,listCutOffSqrd; 	//for updating nebz list.
int nebListCounter; 		//for counting how many steps have took place for nebz list updating.
/*------------------------------------FOR CELL SUBDIVISION ------------------------------------------------*/
double cellCutOff;
int numInCell[MAX_CELLS];
int cells[MAX_CELLS][MAX_BUDS_PER_CELL];
const int shift[4][2] = {{-1,0},{-1,1},{0,1},{1,1}}; //{yShift,xShift};
/*--------------------- FOR LOGARITHAMIC STORAGE -----------------------------------------------------------*/
int *nextStepIndex;
int serial;
int for_ran1;
/*------------------------------ FOR MINIMIZER -------------------------------------------------------------*/

/*----------------------------PARAMETERS FOR PARALLELLIZATION--------------------------------------------------*/
int NUM_TASK, TASK_ID, jobNo;
/**-------------------------------------------- MAIN PROGRAM ----------------------------------------------**/
//initialize with a negative integer. then just give ONE positive one. 
double ran1(int *idum){
    int j;
    int k;
    static int iy=0;
    static int iv[NTAB];
    double temp;

    if ( *idum <= 0 || !iy) { 
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        for (j=NTAB+7;j>=0;j--) { 
            k = *idum/IQ;
            *idum=IA*(*idum-k*IQ)-IR*k;
            if (*idum < 0) *idum += IM;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=*idum/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k; 
    if (*idum < 0) *idum += IM; 
    j=iy/NDIV; 
    iy=iv[j]; 
    iv[j] = *idum;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}

void shellSort(int n, int *a){
    int i,j,inc;
    int v;
    inc=1; //Determine the starting increment.
    do {
        inc *= 3;
        inc++;
    } while (inc <= n);
    do { //Loop over the partial sorts.
        inc /= 3;
        for (i=inc+1;i<=n;i++) { //Outer loop of straight insertion.
            v=a[i];
            j=i;
            while (a[j-inc] > v) { //Inner loop of straight insertion.
                a[j]=a[j-inc];
                j -= inc;
                if (j <= inc) break;
            }
            a[j]=v;
        }
    } while (inc > 1);
    return;
}

double normal(double mean, double variance){
    double v1,v2,rsq;
    do{
        //v1=2.0*UNI-1.0;
        //v2=2.0*UNI-1.0;
        v1=2.0*ran1(&for_ran1)-1.0;
        v2=2.0*ran1(&for_ran1)-1.0;
        rsq =v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    return sqrt(variance)*( v1*sqrt(-2.0*log(rsq)/rsq) ) + mean;
}
int intSqrt(int number){
    int k=1;
    while ( k*k < number )
        k++;
    return k;
}

void initializeGrid(){
    int i,sqrtN,j,temp,tempx,tempy;
    double space,temp1;
    FILE*OP;
    //char radFilename[128];
    //FILE *radFile;
    //op = fopen("rad.dat","wb");
    //sprintf(radFilename,"rad_%.3d.dat",jobNo);
    //radFile = fopen(radFilename,"wb");
    printf("Number of Particle = %d\n",N);
    OP = fopen("INI-CONFIG.dat","wb");
    for (i=0;i<N;i++){
        if(i<=NA) type[i] = 0;
        else type[i] = 0;
    }
    for (i=0; i<N; i++){
        j = (int)(DOUBLE_N*UNI);
        temp = type[i];
        type[i] = type[j];
        type[j] = temp;
    }
    MAXDIA =0.0;
    for(i=0;i<N;i++){
        R[i] = normal(8.5,4.5);
        temp1 = R[i];
        if(temp1 > MAXDIA)
              MAXDIA = temp1; 
        space =  2.0*ran1(&for_ran1)-1.0;
        Ei[i] = 0.001 -space*0.0001; 
        //fprintf(radFile,"%.3f\t%.5f\n",R[i],Ei[i]);
        space =  2.0*ran1(&for_ran1)-1.0;
        eta[i] = 0.005 -space*0.0005;
        //fad[i] = 0.0001;
        fad[i] = 0.0;
    }
    MAXDIA = 2.0*MAXDIA;
    printf("Maxdia = %lf\n",MAXDIA);
    NU = 0.045;
    //fclose(radFile);
    for(i=0; i<N; i++){
        GAMMAM[i] = 6.0*PI*eta[i]*R[i];
    }
    V = 0.0;
    for(i=0;i<N;i++){
        V += PI*R[i]*R[i];
    }
    V = V/DENSITY;
    sqrtN = intSqrt(N);
    L = round(sqrt(V));
    LENGTH = L;
    HALF_LENGTH = 0.5*L;
    printf("Simulation box size L = %.10f\n",L);
    printf("Number DEnsity is = %lf\n",V/(L*L));
    space = 1.0/(double)sqrtN;
    /* for(i=0; i<N; i++){
       rx[i] = space*(double)(i%sqrtN);
       ry[i] = space*(double)(i/sqrtN);
    //px[i] = 2.0*ran1(&for_ran1)-1.0;
    //py[i] = 2.0*ran1(&for_ran1)-1.0;
    } */
    for(i =0;i<N;i++){
        rx[i] = normal(0.5,0.15)*L;
        ry[i] = normal(0.5,0.15)*L;
        //px[i] = 2.0*ran1(&for_ran1)-1.0;
        //py[i] = 2.0*ran1(&for_ran1)-1.0;
    }
    fprintf(OP,"ITEM: TIMESTEP\n");
    fprintf(OP,"%d\n",0);
    fprintf(OP,"ITEM: NUMBER OF ATOMS\n");
    fprintf(OP,"%d\n",N);
    fprintf(OP,"ITEM: ATOMS x y type\n");
    for(i=0;i<N;i++){
        fprintf(OP,"%.10g\t%.10g\t%d\n",(rx[i]),(ry[i]),type[i]);
    }
    fclose(OP);
    initializeSystem();
    return;
}

void saveSnapShot(char *fileName){
    int i;
    FILE *file;
    file = fopen(fileName,"wb");
    fprintf(file,"%f\t%.10f\t%.10f\n",L,u,T);
    for (i=0; i<N; i++){
        fprintf(file,"%.8f\t%.8f\t%d\n",rx[i],ry[i],type[i]);
    }
    fclose(file);
    return;
}

void saveCompleteState(char *outFileName){
    int i;
    FILE *file;

    file = fopen(outFileName,"wb");
    fprintf(file,"%.15g\t%.8g\t%.10g\t%.10f\t%.10f\n",L,u,T,1.1,kinetic); 
    for (i=0; i<N; i++){
        fprintf(file,"%.8g\t%.8g\t%d\t%.8g\t%.8g\n",rx[i],ry[i],type[i]
                ,px[i],py[i]);
    }
    fclose(file);

    return;
}

void saveCompleteRestartState(char *outFileName, int iter){
    int i;
    FILE *file;

    file = fopen(outFileName,"wb");
    fprintf(file,"%.15g\t%.8g\t%.10g\t%d\t%.10f\t1.1\t1.1\n",L,u,T,iter,kinetic); 
    for (i=0; i<N; i++){
        fprintf(file,"%.8g\t%.8g\t%d\t%.8g\t%.8g\n",rx[i],ry[i],type[i]
                ,px[i],py[i]);
    }
    fclose(file);

    return;
}


void readCompleteRestartState(char *inFileName){
    int i;
    double dummy;
    FILE *file;

    //strain = 0.0;//this we need for the minimizer as it is defined for general case
    file = fopen(inFileName,"rb");
    fscanf(file, "%lf",&(L));
    fscanf(file, "%lf",&(u));
    fscanf(file, "%lf",&(T));
    fscanf(file, "%d",&(iterFinish));
    fscanf(file,"%lf",&dummy);
    fscanf(file,"%lf",&dummy);
    fscanf(file,"%lf",&dummy);

    for (i=0;i<N;i++){
        fscanf(file, "%lf",&(rx[i]));
        fscanf(file, "%lf",&(ry[i]));
        fscanf(file, "%d",&(type[i]));
        fscanf(file, "%lf",&(px[i]));
        fscanf(file, "%lf",&(py[i]));
    }
    fclose(file);
    initializeSystem();

    return;
}

void fixDrift(){
    double Px,Py;
    int i;
    Px = 0.0; Py = 0.0; 
    for (i=0;i<N;i++){
        Px += px[i]; Py += py[i];
    }
    Px = Px/DOUBLE_N; Py = Py/DOUBLE_N;
    for (i=0;i<N;i++){
        px[i] -= Px; py[i] -= Py;
    }

    return;
}
void applyPBC(){
    int i;
    for(i=0;i<N;i++){
        if (rx[i] >= LENGTH){
            while(rx[i] >= LENGTH)
                rx[i] = rx[i] - LENGTH;
        }
        else if (rx[i] < 0.0){
            while(rx[i] < 0.0)
                rx[i] = rx[i] + LENGTH;
        }
        if (ry[i] >= LENGTH){
            while(ry[i] >= LENGTH)
                ry[i] = ry[i] - LENGTH;
        }
        else if (ry[i] < 0.0){
            while(ry[i] < 0.0)
                ry[i] = ry[i] + LENGTH;
        }
    }
    return;
}

void adjustStartingTemperature(){
    int i;
    double xi_T;

    kinetic = 0.0;
    for(i=0;i<N;i++){
        kinetic += px[i]*px[i] + py[i]*py[i];
    } 
    kinetic = 0.50*kinetic;
    inst_T = (2.0*kinetic)/(2.0*DOUBLE_N - 3.0);
    xi_T = sqrt(T/inst_T);

    for(i=0;i<N;i++){
        px[i] *= xi_T;
        py[i] *= xi_T;
    }
    return;
}
void initializeSystem(){
    double rcutSq,sr2,sr6,sr12;
    int i,j;
    double temp1,temp2;
    FILE*FP;
    invL = 1.0/L;
    V = L*L;
    cellCutOff = MAXDIA*sqrt(CUTOFF_SQRD) + DR*MAXDIA;
    listCutOffSqrd = cellCutOff*cellCutOff;
    for(i=0;i<N;i++){
        rxUnFolded[i] = rx[i];
        ryUnFolded[i] = ry[i];
    }
    applyPBC();
    updateNebzLists();
    calculateForcesCD();
    temp1 =0.0; temp2 =0.0;
    FP = fopen("FORCES.txt","wb");
    for(i=0;i<N;i++){
        temp1 += fx[i];
        temp2 += fy[i];
        fprintf(FP,"%.10lf\t%.10lf\n",fx[i],fy[i]);
    }
    fclose(FP);
    printf("\nsumFx = %lf sumFy = %lf\n",temp1,temp2);
    nebListCounter = 0;
    maxD = 0.0;
    return;
}
/*********************************** FILE OUTPUT FUNCTIONS ************************************************************/
void updateNebzLists(){
    int i,j,x,y,current,m,mSqrd,ii,jj;
    int a,b,c,k,l,numHere,numThere,target,w,xIndex,yIndex;
    double dx,dy,r2,invCellSize;
    double rxi,ryi;

    nebListCounter = 0;
    maxD = 0.0;
    for (i=0; i<N; i++)
        numOfNebz[i] = 0;
    for (i=0;i<N-1;i++){
        rxi = rx[i];
        ryi = ry[i];
        for (j=i+1; j<N; j++){
            dx = rx[j] - rxi;
            dy = ry[j] - ryi;
            //if ( dx >= HALF_LENGTH )
            //    dx -= LENGTH;
            //else if ( dx < -HALF_LENGTH )
            //    dx += LENGTH;
            //if ( dy >= HALF_LENGTH )
            //    dy -= LENGTH;
            //else if ( dy < -HALF_LENGTH )
            //    dy += LENGTH;
            dx = dx - round(dx/L)*L;
            dy = dy - round(dy/L)*L;
            r2 = ( dx*dx + dy*dy );
            if (r2 < listCutOffSqrd){
                nebz[i][numOfNebz[i]] = j;
                nebz[j][numOfNebz[j]] = i;
                numOfNebz[i]++;
                numOfNebz[j]++;
            }
        }
    }
    return;
}
void calculateForcesCD(){
    int i,j,m,k,typeI,typeJ;
    double r;
    double dx,dy,dz,rxi,ryi,rzi;
    double temp1,temp2;
    double factor,factor2,Lij,F;

    for (i=0;i<N;i++){
        fx[i] = 0.0;
        fy[i] = 0.0;
        gammaMt[i] = 0.0;
    }
    for (i=0;i<N-1;i++){
        rxi = rx[i];
        ryi = ry[i];
        typeI = type[i];
        m = numOfNebz[i];
        for(k=0;k<m;k++){
            j = nebz[i][k];
            if(j>i){
                dx = rx[j] - rxi;
                dy = ry[j] - ryi;
                //if ( dx >= HALF_LENGTH )
                //    dx -= LENGTH;
                //else if ( dx < -HALF_LENGTH )
                //    dx += LENGTH;
                //if ( dy >= HALF_LENGTH )
                //    dy -= LENGTH;
                //else if ( dy < -HALF_LENGTH )
                //    dy += LENGTH;
                dx = dx - round(dx/L)*L;
                dy = dy - round(dy/L)*L;

                typeJ = type[j];
                r = sqrt((dx*dx + dy*dy));
                deltaR = R[i]+R[j] - r;

                if(deltaR > 0.0){
                    sqrtoneByR = sqrt(1.0/R[i] +1.0/R[j]);
                    hijt = deltaR;
                    Aijt = 4.0*r*r*R[i]*R[i] - (r*r - R[j]*R[j] + R[i]*R[i])*(r*r - R[j]*R[j] + R[i]*R[i]);
                    Aijt = sqrt(abs(Aijt));
                    Lij  = Aijt/r;

                    temp1 = pow(hijt,1.5);
                    //temp1 = pow(hijt,1.0);
                    temp2 = 0.75*((1.0 - nui*nui)/Ei[i] + (1.0 - nui*nui)/Ei[j])*sqrtoneByR;
                    factor = temp1/temp2;
                    factor2 = Lij*fad[i];

                    fx[j] += (factor - factor2)*dx/r;
                    fx[i] -= (factor - factor2)*dx/r ;

                    fy[j] += (factor - factor2)*dy/r;
                    fy[i] -= (factor - factor2)*dy/r;
                }
            }
        }
    }

    return;
}
void advanceTimeBD(){
    int i;
    double temp,dx,dy,dz,sum;
    double gfactor,dtbygammam,dtbygammamt,GAMMAMt,gfactorOld;
    double NX,NY;
    for(i=0;i<N;i++){
        GAMMAMt = GAMMAM[i] ;
        dtbygammamt = TIME_STEP/(GAMMAMt);
        //gfactorOld = sqrt(2.0*KB*T/(GAMMAMt));
        //printf("gfactor = %lf\n and NU = %lf\n",gfactorOld,NU);
        gfactor = sqrt(TIME_STEP)*NU;
        NX = normal(0.0,1.0);
        dx = dtbygammamt*(fx[i]) + (gfactor*NX);
        rxUnFolded[i] += dx;
        temp = rx[i] + dx;
        if (temp >= LENGTH)
            rx[i] = temp - LENGTH;
        else if (temp < 0.0)
            rx[i] = temp + LENGTH;
        else
            rx[i] = temp;
        NY = normal(0.0,1.0);
        dy = dtbygammamt*(fy[i]) + (gfactor*NY);
        ryUnFolded[i] += dy;
        temp = ry[i] + dy;
        if (temp >= LENGTH)
            ry[i] = temp - LENGTH;
        else if (temp < 0.0)
            ry[i] = temp + LENGTH;
        else
            ry[i] = temp;
        temp = (dx*dx + dy*dy);


        if (temp > maxD)
            maxD = temp;


        //DRR[i] = sqrt(dx*dx+dy*dy);
    }
    nebListCounter++;
    if ( 2.0*((double)nebListCounter)*sqrt(maxD) > DR*MAXDIA )
        updateNebzLists();
    calculateForcesCD();
    return;
}
void advanceTimeBDBAOAB(){
    int i;
    double temp,dx,dy,dz,sum;
    double gfactor,dtbygammam,dtbygammamt,GAMMAMt;
    for(i=0;i<N;i++){
        GAMMAMt = GAMMAM[i] ;
        dtbygammamt = TIME_STEP/(GAMMAMt);
        gfactor = sqrt(0.5*KB*T*(GAMMAMt)/TIME_STEP);

        dx = dtbygammamt*(fx[i] + gfactor*(wx+ normal(0.0,1.0)));
        rxUnFolded[i] += dx*invL;

        temp = rx[i] + dx*invL;
        if (temp >= LENGTH)
            rx[i] = temp - LENGTH;
        else if (temp < 0.0)
            rx[i] = temp + LENGTH;
        else
            rx[i] = temp;

        dy = dtbygammamt*(fy[i] + gfactor*(wy+normal(0.0,1.0)));
        ryUnFolded[i] += dy*invL; 

        temp = ry[i] + dy*invL;
        if (temp >= LENGTH)
            ry[i] = temp - LENGTH;
        else if (temp < 0.0)
            ry[i] = temp + LENGTH;
        else
            ry[i] = temp;

        temp = (dx*dx + dy*dy);
        if (temp > maxD)
            maxD = temp;
    }
    nebListCounter++;
    if ( 2.0*((double)nebListCounter)*sqrt(maxD) > DR )
        updateNebzLists();
    calculateForcesCD();
    return;
}


void prepareSavingArray(int runLength, int numOfOrigins, double factor){
    unsigned long long int i,j,k;
    unsigned long long int linearInterval;
    unsigned long long int maximalInterval;
    unsigned long long int offset,current,index;
    char nextStepIndexFileName[128];
    FILE *outFile;

    linearInterval = runLength/numOfOrigins;
    maximalInterval = runLength/DIVIDE_TO;

    current = 0;
    for (k=0;k<numOfOrigins;k++){
        nextStepIndex[current] = k*linearInterval;
        current++;
        offset = 1; //the smallest interval
        while (offset < maximalInterval){
            index = k*linearInterval + offset;
            if (index<runLength){
                nextStepIndex[current] = index;
                current++;
            }
            if ( (int)(offset*factor) == offset )
                offset++;
            else
                offset = (int)(offset*factor);
        }
    }

    shellSort(current, nextStepIndex);
    j=0; i=0;
    while (j<current){
        while ( nextStepIndex[j] == nextStepIndex[i] )
            j++;
        i++;
        nextStepIndex[i] = nextStepIndex[j];
    }
    current = i+1;

    sprintf(nextStepIndexFileName,"nextIndex_%d_%.2f.dat",N,T);
    outFile = fopen(nextStepIndexFileName,"wb");
    for (i=0;i<current;i++)
        fprintf(outFile,"%d\n",nextStepIndex[i]);
    return;
}

void equilibrate(double duration){
    unsigned long long  int i,j,steps;
    int t0,t1;      /*for timing purposes...*/
    char restartFileName[128];
    unsigned long long int start;
    double temp1,temp2;
    t0 = time(0);
    steps = (int)(duration/TIME_STEP);
    sprintf(restartFileName,"%.3d/restartEq_%d_%.2f_%.3d.dat",jobNo,N,T,jobNo);

    if(restartEq){
        start = iterFinish;
        printf("Restarting run at T = %g, u = %g\n\n",T,u);
    }
    else{
        start = 0;
        printf("Starting run at T = %g, u = %g\n\n",T,u);
    }
    //FILE*op;
    //op = fopen("disp.dat","wb");

    for (i=start; i<steps; i++){
        //wx = normal(0.0,1.0);
        //wy = normal(0.0,1.0);
        //wz = normal(0.0,1.0);
        //advanceTimeBDBAOAB();
        advanceTimeBD();
        if(!(i%10000)){
            saveCompleteRestartState(restartFileName,i);
            //for(j=0;j<N;j++){
            //    fprintf(op,"%d\t%lf\t%.10g\n",j,2.0*R[j],DRR[j]);
            //}
      }
    }
    //fclose(op);
    temp1 =0.0; temp2 =0.0;
    for(i=0;i<N;i++){
        temp1 += fx[i];
        temp2 += fy[i];
    }
    printf("sumFx = %lf sumFy = %lf\n",temp1,temp2);
    t1 = time(0);
    printf("equilibration at T = %g - Time in seconds: %d\n\n\n",T, t1-t0);
    for(i=0;i<N;i++){
        rxUnFolded[i] = rx[i];
        ryUnFolded[i] = ry[i];
    }
    return;
}
void run(double duration){
    unsigned long long int k,i,steps,between,current;
    int t0,t1;      
    char stuffFileName[128], dataFileName[128],dataEqFileName[128];
    char tempFileName[128];
    FILE *stuffFile, *dataFile, *dataEqFile, *lambdaFile;
    unsigned long long int counter, start, intDummy, numOfFrames;
    double mu,lambda,lambdaExact;
    double rxOrg[N], ryOrg[N];
    double dummy;
    char restartFileName[128],moveDataFile[128];
    FILE *tempFile;
    double temp1,temp2;
    sprintf(tempFileName,"%.3d/tempData.dat",jobNo);
    sprintf(restartFileName,"%.3d/restartRun_%d_%.2f_%.3d.dat",jobNo,N,DENSITY,jobNo);
    sprintf(stuffFileName,"%.3d/stuff_%d_%.2f_%.3d.dat",jobNo,N,DENSITY,jobNo);
    sprintf(dataFileName,"%.3d/data_%d_%.2f_%.3d.dat",jobNo,N,DENSITY,jobNo);
    sprintf(dataEqFileName,"%.3d/dataEq_%d_%.2f_%.3d.dat",jobNo,N,DENSITY,jobNo);

    t0 = time(0);
    steps = (int)(duration/TIME_STEP);
    between = (int)(steps/200);
    if(restartRun){
        printf("Restarting run at T = %g, u = %g\n\n",T,u);
        stuffFile = fopen(stuffFileName,"ab");
        dataEqFile = fopen(dataEqFileName,"ab");
        dataFile = fopen(dataFileName,"rb"); // lets see whether the writing is done properly
        k = 0;
        while ( fscanf(dataFile,"%lf %lf %d",&(dummy),&(dummy),&(intDummy)) != EOF )
            k++;
        numOfFrames = (int)k/DOUBLE_N;
        current = numOfFrames - 1;
        start = nextStepIndex[current] + 1;
        current = current + 1;
        rewind(dataFile);
        tempFile = fopen(tempFileName,"wb");
        for(k=0;k<numOfFrames;k++){
            for(i=0;i<N;i++){
                fscanf(dataFile,"%lf %lf %d",&(rx[i]),&(ry[i]),&(type[i]));
                fprintf(tempFile,"%.10g\t%.10g\t%d\n",rx[i],ry[i],type[i]);
            }
        }
        fclose(dataFile);
        fclose(tempFile);
        sprintf(moveDataFile,"mv %.3d/tempData.dat %.3d/data_%d_%.2f_%.3d.dat",jobNo,jobNo,N,T,jobNo);
        system(moveDataFile);
        dataFile = fopen(dataFileName,"ab"); // lets see whether the writing is done properly
    }
    else{
        printf("Starting run at T = %g, u = %g\n\n",T,u);
        current = 0;
        start = 0;
        stuffFile = fopen(stuffFileName,"wb");
        dataFile = fopen(dataFileName,"wb");
        dataEqFile = fopen(dataEqFileName,"wb");
    }
    counter = 0;
    for (k=start; k<steps; k++){
        //wx = normal(0.0,1.0);
        //wy = normal(0.0,1.0);
        //wz = normal(0.0,1.0);
        //advanceTimeBDBAOAB();
        advanceTimeBD();
        if ( k==nextStepIndex[current] ){
            //if ( !(k%(int)(steps/500))){
            fprintf(stuffFile,"%d\t%.10g\t%.10g\n",k,u,inst_T);
            fprintf(dataFile,"ITEM: TIMESTEP\n");
            fprintf(dataFile,"%d\n",k);
            fprintf(dataFile,"ITEM: NUMBER OF ATOMS\n");
            fprintf(dataFile,"%d\n",N);
            fprintf(dataFile,"ITEM: ATOMS x y radius type\n");
            for (i=0;i<N;i++)
                fprintf(dataFile,"%.10g\t%.10g\t%.4g\t%d\n",rxUnFolded[i],ryUnFolded[i],R[i],type[i]);
            //fprintf(dataFile,"%6d%6d%6d%20.8lf\t%.10g\t%.10g\n",(i+1),(i+1),type[i],R[i],rxUnFolded[i],ryUnFolded[i]);
            current++;
        }
        if(!(k%between)){
            for (i=0;i<N;i++){
                fprintf(dataEqFile,"%.10g\t%.10g\t%.4g\t%d\n",rx[i],ry[i],R[i],type[i]);
                rxOrg[i] = rx[i];
                ryOrg[i] = ry[i];
            }
        }
        if(!(k%100000)){
            saveCompleteRestartState(restartFileName,k);
            fclose(dataFile);
            fclose(dataEqFile);
            fclose(stuffFile);
            stuffFile = fopen(stuffFileName,"ab");
            dataFile = fopen(dataFileName,"ab");
            dataEqFile = fopen(dataEqFileName,"ab");
        }

        if ( !(k%655360) )
            fixDrift();
        }
        temp1 =0.0; temp2 =0.0;
        for(i=0;i<N;i++){
            temp1 += fx[i];
            temp2 += fy[i];
        }
        printf("sumFx = %lf sumFy = %lf\n",temp1,temp2);
        t1 = time(0);
        printf("data collection at T = %g - Time in seconds: %d\n\n\n",T, t1-t0);
        fclose(stuffFile);
        fclose(dataFile);
        fclose(dataEqFile);
        return;
    }
    int main(int argc,char *argv[]){
        unsigned long long int runLength,numOfOrigins,nextStepIndexSize;
        double factor,temp,eqDuration,prodDuration,T_org;
        char snapFileName[128],dir[128];
        int ierr;
        int i;

        ierr = MPI_Init( &argc, &argv );
        ierr = MPI_Comm_size ( MPI_COMM_WORLD, &NUM_TASK );
        ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &TASK_ID );

        if ( TASK_ID == 0 ){
            printf ( "\n" );
            printf ( "  MPI is initialized and running on %i processors\n", NUM_TASK );
        }
        /*___ COnstants For the Model _____*/
        /*                                  *
         * **********************************/  
        //Ei = 0.009;
        //fad = 0.0001;
        //eta = 0.00005; 
        /*********************************/
        for_ran1 = -(int)time(0);
        ran1(&for_ran1); //initialize the random number
        sscanf(argv[1],"%d",&serial);

        jobNo = serial*NUM_TASK + TASK_ID;
        sprintf(dir,"mkdir %.3d",jobNo);
        system(dir);
        srand( jobNo + 100 - for_ran1);
        for_ran1 = -( (int)time(0) + jobNo);
        ran1(&for_ran1); //initialize the random number

        eqDuration =   500000000.0;
        prodDuration = 2000000000.0;
        /*********** initialization ************/
        T = 300.0;
        T_org = T;

        restartEq = 0;
        restartRun = 0;
        if(restartEq != 1 && restartRun != 1){
            initializeGrid();
            //cool(1e-1,0.20);
            printf("energy after cooling u = %.10g\n",u);
            T = T_org;
            //adjustStartingTemperature();
            equilibrate(eqDuration);
            sprintf(snapFileName,"%.3d/snap_%d_%.2f_%.3d.dat",jobNo,N,T,jobNo);
            saveSnapShot(snapFileName);
        }
        else if ( restartEq == 1  && restartRun == 0){
            sprintf(snapFileName,"%.3d/restartEq_%d_%.2f_%.3d.dat",jobNo,N,T,jobNo);
            //sprintf(snapFileName,"%.3d/restartRun_%d_%.2f_%.3d.dat",jobNo,N,T,jobNo);
            readCompleteRestartState(snapFileName);
            //equilibrate(eqDuration);
            sprintf(snapFileName,"%.3d/snap_%d_%.2f_%.3d.dat",jobNo,N,T,jobNo);
            saveSnapShot(snapFileName);
        }
        else if (restartEq == 0 && restartRun == 1){
            sprintf(snapFileName,"%.3d/restartRun_%d_%.2f_%.3d.dat",jobNo,N,T,jobNo);
            readCompleteRestartState(snapFileName);
        }
        runLength = (int)(prodDuration/TIME_STEP);
        factor = 1.4;
        numOfOrigins = 200;
        nextStepIndexSize = numOfOrigins*(10+(int)( log((double)(runLength/DIVIDE_TO))/log(factor) ) );
        nextStepIndex = (int *)malloc(sizeof(int)*nextStepIndexSize);
        prepareSavingArray(runLength, numOfOrigins, factor);
        run(prodDuration);

        sprintf(snapFileName,"%.3d/snap_%d_%.2f_%.3d.dat",jobNo,N,T,jobNo);
        saveCompleteState(snapFileName);

        free(nextStepIndex);
        ierr = MPI_Finalize ( );

        return 0;
    }
