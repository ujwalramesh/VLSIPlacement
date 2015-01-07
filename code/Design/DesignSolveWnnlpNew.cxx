#include <Design.h>

#include <math.h>
#include "wnlib/wnlib.h"
#include "wnlib/wnasrt.h"
#include "wnlib/wnsll.h"
#include "wnlib/wnnlp.h"
#include <typeinfo>

double myDivideNew(double num, double denom) 
{
        if (denom <= 0.0)
        {
                return (WN_FHUGE);
        }
        else
        {
                if (num < 0.0)
                {
                        num = 0.0;
                }
                return (num/denom);
        }
}

void
getNLPCellsToSolveNew(Design &myDesign,vector<Cell *> &cellsToSolve)
{

Cell *cellPtr;
string cellName;

DESIGN_FOR_ALL_CELLS(myDesign,cellName,cellPtr){
        if ((*cellPtr).CellIsTerminal()) continue;
        cellsToSolve.push_back(cellPtr);
}DESIGN_END_FOR;

}

void
nlpInitialValuesNew(vector<Cell *> &cellsToSolve, double* x)
{

Cell *cellPtr;
string cellName;
uint idx=0;

int size = cellsToSolve.size();
cout << "Size of cells: " << size << endl;
//map<uint, Cell *> cellToSolveMap;

//map<uint, Cell *>::iterator cellToSolveMapItr; 
VECTOR_FOR_ALL_ELEMS(cellsToSolve,Cell *, cellPtr){
   //     cout << "populating Inital values for cell: " << (*cellPtr).CellGetName();
        x[idx] = (*cellPtr).CellGetXposDbl();
        x[idx+1] = (*cellPtr).CellGetYposDbl();
  //      cellToSolveMap[idx] = CellPtr;
       cout << "\t Xpos: " << x[idx] << "\t Ypos: " << x[idx+1] << "\t actualXpos: " << (*cellPtr).CellGetXposDbl() << "\t actualYpos: " << (*cellPtr).CellGetYposDbl()<< endl;
        idx = idx+2;
}END_FOR;

//cout << "I value at the end of iteration is: " << idx << endl;

}

/* Objective function here calculates the Log sum exponential of my HPWL and also adds penalty to my objective based on overlap*/

double
wirelengthObjFunc(int size, double *values,ptr myDesign) 
{
  /* objective here is total_wirelength = sum over all nets ( wirelength (net) * net weight) */
  static int objIterationCount=0;
  Cell *cellPtr;
  string cellName;
  uint idx=0;
  double penaltyParameter;
  penaltyParameter = (*((Design*)myDesign)).DesignGetpenaltyParameter();
        //cout << " Void pointer Type: " <<typeid(*((Design*)myDesign)).name() << endl;  
        DESIGN_FOR_ALL_CELLS((*((Design*)myDesign)),cellName,cellPtr){
                //cout << "Changing Xposition for cell: " << (*cellPtr).CellGetName() << " to " << values[idx] << endl;
                 if ((*cellPtr).CellIsTerminal()) continue;     
                (*cellPtr).CellSetXposDbl(values[idx]);
                (*cellPtr).CellSetYposDbl(values[idx+1]);
                //cout << "The new Xposition of the cell is " << (*cellPtr).CellGetXpos() << endl;
                idx = idx+2;

        }END_FOR;
        ulong LseHPWL;
        LseHPWL = (*((Design*)myDesign)).DesignComputeLseHPWL();
   /*Adding penalty function for density constraint as per Will Naylor's patent*/
   double totalDensityPenalty=0; /*variable for adding penalty for my objective function*/
   
   /* Below penalty method is my penalty based on percentage overlap and tried out of intution by taking (percentageOverlap)^3 to make it differentiable*/  
/*   double totalOverlap;
   totalOverlap = (*((Design*)myDesign)).DesignDumpClusterOverlapForPenalty();
   uint lambda = 2; 
   bool loop1=true;
   bool loop2=true;
   bool loop3=true;
   
   totalDensityPenalty= pow(totalOverlap,3);

   if (totalOverlap > 300 ) {
        totalDensityPenalty= 3000*totalDensityPenalty;
        //lambda=lambda*2;
   } else if ((totalOverlap >200) && (totalOverlap <= 300 )) {
           if (loop1){
                   lambda=2;
                   loop1=false;
           }
           totalDensityPenalty= 5000*totalDensityPenalty;
         //  lambda=lambda*2;
   } else if ((totalOverlap >100 ) && (totalOverlap <= 200 )){
           if (loop2){
                   lambda=2;
                   loop2=false;
           }
           totalDensityPenalty= 8000*totalDensityPenalty;
           //lambda=lambda*2;
   } else if ((totalOverlap > 10 ) && (totalOverlap <= 100 )){
           if (loop3){
                   lambda=2;
                   loop3=false;
           }
           totalDensityPenalty= 12000*totalDensityPenalty;
          // lambda=lambda*2;
   } else if ((totalOverlap <= 10 ) && (totalOverlap >= 0 )) {
           totalDensityPenalty =0 ;
   }

*/

/*  The below penalty method is for penalizing my objective based on grod point potentials */
   (*((Design*)myDesign)).DesignUpdateGridPotentials();
   totalDensityPenalty =(*((Design*)myDesign)).DesignComputeTotalDensityPenalty();
    if (penaltyParameter < 0 ) {
            penaltyParameter = -penaltyParameter;
   }
   //static uint lambda;
   //lambda=2;
        double rtv;
       //rtv = LseHPWL+(penaltyParameter*totalDensityPenalty);
       rtv = LseHPWL+(1*totalDensityPenalty);
        cout << "ulong wirelength is: " << LseHPWL << "\t";
        cout << "Density Penalty is : " << penaltyParameter*totalDensityPenalty << "\t";
        cout << "objective Value is : " << rtv << endl;
        //lambda=lambda*2;
    objIterationCount++;
        return rtv;
}


// The below function is the working gradient function before adding the gradient of the penalty.
// Delete the below function once the gradient function including penalty works.
// A git backup is also saved before adding the new function

/*void 
gradientFunc(double *grad,int size,double *values,ptr myDesign)
{
        static int gradIterationCount=0;
        Cell *cellPtr;
        string cellName;
        uint alpha = 500;
        uint idx=0;
        double cellXpos;  
        double cellYpos;
        DESIGN_FOR_ALL_CELLS((*((Design*)myDesign)),cellName,cellPtr){
                if ((*cellPtr).CellIsTerminal()) continue;
                double gradX=0;
                double gradY=0;
                double cellMaxx;
                double cellMinx;
                double cellMaxy;
                double cellMiny;
                cellXpos = (*cellPtr).CellGetXposDbl();
                cellYpos = (*cellPtr).CellGetYposDbl();
                Net *netPtr;
                double sumPinsPos;
                double pinMaxx;
                double pinMinx;
                double pinMaxy;
                double pinMiny;
                double tempDivideX;
                double tempDivideY;
                double pinXpos;
                double pinYpos;
                CELL_FOR_ALL_NETS_NO_DIR((*cellPtr),netPtr){
                        Pin *pinPtr;
                        NET_FOR_ALL_PINS((*netPtr),pinPtr){
                                Cell* cellParentPtr;
                                cellParentPtr = (*pinPtr).PinGetParentCellPtr();
                                pinXpos = pinPtr->xOffset + cellParentPtr->x;
                                pinYpos = pinPtr->yOffset + cellParentPtr->y;
                                tempDivideX= myDivideNew(pinXpos,alpha);
                                tempDivideY= myDivideNew(pinYpos,alpha);
                                pinMaxx += exp(tempDivideX);
                                pinMinx +=myDivideNew(1,exp(tempDivideX));
                                pinMaxy += exp(tempDivideY);
                                pinMiny +=myDivideNew(1,exp(tempDivideY));
                        }NET_END_FOR;
                }CELL_END_FOR;
                double tempDivX;
                double tempDivY;
                tempDivX = myDivideNew(cellXpos,alpha);
                tempDivY = myDivideNew(cellYpos,alpha);
                cellMaxx = exp(tempDivX);
                cellMinx = myDivideNew(1,exp(tempDivX));
                cellMaxy = exp(tempDivY);
                cellMiny = myDivideNew(1,exp(tempDivY));
                double temp1;
                double temp2;
                temp1 = myDivideNew(cellMaxx,pinMaxx);
                temp2 = myDivideNew(cellMinx,pinMinx);
                gradX = (temp1)-(temp2);
                double temp3;
                double temp4;
                temp3 = myDivideNew(cellMaxy,pinMaxy);
                temp4 = myDivideNew(cellMiny,pinMiny);
                gradY = (temp3)-(temp4);
                grad[idx] = gradX;
                grad[idx+1]=gradY;
                //cout << "gradient value for cell: " << (*cellPtr).CellGetName() << "is gradx: " << gradX << " gradY: " <<gradY << endl;
                if (isnan(gradX)){
                     cout << "gradX is nan: Find the debug information below:" << endl;
                     cout << "temp1 and temp2 are: "<<temp1 <<"\t" <<temp2<<endl;
                     cout << "cellMaxX and CellMinX are: "<<cellMaxx<<"\t"<<cellMinx<<endl;
                     cout << "pinMaxx and PinMinx are: " <<pinMaxx<<"\t"<<pinMinx<<endl;
                     cout << "cellXpos is: "<< cellXpos<<" cellName: " << cellName<<endl;  
                //     grad[idx]=0.1;
                }
                if (isnan(gradY)){
                     cout << "gradY is nan: Find the debug information below:" << endl;
                     cout << "temp3 and temp4 are: "<<temp3 <<"\t" <<temp4<<endl;
                     cout << "cellMaxY and CellMinY are: "<<cellMaxy<<"\t"<<cellMiny<<endl;
                     cout << "pinMaxy and PinMiny are: " <<pinMaxy<<"\t"<<pinMiny<<endl;
                     cout << "cellYpos is: "<< cellYpos<<" cellName: " << cellName<<endl;   
                  //   grad[idx+1]=0.1;
                }
                cout <<"iteration Number: "<< gradIterationCount <<"Cell Name: " <<cellName << " CellXpos: "<<cellXpos<< " cellYpos: "<<cellYpos<<" gradX: "<<gradX<<" gradY: "<<gradY<<endl;  
                
                idx=idx+2;
        }DESIGN_END_FOR;
        gradIterationCount++;
}*/


void 
gradientFunc(double *grad,int size,double *values,ptr myDesign)
{
        static int gradIterationCount=0;
        Cell *cellPtr;
        string cellName;
        uint alpha = 500;
        uint idx=0;
        double cellXpos;  
        double cellYpos;
        double densityPenaltyGradient;
        double penaltyParameter;
        penaltyParameter = (*((Design*)myDesign)).DesignGetpenaltyParameter();
        densityPenaltyGradient = (*((Design*)myDesign)).DesignComputeTotalDensityPenaltyGradient();
        cout << "Total Density Penalty Gradient: " << densityPenaltyGradient << endl; 
        DESIGN_FOR_ALL_CELLS((*((Design*)myDesign)),cellName,cellPtr){
                if ((*cellPtr).CellIsTerminal()) continue;
                double gradX=0;
                double gradY=0;
                double cellMaxx;
                double cellMinx;
                double cellMaxy;
                double cellMiny;
                cellXpos = (*cellPtr).CellGetXposDbl();
                cellYpos = (*cellPtr).CellGetYposDbl();
                Net *netPtr;
                double sumPinsPos;
                double pinMaxx;
                double pinMinx;
                double pinMaxy;
                double pinMiny;
                double tempDivideX;
                double tempDivideY;
                double pinXpos;
                double pinYpos;
                CELL_FOR_ALL_NETS_NO_DIR((*cellPtr),netPtr){
                        Pin *pinPtr;
                        NET_FOR_ALL_PINS((*netPtr),pinPtr){
                                Cell* cellParentPtr;
                                cellParentPtr = (*pinPtr).PinGetParentCellPtr();
                                pinXpos = pinPtr->xOffset + cellParentPtr->x;
                                pinYpos = pinPtr->yOffset + cellParentPtr->y;
                                tempDivideX= myDivideNew(pinXpos,alpha);
                                tempDivideY= myDivideNew(pinYpos,alpha);
                                pinMaxx += exp(ceil(tempDivideX));
                                pinMinx +=myDivideNew(1,exp(ceil(tempDivideX)));
                                pinMaxy += exp(ceil(tempDivideY));
                                pinMiny +=myDivideNew(1,exp(ceil(tempDivideY)));
                        }NET_END_FOR;
                }CELL_END_FOR;
                double tempDivX;
                double tempDivY;
                tempDivX = myDivideNew(cellXpos,alpha);
                tempDivY = myDivideNew(cellYpos,alpha);
                cellMaxx = exp(ceil(tempDivX));
                cellMinx = myDivideNew(1,exp(ceil(tempDivX)));
                cellMaxy = exp(ceil(tempDivY));
                cellMiny = myDivideNew(1,exp(ceil(tempDivY)));
                double temp1;
                double temp2;
                temp1 = myDivideNew(cellMaxx,pinMaxx);
                temp2 = myDivideNew(cellMinx,pinMinx);
                gradX = (temp1)-(temp2);
                double temp3;
                double temp4;
                temp3 = myDivideNew(cellMaxy,pinMaxy);
                temp4 = myDivideNew(cellMiny,pinMiny);
                gradY = (temp3)-(temp4);
                // Change to add the gradient for the penalty is here
                double gradPotentialX;
                double gradPotentialY;
                (*((Design*)myDesign)).DesignComputePenaltyGradientforCell(cellPtr,gradPotentialX,gradPotentialY);
               // grad[idx] = gradX + penaltyParameter*gradPotentialX*densityPenaltyGradient*alpha;
                //grad[idx+1]=gradY + penaltyParameter*gradPotentialY*densityPenaltyGradient*alpha;
                grad[idx] = gradX*alpha + penaltyParameter*gradPotentialX;
                grad[idx+1]=gradY*alpha + penaltyParameter*gradPotentialY;
                //grad[idx] = gradX + 0.5*gradPotentialX;
                //grad[idx+1]=gradY + 0.5*gradPotentialY;
                //cout << "gradient value for cell: " << (*cellPtr).CellGetName() << "is gradx: " << gradX << " gradY: " <<gradY << endl;
                if (isnan(grad[idx])){
                     cout << "gradX is nan: Find the debug information below:" << endl;
                     cout << "temp1 and temp2 are: "<<temp1 <<"\t" <<temp2<<endl;
                     cout << "cellMaxX and CellMinX are: "<<cellMaxx<<"\t"<<cellMinx<<endl;
                     cout << "pinMaxx and PinMinx are: " <<pinMaxx<<"\t"<<pinMinx<<endl;
                     cout << "cellXpos is: "<< cellXpos<<" cellName: " << cellName<<endl;  
                //     grad[idx]=0.1;
                }
                if (isnan(grad[idx+1])){
                     cout << "gradY is nan: Find the debug information below:" << endl;
                     cout << "temp3 and temp4 are: "<<temp3 <<"\t" <<temp4<<endl;
                     cout << "cellMaxY and CellMinY are: "<<cellMaxy<<"\t"<<cellMiny<<endl;
                     cout << "pinMaxy and PinMiny are: " <<pinMaxy<<"\t"<<pinMiny<<endl;
                     cout << "cellYpos is: "<< cellYpos<<" cellName: " << cellName<<endl;   
                  //   grad[idx+1]=0.1;
                }
                cout <<"iteration Number: "<< gradIterationCount <<"Cell Name: " <<cellName << " CellXpos: "<<cellXpos<< " cellYpos: "<<cellYpos<<" gradX: "<<grad[idx]<<" gradY: "<<grad[idx+1]<<endl;  
                
                idx=idx+2;
        }DESIGN_END_FOR;
        gradIterationCount++;
}


void
Design::DesignSolveForAllCellsWnnlpNew(void)
{
bool debug = true;
//extern int wn_nlp_verbose=3;
//All Variable Declarations for cluster placement
Cell *clusterCellPtr;
vector<Cell *> clusterCells;
double clusterXpos,clusterYpos;
uint maxx,maxy;
uint constMaxx,constMaxy;
ulong lseXHPWL,lseYHPWL;
uint averageClusterWidth,averageClusterHeight;
uint averageStdCellWidth,acerageStdCellHeight; 
uint numClusters,numRows,numSites;
uint siteWidth,rowHeight;

//All variable Declaration for cell placement
void *cellObj;
ofstream logFile;
vector<Cell *> inputCells;
uint numVars,nodeIdx;
uint inputCellCount;
Cell* allCellsptr;
string cellName;

string DesignPath, DesignName;
string DirName;

//All variable declration for nlp solver X and Y component are solved seperately

wn_sll constraint_list;

wn_linear_constraint_type linear_constraint;
wn_nonlinear_constraint_type nonlinear_constraint;

wn_nonlinear_constraint_type objective;

uint i,j;

int code;

double val_min;

// All initialization
Env &DesignEnv = (*this).DesignEnv;
DesignPath = DesignEnv.EnvGetDesignPath();
DirName = DesignPath + "/.solverData";
DesignName = DesignEnv.EnvGetDesignName();

HyperGraph &myGraph = (*this).DesignGetGraph();

HYPERGRAPH_FOR_ALL_NODES(myGraph, nodeIdx, cellObj) {
    if ((*(Cell*)cellObj).CellIsTerminal()) continue;
        inputCells.push_back((Cell *)cellObj);
        cellsToSolve.push_back((Cell *)cellObj);
} HYPERGRAPH_END_FOR;

inputCellCount = inputCells.size(); 

/*Create Placeable blocks in the design*/
DesignGetBoundingBox(maxx,maxy);

averageStdCellWidth = (uint)DesignGetAverageStdCellWidth();
averageStdCellHeight=(uint)DesignGetAverageStdCellHeight();

//averageClusterWidth = (uint)DesignGetAverageClusterCellWidth()+2*averageStdCellWidth;
//averageClusterHeight = (uint)DesignGetAverageClusterCellHeight()+2*averageStdCellHeight;

averageClusterWidth = (uint)DesignGetAverageClusterCellWidth();
averageClusterHeight = (uint)DesignGetAverageClusterCellHeight();

numClusters = DesignGetNumClusters();
numRows = floor(((double)maxy) / averageClusterHeight);
numSites = ceil(((double)numClusters) / numRows);
siteWidth = floor(((double)maxx) / numSites);
rowHeight = averageClusterHeight;
constMaxx = maxx -  averageClusterWidth;
constMaxy = maxy - averageClusterHeight;

uint Xcenterdie,Ycenterdie;
Xcenterdie=maxx/2;
Ycenterdie=maxy/2;

uint siteNum,rowNum;
/* Place Clusters Constructively*/ 
siteNum = 0;
rowNum = 0;
uint count = 0;
/*DESIGN_FOR_ALL_CLUSTERS((*this), cellName, clusterCellPtr) {
        clusterXpos = siteNum * siteWidth;
        clusterYpos = rowNum * rowHeight;
//        clusterXpos = Xcenterdie;
       //clusterYpos = Ycenterdie;
        (*clusterCellPtr).CellSetXpos(clusterXpos); 
        siteNum++;
        (*clusterCellPtr).CellSetYpos(clusterYpos); 
        if (siteNum == numSites) {
             siteNum = 0;
             rowNum++;
        }   
        clusterCells.push_back(clusterCellPtr);
        count++;
} DESIGN_END_FOR;
*/


/* rameshul Debug*/
if (debug) { 
        ulong lseHPWL;
        cout << "maxx " << maxx << " maxy " << maxy << endl;
        cout << "avgCluster Width " << averageClusterWidth << " avg Cluster Height" << averageClusterHeight << endl;
        /* For loop to display cell information of all input cells*/
      // printAllVisibleCellsInDesign((*this),"nlp_cells"); 
       lseHPWL=DesignComputeLseHPWL();
       cout << "LSE HPWL of Visible cells is : " << lseHPWL << endl;
      // printVisibleCellsineachCluster((*this),"nlp_cluster");
   double percentageOverlap;
   percentageOverlap = (*this).DesignDumpClusterOverlapForPenalty();
   cout << "Percentage Overlap for debug" << percentageOverlap << endl; 
   
}
/*End Debug*/


/* Get all the variables to be solved for from the design class*/

vector<Cell *> cellsToSolve;
getNLPCellsToSolveNew ((*this),cellsToSolve);
numVars= 2*(cellsToSolve.size());
/*rameshul Debug -print the cells returned from the above function  

Cell * cellsToSolvePtr; 
string cellToSolveName;
VECTOR_FOR_ALL_ELEMS(cellsToSolve,Cell *,cellsToSolvePtr){
        cellToSolveName=(*cellsToSolvePtr).CellGetName();
        cout << "Cell Name to Solve: " << cellToSolveName << " Xpos: " << (*cellsToSolvePtr).CellGetXpos() << " Ypos: " << (*cellsToSolvePtr).CellGetYpos() <<endl;
}END_FOR;

End Debug*/

double values[numVars];
/* Populate Initial Values of x and y*/
nlpInitialValuesNew(cellsToSolve,values);
DesignComputepenaltyParameter();
/* rameshul Debug for objective Function
ulong wirelengthBeforeX;
wirelengthBeforeX = (*this).DesignComputeLseXHPWL();
cout << "Wirelength X  before changing  is: " << wirelengthBeforeX << endl;
ulong wirelengthBeforeY;
wirelengthBeforeY = (*this).DesignComputeLseYHPWL();
cout << "Wirelength Y  before changing  is: " << wirelengthBeforeY << endl;
double wirelengthTempX;
wirelengthTempX = wirelengthObjFunc(numVars,values,this);  // was just a check to call function
cout << "Temporary X  Wirelength is: " << wirelengthTempX << endl; 
double wirelengthTempY;
wirelengthTempY = wirelengthObjFuncY(numVars,y,this);  // was just a check to call function
cout << "Temporary Y  Wirelength is: " << wirelengthTempY << endl; 
 End Debug */

/*rameshul debug for gradient function 

double gradientX[numVars];
double gradientY[numVars];

gradientFuncX(gradientX,numVars,x,this);        
gradientFuncY(gradientY,numVars,y,this);        

rameshul End debug */


/****************** Declaration of Objective for X and Y variables *****************************/

wn_make_nonlinear_constraint(&nonlinear_constraint,numVars,WN_EQ_COMPARISON);
nonlinear_constraint->pfunction = &wirelengthObjFunc;
nonlinear_constraint->pgradient = &gradientFunc;
for(i=0;i<numVars;++i)
  {
    (nonlinear_constraint->vars)[i] = i;
  }
nonlinear_constraint->client_data = this;
objective = nonlinear_constraint;

/*wn_make_nonlinear_constraint(&nonlinear_constraintY,numVars,WN_EQ_COMPARISON);
nonlinear_constraintY->pfunction = &wirelengthObjFuncY;
//nonlinear_constraintY->pgradient = &gradientFuncY;
for(iY=0;iY<numVars;++iY)
  {
    (nonlinear_constraintY->vars)[iY] = iY;
  }
nonlinear_constraintY->client_data = this;
objectiveY = nonlinear_constraintY;*


/***************** Declaration of constraints for solving x and y part of the nonlinear problem **************/
constraint_list = NULL;
//constraint 1 - set minimum value of x and y
for (i=0;i<numVars;++i){
        wn_make_linear_constraint(&linear_constraint,1,0.0,WN_GT_COMPARISON);
        (linear_constraint->vars)[0] = i;
        (linear_constraint->weights)[0] = 1.0;
        wn_sllins(&constraint_list,linear_constraint);
}
//constraint 2 - set maximum value of x

for (i=0;i<numVars;i=i+2){
        wn_make_linear_constraint(&linear_constraint,1,constMaxx,WN_LT_COMPARISON);
        (linear_constraint->vars)[0] = i;
        (linear_constraint->weights)[0] = 1.0;
        wn_sllins(&constraint_list,linear_constraint);
}


//constraint 2 - set maximum value of y

for (i=1;i<numVars;i=i+2){
        wn_make_linear_constraint(&linear_constraint,1,constMaxy,WN_LT_COMPARISON);
        (linear_constraint->vars)[0] = i;
        (linear_constraint->weights)[0] = 1.0;
        wn_sllins(&constraint_list,linear_constraint);
}

/*  Trying below the overlap constraints based on LP approach as suggested by Prof.Jeffery Camm */

/*for (i=0;i<numVars;i=i+2){
        for (j=0;j<numVars;j=j+2){
                if (j <= i) continue;
                wn_make_linear_constraint(&linear_constraint,1,(averageClusterWidth+values[j]),WN_GT_COMPARISON);
                (linear_constraint->vars)[0] = i;
                (linear_constraint->weights)[0] = 1.0;
               wn_sllins(&constraint_list,linear_constraint);
                wn_make_linear_constraint(&linear_constraint,1,(-averageClusterWidth-values[j]),WN_GT_COMPARISON);
                (linear_constraint->vars)[0] = i;
                (linear_constraint->weights)[0] = 1.0;
               wn_sllins(&constraint_list,linear_constraint);
        }
}
for (i=1;i<numVars;i=i+2){
        for (j=1;j<numVars;j=j+2){
                if (j <= i) continue;
                wn_make_linear_constraint(&linear_constraint,1,(averageClusterHeight+values[j]),WN_GT_COMPARISON);
                (linear_constraint->vars)[0] = i;
                (linear_constraint->weights)[0] = 1.0;
               wn_sllins(&constraint_list,linear_constraint);
                wn_make_linear_constraint(&linear_constraint,1,(-averageClusterHeight-values[j]),WN_GT_COMPARISON);
                (linear_constraint->vars)[0] = i;
                (linear_constraint->weights)[0] = 1.0;
               wn_sllins(&constraint_list,linear_constraint);
        }
}
*/





/***** JUST a try with the deltaX ****/

//double *delta=NULL;
uint numDeltas=2*numVars;
double delta[numDeltas];

// wn_make_vect(&deltaX,numVars);
for(i=0;i<numDeltas;++i)
 {
   delta[i] = 0.0001;
 }



/**************** Call the conjugate gradient minimizer ******************************/

// For X coordinates
//Commented for compiling
wn_nlp_conj_method(&code,&val_min,values,delta,(wn_nonlinear_constraint_type)objective,constraint_list,numVars,numVars,10,1);

cout << " Final value after minimization is: " << val_min << " code: " <<code <<endl;   
cout << "Final HPWL after minimization is: " <<  (*this).DesignComputeLseHPWL()<<endl;

for(int i=0;i<numVars;++i){

        cout << "Final X: " <<values[i]<<" Final Y: "<<values[i+1]<<endl;
        i=i+1;
}

   double percentageOverlap;
   percentageOverlap = (*this).DesignDumpClusterOverlapForPenalty();
   cout << "Percentage Overlap for debug" << percentageOverlap << endl; 

}
