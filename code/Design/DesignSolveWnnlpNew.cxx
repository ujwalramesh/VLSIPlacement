#include <Design.h>

#include <math.h>
#include "wnlib/wnlib.h"
#include "wnlib/wnasrt.h"
#include "wnlib/wnsll.h"
#include "wnlib/wnnlp.h"
#include <typeinfo>

local double myDivideNew(double num, double denom) 
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
//cout << "Size of cells: " << size << endl;
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


double
wirelengthObjFunc(int size, double *values,ptr myDesign) /* Need to justify size, values and client_data. Just reused it from definition of pfunction*/
{
  /* objective here is total_wirelength = sum over all nets ( wirelength (net) * net weight) */ 
  Cell *cellPtr;
  string cellName;
  uint idx=0;
        //cout << " Void pointer Type: " <<typeid(*((Design*)myDesign)).name() << endl;  
        DESIGN_FOR_ALL_CELLS((*((Design*)myDesign)),cellName,cellPtr){
//                cout << "Changing Xposition for cell: " << (*cellPtr).CellGetName() << " to " << values[idx] << endl; 
                (*cellPtr).CellSetXposDbl(values[idx]);
                (*cellPtr).CellSetYposDbl(values[idx+1]);
  //              cout << "The new Xposition of the cell is " << (*cellPtr).CellGetXpos() << endl;
                idx = idx+2;

        }END_FOR;
        ulong LseHPWL;
        LseHPWL = (*((Design*)myDesign)).DesignComputeLseHPWL();
   /*Adding penalty function for density constraint as per Will Naylor's patent*/
    double totalDensityPenalty=0;
    (*((Design*)myDesign)).DesignUpdateGridPotentials();
    totalDensityPenalty =(*((Design*)myDesign)).DesignComputeTotalDensityPenalty();
        double rtv;
        rtv = LseHPWL+totalDensityPenalty;
   //     cout << "ulong wirelength is: " << LseHPWL <<endl;
    //    cout << "Density Penalty is : " << totalDensityPenalty<<endl;
        return rtv;
}
        
void 
gradientFunc(double *grad,int size,double *values,ptr myDesign)
{

        Cell *cellPtr;
        string cellName;
        uint alpha = 500;
        uint idx=0;
        uint cellXpos;  
        uint cellYpos;
        DESIGN_FOR_ALL_CELLS((*((Design*)myDesign)),cellName,cellPtr){
                double gradX=0;
                double gradY=0;
                double cellMaxx;
                double cellMinx;
                double cellMaxy;
                double cellMiny;
                cellXpos = (*cellPtr).CellGetXpos();
                cellYpos = (*cellPtr).CellGetYpos();
                Net *netPtr;
                double sumPinsPos;
                double pinMaxx;
                double pinMinx;
                double pinMaxy;
                double pinMiny;
                double tempDivideX;
                double tempDivideY;
                uint pinXpos;
                uint pinYpos;
                CELL_FOR_ALL_NETS_NO_DIR((*cellPtr),netPtr){
                        Pin *pinPtr;
                        NET_FOR_ALL_PINS((*netPtr),pinPtr){
                                pinXpos = pinPtr->xOffset + cellXpos;
                                pinYpos = pinPtr->yOffset + cellYpos;
                                tempDivideX= myDivideNew(pinXpos,alpha);
                                tempDivideY= myDivideNew(pinYpos,alpha);
                                pinMaxx += exp(tempDivideX);
                                pinMinx +=1/(exp(tempDivideX));
                                pinMaxy += exp(tempDivideY);
                                pinMiny +=1/(exp(tempDivideY));
                        }NET_END_FOR;
                }CELL_END_FOR;
                double tempDivX;
                double tempDivY;
                tempDivX = myDivideNew(cellXpos,alpha);
                tempDivY = myDivideNew(cellYpos,alpha);
                cellMaxx = exp(tempDivX);
                cellMinx = 1/(exp(tempDivX));
                cellMaxy = exp(tempDivY);
                cellMiny = 1/(exp(tempDivY));
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
                //cout << "gradient value for cell: " << (*cellPtr).CellGetName() << "is" << gradX << endl;
                idx=idx+2;
        }DESIGN_END_FOR;
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
ulong lseXHPWL,lseYHPWL;
uint averageClusterWidth,averageClusterHeight;
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

uint i;

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
averageClusterWidth = (uint)DesignGetAverageClusterCellWidth();
averageClusterHeight = (uint)DesignGetAverageClusterCellHeight();
numClusters = DesignGetNumClusters();
numRows = floor(((double)maxy) / averageClusterHeight);
numSites = ceil(((double)numClusters) / numRows);
siteWidth = floor(((double)maxx) / numSites);
rowHeight = averageClusterHeight;

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
        wn_make_linear_constraint(&linear_constraint,1,maxx,WN_LT_COMPARISON);
        (linear_constraint->vars)[0] = i;
        (linear_constraint->weights)[0] = 1.0;
        wn_sllins(&constraint_list,linear_constraint);
}


//constraint 2 - set maximum value of y

for (i=1;i<numVars;i=i+2){
        wn_make_linear_constraint(&linear_constraint,1,maxy,WN_LT_COMPARISON);
        (linear_constraint->vars)[0] = i;
        (linear_constraint->weights)[0] = 1.0;
        wn_sllins(&constraint_list,linear_constraint);
}



/***** JUST a try with the deltaX ****/

double *delta=NULL;
/*double delta[numVars];

// wn_make_vect(&deltaX,numVars);
for(i=0;i<numVars;++i)
 {
   delta[i] = 0.0001;
 }*/



/**************** Call the conjugate gradient minimizer ******************************/

// For X coordinates
//Commented for compiling
wn_nlp_conj_method(&code,&val_min,values,delta,(wn_nonlinear_constraint_type)objective,constraint_list,numVars,numVars,10,1.0);

cout << " Final value after minimization is: " << val_min << " code: " <<code <<endl;   
cout << "Final HPWL after minimization is: " <<  (*this).DesignComputeLseHPWL()<<endl;

for(int i=0;i<numVars;++i){

        cout << "Final X: " <<values[i]<<" Final Y: "<<values[i+1]<<endl;
        i=i+1;
}

}
