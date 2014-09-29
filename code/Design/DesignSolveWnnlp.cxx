#include <Design.h>

#include <math.h>
#include "wnlib/wnlib.h"
#include "wnlib/wnasrt.h"
#include "wnlib/wnsll.h"
#include "wnlib/wnnlp.h"
#include <typeinfo>

void
getNLPCellsToSolve(Design &myDesign,vector<Cell *> &cellsToSolve)
{

Cell *cellPtr;
string cellName;

DESIGN_FOR_ALL_CELLS(myDesign,cellName,cellPtr){
        cellsToSolve.push_back(cellPtr);
}DESIGN_END_FOR;

}

void
nlpInitialValues(vector<Cell *> &cellsToSolve, double* x, double* y)
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
        y[idx] = (*cellPtr).CellGetYposDbl();
  //      cellToSolveMap[idx] = CellPtr;
       cout << "\t Xpos: " << x[idx] << "\t Ypos: " << y[idx] <<endl;
        idx = idx+1;
}END_FOR;

//cout << "I value at the end of iteration is: " << idx << endl;

}

double
wirelengthObjFuncX(int size, double *values,ptr myDesign) /* Need to justify size, values and client_data. Just reused it from definition of pfunction*/
{
  /* objective here is total_wirelength = sum over all nets ( wirelength (net) * net weight) */ 
  Cell *cellPtr;
  string cellName;
  uint idx=0;
        //cout << " Void pointer Type: " <<typeid(*((Design*)myDesign)).name() << endl;  
        DESIGN_FOR_ALL_CELLS((*((Design*)myDesign)),cellName,cellPtr){
//                cout << "Changing Xposition for cell: " << (*cellPtr).CellGetName() << " to " << values[idx] << endl; 
                (*cellPtr).CellSetXposDbl(values[idx]);
  //              cout << "The new Xposition of the cell is " << (*cellPtr).CellGetXpos() << endl;
                idx = idx+1;

        }END_FOR;
        ulong xHPWL;
        xHPWL = (*((Design*)myDesign)).DesignComputeLseXHPWL();
        double rtv;
        rtv = xHPWL;
        //cout << "ulong wirelength is: " << xHPWL <<endl;
        return rtv;
}
        

double
wirelengthObjFuncY(int size, double *values,ptr myDesign) /* Need to justify size, values and client_data. Just reused it from definition of pfunction*/
{
  /* objective here is total_wirelength = sum over all nets ( wirelength (net) * net weight) */ 
  Cell *cellPtr;
  string cellName;
  uint idx=0;
        //cout << " Void pointer Type: " <<typeid(*((Design*)myDesign)).name() << endl;  
        DESIGN_FOR_ALL_CELLS((*((Design*)myDesign)),cellName,cellPtr){
//                cout << "Changing Yposition for cell: " << (*cellPtr).CellGetName() << " to " << values[idx] << endl; 
                (*cellPtr).CellSetYposDbl(values[idx]);
  //              cout << "The new Yposition of the cell is " << (*cellPtr).CellGetYpos() << endl;
                idx = idx+1;

        }END_FOR;
        ulong yHPWL;
        yHPWL = (*((Design*)myDesign)).DesignComputeLseYHPWL();
        double rtv;
        rtv = yHPWL;
        //cout << "ulong wirelength is: " << yHPWL <<endl;
        return rtv;
}

void 
gradientFuncX(double *grad,int size,double *values,ptr myDesign)
{

        Cell *cellPtr;
        string cellName;
        uint alpha = 500;
        uint idx=0;
        uint cellXpos;  
        DESIGN_FOR_ALL_CELLS((*((Design*)myDesign)),cellName,cellPtr){
                double gradX=0;
                double cellMaxx;
                double cellMinx;
                cellXpos = (*cellPtr).CellGetXpos();
                Net *netPtr;
                double sumPinsPos;
                double pinMaxx;
                double pinMinx;
                uint pinXpos;
                CELL_FOR_ALL_NETS_NO_DIR((*cellPtr),netPtr){
                        Pin *pinPtr;
                        NET_FOR_ALL_PINS((*netPtr),pinPtr){
                                pinXpos = pinPtr->xOffset + cellXpos;
                                pinMaxx += exp(pinXpos/alpha);
                                pinMinx +=1/(exp(pinXpos/alpha));
                        }NET_END_FOR;
                }CELL_END_FOR;
                cellMaxx = exp(cellXpos/alpha);
                cellMinx = 1/(exp(cellXpos/alpha));
                gradX = (cellMaxx/pinMaxx)-(cellMinx/pinMinx);
                grad[idx] = gradX;
                //cout << "gradient value for cell: " << (*cellPtr).CellGetName() << "is" << gradX << endl;
                idx=idx+1;
        }DESIGN_END_FOR;
}


void 
gradientFuncY(double *grad,int size,double *values,ptr myDesign)
{
        Cell *cellPtr;
        string cellName;
        uint alpha = 500;
        uint idx=0;
        uint cellYpos;  
        DESIGN_FOR_ALL_CELLS((*((Design*)myDesign)),cellName,cellPtr){
                double gradY=0;
                double cellMaxy;
                double cellMiny;
                cellYpos = (*cellPtr).CellGetYpos();
                Net *netPtr;
                double sumPinsPos;
                double pinMaxy;
                double pinMiny;
                uint pinYpos;
                CELL_FOR_ALL_NETS_NO_DIR((*cellPtr),netPtr){
                        Pin *pinPtr;
                        NET_FOR_ALL_PINS((*netPtr),pinPtr){
                                pinYpos = pinPtr->yOffset + cellYpos;
                                pinMaxy += exp(pinYpos/alpha);
                                pinMiny +=1/(exp(pinYpos/alpha));
                        }NET_END_FOR;
                }CELL_END_FOR;
                cellMaxy = exp(cellYpos/alpha);
                cellMiny = 1/(exp(cellYpos/alpha));
                gradY = (cellMaxy/pinMaxy)-(cellMiny/pinMiny);
                grad[idx] = gradY;
                //cout << "Y gradient value for cell: " << (*cellPtr).CellGetName() << "is" << gradY << endl;
                idx=idx+1;
        }DESIGN_END_FOR;
}

void
Design::DesignSolveForAllCellsWnnlp(void)
{
bool debug = true;
extern int wn_nlp_verbose=3;
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


wn_sll constraint_listX;
wn_sll constraint_listY;

wn_linear_constraint_type linear_constraintX,linear_constraintY;
wn_nonlinear_constraint_type nonlinear_constraintX,nonlinear_constraintY;

wn_nonlinear_constraint_type objectiveX;
wn_nonlinear_constraint_type objectiveY;

uint iX;
uint iY;

int codeX;
int codeY;

double val_minX;
double val_minY;


double *deltaX=NULL;
double *deltaY=NULL;


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
if (!debug) { 
        ulong lseHPWL;
        cout << "maxx " << maxx << " maxy " << maxy << endl;
        cout << "avgCluster Width " << averageClusterWidth << " avg Cluster Height" << averageClusterHeight << endl;
        /* For loop to display cell information of all input cells*/
       printAllVisibleCellsInDesign((*this),"nlp_cells"); 
       lseHPWL=DesignComputeLseHPWL();
       cout << "LSE HPWL of Visible cells is : " << lseHPWL << endl;
       printVisibleCellsineachCluster((*this),"nlp_cluster");
}
/*End Debug*/


/* Get all the variables to be solved for from the design class*/

vector<Cell *> cellsToSolve;
getNLPCellsToSolve ((*this),cellsToSolve);
numVars= cellsToSolve.size();
/*rameshul Debug -print the cells returned from the above function  

Cell * cellsToSolvePtr; 
string cellToSolveName;
VECTOR_FOR_ALL_ELEMS(cellsToSolve,Cell *,cellsToSolvePtr){
        cellToSolveName=(*cellsToSolvePtr).CellGetName();
        cout << "Cell Name to Solve: " << cellToSolveName << " Xpos: " << (*cellsToSolvePtr).CellGetXpos() << " Ypos: " << (*cellsToSolvePtr).CellGetYpos() <<endl;
}END_FOR;

End Debug*/

double x[numVars];
double y[numVars];
/* Populate Initial Values of x and y*/
nlpInitialValues(cellsToSolve,x,y);

/* rameshul Debug for objective Functioni*/
ulong wirelengthBeforeX;
wirelengthBeforeX = (*this).DesignComputeLseXHPWL();
cout << "Wirelength X  before changing  is: " << wirelengthBeforeX << endl;
ulong wirelengthBeforeY;
wirelengthBeforeY = (*this).DesignComputeLseYHPWL();
cout << "Wirelength Y  before changing  is: " << wirelengthBeforeY << endl;
/*double wirelengthTempX;
wirelengthTempX = wirelengthObjFuncX(numVars,x,this);  // was just a check to call function
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

wn_make_nonlinear_constraint(&nonlinear_constraintX,numVars,WN_EQ_COMPARISON);
nonlinear_constraintX->pfunction = &wirelengthObjFuncX;
nonlinear_constraintX->pgradient = &gradientFuncX;
for(iX=0;iX<numVars;++iX)
  {
    (nonlinear_constraintX->vars)[iX] = iX;
  }
nonlinear_constraintX->client_data = this;
objectiveX = nonlinear_constraintX;

wn_make_nonlinear_constraint(&nonlinear_constraintY,numVars,WN_EQ_COMPARISON);
nonlinear_constraintY->pfunction = &wirelengthObjFuncY;
nonlinear_constraintY->pgradient = &gradientFuncY;
for(iY=0;iY<numVars;++iY)
  {
    (nonlinear_constraintY->vars)[iY] = iY;
  }
nonlinear_constraintY->client_data = this;
objectiveY = nonlinear_constraintY;


/***************** Declaration of constraints for solving x and y part of the nonlinear problem **************/
constraint_listX = NULL;
//constraint 1 - set minimum value of x
for (iX=0;iX<numVars;++iX){
        wn_make_linear_constraint(&linear_constraintX,1,0.0,WN_GT_COMPARISON);
        (linear_constraintX->vars)[0] = iX;
        (linear_constraintX->weights)[0] = 1.0;
        wn_sllins(&constraint_listX,linear_constraintX);
}
//constraint 2 - set maximum value of x

for (iX=0;iX<numVars;++iX){
        wn_make_linear_constraint(&linear_constraintX,1,maxx,WN_LT_COMPARISON);
        (linear_constraintX->vars)[0] = iX;
        (linear_constraintX->weights)[0] = 1.0;
        wn_sllins(&constraint_listX,linear_constraintX);
}

constraint_listY = NULL;
//constraint 1 - set minimum value of y
for (iY=0;iY<numVars;++iY){
        wn_make_linear_constraint(&linear_constraintY,1,0.0,WN_GT_COMPARISON);
        (linear_constraintY->vars)[0] = iY;
        (linear_constraintY->weights)[0] = 1.0;
        wn_sllins(&constraint_listY,linear_constraintY);
}
//constraint 2 - set maximum value of y

for (iY=0;iY<numVars;++iY){
        wn_make_linear_constraint(&linear_constraintY,1,maxy,WN_LT_COMPARISON);
        (linear_constraintY->vars)[0] = iY;
        (linear_constraintY->weights)[0] = 1.0;
        wn_sllins(&constraint_listY,linear_constraintY);
}



/***** JUST a try with the deltaX ****/

/* wn_make_vect(&deltaX,numVars);
for(ix=0;ix<numVars;++ix)
 {
   deltaX[ix] = 0.0001;
 }

 wn_make_vect(&deltaY,numVars);
for(iy=0;iy<numVars;++iy)
 {
   deltaY[iy] = 0.0001;
 }*/



            



/**************** Call the conjugate gradient minimizer ******************************/

// For X coordinates
//Commented for compiling
wn_nlp_conj_method(&codeX,&val_minX,x,deltaX,(wn_nonlinear_constraint_type)objectiveX,constraint_listX,numVars,numVars,10,1.0);

//For Y coordinates
//commented for compiling
wn_nlp_conj_method(&codeY,&val_minY,y,deltaY,(wn_nonlinear_constraint_type)objectiveY,constraint_listY,numVars,numVars,10,1.0);

cout << " Final X value after minimization is: " << val_minX << "\t Final Y value: " << val_minY << " codeX: " <<codeX <<" codeY: "<<codeY<<endl;   


for(int i=0;i<numVars;++i){

        cout << "Final X: " <<x[i]<<" Final Y: "<<y[i]<<endl;
}

}
