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
nlpInitialValuesNew(vector<Cell *> &cellsToSolve, double* x,double* y)
{

        Cell *cellPtr;
        string cellName;
        uint idx=0;

        int size = cellsToSolve.size();
        cout << "Size of cells: " << size << endl;

        VECTOR_FOR_ALL_ELEMS(cellsToSolve,Cell *, cellPtr){
                x[idx] = (*cellPtr).CellGetXposDbl();
                y[idx] = (*cellPtr).CellGetYposDbl();
                cout << "\t Xpos: " << x[idx] << "\t Ypos: " << y[idx] << "\t actualXpos: " << (*cellPtr).CellGetXposDbl() << "\t actualYpos: ";
                cout << (*cellPtr).CellGetYposDbl()<< endl;
                idx = idx+1;
        }END_FOR;

}

double
wirelengthObjFuncX(int size, double *values,ptr myDesign) 
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
                //(*cellPtr).CellSetYposDbl(values[idx+1]);
                //cout << "The new Xposition of the cell is " << (*cellPtr).CellGetXpos() << endl;
                idx = idx+1;

        }END_FOR;
        ulong LseHPWL;
        LseHPWL = (*((Design*)myDesign)).DesignComputeLseHPWL();
   /*Adding penalty function for density constraint as per Will Naylor's patent*/
   double totalDensityPenalty=0; /*variable for adding penalty for my objective function*/

/*  The below penalty method is for penalizing my objective based on grod point potentials */
   (*((Design*)myDesign)).DesignUpdateGridPotentials();
   totalDensityPenalty =(*((Design*)myDesign)).DesignComputeTotalDensityPenalty();
    if (penaltyParameter < 0 ) {
            penaltyParameter = -penaltyParameter;
   }
   //static uint lambda;
   //lambda=2;
        double rtv;
       rtv = LseHPWL+(penaltyParameter*totalDensityPenalty);
        cout << "ulong wirelength X is: " << LseHPWL << "\t";
        cout << "Density Penalty X is : " << penaltyParameter*totalDensityPenalty << "\t";
        cout << "objective Value X is : " << rtv << endl;
        //lambda=lambda*2;
    objIterationCount++;
        return rtv;
}

double
wirelengthObjFuncY(int size, double *values,ptr myDesign) 
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
                //(*cellPtr).CellSetXposDbl(values[idx]);
                (*cellPtr).CellSetYposDbl(values[idx]);
                //cout << "The new Xposition of the cell is " << (*cellPtr).CellGetXpos() << endl;
                idx = idx+1;

        }END_FOR;
        ulong LseHPWL;
        LseHPWL = (*((Design*)myDesign)).DesignComputeLseHPWL();
   /*Adding penalty function for density constraint as per Will Naylor's patent*/
   double totalDensityPenalty=0; /*variable for adding penalty for my objective function*/

/*  The below penalty method is for penalizing my objective based on grod point potentials */
   (*((Design*)myDesign)).DesignUpdateGridPotentials();
   totalDensityPenalty =(*((Design*)myDesign)).DesignComputeTotalDensityPenalty();
    if (penaltyParameter < 0 ) {
            penaltyParameter = -penaltyParameter;
   }
   //static uint lambda;
   //lambda=2;
        double rtv;
       rtv = LseHPWL+(penaltyParameter*totalDensityPenalty);
        cout << "ulong wirelength Y is: " << LseHPWL << "\t";
        cout << "Density Penalty Y is : " << penaltyParameter*totalDensityPenalty << "\t";
        cout << "objective Value Y is : " << rtv << endl;
        //lambda=lambda*2;
    objIterationCount++;
        return rtv;
}

void
gradientFuncX(double *grad,int size,double *values,ptr myDesign)
{
        static int gradIterationCount=0;
        Cell *cellPtr;
        string cellName;
        uint alpha = 500;
        uint idx=0;
        double cellXpos;
        //double cellYpos;
        double densityPenaltyGradient;
        double penaltyParameter;
        penaltyParameter = (*((Design*)myDesign)).DesignGetpenaltyParameter();
        densityPenaltyGradient = (*((Design*)myDesign)).DesignComputeTotalDensityPenaltyGradient();
        cout << "Total Density Penalty Gradient: " << densityPenaltyGradient << endl;
        DESIGN_FOR_ALL_CELLS((*((Design*)myDesign)),cellName,cellPtr){
                if ((*cellPtr).CellIsTerminal()) continue;
                double gradX=0;
                //double gradY=0;
                double cellMaxx;
                double cellMinx;
                //double cellMaxy;
                //double cellMiny;
                cellXpos = (*cellPtr).CellGetXposDbl();
                //cellYpos = (*cellPtr).CellGetYposDbl();
                Net *netPtr;
                double sumPinsPos;
                double pinMaxx;
                double pinMinx;
                //double pinMaxy;
                //double pinMiny;
                 double tempDivideX;
               // double tempDivideY;
                double pinXpos;
                //double pinYpos;
                CELL_FOR_ALL_NETS_NO_DIR((*cellPtr),netPtr){
                        Pin *pinPtr;
                        NET_FOR_ALL_PINS((*netPtr),pinPtr){
                                Cell* cellParentPtr;
                                cellParentPtr = (*pinPtr).PinGetParentCellPtr();
                                pinXpos = pinPtr->xOffset + cellParentPtr->x;
                                //pinYpos = pinPtr->yOffset + cellParentPtr->y;
                                tempDivideX= myDivideNew(pinXpos,alpha);
                                //tempDivideY= myDivideNew(pinYpos,alpha);
                                pinMaxx += exp(tempDivideX);
                                pinMinx +=myDivideNew(1,exp(tempDivideX));
                                //pinMaxy += exp(ceil(tempDivideY));
                                //pinMiny +=myDivideNew(1,exp(ceil(tempDivideY)));
                        }NET_END_FOR;
                }CELL_END_FOR;
                double tempDivX;
                //double tempDivY;
                tempDivX = myDivideNew(cellXpos,alpha);
                //tempDivY = myDivideNew(cellYpos,alpha);
                cellMaxx = exp(tempDivX);
                cellMinx = myDivideNew(1,exp(ceil(tempDivX)));
             //   cellMaxy = exp(ceil(tempDivY));
              //  cellMiny = myDivideNew(1,exp(ceil(tempDivY)));
                double temp1;
                double temp2;
                temp1 = myDivideNew(cellMaxx,pinMaxx);
                temp2 = myDivideNew(cellMinx,pinMinx);
                gradX = (temp1)-(temp2);
             //   double temp3;
             //   double temp4;
             //   temp3 = myDivideNew(cellMaxy,pinMaxy);
             //   temp4 = myDivideNew(cellMiny,pinMiny);
             //   gradY = (temp3)-(temp4);
                // Change to add the gradient for the penalty is here
                double gradPotentialX;
              //  double gradPotentialY;
                (*((Design*)myDesign)).DesignComputePenaltyGradientforCellX(cellPtr,gradPotentialX);
                grad[idx] = gradX + penaltyParameter*gradPotentialX*densityPenaltyGradient;
                //grad[idx] = gradX + penaltyParameter*gradPotentialX;
                //grad[idx+1]=gradY + penaltyParameter*gradPotentialY*densityPenaltyGradient;
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
        cout <<"iteration Number: "<< gradIterationCount <<"Cell Name: " <<cellName << " CellXpos: "<<cellXpos<< " gradX: "<<grad[idx] << endl;

                idx=idx+1;
        }DESIGN_END_FOR;
        gradIterationCount++;
}


                           
void
gradientFuncY(double *grad,int size,double *values,ptr myDesign)
{
        static int gradIterationCount=0;
        Cell *cellPtr;
        string cellName;
        uint alpha = 500;
        uint idx=0;
        //double cellXpos;
        double cellYpos;
        double densityPenaltyGradient;
        double penaltyParameter;
        penaltyParameter = (*((Design*)myDesign)).DesignGetpenaltyParameter();
        densityPenaltyGradient = (*((Design*)myDesign)).DesignComputeTotalDensityPenaltyGradient();
        cout << "Total Density Penalty Gradient: " << densityPenaltyGradient << endl;
        DESIGN_FOR_ALL_CELLS((*((Design*)myDesign)),cellName,cellPtr){
                if ((*cellPtr).CellIsTerminal()) continue;
               // double gradX=0;
                double gradY=0;
                //double cellMaxx;
                //double cellMinx;
                double cellMaxy;
                double cellMiny;
                //cellXpos = (*cellPtr).CellGetXposDbl();
                cellYpos = (*cellPtr).CellGetYposDbl();
                Net *netPtr;
                double sumPinsPos;
                //double pinMaxx;
                //double pinMinx;
                double pinMaxy;
                double pinMiny;
                // double tempDivideX;
                double tempDivideY;
               // double pinXpos;
                double pinYpos;
                CELL_FOR_ALL_NETS_NO_DIR((*cellPtr),netPtr){
                        Pin *pinPtr;
                        NET_FOR_ALL_PINS((*netPtr),pinPtr){
                                Cell* cellParentPtr;
                                cellParentPtr = (*pinPtr).PinGetParentCellPtr();
                               // pinXpos = pinPtr->xOffset + cellParentPtr->x;
                                pinYpos = pinPtr->yOffset + cellParentPtr->y;
                                //tempDivideX= myDivideNew(pinXpos,alpha);
                                tempDivideY= myDivideNew(pinYpos,alpha);
                                //pinMaxx += exp(tempDivideX);
                                //pinMinx +=myDivideNew(1,exp(tempDivideX));
                                pinMaxy += exp(tempDivideY);
                                pinMiny +=myDivideNew(1,exp(tempDivideY));
                        }NET_END_FOR;
                }CELL_END_FOR;
                //double tempDivX;
                double tempDivY;
                //tempDivX = myDivideNew(cellXpos,alpha);
                tempDivY = myDivideNew(cellYpos,alpha);
               // cellMaxx = exp(tempDivX);
               // cellMinx = myDivideNew(1,exp(ceil(tempDivX)));
                cellMaxy = exp(tempDivY);
                cellMiny = myDivideNew(1,exp(tempDivY));
                //double temp1;
                //double temp2;
                //temp1 = myDivideNew(cellMaxx,pinMaxx);
                //temp2 = myDivideNew(cellMinx,pinMinx);
                //gradX = (temp1)-(temp2);
                double temp3;
                double temp4;
                temp3 = myDivideNew(cellMaxy,pinMaxy);
                temp4 = myDivideNew(cellMiny,pinMiny);
                gradY = (temp3)-(temp4);
                // Change to add the gradient for the penalty is here
             //   double gradPotentialX;
                double gradPotentialY;
                (*((Design*)myDesign)).DesignComputePenaltyGradientforCellY(cellPtr,gradPotentialY);
              //  grad[idx] = gradX + penaltyParameter*gradPotentialX*densityPenaltyGradient;
                grad[idx]=gradY + penaltyParameter*gradPotentialY*densityPenaltyGradient;
               // grad[idx]=gradY + penaltyParameter*gradPotentialY;
                //grad[idx] = gradX + 0.5*gradPotentialX;
                //grad[idx+1]=gradY + 0.5*gradPotentialY;
                //cout << "gradient value for cell: " << (*cellPtr).CellGetName() << "is gradx: " << gradX << " gradY: " <<gradY << endl;
                if (isnan(grad[idx])){
                     cout << "gradY is nan: Find the debug information below:" << endl;
                     cout << "temp3 and temp4 are: "<<temp3 <<"\t" <<temp4<<endl;
                     cout << "cellMaxy and CellMiny are: "<<cellMaxy<<"\t"<<cellMiny<<endl;
                     cout << "pinMaxy and PinMiny are: " <<pinMaxy<<"\t"<<pinMiny<<endl;
                     cout << "cellypos is: "<< cellYpos<<" cellName: " << cellName<<endl;
                //     grad[idx]=0.1;
                }
        cout <<"iteration Number: "<< gradIterationCount <<"Cell Name: " <<cellName << " CellYpos: "<<cellYpos<< " gradY: "<<grad[idx]<<endl;
               idx=idx+1;
        }DESIGN_END_FOR;
        gradIterationCount++;
}



void
Design::DesignSolveForAllCellsWnnlpNew(void)
{
        bool debug = true;
        Cell *clusterCellPtr;
        vector<Cell *> clusterCells;
        //All Variable Declarations for cluster placement
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
        //siteWidth =3000;
        //rowHeight = 3000;
        /*DESIGN_FOR_ALL_CLUSTERS((*this), cellName, clusterCellPtr) {
                clusterXpos = siteNum * siteWidth;
                clusterYpos = rowNum * rowHeight;
                //clusterXpos = Xcenterdie;
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
        } DESIGN_END_FOR;*/
        


        //Distort cluster position as the clusters are having the same X, Y , height and width due to which the gradient does not change at all. 
        // Distort by a factor for 3000 which is averagewidth/1000
        uint distortFactor = 0;
        uint smallDistort = 0;
        /*DESIGN_FOR_ALL_CLUSTERS((*this), cellName, clusterCellPtr) {
                clusterXpos = distortFactor ;
                clusterYpos = smallDistort ;
                //clusterXpos = Xcenterdie;
                //clusterYpos = Ycenterdie;
                (*clusterCellPtr).CellSetXpos(clusterXpos); 
                distortFactor = distortFactor + 3000;
                smallDistort = smallDistort + 300;
                //siteNum++;
                //(*clusterCellPtr).CellSetYpos(clusterYpos); 
                //if (siteNum == numSites) {
                //        siteNum = 0;
                 //       rowNum++;
                //}   
                clusterCells.push_back(clusterCellPtr);
                count++;
        } DESIGN_END_FOR;*/


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
                // The below debug info is required for debugging same gradient value in every iteration of solver

                Cell *cellPtr;
                string cellName;

                /*DESIGN_FOR_ALL_CELLS((*this),cellName,cellPtr){
                        if ((*cellPtr).CellIsTerminal()) continue;
                        cout << "Cell Name is: " << cellName ;
                        cout << " cellXpos is: " << (*cellPtr).CellGetXposDbl(); 
                        cout << " cellYpos is: " << (*cellPtr).CellGetYposDbl();
                        cout << " cellHeight is: " << (*cellPtr).CellGetHeight();
                        cout << " cellWidth is: " << (*cellPtr).CellGetWidth() << endl;
                 }DESIGN_END_FOR;*/


        }       
        /*End Debug*/
        
        // Populate the list  of the cluster cells for which the nonlinear placer has to be invoked 
        vector<Cell *> cellsToSolve;
        getNLPCellsToSolveNew ((*this),cellsToSolve);
        numVars= (cellsToSolve.size());

        double valuesX[numVars];
        double valuesY[numVars];

        /// Initialize  objective function for the X and Y values to be solved
        nlpInitialValuesNew(cellsToSolve,valuesX,valuesY);
        DesignComputepenaltyParameter();

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

        //Make the bounds constraints on the X and Y variable

        constraint_listX = NULL;
        //constraint 1 - set minimum value of x and y
        for (iX=0;iX<numVars;++iX){
                wn_make_linear_constraint(&linear_constraintX,1,0.0,WN_GT_COMPARISON);
                (linear_constraintX->vars)[0] = iX;
                (linear_constraintX->weights)[0] = 1.0;
                wn_sllins(&constraint_listX,linear_constraintX);
        }
        //constraint 2 - set maximum value of x

        for (iX=0;iX<numVars;++iX){
                wn_make_linear_constraint(&linear_constraintX,1,constMaxx,WN_LT_COMPARISON);
                (linear_constraintX->vars)[0] = iX;
                (linear_constraintX->weights)[0] = 1.0;
                wn_sllins(&constraint_listX,linear_constraintX);
        }


        constraint_listY = NULL;
        //constraint 1 - set minimum value of x and y
        for (iY=0;iY<numVars;++iY){
                wn_make_linear_constraint(&linear_constraintY,1,0.0,WN_GT_COMPARISON);
                (linear_constraintY->vars)[0] = iY;
                (linear_constraintY->weights)[0] = 1.0;
                wn_sllins(&constraint_listY,linear_constraintY);
        }
        //constraint 2 - set maximum value of x

        for (iY=0;iY<numVars;++iY){
                wn_make_linear_constraint(&linear_constraintY,1,constMaxy,WN_LT_COMPARISON);
                (linear_constraintY->vars)[0] = iY;
                (linear_constraintY->weights)[0] = 1.0;
                wn_sllins(&constraint_listY,linear_constraintY);
        }


        double *deltaX =NULL;
        double *deltaY =NULL;

        /**************** Call the conjugate gradient minimizer ******************************/

        // For X coordinates
        //Commented for compiling
        wn_nlp_conj_method(&codeX,&val_minX,valuesX,deltaX,(wn_nonlinear_constraint_type)objectiveX,constraint_listX,numVars,numVars,10,1);
        cout << " Final value after X minimization is: " << val_minX << " codeX: " <<codeX <<endl;

        wn_nlp_conj_method(&codeY,&val_minY,valuesY,deltaY,(wn_nonlinear_constraint_type)objectiveY,constraint_listY,numVars,numVars,10,1);
        cout << " Final value after Y minimization is: " << val_minY << " codeY: " <<codeY <<endl;


        cout << "Final HPWL after minimization is: " <<  (*this).DesignComputeLseHPWL()<<endl;

        for(int i=0;i<numVars;++i){

                cout << "Final X: " <<valuesX[i]<<" Final Y: "<<valuesY[i]<<endl;
         //       i=i+1;
        }

   double percentageOverlap;
   percentageOverlap = (*this).DesignDumpClusterOverlapForPenalty();
   cout << "Percentage Overlap for debug" << percentageOverlap << endl;

}


