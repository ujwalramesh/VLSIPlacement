#include "Design.h"

/*local double myDivideNew(double num, double denom) 
{
       if (denom <= 0.0)
        {   
           return (HUGE_VAL);
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



local void
getNLPCellsToSolveNew(Design &myDesign,vector<Cell *> &cellsToSolve)
{

        Cell *cellPtr;
        string cellName;

        DESIGN_FOR_ALL_CELLS(myDesign,cellName,cellPtr){
                        if ((*cellPtr).CellIsTerminal()) continue;
                                cellsToSolve.push_back(cellPtr);
        }DESIGN_END_FOR;

}*/



/*The below function returns the variable types involved in the problem*/
bool Design::get_variables_types(Index n, VariableType* var_types){
        std::vector<Cell *> cellsToSolve;
        getNLPCellsToSolveNew ((*this),cellsToSolve);

        int numVars = 2*(cellsToSolve.size());
          
        for (int i=0;i<numVars;i++){
                var_types[i]=CONTINUOUS;
      
        }
       
return true;
}


/* The below function returns the linearity of the variables involved in the problem*/
bool Design::get_variables_linearity(Index n, Ipopt::TNLP::LinearityType* var_types){
        
        std::vector<Cell *> cellsToSolve;
        getNLPCellsToSolveNew ((*this),cellsToSolve);

        int numVars = 2*(cellsToSolve.size());
        for (int i=0;i<numVars;i++){
                var_types[i]=Ipopt::TNLP::NON_LINEAR;
        }
        
return true;
}

bool Design::get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_types){
//        const_types[0]=Ipopt::TNLP::LINEAR;
      //  const_types[1]=Ipopt::TNLP::LINEAR;

        return true;}

bool Design::get_nlp_info(Index& n, Index&m, Index& nnz_jac_g,Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style){

std::vector<Cell *> cellsToSolve;
getNLPCellsToSolveNew ((*this),cellsToSolve);
int size = 2*(cellsToSolve.size());
n = size;
m=0;

nnz_jac_g = 0;
nnz_h_lag = 0;

index_style = TNLP::C_STYLE;

        
return true;
}



/* The below function will be used to get the bounds for variables and constraints  - Initially no constriants are coded*/  
bool Design::get_bounds_info(Index n, Number* x_l, Number* x_u,Index m, Number* g_l, Number* g_u){
        
        uint maxx,maxy;

        DesignGetBoundingBox(maxx,maxy);

        uint averageClusterWidth,averageClusterHeight;

        averageClusterWidth = (uint)DesignGetAverageClusterCellWidth();
        averageClusterHeight = (uint)DesignGetAverageClusterCellHeight();

        uint XupperBound = maxx-averageClusterWidth;
        uint YupperBound = maxy-averageClusterHeight;

        uint idx=0;
        // Lower bound for x and y variable is 0 but the upperbound is maxx/maxy - averageClusterWidth/averageClusterHeight 
        for (idx=0;idx<n;idx++){
                x_l[idx]=0; 
                x_l[idx+1]=0;
                x_u[idx]=XupperBound;
                x_u[idx+1]=YupperBound;
                idx=idx+1;
        }
   /*Initially the number of constraints is 0. Will have to add it once constraints are modeled*/
             //   g_l[0]=10;
                //g_u[1]=10000;
           return true;
}



/* The below function is used to obtain the intial solution for the NLP problem*/
bool Design::get_starting_point(Index n, bool init_x, Number* x,bool init_z, Number* z_L, Number* z_U,
                                Index m, bool init_lambda,Number* lambda){
        assert(init_x);
        assert(!init_lambda);
        assert(!init_z); 
        Cell *cellPtr;
        string cellName;
        uint idx=0;
        
        std::vector<Cell *> cellsToSolve;
        getNLPCellsToSolveNew ((*this),cellsToSolve);
        
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


return true;
}
   
/* The below function computes the value of the objective*/

bool Design::eval_f(Index n, const Number* x, bool new_x, Number& obj_value){
        
  /* objective here is total_wirelength = sum over all nets ( wirelength (net) * net weight) */ 
  Cell *cellPtr;
  string cellName;
  uint idx=0;
        //cout << " Void pointer Type: " <<typeid(*((Design*)myDesign)).name() << endl;  
        DESIGN_FOR_ALL_CELLS((*this),cellName,cellPtr){
//                cout << "Changing Xposition for cell: " << (*cellPtr).CellGetName() << " to " << values[idx] << endl;
                   if ((*cellPtr).CellIsTerminal()) continue;     
                (*cellPtr).CellSetXposDbl(x[idx]);
                (*cellPtr).CellSetYposDbl(x[idx+1]);
  //              cout << "The new Xposition of the cell is " << (*cellPtr).CellGetXpos() << endl;
                idx = idx+2;

        }END_FOR;
        ulong LseHPWL;
        double LseHPWLconv;
        LseHPWL = (*this).DesignComputeLseHPWL();
        LseHPWLconv=LseHPWL; 
        double totalDensityPenalty=0;
        //Penalty Method 1  : Used Grid / Bin penalty method
       (*this).DesignUpdateGridPotentials();
       totalDensityPenalty =(*this).DesignComputeTotalDensityPenalty();


        //Penalty Method 2 : using penalty based on percentage overlap and tried out of intution by taking (percentageOverlap)^3 to make it differentiable
       /* double totalOverlap;
        totalOverlap = (*this).DesignDumpClusterOverlapForPenalty();
        uint lambda = 22; 
        bool loop1=true;
        bool loop2=true;
        bool loop3=true;


        totalDensityPenalty= lambda*pow(totalOverlap,3);
        //lambda=lambda;


        if (totalOverlap > 300 ) {
                totalDensityPenalty= 64*totalDensityPenalty;
        } else if ((totalOverlap >250) && (totalOverlap <= 300 )) {
                totalDensityPenalty= 32*totalDensityPenalty;
        } else if ((totalOverlap >200 ) && (totalOverlap <= 250 )){
                totalDensityPenalty= 16*totalDensityPenalty;
        } else if ((totalOverlap > 150 ) && (totalOverlap <= 200 )){
                totalDensityPenalty= 8*totalDensityPenalty;
        } else if ((totalOverlap > 100 ) && (totalOverlap <= 150 )){
                totalDensityPenalty= 4*totalDensityPenalty;
        } else if ((totalOverlap > 50 ) && (totalOverlap <= 100 )){
                totalDensityPenalty= 2*totalDensityPenalty;
        } else if ((totalOverlap > 10 ) && (totalOverlap <= 50 )){
                totalDensityPenalty= totalDensityPenalty;
        } else if ((totalOverlap <= 10 ) && (totalOverlap >= 0 )) {
                totalDensityPenalty = 0 ;
        }*/
        
        
       double penaltyParameter;
         penaltyParameter = (*this).DesignGetpenaltyParameter();

            if (penaltyParameter < 0 ) {
                  penaltyParameter = -penaltyParameter;
            }
        
        //obj_value = LseHPWLconv+(penaltyParameter*totalDensityPenalty);
        obj_value = LseHPWLconv;
   /*Adding penalty function for density constraint as per Will Naylor's patent*/
      cout << "ulong wirelength is: " << LseHPWL << "\t";
      cout << "Density Penalty is : " << penaltyParameter*totalDensityPenalty << "\t";
      cout << "Density Penalty is : " << obj_value<<endl;
      
return true;
}

bool Design::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f){
        Cell *cellPtr;
        string cellName;
        uint alpha = 500;
        uint idx=0;
        double cellXpos;
        double cellYpos;
        double densityPenaltyGradient;
        double penaltyParameter;
        penaltyParameter = (*this).DesignGetpenaltyParameter();
        densityPenaltyGradient = (*this).DesignComputeTotalDensityPenaltyGradient();
        cout << "Total Density Penalty Gradient: " << densityPenaltyGradient << endl;
        DESIGN_FOR_ALL_CELLS((*this),cellName,cellPtr){
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
                double gradPotentialX;
                double gradPotentialY;
                (*this).DesignComputePenaltyGradientforCell(cellPtr,gradPotentialX,gradPotentialY);
                //grad_f[idx] = gradX + penaltyParameter*gradPotentialX*densityPenaltyGradient;
                //grad_f[idx+1]=gradY + penaltyParameter*gradPotentialY*densityPenaltyGradient;
                grad_f[idx] = gradX;
                grad_f[idx+1]=gradY;
                cout <<"Cell Name: " <<cellName << " CellXpos: "<<cellXpos<< " cellYpos: "<<cellYpos<<" gradX: "<<grad_f[idx]<<" gradY: "<<grad_f[idx+1]<<endl;
                idx=idx+2;
        }DESIGN_END_FOR;
        
return true;
}


bool Design::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g){
     //   g[0]=x[0]-x[2];
      //  g[1]=-x[0]+x[2];
        //g=NULL; 
        return true;}

bool Design::eval_jac_g(Index n, const Number* x, bool new_x,
                        Index m, Index nele_jac, Index* iRow, Index *jCol,
                        Number* values){
      /*  if (values == NULL) {
                iRow[0]=0;
                jCol[0]=0;
                iRow[1]=0;
                iRow[1]=1;
        //        iRow[2]=1;
         //       jCol[2]=0;
          //      iRow[3]=1;
          //      iRow[3]=1;
        } else {
                values[0]=-x[2];
                values[1]=x[0];
           //     values[2]=x[2];
           //     values[3]=-x[0];
        }*/
        //iRow=NULL;
        //jCol=NULL;
        //values=NULL;


       
        return true;}

bool Design::eval_h(Index n, const Number* x, bool new_x,
                        Number obj_factor, Index m, const Number* lambda,
                        bool new_lambda, Index nele_hess, Index* iRow,
                        Index* jCol, Number* values){
        //iRow=NULL;
        //jCol=NULL;
        //values=NULL;
        //lambda=NULL;
        
        return true;}



Index Design::get_number_of_nonlinear_variables(void){

        std::vector<Cell *> cellsToSolve;
        getNLPCellsToSolveNew ((*this),cellsToSolve);
        Index size = 2*cellsToSolve.size();
        return size;
}


bool Design::get_list_of_nonlinear_variables(Index num_nonlin_vars,Index* pos_nonlin_vars) {


     
        std::vector<Cell *> cellsToSolve;
        getNLPCellsToSolveNew ((*this),cellsToSolve);
        Index size = 2*cellsToSolve.size();
        assert (num_nonlin_vars == size);
     
        for (int i=0;i<size;i++){
                  pos_nonlin_vars[i]=i;
        }   


        return true;
}
       

void Design::finalize_solution(TMINLP::SolverReturn status,
                                Index n, const Number* x, Number obj_value){
        
 std::cout<<"Problem status: "<<status<<std::endl;
   std::cout<<"Objective value: "<<obj_value<<std::endl;
     if(x != NULL){
        std::cout<<"Solution:"<<std::endl;
        for(int i = 0 ; i < n ; i++){
        std::cout<<"x["<<i<<"] = "<<x[i];
        i=i+1;
        std::cout<<"y["<<i<<"] = "<<x[i]<<endl;
              
        }
        std::cout<<std::endl;
        }

        
}

