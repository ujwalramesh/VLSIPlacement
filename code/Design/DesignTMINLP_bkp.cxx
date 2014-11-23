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
        for (int idx=0;idx<m;idx++){
                const_types[idx]=Ipopt::TNLP::LINEAR;
        }
        
        return true;}

bool Design::get_nlp_info(Index& n, Index&m, Index& nnz_jac_g,Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style){

std::vector<Cell *> cellsToSolve;
getNLPCellsToSolveNew ((*this),cellsToSolve);
int size = 2*(cellsToSolve.size());
n = size;
m=(((size/2)*((size/2)+1))/2)*4;

nnz_jac_g = m*2;
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
        idx=0;
        int numOneDConst = m/2;
   for (idx=0;idx<numOneDConst;idx++){
           g_l[idx]=averageClusterWidth;
   }
   idx=0;
   for (idx=numOneDConst;idx<m;idx++){
           g_l[idx]=averageClusterHeight;
   }
   
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
        //(*this).DesignUpdateGridPotentials();
        //totalDensityPenalty =(*this).DesignComputeTotalDensityPenalty();
        obj_value = LseHPWLconv+totalDensityPenalty;
   /*Adding penalty function for density constraint as per Will Naylor's patent*/
        
return true;
}

bool Design::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f){
        Cell *cellPtr;
        string cellName;
        uint alpha = 500;
        uint idx=0;
        double cellXpos;
        double cellYpos;
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
                grad_f[idx] = gradX;
                grad_f[idx+1]=gradY;
                cout <<"Cell Name: " <<cellName << " CellXpos: "<<cellXpos<< " cellYpos: "<<cellYpos<<" gradX: "<<gradX<<" gradY: "<<gradY<<endl;
                idx=idx+2;
        }DESIGN_END_FOR;
        
return true;
}


bool Design::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g){
        int idx=0;
        std::vector<Cell *> cellsToSolve;
        getNLPCellsToSolveNew ((*this),cellsToSolve);
        int size = 2*cellsToSolve.size();
        int idx2=0
        for (int i=0;i<size;i=i+2) {
                      for (int j=0;j<size;j=j+2){
                              if (j <= i) continue;
                              g[idx] = x[i]-x[j] + 152471 * x[idx2+size];
                              idx=idx+1;
                              g[idx] =x[i]-x[j] + 152471 * x[idx2+size];
                              idx=idx+1;
                              idx2=idx2+1;
                      }
        }
        for (int i=1;i<size;i=i+2) {
                      for (int j=1;j<size;j=j+2){
                              if (j <= i) continue;
                              g[idx] = x[i]-x[j] + 152471 * x[idx2+size];
                              idx=idx+1;
                              g[idx] =x[i]-x[j] 152471 * x[idx2+size];
                              idx=idx+1;
                              idx2=idx2+1;
                      }
         }
              

        return true;}

bool Design::eval_jac_g(Index n, const Number* x, bool new_x,
                        Index m, Index nele_jac, Index* iRow, Index *jCol,
                        Number* values){
        if (values==NULL){
                int idj =0;
                for (int idx=0;idx<m;idx++){
                        iRow[idj]=idx;
                        jCol[idj]=0;
                        idj=idj+1;
                        iRow[idj]=idx;
                        jCol[idj]=1;
                        idj=idj+1;
                }
       } else {

       int idx=0;  
        for (int i=0;i<n;i=i+2) {
                      for (int j=0;j<n;j=j+2){
                              if (j <= i) continue;
                              values[idx] = x[i];
                              idx=idx+1;
                              values[idx] = -x[j];
                              idx=idx+1;
                              values[idx] =-x[i];
                              idx=idx+1;
                              values[idx] = x[j];
                              idx=idx+1;
                      }
        }
        for (int i=1;i<n;i=i+2) {
                      for (int j=1;j<n;j=j+2){
                              if (j <= i) continue;
                              values[idx] = x[i];
                              idx=idx+1;
                              values[idx] = -x[j];
                              idx=idx+1;
                              values[idx] =-x[i];
                              idx=idx+1;
                              values[idx] = x[j];
                              idx=idx+1;
                      }
         }




       }
        return true;}

bool Design::eval_h(Index n, const Number* x, bool new_x,
                        Number obj_factor, Index m, const Number* lambda,
                        bool new_lambda, Index nele_hess, Index* iRow,
                        Index* jCol, Number* values){
        iRow=NULL;
        jCol=NULL;
        values=NULL;
        lambda=NULL;
        
        return false;}

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

