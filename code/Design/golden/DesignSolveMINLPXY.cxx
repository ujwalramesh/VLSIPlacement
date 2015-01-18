#include <Design.h>

#include <math.h> 

#include <iomanip>
#include <fstream>

#include "CoinPragma.hpp"
#include "CoinTime.hpp"
#include "CoinError.hpp"

#include "BonOsiTMINLPInterface.hpp"
#include "BonIpoptSolver.hpp"
#include "BonCbc.hpp"
#include "BonBonminSetup.hpp"

#include "BonOACutGenerator2.hpp"
#include "BonEcpCuts.hpp"
#include "BonOaNlpOptim.hpp"
#include "BonOaDecBase.hpp"
#define REDIRECT

/*local void
getNLPCellsToSolveNew(Design &myDesign,vector<Cell *> &cellsToSolve)
{

        Cell *cellPtr;
        string cellName;
        DESIGN_FOR_ALL_CELLS(myDesign,cellName,cellPtr){
             if ((*cellPtr).CellIsTerminal()) continue;
                cellsToSolve.push_back(cellPtr);
       }DESIGN_END_FOR;
}*/

void
Design::DesignSolveForAllCellsMINLP(void)
{
using namespace boost;
using namespace Ipopt;
using namespace Bonmin;

/* WIll start defining problem by passing values to TMINLP virtual functions*/ 
bool debug = true;

/* Variable declarations for cluster cell palcement*/
double scaleFactor = (*this).DesignComputeScalingFactor();
cout << " The computed scaling factor is " << scaleFactor << endl;

Cell *clusterCellPtr;
std::vector<Cell *> clusterCells;
double clusterXpos,clusterYpos;
uint maxx,maxy;
uint constMaxx,constMaxy;
ulong lseXHPWL,lseYHPWL;
uint averageClusterWidth,averageClusterHeight;
uint numClusters,numRows,numSites;
uint siteWidth,rowHeight;

/*Create placeable blocks in the design*/
DesignGetBoundingBox(maxx,maxy);
averageClusterWidth = (uint)DesignGetAverageClusterCellWidth();
averageClusterHeight = (uint)DesignGetAverageClusterCellHeight();

numClusters = DesignGetNumClusters();
numRows = floor(((double)maxy) / averageClusterHeight);
numSites = ceil(((double)numClusters) / numRows);
siteWidth = floor(((double)maxx) / numSites);
rowHeight = averageClusterHeight;
constMaxx = maxx -  averageClusterWidth;
constMaxy = maxy - averageClusterHeight;

uint siteNum,rowNum;
/* Place clusters Constructively*/
siteNum = 0;
rowNum = 0;
uint count = 0;
string cellName;

// Constructive placement commented to check the working of Bonmin
uint Xcenterdie,Ycenterdie;
Xcenterdie=maxx/2;
Ycenterdie=maxy/2;

/*
DESIGN_FOR_ALL_CLUSTERS((*this), cellName, clusterCellPtr) {
                clusterXpos = siteNum * siteWidth;
                clusterYpos = rowNum * rowHeight;
                clusterXpos = Xcenterdie;
                clusterYpos = Ycenterdie;
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


siteNum=0;
rowNum =0;
count =0;

// Trying to distort cell positions intially so as to avoid redundant constraints
uint siteWidthDistort, rowHeightDistort; 
siteWidthDistort = 3000;
rowHeightDistort = 3000;
DESIGN_FOR_ALL_CLUSTERS((*this), cellName, clusterCellPtr) {
                clusterXpos = Xcenterdie+siteNum * siteWidthDistort;
              clusterYpos = Ycenterdie+rowNum * rowHeightDistort;
            //    clusterXpos = Xcenterdie;
           //     clusterYpos = Ycenterdie;
              //  clusterXpos = siteNum * siteWidthDistort;
              //  clusterYpos = rowNum * rowHeightDistort;
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

std::vector<Cell *> cellsToSolve; // Now coded in get_variable_types 
getNLPCellsToSolveNew ((*this),cellsToSolve);
Index numVars= 2*(cellsToSolve.size());
Number x[numVars];

DesignComputepenaltyParameter();

if (debug) {
        cout << "maxx " << maxx << " maxy " << maxy << endl;
        get_starting_point(numVars,true,x,false,NULL,NULL,2*numVars,false,NULL);
        //Number wirelength;
        //eval_f(numVars,x,false,wirelength);
        //cout << "wirelength minlp is " << wirelength<<endl;
}

/* Get All the variables to be solved for from the design class*/


/*Populate Initial Values for the variables*/



/* Call the Bonmin functions ---- Prototype taken from /home/rameshul/Bonmin-1.7/Bonmin/examples/CppExample/MyBonmin.cpp ------ */
//shared_ptr<Design> p = (*this).returnSharedPointerFromThis();
SmartPtr<Design> tminlp = this;
//tminlp = this;
#ifdef REDIRECT
        FILE * fpX = fopen("bonminXLog.out","w");
        CoinMessageHandler handlerX(fpX);
        BonminSetup bonminX(&handlerX);
#else
        BonminSetup bonminX;
#endif
bonminX.initializeOptionsAndJournalist();
 bonminX.roptions()->AddStringOption2("print_solution","Do we print the solution or not?",
                                  "yes",
                                  "no", "No, we don't.",
                                 "yes", "Yes, we do.",
                                  "A longer comment can be put here");
    // Here we can change the default value of some Bonmin or Ipopt option
     bonminX.options()->SetNumericValue("bonminX.time_limit", 5); //changes bonminX's time limit
     bonminX.options()->SetStringValue("mu_oracle","loqo");
     //bonminX.options()->SetStringValue("mu_strategy","monotone");
     bonminX.options()->SetStringValue("print_info_string","yes");
     bonminX.options()->SetStringValue("print_user_options","yes");
     // the below option is used because the iteration had a y tag and based on http://www.gams.com/dd/docs/solvers/ipopt.pdf documents page 12 - this option is set   
          //bonminX.options()->SetStringValue("least_square_init_duals","yes");
          //bonminX.options()->SetStringValue("least_square_init_primal","yes");
          //bonminX.options()->SetStringValue("accept_every_trial_step","yes");
     //bonminX.options()->SetStringValue("check_derivatives_for_naninf","yes");

//########################################################################################################################################### 
    //     Big  M variables are not being treated as integers. Changing the below option to check if they are being treated as integers     
    // bonminX.options()->SetIntegerValue("bonminX.integer_tolerance",1);
    // Did not work. Refer to mail from Stefan why you should not use it
//############################################################################################################################################     
     
     bonminX.options()->SetIntegerValue("bonminX.bb_log_level",3);
     bonminX.options()->SetIntegerValue("print_level",7);
     //bonminX.options()->SetIntegerValue("obj_scaling_factor",100000);
    //bonminX.options()->SetIntegerValue("bound_relax_factor",0);
     //bonminX.options()->SetIntegerValue("max_iter",2147483647);
     //bonminX.options()->SetIntegerValue("max_iter",6000);
     //bonminX.options()->SetIntegerValue("max_iter",110);

//#################This is just a try. Delete If it does not work####//
    // bonminX.options()->SetNumericValue("acceptable_constr_viol_tol",1);
///////////////////////////////////////////////////
     bonminX.options()->SetStringValue("hessian_approximation","limited-memory");
    // bonminX.options()->SetStringValue("derivative_test","first-order");
     bonminX.options()->SetStringValue("output_file","IpoptX.log");
     //Here we read several option files
     bonminX.readOptionsFile("MybonminX.opt");
    // bonminX.readOptionsFile();// This reads the default file "bonminX.opt"

     // Options can also be set by using a string with a format similar to the bonminX.opt file
     //bonminX.readOptionsString("bonminX.algorithm B-OA\n");
     bonminX.readOptionsString("bonminX.algorithm B-BB\n");
     // Now we can obtain the value of the new option
     int printSolutionX;
     bonminX.options()->GetEnumValue("print_solution", printSolutionX,"");
      if(printSolutionX == 1){
     //      tminlp->printSolutionAtEndOfAlgorithm();
      }
    //Now initialize from tminlp
    bonminX.initialize(GetRawPtr(tminlp),false);
    //bonminX.initialize(this,false);
    //Set up done, now let's branch and bound
    try {
       Bab bbX;
      // OaDecompositionBase bb;
       bbX(bonminX);//process parameter file using Ipopt and do branch and bound using Cbc
     }
   catch(TNLPSolver::UnsolvedError *E) {
//There has been a failure to solve a problem with Ipopt.
    std::cerr<<"Ipopt has failed to solve a problem in X"<<std::endl;
      }
    catch(OsiTMINLPInterface::SimpleError &E) {
        std::cerr<<E.className()<<"::"<<E.methodName()
                        <<std::endl
                        <<E.message()<<std::endl;
    }
    catch(CoinError &E) {
    std::cerr<<E.className()<<"::"<<E.methodName()
        <<std::endl
        <<E.message()<<std::endl;
    }
//tminlp=NULL;
   DesignSolveMINLPinY = true;     
#ifdef REDIRECT
        FILE * fpY = fopen("bonminYLog.out","w");
        CoinMessageHandler handlerY(fpY);
        BonminSetup bonminY(&handlerY);
#else
        BonminSetup bonminY;
#endif
bonminY.initializeOptionsAndJournalist();
 bonminY.roptions()->AddStringOption2("print_solution","Do we print the solution or not?",
                                  "yes",
                                  "no", "No, we don't.",
                                 "yes", "Yes, we do.",
                                  "A longer comment can be put here");
    // Here we can change the default value of some Bonmin or Ipopt option
     bonminY.options()->SetNumericValue("bonminY.time_limit", 5); //changes bonminY's time limit
     bonminY.options()->SetStringValue("mu_oracle","loqo");
     //bonminY.options()->SetStringValue("mu_strategy","monotone");
     bonminY.options()->SetStringValue("print_info_string","yes");
     bonminY.options()->SetStringValue("print_user_options","yes");
     // the below option is used because the iteration had a y tag and based on http://www.gams.com/dd/docs/solvers/ipopt.pdf documents page 12 - this option is set   
          //bonminY.options()->SetStringValue("least_square_init_duals","yes");
          //bonminY.options()->SetStringValue("least_square_init_primal","yes");
          //bonminY.options()->SetStringValue("accept_every_trial_step","yes");
     //bonminY.options()->SetStringValue("check_derivatives_for_naninf","yes");

//########################################################################################################################################### 
    //     Big  M variables are not being treated as integers. Changing the below option to check if they are being treated as integers     
    // bonminY.options()->SetIntegerValue("bonminY.integer_tolerance",1);
    //     Did not work. Refer to mail from Stefan why you should not use it 
//############################################################################################################################################     
     
     bonminY.options()->SetIntegerValue("bonminY.bb_log_level",3);
     bonminY.options()->SetIntegerValue("print_level",7);
     //bonminY.options()->SetIntegerValue("obj_scaling_factor",100000);
    //bonminY.options()->SetIntegerValue("bound_relax_factor",0);
     //bonminY.options()->SetIntegerValue("max_iter",2147483647);
     //bonminY.options()->SetIntegerValue("max_iter",6000);
     //bonminY.options()->SetIntegerValue("max_iter",110);

//#################This is just a try. Delete If it does not work####//
    // bonminY.options()->SetNumericValue("acceptable_constr_viol_tol",1);
///////////////////////////////////////////////////
     bonminY.options()->SetStringValue("hessian_approximation","limited-memory");
    // bonminY.options()->SetStringValue("derivative_test","first-order");
     bonminY.options()->SetStringValue("output_file","IpoptY.log");
     //Here we read several option files
     bonminY.readOptionsFile("MybonminY.opt");
    // bonminY.readOptionsFile();// This reads the default file "bonminY.opt"

     // Options can also be set by using a string with a format similar to the bonminY.opt file
     //bonminY.readOptionsString("bonminY.algorithm B-OA\n");
     bonminY.readOptionsString("bonminY.algorithm B-BB\n");
     // Now we can obtain the value of the new option
     int printSolutionY;
     bonminY.options()->GetEnumValue("print_solution", printSolutionY,"");
      if(printSolutionY == 1){
     //      tminlp->printSolutionAtEndOfAlgorithm();
      }
    //Now initialize from tminlp
    bonminY.initialize(GetRawPtr(tminlp),false);
    //bonminY.initialize(this,false);
    //Set up done, now let's branch and bound
   try {
           Bab bbY;
           bbY(bonminY);
   }
   catch(TNLPSolver::UnsolvedError *E) {
           std::cerr<<"Ipopt has failed to solve a problem in Y "<<std::endl;
   }
   catch(OsiTMINLPInterface::SimpleError &E) {
        std::cerr<<E.className()<<"::"<<E.methodName()
                        <<std::endl
                        <<E.message()<<std::endl;
    }
    catch(CoinError &E) {
    std::cerr<<E.className()<<"::"<<E.methodName()
        <<std::endl
        <<E.message()<<std::endl;
    }


}
