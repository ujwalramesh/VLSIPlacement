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
double maxx,maxy;
double constMaxx,constMaxy;
ulong lseXHPWL,lseYHPWL;
double averageClusterWidth,averageClusterHeight;
double numClusters,numRows,numSites;
double siteWidth,rowHeight;

/*Create placeable blocks in the design*/
DesignGetBoundingBox(maxx,maxy);
averageClusterWidth = (double)DesignGetAverageClusterCellWidth();
averageClusterHeight = (double)DesignGetAverageClusterCellHeight();

maxx=maxx/scaleFactor;
maxy=maxy/scaleFactor;
averageClusterWidth = averageClusterWidth/scaleFactor;
averageClusterHeight = averageClusterHeight/scaleFactor;


numClusters = DesignGetNumClusters();
numRows = floor(((double)maxy) / averageClusterHeight);
numSites = ceil(((double)numClusters) / numRows);
siteWidth = ((double)maxx) / numSites;
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
double Xcenterdie,Ycenterdie;
Xcenterdie=maxx/2;
Ycenterdie=maxy/2;


/*DESIGN_FOR_ALL_CLUSTERS((*this), cellName, clusterCellPtr) {
           //     clusterXpos = siteNum * siteWidth;
            //    clusterYpos = rowNum * rowHeight;
                clusterXpos = Xcenterdie;
                clusterYpos = Ycenterdie;
                (*clusterCellPtr).CellSetXpos(clusterXpos);
           //     siteNum++;
                (*clusterCellPtr).CellSetYpos(clusterYpos);
           //      if (siteNum == numSites) {
            //                siteNum = 0;
             //               rowNum++;
              //  }
                clusterCells.push_back(clusterCellPtr);
                count++;
} DESIGN_END_FOR;
*/


siteNum=0;
rowNum =0;
count =0;

// Trying to distort cell positions intially so as to avoid redundant constraints
double siteWidthDistort, rowHeightDistort; 
siteWidthDistort = 3000/scaleFactor;
rowHeightDistort = 3000/scaleFactor;
DESIGN_FOR_ALL_CLUSTERS((*this), cellName, clusterCellPtr) {
              //  clusterXpos = Xcenterdie+siteNum * siteWidthDistort;
              // clusterYpos = Ycenterdie+rowNum * rowHeightDistort;
                clusterXpos = Xcenterdie;
                clusterYpos = Ycenterdie;
                //clusterXpos = siteNum * siteWidthDistort;
                //clusterYpos = rowNum * rowHeightDistort;
                (*clusterCellPtr).CellSetXposDbl(clusterXpos*scaleFactor);
                siteNum++;
                (*clusterCellPtr).CellSetYposDbl(clusterYpos*scaleFactor);
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
        FILE * fp = fopen("bonminLog.out","w");
        CoinMessageHandler handler(fp);
        BonminSetup bonmin(&handler);
#else
        BonminSetup bonmin;
#endif
bonmin.initializeOptionsAndJournalist();
 bonmin.roptions()->AddStringOption2("print_solution","Do we print the solution or not?",
                                  "yes",
                                  "no", "No, we don't.",
                                 "yes", "Yes, we do.",
                                  "A longer comment can be put here");
    // Here we can change the default value of some Bonmin or Ipopt option
     bonmin.options()->SetNumericValue("bonmin.time_limit", 5); //changes bonmin's time limit
     bonmin.options()->SetStringValue("mu_oracle","loqo");
     bonmin.options()->SetStringValue("print_info_string","yes");
     bonmin.options()->SetStringValue("print_user_options","yes");
     //bonmin.options()->SetStringValue("check_derivatives_for_naninf","yes");
     bonmin.options()->SetIntegerValue("bonmin.bb_log_level",3);
     bonmin.options()->SetIntegerValue("print_level",6);
     //bonmin.options()->SetIntegerValue("obj_scaling_factor",100000);
     //bonmin.options()->SetIntegerValue("bound_relax_factor",10);
    // bonmin.options()->SetIntegerValue("max_iter",2147483647);
     //bonmin.options()->SetIntegerValue("max_iter",6000);
     //bonmin.options()->SetIntegerValue("max_iter",110);
     bonmin.options()->SetStringValue("hessian_approximation","limited-memory");
     //bonmin.options()->SetStringValue("derivative_test","first-order");
     bonmin.options()->SetStringValue("output_file","Ipopt.log");
     //Here we read several option files
     bonmin.readOptionsFile("Mybonmin.opt");
    // bonmin.readOptionsFile();// This reads the default file "bonmin.opt"

     // Options can also be set by using a string with a format similar to the bonmin.opt file
     bonmin.readOptionsString("bonmin.algorithm B-BB\n");
     // Now we can obtain the value of the new option
     int printSolution;
     bonmin.options()->GetEnumValue("print_solution", printSolution,"");
      if(printSolution == 1){
     //      tminlp->printSolutionAtEndOfAlgorithm();
      }
    //Now initialize from tminlp
    bonmin.initialize(GetRawPtr(tminlp),false);
    //bonmin.initialize(this,false);
    //Set up done, now let's branch and bound
    try {
       Bab bb;
       bb(bonmin);//process parameter file using Ipopt and do branch and bound using Cbc
     }
   catch(TNLPSolver::UnsolvedError *E) {
//There has been a failure to solve a problem with Ipopt.
    std::cerr<<"Ipopt has failed to solve a problem"<<std::endl;
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
}


