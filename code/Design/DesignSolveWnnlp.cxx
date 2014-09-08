#include <Design.h>

#include <wnlib.h>
#include <wnasrt.h>
#include <wnsll.h>
#include <wnnlp.h>
void
getObjFuncInWnnlpFormat()
{

}

double logSumExpWireLength ( )
{

}
double wirelengthObjFunc(int size, double *values,ptr client_data) /* Need to justify size, values and client_data. Just reused it from definition of pfunction*/
{
  /* objective here is total_wirelength = sum over all nets ( wirelength (net) * net weight) */ 
  /* First step is to get list of all nets, second is to obtain net weights and then form wirelength function using log sum exponential wirelength model*/
        double totalWirelength;
        double num_nets; 

        for(int i=0;i<num_nets;i++)
        {
          totalWirelength =+ logSumExpWireLength ( );       
        
        } 
                
        return totalWirelength; 

}
        

void
Design::DesignSolveForAllCellsWnnlp(void)
{

//All Variable Declarations


string DesignPath, DesignName;
string DirName;
// All initialization
Env &DesignEnv = (*this).DesignEnv;
DesignPath = DesignEnv.EnvGetDesignPath();
DirName = DesignPath + "/.solverData";
DesignName = DesignEnv.EnvGetDesignName();

HyperGraph &myGraph = (*this).DesignGetGraph();

}
