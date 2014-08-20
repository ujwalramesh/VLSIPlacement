#include <Design.h>

#include <wnlib.h>
#include <wnasrt.h>
#include <wnsll.h>
#include <wnnlp.h>
void
getObjFuncInWnnlpFormat()
{

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
