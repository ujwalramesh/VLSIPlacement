# include <common.h>

# ifndef ENV_H
# define ENV_H

/***************************************************
   This class is meant to store specifics
   about the environment of the tool and the 
   configurations that the tool is currently being 
   run in. An instantiation of the class can 
   be made in the class Design.
***************************************************/

/***************************************************
  SOME TYPE DEFINITIONS FOR THE ENVIRONMENT
***************************************************/
/* Define the solver type to be used in the design */
typedef enum {
  ENV_SOLVER_QUADRATIC_MOSEK=0,
  ENV_SOLVER_QUADRATICXY_MOSEK,
  ENV_SOLVER_QUADRATIC_CONJ_GRAD,
  ENV_SOLVER_QUADRATICXY_CONJ_GRAD,
  ENV_SOLVER_NON_LINEAR,
  ENV_SOLVER_NON_LINEAR_CONJ_GRAD,
  NUM_ENV_SOLVERS
} EnvSolverType;
# define DEFAULT_ENV_SOLVER_TYPE ENV_SOLVER_QUADRATIC_MOSEK

/* Define the net model used for hypernets of the placement tool  */
typedef enum {
  ENV_CLIQUE_MODEL = 0,
  ENV_STAR_MODEL,
  ENV_HYBRID_MODEL,
  ENV_VARIABLE_HYBRID_MODEL,
  NUM_NET_MODELS
} EnvNetModel;
# define DEFAULT_ENV_NET_MODEL ENV_CLIQUE_MODEL

/* Define the operating mode of the placement tool  */
typedef enum {
  ENV_MODE_OPTIMAL = 0,
  ENV_MODE_TEST,
  NUM_MODES
} EnvToolMode;
# define DEFAULT_ENV_TOOL_MODE ENV_MODE_OPTIMAL

/* Define the legalizer used in placement tool  */
typedef enum {
  ENV_SIMPLE_LEGALIZER = 0,
  ENV_BIN_BASED_LEGALIZER,
  NUM_LEGALIZERS
} EnvLegalizer;
# define DEFAULT_ENV_LEGALIZER ENV_BIN_BASED_LEGALIZER

/* Define the type of optimization desired in the placement tool  */
typedef enum {
  ENV_OPTIMIZE_TIMING = 0,
  ENV_OPTIMIZE_WIRELENGTH
} EnvOptType;
# define DEFAULT_ENV_OPT_TYPE ENV_OPTIMIZE_TIMING

/* Define the type of clustering desired in the placement tool  */
typedef enum {
  ENV_NO_CLUSTERING = 0,
  ENV_FIRST_CHOICE_CLUSTERING,
  ENV_BEST_CHOICE_CLUSTERING,
  ENV_NET_CLUSTERING,
  ENV_EDGE_SEPARABILITY_BASED_CLUSTERING,
  ENV_TIMING_DRIVEN_CLUSTERING,
  ENV_GLOBAL_PLACEMENT_TIMING_DRIVEN_CLUSTERING,
  ENV_NUM_CLUSTERING_TYPES,
} EnvClusterType;
# define DEFAULT_ENV_CLUSTER_TYPE ENV_NO_CLUSTERING

/* Define the type of shape selection desired in the placement tool  */
typedef enum {
  ENV_NO_SHAPE_SELECTION = 0,
  ENV_STD_SHAPE_SELECTION_ALGORITHM,
  ENV_NUM_SHAPE_SELECTION_TYPES,
} EnvShapeSelectionType;
# define DEFAULT_ENV_SHAPE_SELECTION_TYPE ENV_NO_SHAPE_SELECTION

/***************************************************
 EVERY VARIABLE DEFINED BELOW HAS TO BE DOCUMENTED!
 SET AND GET FUNCTIONS ARE ADDED TO THE ENV CLASS
 TO QUERY AND SET THE VARIABLES DURING DIFFERENT
 PHASES OF OPTIMIZATION
***************************************************/
class Env {
 private:
  /***************************************/
  /* STRING VARIABLES                    */
  /***************************************/
  /* Variable to define the design name */
  string DesignName;
  /* Variable to define the design path */
  string DesignPath;
  
  /***************************************/
  /* INTEGER VARIABLES                   */
  /***************************************/
  /* Variable to define the different types of cells */
  /* Total number of objects read from the nodes file */
  uint NumCellsFromFile;
  /* Total number of ports that have been read */
  uint NumPorts;
  /* Total number of terminals that have been read */
  uint NumTerminals;
  /* Total number of fixed cells that have been read */
  uint NumFixedCells;
  /* Total number of movable cells that have to be solved 
     for */
  uint NumMovableCells;
  /* Total number of standard cells */
  uint NumStdCells;
  /* Total number of movable macros */
  uint NumMacroCells;
  /* Total number of fixed macros */
  uint NumFixedMacros;
  /* Total number of nets that are read */
  uint NumNets;

  /***************************************/
  /* BOOLEAN VARIABLES                   */
  /***************************************/
  /* Flag to determine the use of variable bounds in the solver */
  bool UseVarBounds;
  /* Flag to determine if the tool has to dump the hessian and 
     the linear matrix for testing with solvers */
  bool DumpHessian;

  /***************************************/
  /* FLOAT/DOUBLE VARIABLES              */
  /***************************************/

  /***************************************/
  /* CUSTOM ENUM VARIABLES               */
  /* ALL VALUES ARE DEFAULTED TO DEFAULT */
  /* VALUES DEFINED ABOVE                */
  /***************************************/
  /* Decide what the type of the solver to use */
  EnvSolverType SolverType;
  /* Decide what mode the tool has to be used in */
  EnvNetModel NetModel;
  /* Decide what mode the tool has to be used in */
  EnvToolMode ToolMode;
  /* Decide what type of legalizer has to be used */
  EnvLegalizer Legalizer;
  /* Decide what type of optimization cost the tool 
     has to consider */
  EnvOptType OptType;
  /* Decide what the clustering algorithm type should be */
  EnvClusterType ClusterType;
  /* Decide what the shape selection algorithm type should be */
  EnvShapeSelectionType ShapeSelectionType;

 public:
  /****************************************/
  /* SECTION SHOULD ONLY CONTAIN          */
  /* CONSTRUCTORS DESTRUCTORS AND GET AND */
  /* SET FUNCTIONS                        */
  /****************************************/
  
  /* Constructors */
  Env();
  /* Destructors */
  ~Env();
  
  /* Set and get functions in pairs since its */
  /* easy to specify it that way              */
  void EnvSetDesignName(string);
  string EnvGetDesignName(void);

  void EnvSetDesignPath(string);
  string EnvGetDesignPath(void);
  
  void EnvSetNumCellsFromFile(uint);
  uint EnvGetNumCellsFromFile(void);
  
  void EnvSetNumPorts(uint);
  uint EnvGetNumPorts(void);

  void EnvSetNumTerminals(uint);
  uint EnvGetNumTerminals(void);
  
  void EnvSetNumFixedCells(uint);
  uint EnvGetNumFixedCells(void);

  void EnvSetNumMovableCells(uint);
  uint EnvGetNumMovableCells(void);

  void EnvSetNumStdCells(uint);
  uint EnvGetNumStdCells(void);

  void EnvSetNumMacroCells(uint);
  uint EnvGetNumMacroCells(void);

  void EnvSetNumFixedMacros(uint);
  uint EnvGetNumFixedMacros(void);

  void EnvSetNumNets(uint);
  uint EnvGetNumNets(void);

  void EnvSetUseVarBounds(bool);
  bool EnvGetUseVarBounds(void);

  void EnvSetDumpHessian(bool);
  bool EnvGetDumpHessian(void);

  void EnvSetSolverType(EnvSolverType);
  EnvSolverType EnvGetSolverType(void);

  void EnvSetNetModel(EnvNetModel);
  EnvNetModel EnvGetNetModel(void);
  
  void EnvSetToolMode(EnvToolMode);
  EnvToolMode EnvGetToolMode(void);
  
  void EnvSetLegalizer(EnvLegalizer);
  EnvLegalizer EnvGetLegalizer(void);

  void EnvSetOptType(EnvOptType);
  EnvOptType EnvGetOptType(void);
  
  void EnvSetClusterType(EnvClusterType);
  EnvClusterType EnvGetClusterType(void);

  void EnvSetShapeSelectionType(EnvShapeSelectionType);
  EnvShapeSelectionType EnvGetShapeSelectionType(void);
};

# define DEFAULT_DESIGN_NAME ""
# define DEFAULT_DESIGN_PATH ""

# endif