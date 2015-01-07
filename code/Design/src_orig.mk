MODULE = Design
SRCFILES:=  DesignMain DesignUtils DesignRead DesignGraph  DesignAnalysis DesignWrite \
	    DesignProperties DesignSpread DesignDebug \
	    DesignSolveFastConjGrad DesignPlace DesignCluster DesignClusterStrategy \
	    DesignClusterBestChoice DesignClusterNetCluster DesignClusterKWay DesignClusterLarge \
	    DesignWriteCluster DesignSolveForceDirected DesignDump DesignSolveWnnlpNew \
            DesignGrid DesignSolveMINLP DesignTMINLP_LPtest2 

HFILES:= Design DesignIter PriorityQueue HyperGraph Cell Pin Net Env 
