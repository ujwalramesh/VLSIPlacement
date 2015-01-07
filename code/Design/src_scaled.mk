MODULE = Design
SRCFILES:=  DesignMain DesignUtils DesignRead DesignGraph  DesignAnalysis DesignWrite \
	    DesignProperties DesignSpread DesignDebug \
	    DesignSolveFastConjGrad DesignPlace DesignCluster DesignClusterStrategy \
	    DesignClusterBestChoice DesignClusterNetCluster DesignClusterKWay DesignClusterLarge \
	    DesignWriteCluster DesignSolveForceDirected DesignDump DesignSolveWnnlpNew \
            DesignGrid DesignSolveMINLPscaled DesignTMINLP_LPtest2Scaled 

HFILES:= Design DesignIter PriorityQueue HyperGraph Cell Pin Net Env 
