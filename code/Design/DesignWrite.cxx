# include <Design.h>

void
DesignWriteHeaderFile(Design &myDesign, ofstream &opFile)
{
  Env &DesignEnv = myDesign.DesignGetEnv();
  opFile << "# Created  :  " << getCurrentTime() << endl;		
  opFile << "# User     :  " << getUserName() << "@" << getHostName() << endl;
  opFile << "# Platform :  " << getPlatformString() << endl;
}

void DesignWriteNets(Design &myDesign, string fname) 
{
  string netName, designName, fileName;
  string pinDir;
  string cellName, netDegree;
  Cell *cellPtr;
  Pin *pinPtr;
  Net *netPtr;
  ofstream opFile;
  uint numPins, numNetPins;
  uint numCellsIterated;
  uint pinXOffset, pinYOffset;

  /* Compute the number of nets and pins */  
  numPins = 0;
  numCellsIterated = 0;
  DESIGN_FOR_ALL_CELLS(myDesign, cellName, cellPtr) {
    numPins += (*cellPtr).CellGetNumPins();
    numCellsIterated++;
  } DESIGN_END_FOR;

  cout << "Read " << numPins << " pins from " 
       << numCellsIterated << " cells" << endl;
  cout << "Writing nets file.." << endl;

  if (fname == "") {
    designName = myDesign.DesignGetName();
    fileName = designName + ".nets";
  } else {
    fileName = fname + ".nets";
  }

  opFile.open(fileName.data(), ifstream::out);
  opFile << "UCLA nets 1.0" << endl; 
  DesignWriteHeaderFile(myDesign, opFile);
  opFile << endl << "NumNets : " << myDesign.DesignGetNumNets() << endl;
  opFile << "NumPins : " << numPins << endl;
  
  opFile << endl;

  DESIGN_FOR_ALL_NETS(myDesign, netName, netPtr) {
    numNetPins = 0;
    string netString;
    NET_FOR_ALL_PINS((*netPtr), pinPtr) {
      if ((*pinPtr).PinGetDirection() == PIN_DIR_INPUT) {
	pinDir = "I";
      } else {
	pinDir = "O";
      }
      Cell &cellOfPin = (*pinPtr).PinGetParentCell();
      cellName = cellOfPin.CellGetName();
      pinXOffset = (*pinPtr).PinGetXOffset();
      pinYOffset = (*pinPtr).PinGetYOffset();
      pinXOffset -= (cellOfPin.CellGetWidth() / 2);
      pinYOffset -= (cellOfPin.CellGetHeight() / 2);
      netString += "\t" + cellName + "  " + pinDir + "  : " +
	getStrFromInt(pinXOffset) + "  " + getStrFromInt(pinYOffset) + " \n";
      numNetPins++;
    } NET_END_FOR;
    opFile << "NetDegree : " << numNetPins << "  " << netName << endl;
    opFile << netString;
  } DESIGN_END_FOR;

  opFile.close();
}

void DesignWriteNodes(Design &myDesign, string fname) 
{
  string cellName, designName, fileName;
  Cell *cellPtr;
  uint cellCount;
  ofstream opFile;

  cout << "Writing nodes file.." << endl;
  if (fname == "") {
    designName = myDesign.DesignGetName();
    fileName = designName + ".nodes";
  } else {
    fileName = fname + ".nodes";
  }
  cellCount = 0;
  DESIGN_FOR_ALL_CELLS(myDesign, cellName, cellPtr) {
    cellCount++;
  } DESIGN_END_FOR;
  
  opFile.open(fileName.data(), ifstream::out);
  opFile << "UCLA nodes 1.0" << endl;
  DesignWriteHeaderFile(myDesign, opFile);
  opFile << endl << "NumNodes :              " << cellCount << endl;
  opFile << "NumTerminals :           " << myDesign.DesignGetNumTerminalCells() << endl;

  DESIGN_FOR_ALL_CELLS(myDesign, cellName, cellPtr) {
    opFile << "        " << cellName << " " << (*cellPtr).CellGetWidth();
    opFile << "      " << (*cellPtr).CellGetHeight(); 
    if ((*cellPtr).CellIsTerminal()) {
      opFile << "     " << "terminal";
    }
    opFile << endl;
  } DESIGN_END_FOR;

  opFile.close();
}

void
DesignWritePlacement(Design &myDesign, string fname) 
{
  string designName, fileName;
  string cellName;
  Cell *cellPtr;
  ofstream opFile;

  if (fname == "") {
    designName = myDesign.DesignGetName();
    fileName = designName + ".pl";
  } else {
    fileName = fname + ".pl";
  }

  _STEP_BEGIN("Writing placement for current design");
  opFile.open(fileName.data(), ifstream::out);

  opFile << "UCLA pl 1.0" << endl;
  DesignWriteHeaderFile(myDesign, opFile);
  opFile << endl;

  DESIGN_FOR_ALL_CELLS(myDesign, cellName, cellPtr) {
    uint cellXpos;
    uint cellYpos;
    cellXpos = (*cellPtr).CellGetXpos();
    cellYpos = (*cellPtr).CellGetYpos();
    opFile << cellName << "\t" << cellXpos << "\t" << cellYpos 
	   << "\t:\t" << getStrForOrientation((*cellPtr).CellGetOrientation());
    if ((*cellPtr).CellIsTerminal() && !(*cellPtr).CellIsPort()) {
      opFile << "\t" << "/FIXED";
    }
    opFile << endl;
  } DESIGN_END_FOR;

  _STEP_END("Writing placement for current design");
  
  opFile.close();
}

void
DesignWritePlacementFP(Design &myDesign, string fname) 
{
  string designName, fileName;
  string cellName;
  Cell *cellPtr;
  ofstream opFile;

  if (fname == "") {
    designName = myDesign.DesignGetName();
    fileName = designName + ".pl";
  } else {
    fileName = fname + ".pl";
  }

  _STEP_BEGIN("Writing placement for current design");
  opFile.open(fileName.data(), ifstream::out);

  opFile << "UCLA pl 1.0" << endl;
  DesignWriteHeaderFile(myDesign, opFile);
  opFile << endl;

  DESIGN_FOR_ALL_CELLS(myDesign, cellName, cellPtr) {
    uint cellXpos;
    uint cellYpos;
    cellXpos = (*cellPtr).CellGetXpos();
    cellYpos = (*cellPtr).CellGetYpos();
    opFile << cellName << "\t" << cellXpos << "\t" << cellYpos 
	   << "\t:\t" << getStrForOrientation((*cellPtr).CellGetOrientation());
    if ((*cellPtr).CellIsTerminal()) {
      opFile << "\t" << "/FIXED";
    }
    opFile << endl;
  } DESIGN_END_FOR;

  _STEP_END("Writing placement for current design");
  opFile.close();
}

void
DesignWriteScl(Design &myDesign, string fname)
{
  map<uint, uint> subRows;
  string designName, fileName;
  string cellName;
  PhysRow *rowPtr;
  ofstream opFile;
  uint rowIdx;
  uint subRowOrigin, numSites;

  if (fname == "") {
    designName = myDesign.DesignGetName();
    fileName = designName + ".scl";
  } else {
    fileName = fname + ".scl";
  }

  _STEP_BEGIN("Writing rows for current design");
  opFile.open(fileName.data(), ifstream::out);

  opFile << "UCLA scl 1.0" << endl;
  DesignWriteHeaderFile(myDesign, opFile);
  opFile << endl;
  opFile << "NumRows : " << myDesign.DesignGetNumPhysRows() << endl;
  opFile << endl;

  DESIGN_FOR_ALL_ROWS(myDesign, rowIdx, rowPtr) {
    opFile << "CoreRow " << (*rowPtr).PhysRowGetTypeStr() << endl;
    opFile << "  Coordinate    :   " << (*rowPtr).PhysRowGetCoordinate() << endl;
    opFile << "  Height        :   " << (*rowPtr).PhysRowGetHeight() << endl;
    opFile << "  Sitewidth     :   " << (*rowPtr).PhysRowGetSiteWidth() << endl;
    opFile << "  Sitespacing   :   " << (*rowPtr).PhysRowGetSiteSpacing() << endl;
    opFile << "  Siteorient    :   " << (*rowPtr).PhysRowGetSiteOrientationStr() << endl;
    opFile << "  Sitesymmetry  :   " << (*rowPtr).PhysRowGetSiteSymmetryStr() << endl;
    subRows = (*rowPtr).PhysRowGetSubRows();
    MAP_FOR_ALL_ELEMS(subRows, uint, uint, subRowOrigin, numSites) {
      opFile << "  SubrowOrigin  :   " << subRowOrigin 
	     << "  NumSites  :  " << numSites << endl;
    } END_FOR;
    opFile << "End" << endl;
  } DESIGN_END_FOR;

  _STEP_END("Writing rows for current design");

  opFile.close();
}

void
DesignWriteWtsFile(Design &myDesign, string fname) 
{
  string designName, fileName;
  ofstream opFile;

  _STEP_BEGIN("Writing weights file");

  if (fname == "") {
    designName = myDesign.DesignGetName();
    fileName = designName + ".wts";
  } else {
    fileName = fname + ".wts";
  }

  opFile.open(fileName.data(), ifstream::out);
  opFile << "UCLA wts 1.0" << endl;
  DesignWriteHeaderFile(myDesign, opFile);

  _STEP_END("Writing weights file");

  opFile.close();
}

void
DesignWriteAuxFile(Design &myDesign, string fname)
{
  string designName, fileName;
  ofstream opFile;

  _STEP_BEGIN("Writing aux file");

  if (fname == "") {
    designName = myDesign.DesignGetName();
    fileName = designName + ".aux";
  } else {
    designName = fname;
    fileName = fname + ".aux";
  }

  opFile.open(fileName.data(), ifstream::out);
  opFile << "RowBasedPlacement : " 
	 << (designName + ".nodes ") 
	 << (designName + ".nets ")
	 << (designName + ".wts ")
	 << (designName + ".pl ")
	 << (designName + ".scl ")
	 << endl;
  _STEP_END("Writing aux file");

  opFile.close();
}

void DesignWriteBookShelfOutput(Design &myDesign, string opBenchName)
{
  /* Write the nets, nodes and pl file */
  DesignWriteNets(myDesign, opBenchName);
  DesignWriteNodes(myDesign, opBenchName);
  DesignWritePlacement(myDesign, opBenchName);
  DesignWriteScl(myDesign, opBenchName);

  /* Create the wts file */
  DesignWriteWtsFile(myDesign, opBenchName);
  /* Create the aux file */
  DesignWriteAuxFile(myDesign, opBenchName);
}

void DesignWriteBookShelfOutput(Design &myDesign, string opBenchName,
				bool forMacroPlacer)
{
  /* Write the nets, nodes and pl file */
  DesignWriteNets(myDesign, opBenchName);
  DesignWriteNodes(myDesign, opBenchName);
  if (forMacroPlacer) {
    DesignWritePlacementFP(myDesign, opBenchName);
  } else {
    DesignWritePlacement(myDesign, opBenchName);
  }
  DesignWriteScl(myDesign, opBenchName);

  /* Create the wts file */
  DesignWriteWtsFile(myDesign, opBenchName);
  /* Create the aux file */
  DesignWriteAuxFile(myDesign, opBenchName);
}

void DesignWriteOutputPlacement(Design &myDesign)
{
  DesignWritePlacement(myDesign, (myDesign.DesignGetName() + ".ourplacer"));
}

void DesignWriteOutputPlacement(Design &myDesign, string outputFileName)
{
  DesignWritePlacement(myDesign, outputFileName);
}

void 
printAllVisibleCellsInDesign (Design& myDesign,string fname)
{
   string fileName;
   string cellName;
   Cell *cellPtr;
   ofstream opFile;
   fileName = fname + "_debug.pl";
   opFile.open(fileName.data(), ifstream::out);
   opFile << "UCLA pl 1.0" << endl;
   DesignWriteHeaderFile(myDesign, opFile);
   opFile << endl;
   DESIGN_FOR_ALL_CELLS(myDesign, cellName, cellPtr) {
           uint cellXpos;
           uint cellYpos;
           uint clusterLevel;
           clusterLevel = (*cellPtr).CellGetClusterLevel();    
           cellXpos = (*cellPtr).CellGetXpos();
           cellYpos = (*cellPtr).CellGetYpos();
           opFile << cellName << "\t" << cellXpos << "\t" << cellYpos 
           << "\t:\t" << getStrForOrientation((*cellPtr).CellGetOrientation());
          if ((*cellPtr).CellIsTerminal() && !(*cellPtr).CellIsPort()) {
               opFile << "\t" << "/FIXED";
          }   
          if ((*cellPtr).CellIsCluster()){
          opFile << "\t" << "/cluster";
          }   
          opFile << "\t CL :" <<clusterLevel;    
          opFile << endl;
   } DESIGN_END_FOR;
  opFile.close();
}

void 
printVisibleCellsineachCluster ( Design& myDesign, string fname)
{
        string fileName;
        fileName = fname + "_debug.pl";
        ofstream opFile;
        opFile.open(fileName.data(), ifstream::out);
        opFile << "UCLA pl 1.0" << endl;
        DesignWriteHeaderFile(myDesign, opFile);
        opFile << endl;
  
  
 /* failed attempt to get all the cells which are visible in a cluster
        string clusterName;
        string cellName;
        Cluster *clusterPtr;
        Cell *clusterCellPtr;
        DESIGN_FOR_ALL_CLUSTERS(myDesign,clusterName,clusterCellPtr) {
                 clusterPtr = (Cluster *)CellGetCluster(clusterCellPtr);
                 opFile << "cluster Name: " << clusterName << endl; 
                 vector<Cell *> &cellsInCluster = (*clusterPtr).ClusterGetCellsOfCluster();
                 Cell *cellPtr;
                 VECTOR_FOR_ALL_ELEMS (cellsInCluster,Cell*,cellPtr){
                        //if ((*cellPtr).CellIsHidden()) continue;
                        cellName=(*cellPtr).CellGetName();
                        opFile << "\t" << cellName << endl;
                 }END_FOR;
        }DESIGN_END_FOR;*/


// New attempt to get the cells which are visible to top level nets
        string netName;
        Net *netPtr;
        DESIGN_FOR_ALL_NETS(myDesign,netName,netPtr) {
                opFile << "Net Name: " << netName << endl;
                Cell *cellPtr;
                string cellName;
                Cell *parentClusterPtr;
                string parentClusterName;
                NET_FOR_ALL_CELLS((*netPtr),cellPtr){
                        cellName = (*cellPtr).CellGetName();
                        opFile << "\t" << cellName;
                        if (((*cellPtr).CellIsClusterChild())){
                                parentClusterPtr = (*cellPtr).CellGetParentCluster();
                                parentClusterName = (*parentClusterPtr).CellGetName();
                                opFile << "\t" << "parentCluster: " << parentClusterName << endl;

                        }else{
                        opFile << "\t" << "Not Cluster child \t" << endl;
                        }
                }NET_END_FOR;
        }DESIGN_END_FOR;

// The above Worked to display all the cells connected to top nets. The below subroutine is to compare cell list generated by nets of cluster//

//The below sub routine was taken from "createFDPNetlist". however the net,cell information in the routine above and the one below is not matching
// need to check this concept with Nakul or Tuhin


        string fileName2;
        fileName2 = fname + "_debug2.pl";
        ofstream opFile2;
        opFile2.open(fileName2.data(), ifstream::out);
        opFile2 << "UCLA pl 1.0" << endl;
        DesignWriteHeaderFile(myDesign, opFile2);
        opFile2 << endl;
        Cell *clusterCellPtr;
        Cell *cellPtr2;
        Net *netPtr2;
        string clusterCellName;
        DESIGN_FOR_ALL_CLUSTERS(myDesign,clusterCellName,clusterCellPtr){
                string netName2;
                map<Net *,bool> visitedNets;
                CELL_FOR_ALL_NETS_NO_DIR((*clusterCellPtr),netPtr2){
                        _KEY_EXISTS(visitedNets,netPtr2){
                                continue;
                        } else {
                                visitedNets[netPtr2] = true;
                        }
                        string cellName2;
                        netName2 = (*netPtr2).NetGetName();
                        opFile2 << "NetName: " << netName2 << endl ;
                        NET_FOR_ALL_CELLS((*netPtr2),cellPtr2){
                                if (cellPtr2 == clusterCellPtr) continue;
                                cellName2 = (*cellPtr2).CellGetName();
                                opFile2 << "\t" << cellName2 << endl;
                        }NET_END_FOR;
                }CELL_END_FOR;
        }DESIGN_END_FOR;

        string filename3;
        filename3 = fname + "_debug3.pl";
        ofstream opFile3;
        opFile3.open(filename3.data(),ifstream::out);
        opFile3 <<"UCLS pl 1.0" << endl;
        DesignWriteHeaderFile(myDesign, opFile3);
        opFile3 << endl;
        Cell *cellPtr3;
        Net *netPtr3;
        string netName3;
       
       //Above Iteration only gives Cluster cell names. Need to try with pins  
        DESIGN_FOR_ALL_NETS(myDesign,netName3,netPtr3) {
                opFile3 << "Net Name: " << netName3 << endl;
                Cell *cellPtr3;
                string cellName3;
                Pin *pinPtr3;
                NET_FOR_ALL_PINS ((*netPtr3),pinPtr3){
                        cellPtr3 = (*pinPtr3).PinGetParentCellPtr ();
                        cellName3 = (*cellPtr3).CellGetName();
                        opFile3 << "\tparentCell :" << cellName3 << endl;
                }NET_END_FOR;
        }DESIGN_END_FOR;
        // The pins after clustering gives Cluster names only. need to figure out a way to get Cell names 
}


        
