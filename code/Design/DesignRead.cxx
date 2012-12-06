# include <Design.h>
# include <assert.h>

void
Design::DesignReadDesign(string DesignPath, string DesignName)
{
  string fullFileName, fileName;
  string listOfFiles;
  string fullDesignPath;
  bool rowBasedPlacement;
  vector<string> fileList;
  
  fullDesignPath = DesignPath + DIR_SEP + DesignName + DIR_SEP;
  fullFileName = fullDesignPath + DesignName + DESIGN_AUX_FILE_EXT;

  /* Open the aux fle and read all the other files that need to be
     read */
  ifstream DesignFile(fullFileName.data(), ios_base::in);
  getline(DesignFile, listOfFiles);
  if (listOfFiles == "") {
    return;
  }

  /* Collect all files to be read from the aux file */
  istringstream stream(listOfFiles, istringstream::in);
  do {
    stream >> fileName;
    if (fileName == ROW_BASED_PLACEMENT) {
      RowBasedPlacement = true;
    } else if (fileName.find('.', 0) == string::npos) {
      continue;
    }
    fileList.push_back(fileName);
  } while (stream);
  
  VECTOR_FOR_ALL_ELEMS(fileList, string, fileName) {
    int dotpos = fileName.find('.', 0);
    _ASSERT("Filename of unexpected type", dotpos < 0);

    string ext = fileName.substr(dotpos);
    if (ext == DESIGN_CELL_FILE_EXT) {
      DesignCellFileName = fullDesignPath + fileName;
    } else if (ext == DESIGN_NET_FILE_EXT) {
      DesignNetFileName = fullDesignPath + fileName;
    } else if (ext == DESIGN_NET_WTS_FILE_EXT) {
      DesignNetWtsFileName = fullDesignPath + fileName;
    } else if (ext == DESIGN_SCL_FILE_EXT) {
      DesignSclFileName = fullDesignPath + fileName;
    } else if (ext == DESIGN_PL_FILE_EXT) {
      DesignPlFileName = fullDesignPath + fileName;
    } else {
      _ASSERT("Unknown file type found in aux file", 1);
    }
  } END_FOR;

  DesignReadCells();
  DesignReadNets();
  DesignReadRows();
  DesignReadCellPlacement();
  
  DesignFile.close();
}

void
Design::DesignFileReadOneNode(ifstream& file)
{
  Cell *newCell;
  string NodeName, terminalString;
  string line;
  unsigned int height, width;
  bool terminalCell;

  while (!file.eof()) {
    getline(file, line);
    if (line == "") {
      continue;
    }

    istringstream stream(line, istringstream::in);

    stream >> NodeName; stream >> height; stream >> width;
    if (!stream.eof()) stream >> terminalString;

    if (terminalString == NODE_TERMINAL_KEYWORD) {
      terminalCell = true;
    }
    newCell = new Cell(height, width, NodeName, terminalCell);
    DesignAddOneCellToDesignDB(newCell);
    break;
  }
}


void
Design::DesignFileReadNodes(ifstream& file)
{
  string Property, Value;
  unsigned int numNodes, numTerminals;
  int idx, cellCount;

  /* Skip empty lines */
  
  for (idx = 0; idx < NUM_NODE_PROPERTIES; idx++) {
    DesignProcessProperty(file, Property, Value);
    if (Property == NUM_NODES_PROPERTY) {
      numNodes = atoi(Value.data());
    } else if (Property == NUM_TERMINALS_PROPERTY) {
      numTerminals = atoi(Value.data());
    }
  }

  while (!file.eof()) {
    DesignFileReadOneNode(file);
  }
}


void
Design::DesignReadCells()
{
  int NumCells;
  string Msg;

  DesignOpenFile(DesignCellFileName);
  DesignFileReadHeader(DesignFile);
  DesignFileReadNodes(DesignFile);

  Msg += "Loaded " + getStrFromInt(DesignCells.size()) + " cells";

  common_message(Msg);
}


void 
Design::DesignFileReadPins(ifstream &file, unsigned int netDegree,
			   Net &newNet) 
{
  string nodeName, pinDir, garbage;
  string PinName, Msg, line;
  string NetName;
  unsigned int pinNum;
  float xoffset, yoffset;
  char dir;

  pinNum = 1;
  NetName = newNet.NetGetName();
  for (int i=0; i < netDegree; i++, pinNum++) {
    if (file.eof()) {
      break;
    }
    getline(file, line);
    
    Cell *node;
    istringstream stream(line, istringstream::in);
    stream >> nodeName; stream >> pinDir; 
    stream >> garbage;  stream >> xoffset; 
    stream >> yoffset;

    if (pinDir == PIN_DIR_INPUT_STRING) dir = PIN_DIR_INPUT;
    else if (pinDir == PIN_DIR_OUTPUT_STRING) dir = PIN_DIR_OUTPUT;

    PinName = nodeName + "_" + NetName + "_" + getStrFromInt(pinNum);
    Pin *newPin = new Pin(0, (int)xoffset, (int)yoffset, dir, *node, PinName);

    Msg = "Created Pin " + PinName + " on cell " + nodeName;
    common_message(Msg);                                   
    (*newPin).Connect(newNet);
    node = DesignGetNode(nodeName);
    (*newPin).PinSetParentCell(*node);
    (*node).CellAddPin(newPin);
    newNet.NetAddPin(*newPin);
  }

  Msg = "Read and created " + getStrFromInt((pinNum - 1)) + " pins";
  common_message(Msg);                                   
}


void
Design::DesignFileReadOneNet(ifstream &file)
{
  string line, garbage;
  string NetDegree, NetName;
  unsigned int netDegree;
  Cell *node;
  Net *newNet;

  while (!file.eof()) {
    getline(file, line);
    if (line == "") {
      continue;
    }

    istringstream stream(line, istringstream::in);
    stream >> NetDegree;

    if (NetDegree != NET_DEGREE_PROPERTY) {
      string Msg = "Error in benchmark file.";
      common_error(Msg);
    }

    stream >> garbage; stream >> netDegree;
    stream >> NetName;
    Net *newNet;
    string Msg;

    newNet = new Net(0, NetName);
    DesignNets[NetName] = newNet;

    Msg = "Created Net " + NetName + " of Degree " + getStrFromInt(netDegree);
    DesignFileReadPins(file, netDegree, *newNet);
    DesignAddOneNetToDesignDB(newNet);
    break;
  }
}


void
Design::DesignFileReadNets(ifstream& file)
{
  string Property, Value;
  unsigned int numNets, numPins, numNodes;
  int idx, netCount;

  for (idx = 0; idx < NUM_NET_PROPERTIES; idx++) {
    DesignProcessProperty(file, Property, Value);
    if (Property == NUM_NETS_PROPERTY) {
      numNodes = atoi(Value.data());
    } else if (Property == NUM_PINS_PROPERTY) {
      numPins = atoi(Value.data());
    }
  }

  while (!file.eof()) {
    DesignFileReadOneNet(file);
  }
}


void
Design::DesignReadNets()
{
  int NumCells;
  string Msg;

  DesignOpenFile(DesignNetFileName);
  DesignFileReadHeader(DesignFile);
  DesignFileReadRows(DesignFile);

  Msg += "Loaded " + getStrFromInt(DesignNets.size()) + " nets";

  common_message(Msg);
}

void
Design::DesignFileReadOneRow(ifstream &file)
{
  string rowProperty, garbage, line;
  bool rowBegin;
  int rowCoordinate, subRowOrigin;
  unsigned int height, siteWidth, siteSpacing, numSites;
  rowOrientation rowType;
  siteOrientation siteOrient;
  siteSymmetry symmetry;
  vector<int> subRowOrigins;
  vector<unsigned int> numSitesForSubRows;

  rowBegin = false;
  while (!file.eof()) {
    getline(file, line);
    if (line == "") {
      continue;
    }
    istringstream stream(line, istringstream::in);
    stream >> rowProperty;
    if (rowProperty == ROW_BEGIN_KEYWORD) {
      if (rowBegin == true) _ASSERT_TRUE("Benchmark error");
      rowBegin = true;
      stream >> rowProperty;
      rowType = PhysRowGetRowTypeFromStr(rowProperty);
    } else if (rowProperty == ROW_COORDINATE) {
      if (rowBegin == false) _ASSERT_TRUE("Benchmark error");
      stream >> garbage;
      if (garbage != ":") _ASSERT_TRUE("Benchmark error: Missing ':'");
      stream >> rowCoordinate;
    } else if (rowProperty == ROW_HEIGHT) {
      if (rowBegin == false) _ASSERT_TRUE("Benchmark error");
      stream >> garbage;
      if (garbage != ":") _ASSERT_TRUE("Benchmark error: Missing ':'");
      stream >> height;
    } else if (rowProperty == ROW_SITE_WIDTH) {
      if (rowBegin == false) _ASSERT_TRUE("Benchmark error");
      stream >> garbage;
      if (garbage != ":") _ASSERT_TRUE("Benchmark error: Missing ':'");
      stream >> siteWidth;
    } else if (rowProperty == ROW_SITE_SPACING) {
      if (rowBegin == false) _ASSERT_TRUE("Benchmark error: Row begin unspecified");
      stream >> garbage;
      if (garbage != ":") _ASSERT_TRUE("Benchmark error: Missing ':'");
      stream >> siteSpacing;
    } else if (rowProperty == ROW_SITE_ORIENTATION) {
      if (rowBegin == false) _ASSERT_TRUE("Benchmark error: Row begin unspecified");
      stream >> garbage;
      if (garbage != ":") _ASSERT_TRUE("Benchmark error: Missing ':'");
      stream >> rowProperty;
      siteOrient = PhysRowGetSiteOrientationFromStr(rowProperty);
    } else if (rowProperty == ROW_SITE_SYMMETRY) {
      if (rowBegin == false) _ASSERT_TRUE("Benchmark error: Row begin unspecified");
      stream >> garbage;
      if (garbage != ":") _ASSERT_TRUE("Benchmark error: Missing ':'");
      stream >> rowProperty;
      symmetry = PhysRowGetSiteSymmetryFromStr(rowProperty);
    } else if (rowProperty == SUBROW_ORIGIN) {
      if (rowBegin == false) _ASSERT_TRUE("Benchmark error: Row begin unspecified");
      stream >> garbage;
      if (garbage != ":") _ASSERT_TRUE("Benchmark error: Missing ':'");
      stream >> subRowOrigin;
      stream >> garbage;
      if (garbage != ":") _ASSERT_TRUE("Benchmark error: Missing ':'");
      stream >> numSites;
      subRowOrigins.push_back(subRowOrigin);
      numSitesForSubRows.push_back(numSites);
    } else if (rowProperty == ROW_END_KEYWORD) {
      if (rowBegin == false) _ASSERT_TRUE("Benchmark error: Row begin unspecified");
      rowBegin = false;
    } else {
      _ASSERT_TRUE("Benchmark error: Unknown keyword");
    }
    break;
  }
  PhysRow *row;
  row = new PhysRow(rowType, rowCoordinate, height, siteWidth, siteSpacing,
		    siteOrient, symmetry);

  VECTOR_FOR_ALL_ELEMS_DOUBLE(subRowOrigins, int, subRowOrigin,
			      numSitesForSubRows, unsigned int, numSites) {
    (*row).PhysRowAddSubRow(subRowOrigin, numSites);
  } END_FOR;
  DesignAddOnePhysRowToDesignDB(row);
}

void
Design::DesignFileReadRows(ifstream& file)
{
  string Property, Value;
  unsigned int numRows
  int idx, netCount;

  for (idx = 0; idx < NUM_ROWS_PROPERTIES; idx++) {
    do {
      DesignProcessProperty(file, Property, Value);
    } while (Property != "" && !file.eof());
    if (Property == NUM_ROWS_PROPERTY) {
      numRows = atoi(Value.data());
    }
  }
  for (idx = 0; idx < numRows; idx++) {
    DesignFileReadOneRow(file);
  }
}

void
Design::DesignReadRows()
{
  int NumRows;
  string Msg;

  DesignOpenFile(DesignSclFileName);
  DesignFileReadHeader(DesignFile);
  DesignFileReadPhysRows(DesignFile);

  Msg += "Loaded " + getStrFromInt(DesignPhysRows.size()) + " rows";

  common_message(Msg);
}

void 
Design::DesignReadCellPlacement()
{

}


