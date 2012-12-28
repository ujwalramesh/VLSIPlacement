# include <Design.h>

using namespace std;

void
Design::DesignProcessProperty(ifstream& file, string& Property, 
			      string& Value)
{
  string line;
  unsigned int numNodes, numTerms;
  int idx1, idx2, idx3;
  bool result;

  if (!file.eof()) getline(file, line);

  while (line == "" && !file.eof()) 
    getline(file, line);

  idx1 = 0;
  while (line[idx1] == ' ' || line[idx1] == '\t') 
    idx1++;

  idx2 = idx1;
  while (line[idx2] != ' ' && line[idx2] != '\t') 
    idx2++;

  idx3 = idx2;
  while (!(line[idx3] >= ZERO && line[idx3] <= NINE)) 
    idx3++;
  
  Property = line.substr(idx1, (idx2 - idx1));
  Value = line.substr(idx3);
}


void
Design::DesignFileReadHeader(ifstream& file)
{
  string buffer;
  int i;

  i = 0;
  while (!file.eof() && i < 4) {
    getline(file, buffer);
    i++;
  }
}

bool
DesignCellIsStdCell(Design &myDesign, Cell &thisCell)
{
  map<unsigned int, unsigned int>rowHeights;
  int rowHeight, cellHeight;
  bool stdCell;

  rowHeights = myDesign.DesignGetRowHeights();
  rowHeight = myDesign.DesignGetSingleRowHeight();
  cellHeight = thisCell.CellGetHeight();
  stdCell = false;

  if (rowHeight == -1) {
    if (rowHeights.find(cellHeight) != rowHeights.end()) {
      stdCell = true;
    }
  } else if (cellHeight == rowHeight) {
    stdCell = true;
  }

  return (stdCell);
}
