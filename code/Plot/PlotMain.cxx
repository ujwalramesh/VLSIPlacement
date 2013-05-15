# include <Plot.h>

Plot::Plot()
{
}

Plot::~Plot()
{
}

Plot::Plot(string Title)
{
  PlotSetTitle(Title);
}

Plot::Plot(string Title, string plotFileName)
{
  PlotSetTitle(Title);
  PlotSetFileName(plotFileName);
}

Plot::Plot(double minx, double miny, double maxx, double maxy)
{
  PlotSetBoundary(minx, miny, maxx, maxy);
}

void
Plot::PlotSetTitle(string Title)
{
  this->plotTitle = Title;
}

void
Plot::PlotSetFileName(string plotFileName)
{
  this->plotFileName = plotFileName;
}

void
Plot::PlotSetBoundary(double minx, double miny, double maxx,
		      double maxy)
{
  this->minx = minx;
  this->miny = miny;
  this->maxx = maxx;
  this->maxy = maxy;
}

void
Plot::PlotSetBoundary(Design& myDesign)
{
  uint minx, miny, maxx, maxy;
  
  minx = 0; miny = 0;
  myDesign.DesignGetBoundingBox(maxx, maxy);
  PlotSetBoundary(minx, miny, maxx, maxy);
}

string
Plot::PlotGetTitle(void)
{
  return (this->plotTitle);
}

string
Plot::PlotGetFileName(void)
{
  return (this->plotFileName);
}

void
Plot::PlotGetBoundary(double &minx, double &miny, double &maxx,
		      double &maxy)
{
  minx = this->minx;
  miny = this->miny;
  maxx = this->maxx;
  maxy = this->maxy;
}

void
Plot::PlotAddLine(double x1, double y1, double x2, double y2)
{
  Line newLine(x1, y1, x2, y2);
  
  lines.push_back(newLine);
}

void
Plot::PlotAddRow(PhysRow &thisRow)
{
  vector<int> bbox;
  double left, bot, right, top;
  
  thisRow.PhysRowGetBoundingBox(bbox);
  
  left = bbox[0]; bot = bbox[1]; 
  right = bbox[2]; top = bbox[3];
  
  Rect RowRect(left, bot, right, top);
  
  normalRects.push_back(RowRect);
}

void
Plot::PlotAddCell(Cell &thisCell)
{
  double left, bot, right, top;
  
  left = thisCell.CellGetXpos();
  bot = thisCell.CellGetYpos();
  right = left + thisCell.CellGetWidth();
  top = bot + thisCell.CellGetHeight();

  Rect CellRect(left, bot, right, top);
  
  normalRects.push_back(CellRect);
}

void
Plot::PlotAddCells(vector<Cell *> &cellsToPlot)
{
  Cell *cellPtr;

  VECTOR_FOR_ALL_ELEMS(cellsToPlot, Cell *, cellPtr) {
    PlotAddCell(*cellPtr);
  } END_FOR;
}

void
Plot::PlotHighlightCell(Cell &thisCell)
{
  double left, bot, right, top;
  
  left = thisCell.CellGetXpos();
  bot = thisCell.CellGetYpos();
  right = left + thisCell.CellGetWidth();
  top = bot + thisCell.CellGetHeight();

  Rect CellRect(left, bot, right, top);
  
  highlightedRects.push_back(CellRect);
}

void
Plot::PlotAddPort(Cell &thisCell)
{
  double left, bot, right, top;
  
  left = thisCell.CellGetXpos();
  bot = thisCell.CellGetYpos();
  right = left + thisCell.CellGetWidth();
  top = bot + thisCell.CellGetHeight();

  Rect CellRect(left, bot, right, top);
  
  ports.push_back(CellRect);
}

void
Plot::PlotAddBin(Bin &thisBin)
{
  double left, bot, right, top;
  string rectLabel;

  left = thisBin.BinGetLeft();
  bot = thisBin.BinGetBot();
  right = thisBin.BinGetRight();
  top = thisBin.BinGetTop();
  
  Rect BinRect(left, bot, right, top);
  BinRect.RectSetLabel(getStrFromInt(thisBin.BinGetIdx()));
  BinRect.RectSetLabelPos(RIGHT_TOP);

  bins.push_back(BinRect);
}

void
Plot::PlotAddStretchedBin(Bin &thisBin)
{
  double left, bot, right, top;

  left = thisBin.BinGetLeft();
  bot = thisBin.BinGetBot();
  right = thisBin.BinGetNewRight();
  top = thisBin.BinGetNewTop();
  
  Rect BinRect(left, bot, right, top);
  
  stretchedBins.push_back(BinRect);
}

void
Plot::PlotWriteOutput(void)
{
  if (!plotOpFile.is_open()) {
    plotOpFile.open(plotFileName.data(), ios_base::out);
  }
  
  /* Start writing output into the gnuplot plot file */
  plotOpFile << "########################################################" << endl;
  plotOpFile << "#  GNUPLOT FILE WRITTEN BY OURPLACER                   #" << endl;
  plotOpFile << "########################################################" << endl;
  plotOpFile << endl;
  plotOpFile << "set title \"" << plotTitle << "\" font \"Times, 18\"" << endl;
  plotOpFile << "set size ratio 1" << endl;
  plotOpFile << "set nokey" << endl;
  plotOpFile << endl;
  
  plotOpFile << "########################################################" << endl;
  plotOpFile << "#  NORMAL RECTANGLES WILL BE PLOTTED BELOW             #" << endl;
  plotOpFile << "########################################################" << endl;
  
  Rect myRect;
  VECTOR_FOR_ALL_ELEMS(normalRects, Rect, myRect) {
    myRect.RectSetStyle(NORMAL_RECT);
    myRect.RectWriteOutput(plotOpFile);
  } END_FOR;

  VECTOR_FOR_ALL_ELEMS(highlightedRects, Rect, myRect) {
    myRect.RectSetStyle(HIGHLIGHTED_RECT);
    myRect.RectWriteOutput(plotOpFile);
  } END_FOR;

  VECTOR_FOR_ALL_ELEMS(ports, Rect, myRect) {
    myRect.RectSetStyle(PORT_RECT);
    myRect.RectWriteOutput(plotOpFile);
  } END_FOR;

  VECTOR_FOR_ALL_ELEMS(bins, Rect, myRect) {
    myRect.RectSetStyle(BIN_RECT);
    myRect.RectWriteOutput(plotOpFile);
  } END_FOR;
  
  VECTOR_FOR_ALL_ELEMS(stretchedBins, Rect, myRect) {
    myRect.RectSetStyle(STRETCHED_BIN_RECT);
    myRect.RectWriteOutput(plotOpFile);
  } END_FOR;

  Line myLine;
  VECTOR_FOR_ALL_ELEMS(lines, Line, myLine) {
    myLine.LineSetStyle(NORMAL_LINE);
    myLine.LineWriteOutput(plotOpFile);
  } END_FOR;

  VECTOR_FOR_ALL_ELEMS(highlightedLines, Line, myLine) {
    myLine.LineSetStyle(HIGHLIGHTED_LINE);
    myLine.LineWriteOutput(plotOpFile);
  } END_FOR;

  
  plotOpFile << "plot[" << minx << ":" << maxx << "][" << miny << ":" << maxy << "]" 
	     << " '-' w l lt 3" << endl;
  plotOpFile << "# draw the outer chip boundary #" << endl;
  plotOpFile << minx << " " << miny << endl;
  plotOpFile << minx << " " << maxy << endl;
  plotOpFile << maxx << " " << maxy << endl;
  plotOpFile << endl;
  plotOpFile << minx << " " << miny << endl;
  plotOpFile << maxx << " " << miny << endl;
  plotOpFile << maxx << " " << maxy << endl;
  
  plotOpFile << endl << "EOF" << endl;

  plotOpFile << endl << "pause -1 'Press any key to exit'" << endl;
  plotOpFile.close();
}

string 
LineGetStyleString(Line &thisLine)
{
  return (" lt 2");
}

string 
RectGetStyleString(Rect &thisRect)
{
  lineStyle rectStyle;
  string styleString;

  rectStyle = thisRect.RectGetStyle();
  switch (rectStyle) {
  case NORMAL_RECT: 
    styleString = "fs empty border lc rgb \"blue\"";
    break;
  case HIGHLIGHTED_RECT: 
    styleString = "fs empty border lc rgb \"red\"";
    break;
  case BIN_RECT: 
    styleString = "fs empty border lc rgb \"gray\"";
    break;
  case STRETCHED_BIN_RECT: 
    styleString = "fs empty border lc rgb \"green\"";
    break;
  default: 
    break;
  };
  
  return (styleString);
}
