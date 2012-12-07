# include <PhysRow.h>

void
PhysRow::PhysRowSetCoordinate(int coordinate)
{
  this->coordinate = coordinate;
}

void
PhysRow::PhysRowSetType(rowOrientation rowType)
{
  this->rowType = rowType;
}

void
PhysRow::PhysRowSetHeight(unsigned int height)
{
  this->height = height;
}

void
PhysRow::PhysRowSetSiteWidth(unsigned int siteWidth)
{
  this->siteWidth = siteWidth;
}

void
PhysRow::PhysRowSetSiteSpacing(unsigned int siteSpacing)
{
  this->siteWidth = siteSpacing;
}

void
PhysRow::PhysRowSetSiteOrientation(siteOrientation orient)
{
  this->siteOrient = orient;
}

void
PhysRow::PhysRowSetSiteSymmetry(siteSymmetry symmetry)
{
  this->siteSym = symmetry;
}

void
PhysRow::PhysRowSetSubRows(map<unsigned int, unsigned int> subRows)
{
  this->subRows = subRows;
}

void
PhysRow::PhysRowSetNumSubRows(unsigned int numSubRows)
{
  this->numSubRows = numSubRows;
}

void
PhysRow::PhysRowSetNumSites(unsigned int numSites) 
{
  this->numSites = numSites;
}

void
PhysRow::PhysRowIncrementSubRows(void)
{
  this->numSubRows++;
}

int
PhysRow::PhysRowGetCoordinate(void)
{
  return coordinate;
}

siteOrientation
PhysRow::PhysRowGetSiteOrientation(void)
{
  return (siteOrient);
}

siteSymmetry
PhysRow::PhysRowGetSiteSymmetry(void)
{
  return siteSym;
}

rowOrientation
PhysRow::PhysRowGetType(void)
{
  return rowType;
}

unsigned int
PhysRow::PhysRowGetHeight(void)
{
  return height;
}

unsigned int
PhysRow::PhysRowGetSiteWidth(void)
{
  return siteWidth;
}

unsigned int 
PhysRow::PhysRowGetSiteSpacing(void)
{
  return siteSpacing;
}

unsigned int
PhysRow::PhysRowGetNumSubRows(void)
{
  return numSubRows;
}

map<unsigned int, unsigned int>
PhysRow::PhysRowGetSubRows(void)
{
  map<unsigned int, unsigned int> subRows;
  
  subRows = this->subRows;
  
  return subRows;
}

void
PhysRow::PhysRowAddSubRow(unsigned int rowOrigin, unsigned int numSites) 
{
  this->subRows[rowOrigin] = numSites;
  this->numSites += numSites;
}

/* Constructors begin */
PhysRow::PhysRow(rowOrientation orient)
{
  PhysRowSetType(orient);
  PhysRowSetCoordinate(DEFAULT_COORDINATE);
  PhysRowSetSiteOrientation(DEFAULT_SITE_ORIENTATION);
  PhysRowSetSiteSymmetry(DEFAULT_SITE_SYMMETRY);
  PhysRowSetHeight(DEFAULT_SITE_HEIGHT);
  PhysRowSetSiteSpacing(DEFAULT_SITE_SPACING);
  PhysRowSetSiteWidth(DEFAULT_SITE_WIDTH);
  PhysRowSetNumSites(DEFAULT_NUM_SITES);
  PhysRowSetNumSubRows(DEFAULT_NUM_SUBROWS);
  PhysRowSetSubRows(DEFAULT_SUBROWS);
}

PhysRow::PhysRow(rowOrientation orient, int coordinate)
{
  PhysRowSetType(orient);
  PhysRowSetCoordinate(coordinate);
  PhysRowSetSiteOrientation(DEFAULT_SITE_ORIENTATION);
  PhysRowSetSiteSymmetry(DEFAULT_SITE_SYMMETRY);
  PhysRowSetHeight(DEFAULT_SITE_HEIGHT);
  PhysRowSetSiteSpacing(DEFAULT_SITE_SPACING);
  PhysRowSetSiteWidth(DEFAULT_SITE_WIDTH);
  PhysRowSetNumSites(DEFAULT_NUM_SITES);
  PhysRowSetNumSubRows(DEFAULT_NUM_SUBROWS);
  PhysRowSetSubRows(DEFAULT_SUBROWS);
}

PhysRow::PhysRow(rowOrientation orient, unsigned int height)
{
  PhysRowSetType(orient);
  PhysRowSetCoordinate(DEFAULT_COORDINATE);
  PhysRowSetSiteOrientation(DEFAULT_SITE_ORIENTATION);
  PhysRowSetSiteSymmetry(DEFAULT_SITE_SYMMETRY);
  PhysRowSetHeight(DEFAULT_SITE_HEIGHT);
  PhysRowSetSiteSpacing(DEFAULT_SITE_SPACING);
  PhysRowSetSiteWidth(DEFAULT_SITE_WIDTH);
  PhysRowSetNumSites(DEFAULT_NUM_SITES);
  PhysRowSetNumSubRows(DEFAULT_NUM_SUBROWS);
  PhysRowSetSubRows(DEFAULT_SUBROWS);
}

PhysRow::PhysRow(rowOrientation orient, int coordinate, unsigned int height)
{
  PhysRowSetType(orient);
  PhysRowSetCoordinate(DEFAULT_COORDINATE);
  PhysRowSetSiteOrientation(DEFAULT_SITE_ORIENTATION);
  PhysRowSetSiteSymmetry(DEFAULT_SITE_SYMMETRY);
  PhysRowSetHeight(DEFAULT_SITE_HEIGHT);
  PhysRowSetSiteSpacing(DEFAULT_SITE_SPACING);
  PhysRowSetSiteWidth(DEFAULT_SITE_WIDTH);
  PhysRowSetNumSites(DEFAULT_NUM_SITES);
  PhysRowSetNumSubRows(DEFAULT_NUM_SUBROWS);
  PhysRowSetSubRows(DEFAULT_SUBROWS);
}

PhysRow::PhysRow(rowOrientation orient, int coordinate, unsigned int height,
		 unsigned int siteWidth, unsigned int siteSpacing) 
{
  PhysRowSetType(orient);
  PhysRowSetCoordinate(DEFAULT_COORDINATE);
  PhysRowSetSiteOrientation(DEFAULT_SITE_ORIENTATION);
  PhysRowSetSiteSymmetry(DEFAULT_SITE_SYMMETRY);
  PhysRowSetHeight(DEFAULT_SITE_HEIGHT);
  PhysRowSetSiteSpacing(DEFAULT_SITE_SPACING);
  PhysRowSetSiteWidth(DEFAULT_SITE_WIDTH);
  PhysRowSetNumSites(DEFAULT_NUM_SITES);
  PhysRowSetNumSubRows(DEFAULT_NUM_SUBROWS);
  PhysRowSetSubRows(DEFAULT_SUBROWS);
}

PhysRow::PhysRow(rowOrientation orient, int coordinate, unsigned int height,
		 unsigned int siteWidth, unsigned int siteSpacing, 
		 siteOrientation siteOrient, siteSymmetry symmetry)
{
  PhysRowSetType(orient);
  PhysRowSetCoordinate(coordinate);
  PhysRowSetSiteOrientation(siteOrient);
  PhysRowSetSiteSymmetry(symmetry);
  PhysRowSetHeight(height);
  PhysRowSetSiteSpacing(siteSpacing);
  PhysRowSetSiteWidth(siteWidth);
  PhysRowSetNumSites(DEFAULT_NUM_SITES);
  PhysRowSetNumSubRows(DEFAULT_NUM_SUBROWS);
  PhysRowSetSubRows(DEFAULT_SUBROWS);
}

PhysRow::PhysRow(rowOrientation orient, int coordinate, unsigned int height,
		 unsigned int siteWidth, unsigned int siteSpacing, 
		 unsigned int numSubRows, map<unsigned int, unsigned int> subRows)
{
  PhysRowSetType(orient);
  PhysRowSetCoordinate(coordinate);
  PhysRowSetSiteOrientation(DEFAULT_SITE_ORIENTATION);
  PhysRowSetSiteSymmetry(DEFAULT_SITE_SYMMETRY);
  PhysRowSetHeight(height);
  PhysRowSetSiteSpacing(siteSpacing);
  PhysRowSetSiteWidth(siteWidth);
  PhysRowSetNumSites(DEFAULT_NUM_SITES);
  PhysRowSetNumSubRows(numSubRows);
  PhysRowSetSubRows(subRows);  
}

PhysRow::PhysRow(rowOrientation orient, int coordinate, unsigned int height,
		 map<unsigned int, unsigned int> subRows)
{
  PhysRowSetType(orient);
  PhysRowSetCoordinate(coordinate);
  PhysRowSetSiteOrientation(DEFAULT_SITE_ORIENTATION);
  PhysRowSetSiteSymmetry(DEFAULT_SITE_SYMMETRY);
  PhysRowSetHeight(height);
  PhysRowSetSiteSpacing(DEFAULT_SITE_SPACING);
  PhysRowSetSiteWidth(DEFAULT_SITE_WIDTH);
  PhysRowSetNumSites(DEFAULT_NUM_SITES);
  PhysRowSetNumSubRows(DEFAULT_NUM_SUBROWS);
  PhysRowSetSubRows(subRows);
}

rowOrientation 
PhysRowGetRowTypeFromStr(string rowType)
{
  rowOrientation retVal;

  if (rowType == "Horizontal") {
    retVal = HORIZONTAL;
  } else if (rowType == "Vertical") {
    retVal = VERTICAL;
  }
  
  return (retVal);
}	    

siteOrientation 
PhysRowGetSiteOrientationFromStr(string siteOrient)
{
  siteOrientation retVal;
  
  if (siteOrient == "N" || siteOrient == "1") retVal = N;
  else if (siteOrient == "E") retVal = E;
  else if (siteOrient == "S") retVal = S;
  else if (siteOrient == "W") retVal = W;
  else if (siteOrient == "FN") retVal = FN;
  else if (siteOrient == "FE") retVal = FE;
  else if (siteOrient == "FS") retVal = FS;
  else if (siteOrient == "FW") retVal = FW;
  
  return (retVal);
}

siteSymmetry 
PhysRowGetSiteSymmetryFromStr(string symmetry)
{
  siteSymmetry retVal;
  
  if (symmetry == "1" || "Y") retVal = YES_SYMMETRY;
  else retVal = NO_SYMMETRY;
  
  return retVal;
}
