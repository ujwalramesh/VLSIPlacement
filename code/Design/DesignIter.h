# include <Design.h>

# ifndef DESIGN_ITER_H
# define DESIGN_ITER_H

/* Macros for iteration over elements of the design */

# define DESIGN_FOR_ALL_CELLS(Design, CellName, CellPtr) \
  { \
  map<string, Cell*> DesignCells = Design.DesignGetCells(); \
  map<string, Cell*>::iterator mapIter; \
  for (mapIter = DesignCells.begin(); mapIter != DesignCells.end(); mapIter++) { \
    CellName = mapIter->first; \
    CellPtr = mapIter->second;

# define DESIGN_FOR_ALL_NETS(Design, NetName, NetPtr) \
  { \
  map<string, Net*> DesignNets = Design.DesignGetNets(); \
  map<string, Net*>::iterator mapIter; \
  for (mapIter = DesignNets.begin(); mapIter != DesignNets.end(); mapIter++) { \
    NetName = (string)mapIter->first; \
    NetPtr = (Net *)mapIter->second;

# define DESIGN_END_FOR }}
  
# endif
