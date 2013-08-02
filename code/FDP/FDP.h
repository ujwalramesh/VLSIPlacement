# ifndef FDPLACE_H
# define FDPLACE_H

# include <common.h>
# include <PhysRow.h>
# include <Cell.h>
# include <Pin.h>
# include <Design.h>
# include <Flags.h>
# include <HyperGraph.h>
# include <HyperGraphIter.h>

using namespace std;

/* Default definitions */

# define debugPrint 0

class FDPSite{
 private:
  Cell *siteCell;
  bool siteIsLocked;
  int siteX;
  int siteY;
  int siteNum;
  int rowNum;
  
 public:
  /* Constructor and Destructor */
  FDPSite();
  
  /* Set Functions */
  void FDPSiteSetCell(Cell*);
  void FDPSiteSetIsLocked(bool);
  void FDPSiteSetXpos(int);
  void FDPSiteSetYpos(int);
  void FDPSiteSetSiteNum(int);
  void FDPSiteSetRowNum(int);
  
  /* Get functions */
  Cell* FDPSiteGetCell(void);
  bool FDPSiteGetIsLocked(void);
  bool FDPSiteIsOccupied(void);
  bool FDPSiteIsVacant(void);
  int FDPSiteGetXpos(void);
  int FDPSiteGetYpos(void);
  int FDPSiteGetSiteNum(void);
  int FDPSiteGetRowNum(void);

  /* Other Function */
  void FDPSiteRemoveCell(void);
};

extern void FDPTopLevel(Design &, vector<Cell*> &, uint, uint, uint, uint,
			uint, uint);


# endif

