# include <Design.h>


void 
Design::DesignCreateGridPoints(uint numGridPoints)
{
        bool debug = true;
        uint maxx,maxy;
        DesignGetBoundingBox(maxx,maxy);
        uint intervalX,intervalY;
        intervalX = maxx/numGridPoints;
        intervalY = maxy/numGridPoints;
        char *gridName = "grid";
        int suffix = 1;
       // char buf[numGridPoints*numGridPoints];

        /* For loop to create grid points*/
        for (int i = 1;i <= numGridPoints;i++)
        {
                for(int j = 1;j <= numGridPoints;j++)
                {
               /* Create Grid object and initialize the X and Y values based on the interval and potential to 0*/ 
                Grid *GridObj = new Grid;
                char buf[BUFSIZ];
                sprintf (buf,"%s%d",gridName,suffix);
                (*GridObj).GridSetName(buf);
                (*GridObj).GridSetgridX(i*intervalX);
                (*GridObj).GridSetgridY(j*intervalY);
                suffix +=1;
                DesignAddOneGridToDesignDB(GridObj);
                }

        }

}


void 
computeGridPotential (double &distance,uint &radius,double &computedPotential)
{

        if (distance >= radius){
                computedPotential = 0;
        } else if ((distance >= 0) && (distance < (radius/2))){
                computedPotential = (1-2*pow((distance/radius),2));
        }else if ((distance >= (radius/2)) && (distance < radius)){
                computedPotential = 2*pow(((distance-radius)/radius),2);
        }

}



void
Design::DesignUpdateGridPotentials()
{

double cellXpos,cellYpos;
int cellHeight,cellWidth;
double cellXright,cellYtop;
double cellXcenter,cellYcenter;
uint radius;
double distanceX,distanceY;
string cellName;
Cell *cellPtr;

DESIGN_FOR_ALL_CELLS((*this),cellName,cellPtr){

        cellXpos = (*cellPtr).CellGetXposDbl();
        cellYpos = (*cellPtr).CellGetYposDbl();
        cellHeight = (*cellPtr).CellGetHeight();
        cellWidth = (*cellPtr).CellGetWidth();
        cellXcenter = (cellXpos+cellWidth)/2;
        cellYcenter = (cellYpos+cellHeight)/2;
        cellXright = cellXpos + cellWidth;
        cellYtop = cellYpos + cellHeight;
        radius = (cellHeight+cellWidth)/2;
        string gridName;
        Grid *gridPtr;
        double gridXpos,gridYpos;
        double distanceX,distanceY;
        double gridPotentialX,gridPotentialY;
        double gridPotential;
        double currentPotential;
        MAP_FOR_ALL_ELEMS(DesignGridPoints,string, Grid *,gridName,gridPtr){
                gridXpos = (*gridPtr).GridGetgridX();
                gridYpos = (*gridPtr).GridGetgridY();
                if((gridXpos>=cellXpos) && (gridXpos <= cellXright) && (gridYpos>=cellYpos) && (gridYpos <= cellYtop)){
                        distanceX= fabs (gridXpos - cellXcenter); 
                        distanceY= fabs (gridYpos - cellYcenter);
                        computeGridPotential  (distanceX,radius, gridPotentialX); 
                        computeGridPotential  (distanceY,radius, gridPotentialY);
                        gridPotential = gridPotentialX*gridPotentialY; 
                        currentPotential = (*gridPtr).GridGetgridPotential();
                        currentPotential = currentPotential + gridPotential; 
                        (*gridPtr).GridSetgridPotential(currentPotential);
                }

        }END_FOR;

}DESIGN_END_FOR;

}


double 
Design::DesignComputeTotalDensityPenalty(void)
{
       double rtv;
       double averagePotential;
       uint maxx,maxy;
       DesignGetBoundingBox(maxx,maxy);
       double numGridPoints = 100;
       averagePotential = (maxx*maxy)/numGridPoints;
       double totalPenalty;
       string gridName;
       Grid *gridPtr;
       double gridPotential;
       double gridPointPenalty;
       MAP_FOR_ALL_ELEMS(DesignGridPoints,string,Grid *,gridName,gridPtr){
                 gridPotential=(*gridPtr).GridGetgridPotential();
                 gridPointPenalty = pow((gridPotential - averagePotential),2);
                 totalPenalty = totalPenalty + gridPointPenalty ;
       }END_FOR;

        rtv = totalPenalty;
        return rtv;
}
