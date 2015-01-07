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
computeGridPotentialGradient (double &distance,uint &radius,double &computedPotential)
{

        if (distance >= radius){
                computedPotential = 0;
        } else if ((distance >= 0) && (distance < (radius/2))){
                computedPotential = 4*(distance)/(radius*radius);
        }else if ((distance >= (radius/2)) && (distance < radius)){
                computedPotential = 4*((radius-distance)/(radius*radius));
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
double normalizationFactor; 
uint maxx,maxy;
DesignGetBoundingBox(maxx,maxy);

// the below section of code calculates the number of gridPoints influenced by a cell 
Env &DesignEnv = this->DesignEnv;
uint numGridPoints = DesignEnv.EnvGetNumGridPoints();
double numGridPointsInfluenced;
double gridPointXInterval;
double gridPointYInterval;
gridPointXInterval = maxx/numGridPoints;
gridPointYInterval = maxy/numGridPoints;
// Iitialize all cell potentials to zero before starting to update the grid potential

Grid *gridPtrInt;
string gridNameInt;

MAP_FOR_ALL_ELEMS(DesignGridPoints,string, Grid *,gridNameInt,gridPtrInt){
                (*gridPtrInt).GridSetgridPotential(0);
}END_FOR;


DESIGN_FOR_ALL_CELLS((*this),cellName,cellPtr){
        if ((*cellPtr).CellIsTerminal()) continue;
        cellXpos = (*cellPtr).CellGetXposDbl();
        cellYpos = (*cellPtr).CellGetYposDbl();
        cellHeight = (*cellPtr).CellGetHeight();
        cellWidth = (*cellPtr).CellGetWidth();
        cellXcenter = cellXpos+(cellWidth/2);
        cellYcenter = cellYpos+(cellHeight/2);
        cellXright = cellXpos + cellWidth;
        cellYtop = cellYpos + cellHeight;
        radius = (cellHeight+cellWidth)/4;
        //normalization factor will ideally be area/numGridPoints influenced by cell 
        numGridPointsInfluenced = (cellWidth/gridPointXInterval) * (cellHeight/gridPointYInterval);
        normalizationFactor= (cellHeight*cellWidth)/numGridPointsInfluenced;
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
                        gridPotential = normalizationFactor*gridPotentialX*gridPotentialY; 
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
       Env &DesignEnv = this->DesignEnv;
       uint numGridPoints = DesignEnv.EnvGetNumGridPoints();
       double averagePotential;
       uint maxx,maxy;
       DesignGetBoundingBox(maxx,maxy);
       //numGridPoints in one axis. total nuber of grid points will be equal to numGridPoints*numGridPoints  
      // double numGridPoints = 100;
       uint averageClusterWidth,averageClusterHeight;
       averageClusterWidth = (uint)DesignGetAverageClusterCellWidth();
       averageClusterHeight = (uint)DesignGetAverageClusterCellHeight()+2*averageStdCellHeight;
       //averagePotential = (averageClusterWidth*averageClusterHeight)/(numGridPoints*numGridPoints);
       // Above average potential was just a try . The below average potential is after considering the normalization parameter
       averagePotential = (maxx*maxy)/(numGridPoints*numGridPoints);
       double totalPenalty=0;
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

double 
Design::DesignComputeTotalDensityPenaltyGradient(void){

       double rtv;
       Env &DesignEnv = this->DesignEnv;
       uint numGridPoints = DesignEnv.EnvGetNumGridPoints();
       double averagePotential;
       uint maxx,maxy;
       DesignGetBoundingBox(maxx,maxy);
       //numGridPoints in one axis. total nuber of grid points will be equal to numGridPoints*numGridPoints  
      // double numGridPoints = 100;
       uint averageClusterWidth,averageClusterHeight;
       averageClusterWidth = (uint)DesignGetAverageClusterCellWidth();
       averageClusterHeight = (uint)DesignGetAverageClusterCellHeight()+2*averageStdCellHeight;
       //averagePotential = (averageClusterWidth*averageClusterHeight)/(numGridPoints*numGridPoints);
       // Above average potential was just a try . The below average potential is after considering the normalization parameter
       averagePotential = (maxx*maxy)/(numGridPoints*numGridPoints);
       double totalPenalty=0;
       string gridName;
       Grid *gridPtr;
       double gridPotential;
       double gridPointPenalty;
       MAP_FOR_ALL_ELEMS(DesignGridPoints,string,Grid *,gridName,gridPtr){
                 gridPotential=(*gridPtr).GridGetgridPotential();
                 gridPointPenalty = (gridPotential - averagePotential);
                 totalPenalty = totalPenalty + gridPointPenalty ;
       }END_FOR;

        rtv = totalPenalty;
        return rtv;

}

// Below functions are used to update the gradient of the peanlty function for each iteration of the solver.
// The object for the cell need to be passed. This object has the new x and y positions for the cells
// Based on the influence this cell has over grid points, gradient is computed
void
Design::DesignComputePenaltyGradientforCell(Cell *cellPtr,double &gradPenaltyX,double &gradPenaltyY){

        double cellXpos,cellYpos;
        int cellHeight,cellWidth;
        double cellXright,cellYtop;
        double cellXcenter,cellYcenter;
        uint radius;
        double normalizationFactor; 
        uint maxx,maxy;
        DesignGetBoundingBox(maxx,maxy);


        Env &DesignEnv = this->DesignEnv;
        uint numGridPoints = DesignEnv.EnvGetNumGridPoints();
        double numGridPointsInfluenced;
        double gridPointXInterval;
        double gridPointYInterval;
        gridPointXInterval = maxx/numGridPoints;
        gridPointYInterval = maxy/numGridPoints;
        
        cellXpos = (*cellPtr).CellGetXposDbl();
        cellYpos = (*cellPtr).CellGetYposDbl();
        cellHeight = (*cellPtr).CellGetHeight();
        cellWidth = (*cellPtr).CellGetWidth();
        cellXcenter = cellXpos+(cellWidth/2);
        cellYcenter = cellYpos+(cellHeight/2);
        cellXright = cellXpos + cellWidth;
        cellYtop = cellYpos + cellHeight;
        radius = (cellHeight+cellWidth)/4;
        numGridPointsInfluenced = (cellWidth/gridPointXInterval) * (cellHeight/gridPointYInterval);
        normalizationFactor= (cellHeight*cellWidth)/numGridPointsInfluenced;
                
       double averagePotential;
       averagePotential = (maxx*maxy)/(numGridPoints*numGridPoints);
       
        string gridName;
        Grid *gridPtr;
        double gridXpos,gridYpos;
        double distanceX,distanceY;
        double gridPotentialX,gridPotentialY;
        double gridPotential;
        double currentPotential;
        gradPenaltyX =0;
        gradPenaltyY=0;

        MAP_FOR_ALL_ELEMS(DesignGridPoints,string, Grid *,gridName,gridPtr){
                gridXpos = (*gridPtr).GridGetgridX();
                gridYpos = (*gridPtr).GridGetgridY();
                if((gridXpos>=cellXpos) && (gridXpos <= cellXright) && (gridYpos>=cellYpos) && (gridYpos <= cellYtop)){
                        distanceX= fabs (gridXpos - cellXcenter); 
                        distanceY= fabs (gridYpos - cellYcenter);
                        computeGridPotentialGradient  (distanceX,radius, gridPotentialX); 
                        computeGridPotentialGradient  (distanceY,radius, gridPotentialY);
                        gradPenaltyX = gradPenaltyX + normalizationFactor*gridPotentialX; 
                        gradPenaltyY = gradPenaltyY + normalizationFactor*gridPotentialY; 
                                
                }
        }END_FOR;

}


void
Design::DesignComputePenaltyGradientforCellX(Cell *cellPtr,double &gradPenaltyX){

        double cellXpos,cellYpos;
        int cellHeight,cellWidth;
        double cellXright,cellYtop;
        double cellXcenter,cellYcenter;
        uint radius;
        double normalizationFactor; 
        uint maxx,maxy;
        DesignGetBoundingBox(maxx,maxy);


        Env &DesignEnv = this->DesignEnv;
        uint numGridPoints = DesignEnv.EnvGetNumGridPoints();
        double numGridPointsInfluenced;
        double gridPointXInterval;
        double gridPointYInterval;
        gridPointXInterval = maxx/numGridPoints;
        gridPointYInterval = maxy/numGridPoints;
        
        cellXpos = (*cellPtr).CellGetXposDbl();
        cellYpos = (*cellPtr).CellGetYposDbl();
        cellHeight = (*cellPtr).CellGetHeight();
        cellWidth = (*cellPtr).CellGetWidth();
        cellXcenter = cellXpos+(cellWidth/2);
        cellYcenter = cellYpos+(cellHeight/2);
        cellXright = cellXpos + cellWidth;
        cellYtop = cellYpos + cellHeight;
        radius = (cellHeight+cellWidth)/4;
        numGridPointsInfluenced = (cellWidth/gridPointXInterval) * (cellHeight/gridPointYInterval);
        normalizationFactor= (cellHeight*cellWidth)/numGridPointsInfluenced;
                
       double averagePotential;
       averagePotential = (maxx*maxy)/(numGridPoints*numGridPoints);
       
        string gridName;
        Grid *gridPtr;
        double gridXpos,gridYpos;
        double distanceX,distanceY;
        double gridPotentialX,gridPotentialY;
        double gridPotential;
        double currentPotential;
        gradPenaltyX =0;
        //gradPenaltyY=0;

        MAP_FOR_ALL_ELEMS(DesignGridPoints,string, Grid *,gridName,gridPtr){
                gridXpos = (*gridPtr).GridGetgridX();
                gridYpos = (*gridPtr).GridGetgridY();
                if((gridXpos>=cellXpos) && (gridXpos <= cellXright) && (gridYpos>=cellYpos) && (gridYpos <= cellYtop)){
                        distanceX= fabs (gridXpos - cellXcenter); 
                        distanceY= fabs (gridYpos - cellYcenter);
                        computeGridPotentialGradient  (distanceX,radius, gridPotentialX); 
                        //computeGridPotentialGradient  (distanceY,radius, gridPotentialY);
                        gradPenaltyX = gradPenaltyX + normalizationFactor*gridPotentialX; 
                        //gradPenaltyY = gradPenaltyY + normalizationFactor*gridPotentialY; 
                                
                }
        }END_FOR;

}


void
Design::DesignComputePenaltyGradientforCellY(Cell *cellPtr,double &gradPenaltyY){

        double cellXpos,cellYpos;
        int cellHeight,cellWidth;
        double cellXright,cellYtop;
        double cellXcenter,cellYcenter;
        uint radius;
        double normalizationFactor; 
        uint maxx,maxy;
        DesignGetBoundingBox(maxx,maxy);


        Env &DesignEnv = this->DesignEnv;
        uint numGridPoints = DesignEnv.EnvGetNumGridPoints();
        double numGridPointsInfluenced;
        double gridPointXInterval;
        double gridPointYInterval;
        gridPointXInterval = maxx/numGridPoints;
        gridPointYInterval = maxy/numGridPoints;
        
        cellXpos = (*cellPtr).CellGetXposDbl();
        cellYpos = (*cellPtr).CellGetYposDbl();
        cellHeight = (*cellPtr).CellGetHeight();
        cellWidth = (*cellPtr).CellGetWidth();
        cellXcenter = cellXpos+(cellWidth/2);
        cellYcenter = cellYpos+(cellHeight/2);
        cellXright = cellXpos + cellWidth;
        cellYtop = cellYpos + cellHeight;
        radius = (cellHeight+cellWidth)/4;
        numGridPointsInfluenced = (cellWidth/gridPointXInterval) * (cellHeight/gridPointYInterval);
        normalizationFactor= (cellHeight*cellWidth)/numGridPointsInfluenced;
                
       double averagePotential;
       averagePotential = (maxx*maxy)/(numGridPoints*numGridPoints);
       
        string gridName;
        Grid *gridPtr;
        double gridXpos,gridYpos;
        double distanceX,distanceY;
        double gridPotentialX,gridPotentialY;
        double gridPotential;
        double currentPotential;
        //gradPenaltyX =0;
        gradPenaltyY=0;

        MAP_FOR_ALL_ELEMS(DesignGridPoints,string, Grid *,gridName,gridPtr){
                gridXpos = (*gridPtr).GridGetgridX();
                gridYpos = (*gridPtr).GridGetgridY();
                if((gridXpos>=cellXpos) && (gridXpos <= cellXright) && (gridYpos>=cellYpos) && (gridYpos <= cellYtop)){
                        distanceX= fabs (gridXpos - cellXcenter); 
                        distanceY= fabs (gridYpos - cellYcenter);
                        computeGridPotentialGradient  (distanceX,radius, gridPotentialX); 
                        computeGridPotentialGradient  (distanceY,radius, gridPotentialY);
               //         gradPenaltyX = gradPenaltyX + normalizationFactor*gridPotentialX; 
                        gradPenaltyY = gradPenaltyY + normalizationFactor*gridPotentialY; 
                                
                }
        }END_FOR;

}

void
Design::DesignComputepenaltyParameter(void){
        static int gradIterationCount=0;
        Cell *cellPtr;
        string cellName;
        uint alpha = 500;
        uint idx=0;
        double cellXpos;  
        double cellYpos;
        double WLgradientTotal;
        double penaltyGradientTotal;
        double penaltyParameter;
        DESIGN_FOR_ALL_CELLS((*this),cellName,cellPtr){
                if ((*cellPtr).CellIsTerminal()) continue;
                double gradX=0;
                double gradY=0;
                double cellMaxx;
                double cellMinx;
                double cellMaxy;
                double cellMiny;
                cellXpos = (*cellPtr).CellGetXposDbl();
                cellYpos = (*cellPtr).CellGetYposDbl();
                Net *netPtr;
                double sumPinsPos;
                double pinMaxx;
                double pinMinx;
                double pinMaxy;
                double pinMiny;
                double tempDivideX;
                double tempDivideY;
                double pinXpos;
                double pinYpos;
                CELL_FOR_ALL_NETS_NO_DIR((*cellPtr),netPtr){
                      Pin *pinPtr;
                      NET_FOR_ALL_PINS((*netPtr),pinPtr){
                             Cell* cellParentPtr;
                             cellParentPtr = (*pinPtr).PinGetParentCellPtr();
                             pinXpos = pinPtr->xOffset + cellParentPtr->x;
                             pinYpos = pinPtr->yOffset + cellParentPtr->y;
                             tempDivideX= myDivideNew(pinXpos,alpha);
                             tempDivideY= myDivideNew(pinYpos,alpha);
                             pinMaxx += exp(tempDivideX);
                             pinMinx +=myDivideNew(1,exp(tempDivideX));
                             pinMaxy += exp(tempDivideY);
                             pinMiny +=myDivideNew(1,exp(tempDivideY));
                      }NET_END_FOR;
                }CELL_END_FOR;
                double tempDivX;
                double tempDivY;
                tempDivX = myDivideNew(cellXpos,alpha);
                tempDivY = myDivideNew(cellYpos,alpha);
                cellMaxx = exp(tempDivX);
                cellMinx = myDivideNew(1,exp(tempDivX));
                cellMaxy = exp(tempDivY);
                cellMiny = myDivideNew(1,exp(tempDivY));
                double temp1;
                double temp2;
                temp1 = myDivideNew(cellMaxx,pinMaxx);
                temp2 = myDivideNew(cellMinx,pinMinx);
                gradX = (temp1)-(temp2);
                double temp3;
                double temp4;
                temp3 = myDivideNew(cellMaxy,pinMaxy);
                temp4 = myDivideNew(cellMiny,pinMiny);
                gradY = (temp3)-(temp4);
                // Change to add the gradient for the penalty is here
                double gradPotentialX;
                double gradPotentialY;
                DesignComputePenaltyGradientforCell(cellPtr,gradPotentialX,gradPotentialY);
                WLgradientTotal = WLgradientTotal + gradX + gradY;
                penaltyGradientTotal = penaltyGradientTotal + gradPotentialX + gradPotentialY;
                
  }DESIGN_END_FOR;
  double totalDensityPenaltyGradient;
  totalDensityPenaltyGradient = DesignComputeTotalDensityPenaltyGradient();
  penaltyParameter= (WLgradientTotal/penaltyGradientTotal);
  /*if (penaltyParameter < 0 ) {
          penaltyParameter = -penaltyParameter;
  }*/
  cout << "WLgradientTotal is: " << WLgradientTotal << " penaltyGradientTotal is: " << penaltyGradientTotal ; 
  cout << " penaltyParameterInitial Value is: " << penaltyParameter << endl;
  DesignSetpenaltyParameter(penaltyParameter);
  gradIterationCount++;

}
