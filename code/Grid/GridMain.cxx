# include <Grid.h>

Grid::Grid()
{
     GridSetgridX(0);
     GridSetgridY(0);
     GridSetgridPotential(0);
        
}

void
Grid::GridSetName(const string& Name)
{
      gridName = Name;
}


void
Grid::GridSetgridX(double X)
{
      this->gridX = X;
}

void
Grid::GridSetgridY(double Y)
{
      this->gridY = Y ; 
}

void
Grid::GridSetgridPotential(double pot)
{
     this->gridPotential = pot;
}


string  
Grid::GridGetName(void)
{
     return gridName;
}

double 
Grid::GridGetgridX(void)
{
   return this->gridX;
}

double
Grid::GridGetgridY(void)
{
    return this->gridY;
}

double
Grid::GridGetgridPotential(void)
{
                return this->gridPotential;
}

Grid::~Grid()
{}













