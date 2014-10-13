# ifndef GRID_H
# define GRID_H

# include <common.h>
# include <Flags.h>

class Grid {
   private:
        double gridX;
        double gridY;
        double gridPotential;

   public: 
        string gridName;
        
        
   /* Constructor*/     
        Grid();
   
   /* Set Functions*/
        void GridSetName(const string&);
        void GridSetgridX(double);
        void GridSetgridY(double);
        void GridSetgridPotential(double);

   /* Get Functions*/
        string GridGetName(void);
        double GridGetgridX(void);
        double GridGetgridY(void);
        double GridGetgridPotential(void);

    /* Destructor*/
    ~Grid();
};

# endif
