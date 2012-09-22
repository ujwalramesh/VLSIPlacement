# include <Pin.h>
# include <common.h>

Pin::Pin() 
{
  Id = 0;
  Name = NIL(string);
  xOffset = 0;
  yOffset = 0;
  dir = PIN_DIR_INPUT;
}
Pin::Pin(int id) 
{
  Id = id;
  ParentCell = NULL;
  ConnectedNet = NULL;
  xOffset = 0;
  yOffset = 0;
  dir = PIN_DIR_INPUT;
}

Pin::Pin(int id, int xoffset, int yoffset) 
{
  Id = id;
  xOffset = xoffset;
  yOffset = yoffset;
  ParentCell = NULL;
  ConnectedNet = NULL;
  dir = PIN_DIR_INPUT;
}

Pin::Pin(int id, int xoffset, int yoffset, char direction) 
{
  Id = id;
  xOffset = xoffset;
  yOffset = yoffset;
  ParentCell = NULL;
  ConnectedNet = NULL;
}

Pin::Pin(int id, const string& name) 
{
  Id = id;
  Name = name;
  xOffset = 0;
  yOffset = 0;
  dir = PIN_DIR_INPUT;
  ConnectedNet = NULL;
  ParentCell = NULL;
}

Pin::Pin(int id, int xoffset, int yoffset, const string& name) 
{
  Id = id;
  Name = name;
  ConnectedNet = NULL;
  ParentCell = NULL;
  xOffset = xoffset;
  yOffset = yoffset;
  dir = PIN_DIR_INPUT;
}

Pin::Pin(int id, int xoffset, int yoffset, char direction, const string& name)
{
  Id = id;
  Name = name;
  ConnectedNet = NULL;
  ParentCell = NULL;
  xOffset = xoffset;
  yOffset = yoffset;
  dir = direction;
}

Pin::Pin(int id, const Cell& parentCell) 
{
  Id = id;
  ParentCell = (Cell *) &parentCell;
  ConnectedNet = NULL;
  xOffset = 0;
  yOffset = 0;
  dir = PIN_DIR_INPUT;
}

Pin::Pin(int id, int xoffset, int yoffset, const Cell& parentCell) 
{
  Id = id;
  ParentCell = (Cell *) &parentCell;
  ConnectedNet = NULL;
  xOffset = xoffset;
  yOffset = yoffset;
  dir = PIN_DIR_INPUT;
}

Pin::Pin(int id, int xoffset, int yoffset, char direction, const Cell& parentCell) 
{
  Id = id;
  ParentCell = (Cell *) &parentCell;
  ConnectedNet = NULL;
  xOffset = xoffset;
  yOffset = yoffset;
  dir = direction;
}

Pin::Pin(int id, const Cell& parentCell, const string& name) 
{
  Id = id;
  ParentCell = (Cell *) &parentCell;
  Name = name;
  ConnectedNet = NULL;
  xOffset = 0;
  yOffset = 0;
  dir = PIN_DIR_INPUT;
}

Pin::Pin(int id, int xoffset, int yoffset, const Cell& parentCell, 
	 const string& name) 
{
  Id = id;
  ParentCell = (Cell *) &parentCell;
  Name = name;
  ConnectedNet = NULL;
  xOffset = xoffset;
  yOffset = yoffset;
  dir = PIN_DIR_INPUT;
}

Pin::Pin(int id, int xoffset, int yoffset, char direction, 
	 const Cell& parentCell, const string& name) 
{
  Id = id;
  ParentCell = (Cell *) &parentCell;
  Name = name;
  ConnectedNet = NULL;
  xOffset = xoffset;
  yOffset = yoffset;
  dir = direction;
}

Pin::Pin(int id, const Cell& parentCell, const Net& connectedNet)
{
  Id = id;
  ParentCell = (Cell *) &parentCell;
  ConnectedNet = (Net *) &connectedNet;
  xOffset = 0;
  yOffset = 0;
  dir = PIN_DIR_INPUT;
}

Pin::Pin(int id, int xoffset, int yoffset, const Cell& parentCell, 
	 const Net& connectedNet)
{
  Id = id;
  ParentCell = (Cell *) &parentCell;
  ConnectedNet = (Net *) &connectedNet;
  xOffset = xoffset;
  yOffset = yoffset;
  dir = PIN_DIR_INPUT;
}

Pin::Pin(int id, int xoffset, int yoffset, char direction, 
	 const Cell& parentCell, const Net& connectedNet)
{
  Id = id;
  ParentCell = (Cell *) &parentCell;
  ConnectedNet = (Net *) &connectedNet;
  xOffset = xoffset;
  yOffset = yoffset;
  dir = direction;
}

/* Set functions */
void
Pin::PinSetId(int id)
{
  Id = id;
}

void 
Pin::PinSetName(string pinName)
{
  Name = pinName;
}
  
void 
Pin::PinSetParentCell(const Cell& parentCell)
{
  ParentCell = (Cell *) &parentCell;
}

void
Pin::PinSetXOffset(int xoffset)
{
  xOffset = xoffset;
}

void 
Pin::PinSetYOffset(int yoffset)
{
  yOffset = yoffset;
}

void
Pin::PinSetDirection(char direction)
{
  dir = direction;
}

/* Get functions */
int 
Pin::PinGetId(void)
{
  return (Id);
}

string 
Pin::PinGetName(void) const
{
  return (Name);
}

int
Pin::PinGetXOffset(void)
{
  return (xOffset);
}

int 
Pin::PinGetYOffset(void)
{
  return (yOffset);
}

char
Pin::PinGetDirection(void)
{
  return (dir);
}

Cell& 
Pin::PinGetParentCell(void)
{
  Cell &parentCell = *ParentCell;
  
  return (parentCell);
}

Net& 
Pin::PinGetNet(void)
{
  Net &connectedNet = *ConnectedNet;
  
  return (connectedNet);
}

/* Other functions */
void 
Pin::Connect(const Net& netToConnect)
{
  ConnectedNet = (Net *) &netToConnect;
}

Net& 
Pin::Disconnect(void)
{
  ConnectedNet = NULL;
}

void 
Pin::PinGetXposYpos(int *xpos, int *ypos)
{
  *xpos = xOffset;
  *ypos = yOffset;
}

