# include <common.h>
# include <stdlib.h>

void 
common_error(std::string ErrorString)
{
# ifdef COMMON_DEBUG  
  std::cout << std::endl << "Error: " << ErrorString << std::endl << "Quitting.." << std::endl;
  exit(0);
# endif
}

void 
common_message(std::string MessageString)
{
# ifdef COMMON_DEBUG
  //  std::cout << std::endl << "Message from program:" << MessageString << std::endl;
# endif 
}


string 
getStrFromInt(int thisInt)
{
  stringstream oss;
  
  oss << thisInt;
  
  return oss.str();
}

