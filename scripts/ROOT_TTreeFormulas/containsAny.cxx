#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>

bool containsAny( const std::vector<std::string>& inList = std::vector<string>(), const std::string& checkStr = "" )
{

  std::vector<std::string> checkList;

  std::string token;
  std::istringstream ss(checkStr);
  while ( std::getline(ss, token, ',') ) {
    checkList.push_back(token);
  }

  for ( auto itr : checkList ) {
    if ( std::find(inList.begin(), inList.end(), itr) != inList.end() ) { return true; }
  }

  return false;
}
