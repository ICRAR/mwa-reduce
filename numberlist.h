#ifndef NUMBERLIST_H
#define NUMBERLIST_H

#include <cstdlib>
#include <string>

#include <aocommon/uvector.h>

class NumberList {
 public:
  static aocommon::UVector<int> ParseIntList(const std::string& str) {
    aocommon::UVector<int> list;
    std::string temp = str;
    size_t pos = temp.find(",");
    while (pos != std::string::npos) {
      std::string idStr = temp.substr(0, pos);
      temp = temp.substr(pos + 1);
      int num = atoi(idStr.c_str());
      list.push_back(num);
      pos = temp.find(",");
    }
    int num = atoi(temp.c_str());
    list.push_back(num);
    return list;
  }

  static aocommon::UVector<double> ParseDoubleList(const std::string& str) {
    aocommon::UVector<double> list;
    std::string temp = str;
    size_t pos = temp.find(",");
    while (pos != std::string::npos) {
      std::string idStr = temp.substr(0, pos);
      temp = temp.substr(pos + 1);
      double num = atof(idStr.c_str());
      list.push_back(num);
      pos = temp.find(",");
    }
    double num = atof(temp.c_str());
    list.push_back(num);
    return list;
  }
};

#endif
