//
// Created by lurker on 5/3/17.
//

#ifndef LEVELSET_CONFIG_H
#define LEVELSET_CONFIG_H

#include "utils.h"
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>


class Config {
public:
    Config();
    ~Config();
    std::map<std::string, std::string> options;
    void parse(std::istream &cfgFile);
};


#endif //LEVELSET_CONFIG_H
