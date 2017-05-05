//
// Created by lurker on 5/3/17.
//

#include "Config.h"

Config::Config() {
    options.clear();
}

Config::~Config() {
    options.clear();
}

void Config::parse(std::istream &cfgFile) {
    for (std::string line; std::getline(cfgFile, line); ) {
        std::istringstream iss(line);
        std::string id, eq, val;


        bool error = false;

        if (!(iss >> id)) {error = true;}
        else if(id[0] =='#') { continue;}
        else if (!(iss >> eq >> val >> std::ws) || eq != "=" || iss.get() != EOF) {
            error = true;
        }

        if (line.size() == 0) error = false;

        if (error) {
            std::cout << "Incomplete configuration! Did you forget the space?" << std::endl;
        }
        else {
            /*
             * eliminate empty string key
             */
            if (line.size()) options[id] = val;
        }

    }
}

void Config::print(){
    std::cout << std::setw(15) << "=============== OPTIONS ===============" << std::endl;
    for (auto key:options) {
        std::cout <<std::setw(15) <<  key.first << " = " << key.second << std::endl;
    }
    std::cout << std::setw(15) << "=======================================" << std::endl;
}

