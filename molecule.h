//
// Created by lurker on 2/10/17.
//

#ifndef MOLECULE_H
#define MOLECULE_H

#include "ls_point.h"

class Molecule {
public:
    vector<ls_point> centers;
    vector<scalar_t > radii;
    vector<scalar_t > charges;
    vector<vector<index_t >> adjacent;

    ls_point center;
    scalar_t radius;
    index_t N;

    void load(std::string fileName) {
        std::ifstream pdbFile;
        std::string line;

        pdbFile.open(fileName);

        if (pdbFile.is_open()) {
            // read each line
            while (getline(pdbFile, line)) {
                if (line.find("ATOM") == 0) {
                    std::stringstream ss(line);
                    std::string buf;
                    vector<std::string> tokens;
                    while (ss >> buf) {
                        tokens.push_back(buf);
                    }

                    centers.push_back(ls_point(
                            std::stod(tokens[5]),
                            std::stod(tokens[6]),
                            std::stod(tokens[7])
                    ));

                    radii.push_back(std::stof(tokens[9]));
                    charges.push_back(std::stof(tokens[8]));
                }
            }
            pdbFile.close();
            N = (index_t) centers.size();
            getCenter();
        }
    }
    void getAdjacent(scalar_t probe) {
        adjacent.resize(N);
        for (int i = 0; i < N; ++i) {
            scalar_t ri = radii[i] + probe;
            for (int j = 0; j < i; ++j) {
                // check if intersected.
                scalar_t  d = norm(centers[i] - centers[j]);
                scalar_t  rj = radii[j] + probe;
                if (d < ri + rj) {
                    adjacent[i].push_back(j);
                    adjacent[j].push_back(i);

                }
            }
        }
    }
    void getCenter() {
        assert(N > 0);
        ls_point minP = {
                centers[0].data[0] - radii[0],
                centers[0].data[1] - radii[0],
                centers[0].data[2] - radii[0]};
        ls_point maxP = {
                centers[0].data[0] + radii[0],
                centers[0].data[1] + radii[0],
                centers[0].data[2] + radii[0]};

        for (int i = 1; i < N; ++i) {
            ls_point current_max_P = {
                    centers[i].data[0] + radii[i],
                    centers[i].data[1] + radii[i],
                    centers[i].data[2] + radii[i]
            };

            ls_point current_min_P = {
                    centers[i].data[0] - radii[i],
                    centers[i].data[1] - radii[i],
                    centers[i].data[2] - radii[i]
            };


            for (int j = 0; j < 3; ++j) {
                if (maxP.data[j] <= current_max_P.data[j]) {maxP.data[j] = current_max_P.data[j];}

                if (minP.data[j] >= current_min_P.data[j]) {minP.data[j] = current_min_P.data[j];}
            }

        }

        this->center = (minP + maxP) * 0.5;
        this->radius = std::max(
                0.5 * (maxP.data[0] - minP.data[0]),
                std::max(
                        0.5 * (maxP.data[1] - minP.data[1]),
                        0.5 * (maxP.data[2] - minP.data[2]))
        );
    }


    scalar_t centralize(float range) {
        scalar_t s = range / radius;

        for (index_t i = 0; i < (index_t)centers.size(); ++i) {
            centers[i].data[0] = s * (centers[i].data[0] - center.data[0]);
            centers[i].data[1] = s * (centers[i].data[1] - center.data[1]);
            centers[i].data[2] = s * (centers[i].data[2] - center.data[2]);
            radii[i] = s * radii[i];
        }

        std::cout << std::setw(15)<< "ATOMS" << " " << std::setw(8) << N  <<std::endl;


        return s;
    }
};


#endif //MOLECULE_H
