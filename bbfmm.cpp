//
// Created by lurker on 3/22/17.
//

#include "bbfmm.h"

namespace bbfmm {
    void
    tree::populate(vector<point> &_source, vector<point> &_target, index_t _nSource, index_t _nTarget, index_t _rank,
                   index_t _maxLevel) {
        this->sourceTree = _source;
        this->targetTree = _target;
        this->nSource = _nSource;
        this->nTarget = _nTarget;
        this->maxLevel = 0;
        this->rank = _rank;

        getCenterRadius(_source);
        this->root = 0;

        this->dict.push_back(node(0, 0));
        this->maxId = root;

        dict[root].nSource = nSource;
        dict[root].nTarget = nTarget;
        dict[root].center = center;
        dict[root].radius = radius;
        dict[root].sourceIndex.resize((unsigned long) nSource);
        dict[root].targetIndex.resize((unsigned long) nTarget);

        for (index_t i = 0; i < nSource; ++i) {
            dict[root].sourceIndex[i] = i;
        }
        for (index_t i = 0; i < nTarget; ++i) {
            dict[root].targetIndex[i] = i;
        }


        RUN("initialization", assignChildren(root, _maxLevel));
        RUN("assign lists", buildTree());
    }

    void tree::getCenterRadius(vector<point> &_source) {
        assert(_source.size() > 0);
        scalar_t x_max = _source[0].x;
        scalar_t x_min = _source[0].x;
        scalar_t y_max = _source[0].y;
        scalar_t y_min = _source[0].y;
        scalar_t z_max = _source[0].z;
        scalar_t z_min = _source[0].z;
        for (size_t i = 0; i < _source.size(); ++i) {
            x_max = std::max(x_max, _source[i].x);
            y_max = std::max(y_max, _source[i].y);
            z_max = std::max(z_max, _source[i].z);
            x_min = std::min(x_min, _source[i].x);
            y_min = std::min(y_min, _source[i].y);
            z_min = std::min(z_min, _source[i].z);
        }
        this->center.x = (x_max + x_min) / 2.0;
        this->center.y = (y_max + y_min) / 2.0;
        this->center.z = (z_max + z_min) / 2.0;
        this->radius.x = (x_max - x_min) / 2.0;
        this->radius.y = (y_max - y_min) / 2.0;
        this->radius.z = (z_max - z_min) / 2.0;
    }

    void tree::assignChildren(index_t _id, index_t _maxLevel) {
        /*
         * when assigning children nodes, the points are not assigned due to storage.
         *
         * Now the limitation of nodes is around 2^24.
         */
        assert(root != -1); // check tree is non-empty

        // check source
        if (dict[_id].nSource == 0) {
            dict[_id].isLeaf = true;
            dict[_id].isEmpty = true;
        } else {
            // divide
            if ((dict[_id].nSource <= rank) || (dict[_id].nLevel == _maxLevel)) {
                dict[_id].isLeaf = true;
                if (maxLevel < dict[_id].nLevel) {
                    maxLevel = dict[_id].nLevel;
                }
            } else {
                // not a leaf
                for (index_t i = 0; i < 8; ++i) {
                    maxId += 1;
                    dict[_id].child[i] = maxId;
                    dict.push_back(node(dict[_id].nLevel + 1, i));
                    dict[maxId].parent = _id;
                    dict[maxId].center.x = dict[_id].center.x + ((i & 1) - 0.5) * dict[_id].radius.x;
                    dict[maxId].center.y = dict[_id].center.y + (((i >> 1) & 1) - 0.5) * dict[_id].radius.y;
                    dict[maxId].center.z = dict[_id].center.z + ((i >> 2) - 0.5) * dict[_id].radius.z;
                    dict[maxId].radius.x = dict[_id].radius.x * 0.5;
                    dict[maxId].radius.y = dict[_id].radius.y * 0.5;
                    dict[maxId].radius.z = dict[_id].radius.z * 0.5;
                    dict[maxId].nSource = 0;
                    dict[maxId].nTarget = 0;
                }

                /*
                 * can be accelerated by **reduce**
                 */
                for (index_t i = 0; i < dict[_id].nSource; ++i) {
                    index_t index = dict[_id].sourceIndex[i];
                    index_t z_bit = sourceTree[index].z < dict[_id].center.z ? 0 : 1;
                    index_t y_bit = sourceTree[index].y < dict[_id].center.y ? 0 : 1;
                    index_t x_bit = sourceTree[index].x < dict[_id].center.x ? 0 : 1;
                    index_t childIndex = 4 * z_bit + 2 * y_bit + x_bit;

                    index_t childId = dict[_id].child[childIndex];
                    dict[childId].sourceIndex.push_back(index);
                    dict[childId].nSource += 1;
                }

                /*
                 * can be accelerated by **reduce**
                 */
                for (index_t i = 0; i < dict[_id].nTarget; ++i) {
                    index_t index = dict[_id].targetIndex[i];
                    index_t z_bit = targetTree[index].z < dict[_id].center.z ? 0 : 1;
                    index_t y_bit = targetTree[index].y < dict[_id].center.y ? 0 : 1;
                    index_t x_bit = targetTree[index].x < dict[_id].center.x ? 0 : 1;
                    index_t childIndex = 4 * z_bit + 2 * y_bit + x_bit;

                    index_t childId = dict[_id].child[childIndex];
                    dict[childId].targetIndex.push_back(index);
                    dict[childId].nTarget += 1;
                }

                for (index_t i = 0; i < 8; ++i) {
                    assignChildren(dict[_id].child[i], _maxLevel);
                }
            }
        }
    }


    void tree::buildTree() {
        point min_p(dict[root].center.x - dict[root].radius.x,
                    dict[root].center.y - dict[root].radius.y,
                    dict[root].center.z - dict[root].radius.z);
        point max_p(dict[root].center.x + dict[root].radius.x,
                    dict[root].center.y + dict[root].radius.y,
                    dict[root].center.z + dict[root].radius.z);
        index_t i;
#ifdef RUN_OMP
#pragma omp parallel for private(i) shared(min_p, max_p) schedule(dynamic)
#endif
        for (i = 0; i < dict.size(); ++i) {
            buildNode(i, min_p, max_p);
        }
    }


    void tree::buildNode(index_t _id, point &min_p, point &max_p) {
        node &n = dict[_id];
        n.uList.clear();
        n.vList.clear();
        n.wList.clear();
        n.xList.clear();

        // not root
        if (n.parent != -1) {
            node &pn = dict[n.parent];
            scalar_t dx = n.radius.x;
            scalar_t dy = n.radius.y;
            scalar_t dz = n.radius.z;
            scalar_t xs = pn.center.x - dx;
            scalar_t ys = pn.center.y - dy;
            scalar_t zs = pn.center.z - dz;

            point cur;

            for (index_t x_id = -2; x_id < 4; x_id++) {
                for (index_t y_id = -2; y_id < 4; y_id++) {
                    for (index_t z_id = -2; z_id < 4; z_id++) {
                        cur.x = xs + 2 * x_id * dx;
                        cur.y = ys + 2 * y_id * dy;
                        cur.z = zs + 2 * z_id * dz;

                        // check box and not itself.
                        if (cur <= max_p && cur >= min_p && !(cur == n.center)) {
                            //find node.
                            index_t curId = findNode(0, cur);
                            bool adj = isAdjacent(_id, curId);
                            node &curNode = dict[curId];

                            if (curNode.nLevel < n.nLevel) {
                                if (adj) {
                                    if (curNode.isLeaf) {
                                        n.uList.insert(curId);
                                    }
                                } else {
                                    n.xList.insert(curId);
                                }
                            }

                            if (curNode.nLevel == n.nLevel) {
                                if (!adj) {
                                    n.vList.insert(curId);
                                } else {
                                    if (n.isLeaf) {
                                        std::queue<index_t> rest;
                                        rest.push(curId);
                                        while (!rest.empty()) {
                                            index_t frontId = rest.front();
                                            rest.pop();
                                            node &frontNode = dict[frontId];
                                            if (!isAdjacent(frontId, _id)) {
                                                n.wList.insert(frontId);
                                            } else {
                                                if (frontNode.isLeaf) {
                                                    n.uList.insert(frontId);
                                                } else {
                                                    for (index_t i = 0; i < 8; ++i) {
                                                        rest.push(frontNode.child[i]);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        if (n.isLeaf) {
            n.uList.insert(_id);
        }

        n.nUList = (index_t) n.uList.size();
        n.nWList = (index_t) n.wList.size();
        n.nVList = (index_t) n.vList.size();
        n.nXList = (index_t) n.xList.size();
    }


    index_t tree::findNode(index_t _id, point &p) {
        node &n = dict[_id];
        if (n.center == p) return _id;
        else {
            if (n.isLeaf) {
                return _id;
            } else {
                index_t x_bit = n.center.x > p.x ? 0 : 1;
                index_t y_bit = n.center.y > p.y ? 0 : 1;
                index_t z_bit = n.center.z > p.z ? 0 : 1;
                index_t id = 4 * z_bit + 2 * y_bit + x_bit;
                return findNode(n.child[id], p);
            }
        }
    }

    bool tree::isAdjacent(index_t _aId, index_t _bId) {
        node &nA = dict[_aId];
        node &nB = dict[_bId];
        scalar_t diff_x = fabs(nA.center.x - nB.center.x);
        scalar_t diff_y = fabs(nA.center.y - nB.center.y);
        scalar_t diff_z = fabs(nA.center.z - nB.center.z);
        scalar_t r_x = fabs(nA.radius.x + nB.radius.x);
        scalar_t r_y = fabs(nA.radius.y + nB.radius.y);
        scalar_t r_z = fabs(nA.radius.z + nB.radius.z);

        bool rdx = r_x >= diff_x - EPS;
        bool rdy = r_y >= diff_y - EPS;
        bool rdz = r_z >= diff_z - EPS;

        bool x_adj = (fabs(diff_x - r_x) < EPS) && (rdy && rdz);
        bool y_adj = (fabs(diff_y - r_y) < EPS) && (rdx && rdz);
        bool z_adj = (fabs(diff_z - r_z) < EPS) && (rdy && rdx);

        return x_adj || y_adj || z_adj;

    }


    void tree::output(std::string file) {
        std::ofstream file_stream(file);
        if (file_stream.is_open()) {

            for (size_t i = 0; i < dict.size(); ++i) {
                file_stream << dict[i].center.x << " "
                            << dict[i].center.y << " "
                            << dict[i].center.z << " "
                            << dict[i].radius.x << " "
                            << dict[i].radius.y << " "
                            << dict[i].radius.z << " "
                            << dict[i].nVList << " " << dict[i].nXList << " " << dict[i].nUList << " "
                            << dict[i].nWList << " " << dict[i].isLeaf << " " << dict[i].nSource << "\n";
            }

            file_stream.close();
        } else {
            std::cout << "cannot open file: " << file << std::endl;
        }
    }
}

