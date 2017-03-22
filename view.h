//
// Created by lurker on 2/13/17.
//

#ifndef IBIM_OCTREE_VIEW_H
#define IBIM_OCTREE_VIEW_H

#include "levelset.h"
#include <GL/glut.h>

#define TITLE "VIEW"
#define WIDTH 800
#define HEIGHT 800

class view {
private:
    Grid* phi;
    levelset* ls;
    Surface* surf;
    bool cached;
    char tag;
    vector<point> cachePoints;
    vector<scalar_t> cacheRadius;

    void _cube(point& p, scalar_t r);
    void _key(unsigned char k, int x, int y);
    void _display();
    void _levelset();

    view(char flag = 0) {
        rotate_x = 0.;
        rotate_y = 0.;
        zoom     = 0.003;
        rotate_step = 5.0;
        zoom_step = 0.1;
        cached = false;
        tag = flag;

    }

    ~view() {
#ifdef DEBUG
        std::cout << "destructing view." << std::endl;
#endif
        phi = nullptr; ls = nullptr;
        surf = nullptr;
    }


protected:
    scalar_t rotate_x;
    scalar_t rotate_y;
    scalar_t zoom;
    scalar_t rotate_step;
    scalar_t zoom_step;

public:


    static view& getInstance(char flag = 2) {
        static view viewInstance(flag);
        return viewInstance;
    }

    void loadLevelSet(Grid& _phi, levelset& _ls, Surface& _surf) {
        phi = &_phi;
        ls =&_ls;
        surf = &_surf;
    }

    void run();

    static void key_callback(unsigned char k, int x, int y) {
        getInstance()._key(k, x, y);
    }

    static void reshape_callback(int width, int height) {
        if (height <= 0) height = 1;

        glViewport(0, 0, width, height);
        glOrtho(0., width, 0., height, 0, 50000);
    }

    static void display_callback() {
        getInstance()._display();

    }
};


#endif //IBIM_OCTREE_VIEW_H
