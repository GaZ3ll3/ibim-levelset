//
// Created by lurker on 2/13/17.
//

#include "view.h"

void view::_cube(point& p, scalar_t r) {
    GLfloat front  = (GLfloat) (p.data[0] - r);
    GLfloat back   = (GLfloat) (p.data[0] + r);
    GLfloat left   = (GLfloat) (p.data[1] - r);
    GLfloat right  = (GLfloat) (p.data[1] + r);
    GLfloat bottom = (GLfloat) (p.data[2] - r);
    GLfloat top    = (GLfloat) (p.data[2] + r);

    glPushMatrix();

    glBegin(GL_POLYGON);
    glColor3f(0.0, 1.0, 0.0);
    glVertex3f(right, bottom, front);
    glVertex3f(right, bottom, front);
    glVertex3f(right, top,    front);
    glVertex3f(left,  top,    front);
    glVertex3f(left,  bottom, front);
    glEnd();

    glBegin(GL_POLYGON);
    glColor3f(1.0, 0.0, 0.0);
    glVertex3f(right, bottom, back);
    glVertex3f(right, top,    back);
    glVertex3f(left,  top,    back);
    glVertex3f(left,  bottom, back);
    glEnd();

    glBegin(GL_POLYGON);
    glColor3f(0.0, 0.0, 1.0);
    glVertex3f(right, bottom, front);
    glVertex3f(right, top,    front);
    glVertex3f(right, top,    back);
    glVertex3f(right, bottom, back);
    glEnd();

    glBegin(GL_POLYGON);
    glColor3f(1.0, 1.0, 0.0);
    glVertex3f(left, bottom, back);
    glVertex3f(left, top,    back);
    glVertex3f(left, top,    front);
    glVertex3f(left, bottom, front);
    glEnd();

    glBegin(GL_POLYGON);
    glColor3f(0.0, 1.0, 1.0);
    glVertex3f(right, top, back);
    glVertex3f(right, top, front);
    glVertex3f(left,  top, front);
    glVertex3f(left,  top, back);
    glEnd();

    glBegin(GL_POLYGON);
    glColor3f(1.0, 0.0, 1.0);
    glVertex3f(right, bottom, front);
    glVertex3f(right, bottom, back);
    glVertex3f(left,  bottom, back);
    glVertex3f(left,  bottom, front);
    glEnd();

    glPopMatrix();
}

void view::_key(unsigned char k, int x, int y) {
    if (k == 'd') rotate_y += rotate_step;
    else if (k == 'a') rotate_y -= rotate_step;
    else if (k == 'w') rotate_x += rotate_step;
    else if (k == 's') rotate_x -= rotate_step;
    else if (k == 'q') zoom *= (1.0 + zoom_step);
    else if (k == 'e') zoom *= (1.0 - zoom_step);

    glutPostRedisplay();
}

void view::_display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();
    glRotatef((GLfloat) rotate_x, 1.0, 0.0, 0.0);
    glRotatef((GLfloat) rotate_y, 0.0, 1.0, 0.0);
    glScaled(zoom, zoom, zoom);


    _levelset();


    glFlush();
    glutSwapBuffers();

}


void view::_levelset() {
    int iX, iY, iZ;
    int lx = ls->Nx;
    int ly = ls->Ny;
    int lz = ls->Nz;
    double dx = ls->dx;
    double sx = ls->sx;
    double sy = ls->sy;
    double sz = ls->sz;
    int count = 0;

    if (!cached) {
        for (iX =1; iX < lx-1; ++iX) {
            for (iY = 1; iY < ly-1; ++iY) {
                for (iZ = 1; iZ < lz-1; ++iZ) {
                    if (fabs(phi->get(iX, iY, iZ)) < 3 * dx) {
                        count++;
                        point P = {
                                sx + iX * dx,
                                sy + iY * dx,
                                sz + iZ * dx
                        };
                        scalar_t r = dx / 2.0;
                        cachePoints.push_back(P);
                        cacheRadius.push_back(r);
                    }
                }
            }
        }
        cached = true;
        std::cout << "boundary size : " << cachePoints.size() << std::endl;
    }
    else {
        for (index_t i = 0; i < (index_t )cachePoints.size(); ++i) {
//            if (cachePoints[i].data[1] > -20)
            _cube(cachePoints[i], cacheRadius[i]);
        }
    }


}


void view::run() {
    int argc = 1;
    char* argv[1]; argv[0] = strdup(TITLE);
    glutInit(&argc, argv);

    glutInitWindowSize(WIDTH, HEIGHT);
    glutInitWindowPosition((glutGet(GLUT_SCREEN_WIDTH) - WIDTH)/2,
                           (glutGet(GLUT_SCREEN_HEIGHT) - HEIGHT)/2);

    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

    glutCreateWindow(TITLE);
    glEnable(GL_DEPTH_TEST);

    glutKeyboardFunc(view::key_callback);
    glutDisplayFunc(view::display_callback);
    glutReshapeFunc(view::reshape_callback);

    glutMainLoop();

    return;
}