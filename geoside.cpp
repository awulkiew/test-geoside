#include <GL/glut.h>

#include <boost/geometry.hpp>
#include <boost/geometry/extensions/gis/geographic/strategies/vincenty.hpp>

namespace bg = boost::geometry;
namespace bgm = bg::model;
namespace bgd = bg::detail;

template <typename T>
void normalize(T & a)
{
    while ( a > bg::math::pi<T>() )
        a -= bg::math::pi<T>() * 2;
    while ( a <= -bg::math::pi<T>() )
        a += bg::math::pi<T>() * 2;
}

void render_scene(void)
{
    glClear(GL_COLOR_BUFFER_BIT);

    typedef bgm::point<double, 2, bg::cs::geographic<bg::degree> > geo_point;
    typedef bgm::point<double, 2, bg::cs::spherical_equatorial<bg::degree> > sph_point;
    
    bg::srs::spheroid<double> sph(1, 0.75);

    bgd::vincenty_inverse<double> vi(-51 * bg::math::d2r,
                                     -51 * bg::math::d2r,
                                     51 * bg::math::d2r,
                                     51 * bg::math::d2r,
                                     sph);
    double fwd = vi.azimuth12();
    double bck = vi.azimuth21();

    double step = 0.1;
    for ( double x = -50 ; x <= 50 ; x += step )
    {
        for ( double y = -50 ; y <= 50 ; y += step )
        {
            bgd::vincenty_inverse<double> vi2(-51 * bg::math::d2r,
                                              -51 * bg::math::d2r,
                                              x * bg::math::d2r,
                                              y * bg::math::d2r,
                                              sph);

            bg::strategy::side::spherical_side_formula<double> ssf;
            int ss = ssf.apply(sph_point(-51, -51), sph_point(51, 51), sph_point(x, y));

            double fwd2 = vi2.azimuth12();
            double bck2 = vi2.azimuth21();

            normalize(fwd);
            normalize(fwd2);

            if ( fwd2 > fwd ) // right
            {
                if ( ss < 0 ) // right
                    glColor3f(0.25, 0, 0);
                else
                    glColor3f(1, 0, 1);
            }
            else // left
            {
                if ( ss > 0 ) // left
                    glColor3f(0, 0.25, 0);
                else
                    glColor3f(0, 1, 1);
            }
            
            glBegin(GL_QUADS);
            glVertex3f(x-step/2, y-step/2, 0);
            glVertex3f(x+step/2, y-step/2, 0);
            glVertex3f(x+step/2, y+step/2, 0);
            glVertex3f(x-step/2, y+step/2, 0);
            glEnd();
        }
    }

    glFlush();
}

void resize(int w, int h)
{
    if ( h == 0 )
        h = 1;

    //float ratio = float(w) / h;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glViewport(0, 0, w, h);

    //gluPerspective(45, ratio, 1, 1000);
    glOrtho(-70, 70, -70, 70, -10, 100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    /*gluLookAt(
        120.0f, 120.0f, 120.0f, 
        50.0f, 50.0f, -1.0f,
        0.0f, 1.0f, 0.0f);*/
    gluLookAt(
        0.0f, 0.0f, 100.0f, 
        0.0f, 0.0f, -1.0f,
        0.0f, 1.0f, 0.0f);

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glLineWidth(1.5f);
}

void mouse(int button, int state, int /*x*/, int /*y*/)
{
}

void keyboard(unsigned char key, int /*x*/, int /*y*/)
{
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_SINGLE | GLUT_RGBA);
    glutInitWindowPosition(100,100);
    glutInitWindowSize(800, 800);
    glutCreateWindow("boost::geometry::index::rtree GLUT test");

    glutDisplayFunc(render_scene);
    glutReshapeFunc(resize);
    glutMouseFunc(mouse);
    glutKeyboardFunc(keyboard);

    glutMainLoop();

    return 0;
}