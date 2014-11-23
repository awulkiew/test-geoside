#include <GL/glut.h>

#include <boost/geometry.hpp>
#include <boost/geometry/extensions/gis/geographic/strategies/vincenty.hpp>

#include <iostream>
#include <vector>

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

struct point_info
{
    double lon, lat, azi_s, azi_p;
    int sph_s;
};

std::vector<point_info> points;

double lon_s1 = -51;
double lat_s1 = -51;
double lon_s2 = 51;
double lat_s2 = 51;
double lon_min = -50;
double lat_min = -50;
double lon_max = 50;
double lat_max = 50;
double step = 0.1;
double a = 1;
double b = 0.75;

double disp_x_min = -60;
double disp_y_min = -60;
double disp_x_max = 60;
double disp_y_max = 60;

void fill_points()
{
    typedef bgm::point<double, 2, bg::cs::geographic<bg::degree> > geo_point;
    typedef bgm::point<double, 2, bg::cs::spherical_equatorial<bg::degree> > sph_point;

    bg::srs::spheroid<double> sph(a, b);

    bgd::vincenty_inverse<double> vi(lon_s1 * bg::math::d2r,
                                     lat_s1 * bg::math::d2r,
                                     lon_s2 * bg::math::d2r,
                                     lat_s2 * bg::math::d2r,
                                     sph);
    double fwd = vi.azimuth12();
    normalize(fwd);

    for ( double x = lon_min ; x <= lon_max ; x += step )
    {
        for ( double y = lat_min ; y <= lat_max ; y += step )
        {
            bgd::vincenty_inverse<double> vi2(lon_s1 * bg::math::d2r,
                                              lon_s1 * bg::math::d2r,
                                              x * bg::math::d2r,
                                              y * bg::math::d2r,
                                              sph);

            double fwd2 = vi2.azimuth12();
            normalize(fwd2);

            bg::strategy::side::spherical_side_formula<double> ssf;
            int ss = ssf.apply(sph_point(lon_s1, lat_s1), sph_point(lon_s2, lat_s2), sph_point(x, y));

            point_info pi;
            pi.lon = x;
            pi.lat = y;
            pi.azi_s = fwd;
            pi.azi_p = fwd2;
            pi.sph_s = ss;

            points.push_back(pi);
        }
    }
}

void render_scene(void)
{
    glClear(GL_COLOR_BUFFER_BIT);

    double disp_width = disp_x_max - disp_x_min;
    double disp_height = disp_y_max - disp_y_min;
    double lon_width = lon_max - lon_min;
    double lat_height = lat_max - lat_min;

    double scale_x = disp_width / lon_width;
    double scale_y = disp_height / lat_height;

    size_t i = 0;
    for ( double x = lon_min ; x <= lon_max ; x += step )
    {
        for ( double y = lat_min ; y <= lat_max ; y += step )
        {
            point_info pi = points[i++];

            if ( pi.azi_p > pi.azi_s ) // right
            {
                if ( pi.sph_s < 0 ) // right
                    glColor3f(0.25, 0, 0);
                else
                    glColor3f(1, 0, 1);
            }
            else // left
            {
                if ( pi.sph_s > 0 ) // left
                    glColor3f(0, 0.25, 0);
                else
                    glColor3f(0, 1, 1);
            }
            
            double xt = (x - lon_min) * scale_x + disp_x_min;
            double yt = (y - lat_min) * scale_y + disp_y_min;
            double step_h_x = step / 2 * scale_x;
            double step_h_y = step / 2 * scale_y;

            glBegin(GL_QUADS);
            glVertex3f(xt-step_h_x, yt-step_h_y, 0);
            glVertex3f(xt+step_h_x, yt-step_h_y, 0);
            glVertex3f(xt+step_h_x, yt+step_h_y, 0);
            glVertex3f(xt-step_h_x, yt+step_h_y, 0);
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
    std::cout << "filling points" << std::endl;
    fill_points();
    std::cout << "visualizing" << std::endl;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_SINGLE | GLUT_RGBA);
    glutInitWindowPosition(100,100);
    glutInitWindowSize(800, 800);
    glutCreateWindow("test-geoside");

    glutDisplayFunc(render_scene);
    glutReshapeFunc(resize);
    glutMouseFunc(mouse);
    glutKeyboardFunc(keyboard);

    glutMainLoop();

    return 0;
}
