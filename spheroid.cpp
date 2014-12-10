#include <GL/glut.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/register/point.hpp>

#include <boost/geometry/strategies/geographic/ssf.hpp>

#include <iostream>
#include <vector>

namespace bg = boost::geometry;
namespace bgm = bg::model;
namespace bgd = bg::detail;

typedef bgm::point<double, 2, bg::cs::geographic<bg::degree> > point_geo;
typedef bgm::point<double, 2, bg::cs::spherical_equatorial<bg::degree> > point_sph;
typedef bgm::point<double, 3, bg::cs::cartesian> point_3d;
typedef bg::srs::spheroid<double> spheroid;

double a = 1;
double b = 0.5;
spheroid sph(a, b);

double pi = bg::math::pi<double>();

struct color
{
    color()
        : r(0), g(0), b(0), a(1)
    {}
    color(float c1, float c2, float c3, float c4 = 1.0f)
        : r(c1), g(c2), b(c3), a(c4)
    {}

    float r, g, b, a;
};

template <typename P>
void convert(P const& p, point_3d & res)
{
    typedef typename bg::coordinate_type<P>::type ct;
    ct lon = bg::get_as_radian<0>(p);
    ct lat = bg::get_as_radian<1>(p);
    ct cos_lat = cos(lat);
    bg::set<0>(res, bg::get_radius<0>(sph) * cos_lat * cos(lon));
    bg::set<1>(res, bg::get_radius<0>(sph) * cos_lat * sin(lon));
    bg::set<2>(res, bg::get_radius<2>(sph) * sin(lat));
}

void convert(point_3d const& p, point_3d & res)
{
    res = p;
}

template <typename Res, typename P>
Res pcast(P const& p)
{
    Res res;
    convert(p, res);
    return res;
}

point_geo projected_to_equator(point_geo const& p)
{
    point_geo res = p;
    bg::set<1>(res, 0);
    return res;
}

point_3d projected_to_xy(point_geo const& p)
{
    point_3d res;
    convert(p, res);
    bg::set<2>(res, 0);
    return res;
}

point_3d projected_to_xy_geod(point_geo const& p)
{
    double r = 0;
    double lon = bg::get_as_radian<0>(p);
    double lat = bg::get_as_radian<1>(p);

    point_3d p3d;
    convert(p, p3d);

    // lat == 0
    if ( bg::math::equals(bg::get<2>(p3d), 0) )
        r = bg::get_radius<0>(sph);
    // |lat| == pi/2
    else if ( bg::math::equals(bg::get<0>(p3d), 0) && bg::math::equals(bg::get<1>(p3d), 0) )
        r = 0;
    // sqrt(xx+yy) - |z|/tan(lat)
    else
        r = bg::math::sqrt(bg::get<0>(p3d) * bg::get<0>(p3d) + bg::get<1>(p3d) * bg::get<1>(p3d))
            - bg::math::abs(bg::get<2>(p3d)) / tan(lat);

    point_3d res;
    bg::set<0>(res, r * cos(lon));
    bg::set<1>(res, r * sin(lon));
    bg::set<2>(res, 0);

    return res;
}

void draw_meridian(double lon, double step = pi/32)
{
    glBegin(GL_LINE_STRIP);
    for (double lat = -pi/2 ; lat < pi/2+step ; lat += step)
    {
        double cos_lat = cos(lat);
        double x = a * cos_lat * cos(lon);
        double y = a * cos_lat * sin(lon);
        double z = b * sin(lat);
        glVertex3f(x, y, z);
    }
    glEnd();
}

void draw_parallel(double lat, double step = pi/32)
{
    glBegin(GL_LINE_STRIP);
    for (double lon = 0 ; lon < 2*pi+step ; lon += step)
    {
        double cos_lat = cos(lat);
        double x = a * cos_lat * cos(lon);
        double y = a * cos_lat * sin(lon);
        double z = b * sin(lat);
        glVertex3f(x, y, z);
    }
    glEnd();
}

void draw_sphere(double x, double y, double z, double r)
{
    glPushMatrix();
    glTranslatef(x, y, z);
    GLUquadricObj * sph = gluNewQuadric();
    gluSphere(sph, r, 8, 8);
    gluDeleteQuadric(sph);
    glPopMatrix();
}

void draw_point(point_3d const& p)
{
    draw_sphere(bg::get<0>(p), bg::get<1>(p), bg::get<2>(p), 0.03);
}

template <typename P>
void draw_point(P const& p)
{
    point_3d p3d;
    convert(p, p3d);
    draw_point(p3d);
}

void draw_line(point_3d const& p1, point_3d const& p2)
{
    glBegin(GL_LINES);
    glVertex3f(bg::get<0>(p1), bg::get<1>(p1), bg::get<2>(p1));
    glVertex3f(bg::get<0>(p2), bg::get<1>(p2), bg::get<2>(p2));
    glEnd();
}

template <typename P1, typename P2>
void draw_line(P1 const& p1, P2 const& p2)
{
    point_3d p3d1, p3d2;
    convert(p1, p3d1);
    convert(p2, p3d2);
    draw_line(p3d1, p3d2);
}

void draw_axes_rgb(double len)
{
    glBegin(GL_LINES);
    glColor3f(1, 0, 0);
    glVertex3f(0, 0, 0);
    glVertex3f(len, 0, 0);
    glColor3f(0, 1, 0);
    glVertex3f(0, 0, 0);
    glVertex3f(0, len, 0);
    glColor3f(0, 0, 1);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, len);
    glEnd();
}

void draw_model()
{
    draw_axes_rgb(1.5);

    glColor3f(0.5, 0.5, 0.5);
    for (double lon = 0 ; lon < 2*pi ; lon += pi/8)
        draw_meridian(lon);

    for (double lat = -pi/2+pi/8 ; lat < pi/2 ; lat += pi/8)
        draw_parallel(lat);

    glColor3f(1, 1, 1);
    draw_point(point_3d(0, 0, 0));
    glColor3f(1, 0, 0);
    draw_point(point_3d(a, 0, 0));
    glColor3f(0, 1, 0);
    draw_point(point_3d(0, a, 0));
    glColor3f(0, 0, 1);
    draw_point(point_3d(0, 0, b));
}

void draw_point_adv(point_geo const& p, color const& c)
{
    point_3d p_xy = projected_to_xy(p);
    point_3d p_xy_geod = projected_to_xy_geod(p);

    glColor4f(c.r, c.g, c.b, c.a);
    draw_point(p);
    draw_point(p_xy_geod);

    draw_line(point_3d(0, 0, 0), p_xy);
    draw_line(p_xy_geod, p);
}

void render_scene(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    draw_model();

    draw_point_adv(point_geo(60, 30), color(1, 0.5, 0));
    draw_point_adv(point_geo(30, 10), color(1, 1, 0));

//    glColor4f(0, 0, 0, 0.5);
//    glBegin(GL_QUADS);
//    glVertex3f(-1, -1, -0.001);
//    glVertex3f(1, -1, -0.001);
//    glVertex3f(1, 1, -0.001);
//    glVertex3f(-1, 1, -0.001);
//    glEnd();

    glFlush();
}

void resize(int w, int h)
{
    if ( h == 0 )
        h = 1;

    float ratio = float(w) / h;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glViewport(0, 0, w, h);

    gluPerspective(45, ratio, 0.1, 10);
    //glOrtho(-5, 5, -5, 5,-5, 5);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(
        3.0f, 3.0f, 3.0f,
        0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 1.0f);

    glClearDepth(1.0f);
    glDepthFunc(GL_LESS);
    glEnable(GL_DEPTH_TEST);

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable( GL_BLEND );

    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

    glLineWidth(1.5f);
}

void mouse(int /*button*/, int /*state*/, int /*x*/, int /*y*/)
{
}

void keyboard(unsigned char /*key*/, int /*x*/, int /*y*/)
{
}

int main(int argc, char **argv)
{
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
