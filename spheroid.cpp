#include <GL/glut.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/register/point.hpp>

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
    bg::set<0>(res, a * cos_lat * cos(lon));
    bg::set<1>(res, a * cos_lat * sin(lon));
    bg::set<2>(res, b * sin(lat));
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

        double eq = (x*x+y*y)/(a*a) + z*z/(b*b);
        assert(bg::math::equals(eq, 1.0));

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

        double eq = (x*x+y*y)/(a*a) + z*z/(b*b);
        assert(bg::math::equals(eq, 1.0));

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
    draw_sphere(bg::get<0>(p), bg::get<1>(p), bg::get<2>(p), 0.02);
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

point_3d operator+(point_3d const& p1, point_3d const& p2) { point_3d res = p1; bg::add_point(res, p2); return res; }
point_3d operator-(point_3d const& p1, point_3d const& p2) { point_3d res = p1; bg::subtract_point(res, p2); return res; }
point_3d operator*(point_3d const& p1, double v) { point_3d res = p1; bg::multiply_value(res, v); return res; }
point_3d operator/(point_3d const& p1, double v) { point_3d res = p1; bg::divide_value(res, v); return res; }

float yaw = 0;
float pitch = 0;
float zoom = 0;

void render_scene()
{
    glPushMatrix();

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glTranslatef(zoom, zoom, zoom);
    glRotatef(pitch, -0.7, 0.7, 0);
    glRotatef(yaw, 0, 0, 1);

    draw_model();


    point_geo p1(60, 30);
    point_geo p2(30, 10);

    draw_point_adv(p1, color(1, 0.5, 0));
    draw_point_adv(p2, color(1, 1, 0));

    point_3d p1_3d = pcast<point_3d>(p1);
    point_3d p2_3d = pcast<point_3d>(p2);
    point_3d p01_3d = projected_to_xy_geod(p1);
    point_3d p02_3d = projected_to_xy_geod(p2);

    point_3d v_surface = p2_3d - p1_3d;
    point_3d v_equator = p02_3d - p01_3d;

    glColor3f(1, 1, 1);
    draw_line(p1_3d, p2_3d);
    draw_line(p01_3d, p02_3d);
    
    std::vector<point_3d> curve;

    double f = 0;
    for ( int i = 0 ; i <= 10 ; ++i, f += 0.1 )
    {
        point_3d pe = p01_3d + v_equator * f;
        point_3d ps = p1_3d + v_surface * f;
        
        glColor3f(1, 0.5+0.5*f, 0);
        draw_line(pe, ps);

        point_3d d = ps - pe;
        
        //(x*x+y*y) / a + z*z/b = 1
        // x = o.x + d.x * t
        // y = o.y + d.y * t
        // z = o.z + d.z * t        
        double ox = bg::get<0>(pe);
        double oy = bg::get<1>(pe);
        double oz = bg::get<2>(pe);
        double dx = bg::get<0>(d);
        double dy = bg::get<1>(d);
        double dz = bg::get<2>(d);

        double a_sqr = a*a;
        double b_sqr = b*b;
        double param_a = (dx*dx+dy*dy)/a_sqr+dz*dz/b_sqr;
        double param_b = 2*((ox*dx+oy*dy)/a_sqr+oz*dz/b_sqr);
        double param_c = (ox*ox+oy*oy)/a_sqr+oz*oz/b_sqr-1;

        double delta = param_b*param_b-4*param_a*param_c;
        double t = delta >= 0 ?
                   (-param_b+sqrt(delta)) / (2*param_a) :
                   0.0;

        point_3d p_curve = pe + d * t;
        
        curve.push_back(p_curve);
    }

    f = 0.1;
    for ( int i = 1 ; i <= 10 ; ++i, f += 0.1 )
    {
        glColor3f(1, 0.5+0.5*f, 1);
        draw_line(curve[i-1], curve[i]);
    }

    glPopMatrix();
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

    glLineWidth(2);
}

bool left_down = false;
bool right_down = false;
int down_x;
int down_y;
float prev_yaw;
float prev_pitch;
float prev_zoom;
void mouse(int button, int state, int x, int y)
{
    bool down = state == 0;
    if ( button == 0 )
        left_down = down;
    if ( button == 2 )
        right_down = down;

    down_x = x;
    down_y = y;
    prev_yaw = yaw;
    prev_pitch = pitch;
    prev_zoom = zoom;
}
void mouse_move(int x, int y)
{
    if ( left_down )
    {
        yaw = prev_yaw + (x - down_x) / 10.0f;
        pitch = prev_pitch + (y - down_y) / 10.0f;
    }
    
    if ( right_down )
    {
        zoom = prev_zoom + (y - down_y) / 100.0f;
    }
}

void keyboard(unsigned char /*key*/, int /*x*/, int /*y*/)
{
}

void idle_fun()
{
    glutPostRedisplay();
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_SINGLE | GLUT_RGBA);
    glutInitWindowPosition(100,100);
    glutInitWindowSize(800, 800);
    glutCreateWindow("test-geoside");

    glutDisplayFunc(render_scene);
    glutIdleFunc(idle_fun);
    glutReshapeFunc(resize);
    glutMouseFunc(mouse);
    glutMotionFunc(mouse_move);
    glutKeyboardFunc(keyboard);

    glutMainLoop();

    return 0;
}
