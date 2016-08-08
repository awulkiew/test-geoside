#include <GL/glut.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/linestring.hpp>

#include <boost/geometry/strategies/geographic/mapping_ssf.hpp>

#include <boost/array.hpp>

#include <iostream>
#include <vector>

namespace bg = boost::geometry;
namespace bgm = bg::model;
namespace bgd = bg::detail;
namespace bgf = bg::formula;

double const d2r = bg::math::d2r<double>();
double const r2d = bg::math::r2d<double>();

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
    double lon, lat, azi_s, azi_p, azi_s_bck, azi_p_bck;
    int sph_s, car_s;
};
// for simplify
BOOST_GEOMETRY_REGISTER_POINT_2D(point_info, double, bg::cs::cartesian, lon, lat)

std::vector<point_info> points;

double lon_s1 = -31;
double lat_s1 = -31;
double lon_s2 = 31;
double lat_s2 = 31;
double step = 0.1;
double lon_min = lon_s1 - 30;
double lat_min = lat_s1 - 30;
double lon_max = lon_s2 + 30;
double lat_max = lat_s2 + 30;
double a = 1;
double b = 1 - 1.0/300;
bg::srs::spheroid<double> spheroid(a, b);

bg::strategy::distance::vincenty<bg::srs::spheroid<double> > vincenty(spheroid);
bg::strategy::distance::thomas<bg::srs::spheroid<double> > thomas(spheroid);
bg::strategy::distance::andoyer<bg::srs::spheroid<double> > andoyer(spheroid);
bg::strategy::distance::haversine<double> haversine_mean((2*a+b)/3);
bg::strategy::distance::haversine<double> haversine_max(a);

typedef bgf::vincenty_inverse<double, true, true, true> vincenty_inv;

double simplify_distance = 0.066;

double disp_x_min = -60;
double disp_y_min = -60;
double disp_x_max = 60;
double disp_y_max = 60;

void fill_points()
{
    typedef bgm::point<double, 2, bg::cs::geographic<bg::degree> > geo_point;
    typedef bgm::point<double, 2, bg::cs::spherical_equatorial<bg::degree> > sph_point;
    typedef bgm::point<double, 2, bg::cs::cartesian> car_point;

    bgf::result_inverse<double>
        vi = vincenty_inv::apply(lon_s1 * d2r,
                                 lat_s1 * d2r,
                                 lon_s2 * d2r,
                                 lat_s2 * d2r,
                                 spheroid);
    double fwd = vi.azimuth;
    normalize(fwd);
    double bck = vi.reverse_azimuth;
    normalize(bck);

    for ( double y = lat_min ; y <= lat_max ; y += step )
    {
        for ( double x = lon_min ; x <= lon_max ; x += step )    
        {
            bgf::result_inverse<double>
                vi2 = vincenty_inv::apply(lon_s1 * d2r,
                                          lat_s1 * d2r,
                                          x * d2r,
                                          y * d2r,
                                          spheroid);

            double fwd2 = vi2.azimuth;
            normalize(fwd2);
            double bck2 = vi2.reverse_azimuth;
            normalize(bck2);

            boost::array<int, 4> ss;

            bg::strategy::side::spherical_side_formula<double> ssf;
            ss[0] = ssf.apply(sph_point(lon_s1, lat_s1), sph_point(lon_s2, lat_s2), sph_point(x, y));

            bg::strategy::side::mapping_spherical_side_formula
                <
                    bg::srs::spheroid<double>,
                    bg::strategy::side::mapping_geodetic
                > ssf1(spheroid);
            ss[1] = ssf1.apply(geo_point(lon_s1, lat_s1), geo_point(lon_s2, lat_s2), geo_point(x, y));

            bg::strategy::side::mapping_spherical_side_formula
                <
                    bg::srs::spheroid<double>,
                    bg::strategy::side::mapping_reduced
                > ssf2(spheroid);
            ss[2] = ssf2.apply(geo_point(lon_s1, lat_s1), geo_point(lon_s2, lat_s2), geo_point(x, y));

            bg::strategy::side::mapping_spherical_side_formula
                <
                    bg::srs::spheroid<double>,
                    bg::strategy::side::mapping_geocentric
                > ssf3(spheroid);
            ss[3] = ssf3.apply(geo_point(lon_s1, lat_s1), geo_point(lon_s2, lat_s2), geo_point(x, y));

            for(int i = 0 ; i < 4 ; ++i)
            {
                for(int j = 1 ; j < 4 ; ++j)
                {
                    if ( ss[i] != ss[j] )
                    {
                        std::cout << "different SSF res for: " << i << " and " << j << std::endl;
                    }
                }
            }

            bg::strategy::side::side_by_triangle<> sbt;
            int cs = sbt.apply(car_point(lon_s1, lat_s1), car_point(lon_s2, lat_s2), car_point(x, y));

            point_info pi;
            pi.lon = x;
            pi.lat = y;
            pi.azi_s = fwd;
            pi.azi_p = fwd2;
            pi.azi_s_bck = bck;
            pi.azi_p_bck = bck2;
            pi.sph_s = ss[0];
            pi.car_s = cs;

            points.push_back(pi);

            boost::ignore_unused(sbt, ssf);
        }
    }
}

bg::model::linestring<point_info> path_geo;
bg::model::linestring<point_info> path_sph;
bg::model::linestring<point_info> path_car;

void measure_paths()
{
    double lin_a = (lat_s2 - lat_s1) / (lon_s2 - lon_s1);
    double lin_b = lat_s1 - lin_a * lon_s1;

    size_t i = 0;
    for ( double y = lat_min ; y <= lat_max ; y += step )
    {
        bool last_left_geo = true;
        bool last_left_sph = true;
        bool last_left_car = true;

        for ( double x = lon_min ; x <= lon_max ; x += step )    
        {
            point_info pi = points[i++];
            
            if ( x < lon_s1 || x > lon_s2 || y < lat_s1 || y > lat_s2 )
            {
                continue;
            }

            bool is_s1 = lon_s1 == pi.lon && lat_s1 == pi.lat;
            bool is_s2 = lon_s2 == pi.lon && lat_s2 == pi.lat;

            if ( !is_s1 && !is_s2 )
            {
                pi.lon -= step/2; // the coordinate between left and right point
            }

            bool is_geo_right = ::sin(pi.azi_p-pi.azi_s) >= 0;
            //bool is_geo_bck_right = ::sin(pi.azi_p_bck-pi.azi_s_bck) <= 0;
            bool is_sph_right = pi.sph_s < 0;
            bool is_car_right = pi.car_s < 0;

            if ( is_geo_right || is_s1 ) // right
            {
                if ( last_left_geo )
                {
                    path_geo.push_back(pi);
                    last_left_geo = false;
                }
            }

            if ( is_sph_right || is_s1 ) // right
            {
                if ( last_left_sph )
                {
                    path_sph.push_back(pi);
                    last_left_sph = false;
                }
            }

            if ( is_car_right || is_s1 )
            {
                // sanity check
                if ( !(y < lin_a * x + lin_b) )
                {
                    std::cerr << "Error: cartesian right not compatible with line equation." << std::endl;
                }

                if ( last_left_car )
                {
                    path_car.push_back(pi);
                    last_left_car = false;
                }
            }
        }
    }

    bg::model::linestring<point_info> temp;
    bg::simplify(path_geo, temp, simplify_distance);
    path_geo = temp;
    bg::simplify(path_sph, temp, simplify_distance);
    path_sph = temp;
    // simplify cartesian path but make sure that it contain similar number of points
    std::size_t count = path_car.size() / (std::min)(path_geo.size(), path_sph.size());
    temp.clear();
    for ( size_t i = 0 ; i < path_car.size() ; i += count )
    {
        temp.push_back(path_car[i]);
    }
    if ( temp.back().lon != path_car.back().lon || temp.back().lat != path_car.back().lat )
    {
        temp.push_back(path_car.back());
    }
    path_car = temp;

    typedef bgm::point<double, 2, bg::cs::geographic<bg::degree> > geo_point;
    typedef bgm::point<double, 2, bg::cs::spherical_equatorial<bg::degree> > sph_point;

    double distance_geo1 = 0;
    double distance_geo2 = 0;
    double distance_geo3 = 0;
    double distance_geo4 = 0;
    for ( size_t i = 1 ; i < path_geo.size() ; ++i )
    {
        point_info const& pi1 = path_geo[i-1];
        point_info const& pi2 = path_geo[i];
        distance_geo1 += bg::distance(geo_point(pi1.lon, pi1.lat),
                                      geo_point(pi2.lon, pi2.lat),
                                      vincenty);
        distance_geo2 += bg::distance(geo_point(pi1.lon, pi1.lat),
                                      geo_point(pi2.lon, pi2.lat),
                                      andoyer);
        distance_geo3 += bg::distance(sph_point(pi1.lon, pi1.lat),
                                      sph_point(pi2.lon, pi2.lat),
                                      haversine_mean);
        distance_geo4 += bg::distance(sph_point(pi1.lon, pi1.lat),
                                      sph_point(pi2.lon, pi2.lat),
                                      haversine_max);
    }

    double distance_sph1 = 0;
    double distance_sph2 = 0;
    double distance_sph3 = 0;
    double distance_sph4 = 0;
    for ( size_t i = 1 ; i < path_sph.size() ; ++i )
    {
        point_info const& pi1 = path_sph[i-1];
        point_info const& pi2 = path_sph[i];
        distance_sph1 += bg::distance(geo_point(pi1.lon, pi1.lat),
                                      geo_point(pi2.lon, pi2.lat),
                                      vincenty);
        distance_sph2 += bg::distance(geo_point(pi1.lon, pi1.lat),
                                      geo_point(pi2.lon, pi2.lat),
                                      andoyer);
        distance_sph3 += bg::distance(sph_point(pi1.lon, pi1.lat),
                                      sph_point(pi2.lon, pi2.lat),
                                      haversine_mean);
        distance_sph4 += bg::distance(sph_point(pi1.lon, pi1.lat),
                                      sph_point(pi2.lon, pi2.lat),
                                      haversine_max);
    }

    double distance_car1 = 0;
    double distance_car2 = 0;
    double distance_car3 = 0;
    double distance_car4 = 0;
    for ( size_t i = 1 ; i < path_car.size() ; ++i )
    {
        point_info const& pi1 = path_car[i-1];
        point_info const& pi2 = path_car[i];
        distance_car1 += bg::distance(geo_point(pi1.lon, pi1.lat),
                                      geo_point(pi2.lon, pi2.lat),
                                      vincenty);
        distance_car2 += bg::distance(geo_point(pi1.lon, pi1.lat),
                                      geo_point(pi2.lon, pi2.lat),
                                      andoyer);
        distance_car3 += bg::distance(sph_point(pi1.lon, pi1.lat),
                                      sph_point(pi2.lon, pi2.lat),
                                      haversine_mean);
        distance_car4 += bg::distance(sph_point(pi1.lon, pi1.lat),
                                      sph_point(pi2.lon, pi2.lat),
                                      haversine_max);
    }

    std::cout << "(V) - vincenty, (T) - thomas, (A) - andoyer, (H) - haversine (mean radius), (M) - haversine (max radius)" << std::endl;

    std::cout << "distance (V) = " << std::setprecision(32)
              << bg::distance(geo_point(lon_s1, lat_s1), geo_point(lon_s2, lat_s2), vincenty) << std::endl;
    std::cout << "distance (T) = " << std::setprecision(32)
              << bg::distance(geo_point(lon_s1, lat_s1), geo_point(lon_s2, lat_s2), thomas) << std::endl;
    std::cout << "distance (A) = " << std::setprecision(32)
              << bg::distance(geo_point(lon_s1, lat_s1), geo_point(lon_s2, lat_s2), andoyer) << std::endl;
    std::cout << "distance (H) = " << std::setprecision(32)
              << bg::distance(sph_point(lon_s1, lat_s1), sph_point(lon_s2, lat_s2), haversine_mean) << std::endl;
    std::cout << "distance (M) = " << std::setprecision(32)
              << bg::distance(sph_point(lon_s1, lat_s1), sph_point(lon_s2, lat_s2), haversine_max) << std::endl;

    std::cout << "length geo (V) = " << std::setprecision(32) << distance_geo1 << std::endl;
    std::cout << "length sph (V) = " << std::setprecision(32) << distance_sph1 << std::endl;
    std::cout << "length car (V) = " << std::setprecision(32) << distance_car1 << std::endl;
    std::cout << "length geo (A) = " << std::setprecision(32) << distance_geo2 << std::endl;
    std::cout << "length sph (A) = " << std::setprecision(32) << distance_sph2 << std::endl;
    std::cout << "length car (A) = " << std::setprecision(32) << distance_car2 << std::endl;
    std::cout << "length geo (H) = " << std::setprecision(32) << distance_geo3 << std::endl;
    std::cout << "length sph (H) = " << std::setprecision(32) << distance_sph3 << std::endl;
    std::cout << "length car (H) = " << std::setprecision(32) << distance_car3 << std::endl;
    std::cout << "length geo (M) = " << std::setprecision(32) << distance_geo4 << std::endl;
    std::cout << "length sph (M) = " << std::setprecision(32) << distance_sph4 << std::endl;
    std::cout << "length car (M) = " << std::setprecision(32) << distance_car4 << std::endl;
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
    for ( double y = lat_min ; y <= lat_max ; y += step )
    {
        for ( double x = lon_min ; x <= lon_max ; x += step )    
        {
            point_info pi = points[i++];

            bool is_geo_right = ::sin(pi.azi_p-pi.azi_s) >= 0;
            //bool is_geo_bck_right = ::sin(pi.azi_p_bck-pi.azi_s_bck) <= 0;
            bool is_sph_right = pi.sph_s < 0;

            float r = 0.25f, g = 0.25f, b = 0.25f;

            r += is_geo_right ? 0.5f : 0.0f;
            g += is_sph_right ? 0.5f : 0.0f;
            //b += is_geo_bck_right ? 0.5f : 0.0f;

            glColor3f(r, g, b);
            
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

    glColor3f(1, 1, 1);
    glBegin(GL_LINE_STRIP);
    for ( size_t i = 0 ; i < path_geo.size() ; ++i )
    {
        point_info p = path_geo[i];

        double xt = (p.lon - lon_min) * scale_x + disp_x_min;
        double yt = (p.lat - lat_min) * scale_y + disp_y_min;
        
        glVertex3f(xt, yt, 0);
    }
    glEnd();

    glColor3f(0, 0, 1);
    glBegin(GL_LINE_STRIP);
    for ( size_t i = 0 ; i < path_sph.size() ; ++i )
    {
        point_info p = path_sph[i];

        double xt = (p.lon - lon_min) * scale_x + disp_x_min;
        double yt = (p.lat - lat_min) * scale_y + disp_y_min;

        glVertex3f(xt, yt, 0);
    }
    glEnd();

    glColor3f(1, 0, 0);
    glBegin(GL_LINE_STRIP);
    for ( size_t i = 0 ; i < path_car.size() ; ++i )
    {
        point_info p = path_car[i];

        double xt = (p.lon - lon_min) * scale_x + disp_x_min;
        double yt = (p.lat - lat_min) * scale_y + disp_y_min;

        glVertex3f(xt, yt, 0);
    }
    glEnd();

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

void mouse(int /*button*/, int /*state*/, int /*x*/, int /*y*/)
{
}

void keyboard(unsigned char /*key*/, int /*x*/, int /*y*/)
{
}

int main(int argc, char **argv)
{
    std::cout << "filling points" << std::endl;
    fill_points();
    std::cout << "measuring paths" << std::endl;
    measure_paths();
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
