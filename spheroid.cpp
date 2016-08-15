#include <GL/glut.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/arithmetic/cross_product.hpp>

#include <boost/geometry/formulas/gnomonic_intersection.hpp>
#include <boost/geometry/formulas/sjoberg_intersection.hpp>
#include <boost/geometry/formulas/vincenty_direct.hpp>
#include <boost/geometry/formulas/vincenty_inverse.hpp>
#include <boost/geometry/formulas/andoyer_inverse.hpp>
#include <boost/geometry/formulas/thomas_direct.hpp>
#include <boost/geometry/formulas/thomas_inverse.hpp>
#include <boost/geometry/formulas/geographic.hpp>
#include <boost/geometry/formulas/spherical.hpp>

#include <iostream>
#include <vector>

namespace bg = boost::geometry;
namespace bgm = bg::model;
namespace bgd = bg::detail;

typedef bgm::point<double, 2, bg::cs::geographic<bg::degree> > point_geo;
typedef bgm::point<double, 2, bg::cs::spherical_equatorial<bg::degree> > point_sph;
typedef bgm::point<double, 3, bg::cs::cartesian> point_3d;
typedef bg::srs::spheroid<double> spheroid;

double a = 1; //6378137.0;
double b = 0.75; //6356752.314245;
spheroid sph(a, b);

double f_earth = 1.0 / 298.257223563;
double b_earth = a - f_earth * a;

double pi = bg::math::pi<double>();
double d2r = bg::math::d2r<double>();
double r2d = bg::math::r2d<double>();

double flattening(double a, double b) { return (a - b) / a; }

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

color operator+(color const& l, color const& r) { return color(l.r+r.r, l.g+r.g, l.b+r.b, l.a+r.a); }
color operator-(color const& l, color const& r) { return color(l.r-r.r, l.g-r.g, l.b-r.b, l.a-r.a); }
color operator*(color const& l, float f) { return color(l.r*f, l.g*f, l.b*f, l.a*f); }

point_3d operator+(point_3d const& p1, point_3d const& p2) { point_3d res = p1; bg::add_point(res, p2); return res; }
point_3d operator-(point_3d const& p1, point_3d const& p2) { point_3d res = p1; bg::subtract_point(res, p2); return res; }
point_3d operator*(point_3d const& p1, double v) { point_3d res = p1; bg::multiply_value(res, v); return res; }
point_3d operator/(point_3d const& p1, double v) { point_3d res = p1; bg::divide_value(res, v); return res; }

void convert(point_geo const& p, point_3d & res)
{
    res = bg::formula::geo_to_cart3d<point_3d>(p, sph);
    bg::multiply_value(res, a);
}

void convert(point_sph const& p, point_3d & res)
{
    res = bg::formula::sph_to_cart3d<point_3d>(p);
    bg::multiply_value(res, a);
}

void convert(point_3d const& p, point_sph & res)
{
    double r = sqrt(bg::dot_product(p, p));
    bg::set_from_radian<0>(res, atan2(bg::get<1>(p), bg::get<0>(p)));
    bg::set_from_radian<1>(res, asin(bg::get<2>(p) / r));
}

void convert(point_geo const& p, point_sph & res)
{
    bg::set<0>(res, bg::get<0>(p));
    bg::set<1>(res, bg::get<1>(p));
}

void convert(point_sph const& p, point_geo & res)
{
    bg::set<0>(res, bg::get<0>(p));
    bg::set<1>(res, bg::get<1>(p));
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

enum mapping_type { mapping_geodetic, mapping_geocentric, mapping_reduced };

void convert(point_geo const& p, point_sph & res, mapping_type mapping)
{
    bg::set<0>(res, bg::get<0>(p));

    if ( mapping == mapping_geodetic )
    {
        bg::set<1>(res, bg::get<1>(p));
    }
    else
    {
        double b_a = bg::get_radius<2>(sph) / bg::get_radius<0>(sph);
        
        if ( mapping == mapping_geocentric )
            bg::set_from_radian<1>(res, atan(b_a * b_a * tan(bg::get_as_radian<1>(p))));
        else
            bg::set_from_radian<1>(res, atan(b_a * tan(bg::get_as_radian<1>(p))));
    }
}

void convert(point_sph const& p, point_geo & res, mapping_type mapping)
{
    bg::set<0>(res, bg::get<0>(p));
    
    if ( mapping == mapping_geodetic )
    {
        bg::set<1>(res, bg::get<1>(p));
    }
    else
    {
        double a_b = bg::get_radius<0>(sph) / bg::get_radius<2>(sph);

        if ( mapping == mapping_geocentric )
            bg::set_from_radian<1>(res, atan(a_b * a_b * tan(bg::get_as_radian<1>(p))));
        else
            bg::set_from_radian<1>(res, atan(a_b * tan(bg::get_as_radian<1>(p))));
    }
}

point_geo projected_to_equator(point_geo const& p)
{
    point_geo res = p;
    bg::set<1>(res, 0);
    return res;
}

point_3d projected_to_xy_vert(point_geo const& p)
{
    point_3d res;
    ::convert(p, res);
    bg::set<2>(res, 0);
    return res;
}

point_3d projected_to_xy_geod(point_geo const& p)
{
    point_3d p3d;
    ::convert(p, p3d);
    return bg::formula::projected_to_xy(p3d, sph);
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

        //double eq = (x*x+y*y)/(a*a) + z*z/(b*b);
        //assert(bg::math::equals(eq, 1.0));

        glVertex3f(x, y, z);
    }
    glEnd();
}

void draw_parallel(double lat, double step = pi/32)
{
    lat = atan(b/a*tan(lat)); // geodesic lat to geocentric lat

    glBegin(GL_LINE_STRIP);
    for (double lon = 0 ; lon < 2*pi+step ; lon += step)
    {
        double cos_lat = cos(lat);
        double x = a * cos_lat * cos(lon);
        double y = a * cos_lat * sin(lon);
        double z = b * sin(lat);

        //double eq = (x*x+y*y)/(a*a) + z*z/(b*b);
        //assert(bg::math::equals(eq, 1.0));

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

void draw_point(point_3d const& p, double r = 0.01)
{
    draw_sphere(bg::get<0>(p), bg::get<1>(p), bg::get<2>(p), r);
}

template <typename P>
void draw_point(P const& p)
{
    point_3d p3d;
    ::convert(p, p3d);
    draw_point(p3d);
}

void draw_line(point_3d const& p1, point_3d const& p2)
{
    glBegin(GL_LINES);
    glVertex3f(bg::get<0>(p1), bg::get<1>(p1), bg::get<2>(p1));
    glVertex3f(bg::get<0>(p2), bg::get<1>(p2), bg::get<2>(p2));
    glEnd();
}

void draw_triangle(point_3d const& p1, point_3d const& p2, point_3d const& p3)
{
    glBegin(GL_TRIANGLES);
    glVertex3f(bg::get<0>(p1), bg::get<1>(p1), bg::get<2>(p1));
    glVertex3f(bg::get<0>(p2), bg::get<1>(p2), bg::get<2>(p2));
    glVertex3f(bg::get<0>(p3), bg::get<1>(p3), bg::get<2>(p3));
    glEnd();
}

template <typename P1, typename P2>
void draw_line(P1 const& p1, P2 const& p2)
{
    point_3d p3d1, p3d2;
    ::convert(p1, p3d1);
    ::convert(p2, p3d2);
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
    point_3d p_xy = projected_to_xy_vert(p);
    point_3d p_xy_geod = projected_to_xy_geod(p);

    glColor4f(c.r, c.g, c.b, c.a);
    draw_point(p);
    draw_point(p_xy_geod);

    draw_line(point_3d(0, 0, 0), p_xy);
    draw_line(p_xy_geod, p);
}

void normalize(point_3d & p)
{
    double l_sqr = bg::dot_product(p, p);
    if ( !bg::math::equals(l_sqr, 0) )
    {
        bg::divide_value(p, sqrt(l_sqr));
    }
}

double degree(double a, double m, double s)
{
    double sign = a >= 0 ? 1 : -1;
    return sign * ( fabs(a) + fabs(m)/60 + fabs(s)/3600 );
}

std::pair<double, double> andoyer_inverse(point_geo const& p1, point_geo const& p2, spheroid const& sph)
{
    double lon1 = bg::get_as_radian<0>(p1);
    double lat1 = bg::get_as_radian<1>(p1);
    double lon2 = bg::get_as_radian<0>(p2);
    double lat2 = bg::get_as_radian<1>(p2);

    typedef bg::formula::andoyer_inverse<double, true, true> andoyer_t;
    andoyer_t::result_type ai = andoyer_t::apply(lon1, lat1, lon2, lat2, sph);
    return std::make_pair(ai.distance, ai.azimuth);
}

std::pair<double, double> thomas_inverse(point_geo const& p1, point_geo const& p2, spheroid const& sph)
{
    double lon1 = bg::get_as_radian<0>(p1);
    double lat1 = bg::get_as_radian<1>(p1);
    double lon2 = bg::get_as_radian<0>(p2);
    double lat2 = bg::get_as_radian<1>(p2);

    typedef bg::formula::thomas_inverse<double, true, true> thomas_t;
    thomas_t::result_type ti = thomas_t::apply(lon1, lat1, lon2, lat2, sph);
    return std::make_pair(ti.distance, ti.azimuth);
}

double bearing(double rad1, double rad2)
{
    double d = rad1 - rad2;
    while ( d > bg::math::pi<double>() )
        d -= 2 * bg::math::pi<double>();
    while ( d < -bg::math::pi<double>() )
        d += 2 * bg::math::pi<double>();
    return bg::math::abs(d);
}

struct scene_data
{
    enum mode_type { mode_navigation = 0, mode_segment, mode_intersection } mode;
    bool enable_experimental;
    bool enable_mapping_geodetic;
    bool enable_mapping_geocentric;
    bool enable_mapping_reduced;
    bool enable_great_ellipse;
    bool enable_vincenty;
    bool enable_andoyer;
    bool enable_thomas;

    std::vector<point_3d> curve_experimental;
    std::vector<point_3d> curve_mapped_geodetic;
    std::vector<point_3d> curve_mapped_geocentric;
    std::vector<point_3d> curve_mapped_reduced;
    std::vector<point_3d> curve_great_ellipse;
    std::vector<point_3d> curve_vincenty;
    std::vector<point_3d> curve_andoyer;
    std::vector<point_3d> curve_thomas;

    point_3d p1_s; // cartesian point on surface
    point_3d p2_s; // cartesian point on surface
    point_3d p1_xy; // corresponding cartesian point on XY
    point_3d p2_xy; // corresponding cartesian point on XY
    point_3d v_s; // vector between the surface points
    point_3d v_xy; // vector between the XY points

    point_3d q1_s; // cartesian point on surface
    point_3d q2_s; // cartesian point on surface
    point_3d q1_xy; // corresponding cartesian point on XY
    point_3d q2_xy; // corresponding cartesian point on XY
    point_3d w_s; // vector between the surface points
    point_3d w_xy; // vector between the XY points

    point_3d p_east;
    point_3d p_north;
    point_3d q_east;
    point_3d q_north;

    point_3d v_azimuth_vincenty;
    point_3d v_azimuth_andoyer;
    point_3d v_azimuth_thomas;
    point_3d v_azimuth_great_ellipse;

    std::vector<point_3d> curve2_experimental;
    std::vector<point_3d> curve2_vincenty;
    
    point_geo i_experimental;
    point_geo i_vincenty_gnomonic;
    point_geo i_vincenty_sjoberg;

    point_3d i_experimental_s;
    point_3d i_vincenty_gnomonic_s;
    point_3d i_vincenty_sjoberg_s;

    double azimuth_vincenty;
    double azimuth_andoyer;
    double azimuth_thomas;
    double azimuth_great_ellipse;
    double azimuth_mapping_geodetic;
    double azimuth_mapping_geocentric;
    double azimuth_mapping_reduced;
    double azimuth_experimental;

    scene_data()
    {
        mode = mode_navigation;
        enable_experimental = true;
        enable_mapping_geodetic = true;
        enable_mapping_geocentric = true;
        enable_mapping_reduced = true;
        enable_great_ellipse = true;
        enable_vincenty = true;
        enable_andoyer = true;
        enable_thomas = true;

        clear();
    }

    void clear()
    {
        azimuth_vincenty = 0;
        azimuth_andoyer = 0;
        azimuth_thomas = 0;
        azimuth_great_ellipse = 0;
        azimuth_mapping_geodetic = 0;
        azimuth_mapping_geocentric = 0;
        azimuth_mapping_reduced = 0;
        azimuth_experimental = 0;

        curve_experimental.clear();
        curve_mapped_geodetic.clear();
        curve_mapped_geocentric.clear();
        curve_mapped_reduced.clear();
        curve_great_ellipse.clear();
        curve_vincenty.clear();
        curve_andoyer.clear();
        curve_thomas.clear();

        curve2_experimental.clear();
        curve2_vincenty.clear();
    }

    void recalculate_loc(point_geo const& p1, point_geo const& p2,
                         point_geo const& q1, point_geo const& q2)
    {
        // cartesian points
        bg::formula::geo_to_cart3d(p1, p1_s, p_north, p_east, sph);
        p2_s = bg::formula::geo_to_cart3d<point_3d>(p2, sph);
        p1_xy = bg::formula::projected_to_xy(p1_s, sph);
        p2_xy = bg::formula::projected_to_xy(p2_s, sph);
        // cartesian vectors
        v_s = p2_s - p1_s;
        v_xy = p2_xy - p1_xy;

        if (mode == mode_intersection)
        {
            // cartesian points
            bg::formula::geo_to_cart3d(q1, q1_s, q_north, q_east, sph);
            q2_s = bg::formula::geo_to_cart3d<point_3d>(q2, sph);
            q1_xy = bg::formula::projected_to_xy(q1_s, sph);
            q2_xy = bg::formula::projected_to_xy(q2_s, sph);
            // cartesian vectors
            w_s = q2_s - q1_s;
            w_xy = q2_xy - q1_xy;
        }
    }

    point_3d make_p_azimuth(double azimuth)
    {
        return p_north * cos(azimuth) + p_east * sin(azimuth);
    }

    void recalculate(point_geo const& p1, point_geo const& p2,
                     point_geo const& q1, point_geo const& q2)
    {
        clear();

        recalculate_loc(p1, p2, q1, q2);

        int count = 50;
        double f_step = 1.0 / count;

        if (mode == mode_navigation || mode == mode_segment)
        {
            double dist_step = bg::formula::vincenty_inverse<double, true, false>::apply(
                bg::get<0>(p1) * d2r, bg::get<1>(p1) * d2r, bg::get<0>(p2) * d2r, bg::get<1>(p2) * d2r, sph
            ).distance / count;
        
            if (mode == mode_navigation)
            {
                if (enable_experimental)
                {
                    curve_experimental.clear();
                    fill_navigation_curve< experimental_inverse<double> >(p1, p2, dist_step, curve_experimental, azimuth_experimental);
                }

                if (enable_mapping_geodetic)
                {
                    fill_navigation_curve< great_circle_inverse<double, mapping_geodetic> >(p1, p2, dist_step, curve_mapped_geodetic, azimuth_mapping_geodetic);
                }

                if (enable_mapping_geocentric)
                {
                    fill_navigation_curve< great_circle_inverse<double, mapping_geocentric> >(p1, p2, dist_step, curve_mapped_geocentric, azimuth_mapping_geocentric);
                }

                if (enable_mapping_reduced)
                {
                    fill_navigation_curve< great_circle_inverse<double, mapping_reduced> >(p1, p2, dist_step, curve_mapped_reduced, azimuth_mapping_reduced);
                }

                if (enable_great_ellipse)
                {
                    fill_navigation_curve< great_ellipse_inverse<double> >(p1, p2, dist_step, curve_great_ellipse, azimuth_great_ellipse);
                    v_azimuth_great_ellipse = make_p_azimuth(azimuth_great_ellipse);
                }

                if (enable_vincenty)
                {
                    fill_navigation_curve< bg::formula::vincenty_inverse<double, true, true> >(p1, p2, dist_step, curve_vincenty, azimuth_vincenty);
                    v_azimuth_vincenty = make_p_azimuth(azimuth_vincenty);
                }

                if (enable_andoyer)
                {
                    fill_navigation_curve< bg::formula::andoyer_inverse<double, true, true> >(p1, p2, dist_step, curve_andoyer, azimuth_andoyer);
                    v_azimuth_andoyer = make_p_azimuth(azimuth_andoyer);
                }

                if (enable_thomas)
                {
                    fill_navigation_curve< bg::formula::thomas_inverse<double, true, true> >(p1, p2, dist_step, curve_thomas, azimuth_thomas);
                    v_azimuth_thomas = make_p_azimuth(azimuth_thomas);
                }
            }
            else if (mode == mode_segment)
            {
                double lon1 = bg::get_as_radian<0>(p1);
                double lat1 = bg::get_as_radian<1>(p1);
                double lon2 = bg::get_as_radian<0>(p2);
                double lat2 = bg::get_as_radian<1>(p2);

                double azimuth_vincenty = 0;
                double azimuth_thomas = 0;

                if (enable_vincenty)
                {
                    azimuth_vincenty = bg::formula::vincenty_inverse<double, false, true>::apply(
                        lon1, lat1, lon2, lat2, sph).azimuth;
                }

                if (enable_thomas)
                {
                    azimuth_thomas = bg::formula::thomas_inverse<double, false, true>::apply(
                        lon1, lat1, lon2, lat2, sph).azimuth;
                }


                double f = 0;
                double d = 0;
                for (int i = 0; i <= count; ++i, f += f_step, d += dist_step)
                {
                    point_3d ps = p1_s + v_s * f;

                    // experimental method
                    if (enable_experimental)
                    {
                        point_3d o = p1_xy + v_xy * 0.5;
                        point_3d d = ps - o;
                        point_3d p_curve = o, dummy;
                        bg::formula::projected_to_surface(o, d, p_curve, dummy, sph);
                        curve_experimental.push_back(p_curve);
                    }

                    if (enable_great_ellipse)
                    {
                        point_3d p_curve = bg::formula::projected_to_surface(ps, sph);
                        curve_great_ellipse.push_back(p_curve);
                    }

                    if (enable_vincenty)
                    {
                        typedef bg::formula::vincenty_direct<double> direct_t;
                        direct_t::result_type res = direct_t::apply(lon1, lat1, d, azimuth_vincenty, sph);
                        point_geo p;
                        bg::set_from_radian<0>(p, res.lon2);
                        bg::set_from_radian<1>(p, res.lat2);
                        point_3d p_curve = pcast<point_3d>(p);
                        curve_vincenty.push_back(p_curve);
                    }

                    if (enable_thomas)
                    {
                        typedef bg::formula::thomas_direct<double> direct_t;
                        direct_t::result_type res = direct_t::apply(lon1, lat1, d, azimuth_thomas, sph);
                        point_geo p;
                        bg::set_from_radian<0>(p, res.lon2);
                        bg::set_from_radian<1>(p, res.lat2);
                        point_3d p_curve = pcast<point_3d>(p);
                        curve_thomas.push_back(p_curve);
                    }

                    // mapped to sphere
                    if (enable_mapping_geodetic)
                    {
                        push_to_mapped(p1, p2, f, curve_mapped_geodetic, mapping_geodetic);
                    }

                    if (enable_mapping_geocentric)
                    {
                        push_to_mapped(p1, p2, f, curve_mapped_geocentric, mapping_geocentric);
                    }

                    if (enable_mapping_reduced)
                    {
                        push_to_mapped(p1, p2, f, curve_mapped_reduced, mapping_reduced);
                    }
                }
            }
        }
        else if (mode == mode_intersection)
        {
            double plon1 = bg::get_as_radian<0>(p1);
            double plat1 = bg::get_as_radian<1>(p1);
            double plon2 = bg::get_as_radian<0>(p2);
            double plat2 = bg::get_as_radian<1>(p2);
            double qlon1 = bg::get_as_radian<0>(q1);
            double qlat1 = bg::get_as_radian<1>(q1);
            double qlon2 = bg::get_as_radian<0>(q2);
            double qlat2 = bg::get_as_radian<1>(q2);

            typedef bg::formula::vincenty_inverse<double, true, true> inverse_t;
            inverse_t::result_type ires = inverse_t::apply(plon1, plat1, plon2, plat2, sph);
            double dist_step1 = ires.distance / count;
            double azimuth1 = ires.azimuth;

            ires = inverse_t::apply(qlon1, qlat1, qlon2, qlat2, sph);
            double dist_step2 = ires.distance / count;
            double azimuth2 = ires.azimuth;

            double f = 0;
            double d1 = 0;
            double d2 = 0;
            for (int i = 0; i <= count; ++i, f += f_step, d1 += dist_step1, d2 += dist_step2)
            {
                typedef bg::formula::vincenty_direct<double> direct_t;
                direct_t::result_type res = direct_t::apply(plon1, plat1, d1, azimuth1, sph);
                point_geo p;
                bg::set_from_radian<0>(p, res.lon2);
                bg::set_from_radian<1>(p, res.lat2);
                curve_vincenty.push_back(pcast<point_3d>(p));

                res = direct_t::apply(qlon1, qlat1, d2, azimuth2, sph);
                bg::set_from_radian<0>(p, res.lon2);
                bg::set_from_radian<1>(p, res.lat2);
                curve2_vincenty.push_back(pcast<point_3d>(p));

                {
                    point_3d ps = p1_s + v_s * f;
                    point_3d o = p1_xy + v_xy * 0.5;
                    point_3d d = ps - o;
                    point_3d p_curve = o, dummy;
                    bg::formula::projected_to_surface(o, d, p_curve, dummy, sph);
                    curve_experimental.push_back(p_curve);
                }

                {
                    point_3d qs = q1_s + w_s * f;
                    point_3d o = q1_xy + w_xy * 0.5;
                    point_3d d = qs - o;
                    point_3d p_curve = o, dummy;
                    bg::formula::projected_to_surface(o, d, p_curve, dummy, sph);
                    curve2_experimental.push_back(p_curve);
                }
            }

            bg::formula::elliptic_intersection(p1_s, p2_s, q1_s, q2_s, i_experimental_s, sph);
            i_experimental = bg::formula::cart3d_to_geo<point_geo>(i_experimental_s, sph);

            double i_lon = 0, i_lat = 0;
            bg::formula::gnomonic_intersection<double, bg::formula::vincenty_inverse, bg::formula::vincenty_direct>
                ::apply(plon1, plat1, plon2, plat2, qlon1, qlat1, qlon2, qlat2, i_lon, i_lat, sph);
            bg::set_from_radian<0>(i_vincenty_gnomonic, i_lon);
            bg::set_from_radian<1>(i_vincenty_gnomonic, i_lat);
            i_vincenty_gnomonic_s = bg::formula::geo_to_cart3d<point_3d>(i_vincenty_gnomonic, sph);

            bg::formula::sjoberg_intersection<double, bg::formula::vincenty_inverse, 4>
                ::apply(plon1, plat1, plon2, plat2,
                        qlon1, qlat1, qlon2, qlat2,
                        i_lon, i_lat,
                        sph);
            bg::set_from_radian<0>(i_vincenty_sjoberg, i_lon);
            bg::set_from_radian<1>(i_vincenty_sjoberg, i_lat);
            i_vincenty_sjoberg_s = bg::formula::geo_to_cart3d<point_3d>(i_vincenty_sjoberg, sph);
        }
    }

    template<typename T, mapping_type Mapping>
    struct great_circle_inverse
    {
        struct result_type
        {
            T azimuth;
        };

        template <typename Spheroid>
        static inline result_type apply(T const& lon1, T const& lat1, T const& lon2, T const& lat2, Spheroid const& spheroid)
        {
            point_geo p1(lon1 * r2d, lat1 * r2d);
            point_geo p2(lon2 * r2d, lat2 * r2d);

            point_sph p1_m, p2_m;
            ::convert(p1, p1_m, Mapping);
            ::convert(p2, p2_m, Mapping);

            result_type res;
            res.azimuth = bg::detail::azimuth<T>(p1_m, p2_m);
            return res;
        }
    };

    template <typename T>
    struct great_ellipse_inverse
    {
        struct result_type
        {
            T azimuth;
        };

        template <typename Spheroid>
        static inline result_type apply(T const& lon1, T const& lat1, T const& lon2, T const& lat2, Spheroid const& spheroid)
        {
            point_geo p1(lon1 * r2d, lat1 * r2d);
            point_geo p2(lon2 * r2d, lat2 * r2d);

            point_3d p1_s, north, east;
            bg::formula::geo_to_cart3d(p1, p1_s, north, east, sph);
            point_3d p2_s = bg::formula::geo_to_cart3d<point_3d>(p2, sph);

            point_3d n = bg::cross_product(p1_s, p2_s);
            point_3d aziv = bg::cross_product(n, p1_s);
            normalize(aziv);
            
            result_type res;
            res.azimuth = atan2(bg::dot_product(aziv, east), bg::dot_product(aziv, north));
            return res;
        }
    };

    template <typename T>
    struct experimental_inverse
    {
        struct result_type
        {
            T azimuth;
        };

        template <typename Spheroid>
        static inline result_type apply(T const& lon1, T const& lat1, T const& lon2, T const& lat2, Spheroid const& spheroid)
        {
            point_geo p1(lon1 * r2d, lat1 * r2d);
            point_geo p2(lon2 * r2d, lat2 * r2d);

            point_3d p1_s, north, east;
            bg::formula::geo_to_cart3d(p1, p1_s, north, east, sph);
            
            point_3d p2_s = bg::formula::geo_to_cart3d<point_3d>(p2, sph);
            point_3d p1_xy = bg::formula::projected_to_xy(p1_s, sph);
            point_3d p2_xy = bg::formula::projected_to_xy(p2_s, sph);

            //point_3d v12 = p2v - p1v;
            point_3d v12_xy = p2_xy - p1_xy;

            point_3d origin_xy = p1_xy + v12_xy * 0.5;
            
            point_3d op1 = p1_s - origin_xy;
            point_3d op2 = p2_s - origin_xy;

            point_3d n = bg::cross_product(op1, op2);
            point_3d aziv = bg::cross_product(n, op1);
            normalize(aziv);
            
            result_type res;
            res.azimuth = atan2(bg::dot_product(aziv, east), bg::dot_product(aziv, north));
            return res;
        }
    };

    template <typename Inverse>
    static inline void fill_navigation_curve(point_geo const& p1, point_geo const& p2, double dist_step, std::vector<point_3d> & curve, double & azimuth)
    {
        int max_count = 100;

        double lon = bg::get_as_radian<0>(p1);
        double lat = bg::get_as_radian<1>(p1);
        double azi = 0;

        curve.push_back(pcast<point_3d>(p1));

        {
            typename Inverse::result_type inv = Inverse::apply(
                        lon,
                        lat,
                        bg::get_as_radian<0>(p2),
                        bg::get_as_radian<1>(p2),
                        sph);

            azimuth = azi = inv.azimuth;
        }

        for ( int i = 0 ; i <= max_count ; ++i )
        {
            // calculate the position of p2 for given d and azimuth
            typedef bg::formula::vincenty_direct<double> direct_t;
            direct_t::result_type vd = direct_t::apply(
                        lon,
                        lat,
                        dist_step,
                        azi,
                        sph);

            lon = vd.lon2;
            lat = vd.lat2;

            typename Inverse::result_type inv = Inverse::apply(
                        lon,
                        lat,
                        bg::get_as_radian<0>(p2),
                        bg::get_as_radian<1>(p2),
                        sph);

            double ba = bearing(azi, inv.azimuth);
            if (ba > bg::math::pi<double>() / 2)
            {
                curve.push_back(pcast<point_3d>(p2));
                return;
            }
            else
            {
                azi = inv.azimuth;

                point_geo p;
                bg::set_from_radian<0>(p, lon);
                bg::set_from_radian<1>(p, lat);

                curve.push_back(pcast<point_3d>(p));
            }
        }
    }

    void push_to_mapped(point_geo const& p1,
                        point_geo const& p2,
                        double const& f,
                        std::vector<point_3d> & curve,
                        mapping_type mapping)
    {
        // geographic mapped to spherical
        point_sph p1_m, p2_m;
        ::convert(p1, p1_m, mapping);
        ::convert(p2, p2_m, mapping);
        
        // spherical mapped to cartesian 3d
        point_3d p1_mc = pcast<point_3d>(p1_m);
        point_3d p2_mc = pcast<point_3d>(p2_m);
        point_3d v_mc = p2_mc - p1_mc;

        // vector attached at 0.0
        point_3d d_mc = p1_mc + v_mc * f;

        //x*x+y*y+z*z = a*a
        double dx = bg::get<0>(d_mc);
        double dy = bg::get<1>(d_mc);
        double dz = bg::get<2>(d_mc);

        double param_a = dx*dx+dy*dy+dz*dz;
        double param_c = -a*a;

        double delta = -4*param_a*param_c;
        double t = delta >= 0 ?
                   sqrt(delta) / (2*param_a) :
                   0.0;

        // 3d cartesian point on a surface of a sphere
        point_3d p_curve_mc = d_mc * t;
        // spherical point
        point_sph p_curve_m = pcast<point_sph>(p_curve_mc);
        // geographical point
        point_geo p_curve;
        ::convert(p_curve_m, p_curve, mapping);
        // 3d cartesian point on a surface of a spheroid
        point_3d p_curve_3d = pcast<point_3d>(p_curve);

        curve.push_back(p_curve_3d);
    }

    void cycle_mode()
    {
        mode = mode_type((mode + 1) % 3);
    }

    void print_settings() const
    {
        std::cout << "SETTINGS\n"
                  << "mode:               " << mode_str() << '\n'
                  << "experimental:       " << flag_str(enable_experimental) << '\n'
                  << "mapping_geodetic:   " << flag_str(enable_mapping_geodetic) << '\n'
                  << "mapping_geocentric: " << flag_str(enable_mapping_geocentric) << '\n'
                  << "mapping_reduced:    " << flag_str(enable_mapping_reduced) << '\n'
                  << "great_ellipse:      " << flag_str(enable_great_ellipse) << '\n'
                  << "vincenty:           " << flag_str(enable_vincenty) << '\n'
                  << "andoyer:            " << flag_str(enable_andoyer) << '\n'
                  << "thomas:             " << flag_str(enable_thomas) << '\n';
        std::cout.flush();
    }

    const char * flag_str(bool f) const
    {
        return f ? "on" : "off";
    }

    const char * mode_str() const
    {
        return mode == mode_navigation ? "navigation" :
               mode == mode_segment ? "segment" :
               mode == mode_intersection ? "intersection" :
               "unknown";
    }

    void print_curve_lengths() const
    {
        std::cout << std::setprecision(16);
        std::cout << "LENGTHS\n"
                  << "vincenty:           " << curve_length(curve_vincenty) << '\n'
                  << "thomas:             " << curve_length(curve_thomas) << '\n'
                  << "andoyer:            " << curve_length(curve_andoyer) << '\n'
                  << "experimental:       " << curve_length(curve_experimental) << '\n'
                  << "great_ellipse:      " << curve_length(curve_great_ellipse) << '\n'
                  << "mapping_geodetic:   " << curve_length(curve_mapped_geodetic) << '\n'
                  << "mapping_geocentric: " << curve_length(curve_mapped_geocentric) << '\n'
                  << "mapping_reduced:    " << curve_length(curve_mapped_reduced) << '\n';
        std::cout.flush();
    }

    static double curve_length(std::vector<point_3d> const& curve)
    {
        double result = 0.0;

        size_t count = curve.size();
        if ( count == 0 )
            return result;

        for ( size_t i = 1 ; i < count ; ++i )
        {
            point_3d v_seg = curve[i] - curve[i-1];

            result += sqrt(bg::dot_product(v_seg, v_seg));
        }

        return result;
    }

} data;

void draw_curve(std::vector<point_3d> const& curve, color const& color_first, color const& color_last)
{
    size_t count = curve.size();

    if ( count == 0 )
        return;

    double f = 0;
    double f_step = 1.0 / count;
    for ( size_t i = 1 ; i < count ; ++i, f += f_step )
    {
        color col = color_first + (color_last - color_first) * f;
        glColor3f(col.r, col.g, col.b);
        draw_line(curve[i-1], curve[i]);
        
        /*if (i == 1)
            draw_point(curve[0], 0.005);
        draw_point(curve[i], 0.005);*/
    }
}

point_geo p1(-30, -30);
point_geo p2(30, 30);
point_geo q1(-30, 30);
point_geo q2(30, -30);

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

    glLineWidth(1);

    draw_model();

    glColor3f(1, 0.5, 0);
    draw_parallel(bg::get<1>(p1) * d2r);
    glColor3f(1, 1, 0);
    draw_parallel(bg::get<1>(p2) * d2r);

    glLineWidth(2);

    draw_point_adv(p1, color(1, 0.5, 0));
    draw_point_adv(p2, color(1, 1, 0));

    // loc
    glColor3f(1, 0, 0);
    draw_line(data.p1_s, data.p1_s + data.p_north*0.2);
    glColor3f(0, 1, 0);
    draw_line(data.p1_s, data.p1_s + data.p_east*0.2);

    if (data.mode == scene_data::mode_navigation || data.mode == scene_data::mode_segment)
    {
        glLineWidth(3);

        if (data.enable_experimental) // orange -> yellow
            draw_curve(data.curve_experimental, color(1, 0.5, 0), color(1, 1, 0));

        if (data.enable_mapping_geodetic) // red
            draw_curve(data.curve_mapped_geodetic, color(0.5, 0, 0), color(1, 0, 0));
        if (data.enable_mapping_geocentric) // violet
            draw_curve(data.curve_mapped_geocentric, color(0.25, 0, 0.5), color(0.5, 0, 1));
        if (data.enable_mapping_reduced) // blue
            draw_curve(data.curve_mapped_reduced, color(0, 0.25, 0.5), color(0, 0.5, 1));

        if (data.enable_great_ellipse) // green
        {
            draw_curve(data.curve_great_ellipse, color(0, 0.5, 0), color(0, 1, 0));
        }

        if (data.enable_vincenty) // gray->white
        {
            draw_curve(data.curve_vincenty, color(0.75, 0.75, 0.75), color(1, 1, 1));
        }
        if (data.enable_andoyer) // magenta
        {
            draw_curve(data.curve_andoyer, color(0.75, 0, 0.75), color(1, 0, 1));
        }
        if (data.enable_thomas) // cyan
        {
            draw_curve(data.curve_thomas, color(0, 0.75, 0.75), color(0, 1, 1));
        }
    }

    if (data.mode == scene_data::mode_intersection)
    {
        glLineWidth(2);

        // loc
        glColor3f(1, 0, 0);
        draw_line(data.q1_s, data.q1_s + data.q_north*0.2);
        glColor3f(0, 1, 0);
        draw_line(data.q1_s, data.q1_s + data.q_east*0.2);

        glLineWidth(1);

        glColor3f(0, 0.5, 1);
        draw_parallel(bg::get<1>(q1) * d2r);
        glColor3f(0, 1, 1);
        draw_parallel(bg::get<1>(q2) * d2r);

        glLineWidth(2);

        draw_point_adv(q1, color(0, 0.5, 1));
        draw_point_adv(q2, color(0, 1, 1));

        glLineWidth(3);

        draw_curve(data.curve_experimental, color(1, 0.5, 0), color(1, 1, 0));
        draw_curve(data.curve2_experimental, color(0, 0.5, 1), color(0, 1, 1));
        draw_curve(data.curve_vincenty, color(1, 1, 1), color(1, 1, 1));
        draw_curve(data.curve2_vincenty, color(1, 1, 1), color(1, 1, 1));

        glColor3f(1, 1, 0);
        draw_point(data.i_experimental_s);
        glColor3f(1, 1, 1);
        draw_point(data.i_vincenty_gnomonic_s);
        glColor3f(1, 0, 1);
        draw_point(data.i_vincenty_sjoberg_s);
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

void print_help()
{
    std::cout << "UI\n"
              << "Mouse - navigation" << '\n'
              << "WSAD  - move p1" << '\n'
              << "IKJL  - move p2" << '\n'
              << "?     - display help" << '\n'
              << ",     - decrease b, increase flattening" << '\n'
              << ".     - increase b, decrease flattening" << '\n'
              << "m     - mode switch" << '\n'
              << "1     - experimental curve on/off" << '\n'
              << "2     - geodetic curve on/off" << '\n'
              << "3     - geocentric curve on/off" << '\n'
              << "4     - reduced curve on/off" << '\n'
              << "5     - great ellipse curve on/off" << '\n'
              << "6     - vincenty curve on/off" << '\n'
              << "7     - andoyer curve on/off" << '\n';
    std::cout.flush();
}

void print_geometry()
{
    std::cout << std::setprecision(16);

    std::cout << "GEOMETRY\n";
    std::cout << "p1:         (" << bg::get<0>(p1) << ", " << bg::get<1>(p1) << ")\n";
    std::cout << "p2:         (" << bg::get<0>(p2) << ", " << bg::get<1>(p2) << ")\n";
    std::cout << "flattening: " << flattening(a, b) << std::endl;
}

void print_distances_and_azimuths()
{
    std::pair<double, double> ai_res = andoyer_inverse(p1, p2, sph);
    std::pair<double, double> ti_res = thomas_inverse(p1, p2, sph);

    std::cout << std::setprecision(16);
    std::cout << "DISTANCES\n"
              << "vincenty:       " << bg::distance(p1, p2, bg::strategy::distance::vincenty<spheroid>(sph)) << '\n'
              << "thomas:         " << ti_res.first << '\n'
              << "andoyer:        " << ai_res.first << '\n'
              << "haversine:      " << bg::distance(p1, p2, bg::strategy::distance::haversine<double>((2*a+b)/3)) << '\n';
              
    std::cout << "AZIMUTHS\n"
              << "vincenty:       " << bg::detail::azimuth<double>(p1, p2, sph) << '\n'
              << "thomas:         " << ti_res.second << '\n'
              << "andoyer:        " << ai_res.second << '\n'
              << "experimental:   " << data.azimuth_experimental << '\n'
              << "great ellipse:  " << data.azimuth_great_ellipse << '\n'
              << "map geocentric: " << data.azimuth_mapping_geocentric << '\n'
              << "map geodetic:   " << data.azimuth_mapping_geodetic << '\n'
              << "map reduced:    " << data.azimuth_mapping_reduced << '\n';
                  
    if (data.mode == scene_data::mode_intersection)
    {
        std::cout << "INTERSECTIONS\n"
                  << "vincenty gnomonic: " << bg::get<0>(data.i_vincenty_gnomonic) << ", " << bg::get<1>(data.i_vincenty_gnomonic) << '\n'
                  << "vincenty sjoberg:  " << bg::get<0>(data.i_vincenty_sjoberg) << ", " << bg::get<1>(data.i_vincenty_sjoberg) << '\n'
                  << "experimental:      " << bg::get<0>(data.i_experimental) << ", " << bg::get<1>(data.i_experimental) << '\n';
    }

    std::cout.flush();
}

void move_lat(point_geo & p, double diff)
{
    double l = bg::get<1>(p) + diff;
    if ( l > 90 ) l = 90;
    if ( l < -90 ) l = -90;
    bg::set<1>(p, l);
}

void move_lon(point_geo & p, double diff)
{
    double l = bg::get<0>(p) + diff;
    if ( l > 180 ) l -= 360;
    if ( l < -180 ) l += 360;
    bg::set<0>(p, l);
}

void keyboard(unsigned char key, int /*x*/, int /*y*/)
{
    static const double b_step = 0.05;

    bool refresh = true;

    if ( key == '/' )
    {
        print_help();
    }
    else if ( key == '.')
    {
        if (b == b_earth)
        {
            b = a;
        }
        else
        {
            double new_b = b + b_step;
            if (new_b > 2 * a - b_step)
                new_b = 2 * a - b_step;

            double f = flattening(a, b);
            double new_f = flattening(a, new_b);
            if (f > f_earth && f_earth > new_f)
                new_b = b_earth;

            b = new_b;
        }

        sph = spheroid(a, b);
    }
    else if ( key == ',' )
    {
        b -= b_step;
        if ( b < b_step )
            b = b_step;
        sph = spheroid(a, b);
    }
    else if ( key == 'm' )
    {
        data.cycle_mode();
    }
    else if ( key == '1' )
    {
        data.enable_experimental = !data.enable_experimental;
    }
    else if ( key == '2' )
    {
        data.enable_mapping_geodetic = !data.enable_mapping_geodetic;
    }
    else if ( key == '3' )
    {
        data.enable_mapping_geocentric = !data.enable_mapping_geocentric;
    }
    else if ( key == '4' )
    {
        data.enable_mapping_reduced = !data.enable_mapping_reduced;
    }
    else if ( key == '5' )
    {
        data.enable_great_ellipse = !data.enable_great_ellipse;
    }
    else if ( key == '6' )
    {
        data.enable_vincenty = !data.enable_vincenty;
    }
    else if ( key == '7' )
    {
        data.enable_andoyer = !data.enable_andoyer;
    }
    else if ( key == '8' )
    {
        data.enable_thomas = !data.enable_thomas;
    }
    // moving of the p1
    else if ( key == 'w' )
        move_lat(p1, 1);
    else if ( key == 's' )
        move_lat(p1, -1);
    else if ( key == 'a' )
        move_lon(p1, -1);
    else if ( key == 'd' )
        move_lon(p1, 1);
    // moving of the q1
    else if (key == 't')
        move_lat(q1, 1);
    else if (key == 'g')
        move_lat(q1, -1);
    else if (key == 'f')
        move_lon(q1, -1);
    else if (key == 'h')
        move_lon(q1, 1);
    // moving of the q2
    else if (key == 'i')
        move_lat(q2, 1);
    else if (key == 'k')
        move_lat(q2, -1);
    else if (key == 'j')
        move_lon(q2, -1);
    else if (key == 'l')
        move_lon(q2, 1);
    // other key
    else
        refresh = false;

    if ( refresh )
    {
        print_geometry();
        data.print_settings();

        data.recalculate(p1, p2, q1, q2);
        data.print_curve_lengths();
        print_distances_and_azimuths();
    }
}

void special_input(int key, int x, int y)
{
    bool refresh = true;

    // moving of the p2
    switch (key)
    {
    case GLUT_KEY_UP:
        move_lat(p2, 1);
        break;
    case GLUT_KEY_DOWN:
        move_lat(p2, -1);
        break;
    case GLUT_KEY_LEFT:
        move_lon(p2, -1);
        break;
    case GLUT_KEY_RIGHT:
        move_lon(p2, 1);
        break;
    default:
        refresh = false;
    }

    if (refresh)
    {
        print_geometry();
        data.print_settings();

        data.recalculate(p1, p2, q1, q2);
        data.print_curve_lengths();
        print_distances_and_azimuths();
    }
}

void idle_fun()
{
    glutPostRedisplay();
}

int main(int argc, char **argv)
{
    //test
    {
        point_geo p1(degree(0, 0, 0), degree(0, 0, 0));
        point_geo p2(degree(0, 0, 0), degree(10, 0, 0));
        std::cout << "TEST DISTANCES\n" << std::setprecision(16)
                  << "vincenty:     " << bg::distance(p1, p2, bg::strategy::distance::vincenty<spheroid>(sph)) << '\n'
                  << "andoyer: " << bg::distance(p1, p2, bg::strategy::distance::andoyer<spheroid>(sph)) << '\n'
                  << "thomas: " << bg::distance(p1, p2, bg::strategy::distance::thomas<spheroid>(sph)) << '\n'
                  << "haversine:    " << bg::distance(p1, p2, bg::strategy::distance::haversine<double>((2*a+b)/3)) << '\n';
    }

    print_help();
    print_geometry();
    data.print_settings();

    data.recalculate(p1, p2, q1, q2);
    data.print_curve_lengths();
    print_distances_and_azimuths();
    
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
    glutSpecialFunc(special_input);

    glutMainLoop();

    return 0;
}
