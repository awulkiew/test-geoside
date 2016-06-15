#include <GL/glut.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/extensions/arithmetic/cross_product.hpp>

#include <boost/geometry/algorithms/detail/vincenty_direct.hpp>
#include <boost/geometry/algorithms/detail/vincenty_inverse.hpp>
#include <boost/geometry/algorithms/detail/andoyer_inverse.hpp>
#include <boost/geometry/algorithms/detail/thomas_inverse.hpp>

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

color operator+(color const& l, color const& r) { return color(l.r+r.r, l.g+r.g, l.b+r.b, l.a+r.a); }
color operator-(color const& l, color const& r) { return color(l.r-r.r, l.g-r.g, l.b-r.b, l.a-r.a); }
color operator*(color const& l, float f) { return color(l.r*f, l.g*f, l.b*f, l.a*f); }

point_3d operator+(point_3d const& p1, point_3d const& p2) { point_3d res = p1; bg::add_point(res, p2); return res; }
point_3d operator-(point_3d const& p1, point_3d const& p2) { point_3d res = p1; bg::subtract_point(res, p2); return res; }
point_3d operator*(point_3d const& p1, double v) { point_3d res = p1; bg::multiply_value(res, v); return res; }
point_3d operator/(point_3d const& p1, double v) { point_3d res = p1; bg::divide_value(res, v); return res; }

void convert(point_geo const& p, point_3d & res)
{
    double lon = bg::get_as_radian<0>(p);
    double lat = bg::get_as_radian<1>(p);
    lat = atan(b/a*tan(lat)); // geodesic lat to geocentric lat
    double cos_lat = cos(lat);
    bg::set<0>(res, a * cos_lat * cos(lon));
    bg::set<1>(res, a * cos_lat * sin(lon));
    bg::set<2>(res, b * sin(lat));
}

void convert(point_sph const& p, point_3d & res)
{
    double lon = bg::get_as_radian<0>(p);
    double lat = bg::get_as_radian<1>(p);
    double cos_lat = cos(lat);
    bg::set<0>(res, a * cos_lat * cos(lon));
    bg::set<1>(res, a * cos_lat * sin(lon));
    bg::set<2>(res, a * sin(lat));
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

void convert(point_3d const& p, point_3d & res)
{
    res = p;
}

template <typename Res, typename P>
Res pcast(P const& p)
{
    Res res;
    ::convert(p, res);
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
    ::convert(p, res);
    bg::set<2>(res, 0);
    return res;
}

point_3d projected_to_xy_geod(point_geo const& p)
{
    double r = 0;
    double lon = bg::get_as_radian<0>(p);
    double lat = bg::get_as_radian<1>(p);

    point_3d p3d;
    ::convert(p, p3d);

    // lat == 0
    if ( bg::math::equals(bg::get<2>(p3d), 0) )
        r = bg::get_radius<0>(sph);
    // |lat| == pi/2
    else if ( bg::math::equals(bg::get<0>(p3d), 0) && bg::math::equals(bg::get<1>(p3d), 0) )
        r = 0;
    // sqrt(xx+yy) - |z|/tan(lat)
    else
        r = bg::math::sqrt(bg::get<0>(p3d) * bg::get<0>(p3d) + bg::get<1>(p3d) * bg::get<1>(p3d))
            - bg::math::abs(bg::get<2>(p3d) / tan(lat));

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

void draw_point(point_3d const& p)
{
    draw_sphere(bg::get<0>(p), bg::get<1>(p), bg::get<2>(p), 0.02);
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
    point_3d p_xy = projected_to_xy(p);
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

    typedef bg::detail::andoyer_inverse<double, true, true> andoyer_t;
    andoyer_t andoyer;
    andoyer_t::result_type ai = andoyer.apply(lon1, lat1, lon2, lat2, sph);
    return std::make_pair(ai.distance, ai.azimuth);
}

std::pair<double, double> andoyer_inverse_2nd(point_geo const& p1, point_geo const& p2, spheroid const& sph)
{
    double lon1 = bg::get_as_radian<0>(p1);
    double lat1 = bg::get_as_radian<1>(p1);
    double lon2 = bg::get_as_radian<0>(p2);
    double lat2 = bg::get_as_radian<1>(p2);

    typedef bg::detail::thomas_inverse<double, true, true> thomas_t;
    thomas_t thomas;
    thomas_t::result_type ti = thomas.apply(lon1, lat1, lon2, lat2, sph);
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
    enum method_type { method_mean_point = 0, method_interpolate, method_interpolate_vertically, method_nearest } method;
    bool enable_experimental;
    bool enable_mapping_geodetic;
    bool enable_mapping_geocentric;
    bool enable_mapping_reduced;
    bool enable_great_ellipse;
    bool enable_vincenty;
    bool enable_andoyer;
    bool enable_thomas;

    std::vector<point_3d> curve_experimental;
    std::vector<std::pair<point_3d, point_3d> > lines_experimental;
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

    scene_data()
    {
        method = method_mean_point;
        enable_experimental = true;
        enable_mapping_geodetic = true;
        enable_mapping_geocentric = true;
        enable_mapping_reduced = true;
        enable_great_ellipse = true;
        enable_vincenty = true;
        enable_andoyer = true;
        enable_thomas = true;
    }

    void clear()
    {
        curve_experimental.clear();
        lines_experimental.clear();
        curve_mapped_geodetic.clear();
        curve_mapped_geocentric.clear();
        curve_mapped_reduced.clear();
        curve_great_ellipse.clear();
        curve_vincenty.clear();
        curve_andoyer.clear();
        curve_thomas.clear();
    }

    void recalculate(point_geo const& p1, point_geo const& p2)
    {
        clear();

        // cartesian points
        p1_s = pcast<point_3d>(p1);
        p2_s = pcast<point_3d>(p2);
        p1_xy = projected_to_xy_geod(p1);
        p2_xy = projected_to_xy_geod(p2);
        // cartesian vectors
        v_s = p2_s - p1_s;
        v_xy = p2_xy - p1_xy;

        double f = 0;
        int count = 50;
        double f_step = 1.0 / count;
        for ( int i = 0 ; i <= count ; ++i, f += f_step )
        {
            point_3d ps = p1_s + v_s * f;

            // experimental method
            if ( enable_experimental )
            {
                point_3d pxy;

                if ( method == method_interpolate )
                {
                    pxy = p1_xy + v_xy * f;
                }
                else if ( method == method_nearest )
                {
                    double l_sqr = bg::dot_product(v_xy, v_xy);
                    if ( !bg::math::equals(l_sqr, 0) )
                    {
                        double nt = -bg::dot_product(p1_xy - ps, p2_xy - p1_xy) / l_sqr;
                        pxy = p1_xy + v_xy * nt;
                    }
                    else
                    {
                        pxy = point_3d(0, 0, 0);
                    }
                }
                else if ( method == method_mean_point )
                {
                    pxy = p1_xy + v_xy * 0.5;
                }
                else if ( method == method_interpolate_vertically )
                {
                    point_3d v1 = p1_xy - p1_s;
                    point_3d v2 = p2_xy - p2_s;
                    // 0=ox+dx*t
                    // 0=yx+dy*t
                    // z=oz+dz*t
                    double t1 = - bg::get<0>(p1_s) / bg::get<0>(v1);
                    double t2 = - bg::get<0>(p2_s) / bg::get<0>(v2);
                    point_3d p1_v = p1_s + v1 * t1;
                    point_3d p2_v = p2_s + v2 * t2;

                    // for display purposes
                    p1_xy = p1_v;
                    p2_xy = p2_v;
                    v_xy = p2_xy - p1_xy;

                    pxy = p1_xy + v_xy * f;
                }
                else
                {
                    BOOST_ASSERT(false);
                }

                // segments between the lines
                lines_experimental.push_back(std::make_pair(pxy, ps));

                // vectors corresponding to segments
                point_3d d = ps - pxy;

                // calculate the point of intersection of a ray and spheroid's surface
                //(x*x+y*y)/(a*a) + z*z/(b*b) = 1
                // x = o.x + d.x * t
                // y = o.y + d.y * t
                // z = o.z + d.z * t        
                double ox = bg::get<0>(pxy);
                double oy = bg::get<1>(pxy);
                double oz = bg::get<2>(pxy);
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

                point_3d p_curve = pxy + d * t;

                curve_experimental.push_back(p_curve);
            }

            // great ellipse
            if ( enable_great_ellipse )
            {
                double dx = bg::get<0>(ps);
                double dy = bg::get<1>(ps);
                double dz = bg::get<2>(ps);

                double a_sqr = a*a;
                double b_sqr = b*b;
                double param_a = (dx*dx+dy*dy)/a_sqr+dz*dz/b_sqr;
                double param_c = -1;

                double delta = -4*param_a*param_c;
                double t = delta >= 0 ?
                           sqrt(delta) / (2*param_a) :
                           0.0;

                point_3d p_curve = ps * t;

                curve_great_ellipse.push_back(p_curve);
            }

            // mapped to sphere
            if ( enable_mapping_geodetic )
            {
                recalculate_mapped(p1, p2, f, curve_mapped_geodetic, mapping_geodetic);
            }
            
            if ( enable_mapping_geocentric )
            {
                recalculate_mapped(p1, p2, f, curve_mapped_geocentric, mapping_geocentric);
            }

            if ( enable_mapping_reduced )
            {
                recalculate_mapped(p1, p2, f, curve_mapped_reduced, mapping_reduced);
            }
        }

        // vincenty
        if ( enable_vincenty )
        {
            recalculate_curve< bg::detail::vincenty_inverse<double, true, true> >(p1, p2, curve_vincenty);

            //recalculate_azimuth(azi);
            recalculate_azimuth(p1, andoyer_inverse(p1, p2, sph).second);
        }

        if ( enable_andoyer )
        {
            recalculate_curve< bg::detail::andoyer_inverse<double, true, true> >(p1, p2, curve_andoyer);
        }

        if ( enable_thomas )
        {
            recalculate_curve< bg::detail::thomas_inverse<double, true, true> >(p1, p2, curve_thomas);
        }
    }

    template <typename Inverse>
    static inline void recalculate_curve(point_geo const& p1, point_geo const& p2, std::vector<point_3d> & curve)
    {
        double f_step = 1.0 / 50;
        int max_count = 100;

        // calculate the azimuth and distance from p1 to p2
        Inverse inverse;
        typename Inverse::result_type inv = inverse.apply(
                    bg::get_as_radian<0>(p1),
                    bg::get_as_radian<1>(p1),
                    bg::get_as_radian<0>(p2),
                    bg::get_as_radian<1>(p2),
                    sph);

        double azi = inv.azimuth;
        double dist = inv.distance;

        double lon = bg::get_as_radian<0>(p1);
        double lat = bg::get_as_radian<1>(p1);
        for ( int i = 0 ; i <= max_count ; ++i )
        {
            typename Inverse::result_type inv = inverse.apply(
                        lon,
                        lat,
                        bg::get_as_radian<0>(p2),
                        bg::get_as_radian<1>(p2),
                        sph);

            // calculate the position of p2 for given d and azimuth
            typedef bg::detail::vincenty_direct<double> direct_t;
            direct_t direct;
            direct_t::result_type vd = direct.apply(
                        lon,
                        lat,
                        f_step * dist,
                        inv.azimuth,
                        sph);

            lon = vd.lon2;
            lat = vd.lat2;

            {
                typename Inverse::result_type inv = inverse.apply(
                            lon,
                            lat,
                            bg::get_as_radian<0>(p2),
                            bg::get_as_radian<1>(p2),
                            sph);
                double azi_curr = inv.azimuth;
                double ba = bearing(azi, azi_curr);
                if ( ba > bg::math::pi<double>() / 2 )
                {
                    curve.push_back(pcast<point_3d>(p2));
                    return;
                }
                azi = azi_curr;
            }

            point_geo p;
            bg::set_from_radian<0>(p, lon);
            bg::set_from_radian<1>(p, lat);

            curve.push_back(pcast<point_3d>(p));
        }
    }

    point_3d loc_east;
    point_3d loc_north;
    point_3d v_azimuth;

    void recalculate_azimuth(point_geo const& p1, double azimuth)
    {
        if ( bg::math::equals(bg::get<0>(p1_s), 0.0)
          && bg::math::equals(bg::get<1>(p1_s), 0.0) )
        {
            double lon = bg::get_as_radian<0>(p1);

            bg::set<0>(loc_north, cos(lon));
            bg::set<1>(loc_north, sin(lon));
            bg::set<2>(loc_north, 0);

            bg::set<0>(loc_east, -sin(lon));
            bg::set<1>(loc_east, cos(lon));
            bg::set<2>(loc_east, 0);
        }
        else
        {
            point_3d v1_perp = p1_s - p1_xy;
            loc_east = bg::cross_product(point_3d(0, 0, 1), v1_perp);
            normalize(loc_east);
            loc_north = bg::cross_product(v1_perp, loc_east);
            normalize(loc_north);
        }
        v_azimuth = loc_north * cos(azimuth) + loc_east * sin(azimuth);
    }

    void recalculate_mapped(point_geo const& p1,
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

    void cycle_method()
    {
        method = method_type((method + 1) % 4);
    }

    void print_settings() const
    {
        std::cout << "SETTINGS\n"
                  << "method:             " << method_str() << '\n'
                  << "experimental:       " << flag_str(enable_experimental) << '\n'
                  << "mapping_geodetic:   " << flag_str(enable_mapping_geodetic) << '\n'
                  << "mapping_geocentric: " << flag_str(enable_mapping_geocentric) << '\n'
                  << "mapping_reduced:    " << flag_str(enable_mapping_reduced) << '\n'
                  << "great_ellipse:      " << flag_str(enable_great_ellipse) << '\n'
                  << "vincenty:           " << flag_str(enable_vincenty) << '\n'
                  << "andoyer:            " << flag_str(enable_andoyer) << '\n'
                  << "thomas:            " << flag_str(enable_thomas) << '\n';
        std::cout.flush();
    }

    const char * flag_str(bool f) const
    {
        return f ? "on" : "off";
    }

    const char * method_str() const
    {
        return method == method_interpolate ? "interpolate" :
               method == method_nearest ? "nearest" :
               method == method_mean_point ? "mean point" :
               "interpolate vertically";
    }

    void print_curve_lengths() const
    {
        std::cout << std::setprecision(8);
        std::cout << "LENGTHS\n"
                  << "experimental:       " << curve_length(curve_experimental) << '\n'
                  << "mapping_geodetic:   " << curve_length(curve_mapped_geodetic) << '\n'
                  << "mapping_geocentric: " << curve_length(curve_mapped_geocentric) << '\n'
                  << "mapping_reduced:    " << curve_length(curve_mapped_reduced) << '\n'
                  << "great_ellipse:      " << curve_length(curve_great_ellipse) << '\n'
                  << "vincenty:           " << curve_length(curve_vincenty) << '\n'
                  << "andoyer:            " << curve_length(curve_andoyer) << '\n'
                  << "thomas:             " << curve_length(curve_thomas) << '\n';
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
    }
}

point_geo p1(-60, -60);
point_geo p2(60, 60);

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

    glLineWidth(2);

    draw_point_adv(p1, color(1, 0.5, 0));
    draw_point_adv(p2, color(1, 1, 0));

    // draw experimental lines
    if ( data.enable_experimental )
    {
        glColor3f(1, 1, 1);
        draw_line(data.p1_s, data.p2_s);
        draw_line(data.p1_xy, data.p2_xy);

        size_t count = data.lines_experimental.size();
        double f = 0;
        double f_step = 1.0 / count;
        for ( size_t i = 0 ; i < count ; ++i, f += f_step )
        {
            glColor3f(1, 0.5+0.5*f, 0); // orange -> yellow
            draw_line(data.lines_experimental[i].first, data.lines_experimental[i].second);
        }
    }

    glLineWidth(3);
    
    if ( data.enable_experimental ) // orange -> yellow
        draw_curve(data.curve_experimental, color(1, 0.5, 0), color(1, 1, 0));

    if ( data.enable_mapping_geodetic ) // red
        draw_curve(data.curve_mapped_geodetic, color(0.5, 0, 0), color(1, 0, 0));
    if ( data.enable_mapping_geocentric ) // green
        draw_curve(data.curve_mapped_geocentric, color(0, 0.5, 0), color(0, 1, 0));
    if ( data.enable_mapping_reduced ) // blue
        draw_curve(data.curve_mapped_reduced, color(0, 0.25, 0.5), color(0, 0.5, 1));

    if ( data.enable_great_ellipse ) // pink
        draw_curve(data.curve_great_ellipse, color(0.5, 0, 0.5), color(1, 0, 1));

    if ( data.enable_vincenty ) // gray->white
        draw_curve(data.curve_vincenty, color(0.75, 0.75, 0.75), color(1, 1, 1));
    if ( data.enable_andoyer ) // magenta
        draw_curve(data.curve_andoyer, color(0.75, 0, 0.75), color(1, 0, 1));
    if ( data.enable_thomas ) // cyan
        draw_curve(data.curve_thomas, color(0, 0.75, 0.75), color(0, 1, 1));

    // azimuth
    glColor3f(1, 0, 0);
    draw_line(data.p1_s, data.p1_s + data.loc_north*0.2);
    glColor3f(0, 1, 0);
    draw_line(data.p1_s, data.p1_s + data.loc_east*0.2);
    glColor3f(1, 1, 1);
    draw_line(data.p1_s, data.p1_s + data.v_azimuth*0.3);

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
              << "h     - display help" << '\n'
              << ",     - decrease b, increase flattening" << '\n'
              << ".     - increase b, decrease flattening" << '\n'
              << "m     - experimental method switch" << '\n'
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
    std::cout << "GEOMETRY\n";
    std::cout << "p1:         (" << bg::get<0>(p1) << ", " << bg::get<1>(p1) << ")\n";
    std::cout << "p2:         (" << bg::get<0>(p2) << ", " << bg::get<1>(p2) << ")\n";
    std::cout << "flattening: " << (a - b) / a << std::endl;
}

void print_distances_and_azimuths()
{
    std::pair<double, double> ai_res = andoyer_inverse(p1, p2, sph);
    std::pair<double, double> ai2_res = andoyer_inverse_2nd(p1, p2, sph);

    std::cout << std::setprecision(8);
    std::cout << "DISTANCES\n"
              << "vincenty:     " << bg::distance(p1, p2, bg::strategy::distance::vincenty<spheroid>(sph)) << '\n'
              << "andoyer orig: " << bg::distance(p1, p2, bg::strategy::distance::andoyer<spheroid>(sph)) << '\n'
              << "haversine:    " << bg::distance(p1, p2, bg::strategy::distance::haversine<double>((2*a+b)/3)) << '\n'
              << "andoyer:      " << ai_res.first << '\n'
              << "thomas:       " << ai2_res.first << '\n';
    std::cout << "AZIMUTHS\n"
              << "vincenty:     " << bg::detail::azimuth<double>(p1, p2, sph) << '\n'
              << "andoyer:      " << ai_res.second << '\n'
              << "thomas:       " << ai2_res.second << '\n';
                  
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

    if ( key == 'h' )
    {
        print_help();
    }
    else if ( key == '.')
    {
        b += b_step;
        if ( b > 2*a-b_step )
            b = 2*a-b_step;
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
        data.cycle_method();
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
    // moving of the p2
    else if ( key == 'i' )
        move_lat(p2, 1);
    else if ( key == 'k' )
        move_lat(p2, -1);
    else if ( key == 'j' )
        move_lon(p2, -1);
    else if ( key == 'l' )
        move_lon(p2, 1);
    // other key
    else
        refresh = false;

    if ( refresh )
    {
        print_geometry();
        data.print_settings();

        data.recalculate(p1, p2);
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
        std::cout << "TEST DISTANCES\n" << std::setprecision(32)
                  << "vincenty:     " << bg::distance(p1, p2, bg::strategy::distance::vincenty<spheroid>(sph)) << '\n'
                  << "andoyer: " << bg::distance(p1, p2, bg::strategy::distance::andoyer<spheroid>(sph)) << '\n'
                  << "thomas: " << bg::distance(p1, p2, bg::strategy::distance::thomas<spheroid>(sph)) << '\n'
                  << "haversine:    " << bg::distance(p1, p2, bg::strategy::distance::haversine<double>((2*a+b)/3)) << '\n';
    }

    print_help();
    print_geometry();
    data.print_settings();

    data.recalculate(p1, p2);
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

    glutMainLoop();

    return 0;
}
