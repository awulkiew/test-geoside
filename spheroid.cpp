#include <GL/glut.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/arithmetic/cross_product.hpp>

#include <boost/geometry/formulas/vincenty_direct.hpp>
#include <boost/geometry/formulas/vincenty_inverse.hpp>
#include <boost/geometry/formulas/andoyer_inverse.hpp>
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
    enum method_type { method_mean_point = 0, method_interpolate, method_nearest, method_interpolate_vertically } method;
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

    point_3d loc_east;
    point_3d loc_north;
    point_3d v_azimuth_vincenty;
    point_3d v_azimuth_andoyer;
    point_3d v_azimuth_thomas;
    point_3d v_azimuth_great_ellipse;

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
        method = method_mean_point;
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
        lines_experimental.clear();
        curve_mapped_geodetic.clear();
        curve_mapped_geocentric.clear();
        curve_mapped_reduced.clear();
        curve_great_ellipse.clear();
        curve_vincenty.clear();
        curve_andoyer.clear();
        curve_thomas.clear();
    }

    static void calculate_north_east(point_geo const& p, point_3d & north, point_3d & east)
    {
        point_3d p_s;
        bg::formula::geo_to_cart3d(p, p_s, north, east, sph);

        return;

        p_s = pcast<point_3d>(p);

        // pole
        if (bg::math::equals(bg::get<0>(p_s), 0.0)
            && bg::math::equals(bg::get<1>(p_s), 0.0))
        {
            double lon = bg::get_as_radian<0>(p);

            bg::set<0>(north, cos(lon));
            bg::set<1>(north, sin(lon));
            bg::set<2>(north, 0);

            bg::set<0>(east, -sin(lon));
            bg::set<1>(east, cos(lon));
            bg::set<2>(east, 0);
        }
        else
        {
            point_3d p_xy = projected_to_xy_geod(p);
            point_3d v_perp = p_s - p_xy;
            east = bg::cross_product(point_3d(0, 0, 1), v_perp);
            normalize(east);
            north = bg::cross_product(v_perp, east);
            normalize(north);
        }
    }

    void recalculate_loc(point_geo const& p1, point_geo const& p2)
    {
        // cartesian points
        p1_s = pcast<point_3d>(p1);
        p2_s = pcast<point_3d>(p2);
        p1_xy = projected_to_xy_geod(p1);
        p2_xy = projected_to_xy_geod(p2);
        // cartesian vectors
        v_s = p2_s - p1_s;
        v_xy = p2_xy - p1_xy;

        calculate_north_east(p1, loc_north, loc_east);
    }

    point_3d make_v_azimuth(double azimuth)
    {
        return loc_north * cos(azimuth) + loc_east * sin(azimuth);
    }

    void recalculate(point_geo const& p1, point_geo const& p2)
    {
        clear();

        recalculate_loc(p1, p2);

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

            // mapped to sphere
            /*if ( enable_mapping_geodetic )
            {
                push_to_mapped(p1, p2, f, curve_mapped_geodetic, mapping_geodetic);
            }
            
            if ( enable_mapping_geocentric )
            {
                push_to_mapped(p1, p2, f, curve_mapped_geocentric, mapping_geocentric);
            }

            if ( enable_mapping_reduced )
            {
                push_to_mapped(p1, p2, f, curve_mapped_reduced, mapping_reduced);
            }*/
        }

        double dist_step = bg::formula::vincenty_inverse<double, true, false>::apply(
            bg::get<0>(p1) * d2r, bg::get<1>(p1) * d2r, bg::get<0>(p2) * d2r, bg::get<1>(p2) * d2r, sph
        ).distance / 50;

        if (enable_experimental)
        {
            curve_experimental.clear();
            
            if (method == method_mean_point)
                recalculate_curve< experimental_inverse<double, method_mean_point> >(p1, p2, dist_step, curve_experimental, azimuth_experimental);
            else if (method == method_interpolate)
                recalculate_curve< experimental_inverse<double, method_interpolate> >(p1, p2, dist_step, curve_experimental, azimuth_experimental);
            else if (method == method_nearest)
                recalculate_curve< experimental_inverse<double, method_nearest> >(p1, p2, dist_step, curve_experimental, azimuth_experimental);
            else if (method == method_interpolate_vertically)
                recalculate_curve< experimental_inverse<double, method_interpolate_vertically> >(p1, p2, dist_step, curve_experimental, azimuth_experimental);
        }

        if (enable_mapping_geodetic)
        {
            recalculate_curve< great_circle_inverse<double, mapping_geodetic> >(p1, p2, dist_step, curve_mapped_geodetic, azimuth_mapping_geodetic);
        }

        if (enable_mapping_geocentric)
        {
            recalculate_curve< great_circle_inverse<double, mapping_geocentric> >(p1, p2, dist_step, curve_mapped_geocentric, azimuth_mapping_geocentric);
        }

        if (enable_mapping_reduced)
        {
            recalculate_curve< great_circle_inverse<double, mapping_reduced> >(p1, p2, dist_step, curve_mapped_reduced, azimuth_mapping_reduced);
        }

        if (enable_great_ellipse)
        {
            //recalculate_great_ellipse(p1, p2, curve_great_ellipse, azimuth);
            recalculate_curve< great_ellipse_inverse<double> >(p1, p2, dist_step, curve_great_ellipse, azimuth_great_ellipse);
            v_azimuth_great_ellipse = make_v_azimuth(azimuth_great_ellipse);
        }

        if ( enable_vincenty )
        {
            recalculate_curve< bg::formula::vincenty_inverse<double, true, true> >(p1, p2, dist_step, curve_vincenty, azimuth_vincenty);
            v_azimuth_vincenty = make_v_azimuth(azimuth_vincenty);
        }

        if ( enable_andoyer )
        {
            recalculate_curve< bg::formula::andoyer_inverse<double, true, true> >(p1, p2, dist_step, curve_andoyer, azimuth_andoyer);
            v_azimuth_andoyer = make_v_azimuth(azimuth_andoyer);
        }

        if ( enable_thomas )
        {
            recalculate_curve< bg::formula::thomas_inverse<double, true, true> >(p1, p2, dist_step, curve_thomas, azimuth_thomas);
            v_azimuth_thomas = make_v_azimuth(azimuth_thomas);
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
            res.azimuth = acos(bg::dot_product(aziv, north));            
            if (bg::dot_product(aziv, east) < 0.0)
            {
                res.azimuth = -res.azimuth;
            }

            return res;
        }
    };

    template <typename T, method_type Method>
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
            
            point_3d p2_s = ::pcast<point_3d>(p2);
            point_3d p1_xy = projected_to_xy_geod(p1);
            point_3d p2_xy = projected_to_xy_geod(p2);

            //point_3d v12 = p2v - p1v;
            point_3d v12_xy = p2_xy - p1_xy;

            point_3d origin_xy(0, 0, 0);

            if (Method == method_mean_point)
            {
                origin_xy = p1_xy + v12_xy * 0.5;
            }
            else if (Method == method_interpolate)
            {
                origin_xy = p1_xy;
            }
            else if (Method == method_nearest)
            {
                double l_sqr = bg::dot_product(v12_xy, v12_xy);
                if (!bg::math::equals(l_sqr, 0))
                {
                    double nt = -bg::dot_product(p1_xy - p1_s, p2_xy - p1_xy) / l_sqr;
                    origin_xy = p1_xy + v12_xy * nt;
                }
                else
                {
                    origin_xy = point_3d(0, 0, 0);
                }
            }
            else if (Method == method_interpolate_vertically)
            {
                point_3d v1 = p1_xy - p1_s;
                //point_3d v2 = p2_xy - p2_s;
                // 0=ox+dx*t
                // 0=yx+dy*t
                // z=oz+dz*t
                double t1 = -bg::get<0>(p1_s) / bg::get<0>(v1);
                //double t2 = -bg::get<0>(p2_s) / bg::get<0>(v2);
                point_3d p1_v = p1_s + v1 * t1;
                //point_3d p2_v = p2_s + v2 * t2;

                origin_xy = p1_v;
            }

            point_3d op1 = p1_s - origin_xy;
            point_3d op2 = p2_s - origin_xy;

            point_3d n = bg::cross_product(op1, op2);
            point_3d aziv = bg::cross_product(n, op1);
            normalize(aziv);
            
            result_type res;
            res.azimuth = acos(bg::dot_product(aziv, north));
            if (bg::dot_product(aziv, east) < 0.0)
            {
                res.azimuth = -res.azimuth;
            }

            return res;
        }
    };

    static void recalculate_great_ellipse(point_geo const& p1, point_geo const& p2, std::vector<point_3d> & curve, double & azimuth)
    {
        point_3d p1_s = ::pcast<point_3d>(p1);
        point_3d p2_s = ::pcast<point_3d>(p2);
        point_3d v_s = p2_s - p1_s;

        point_3d n = bg::cross_product(p1_s, p2_s);
        point_3d aziv = bg::cross_product(n, p1_s);
        normalize(aziv);

        point_3d north, east;
        calculate_north_east(p1, north, east);

        azimuth = acos(bg::dot_product(aziv, north));
        if (bg::dot_product(aziv, east) < 0.0)
        {
            azimuth = -azimuth;
        }

        double f = 0;
        int count = 50;
        double f_step = 1.0 / count;
        for (int i = 0; i <= count; ++i, f += f_step)
        {
            point_3d ps = p1_s + v_s * f;

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

            curve.push_back(p_curve);
        }
    }

    template <typename Inverse>
    static inline void recalculate_curve(point_geo const& p1, point_geo const& p2, double dist_step, std::vector<point_3d> & curve, double & azimuth)
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

point_geo p1(-60, 0);
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

    glColor3f(1, 0.5, 0);
    draw_parallel(bg::get<1>(p1) * d2r);
    glColor3f(1, 1, 0);
    draw_parallel(bg::get<1>(p2) * d2r);

    glLineWidth(2);

    draw_point_adv(p1, color(1, 0.5, 0));
    draw_point_adv(p2, color(1, 1, 0));

    // draw experimental lines
    if ( data.enable_experimental )
    {
        /*
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
        }*/
    }

    glLineWidth(3);
    
    if ( data.enable_experimental ) // orange -> yellow
        draw_curve(data.curve_experimental, color(1, 0.5, 0), color(1, 1, 0));

    if ( data.enable_mapping_geodetic ) // red
        draw_curve(data.curve_mapped_geodetic, color(0.5, 0, 0), color(1, 0, 0));
    if ( data.enable_mapping_geocentric ) // violet
        draw_curve(data.curve_mapped_geocentric, color(0.25, 0, 0.5), color(0.5, 0, 1));
    if ( data.enable_mapping_reduced ) // blue
        draw_curve(data.curve_mapped_reduced, color(0, 0.25, 0.5), color(0, 0.5, 1));

    if ( data.enable_great_ellipse ) // green
    {
        draw_curve(data.curve_great_ellipse, color(0, 0.5, 0), color(0, 1, 0));
        //glColor3f(0, 0.5, 0);
        //draw_line(data.p1_s, data.p1_s + data.v_azimuth_great_ellipse*0.3);
        /*
        size_t count = data.curve_great_ellipse.size();
        double f = 0;
        double f_step = 1.0 / count;
        for (size_t i = 0; i < count; ++i, f += f_step)
        {
            glColor3f(0, 0.5 + i * 0.5 / count, 0);
            draw_line(point_3d(0, 0, 0), data.curve_great_ellipse[i]);
        }
        */
    }

    if (data.enable_vincenty) // gray->white
    {
        draw_curve(data.curve_vincenty, color(0.75, 0.75, 0.75), color(1, 1, 1));
        //glColor3f(0.75, 0.75, 0.75);
        //draw_line(data.p1_s, data.p1_s + data.v_azimuth_vincenty*0.3);
    }
    if (data.enable_andoyer) // magenta
    {
        draw_curve(data.curve_andoyer, color(0.75, 0, 0.75), color(1, 0, 1));
        //glColor3f(0.75, 0, 0.75);
        //draw_line(data.p1_s, data.p1_s + data.v_azimuth_andoyer*0.3);
    }
    if (data.enable_thomas) // cyan
    {
        draw_curve(data.curve_thomas, color(0, 0.75, 0.75), color(0, 1, 1));
        //glColor3f(0, 0.75, 0.75);
        //draw_line(data.p1_s, data.p1_s + data.v_azimuth_thomas*0.3);
    }

    // loc
    glColor3f(1, 0, 0);
    draw_line(data.p1_s, data.p1_s + data.loc_north*0.2);
    glColor3f(0, 1, 0);
    draw_line(data.p1_s, data.p1_s + data.loc_east*0.2);
    
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
        std::cout << "TEST DISTANCES\n" << std::setprecision(16)
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
