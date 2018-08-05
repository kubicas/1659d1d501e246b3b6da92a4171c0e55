#ifndef INTF_EPS__EPS_BASIC_SHAPES_H
#define INTF_EPS__EPS_BASIC_SHAPES_H

#include "eps.h"
#include <ostream>
#include <string>

namespace eps
{

/******************************************************************************
basic shapes
******************************************************************************/

class line_t
    : public fixed_shape_t<2>
{
public:
    enum points_t
    {
        a,
        b,
        number_of_points // must be the last!!
    };
    static_assert(number_of_points == size);
    line_t(iproperties_t const& parent_properties, point_t pa, point_t pb)
        : fixed_shape_t<2>(parent_properties)
    {
        m_[a] = pa;
        m_[b] = pb;
    }
private:
    void draw(std::ostream& stream, graphicsstate_t& graphicsstate) const override
    {
        new_path(stream);
        vect_t v = m_[b] - m_[a];
        float end_angle = eps::to_deg(std::atan2f(v.m_y, v.m_x));
        float begin_angle = end_angle + 180;
        v *= 1.f / eps::abs(v);
        point_t begin = m_[a] + v * linebegin()->length();
        point_t end = m_[b] - v * lineend()->length();
        linebegin()->draw(stream, m_[a], begin_angle, linewidth(),
            linercolor(), linegcolor(), linebcolor());
        moveto(stream, begin);
        lineto(stream, end);
        stroke(stream, graphicsstate, *this);
        lineend()->draw(stream, m_[b], end_angle, linewidth(),
            linercolor(), linegcolor(), linebcolor());
    }
};

class bezier_t
    : public fixed_shape_t<4>
{
public:
    enum points_t
    {
        a,
        ai,
        bi,
        b,
        number_of_points // must be the last!!
    };
    static_assert(number_of_points == size);
    bezier_t(iproperties_t const& parent_properties, point_t pa, point_t pai, point_t pbi, point_t pb)
        : fixed_shape_t<4>(parent_properties)
    {
        m_[a] = pa;
        m_[ai] = pai;
        m_[bi] = pbi;
        m_[b] = pb;
    }
private:
    void draw(std::ostream& stream, graphicsstate_t& graphicsstate) const override
    {
        new_path(stream);
        moveto(stream, m_[a]);
        curveto(stream, m_[ai], m_[bi], m_[b]);
        stroke(stream, graphicsstate, *this);
    }
    area_t bounding_box(float epsilon) override
    {
        return bezier_bounding_box(m_[a], m_[ai], m_[bi], m_[b], epsilon);
    }
};

class arc_t
    : public fixed_shape_t<5>
{
public:
    enum points_t
    {
        a,
        c,
        cx,
        cy,
        b,
        number_of_points // must be the last!!
    };
    static_assert(number_of_points == size);
    arc_t(iproperties_t const& parent_properties, point_t pa, point_t pc, point_t pcx, point_t pcy, point_t pb)
        : fixed_shape_t<5>(parent_properties)
    {
        m_[a] = pa;
        m_[c] = pc;
        m_[cx] = pcx;
        m_[cy] = pcy;
        m_[b] = pb;
    }
private:
    void draw(std::ostream& stream, graphicsstate_t& graphicsstate) const override
    {
        new_path(stream);
        draw_arc(stream, m_[a], m_[c], m_[cx] - m_[c], m_[cy] - m_[c], m_[b], get_epsilon(graphicsstate));
        stroke(stream, graphicsstate, *this);
    }
    area_t bounding_box(float epsilon) override
    {
        return arc_bounding_box(m_[a], m_[c], m_[cx] - m_[c], m_[cy] - m_[c], m_[b], epsilon);
    }
};

class rectangle_t
    : public fixed_shape_t<4>
{
public:
    enum points_t
    {
        tr,
        tl,
        bl,
        br,
        number_of_points // must be the last!!
    };
    static_assert(number_of_points == size);
    rectangle_t(iproperties_t const& parent_properties, float side_x, float side_y)
        : fixed_shape_t<4>(parent_properties)
        , m_fill(false)
    {
        m_[tr] = point_t(+side_x / 2, +side_y / 2);
        m_[tl] = point_t(-side_x / 2, +side_y / 2);
        m_[bl] = point_t(-side_x / 2, -side_y / 2);
        m_[br] = point_t(+side_x / 2, -side_y / 2);
    }
    rectangle_t(iproperties_t const& parent_properties, float side_x)
        : rectangle_t(parent_properties, side_x, side_x)
    {}
    void fill() { m_fill = true; };
private:
    void draw(std::ostream& stream, graphicsstate_t& graphicsstate) const override
    {
        new_path(stream);
        moveto(stream, m_[tr]);
        lineto(stream, m_[tl]);
        lineto(stream, m_[bl]);
        lineto(stream, m_[br]);
        closepath(stream);
        if (m_fill)
        {
            eps::fill(stream, graphicsstate, *this, true);
        }
        else
        {
            stroke(stream, graphicsstate, *this);
        }
    }
    bool m_fill;
};

class ellipse_t
    : public fixed_shape_t<4>
{
public:
    enum points_t
    {
        tr,
        tl,
        bl,
        br,
        number_of_points // must be the last!!
    };
    static_assert(number_of_points == size);
    ellipse_t(iproperties_t const& parent_properties, float radius_x, float radius_y)
        : fixed_shape_t<4>(parent_properties)
        , m_fill(false)
    {
        m_[tr] = point_t( radius_x,  radius_y);
        m_[tl] = point_t(-radius_x,  radius_y);
        m_[bl] = point_t(-radius_x, -radius_y);
        m_[br] = point_t( radius_x, -radius_y);
    }
    ellipse_t(iproperties_t const& parent_properties, float radius_x)
        : ellipse_t(parent_properties, radius_x, radius_x)
    {}
    void fill() { m_fill = true; };
private:
    void draw(std::ostream& stream, graphicsstate_t& graphicsstate) const override
    {
        new_path(stream);
        point_t c = mid(m_[tr], m_[bl]);
        point_t a = mid(m_[tr], m_[br]);
        vect_t rx = (m_[br] - m_[bl]) * 0.5;
        vect_t ry = (m_[tr] - m_[br]) * 0.5;
        draw_ellipse(stream, a, c, rx, ry, get_epsilon(graphicsstate));
        // ??? clospath ????
        if (m_fill)
        {
            eps::fill(stream, graphicsstate, *this, true);
        }
        else
        {
            stroke(stream, graphicsstate, *this);
        }
    }
    area_t bounding_box(float epsilon) override
    {
        point_t c = mid(m_[tr], m_[bl]);
        point_t a = mid(m_[tr], m_[br]);
        vect_t rx = (m_[br] - m_[bl]) * 0.5;
        vect_t ry = (m_[tr] - m_[br]) * 0.5;
        return ellipse_bounding_box(a, c, rx, ry, epsilon);;
    }
    bool m_fill;
};

struct text_t
    : public fixed_shape_t<2>
{
public:
    enum points_t
    {
        a,
        b,
        number_of_points // must be the last!!
    };
    static_assert(number_of_points == size);
    text_t(iproperties_t const& parent_properties, std::string const& text, text_ref_t text_ref=text_ref_t::default_, float scale=1, float rotate=0/*deg*/)
        : fixed_shape_t<2>(parent_properties)
        , m_s(text)
        , m_tr(text_ref)
    {
        m_[a] = point_t(0, 0);
        rotate = to_rad(rotate);
        m_[b] = point_t(scale*std::cos(rotate), scale*std::sin(rotate));
    }
    text_t(iproperties_t const& parent_properties, std::string const& text, float scale, float rotate = 0/*deg*/)
        : text_t(parent_properties, text, text_ref_t::default_, scale, rotate)
    {}
private:
    bool pure_text() override { return true; }
    void draw(std::ostream& stream, graphicsstate_t& graphicsstate) const override
    {
        moveto(stream, m_[a]);
        vect_t v = m_[b] - m_[a];
        float scale = abs(v);
        float rotate = std::atan2(v.m_y, v.m_x);
        showlatex(stream, m_s, m_tr, scale, rotate);
    }
    area_t bounding_box(float epsilon) override
    {
        area_t area = null_bounding_box();
        // ignore m_[b], and assume text size (0,0)
        min_bounding_box(area.m_min, m_[a]);
        max_bounding_box(area.m_max, m_[a]);
        return area;

    }
    std::string m_s;
    text_ref_t m_tr;
};

class section_t
{
public:
    virtual ~section_t() {};
};

class path_t
    : public dynamic_shape_t
{
public:
    path_t(iproperties_t const& parent_properties);
    path_t(path_t const& rhs);
    EPS_API void moveto(point_t p);
    EPS_API void lineto(point_t p);
    EPS_API void curveto(point_t tangent1, point_t tangent2, point_t end);
    EPS_API void arcto(point_t center, point_t x_ax, point_t y_ax, point_t end);
    EPS_API void closepath();
    void fill() { m_fill = true; };
private:
    area_t bounding_box(float epsilon) override;
    void draw(std::ostream& stream, graphicsstate_t& graphicsstate) const override;
    std::vector<std::unique_ptr<section_t>> m_sections;
    bool m_fill;
};

/******************************************************************************
basic line endings
******************************************************************************/

struct lineending_triangle_arrow_t
    : lineending_t
{
    lineending_triangle_arrow_t(std::string const& name, float length, float width, bool fill, bool close, bool through)
        : m_name(name)
        , m_length(length)
        , m_width(width)
        , m_fill(fill)
        , m_close(close)
        , m_through(through)
    {}
    void draw_procedure(std::ostream& stream) const override
    {
        begin_procedure(stream, m_name, { "x", "y", "angle", "lw", "lcr", "lcg", "lcb" });
        stream << "x y translate\n";
        stream << "angle rotate\n";
        stream << "lw setlinewidth\n";
        stream << "lcr lcg lcb setcolor\n";
        new_path(stream);
        moveto(stream, point_t(-m_length, +m_width / 2));
        rlineto(stream, vect_t(+m_length, -m_width / 2));
        rlineto(stream, vect_t(-m_length, -m_width / 2));
        if (m_close)
        {
            closepath(stream);
        }
        if (m_fill)
        {
            stream << "gsave\n";
            stream << "fill\n";
            stream << "grestore\n";
        }
        stream << "stroke\n";
        end_procedure(stream);
    }
    void draw(std::ostream& stream, point_t point, float rotation, float linewidth,
        float linercolor, float linegcolor, float linebcolor) const override
    {
        stream << point << ' ' << rotation << ' ' << linewidth << ' ' 
            << linercolor << ' ' << linegcolor << ' ' << linebcolor << ' ' << m_name << std::endl;
    }
    float length() const override
    {
        return m_through ? 0.f : m_length;
    }
    std::string m_name;
    float const m_length;
    float const m_width;
    bool const m_fill;
    bool const m_close;
    bool const m_through;
};

struct lineending_circle_center_t
    : lineending_t
{
    lineending_circle_center_t(std::string const& name, float radius, bool fill)
        : m_name(name)
        , m_radius(radius)
        , m_fill(fill)
    {}
    void draw_procedure(std::ostream& stream) const override
    {
        begin_procedure(stream, m_name, { "x", "y", "lw", "lcr", "lcg", "lcb" });
        stream << "lw setlinewidth\n";
        stream << "lcr lcg lcb setcolor\n";
        new_path(stream);
        stream << "x y " << m_radius << " 0 360 arc\n";
        if (m_fill)
        {
            stream << "gsave\n";
            stream << "fill\n";
            stream << "grestore\n";
        }
        stream << "stroke\n";
        end_procedure(stream);
    }
    void draw(std::ostream& stream, point_t point, float rotation, float linewidth,
        float linercolor, float linegcolor, float linebcolor) const override
    {
        stream << point << ' ' << linewidth << ' '
            << linercolor << ' ' << linegcolor << ' ' << linebcolor << ' ' << m_name << std::endl;
    }
    float length() const override
    {
        return m_radius;
    }
    std::string m_name;
    float const m_radius;
    bool const m_fill;
};

/******************************************************************************
basic line styles
******************************************************************************/

struct linestyle_symetric_dash_t
    : linestyle_t
{
    explicit linestyle_symetric_dash_t(float length)
        : m_length(length)
    {}
    void draw(std::ostream& stream) const override
    {
        stream << "[" << m_length << "] " << m_length*1.5 << " setdash\n";
    }
    float const m_length;
};

struct linestyle_asymetric_dash_t
    : linestyle_t
{
    explicit linestyle_asymetric_dash_t(float on, float off)
        : m_on(on)
        , m_off(off)
    {}
    void draw(std::ostream& stream) const override
    {
        stream << "[" << m_on << ' ' << m_off << "] " << m_on + m_off*0.5 << " setdash\n";
    }
    float const m_on;
    float const m_off;
};

}; // namespace eps

#endif // INTF_EPS__EPS_BASIC_SHAPES_H