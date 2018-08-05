#ifndef INTF_EPS__EPS_H
#define INTF_EPS__EPS_H

#include <memory>
#include <vector>
#include <algorithm>
#include <sstream>

#ifdef EPS
#define EPS_API __declspec(dllexport)
#else
#define EPS_API __declspec(dllimport)
#endif

#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4244)
#endif

#define THROW(E, C, M) do { std::stringstream excss; throw E(((excss << __FILE__ << '(' << __LINE__ << "): error " C ": " M), excss.str())); } while (false)

namespace eps
{
static long double const pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286l;
static long double const inch = 72.l; // note that a pt in LaTeX is 72.27l. In LaTeX the unit 'bp' is equal to the eps unit
static long double const mm = inch / 25.4l;
static long double const cm = mm * 10.l;
static long double const dm = mm * 100.l;
static long double const m = mm * 1000.l;

inline float to_deg(float rad){ return rad * static_cast<float>(180 / pi); }
inline float to_rad(float deg) { return deg * static_cast<float>(pi / 180); }
}; // namespace eps

constexpr long double operator"" _deg(long double x)
{
    return x * eps::pi / 180;
}

constexpr long double operator"" _bp(long double x)
{
    return x;
}

constexpr long double operator"" _inch(long double x)
{
    return x * eps::inch;
}

constexpr long double operator"" _mm(long double x)
{
    return x * eps::mm;
}

constexpr long double operator"" _cm(long double x)
{
    return x * eps::cm;
}

constexpr long double operator"" _dm(long double x)
{
    return x * eps::dm;
}

constexpr long double operator"" _m(long double x)
{
    return x * eps::m;
}

namespace eps
{

struct rotation_t;
struct transformation_t;

inline float sqr(float a)
{
    return a*a;
}

inline float clip(float a, float min, float max)
{
    a = std::max(a, min);
    a = std::min(a, max);
    return a;
}

inline float wrap_angle(float angle)
{
    return static_cast<float>(angle - 2.f*pi * floor(angle / (2.f*pi)));
}

inline bool in_arc(float v, float a1, float a2, bool possitive)
{
    v = wrap_angle(v - a1);
    a2 = wrap_angle(a2 - a1);
    if (possitive)
    {
        return v <= a2;
    }
    else
    {
        return v >= a2;
    }
}

inline void bc_formula(float b, float c, float epsilon, int& solutions, float& solution1)
{
    if (std::abs(b) <= epsilon)
    {
        solutions = 0;
    }
    else
    {
        solutions = 1;
        solution1 = -c / b;
    }
}

inline void abc_formula(float a, float b, float c, float epsilon, int& solutions, float& solution1, float& solution2)
{
    if (std::abs(a) <= epsilon)
    {
        bc_formula(b, c, epsilon, solutions, solution1);
        return;
    }
    float d = sqr(b) - 4 * a*c;
    if (std::abs(d) <= epsilon)
    {
        d = 0.f;
        solutions = 1;
        solution1 = -0.5f*b / a;
    }
    else if (d < 0)
    {
        solutions = 0;
    }
    else
    {
        float sqrtd = std::sqrt(d);
        solution1 = 0.5f*(-b + sqrtd) / a;
        solution1 = 0.5f*(-b - sqrtd) / a;
        solutions = 2;
    }
}

inline float bezier(float a, float ai, float bi, float b, float t1)
{
    float u1 = 1.f - t1;
    float t2 = sqr(t1);
    float u2 = sqr(u1);
    float t3 = t1*t2;
    float u3 = u1*u2;
    return u3*a + 3*u2*t1*ai + 3*u1*t2*bi + t3*b;
}

struct vect_t
{
    vect_t()
    {}
    vect_t(float x, float y)
        : m_x(x)
        , m_y(y)
    {}
    vect_t& operator+=(vect_t rhs)
    {
        m_x += rhs.m_x;
        m_y += rhs.m_y;
        return *this;
    }
    vect_t& operator-=(vect_t rhs)
    {
        m_x -= rhs.m_x;
        m_y -= rhs.m_y;
        return *this;
    }
    vect_t& operator*=(float rhs)
    {
        m_x *= rhs;
        m_y *= rhs;
        return *this;
    }
    vect_t& operator*=(rotation_t const& rhs);
    vect_t& operator*=(transformation_t const& rhs);
    float m_x;
    float m_y;
};

inline vect_t operator+(vect_t lhs, vect_t rhs)
{
    return lhs += rhs;
}

inline vect_t operator-(vect_t lhs, vect_t rhs)
{
    return lhs -= rhs;
}

inline vect_t operator*(vect_t lhs, float rhs)
{
    return lhs *= rhs;
}

inline float abs(vect_t v)
{
    return sqrt(v.m_x*v.m_x + v.m_y*v.m_y);
}

struct point_t
{
    point_t()
    {}
    point_t(float x, float y)
        : m_x(x)
        , m_y(y)
    {}
    point_t& operator+=(vect_t rhs)
    {
        m_x += rhs.m_x;
        m_y += rhs.m_y;
        return *this;
    }
    point_t& operator-=(vect_t rhs)
    {
        m_x -= rhs.m_x;
        m_y -= rhs.m_y;
        return *this;
    }
    point_t& operator*=(rotation_t const& rhs);
    point_t& operator*=(transformation_t const& rhs);
    operator vect_t() const { return vect_t(m_x,m_y); }
    float m_x;
    float m_y;
};

inline point_t mid(point_t lhs, point_t rhs)
{
    return point_t((lhs.m_x + rhs.m_x) / 2, (lhs.m_y + rhs.m_y) / 2);
}

inline void min_bounding_box(point_t& min, point_t p)
{
    min.m_x = std::min(min.m_x, p.m_x);
    min.m_y = std::min(min.m_y, p.m_y);
}

inline void max_bounding_box(point_t& max, point_t p)
{
    max.m_x = std::max(max.m_x, p.m_x);
    max.m_y = std::max(max.m_y, p.m_y);
}

inline point_t operator+(point_t lhs, vect_t rhs)
{
    return lhs += rhs;
}

inline point_t operator-(point_t lhs, vect_t rhs)
{
    return lhs -= rhs;
}

inline vect_t operator-(point_t lhs, point_t rhs)
{
    return vect_t(lhs.m_x-rhs.m_x, lhs.m_y - rhs.m_y);
}

struct area_t
{
    area_t(point_t min, point_t max)
        : m_min(min)
        , m_max(max)
    {}
    point_t m_min;
    point_t m_max;
};

EPS_API std::ostream& operator<<(std::ostream& o, eps::vect_t v);
EPS_API std::ostream& operator<<(std::ostream& o, eps::point_t v);
EPS_API std::ostream& operator<<(std::ostream& o, eps::area_t a);
EPS_API std::ostream& operator<<(std::ostream& o, eps::rotation_t r);
EPS_API std::ostream& operator<<(std::ostream& o, eps::transformation_t t);
EPS_API area_t null_bounding_box();

struct rotation_t
{
    rotation_t()
    {}
    rotation_t(vect_t x, vect_t y)
        : m_x(x)
        , m_y(y)
    {}
    rotation_t& operator*=(rotation_t const& rhs);
    rotation_t operator~() const;
    vect_t m_x;
    vect_t m_y;
};

inline vect_t operator*(vect_t lhs, rotation_t const& rhs)
{
    return vect_t(
        lhs.m_x * rhs.m_x.m_x + lhs.m_y * rhs.m_y.m_x,
        lhs.m_x * rhs.m_x.m_y + lhs.m_y * rhs.m_y.m_y);
}

inline vect_t& vect_t::operator*=(rotation_t const& rhs)
{
    return *this = *this * rhs;
}

inline point_t operator*(point_t lhs, rotation_t const& rhs)
{
    return point_t(
        lhs.m_x * rhs.m_x.m_x + lhs.m_y * rhs.m_y.m_x,
        lhs.m_x * rhs.m_x.m_y + lhs.m_y * rhs.m_y.m_y);
}

inline point_t& point_t::operator*=(rotation_t const& rhs)
{
    return *this = *this * rhs;
}

inline rotation_t operator*(rotation_t const& lhs, rotation_t const& rhs)
{
    return rotation_t(
        lhs.m_x * rhs.m_x.m_x + lhs.m_y * rhs.m_y.m_x,
        lhs.m_x * rhs.m_x.m_y + lhs.m_y * rhs.m_y.m_y);
}

inline rotation_t& rotation_t::operator*=(rotation_t const& rhs)
{
    return *this = *this * rhs;
}

inline rotation_t rotation_t::operator~() const
{
    float id = 1.f / (m_x.m_x * m_y.m_y - m_y.m_x * m_x.m_y);
    return rotation_t(vect_t(m_y.m_y * id, -m_x.m_y * id), vect_t(-m_y.m_x * id, m_x.m_x * id));
}

struct transformation_t
{
    transformation_t()
    {}
    transformation_t(rotation_t r, vect_t t)
        : m_r(r)
        , m_t(t)
    {}
    transformation_t& operator*=(transformation_t const& rhs);
    transformation_t operator~() const;
    rotation_t m_r;
    vect_t m_t;
};

inline vect_t operator*(vect_t lhs, transformation_t const& rhs)
{
    return lhs * rhs.m_r;
}

inline vect_t& vect_t::operator*=(transformation_t const& rhs)
{
    return *this = *this * rhs;
}

inline point_t operator*(point_t lhs, transformation_t const& rhs)
{
    return lhs * rhs.m_r + rhs.m_t;
}

inline point_t& point_t::operator*=(transformation_t const& rhs)
{
    return *this = *this * rhs;
}

inline transformation_t transformation_t::operator~() const
{
    rotation_t ir = ~m_r;
    vect_t it = m_t * -1;
    return transformation_t(ir, it*ir);
}

inline transformation_t operator*(transformation_t const& lhs, transformation_t const& rhs)
{
    return transformation_t(
        lhs.m_r * rhs.m_r,
        lhs.m_t * rhs.m_r + rhs.m_t);
}

inline transformation_t& transformation_t::operator*=(transformation_t const& rhs)
{
    return *this = *this * rhs;
}

// see https://en.wikipedia.org/wiki/ISO_216#A_series
static long double const a0_short_side(sqrt(sqrt(0.5E12l)));
static long double const a0_long_side(sqrt(sqrt(2.E12l)));
static vect_t const a0(std::floor(a0_short_side)*mm,     std::floor(a0_long_side)*mm);
static vect_t const a1(std::floor(a0_long_side / 2)*mm,  std::floor(a0_short_side)*mm);
static vect_t const a2(std::floor(a0_short_side / 2)*mm, std::floor(a0_long_side / 2)*mm);
static vect_t const a3(std::floor(a0_long_side / 4)*mm,  std::floor(a0_short_side / 2)*mm);
static vect_t const a4(std::floor(a0_short_side / 4)*mm, std::floor(a0_long_side / 4)*mm);
static vect_t const a5(std::floor(a0_long_side / 8)*mm,  std::floor(a0_short_side / 4)*mm);
static vect_t const a6(std::floor(a0_short_side / 8)*mm, std::floor(a0_long_side / 8)*mm);

enum class cap_t
{
    default_ = 0,
    butt = 0,
    round = 1,
    square = 2
};

enum class join_t
{
    default_ = 0,
    miter = 0,
    round = 1,
    bevel = 2
};

enum class ref_t
{
    default_ = 4,
    bl = 0,
    bc = 1,
    br = 2,
    cl = 3,
    cc = 4,
    cr = 5,
    tl = 6,
    tc = 7,
    tr = 8,
    number_of_refs = 9 // must be the last!!
};

enum class text_ref_t
{
    default_ = 9,
    bl =  0,
    bc =  1,
    br =  2,
    cl =  3,
    cc =  4,
    cr =  5,
    tl =  6,
    tc =  7,
    tr =  8,
    Bl =  9,
    Bc = 10,
    Br = 11,
    number_of_text_refs = 12 // must be the last!!
};

constexpr text_ref_t to_text_ref_t(ref_t ref)
{
    return static_cast<text_ref_t>(static_cast<int>(ref));
}

struct lineending_t
{
    lineending_t() : m_used(false) {}
    virtual ~lineending_t() = 0;
    virtual void draw_procedure(std::ostream& stream) const = 0;
    virtual void draw(std::ostream& stream, eps::point_t point, float rotation, float linewidth, 
        float linercolor, float linegcolor, float linebcolor) const = 0;
    virtual float length() const { return 0; }
    bool m_used;
};
inline lineending_t::~lineending_t() = default;
EPS_API lineending_t* lineending_none();

struct linestyle_t
{
    virtual ~linestyle_t() = 0;
    virtual void draw(std::ostream& stream) const = 0;
};
inline linestyle_t::~linestyle_t() = default;
EPS_API linestyle_t* linestyle_none();

struct iproperties_t
{
    virtual ~iproperties_t() = 0;
    virtual void setlinewidth(float width) = 0;
    virtual void setlinegray(float whiteness) = 0;
    virtual void setlinergbcolor(float r, float g, float b) = 0;
    virtual void setfillgray(float whiteness) = 0;
    virtual void setfillrgbcolor(float r, float g, float b) = 0;
    virtual void setlinecap(cap_t cap) = 0;
    virtual void setlinejoin(join_t join) = 0;
    virtual void setmiterlimit(float miterlimit) = 0;
    virtual void setepsilon(float epsilon) = 0;
    virtual void setlineend(lineending_t* lineend) = 0;
    virtual void setlinebegin(lineending_t* linebegin) = 0;
    virtual void setlinestyle(linestyle_t* linestyle) = 0;
    virtual float linewidth() const = 0;
    virtual float linercolor() const = 0;
    virtual float linegcolor() const = 0;
    virtual float linebcolor() const = 0;
    virtual float fillrcolor() const = 0;
    virtual float fillgcolor() const = 0;
    virtual float fillbcolor() const = 0;
    virtual cap_t linecap() const = 0;
    virtual join_t linejoin() const = 0;
    virtual float miterlimit() const = 0;
    virtual float epsilon() const = 0;
    virtual lineending_t* lineend() const = 0;
    virtual lineending_t* linebegin() const = 0;
    virtual linestyle_t* linestyle() const = 0;
};
inline iproperties_t::~iproperties_t() = default;


class graphicsstate_t
    : public iproperties_t
{
public:
    graphicsstate_t()
        : m_linewidth(static_cast<float>(1.0_bp))
        , m_linercolor(0)
        , m_linegcolor(0)
        , m_linebcolor(0)
        , m_fillrcolor(0)
        , m_fillgcolor(0)
        , m_fillbcolor(0)
        , m_linecap(cap_t::default_)
        , m_linejoin(join_t::default_)
        , m_miterlimit(static_cast<float>(10.0_bp))
        , m_epsilon(static_cast<float>(0.25_bp))
        , m_lineend(lineending_none())
        , m_linebegin(lineending_none())
        , m_linestyle(linestyle_none())
    {}
#define SET_FUNCTION_1(TYPE, NAME) \
    void set##NAME(TYPE NAME) override \
    { \
        m_##NAME = NAME; \
    }
#define SET_FUNCTION_3(TYPE, NAME, NAME1, NAME2, NAME3) \
    void set##NAME(TYPE NAME1, TYPE NAME2, TYPE NAME3) override \
    { \
        m_##NAME1 = NAME1; \
        m_##NAME2 = NAME2; \
        m_##NAME3 = NAME3; \
    }
#define SET_FUNCTION_1_USED(TYPE, NAME) \
    void set##NAME(TYPE NAME) override \
    { \
        m_##NAME = NAME; \
        NAME->m_used = true; \
    }
    SET_FUNCTION_1(float, linewidth)
    SET_FUNCTION_3(float, linergbcolor, linercolor, linegcolor, linebcolor)
    SET_FUNCTION_3(float, fillrgbcolor, fillrcolor, fillgcolor, fillbcolor)
    SET_FUNCTION_1(cap_t, linecap)
    SET_FUNCTION_1(join_t, linejoin)
    SET_FUNCTION_1(float, miterlimit)
    SET_FUNCTION_1(float, epsilon)
    SET_FUNCTION_1_USED(lineending_t*, lineend)
    SET_FUNCTION_1_USED(lineending_t*, linebegin)
    SET_FUNCTION_1(linestyle_t*, linestyle)
#undef SET_FUNCTION_1
#undef SET_FUNCTION_3
#define GET_FUNCTION_1(TYPE, NAME) \
    TYPE NAME() const override \
    { \
        return m_##NAME; \
    }
        GET_FUNCTION_1(float, linewidth)
        GET_FUNCTION_1(float, linercolor)
        GET_FUNCTION_1(float, linegcolor)
        GET_FUNCTION_1(float, linebcolor)
        GET_FUNCTION_1(float, fillrcolor)
        GET_FUNCTION_1(float, fillgcolor)
        GET_FUNCTION_1(float, fillbcolor)
        GET_FUNCTION_1(cap_t, linecap)
        GET_FUNCTION_1(join_t, linejoin)
        GET_FUNCTION_1(float, miterlimit)
        GET_FUNCTION_1(float, epsilon)
        GET_FUNCTION_1(lineending_t*, lineend)
        GET_FUNCTION_1(lineending_t*, linebegin)
        GET_FUNCTION_1(linestyle_t*, linestyle)
#undef GET_FUNCTION_1
    bool operator!=(graphicsstate_t const& rhs) const
    {
        if (m_linewidth != rhs.m_linewidth) return true;
        if (m_linercolor != rhs.m_linercolor) return true;
        if (m_linegcolor != rhs.m_linegcolor) return true;
        if (m_linebcolor != rhs.m_linebcolor) return true;
        if (m_fillrcolor != rhs.m_fillrcolor) return true;
        if (m_fillgcolor != rhs.m_fillgcolor) return true;
        if (m_fillbcolor != rhs.m_fillbcolor) return true;
        if (m_linecap != rhs.m_linecap) return true;
        if (m_linejoin != rhs.m_linejoin) return true;
        if (m_miterlimit != rhs.m_miterlimit) return true;
        return false;
    }
private:
#define SET_FUNCTION_1_THROW(TYPE, NAME) \
    void set##NAME(TYPE NAME) override \
    { \
        THROW(std::logic_error, "E0002", << "Don't use the function 'set" #NAME "'"); \
    }
    SET_FUNCTION_1_THROW(float, linegray)
    SET_FUNCTION_1_THROW(float, fillgray)
#undef SET_FUNCTION_1_THROW
    float m_linewidth;
    float m_linercolor;
    float m_linegcolor;
    float m_linebcolor;
    float m_fillrcolor;
    float m_fillgcolor;
    float m_fillbcolor;
    cap_t m_linecap;
    join_t m_linejoin;
    float m_miterlimit;
    float m_epsilon; // accuracy margin for simplifications
    lineending_t* m_lineend;
    lineending_t* m_linebegin;
    linestyle_t* m_linestyle;
};


class properties_override_t
    : public graphicsstate_t
{
public:
    properties_override_t()
        : m_linewidth_overridden(false)
        , m_linergbcolor_overridden(false)
        , m_fillrgbcolor_overridden(false)
        , m_linecap_overridden(false)
        , m_linejoin_overridden(false)
        , m_miterlimit_overridden(false)
        , m_epsilon_overridden(false)
        , m_lineend_overridden(false)
        , m_linebegin_overridden(false)
        , m_linestyle_overridden(false)
        , m_ref_count(0)
    {}
    properties_override_t(graphicsstate_t& other)
        : graphicsstate_t(other)
        , m_linewidth_overridden(false)
        , m_linergbcolor_overridden(false)
        , m_fillrgbcolor_overridden(false)
        , m_linecap_overridden(false)
        , m_linejoin_overridden(false)
        , m_miterlimit_overridden(false)
        , m_epsilon_overridden(false)
        , m_lineend_overridden(false)
        , m_linebegin_overridden(false)
        , m_linestyle_overridden(false)
        , m_ref_count(0)
    {}

#define COMPARE_1_PROPERTY(NAME) \
    if (m_##NAME##_overridden < rhs.m_##NAME##_overridden) return true; \
    if (rhs.m_##NAME##_overridden < m_##NAME##_overridden) return false; \
    if (m_##NAME##_overridden) \
    { \
        if (graphicsstate_t::NAME() < rhs.graphicsstate_t::NAME()) return true; \
        if (rhs.graphicsstate_t::NAME() < graphicsstate_t::NAME()) return false; \
    }
#define COMPARE_3_PROPERTY(NAME, NAME1, NAME2, NAME3) \
    if (m_##NAME##_overridden < rhs.m_##NAME##_overridden) return true; \
    if (rhs.m_##NAME##_overridden < m_##NAME##_overridden) return false; \
    if (m_##NAME##_overridden) \
    { \
        if (graphicsstate_t::NAME1() < rhs.graphicsstate_t::NAME1()) return true; \
        if (rhs.graphicsstate_t::NAME1() < graphicsstate_t::NAME1()) return false; \
        if (graphicsstate_t::NAME2() < rhs.graphicsstate_t::NAME2()) return true; \
        if (rhs.graphicsstate_t::NAME2() < graphicsstate_t::NAME2()) return false; \
        if (graphicsstate_t::NAME3() < rhs.graphicsstate_t::NAME3()) return true; \
        if (rhs.graphicsstate_t::NAME3() < graphicsstate_t::NAME3()) return false; \
    }
    bool operator<(properties_override_t const& rhs) const
    {
        COMPARE_1_PROPERTY(linewidth);
        COMPARE_3_PROPERTY(linergbcolor, linercolor, linegcolor, linebcolor);
        COMPARE_3_PROPERTY(fillrgbcolor, fillrcolor, fillgcolor, fillbcolor);
        COMPARE_1_PROPERTY(linecap);
        COMPARE_1_PROPERTY(linejoin);
        COMPARE_1_PROPERTY(miterlimit);
        COMPARE_1_PROPERTY(epsilon);
        COMPARE_1_PROPERTY(lineend);
        COMPARE_1_PROPERTY(linebegin);
        COMPARE_1_PROPERTY(linestyle);
        return false;
    }
#undef COMPARE_1_PROPERTY
#undef COMPARE_3_PROPERTY
#define SET_FUNCTION_1A(TYPE, NAME) \
    void set##NAME(TYPE NAME) override \
    { \
        graphicsstate_t::set##NAME(NAME); \
        m_##NAME##_overridden = true; \
    }
#define SET_FUNCTION_1B(TYPE, NAME, NAMEB) \
    void set##NAME(TYPE NAME) override \
    { \
        graphicsstate_t::set##NAMEB(NAME, NAME, NAME); \
        m_##NAMEB##_overridden = true; \
    }
#define SET_FUNCTION_3(TYPE, NAME, NAME1, NAME2, NAME3) \
    void set##NAME(TYPE NAME1, TYPE NAME2, TYPE NAME3) override \
    { \
        graphicsstate_t::set##NAME(NAME1, NAME2, NAME3); \
        m_##NAME##_overridden = true; \
    }
    SET_FUNCTION_1A(float, linewidth)
    SET_FUNCTION_1B(float, linegray, linergbcolor)
    SET_FUNCTION_3(float, linergbcolor, linercolor, linegcolor, linebcolor)
    SET_FUNCTION_1B(float, fillgray, fillrgbcolor)
    SET_FUNCTION_3(float, fillrgbcolor, fillrcolor, fillgcolor, fillbcolor)
    SET_FUNCTION_1A(cap_t, linecap)
    SET_FUNCTION_1A(join_t, linejoin)
    SET_FUNCTION_1A(float, miterlimit)
    SET_FUNCTION_1A(float, epsilon)
    SET_FUNCTION_1A(lineending_t*, lineend)
    SET_FUNCTION_1A(lineending_t*, linebegin)
    SET_FUNCTION_1A(linestyle_t*, linestyle)
#undef SET_FUNCTION_1A
#undef SET_FUNCTION_1B
#undef SET_FUNCTION_3
#define GET_FUNCTION_1A(TYPE, NAME1) \
    TYPE NAME1(TYPE DEFAULT) const \
    { \
        return m_##NAME1##_overridden ? graphicsstate_t::NAME1() : DEFAULT; \
    }
#define GET_FUNCTION_1B(TYPE, NAME1, NAME2) \
    TYPE NAME2(TYPE DEFAULT) const \
    { \
        return m_##NAME1##_overridden ? graphicsstate_t::NAME2() : DEFAULT; \
    }
        GET_FUNCTION_1A(float, linewidth)
        GET_FUNCTION_1B(float, linergbcolor, linercolor)
        GET_FUNCTION_1B(float, linergbcolor, linegcolor)
        GET_FUNCTION_1B(float, linergbcolor, linebcolor)
        GET_FUNCTION_1B(float, fillrgbcolor, fillrcolor)
        GET_FUNCTION_1B(float, fillrgbcolor, fillgcolor)
        GET_FUNCTION_1B(float, fillrgbcolor, fillbcolor)
        GET_FUNCTION_1A(cap_t, linecap)
        GET_FUNCTION_1A(join_t, linejoin)
        GET_FUNCTION_1A(float, miterlimit)
        GET_FUNCTION_1A(float, epsilon)
        GET_FUNCTION_1A(lineending_t*, lineend)
        GET_FUNCTION_1A(lineending_t*, linebegin)
        GET_FUNCTION_1A(linestyle_t*, linestyle)
#undef GET_FUNCTION_1A
#undef GET_FUNCTION_1B
private:
    friend class shape_t;
    bool m_linewidth_overridden;
    bool m_linergbcolor_overridden;
    bool m_fillrgbcolor_overridden;
    bool m_linecap_overridden;
    bool m_linejoin_overridden;
    bool m_miterlimit_overridden;
    bool m_epsilon_overridden;
    bool m_lineend_overridden;
    bool m_linebegin_overridden;
    bool m_linestyle_overridden;
    mutable int m_ref_count;
};

EPS_API float get_epsilon(graphicsstate_t& graphicsstate);

// see ftp://ftp.dante.de/tex-archive/info/epslatex.pdf
// see https://www-cdf.fnal.gov/offline/PostScript/PLRM2.pdf
EPS_API void new_path(std::ostream& stream);
EPS_API void closepath(std::ostream& stream);
EPS_API void stroke(std::ostream& stream, graphicsstate_t& graphicsstate, iproperties_t const& properties);
EPS_API void fill(std::ostream& stream, graphicsstate_t& graphicsstate, iproperties_t const& properties, bool and_stroke);
EPS_API void gsave(std::ostream& stream);
EPS_API void grestore(std::ostream& stream);
EPS_API void initgraphics(std::ostream& stream, graphicsstate_t& graphicsstate);
EPS_API void moveto(std::ostream& stream, eps::point_t p);
EPS_API void rmoveto(std::ostream& stream, eps::vect_t v);
EPS_API void lineto(std::ostream& stream, eps::point_t p);
EPS_API void rlineto(std::ostream& stream, eps::vect_t v);
EPS_API void curveto(std::ostream& stream, eps::point_t tangent1, eps::point_t tangent2, eps::point_t end);
EPS_API void rcurveto(std::ostream& stream, eps::vect_t tangent1, eps::vect_t tangent2, eps::vect_t end);
EPS_API void arc(std::ostream& stream, eps::point_t center, float radius, float begin_angle/*deg*/, float end_angle/*deg*/);
EPS_API void arcn(std::ostream& stream, eps::point_t center, float radius, float begin_angle/*deg*/, float end_angle/*deg*/);
EPS_API void arct(std::ostream& stream, eps::point_t tangent, eps::point_t end, float radius);
EPS_API void show(std::ostream& stream, std::string const& text);
EPS_API void showlatex(std::ostream& stream, std::string const& text, text_ref_t rext_ref = text_ref_t::Bl, float scale = 1, float rotate/*deg*/ = 0);
inline  void showlatex(std::ostream& stream, std::string const& text, float scale = 1, float rotate/*deg*/ = 0) { showlatex(stream, text, text_ref_t::Bl, scale, rotate); }
EPS_API void clip(std::ostream& stream);
EPS_API void pushmatrix(std::ostream& stream);
EPS_API void concatmatrix(std::ostream& stream, eps::transformation_t const& transformation);
EPS_API void popmatrix(std::ostream& stream);
EPS_API void scale(std::ostream& stream, float x, float y);
EPS_API void rotate(std::ostream& stream, float angle);
EPS_API void concat(std::ostream& stream, transformation_t t);
EPS_API void begin_procedure(std::ostream& stream, std::string const& name, std::vector<char const*> l);
EPS_API void end_procedure(std::ostream& stream);
EPS_API void call_procedure(std::ostream& stream, std::string const& name);
EPS_API void draw_ellipse(std::ostream& stream, point_t a, point_t c, vect_t rx, vect_t ry, float epsilon);
EPS_API void draw_arc(std::ostream& stream, point_t a, point_t c, vect_t rx, vect_t ry, point_t b, float epsilon);
EPS_API area_t bezier_bounding_box(point_t a, point_t ai, point_t bi, point_t b, float epsilon);
EPS_API area_t ellipse_bounding_box(point_t a, point_t c, vect_t rx, vect_t ry, float epsilon);
EPS_API area_t arc_bounding_box(point_t a, point_t c, vect_t rx, vect_t ry, point_t b, float epsilon);

EPS_API std::string single_line_it(std::string const& in);
EPS_API std::string single_line_listing(std::string const& in, char delim_char = '$');
EPS_API std::string multi_line_listing(float width, std::string const& in, char delim_char = '$');

class shape_t
    : public iproperties_t
{
public:
    EPS_API explicit shape_t(iproperties_t const& parent_properties);
    EPS_API shape_t(shape_t const& other);
    shape_t& operator=(shape_t const& other) = delete;
    EPS_API virtual ~shape_t();
    void rotate(float angle, bool exluding_text = false)
    {
        apply(transformation_t(rotation_t(
            vect_t(+cos(to_rad(angle)), +sin(to_rad(angle))),
            vect_t(-sin(to_rad(angle)), +cos(to_rad(angle))))
            , vect_t(0, 0)), exluding_text);
    }
    void move(vect_t v, bool exluding_text = false)
    {
        apply(transformation_t(rotation_t(
            vect_t(1, 0),
            vect_t(0, 1))
            , v), exluding_text);
    }
    void scale(float s, bool exluding_text = false)
    {
        apply(transformation_t(rotation_t(
            vect_t(s, 0),
            vect_t(0, s))
            , vect_t(0, 0)), exluding_text);
    }
    void scale(float s_x, float s_y, bool exluding_text = false)
    {
        apply(transformation_t(rotation_t(
            vect_t(s_x, 0),
            vect_t(0, s_y))
            , vect_t(0, 0)), exluding_text);
    }
    void vmirror(bool exluding_text = false)
    {
        apply(transformation_t(rotation_t(
            vect_t(-1, 0),
            vect_t(0, 1))
            , vect_t(0, 0)), exluding_text);
    }
    void hmirror(bool exluding_text = false)
    {
        apply(transformation_t(rotation_t(
            vect_t(1, 0),
            vect_t(0, -1))
            , vect_t(0, 0)), exluding_text);
    }

    EPS_API void setlinewidth(float width) override;
    EPS_API void setlinegray(float whiteness) override;
    EPS_API void setlinergbcolor(float r, float g, float b) override;
    EPS_API void setfillgray(float whiteness) override;
    EPS_API void setfillrgbcolor(float r, float g, float b) override;
    EPS_API void setlinecap(cap_t cap) override;
    EPS_API void setlinejoin(join_t join) override;
    EPS_API void setmiterlimit(float miterlimit) override;
    EPS_API void setepsilon(float epsilon) override;
    EPS_API void setlineend(lineending_t* lineend) override;
    EPS_API void setlinebegin(lineending_t* linebegin) override;
    EPS_API void setlinestyle(linestyle_t* linestyle) override;
    EPS_API float linewidth() const override;
    EPS_API float linercolor() const override;
    EPS_API float linegcolor() const override;
    EPS_API float linebcolor() const override;
    EPS_API float fillrcolor() const override;
    EPS_API float fillgcolor() const override;
    EPS_API float fillbcolor() const override;
    //EPS_API float dash() const override;
    //EPS_API float dashspace() const override;
    //EPS_API float dashoffset() const override;
    EPS_API cap_t linecap() const override;
    EPS_API join_t linejoin() const override;
    EPS_API float miterlimit() const override;
    EPS_API float epsilon() const override;
    EPS_API lineending_t* lineend() const override;
    EPS_API lineending_t* linebegin() const override;
    EPS_API linestyle_t* linestyle() const override;
    virtual eps::area_t bounding_box(float epsilon) = 0;
protected:
    virtual bool pure_text() { return false; }
private:
    iproperties_t const& m_parent_properties;
    properties_override_t const* m_pproperties_override;
    virtual void draw(std::ostream& stream, eps::graphicsstate_t& graphicsstate) const = 0;
    virtual void apply(transformation_t const&, bool exluding_text) = 0;
    friend class group_t;

    void inc_ref();
    void dec_ref();
    void get(properties_override_t& properties_override);
    void add(properties_override_t& properties_override);
};

// see https://www-cdf.fnal.gov/offline/PostScript/PLRM2.pdf

template<int S>
class fixed_shape_t // fixed sized point_t array
    : public shape_t
{
protected:
    fixed_shape_t(iproperties_t const& parent_properties)
        : shape_t(parent_properties)
    {}
    point_t m_[S];
    static int const size = S;
    void apply(transformation_t const& t, bool exluding_text) final
    {
        if (!(exluding_text && pure_text()))
        {
            for (eps::point_t& p : m_)
            {
                p *= t;
            }
        }
    }
    eps::area_t bounding_box(float epsilon) override
    {
        eps::area_t area = null_bounding_box();
        for (eps::point_t p : m_)
        {
            min_bounding_box(area.m_min, p);
            max_bounding_box(area.m_max, p);
        }
        return area;
    }
public:
    eps::point_t& operator[](int idx)
    {
        if (idx >= size)
        {
            THROW(std::logic_error, "E0003", << "idx '" << idx << "' >= size '" << size << "'");
        }
        return m_[idx];
    }
};

class dynamic_shape_t // dynamically sized point_t array
    : public shape_t
{
protected:
    explicit dynamic_shape_t(iproperties_t const& parent_properties)
        : shape_t(parent_properties)
    {}
    std::vector<point_t> m_;
    void apply(transformation_t const& t, bool exluding_text) final
    {
        if (!(exluding_text && pure_text()))
        {
            for (point_t& p : m_)
            {
                p *= t;
            }
        }
    }
    area_t bounding_box(float epsilon) override
    {
        area_t area = null_bounding_box();
        for (point_t p : m_)
        {
            min_bounding_box(area.m_min, p);
            max_bounding_box(area.m_max, p);
        }
        return area;
    }
public:
    eps::point_t& operator[](int idx)
    {
        if (idx >= static_cast<int>(m_.size()))
        {
            THROW(std::logic_error, "E0003", << "idx '" << idx << "' >= size '" << static_cast<int>(m_.size()) << "'");
        }
        return m_[idx];
    }
};

class group_t
    : public shape_t
{
public:
    explicit group_t(iproperties_t const& parent_properties)
        : shape_t(parent_properties)
    {}
    EPS_API virtual void add(std::unique_ptr<shape_t>&& o);
    EPS_API virtual area_t bounding_box(float epsilon) override;
protected:
    EPS_API virtual void draw(std::ostream& stream, eps::graphicsstate_t& graphicsstate) const override;
private:
    EPS_API virtual void apply(transformation_t const&, bool excluding_text) override;
//    EPS_API virtual point_t& operator[](int idx);
    std::vector<std::unique_ptr<shape_t>> m_shapes;
friend class canvas_t;
};

class canvas_t
    : public group_t
{
public:
    canvas_t(iproperties_t const& root_properties)
        : group_t(root_properties)
    {}
    virtual void draw() = 0;
};

EPS_API std::unique_ptr<canvas_t> create_canvas(std::string const& filename);

template<typename T, typename U, typename... ARGS>
T* add(std::unique_ptr<U>& pgroup, ARGS&&... args)
{
    std::unique_ptr<T> p(std::make_unique<T>(*pgroup.get(), args...));
    T* rp = p.get();
    pgroup->add(std::move(p));
    return rp;
}

template<typename T, typename P, typename... ARGS>
T* add(P* pgroup, ARGS&&... args)
{
    std::unique_ptr<T> p(std::make_unique<T>(*pgroup, args...));
    T* rp = p.get();
    pgroup->add(std::move(p));
    return rp;
}

EPS_API void handle_exception();

template<typename T>
void main(char const* file_name)
{
    try
    {
        std::unique_ptr<canvas_t> canvas =
            create_canvas(file_name);
        add<T>(canvas);
        canvas->draw();
    }
    catch (...)
    {
        handle_exception();
    }
}

}; // namespace eps

#ifdef _MSC_VER
#pragma warning( pop )
#endif // _MSC_VER

#endif // INTF_EPS__EPS_H
