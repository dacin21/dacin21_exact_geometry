// Released under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007, see the LICENSE file.
// Copyright (C) 2018-2019 Daniel Rutschmann aka. dacin21

#ifndef GEOM_2D_HPP
#define GEOM_2D_HPP

#include "geom_utility.hpp"
#include "adaptive_int.hpp"

namespace dacin::geom{

template<size_t n>
class Point{
public:
    using coord_t = Adaptive_Int<n>;

    coord_t x, y;

    Point() : x(), y() {}
    template<size_t m, typename = enable_if_t<m <= n> >
    Point(Point<m> const&o) : x(o.x), y(o.y) {}
    template<typename T, typename = enable_if_t<is_constructible_v<coord_t, T> > >
    Point(T const&x_, T const&y_) : x(x_), y(y_) {}

    explicit operator std::pair<double, double>() const {
        return std::make_pair(static_cast<double>(x), static_cast<double>(y));
    }
    explicit operator std::pair<long double, long double>() const {
        return std::make_pair(static_cast<long double>(x), static_cast<long double>(y));
    }

    template<size_t m, size_t k = max(n, m)+1>
    Point<k> operator+(Point<m> const&o) const {
        Point<k> ret(x+o.x, y+o.y);
        return ret;
    }
    template<size_t m, size_t k = max(n, m)+1>
    Point<k> operator-(Point<m> const&o) const {
        Point<k> ret(x-o.x, y-o.y);
        return ret;
    }
    template<size_t m, size_t k = n+m>
    Point<k> operator*(Adaptive_Int<m> const&o) const {
        Point<k> ret(x*o, y*o);
        return ret;
    }
    template<typename T, typename = decltype(declval<coord_t>() / declval<T>())>
    Point operator/(T const&o) const {
        Point ret(x/o, y/o);
        return ret;
    }
    template<typename T, typename = decltype(declval<coord_t>() /= declval<T>())>
    Point& operator/=(T const&o) const {
        x/=o;
        y/=o;
        return *this;
    }

    template<size_t m>
    Point& operator+=(Unsafe_Wrapper<Point<m> const&> o){
        x+=o().x;
        x+=o().y;
        return *this;
    }
    template<size_t m>
    Point& operator-=(Unsafe_Wrapper<Point<m> const&> o){
        x+=o().x;
        x+=o().y;
        return *this;
    }
    template<size_t m, size_t k = n+m>
    Point<k> operator*=(Unsafe_Wrapper<Adaptive_Int<m> const&> o) const {
        x*=o;
        y*=o;
        return *this;
    }

    template<size_t m, size_t k = n+m+1>
    Adaptive_Int<k> dot(Point<m> const&o) const {
        Adaptive_Int<k> ret(x*o.x);
        ret+= make_unsafe(y*o.y);
        return ret;
    }
    template<size_t m, size_t k = n+m+1>
    Adaptive_Int<k> cross(Point<m> const&o) const {
        Adaptive_Int<k> ret(x*o.y);
        ret-= make_unsafe(y*o.x);
        return ret;
    }
    template<size_t k = 2*n+1>
    Adaptive_Int<k> norm_sq() const {
        Adaptive_Int<k> ret(x*x);
        ret+= make_unsafe(y*y);
        return ret;
    }

    template<size_t m>
    int comp_angular_180(Point<m> const&o) const {
        return o.cross(*this).sign();
    }
    template<size_t m, size_t k = n+m>
    int comp_angular_360(Point<m> const&o) const {
        bool low = is_nonneg_angle(), o_low = o.is_nonneg_angle();
        return low != o_low ? o_low-low : comp_angular_180(o);
    }
    template<size_t m>
    int comp_lexicographical(Point<m> const&o) const {
        const int c1 = x.comp(o.x);
        return c1 ? c1 : y.comp(o.y);
    }
    Point conj() const {
        return Point(x, -y);
    }

    template<size_t m, size_t k = n+m>
    Point<k> angle_sum(Point<m> const&o) const {
        return Point<k>(x*o.x - y*o.y, x*o.y + y*o.x);
    }
    template<size_t m, size_t k = n+m>
    Point<k> angle_diff(Point<m> const&o) const {
        return Point<k>(x*o.x + y*o.y, y*o.x - x*o.y);
    }
    template<size_t m>
    bool operator==(Point<m> const&o) const {
        return x == o.x && y == o.y;
    }
    template<size_t m>
    bool operator!=(Point<m> const&o) const {
        return !(operator==(o));
    }
    friend std::istream& operator>>(std::istream&in, Point &p){
        in >> p.x >> p.y;
        return in;
    }
    friend std::ostream& operator<<(std::ostream&o, Point const&p){
        return o << "(" << p.x << ", " << p.y << ")";
    }
private:
    bool is_nonneg_angle() const {
        const int A = y.comp(0);
        return A ? A>0 : x > 0;
    }
};

#ifdef DACIN_HASH_HPP
template<size_t n>
struct Dacin_Hash<Point<n> > {
    using pair_t = std::pair<typename Point<n>::coord_t const&, typename Point<n>::coord_t const&>;
    size_t operator()(Point<n> const&val) const {
        static Dacin_Hash<pair_t> h;
        return h(pair_t(val.x, val.y));
    }
};
#endif // DACIN_HASH_HPP

template<size_t n>
int ccw(Point<n> const&a, Point<n> const&b, Point<n> const&c){
    return - (b-a).comp_angular_180(c-a);
}


template<size_t n>
std::vector<Point<n>> convex_hull(std::vector<Point<n> > pts){
    std::sort(pts.begin(), pts.end(), [](Point<n> const&a, Point<n> const&b){return a.comp_lexicographical(b) < 0;});
    pts.erase(std::unique(pts.begin(), pts.end(), [](Point<n> const&a, Point<n> const&b){return a.comp_lexicographical(b) == 0;}), pts.end());

    std::vector<Point<n>> hull;
    for(size_t it=0;it<2;++it){
        const size_t old_size = hull.size();
        for(auto const&e:pts){
            while(hull.size() > old_size+1 && ccw(hull.rbegin()[1], hull.back(), e) <= 0){
                hull.pop_back();
            }
            hull.push_back(e);
        }
        if(hull.size() > 1) hull.pop_back();
        reverse(pts.begin(), pts.end());
    }
    return hull;
}


template<size_t n, size_t m, size_t k = max(n, m)+1>
std::vector<Point<k> > minkowski_sum(std::vector<Point<n>> a, std::vector<Point<m> > b){
    std::rotate(a.begin(), min_element(a.begin(), a.end(), [](Point<n> const&p1, Point<n> const&p2){return p1.comp_lexicographical(p2) < 0;}), a.end());
    std::rotate(b.begin(), min_element(b.begin(), b.end(), [](Point<m> const&p1, Point<m> const&p2){return p1.comp_lexicographical(p2) < 0;}), b.end());
    Point<k> last_dir(0, -1);
    std::vector<Point<k> > ret;
    if(a.size()>1) a.push_back(a.front());
    if(b.size()>1) b.push_back(b.front());
    for(size_t i=0,j=0;i<a.size() && j<b.size();){
        ret.push_back(a[i]+b[j]);
        if(i+1 == a.size()) ++j;
        else if(j+1 == b.size()) ++i;
        else {
            Point<n+1> d_a = a[i+1] - a[i];
            Point<m+1> d_b = b[j+1] - b[j];
            if(d_a.angle_diff(last_dir).comp_angular_360(d_b.angle_diff(last_dir)) < 0){
                last_dir = d_a;
                ++i;
            } else {
                last_dir = d_b;
                ++j;
            }
        }
    }
    if(ret.size() > 1 && ret.back() == ret.front()){
        ret.pop_back();
    }
    return ret;
}

template<size_t n, size_t k = 2*n+3>
Adaptive_Int<k> polygon_area_doubled(std::vector<Point<n> > const&poly){
    Adaptive_Int<k> ret(0);
    for(size_t i=0;i+1<poly.size();++i){
        ret+=make_unsafe(poly[i].cross(poly[i+1]));
    }
    ret+=make_unsafe(poly.back().cross(poly.front()));
    return ret;
}

template<size_t n>
bool segments_intersect(std::pair<Point<n>, Point<n> > const&s, std::pair<Point<n>, Point<n> > const&t){
    using int_t = Adaptive_Int<n>;
    // 1D bounding box test
    auto interval_intersect = [](int_t a, int_t b, int_t c, int_t d){
        if(a>b) swap(a, b);
        if(c>d) swap(c, d);
        return (a<=d && c<=b);
    };
    return interval_intersect(s.first.x, s.second.x, t.first.x, t.second.x)
        && interval_intersect(s.first.y, s.second.y, t.first.y, t.second.y)
        && ccw(s.first, s.second, t.first) * ccw(s.first, s.second, t.second) <= 0
        && ccw(t.first, t.second, s.first) * ccw(t.first, t.second, s.second) <= 0;
}

/// 1: inside, 0: ontop, -1: outside
template<size_t n>
int is_in_circumcircle(Point<n> const&a, Point<n> const&b, Point<n> const&c, Point<n> const&x){
    auto const A = a-x, B = b-x, C = c-x;
    auto const X = A.norm_sq(), Y = B.norm_sq(), Z = C.norm_sq();
    auto det = A.cross(B)*Z + B.cross(C)*X + C.cross(A)*Y;
    return det.sign() * ccw(a, b, c);
}


} // namespace dacin::geom

#endif // GEOM_2D_HPP
