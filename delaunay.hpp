#ifndef DELAUNAY_HPP
#define DELAUNAY_HPP

#include "geom_2d.hpp"

namespace dacin::geom{

template<typename point_t>
struct Delaunay_Face{
    std::array<point_t, 3> corners;
    std::array<Delaunay_Face*, 3> adj {nullptr, nullptr, nullptr};
    std::vector<int> bucket;
    Delaunay_Face(): corners{} {}
    Delaunay_Face(point_t const&a, point_t const&b, point_t const&c) : corners({a, b, c}) {}
    friend std::ostream& operator<<(std::ostream&o, const Delaunay_Face&ff){
        o << "Face: " << &ff << "\n";
        return o << ff.corners[0] << "\n" << ff.corners[1] << "\n" << ff.corners[2] << "\n";
    }
};
/**
 *  Delaunay triangulation in expected O(n log n)
 *
 *  *Warning*: All points have to be contained in the infinite triangle.
 *  Hence n has to be at least 2 + log_2(max_coord)
 *
 *  Does not work with duplicate points, but works in degenerate cases.
 *
 */
template<size_t n>
class Delaunay{
public:
    using point_t = Point<n>;
    using coord_t = typename Point<n>::coord_t;
    using Face = Delaunay_Face<point_t>;

    static const coord_t INF;
    static const point_t inf_n, inf_sw, inf_se;

    static bool is_infinite(point_t const&p){
        return (p.x == -INF || p.x == INF || p.y == INF || p.y == -INF);
    }
    static bool is_infinite(Face const&f){
        return is_infinite(f.corners[0]) || is_infinite(f.corners[1]) || is_infinite(f.corners[2]);
    }

    Delaunay(){}

    std::vector<Face>& triangulate(std::vector<point_t> const&p){
        points = p;
        int N = p.size();
        faces.reserve(3*N);
        //std::random_shuffle(points.begin(), points.end());
		// start with super triangle that contains all points
        locateFace = new (getFreeFace()) Face(inf_n, inf_sw, inf_se);
        locateFace->bucket.resize(N);
        std::iota(locateFace->bucket.begin(), locateFace->bucket.end(), 0);
        vertex_location.resize(N, locateFace);
        // incremental construction
        for(int i=0;i<N;++i){
            Face* place = vertex_location[i];
            split(place, i);
        }
        compress_faces();
        return faces;
    }

private:
    static coord_t get_inf(){
        coord_t ret(1);
        const size_t exp = n-1;
        for(size_t i=1;i<=exp;i<<=1){
            ret*=make_unsafe(ret);
            if(exp&i){
                ret*=make_unsafe(coord_t(2));
            }
        }
        return ret;
    }
    static bool has_to_flip(Face const&f, point_t const&p){
        point_t const&A = f.corners[0], B = f.corners[1], C = f.corners[2];
        // infinite point
        if(is_infinite(A)) return ccw(B, C, p) > 0;
        if(is_infinite(B)) return ccw(C, A, p) > 0;
        if(is_infinite(C)) return ccw(A, B, p) > 0;
        // degenerate face
        if(ccw(A, B, C) == 0) return ccw(A, B, p) + ccw(B, C, p) + ccw(C, A, p) > 0;
        // infinite point is never in finite face
        if(is_infinite(p)) return false;
        return is_in_circumcircle(A, B, C, p) > 0;
    }
    static int get_other_dir(Face*f, int dir, Face*old_f){
        Face*other = f->adj[dir];
        if(other->adj[0] == old_f) return 0;
        if(other->adj[1] == old_f) return 1;
        if(other->adj[2] == old_f) return 2;
        return -1;
    }
    static void link_face(Face*f, int dir, Face*old){
        if(f->adj[dir] == nullptr) return;
        int other_dir = get_other_dir(f, dir, old);
        f->adj[dir]->adj[other_dir] = f;
    }

    Face* locate(point_t const& p, Face* cur){
        asser(cur);
        for(int i=0;i<3;++i){
            if(ccw(cur->corners[i], cur->corners[(i+1)%3], p)<0){
                cur = cur->adj[i];
                i=-1;
            }
        }
        return cur;
    }
    Face* getFreeFace(){
        faces.emplace_back();
        return &(faces.back());
    }
    void link_bucket(Face*f){
        if(!f->bucket.empty()){
            vertex_location[f->bucket.front()] = f;
        }
    }
    void check_flips(Face*f, int dir){
        if(f->adj[dir] == 0) return;
        int other_dir = get_other_dir(f, dir, f);
        Face*o = f->adj[dir];
        if(has_to_flip(*f, o->corners[other_dir])){
            f->corners[(dir+1)%3] = o->corners[other_dir];
            o->corners[(other_dir+1)%3] = f->corners[dir];
            f->adj[dir] = o->adj[(other_dir+2)%3];
            o->adj[(other_dir+2)%3] = f;
            o->adj[other_dir] = f->adj[(dir+2)%3];
            f->adj[(dir+2)%3] = o;
            link_face(f, dir, o);
            link_face(o, other_dir, f);
            std::vector<int> tmp(f->bucket.size() + o->bucket.size());
            std::merge(f->bucket.begin(), f->bucket.end(), o->bucket.begin(), o->bucket.end(), tmp.begin());
            f->bucket.clear();
            o->bucket.clear();
            for(int e:tmp){
                if(ccw(f->corners[dir], f->corners[(dir+1)%3], points[e])>0){
                    f->bucket.push_back(e);
                } else {
                    o->bucket.push_back(e);
                }
            }
            link_bucket(f);
            link_bucket(o);
            check_flips(f, (dir)%3);
            check_flips(f, (dir+1)%3);
            check_flips(o, (other_dir)%3);
            check_flips(o, (other_dir+1)%3);
        }
    }
    void split(Face*a, int point_index){
        Face*b = new (getFreeFace()) Face(a->corners[0], a->corners[1], points[point_index]);
        Face*c = new (getFreeFace()) Face(a->corners[1], a->corners[2], points[point_index]);
        a->corners[1] = points[point_index];
        b->adj = {c, a, a->adj[2]};
        c->adj = {a, b, a->adj[0]};
        a->adj = {c, a->adj[1], b};

        link_face(b, 2, a);
        link_face(c, 2, a);
        link_face(a, 1, a);
        std::vector<int> tmpBuck;
        tmpBuck.swap(a->bucket);
        for(int e:tmpBuck){
            if(e==point_index) continue;
            if(ccw(b->corners[1], b->corners[2], points[e])>=0 && ccw(b->corners[2], b->corners[0], points[e])>=0){
                b->bucket.push_back(e);
            } else if(ccw(c->corners[1], c->corners[2], points[e])>=0 && ccw(c->corners[2], c->corners[0], points[e])>=0){
                c->bucket.push_back(e);
            } else {
                a->bucket.push_back(e);
            }
        }
        link_bucket(a);
        link_bucket(b);
        link_bucket(c);
        check_flips(a, 1);
        check_flips(b, 2);
        check_flips(c, 2);
    }
    void compress_faces(){
        std::vector<Face> retFaces;
        retFaces.reserve(faces.size()+1);
        std::unordered_map<Face*, Face*> decode;
        decode.reserve(faces.size());
        for(auto &e:faces){
            if(!decode.count(&e)){
                retFaces.push_back(e);
                decode[&e] = &(retFaces.back());
            }
        }
        for(auto &e:retFaces){
            for(int i=0;i<3;++i){
                e.adj[i] = decode[e.adj[i]];
            }
        }
        faces.swap(retFaces);
        locateFace = decode[locateFace];
    }

    std::vector<Face> faces;
    std::vector<point_t> points;
    std::vector<Face*> vertex_location;
    Face* locateFace = 0;
};

template<size_t n> const typename Delaunay<n>::coord_t Delaunay<n>::INF = Delaunay<n>::get_inf();
template<size_t n> const typename Delaunay<n>::point_t Delaunay<n>::inf_n(typename Delaunay<n>::coord_t(0), Delaunay<n>::INF);
template<size_t n> const typename Delaunay<n>::point_t Delaunay<n>::inf_sw(-Delaunay<n>::INF, -Delaunay<n>::INF);
template<size_t n> const typename Delaunay<n>::point_t Delaunay<n>::inf_se{Delaunay<n>::INF, -Delaunay<n>::INF};

#ifdef DACIN_HASH_HPP
template<typename T>
struct Dacin_Hash<Delaunay_Face<T> >{
    size_t operator()(Delaunay_Face<T> const&face) const {
        using Point = T;
        using ppp = std::pair<Point const&, Point const&>;
        using ppppp = std::pair<Point const&, ppp>;
        static Dacin_Hash<ppppp> h;
        return h(ppppp(face.corners[0], ppp(face.corners[1], face.corners[2])));
    }

};
#endif // DACIN_HASH_HPP


} // namespace dacin::geom

#endif // DELAUNAY_HPP
