#pragma GCC optimize("O3")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,tune=native") // yandex


#include "../exact_geometry_all.hpp"

using namespace std;
using namespace dacin::geom;

template<bool enabled>
struct Debug{
    template<typename S, typename T = void> struct Tag_Printable : false_type {};
    template<typename S> struct Tag_Printable<S, decltype((void)(cerr << declval<S>()))> : true_type {};
    template<typename S, typename T = void> struct Tag_Iterable: false_type {};
    template<typename S> struct Tag_Iterable<S, decltype((void)(begin(declval<S>()), end(declval<S>())))> : true_type {};

    template<typename T, size_t N> struct Tuple_Printer{
        template<typename S>
        static S& print(S& stream, T const&t){
            return Tuple_Printer<T, N-1>::print(stream, t) << ", " << get<N>(t);
        }
    };
    template<typename T> struct Tuple_Printer<T, 0>{
        template<typename S>
        static S& print(S& stream, T const&t){
            return stream << get<0>(t);
        }
    };

    template<typename T, typename... Args>
    Debug& print(T const&x, true_type, Args...){
        #ifdef LOCAL_RUN
        if(enabled){
            cerr << boolalpha << x;
        }
        #endif // LOCAL_RUN
        return *this;
    }
    template<typename T>
    Debug& print(T const&x, false_type, true_type){
        *this << "[";
        bool first = true;
        for(auto &e:x){
            if(!first) *this << ", ";
            *this << e;
            first = false;
        }
        return *this << "]";
    }
    template<typename S, typename T>
    Debug& print(pair<S, T> const&x, false_type, false_type){
        return *this << "(" << x.first << ", " << x.second << ")";
    }
    template<typename... Args>
    Debug& print(tuple<Args...> const&t, false_type, false_type){
        *this << "(";
        return Tuple_Printer<decltype(t), sizeof...(Args)-1>::print(*this, t) << ")";
    }
    template<typename T>
    Debug& operator<<(T const&x){
        return print(x, Tag_Printable<T>{}, Tag_Iterable<T>{});
    }
};
// Debug<true> debug;
 Debug<false> debug; // disable debug printing
#define named(x) string(#x) << " : " <<  x


// nicer implementation with struct
template<class Segtree_Data>
struct Segment_Tree{
	using T = typename Segtree_Data::node_t;
	int n;
	vector<T>data;
	Segment_Tree(int _n):n(_n), data(2*n, Segtree_Data::node_ne()){
		for(int i=n-1;i>=0;--i) data[i] = Segtree_Data::merge_nodes(data[i<<1], data[i<<1|1]);
	}
	Segment_Tree(vector<T> const&base):n(base.size()), data(2*n, Segtree_Data::node_ne()){
		copy(base.begin(), base.end(), data.begin()+n);
		for(int i=n-1;i>=0;--i) data[i] = Segtree_Data::merge_nodes(data[i<<1], data[i<<1|1]);
	}
	void update(int pos, typename Segtree_Data::update_t const&val){
		for(Segtree_Data::update_node(data[pos+=n], val);pos>>=1;){
			data[pos] = Segtree_Data::merge_nodes(data[pos<<1], data[pos<<1|1]);
		}
	}
	T query(int l, int r)const{
		T retL = Segtree_Data::node_ne(), retR = Segtree_Data::node_ne();
		for(l+=n, r+=n;l<r;l>>=1, r>>=1){
			if(l&1) retL = Segtree_Data::merge_nodes(retL, data[l++]);
			if(r&1) retR = Segtree_Data::merge_nodes(data[--r], retR);
		}
		return Segtree_Data::merge_nodes(retL, retR);
	}
};

constexpr int B = 21;
using Pt = Point<B>;
using Dir = Point<B+3>;


Dir cur_dir;
size_t merges = 0;

template<typename Agg>
struct Radial_Heap{
    using node_t =  tuple<Dir, Agg, int>;
    static node_t node_ne() {return make_tuple(Dir{}, Agg{}, -1);}

    Radial_Heap(int n_) : n(n_), data(2*n+2, node_ne()){}

    node_t merge_nodes(node_t const&l, node_t const&r){
        ++merges;
        if(get<2>(l) == -1 || !get<0>(l)) return r;
        if(get<2>(r) == -1 || !get<0>(r)) return l;
        if(get<2>(r) < get<2>(l)) return r;
        if(get<2>(r) > get<2>(l)) return l;
        return get<0>(l).comp_angular_360(get<0>(r)) <= 0 ? l : r;
    }
    void recalc(int i){
        for(int j=n+i;j>>=1;){
            data[j] = merge_nodes(data[2*j], data[2*j+1]);
        }
    }
    template<bool allow_cur_time>
    void update(int i, Dir dir, Agg agg){
        int rot = cur_rot;
        const int c = dir.comp_angular_360(cur_dir);
        if(c == -1 || (!allow_cur_time && c == 0)){
            ++rot;
        }
        node_t new_val = make_tuple(dir, agg, rot);
        if(data[n+i] != new_val){
            data[n+i] = new_val;
            recalc(i);
        }
    }
    void clear(int i){
        data[i+n] = node_ne();
    }
    bool empty(){
        return get<2>(top()) == -1;
    }
    pair<Dir, Agg> top(){
        cur_rot = max(cur_rot, get<2>(data[1]));
        return make_pair(get<0>(data[1]), get<1>(data[1]));
    }

    int n;
    int cur_rot = 0;
    vector<node_t> data;
};

struct Segtreedata{
    typedef pair<Dir, int> node_t;
    typedef node_t update_t;
    static  node_t node_ne() {return make_pair(Dir{}, -1);}
    static node_t merge_nodes(node_t const&left, node_t const&right){
        ++merges;
        if(left.second == -1) return right;
        if(right.second == -1) return left;
        const int c = (left.first.angle_diff(cur_dir)).comp_angular_360(right.first.angle_diff(cur_dir));
        return c < 0 ? left : c > 0 ? right : left;
    }
    static void update_node(node_t &node, update_t const&update){
        if(update.second != -1) node = update;
    }
};


struct Fen{
    Fen(int n_):n(n_+5), data(n){}
    int q(int x){
        int ret = 0;
        for(++x;x;x-=x&-x){
            ret+=data[x];
        }
        return ret;
    }
    void u(int x, int v){
        for(++x;x<n;x+=x&-x){
            data[x]+=v;
        }
    }
    int n;
    vector<int> data;
};


template<size_t n>
Point<n> ccw90(Point<n> const&p){
    return Point<n>(-p.y, p.x);
}



void solve(){
    int n, k;
    cin >> n >> k;
    vector<Point<B>> v;
    vector<pair<int, int> > v_orig;
    for(int i=0;i<n;++i){
        int x, y;
        cin >> x >> y;
        v.push_back(Point<B>(x, y));
        v_orig.emplace_back(x, y);
    }


    auto run_it = [&](int output_id){
        vector<int> ord(n);
        string side(n, 'L');
        iota(ord.begin(), ord.end(), 0);
        sort(ord.begin(), ord.end(), [&v](int const&a, int const&b){
            return make_pair(v[a].x*v[a].x, v[a].norm_sq()) < make_pair(v[b].x*v[b].x, v[b].norm_sq());
        });
        for(int i=0;i<n;++i){
            side[i] = v[i].x < 0 ? 'L' : 'R';
        }
        cur_dir = Dir(0, 1);
        debug << named(ord) << "\n";

        // 0..n-1 for points, n .. 2n-1 for adjacent swaps
        //Segment_Tree<Segtreedata> st(n + n-1);
        Radial_Heap<int> st(n + n-1);

        Fen fen_x(n), fen_y(n);

        auto fen_up = [&](int i, int f){
            if(side[ord[i]] == 'L') f=-f;
            fen_x.u(i, f*v_orig[ord[i]].first);
            fen_y.u(i, f*v_orig[ord[i]].second);
        };


        // index in p
        auto re_single = [&](int i){
            //if(!v[i].norm_sq()) return;
            if(side[i] == 'L'){
                st.update<true>(i, v[i], i);
            } else {
                st.update<true>(i, -v[i], i);
            }
        };

        // index in ord
        // find event for b to be closer than a
        auto re_adj = [&](int i){
            if(i<0 || i>=n-1) return;
            assert(0 <= i && i+1 < n);
            const int a = ord[i], b = ord[i+1];
            auto pb = v[b], pa = v[a];
            if(pb.cross(pa).sign() == 0){
                auto dir = [&](){
                    Dir d1 = ccw90(pa - pb);
                    Dir d2 = ccw90(pb - pa);
                    if(d1.comp_angular_360(cur_dir) == 0) return d2;
                    if(d2.comp_angular_360(cur_dir) == 0) return d1;
                    if((d1.angle_diff(cur_dir).comp_angular_360(d2.angle_diff(cur_dir))) < 0){
                        return d1;
                    } else {
                        return d2;
                    }
                }();
                st.update<true>(n+i, dir, -1-i);
            } else {
                auto dir = [&](){
                    if(side[a] == 'L' && side[b] == 'L'){
                        return pb - pa;
                    } else if(side[a] == 'L' && side[b] == 'R'){
                        return -pa - pb;
                    } else if(side[a] == 'R' && side[b] == 'L'){
                        return pa + pb;
                    } else if(side[a] == 'R' && side[b] == 'R'){
                        return pa - pb;
                    } else assert(0);
                }();
                st.update<true>(n+i, dir, n+i);
            }
        };

        int ret = 0;
        long double ans = -1e90;
        Dir ret_dir;
        auto cand_dir = [&](Dir const&d, int id){
            if(!d.norm_sq()) return;
            Dir p_sum(fen_x.q(k-1), fen_y.q(k-1));
            auto numer = (long double) p_sum.cross(d);
            auto denom = sqrtl((long double)d.norm_sq());
            auto val = numer / denom;
            debug << d << " " << val << "\n";
            //debug << named(st.data) << "\n";
            debug << ord << "\n";
            if(val > ans){
                ans = val;
                ret = id;
                ret_dir = d;
            }
            /*if(id == output_id){
                cout << fixed << setprecision(20) << val << "\n";
                cout << (long double)d.x << " " << (long double)d.y << "\n";
                for(int i=0;i<k;++i){
                    cout << ord[i]+1 << " ";
                }
                cout << "\n";

            }*/
        };

        for(int i=0;i<n;++i){
            fen_up(i, 1);
        }



        for(int i=0;i<n;++i){
            re_adj(i);
            re_single(i);
        }
        const int STEPS = n*n + 4*n + 20;
        int id = 0;
        for(int it=0;it<STEPS;++it){
            //if(it%n == 0) cerr << it/n << " : " << merges << "\n";
            auto ev = st.top();
            debug << "\n" << named(ev) << "\n";
            debug << named(st.data) << "\n";
            //assert(ev.second != -1);
            int i = ev.second;
            auto ev_dir = ev.first;
            cand_dir(cur_dir, ++id);
            cand_dir(ev_dir, ++id);
            Dir mid ( -fen_y.q(k-1), fen_x.q(k-1));
            if( (mid.angle_diff(cur_dir)).comp_angular_360(ev_dir.angle_diff(cur_dir)) <= 0){
                cand_dir(mid, ++id);
            }

            cur_dir = ev_dir;
            if (i<0){
                i = -1-i;
                re_adj(i);

            } else if(i < n){
                // point switches side
                const int j = find(ord.begin(), ord.end(), i) - ord.begin();
                fen_up(j, -1);
                if(side[i] == 'L'){
                    side[i] = 'R';
                } else {
                    side[i] = 'L';
                }
                fen_up(j, 1);
                re_single(i);
                re_adj(j);
                re_adj(j-1);
            } else {
                i-=n;
                // swap v[i] with v[i+1] after dir
                fen_up(i, -1);
                fen_up(i+1, -1);
                swap(ord[i], ord[i+1]);
                fen_up(i, 1);
                fen_up(i+1, 1);
                for(int j=-1;j<=1;++j){
                    const int k = i+j;
                    if(0<=k && k<n){
                        re_single(ord[k]);
                    }
                    re_adj(k);
                }
            }

        }
        auto sort_dir = ret_dir.ccw90();
        sort(ord.begin(), ord.end(), [&](int const&a, int const&b){
            return v[a].dot(sort_dir).abs() < v[b].dot(sort_dir).abs();
        });
        cout << fixed << setprecision(20) << ans << "\n";
        cout << (long double)ret_dir.x << " " << (long double)ret_dir.y << "\n";
        for(int i=0;i<k;++i){
            cout << ord[i]+1 << " ";
        }
        cout << "\n";

        debug << ans << "\n";
        return ret;
    };
    int opt_id = run_it(-1);
    //run_it(opt_id);
    cerr  << merges << "\n";
}

signed gen(int T){
    mt19937 rng(43151);
    auto get_rand = [&](int64_t l, int64_t r){
        return uniform_int_distribution<int64_t>(l, r)(rng);
    }; (void) get_rand;
    auto get_double = [&](double l, double r){
        return uniform_real_distribution<double>(l, r)(rng);
    };  (void) get_double;
    ofstream o("gen.txt");
    o << T << "\n";
    for(int cas=0;cas<T;++cas){
        /// GEN HERE
        const int n = 2000;
        const int k = get_rand(n/3, 2*n/3);
        o << n << " " << k << "\n";

        for(int i=0;i<n;++i){
            const int x = get_rand(-1e6, 1e6);
            const int y = get_rand(-1e6, 1e6);
            o << x << " " << y << "\n";
        }

        o << "\n";
    }
    o << endl;
    o.close();
    return 0;
}

signed main()
{
    #ifdef LOCAL_RUN
    freopen("examples/in.txt", "r", stdin);
    //freopen("out.txt", "w", stdout);
    cin.exceptions(ios::badbit | ios::eofbit | ios::failbit);
    cin.tie(nullptr); ios_base::sync_with_stdio(false);
    int TTT; cin >> TTT;
    if(TTT < 0) return gen(-TTT);
    while(TTT--){
    #else
    cin.tie(nullptr); ios_base::sync_with_stdio(false);
    #endif // LOCAL_RUN

    solve();

    #ifdef LOCAL_RUN
    cout << flush;
    }
    #endif // LOCAL_RUN
    return 0;
}
