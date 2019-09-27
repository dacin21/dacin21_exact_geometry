// https://codejam.withgoogle.com/2018/challenges/0000000000007766/dashboard/
// Problem E: 'The Cartesian Job'

#include "../exact_geometry_all.hpp"

#include <bits/stdc++.h>
using namespace std;

using coord_t = int64_t;

using dacin::geom::Point;


template<size_t n>
Point<2*n+2> angle(Point<n> const&base, Point<n> const&a, Point<n> const&b){
    return (b-base).angle_diff(a-base);
}

using ll = long long;
using ull = unsigned long long;
using fl = long double;

template<typename S, typename T>
void xmin(S&a, T const&b){if(b<a) a=b;}
template<typename S, typename T>
void xmax(S&a, T const&b){if(b>a) a=b;}

using namespace std;

constexpr size_t BITS = 30;
constexpr size_t LBITS = 62;

const Point<BITS> seg_base(0, 0), seg_top(0, 1000);
const Point<LBITS> mid_angle(-1, 0);

const fl pipi = 2 * 3.1415926535897932384626;


template<size_t n>
fl angfl(Point<n> const&p){
    fl tmp = atan2((double)p.x, (double)p.first);
    if(tmp < -1e-3){
        tmp+=pipi;
    }
    return tmp;
}

template<size_t n>
bool angle_comp(Point<n> const&p, Point<n> const&q){
    return p.comp_angular_360(q) < 0;
}
template<size_t n>
bool angle_eq(Point<n> const&p, Point<n> const&q){
    return p.comp_angular_360(q) == 0;
}

void solve(){
    /// SOLVE HERE
    int n;
    cin >> n;
    vector<pair<Point<BITS>, Point<BITS>> > pts(n);
    for(auto &e:pts){
        cin >> e.first >> e.second;
    }
    // remove degenerate lasers
    pts.erase(remove_if(pts.begin(), pts.end(), [](pair<Point<BITS>, Point<BITS>> const&p){return p.first.x == 0;}), pts.end());
    n = pts.size();
    vector<pair<Point<LBITS>, Point<LBITS>> > inters(n);
    Point<LBITS> angle_left(1, 0), angle_right(-1, 0);
    for(int i=0;i<n;++i){
        // compute time frame
        auto &l = inters[i].first, &r = inters[i].second;
        l = angle(pts[i].first, pts[i].second, seg_base);
        r = angle(pts[i].first, pts[i].second, seg_top);
        if(pts[i].first.x > 0) swap(l, r);
        //cerr << l << " | " << r << "\n";
        //cerr << angfl(l) << " / " << angfl(r) << "\n";
        // normalize to left side;
        if(angle_comp(l, r)){ // inner
            if(!angle_comp(l, r.conj())){
                l = l.conj();
                r = r.conj();
                swap(l, r);
            }
        } else { // outer
            if(angle_comp(l, r.conj()) || (r.y == 0 && r.x > 0)){
                l = l.conj();
                r = r.conj();
                swap(l, r);
            }
        }
        //cerr << angfl(l) << " / " << angfl(r) << "\n";
        // cut outer stuff
        if(angle_comp(mid_angle, l)){
            auto const cand = l.conj();
            if(angle_comp(angle_left, cand)) angle_left = cand;
            l = angle_left;
        }
        if(angle_comp(mid_angle, r)){
            auto const cand = r.conj();
            if(angle_comp(cand, angle_right)) angle_right = cand;
            r = angle_right;
        }
        //cerr << angfl(l) << " / " << angfl(r) << "\n";
        //cerr << angfl(angle_left) << " // " << angfl(angle_right) << "\n";
    }
    for(auto &e:inters){
        if(angle_comp(e.first, angle_left)) e.first = angle_left;
        if(angle_comp(angle_right, e.second)) e.second = angle_right;
    }
    inters.erase(remove_if(inters.begin(), inters.end(), [](pair<Point<LBITS>, Point<LBITS>> const&p){return !angle_comp(p.first, p.second);}), inters.end());
    n = inters.size();

    vector<Point<LBITS>> pool;
    for(auto const&e:inters){
        pool.push_back(e.first);
        pool.push_back(e.second);
    }
    if(!angle_comp(angle_left, angle_right)){
        cout << "0.0\n";
        return;
    }
    pool.push_back(angle_left);
    pool.push_back(angle_right);
    sort(pool.begin(), pool.end(), angle_comp<LBITS>);
    pool.erase(unique(pool.begin(), pool.end(), angle_eq<LBITS>), pool.end());
    assert(angle_eq(pool.front(), angle_left));
    assert(angle_eq(pool.back(), angle_right));
    auto deco = [&](Point<LBITS> const&p){
        return lower_bound(pool.begin(), pool.end(), p, angle_comp<LBITS>) - pool.begin();
    };
    const int m = pool.size();

    /*cerr << "pool: ";
    for(auto const&e:pool) cerr << angfl(e) << " ";
    cerr << "\n";*/

    vector<pair<int, int> > covers(n);
    for(int i=0;i<n;++i){
        covers[i].first = deco(inters[i].first);
        covers[i].second = deco(inters[i].second);
        //cerr << covers[i].first << " - " << covers[i].second << "\n";
    }
    sort(covers.begin(), covers.end());
    vector<fl> dp(m, 0.0), dp2(m, 0.0);
    int last = 0;
    dp[0] = 1.0;
    for(auto const&e:covers){
        fill(dp2.begin(), dp2.end(), fl{0});
        for(int i=last+1;i<m;++i) assert(dp[i] < 1e-9);

        // try same
        if(e.first <= last){
            for(int i=0;i<m;++i){
                dp2[i]+=dp[i]* 0.5;
            }
        }
        // try different
        if(e.second < last){
            for(int i = e.first;i<e.second;++i){
                dp2[e.second]+=dp[i]* 0.5;
            }
            for(int i = e.second;i<m;++i){
                dp2[i]+=dp[i]* 0.5;
            }
        } else {
            for(int i = e.first;i<last;++i){
                dp2[last]+=dp[i]* 0.5;
            }
            for(int i = last;i<m;++i){
                dp2[i]+=dp[i]* 0.5;
            }
            last = e.second;
        }
        dp.swap(dp2);
        // avoid denormal floats
        /*for(auto &e:dp){
            if(e < 1e-30) e = 0.0;
        }*/
        //for(auto const&f:dp) cerr << f << " ";
        //cerr << "\n";
    }
    fl ans = 1.0;
    if(last == m-1){
        ans = fl{1} - dp[m-1];
    }
    cout << ans << "\n";



}

signed gen(int T){
    mt19937 rng(43151);
    auto get_rand = [&](int64_t l, int64_t r){
        return uniform_int_distribution<int64_t>(l, r)(rng);
    }; (void) get_rand;
    auto get_double = [&](double l, double r){
        return uniform_real_distribution<double>(l, r)(rng);
    }; (void) get_double;
    ofstream o("gen.txt");
    o << T << "\n";
    for(int cas=0;cas<T;++cas){
        /// GEN HERE

        o << "\n";
    }
    o << endl;
    o.close();
    return 0;
}

int main()
{
    #ifdef LOCAL_RUN
    freopen("inE.txt", "r", stdin);
    //freopen("outX.txt", "w", stdout);
    #endif // LOCAL_RUN
    cin.tie(0); cout.tie(0);
    ios_base::sync_with_stdio(false);
    int TTT; cin >> TTT;
    cout << fixed << setprecision(15);
    for(int cas = 1;cas<=TTT;++cas){
        cout << "Case #" << cas << ": ";

        solve();

        //cout << "\n";
        #ifdef LOCAL_RUN
        cout << flush;
        #endif // LOCAL_RUN
    }
    return 0;
}


