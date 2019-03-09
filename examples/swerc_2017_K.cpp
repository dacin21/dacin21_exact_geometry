// https://swerc.eu/2017/problems/
// Problem K: 'Blowing candles'

#include "exact_geometry_all.hpp"

#include <bits/stdc++.h>
using namespace std;

constexpr int BITS = 29;

using dacin::geom::Point;

signed main(){
    #ifdef LOCAL_RUN
    freopen("in.txt", "r", stdin);
    #endif // LOCAL_RUN

    int n, r;
    while(cin >> n >> r){
        vector<Point<BITS> > poly(n), p2;
        for(auto &e:poly){
            cin >> e.x >> e.y;
        }
        poly = dacin::geom::convex_hull(poly);
        //for(auto &e:poly) cerr << e << "\n"; cerr << "\n";
        for(auto &e:poly){
            p2.emplace_back(-e.x, -e.y);
        }
        auto diff = dacin::geom::minkowski_sum(poly, p2);
        //for(auto &e:diff) cerr << e << "\n";
        diff.push_back(diff.front());
        double ans = 1e90;
        for(size_t i=0;i+1 < diff.size();++i){
            using Pt = dacin::geom::decay_t<decltype(diff[0])>;
            auto arr = dacin::geom::polygon_area_doubled(vector<Pt>{Pt(0, 0), diff[i], diff[i+1]});
            auto base = (diff[i+1]-diff[i]).norm_sq();
            //cerr << arr << " " << (double)arr << " / " << base << " " << (double)base << "\n";
            double cand = (double)arr/sqrt((double)base);
            ans = min(ans, cand);
        }
        cout << fixed << setprecision(15) << ans << "\n";
    }
}
