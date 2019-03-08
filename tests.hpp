// Released under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007, see the LICENSE file.
// Copyright (C) 2018-2019 Daniel Rutschmann aka. dacin21

#undef _GLIBCXX_DEBUG

#include "exact_geometry_all.hpp"

using namespace std;

namespace dacin::geom{

    template<size_t bits, typename T>
    void test_hulls_square_impl(uint64_t expected = 0, size_t log_lim = bits){
        cerr << "Running test test_hulls_square_impl " << bits << " " << log_lim << "\n";
        mt19937 rng(100531);
        auto get_rand = [&](T l, T r){return uniform_int_distribution<T>(l, r)(rng);};
        const int n = 3000;
        const int ITER = 1500;
        uint64_t ha = 0;
        for(int it=0;it<ITER;++it){
            vector<Point<bits> > p(n);
            const T lim = (T{1}<<(log_lim-1)) - T{1} + (T{1}<<(log_lim-1));
            for(auto &e:p){
                e.x = Adaptive_Int<bits>(get_rand(0, lim));
                e.y = Adaptive_Int<bits>(get_rand(0, lim));
            }
            auto ret = convex_hull(p);
            ha*=12347;
            ha+=ret.size();
        }
        if(expected == 0){
            cerr << bits << " " << typeid(T).name() << " : " << ha << "\n";
        } else {
            assert(ha == expected);
            cerr << "Test ok\n";
        }
    }
    void test_hulls_square(){
        test_hulls_square_impl<30, int>(3813651317602481189ull);
        test_hulls_square_impl<63, int64_t>(3813651317602481189ull, 30);
        test_hulls_square_impl<321, int>(3813651317602481189ull, 30);
        test_hulls_square_impl<31, int>(11956885618327015784ull);
        test_hulls_square_impl<63, int64_t>(11956885618327015784ull, 31);
        test_hulls_square_impl<321, int>(11956885618327015784ull, 31);
        test_hulls_square_impl<32, int64_t>(16151321917681398953ull);
        test_hulls_square_impl<321, int64_t>(16151321917681398953ull, 32);
        test_hulls_square_impl<61, int64_t>(4189760975843510352ull);
        test_hulls_square_impl<321, int64_t>(4189760975843510352ull, 61);
        test_hulls_square_impl<62, int64_t>(5248897970778511319ull);
        test_hulls_square_impl<321, int64_t>(5248897970778511319ull, 62);
        test_hulls_square_impl<63, int64_t>(7546482978952081156ull);
        test_hulls_square_impl<321, int64_t>(7546482978952081156ull, 63);
    }
    template<size_t bits, typename T>
    void test_circumcircle_impl(uint64_t expected = 0, size_t log_lim = bits){
        cerr << "Running test test_circumcircle " << bits << " " << log_lim << "\n";
        mt19937 rng(100531);
        auto get_rand = [&](T l, T r){return uniform_int_distribution<T>(l, r)(rng);};
        const int ITER = 15000;
        uint64_t ha = 0;
        for(int it=0;it<ITER;++it){
            array<Point<bits>, 4> p;
            const T lim = (T{1}<<(log_lim-1)) - T{1} + (T{1}<<(log_lim-1));
            for(auto &e:p){
                e.x = Adaptive_Int<bits>(get_rand(0, lim));
                e.y = Adaptive_Int<bits>(get_rand(0, lim));
            }
            auto ret = (uint32_t)is_in_circumcircle(p[0], p[1], p[2], p[3]);
            ha*=12347;
            ha+=ret;
        }
        if(expected == 0){
            cerr << bits << " " << typeid(T).name() << " : " << ha << "\n";
        } else {
            assert(ha == expected);
            cerr << "Test ok\n";
        }
    }

    void test_circumcircle(){
        test_circumcircle_impl<20, int>(4884053424985750744ull);
        test_circumcircle_impl<321, int>(4884053424985750744ull, 20);
        test_circumcircle_impl<21, int>(9294910920041591446ull);
        test_circumcircle_impl<321, int>(9294910920041591446ull, 21);
        test_circumcircle_impl<22, int>(6899345392128341650ull);
        test_circumcircle_impl<321, int>(6899345392128341650ull, 22);
        test_circumcircle_impl<62, int64_t>(2678531202223027034ull);
        test_circumcircle_impl<321, int64_t>(2678531202223027034ull, 62);
        test_circumcircle_impl<63, int64_t>(5179094771316369856ull);
        test_circumcircle_impl<321, int64_t>(5179094771316369856ull, 63);

    }

} // namespace dacin::geom



void run_tests(){
    cerr << "Running all tests\n";
    //dacin::geom::test_hulls_square();
    dacin::geom::test_circumcircle();

    cerr << "Done with all tests\n";
}
