// Released under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007, see the LICENSE file.
// Copyright (C) 2018-2019 Daniel Rutschmann aka. dacin21

#ifndef GEOM_UTILITY_HPP
#define GEOM_UTILITY_HPP

#include <bits/stdc++.h>

namespace dacin::geom{

using std::integral_constant;
using std::declval;

#if __cplusplus <= 201103L
template<bool B, class T = void>
using enable_if_t = typename std::enable_if<B,T>::type;
template<bool B, class T, class F>
using conditional_t = typename std::conditional<B,T,F>::type;
template<class T>
using decay_t = typename std::decay<T>::type;
template<typename T>
constexpr T max(T const&a, T const&b){return (a < b) ? b : a;}
template<typename T>
constexpr T min(T const&a, T const&b){return (a > b) ? b : a;}
#else
using std::enable_if_t;
using std::conditional_t;
using std::decay_t;
using std::max;
using std::min;
#endif // __cplusplus

#if __cplusplus < 201703L
template<class T, class... Args>
constexpr bool is_constructible_v = std::is_constructible<T, Args...>::value;
template<class T>
constexpr bool is_integral_v = std::is_integral<T>::value;
template<class T, class U>
constexpr bool is_same_v = std::is_same<T, U>::value;
#else
using std::is_constructible_v;
using std::is_integral_v;
using std::is_same_v;
#endif // __cplusplus


template<size_t n>
struct Priority : Priority<n-1> {};
template<>
struct Priority<0>{};


template<typename T>
struct Unsafe_Wrapper {
    T val;
    Unsafe_Wrapper(T val_) : val(val_) {}
    T operator()() const {return val;}
};
template<typename T>
Unsafe_Wrapper<T const&> make_unsafe(T const& val) {
    return Unsafe_Wrapper<T const&>(val);
}

} // namespace dacin::geom

#endif // GEOM_UTILITY_HPP
