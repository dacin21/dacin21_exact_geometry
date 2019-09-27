// Released under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007, see the LICENSE file.
// Copyright (C) 2018-2019 Daniel Rutschmann aka. dacin21

#ifndef GEOM_UTILITY_HPP
#define GEOM_UTILITY_HPP

#ifdef RELEASE
#undef _GLIBCXX_DEBUG
#endif // RELEASE


#ifndef __SIZEOF_INT128__
static_assert(false, "__int128 not found");
#endif // __SIZEOF_INT128__

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

template<typename S, typename... Types>
struct is_one_of : std::false_type{};
template<typename S, typename T, typename... Types>
struct is_one_of<S, T, Types...> : std::integral_constant<bool, is_same_v<S, T> || is_one_of<S, Types...>::value> {};

using uint128_t = unsigned __int128;
using int128_t = __int128;

std::ostream& operator<<(std::ostream&o, uint128_t x){
    if(x < 10){
        return o << x;
    } else {
        return o << x/10 << (int)x%10;
    }
}
std::ostream& operator<<(std::ostream&o, int128_t x){
    if(x < 0){
        return o << '-' << -(uint128_t)x;
    }
    return o << (uint128_t)x;
}

} // namespace dacin::geom

#endif // GEOM_UTILITY_HPP
