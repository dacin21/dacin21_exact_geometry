// Released under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007, see the LICENSE file.
// Copyright (C) 2018-2019 Daniel Rutschmann aka. dacin21

#ifndef ADAPTIVE_INT_HPP
#define ADAPTIVE_INT_HPP

#include "geom_utility.hpp"
#include "bignum_fixedsize_signed.hpp"

namespace dacin::geom{

template<size_t n>
class Adaptive_Int{
public:
    using backend_t = conditional_t< n <= 31, int32_t, conditional_t<n <= 63, int64_t, Bigint_Fixedsize_Signed<std::max<size_t>(1, n/32+1)> > >;

    template<typename T>
    struct is_adaptive_int : std::false_type{};
    template<size_t m>
    struct is_adaptive_int<Adaptive_Int<m> > : std::true_type{};

    template<typename T, typename S>
    struct has_sign : std::false_type{};
    template<typename T>
    struct has_sign<T, decltype(declval<T>().sign())> : std::true_type{};

    template<typename S, typename T, typename SFINAE>
    struct can_comp : std::false_type{};
    template<typename S, typename T>
    struct can_comp<S, T, decltype(declval<S>().comp(declval<T>))> : std::true_type{};

    static Adaptive_Int ZERO;

    Adaptive_Int() : value() {}

    template<typename S, typename = enable_if_t<is_constructible_v<backend_t, S>>>
    explicit Adaptive_Int(S const&o) : value(o) {}

    template<size_t m, typename = enable_if_t<m <= m> >
    explicit Adaptive_Int(Adaptive_Int<m> const&o) : value(o.get_cvalue()) {}

    backend_t& get_value(){return value;}
    const backend_t& get_cvalue() const {return value;}

    template<size_t m, size_t k = std::max(n, m)+1>
    Adaptive_Int<k> operator+(Adaptive_Int<m> const&o) const {
        Adaptive_Int<k> ret(*this);
        ret+= make_unsafe(o);
        return ret;
    }
    template<size_t m, size_t k = std::max(n, m)+1>
    Adaptive_Int<k> operator-(Adaptive_Int<m> const&o) const {
        Adaptive_Int<k> ret(*this);
        ret-= make_unsafe(o);
        return ret;
    }
    template<size_t m, size_t k = n+m>
    Adaptive_Int<k> operator*(Adaptive_Int<m> const&o) const {
        Adaptive_Int<k> ret(*this);
        ret*= make_unsafe(o);
        return ret;
    }
    template<size_t m>
    Adaptive_Int& operator+=(Unsafe_Wrapper<Adaptive_Int<m> const&> o) {
        value+= o().get_cvalue();
        return *this;
    }
    template<size_t m>
    Adaptive_Int& operator-=(Unsafe_Wrapper<Adaptive_Int<m> const&> o) {
        static_assert(m <= n);
        value-= o().get_cvalue();
        return *this;
    }
    template<size_t m>
    Adaptive_Int& operator*=(Unsafe_Wrapper<Adaptive_Int<m> const&> o) {
        value*= o().get_cvalue();
        return *this;
    }

    template<typename T, typename = decltype(declval<backend_t>() / declval<T>())>
    Adaptive_Int operator/(T const&o) const {
        Adaptive_Int ret(value / o);
        return ret;
    }    template<typename T, typename = decltype(declval<backend_t>() / declval<T>())>
    Adaptive_Int& operator/=(T const&o) {
        value/=o;
        return *this;
    }
    template<typename T, typename = decltype(declval<backend_t>() % declval<T>())>
    Adaptive_Int operator%(T const&o) const {
        Adaptive_Int ret(value / o);
        return ret;
    }    template<typename T, typename = decltype(declval<backend_t>() % declval<T>())>
    Adaptive_Int& operator%=(T const&o) {
        value%=o;
        return *this;
    }

    template<size_t k = n+64>
    Adaptive_Int<k> operator<<(size_t const&o) const {
        assert(o <= 64);
        Adaptive_Int<k> ret(*this);
        ret.value<<= o;
        return ret;
    }
    Adaptive_Int operator>>(size_t const&o) const {
        Adaptive_Int ret(*this);
        ret.value>>= o;
        return ret;
    }
    Adaptive_Int& operator>>=(size_t const&o) {
        value>>= o;
        return *this;
    }
    Adaptive_Int& operator<<=(Unsafe_Wrapper<size_t const&> o) {
        value<<= o;
        return *this;
    }

    friend std::istream& operator>>(std::istream&in, Adaptive_Int &val){
        in >> val.value;
        return in;
    }
    friend std::ostream& operator<<(std::ostream&o, Adaptive_Int const&val){
        o << val.value;
        return o;
    }

    Adaptive_Int operator-() const {
        Adaptive_Int ret(-value);
        return ret;
    }
    bool operator!() const {
        return !value;
    }

    explicit operator double() const {
        return static_cast<double>(value);
    }
    explicit operator long double() const {
        return static_cast<long double>(value);
    }

    template<typename T>
    int comp(T const&o) const {
        return comp_impl(o, is_adaptive_int<T>{});
    }
    template<typename SFINAE=backend_t>
    enable_if_t<is_same_v<SFINAE, backend_t> && has_sign<backend_t, int>::value, int> sign() const{
        return value.sign();
    }
    template<typename SFINAE=backend_t>
    enable_if_t<is_same_v<SFINAE, backend_t> && !has_sign<backend_t, int>::value, int> sign() const{
        return comp(ZERO);
    }
    #define DECLARE_COMPARISON_OPERATOR(op)\
    template<typename T, typename = decltype(declval<Adaptive_Int>().comp(declval<T const&>()))>\
    bool operator op (T const&o)const{\
        return comp(o) op 0;\
    }
    DECLARE_COMPARISON_OPERATOR(<);
    DECLARE_COMPARISON_OPERATOR(<=);
    DECLARE_COMPARISON_OPERATOR(>);
    DECLARE_COMPARISON_OPERATOR(>=);
    DECLARE_COMPARISON_OPERATOR(==);
    DECLARE_COMPARISON_OPERATOR(!=);
    #undef DECLARE_COMPARISON_OPERATOR

private:
    template<typename T>
    int comp_impl(T const&o, std::true_type) const {
        return comp_impl_1(o, can_comp<backend_t, typename T::backend_t, int>{});
    }
    template<typename T>
    int comp_impl_1(T const&o, std::true_type) const {
        return value.comp(o.value);
    }
    template<typename T>
    int comp_impl_1(T const&o, std::false_type) const {
        return (value > o.value) - (value < o.value);
    }
    template<typename T>
    int comp_impl(T const&o, std::false_type) const {
        return comp_impl_2(o, can_comp<backend_t, T, int>{});
    }
    template<typename T>
    int comp_impl_2(T const&o, std::true_type) const {
        return value.comp(o);
    }
    template<typename T>
    int comp_impl_2(T const&o, std::false_type) const {
        return (value > o) - (value < o);
    }

    backend_t value;
};
template<size_t n>
class Adaptive_Int<n> Adaptive_Int<n>::ZERO;


#ifdef DACIN_HASH_HPP
template<size_t n>
struct Dacin_Hash<Adaptive_Int<n> > {
    size_t operator()(Adaptive_Int<n> const&val) const {
        static Dacin_Hash<typename Adaptive_Int<n>::backend_t> h;
        return h(val.get_cvalue());
    }
};
#endif // DACIN_HASH_HPP


} // namespace dacin::geom

#endif // ADAPTIVE_INT_HPP
