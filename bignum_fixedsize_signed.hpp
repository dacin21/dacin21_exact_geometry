// Released under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007, see the LICENSE file.
// Copyright (C) 2018-2019 Daniel Rutschmann aka. dacin21

#ifndef BIGNUM_FIXEDSIZE_SIGNED_HPP
#define BIGNUM_FIXEDSIZE_SIGNED_HPP

#include "geom_utility.hpp"

namespace dacin::geom{

template<size_t word_cnt, typename = enable_if_t<word_cnt != 0> >
class Bigint_Fixedsize_Signed{
public:
    std::array<uint64_t, word_cnt> data;

    template<typename T>
    using is_compatible_with = integral_constant<bool, is_integral_v<T> && (sizeof(T) >= 4)>;
    template<typename T>
    struct is_bigint : std::false_type{};
    template<size_t other_word_cnt, typename T>
    struct is_bigint<Bigint_Fixedsize_Signed<other_word_cnt, T> > : std::true_type{};
    template<typename T>
    using enable_if_by_construction_t = enable_if_t<is_constructible_v<Bigint_Fixedsize_Signed, T> && !is_bigint<decay_t<T> >::value>;

private:
    void negate(){
        size_t j = 0;
        while(j < word_cnt && !data[j]) ++j;
        if(j == word_cnt) return;
        --data[j];
        for(;j<word_cnt;++j){
            data[j]=~data[j];
        }
    }
    template<typename T>
    T convert_to_floating() const {
        static_assert(std::is_floating_point<T>::value);
        if(is_negative()){
            T ret = 0;
            for(auto it = data.rbegin(), it_end = data.rend();it!=it_end;++it){
                ret = (ret * (1ull<< 32) * (1ull<< 32)) - ~*it;
            }
            return ret - 1;
        } else {
            T ret = 0;
            for(auto it = data.rbegin(), it_end = data.rend();it!=it_end;++it){
                ret = (ret * (1ull<< 32) * (1ull<< 32)) + *it;
            }
            return ret;
        }
    }

    Bigint_Fixedsize_Signed& impl_mul_inplace(int64_t const&val, std::true_type){
        uint64_t res = val;
        if(val<0){
            negate();
            res = -res;
        }
        return operator*=(res);
    }
    Bigint_Fixedsize_Signed& impl_mul_inplace(uint64_t const&o, std::true_type){
        uint128_t carry = 0;
        for(auto &e:data){
            carry+= e*static_cast<uint128_t>(o);
            e = carry;
            carry>>=64;
        }
        return *this;
    }
    template<size_t other_word_cnt>
    Bigint_Fixedsize_Signed& impl_mul_inplace(Bigint_Fixedsize_Signed<other_word_cnt> const&o, std::false_type){
        static_assert(other_word_cnt <= word_cnt);
        mul(data, data, o.data);
        return *this;
    }
    template<typename T, typename = enable_if_by_construction_t<T> >
    Bigint_Fixedsize_Signed& impl_mul_inplace(T const& val, std::false_type){
        return operator*=(Bigint_Fixedsize_Signed(val));
    }
    template<typename T>
    Bigint_Fixedsize_Signed& impl_mul_inplace(T const&o){
        return impl_mul_inplace(o, is_one_of<decay_t<T>, int32_t, uint32_t, int64_t, uint64_t>{});
    }



    template<size_t n>
    static uint64_t get_pad(std::array<uint64_t, n> const&a){
        static_assert(n != 0);
        return -(a.back()>>63);
    }
    template<size_t n1, size_t n2>
    static int signed_comp(std::array<uint64_t, n1> const&a, std::array<uint64_t, n2> const&b){
        const uint64_t pad_a = get_pad(a), pad_b = get_pad(b);
        if(pad_a != pad_b){
            return pad_a - pad_b;
        }
        size_t j = max(n2, n1)-1;
        for(;j>=n2;--j){
            if(a[j] != pad_b) return pad_a ? -1:1;
        }
        for(;j>=n1;--j){
            if(b[j] != pad_a) return pad_a ? 1:-1;
        }
        for(;j+1;--j){
            if(a[j] != b[j]) return a[j]<b[j] ? -1 : 1;
        }
        return 0;
    }

    template<size_t n1, size_t n2, typename = enable_if_t<n2 <= n1> >
    static void add(std::array<uint64_t, n1> &a, std::array<uint64_t, n2> const&b){
        uint128_t carry = 0;
        for(size_t i=0;i<n2;++i){
            carry += a[i];
            carry += b[i];
            a[i] = carry;
            carry >>=64;
        }
        const uint64_t pad_b = get_pad(b);
        for(size_t i=n2;i<n1;++i){
            carry += a[i];
            carry += pad_b;
            a[i] = carry;
            carry>>=64;
        }
    }
    template<size_t n1, size_t n2, typename = enable_if_t<n2 <= n1> >
    static void sub(std::array<uint64_t, n1> &a, std::array<uint64_t, n2> const&b){
        uint128_t carry = 0;
        for(size_t i=0;i<n2;++i){
            carry += a[i];
            carry -= b[i];
            a[i] = carry;
            carry >>=64;
            if(carry>>63) carry|=~uint128_t{0}<<64;
        }
        for(size_t i=n2;i<n1;++i){
            carry+=a[i];
            a[i] = carry;
            carry>>=64;
            if(carry>>63) carry|=~uint128_t{0}<<64;
        }
    }

    template<size_t n>
    static void right_shift_small(std::array<uint64_t, n> &a, size_t const&c){
        if(!c) return;
        uint64_t carry = get_pad(a), tmp;
        for(size_t i=n-1;i+1;--i){
            carry<<=64-c;
            tmp = a[i];
            a[i]>>=c;
            a[i]|=carry;
            carry = tmp;
        }
    }
    template<size_t n>
    static void right_shift(std::array<uint64_t, n> &a, size_t c){
        right_shift_small(a, c%64);
        c/=64;
        copy(a.begin()+c, a.end(), a.begin());
        fill(a.end()-c, a.end(), 0);
    }
    template<size_t n>
    static void left_shift_small(std::array<uint64_t, n> &a, size_t const&c){
        if(!c) return;
        uint64_t carry = 0, tmp;
        for(size_t i=0;i<n;--i){
            carry>>=64-c;
            tmp = a[i];
            a[i]<<=c;
            a[i]|=carry;
            carry = tmp;
        }
    }
    template<size_t n>
    static void left_shift(std::array<uint64_t, n> &a, size_t c){
        left_shift_small(a, c%64);
        c/=64;
        const uint64_t pad = get_pad(a);
        copy(a.begin(), a.end()-c, a.begin()+c);
        fill(a.begin(), a.begin()+c, pad);
    }

    template<size_t n, size_t m>
    static void mul(std::array<uint64_t, n>&out, std::array<uint64_t, n> const&a, std::array<uint64_t, m> const&b){
        static std::array<uint64_t, n + max(n, m) + 1> tmp;
        std::fill(tmp.begin(), tmp.end(), 0);
        const uint64_t pad_b = get_pad(b);
        for(size_t i=0;i<n;++i){
            uint128_t carry = 0;
            for(size_t j=0;j<m;++j){
                carry+=a[i] * static_cast<uint128_t>(b[j]);
                carry+=tmp[i+j];
                tmp[i+j] = carry;
                carry>>=64;
            }
            if(m<n && pad_b){
                for(size_t j=m;j<n;++j){
                    carry+=a[i] * static_cast<uint128_t>(pad_b);
                    carry+=tmp[i+j];
                    tmp[i+j] = carry;
                    carry>>=64;
                }
                tmp[i+n] = carry;
            } else {
                tmp[i+m] = carry;
            }
        }
        std::copy(tmp.begin(), tmp.begin()+n, out.begin());
    }
    static uint64_t divmod(Bigint_Fixedsize_Signed &a, uint64_t const&d){
        bool nega = false;
        if(a.is_negative()){
            a.negate();
            nega = true;
        }
        assert(d != 0);
        uint128_t carry = 0;
        for(size_t i = word_cnt-1;i+1;--i){
            carry<<=64;
            carry += a.data[i];
            a.data[i] = carry / d;
            carry = carry % d;
        }
        if(nega){
            a.negate();
            if(carry) carry = d-carry;
        }
        return carry;
    }
    static bool print_destructive(std::ostream&o, Bigint_Fixedsize_Signed &v){
        if(!v) return true;
        uint32_t last = divmod(v, 1000000000);
        bool is_first = print_destructive(o, v);
        if(!is_first){
            o << std::setw(9) << std::setfill('0');
        }
        o << last;
        return false;
    }

public:
    bool is_negative()const{
        return data.back()>>63;
    }
    int sign() const{
        if(is_negative()) return -1;
        for(size_t i=0;i<word_cnt;++i) if(data[i]) return 1;
        return 0;
    }
    bool is_positive()const{
        return sign() == 1;
    }
    Bigint_Fixedsize_Signed():data{}{}
    explicit Bigint_Fixedsize_Signed(uint32_t const&val):data{val}{}
    explicit Bigint_Fixedsize_Signed(int32_t const&val):data{static_cast<uint32_t>(val)}{
        if(val < 0) {
            data[0] |= ~uint64_t{}<<32;
            std::fill(data.begin()+1, data.end(), ~uint64_t{});
        }
    }
    explicit Bigint_Fixedsize_Signed(uint64_t const&val):data{val}{}
    explicit Bigint_Fixedsize_Signed(int64_t val):data{static_cast<uint64_t>(val)}{
        if(val<0) std::fill(data.begin()+1, data.end(), ~uint64_t{});
    }
    explicit Bigint_Fixedsize_Signed(__int128 val):data{static_cast<uint64_t>(val), static_cast<uint64_t>(val>>64)}{
        static_assert(word_cnt >= 2);
        if(val<0) std::fill(data.begin()+2, data.end(), ~uint64_t{});
    }
    template<size_t other_word_cnt, typename = enable_if_t<other_word_cnt <= word_cnt> >
    explicit Bigint_Fixedsize_Signed(Bigint_Fixedsize_Signed<other_word_cnt> const&val):data{}{
        std::copy(val.data.begin(), val.data.end(), data.begin());
        if(val.is_negative()) std::fill(data.begin()+other_word_cnt, data.end(), ~uint64_t{});
    }

    template<size_t other_word_cnt, typename = enable_if_t<other_word_cnt <= word_cnt> >
    Bigint_Fixedsize_Signed& operator+=(Bigint_Fixedsize_Signed<other_word_cnt> const&o){
        add(data, o.data);
        return *this;
    }
    template<typename T, typename = enable_if_by_construction_t<T> >
    Bigint_Fixedsize_Signed& operator+=(T const&o){
        return operator+=(Bigint_Fixedsize_Signed(o));
    }
    template<size_t other_word_cnt, typename = enable_if_t<other_word_cnt <= word_cnt> >
    Bigint_Fixedsize_Signed<word_cnt> operator+(Bigint_Fixedsize_Signed<other_word_cnt> const&o)const{
        Bigint_Fixedsize_Signed<word_cnt> ret(*this);
        ret+=o;
        return ret;
    }
    template<size_t other_word_cnt, typename = enable_if_t<word_cnt < other_word_cnt> >
    Bigint_Fixedsize_Signed<other_word_cnt> operator+(Bigint_Fixedsize_Signed<other_word_cnt> const&o)const{
        Bigint_Fixedsize_Signed<other_word_cnt> ret(o);
        ret+=*this;
        return ret;
    }
    template<typename T, typename = enable_if_by_construction_t<T> >
    Bigint_Fixedsize_Signed& operator+(T const&o){
        return operator+(Bigint_Fixedsize_Signed(o));
    }
    template<size_t other_word_cnt, typename = enable_if_t<other_word_cnt <= word_cnt> >
    Bigint_Fixedsize_Signed& operator-=(Bigint_Fixedsize_Signed<other_word_cnt> const&o){
        sub(data, o.data);
        return *this;
    }
    template<typename T, typename = enable_if_by_construction_t<T> >
    Bigint_Fixedsize_Signed& operator-=(T const&o){
        return operator-=(Bigint_Fixedsize_Signed(o));
    }
    template<size_t other_word_cnt, typename = enable_if_t<other_word_cnt <= word_cnt> >
    Bigint_Fixedsize_Signed<word_cnt> operator-(Bigint_Fixedsize_Signed<other_word_cnt> const&o)const{
        Bigint_Fixedsize_Signed<word_cnt> ret(*this);
        ret-=o;
        return ret;
    }
    template<size_t other_word_cnt, typename = enable_if_t<word_cnt < other_word_cnt  > >
    Bigint_Fixedsize_Signed<other_word_cnt> operator-(Bigint_Fixedsize_Signed<other_word_cnt> const&o)const{
        Bigint_Fixedsize_Signed<other_word_cnt> ret(o);
        ret-=*this;
        ret.negate();
        return ret;
    }
    template<typename T, typename = enable_if_by_construction_t<T> >
    Bigint_Fixedsize_Signed& operator-(T const&o){
        return operator-(Bigint_Fixedsize_Signed(o));
    }

    template<typename T>
    Bigint_Fixedsize_Signed& operator*=(T const&o){
        return impl_mul_inplace(o);
    }

    template<size_t other_word_cnt, typename = enable_if_t<other_word_cnt <= word_cnt> >
    Bigint_Fixedsize_Signed operator*(Bigint_Fixedsize_Signed<other_word_cnt> const&o)const{
        Bigint_Fixedsize_Signed ret;
        mul(ret.data, data, o.data);
        return ret;
    }
    template<size_t other_word_cnt, typename = enable_if_t<word_cnt < other_word_cnt> >
    Bigint_Fixedsize_Signed<other_word_cnt>& operator*(Bigint_Fixedsize_Signed<other_word_cnt> const&o){
        Bigint_Fixedsize_Signed<other_word_cnt> ret;
        mul(ret.data, o.data, data);
        return *this;
    }

    Bigint_Fixedsize_Signed operator*(uint64_t const&val)const{
        Bigint_Fixedsize_Signed ret(*this);
        ret*=val;
        return ret;
    }
    Bigint_Fixedsize_Signed operator*(int64_t const&val)const{
        Bigint_Fixedsize_Signed ret(*this);
        ret*=val;
        return ret;
    }
    template<typename T, typename = enable_if_by_construction_t<T>, typename = enable_if_t<!is_one_of<decay_t<T>, int32_t, uint32_t, int64_t, uint64_t>::value> >
    Bigint_Fixedsize_Signed operator*(T const& val)const{
        return operator*(Bigint_Fixedsize_Signed(val));
    }
    Bigint_Fixedsize_Signed& operator/=(uint32_t const&d){
        divmod(*this, d);
        return *this;
    }
    Bigint_Fixedsize_Signed operator/(uint32_t const&d)const{
        Bigint_Fixedsize_Signed ret(*this);
        ret/=d;
        return ret;
    }
    Bigint_Fixedsize_Signed& operator/=(int32_t const&d){
        uint32_t res = d;
        if(d<0){
            negate();
            res = -res;
        }
        return operator/=(res);
    }
    Bigint_Fixedsize_Signed operator/(int32_t const&d)const{
        Bigint_Fixedsize_Signed ret(*this);
        ret/=d;
        return ret;
    }
    uint32_t operator%(uint32_t const&d){
        Bigint_Fixedsize_Signed tmp(*this);
        return divmod(tmp, d);
    }
    int32_t operator%(int32_t const&d){
        bool nega = is_negative();
        uint32_t res = d;
        if(d < 0){
            nega = !nega;
            res=-res;
        }
        uint32_t ret = operator%(res);
        if(nega && ret) ret-=res;
        return *reinterpret_cast<int32_t*>(&ret);
    }

    Bigint_Fixedsize_Signed<word_cnt>& operator<<=(size_t const&s){
        left_shift(data, s);
        return *this;
    }
    Bigint_Fixedsize_Signed<word_cnt> operator<<(size_t const&s)const{
        Bigint_Fixedsize_Signed ret(*this);
        ret<<=s;
        return ret;
    }
    Bigint_Fixedsize_Signed<word_cnt>& operator>>=(size_t const&s){
        right_shift(data, s);
        return *this;
    }
    Bigint_Fixedsize_Signed<word_cnt> operator>>(size_t const&s) const {
        Bigint_Fixedsize_Signed ret(*this);
        ret>>=s;
        return ret;
    }

    Bigint_Fixedsize_Signed<word_cnt> operator-() const {
        Bigint_Fixedsize_Signed ret(*this);
        ret.negate();
        return ret;
    }
    bool operator!() const {
        return sign() == 0;
    }
    explicit operator double() const {
        return convert_to_floating<double>();
    }
    explicit operator long double() const {
        return convert_to_floating<long double>();
    }


    static void print_bin(std::ostream&o, Bigint_Fixedsize_Signed val){
        for(auto it = val.data.rbegin(); it != val.data.rend();++it){
            o << std::bitset<64>(*it);
        }
    }
    friend std::ostream& operator<<(std::ostream&o, Bigint_Fixedsize_Signed val){
        const int s = val.sign();
        if(s == 0) return o << 0;
        if(s < 0){
            uint32_t tail = divmod(val, 10);
            if(!val){
                return o << "-" << (10-tail)%10;
            }
            val.negate();
            o << '-';
            print_destructive(o, val);
            o << (10-tail)%10;
        } else {
            print_destructive(o, val);
        }
        return o;
    }
    template<size_t other_word_cnt>
    int comp(Bigint_Fixedsize_Signed<other_word_cnt> const&o)const{
        return signed_comp(data, o.data);
    }
    template<typename T, typename = enable_if_by_construction_t<T>>
    int comp(T const&o)const{
        return comp(Bigint_Fixedsize_Signed(o));
    }
    #define DECLARE_COMPARISON_OPERATOR(op)\
    template<size_t other_word_cnt>\
    bool operator op (Bigint_Fixedsize_Signed const&o) const {\
        return comp(o) op 0;\
    }\
    template<typename T>\
    bool operator op (T const&o) const {\
        return comp(o) op 0;\
    }\
    template<typename T>\
    bool friend operator op (T const&o, Bigint_Fixedsize_Signed const&me){\
        return 0 op me.comp(o);\
    }
    DECLARE_COMPARISON_OPERATOR(<);
    DECLARE_COMPARISON_OPERATOR(<=);
    DECLARE_COMPARISON_OPERATOR(>);
    DECLARE_COMPARISON_OPERATOR(>=);
    DECLARE_COMPARISON_OPERATOR(==);
    DECLARE_COMPARISON_OPERATOR(!=);
    #undef DECLARE_COMPARISON_OPERATOR
};

#ifdef DACIN_HASH_HPP
template<size_t word_cnt>
struct Dacin_Hash<Bigint_Fixedsize_Signed<word_cnt>>{
    size_t old_hash(std::vector<uint32_t> const&data) const {
        size_t ret = SALT;
        const uint32_t pad = data.back();
        int i = data.size()-1;
        if(pad == 0 || pad == ~0u){
            while(i>0 && data[i] == pad) --i;
        }
        if(pad != 0 && i+1 < (int)data.size()) ++i;
        for(size_t j=0, lim=i+1;j<lim;++j){
            ret = splitmix64(data[i] + ret);
        }
        return ret;
    }
    size_t new_hash(Bigint_Fixedsize_Signed<word_cnt> const&val) const {
        size_t ret = SALT;
        const uint64_t pad = val.data.back();
        int i = val.data.size()-1;
        if(pad == 0 || pad == ~0ull){
            while(i>0 && val.data[i] == pad) --i;
        }
        if(pad != 0 && i+1 < (int)val.data.size()) ++i;
        for(size_t j=0, lim=i+1;j<lim;++j){
            ret = splitmix64(val.data[i] + ret);
        }
        return ret;
    }
    size_t operator() (Bigint_Fixedsize_Signed<word_cnt> const&val) const {
        std::vector<uint32_t> data(2*val.data.size());
        for(int i=0;i<(int)val.data.size();++i){
            data[2*i] = val.data[i];
            data[2*i+1] = val.data[i]>>32;
        }
        return old_hash(data);
    }
};
#endif // DACIN_HASH_HPP

} // namespace dacin::geom

#endif // BIGNUM_FIXEDSIZE_SIGNED_HPP
