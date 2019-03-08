
#ifndef DACIN_HASH_HPP
#define DACIN_HASH_HPP

#include "geom_utility.hpp"

namespace dacin::geom{

const uint64_t SALT = std::chrono::steady_clock::now().time_since_epoch().count();
uint64_t splitmix64(uint64_t x) {
    // http://xorshift.di.unimi.it/splitmix64.c
    x+= 0x9e3779b97f4a7c15;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
    x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
    return x ^ (x >> 31);
}

template<typename T, typename SFINAE = void>
struct Dacin_Hash{};

template<typename T>
struct Dacin_Hash<T, enable_if_t<is_same_v<decltype(declval<std::hash<T> >()(declval<T>())), size_t> > >{
    size_t operator()(T const&t) const {
        static std::hash<T> hasher;
        return splitmix64(hasher(t)+SALT);
    }
};

template<typename S, typename T>
struct Dacin_Hash<std::pair<S, T> >{
    size_t operator()(std::pair<S, T> const&x) const {
        static Dacin_Hash<decay_t<S> > h1;
        static Dacin_Hash<decay_t<T> > h2;
        return splitmix64((h1(x.first)^SALT) + h2(x.second));
    }
};
template<typename T>
struct Dacin_Hash<std::vector<T> >{
    size_t operator()(std::vector<T> const&o) const {
        static Dacin_Hash<decay_t<T> > h;
        size_t ret = 0;
        for(auto &e:o){
            ret = splitmix64(ret+h(e));
        }
        return ret;
    }
};

} // namespace dacin::geom

#endif // DACIN_HASH_HPP
