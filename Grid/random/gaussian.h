// -*- C++ -*-
//===--------------------------- random -----------------------------------===//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

// Peter Boyle: Taken from libc++ in Clang/LLVM.
// Reason is that libstdc++ and clang differ in their return order in the normal_distribution / box mueller type step.
// standardise on one and call it "gaussian_distribution".

#pragma once

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <type_traits>
#include <initializer_list>
#include <limits>
#include <algorithm>
#include <numeric>
#include <vector>
#include <string>
#include <istream>
#include <ostream>
#include <random>

// normal_distribution -> gaussian distribution
namespace Grid {

template<class _RealType = double>
class  gaussian_distribution
{
public:
    // types
    typedef _RealType result_type;

    class param_type
    {
        result_type __mean_;
        result_type __stddev_;
    public:
        typedef gaussian_distribution distribution_type;

        strong_inline
        explicit param_type(result_type __mean = 0, result_type __stddev = 1)
            : __mean_(__mean), __stddev_(__stddev) {}

        strong_inline
        result_type mean() const {return __mean_;}
        strong_inline
        result_type stddev() const {return __stddev_;}

        friend strong_inline
            bool operator==(const param_type& __x, const param_type& __y)
            {return __x.__mean_ == __y.__mean_ && __x.__stddev_ == __y.__stddev_;}
        friend strong_inline
            bool operator!=(const param_type& __x, const param_type& __y)
            {return !(__x == __y);}
    };

private:
    param_type __p_;
    result_type _V_;
    bool _V_hot_;

public:
    // constructors and reset functions
    strong_inline
    explicit gaussian_distribution(result_type __mean = 0, result_type __stddev = 1)
        : __p_(param_type(__mean, __stddev)), _V_hot_(false) {}
    strong_inline
    explicit gaussian_distribution(const param_type& __p)
        : __p_(__p), _V_hot_(false) {}
    strong_inline
    void reset() {_V_hot_ = false;}

    // generating functions
    template<class _URNG>
        strong_inline
        result_type operator()(_URNG& __g)
        {return (*this)(__g, __p_);}
    template<class _URNG> result_type operator()(_URNG& __g, const param_type& __p);

    // property functions
    strong_inline
    result_type mean() const {return __p_.mean();}
    strong_inline
    result_type stddev() const {return __p_.stddev();}

    strong_inline
    param_type param() const {return __p_;}
    strong_inline
    void param(const param_type& __p) {__p_ = __p;}

    strong_inline
    result_type min() const {return -std::numeric_limits<result_type>::infinity();}
    strong_inline
    result_type max() const {return std::numeric_limits<result_type>::infinity();}

    friend strong_inline
        bool operator==(const gaussian_distribution& __x,
                        const gaussian_distribution& __y)
        {return __x.__p_ == __y.__p_ && __x._V_hot_ == __y._V_hot_ &&
                (!__x._V_hot_ || __x._V_ == __y._V_);}
    friend strong_inline
        bool operator!=(const gaussian_distribution& __x,
                        const gaussian_distribution& __y)
        {return !(__x == __y);}

    template <class _CharT, class _Traits, class _RT>
    friend
    std::basic_ostream<_CharT, _Traits>&
    operator<<(std::basic_ostream<_CharT, _Traits>& __os,
               const gaussian_distribution<_RT>& __x);

    template <class _CharT, class _Traits, class _RT>
    friend
    std::basic_istream<_CharT, _Traits>&
    operator>>(std::basic_istream<_CharT, _Traits>& __is,
               gaussian_distribution<_RT>& __x);
};

template <class _RealType>
template<class _URNG>
_RealType
gaussian_distribution<_RealType>::operator()(_URNG& __g, const param_type& __p)
{
    result_type _Up;
    if (_V_hot_)
    {
        _V_hot_ = false;
        _Up = _V_;
    }
    else
    {
        std::uniform_real_distribution<result_type> _Uni(-1, 1);
        result_type __u;
        result_type __v;
        result_type __s;
        do
        {
            __u = _Uni(__g);
            __v = _Uni(__g);
            __s = __u * __u + __v * __v;
        } while (__s > 1 || __s == 0);
        result_type _Fp = std::sqrt(-2 * std::log(__s) / __s);
        _V_ = __v * _Fp;
        _V_hot_ = true;
        _Up = __u * _Fp;
    }
    return _Up * __p.stddev() + __p.mean();
}

template <class _CharT, class _Traits, class _RT>
std::basic_ostream<_CharT, _Traits>&
operator<<(std::basic_ostream<_CharT, _Traits>& __os,
           const gaussian_distribution<_RT>& __x)
{
    auto __save_flags = __os.flags();
    __os.flags(std::ios_base::dec | std::ios_base::left | std::ios_base::fixed |
               std::ios_base::scientific);
    _CharT __sp = __os.widen(' ');
    __os.fill(__sp);
    __os << __x.mean() << __sp << __x.stddev() << __sp << __x._V_hot_;
    if (__x._V_hot_)
        __os << __sp << __x._V_;
    __os.flags(__save_flags);
    return __os;
}

template <class _CharT, class _Traits, class _RT>
std::basic_istream<_CharT, _Traits>&
operator>>(std::basic_istream<_CharT, _Traits>& __is,
           gaussian_distribution<_RT>& __x)
{
    typedef gaussian_distribution<_RT> _Eng;
    typedef typename _Eng::result_type result_type;
    typedef typename _Eng::param_type param_type;
    auto __save_flags = __is.flags();
    __is.flags(std::ios_base::dec | std::ios_base::skipws);
    result_type __mean;
    result_type __stddev;
    result_type _Vp = 0;
    bool _V_hot = false;
    __is >> __mean >> __stddev >> _V_hot;
    if (_V_hot)
        __is >> _Vp;
    if (!__is.fail())
    {
        __x.param(param_type(__mean, __stddev));
        __x._V_hot_ = _V_hot;
        __x._V_ = _Vp;
    }
    __is.flags(__save_flags);
    return __is;
}
}
