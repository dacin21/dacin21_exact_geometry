// Released under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007, see the LICENSE file.
// Copyright (C) 2018-2019 Daniel Rutschmann aka. dacin21

#ifdef RELEASE
#undef _GLIBCXX_DEBUG
#endif // RELEASE

#include "svg_plot.hpp"
#include "geom_utility.hpp"
#include "dacin_hash.hpp"
#include "bignum_fixedsize_signed.hpp"
#include "adaptive_int.hpp"
#include "geom_2d.hpp"
#include "delaunay.hpp"
