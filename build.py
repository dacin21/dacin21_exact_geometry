#!/usr/bin/python3

##
#  Movified version of Buildscript from jgen
#  https://github.com/ifsmirnov/jngen
#  used to move stuff into a single file
##

import re, os

HEADER_REGEX = re.compile('#include "(.*)"')

# This list may contain not all headers directly, but each jngen header
# must be among the dependencies of some file from here.
LIBRARY_HEADERS = [
    "examples/opencup_2019_spb_G.cpp",
]


def posix_path_to_native(posix_path):
    return os.path.join(*posix_path.split('/'))


def extract_header(line):
    res = HEADER_REGEX.match(line)
    if res:
        return res.groups()[0]


def extract_direct_deps(posix_filename):
    dir = os.path.dirname(posix_filename) # check explicitly on win

    res = set()
    with open(posix_path_to_native(posix_filename)) as fin:
        for line in fin.readlines():
            t = extract_header(line)
            if t and not t.endswith("_inl.h"):
                res.add(dir + '/' + t if dir else t)

    return res


deps = {}


def extract_deps(posix_filename):
    posix_filename = os.path.normpath(posix_filename)
    if posix_filename in deps:
        return deps[posix_filename]
    deps[posix_filename] = set((posix_filename,))

    for dep in extract_direct_deps(posix_filename):
        deps[posix_filename].update(extract_deps(dep))
    return deps[posix_filename]

def compress_whitespace(line):
    old_indent = 4
    new_indent = 1 # lower -> smaller file size
    line_rest = line.lstrip()
    if not line_rest: # keep empty lines
        return line
    # line_rest = str(line).lstrip()
    indent = (len(line) - len(line_rest) + old_indent-1)//old_indent
    return ' '*(new_indent*indent) + line_rest

def write_file(filename, stream):
    dir = os.path.dirname(filename) # check explicitly on win
    with open(posix_path_to_native(filename)) as fin:
        for line in fin.readlines():
            include_or_not = HEADER_REGEX.match(line)
            if include_or_not:
                if include_or_not.groups()[0].endswith("_inl.h"):
                    t = include_or_not.groups()[0]
                    write_file(dir + '/' + t if dir else t, stream)
            elif '#pragma once' not in line:
                stream.write(compress_whitespace(line))


headers = set()
for h in LIBRARY_HEADERS:
    headers.update(extract_deps(h))
headers = ['header.h'] + sorted(headers)
deps['header.h'] = set(('header.h',))

print("header: {}".format(headers))

headers_in_order = []
while headers:
    for h in headers:
        if len(deps[h]) == 1:
            headers_in_order.append(h)
            for other in deps:
                deps[other].discard(h)
            del deps[h]
            headers.remove(h)
            break

with open("exact_geometry_merged.cpp", "w") as fout:
    for filename in headers_in_order:
        write_file(filename, fout)

