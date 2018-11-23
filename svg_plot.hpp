// Released under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007, see the LICENSE file.
// Copyright (C) 2018 Daniel Rutschmann aka. dacin21

#ifndef SVG_PLOT_HPP
#define SVG_PLOT_HPP

#include <bits/stdc++.h>

namespace dacin::geom{

class SVG{
public:
    enum class Coords{
        X_RIGHT_Y_DOWN,
        X_RIGHT_Y_UP,
        X_DOWN_Y_RIGHT,
    };
    SVG(std::string const&path, const double x1, const double x2, const double y1, const double y2, Coords coords_ = Coords::X_RIGHT_Y_UP) : outfile(path), coords(coords_), range_x(x1, x2), range_y(y1, y2){
        outfile << std::fixed << std::setprecision(4);
        outfile << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n";
        outfile << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" ";
        outfile << "viewBox=\"";
        auto const A = resolve_coords(x1, y1), B = resolve_coords(x2, y2);
        outfile << std::min(A.first, B.first) << " " << std::min(A.second, B.second) << " " << std::abs(A.first- B.first) << " " << std::abs(A.second- B.second);
        outfile << "\">\n";
    }
    SVG(SVG const&) = delete;
    SVG(SVG const&&) = delete;
    SVG& operator=(SVG const&) = delete;
    SVG& operator=(SVG const&&) = delete;
    ~SVG(){
        outfile << "</svg>\n";
        outfile.close();
    }

    void set_fill_color(int r, int g, int b){
        set_color_string(color_fill, r, g, b);
    }
    void set_stroke_color(int r, int g, int b){
        set_color_string(color_stroke, r, g, b);
    }
    void set_stroke_width(int w){
        stroke_width = w;
    }

    void draw_point(double x, double y, double scale = 0.005){
        double r = (range_x.second-range_x.first + range_y.second-range_y.first) / 2 * scale;
        draw_circle(x, y, r, false, true);
        draw_circle(x, y, r/2, true, false);
    }
    void draw_segment(double x1, double y1, double x2, double y2, bool stroke = true, bool fill = false){
        std::tie(x1,y1) = resolve_coords(x1, y1);
        std::tie(x2,y2) = resolve_coords(x2, y2);
        print_with_arributes("line", stroke, fill, "x1", x1, "y1", y1, "x2", x2, "y2", y2);
    }
    void draw_rectangle(double x, double y, double w, double h, bool stroke = true, bool fill = false){
        assert(w >= 0 && h >= 0);
        double x1, y1, x2, y2;
        std::tie(x1,y1) = resolve_coords(x, y);
        std::tie(x2,y2) = resolve_coords(x+w, y+h);
        print_with_arributes("rect", stroke, fill, "x", std::min(x1, x2), "y", std::min(y1, y2), "width", std::abs(x2-x1), "height", std::abs(y2-y1));
        // needs to deal with coordinate systems
    }
    void draw_circle(double cx, double cy, double r, bool stroke = true, bool fill = false){
        assert(r >= 0);
        std::tie(cx, cy) = resolve_coords(cx, cy);
        print_with_arributes("circle", stroke, fill, "cx", cx, "cy", cy, "r", r);
    }
    void draw_ellipse(double cx, double cy, double rx, double ry, bool stroke = true, bool fill = false){
        assert(rx >= 0 && ry >= 0);
        std::tie(cx, cy) = resolve_coords(cx, cy);
        print_with_arributes("ellipse", stroke, fill, "cx", cx, "cy", cy, "rx", rx, "ry", ry);
    }
    // needs support for std::vector in attributes
    void draw_polygon(std::vector<std::pair<int, int> > v, bool stroke = true, bool fill = false){
        for(auto &e:v) e = resolve_coords(e.first, e.second);
        print_with_arributes("polygon", stroke, fill, "points", v, "stroke-linejoin", "round");
    }
    void draw_polyline(std::vector<std::pair<int, int> > v, bool stroke = true, bool fill = false){
        for(auto &e:v) e = resolve_coords(e.first, e.second);
        print_with_arributes("polyline", stroke, fill, "points", v);
    }

private:
    void print_stroke_width(){
        outfile << " stroke-width=\"" << stroke_width << "px\"";
    }
    void print_stroke(){
        outfile << " stroke=\"" << color_stroke << '"';
    }
    void print_fill(){
        outfile << " fill=\"" << color_fill << '"';
    }
    void print_color(bool stroke, bool fill){
        if(fill) print_fill();
        else outfile << " fill=\"none\"";
        if(stroke){
            print_stroke();
            print_stroke_width();
        } else {
            outfile << " stroke=\"none\"";
        }
    }
    template<typename T, typename = typename std::enable_if<!std::is_class<T>::value>::type>
    void print_value(T const&val){
        outfile << val;
    }
    template<typename S, typename T>
    void print_value(std::pair<S, T> const&val){
        print_value(val.first);
        outfile << ',';
        print_value(val.second);
    }
    template<typename T, typename = decltype(std::begin(std::declval<T>())), typename = decltype(std::end(std::declval<T>()))>
    void print_value(T const&v){
        bool first = true;
        for(auto &e:v){
            if(!first) outfile << ' ';
            print_value(e);
            first = false;
        }
    }
    void print_attributes(){}
    template<typename T, typename ... Args>
    void print_attributes(std::string const&name, T const&val, Args&& ... args){
        outfile << ' ' << name << '=' << '"';
        print_value(val);
        outfile << '"';
        print_attributes(std::forward<Args>(args) ...);
    }
    template<typename ... Args>
    void print_with_arributes(std::string const&name, bool stroke, bool fill, Args&& ... args){
        outfile << "<" << name;
        print_color(stroke, fill);
        print_attributes(std::forward<Args>(args) ...);
        outfile << "/>\n";
    }

    std::pair<double, double> resolve_coords(std::pair<double, double> const&p){
        switch(coords){
            case Coords::X_RIGHT_Y_DOWN:
                return p;
            case Coords::X_RIGHT_Y_UP:
                return std::make_pair(p.first, -p.second);
            case Coords::X_DOWN_Y_RIGHT:
                return std::make_pair(p.second, p.first);
            default:
                assert(0);
                return p;
        }
    }
    std::pair<double, double> resolve_coords(double const&x, double const&y){
        return resolve_coords(std::make_pair(x, y));
    }

    void set_color_string(std::string&s, int r, int g, int b){
        std::stringstream ss;
        ss << '#' << std::setfill('0') << std::hex << std::setw(2) << r << std::setw(2) << g << std::setw(2) << b;
        s = ss.str();
        assert(s.size() == 7);
    }

    std::ofstream outfile;
    const Coords coords;
    const std::pair<double, double> range_x, range_y;
    std::string color_stroke{"#222222"}, color_fill{"#444444"};
    int stroke_width = 2;
};

} // namespace dacin::geom

#endif // SVG_PLOT_HPP
