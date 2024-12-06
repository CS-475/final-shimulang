#ifndef PA4_SHIMULANG_MAIN_MYGCANVAS_HPP
#define PA4_SHIMULANG_MAIN_MYGCANVAS_HPP

#include "include/GCanvas.h"
#include "include/GBitmap.h"
#include "include/GPathBuilder.h"
#include "include/GPath.h"
#include "MyGShader.hpp"
#include "blendModes.h"

#include <stack>
#include <map>
#include <iostream>
#include <algorithm>

[[maybe_unused]] static GPoint BezierQuad(const GPoint& P0, const GPoint& P1, const GPoint& P2, float t) {
    float x = (1 - t) * (1 - t) * P0.x + 2 * (1 - t) * t * P1.x + t * t * P2.x;
    float y = (1 - t) * (1 - t) * P0.y + 2 * (1 - t) * t * P1.y + t * t * P2.y;
    return GPoint{x, y};
}

[[maybe_unused]] static GPoint BezierCubic(const GPoint& P0, const GPoint& P1, const GPoint& P2,const GPoint& P3, float t) {
    float x = (1 - t) * (1 - t) * (1 - t) * P0.x + 3 * (1 - t) * (1 - t) * t * P1.x + 3 * (1 - t) * t * t * P2.x + t * t * t * P3.x;
    float y = (1 - t) * (1 - t) * (1 - t) * P0.y + 3 * (1 - t) * (1 - t) * t * P1.y + 3 * (1 - t) * t * t * P2.y + t * t * t * P3.y;
    return GPoint{x, y};
}

static inline GIRect intersect(const GRect& src_rect, const GIRect& dest_rect){
    GIRect src_rect_int = src_rect.round();
    int top = std::max(src_rect_int.top, dest_rect.top);
    int bottom = std::min(src_rect_int.bottom, dest_rect.bottom);
    int left = std::max(src_rect_int.left, dest_rect.left);
    int right = std::min(src_rect_int.right, dest_rect.right);

    GIRect result = GIRect::LTRB(left, top, right, bottom);

    return result;
}

void drawRow(int x0, int y0, int width, GPixel src, GPaint paint, GBitmap bm){
    // todo make a cpp file and declare drawRow here
    if (width < 1) return;

    unsigned sa = GPixel_GetA(src);
    GBlendMode mode = paint.getBlendMode();
    mode = blendTable(mode, sa);

    GPixel* p0 = bm.getAddr(x0, y0);
    if(mode == GBlendMode::kSrc){
        std::fill_n(p0, width, src);
        return;
    }
    if(mode == GBlendMode::kClear){
        std::fill_n(p0, width, GPixel_PackARGB(0, 0, 0, 0));
        return;
    }
    if(mode == GBlendMode::kDst){
        return;
    }
    for (int x = x0; x < x0 + width; x++){
        GPixel* p = bm.getAddr(x, y0);
        *p = blend(src, *p, mode);
    }
}

class MyGCanvas : public GCanvas{
public:
    MyGCanvas(const GBitmap& bitmap);

public:
    void save() override;
    void restore() override;
    void concat(const GMatrix& matrix) override;
    void clear(const GColor& color) override;

    void drawRect(const GRect& rect, const GPaint& paint) override;
    void drawConvexPolygon(const GPoint[], int count, const GPaint& paint) override;
    void drawPath(const GPath& path, const GPaint& paint) override;
    void drawMesh(const GPoint verts[], const GColor colors[], const GPoint texs[],
                  int count, const int indices[], const GPaint&) override;
    void drawQuad(const GPoint verts[4], const GColor colors[4], const GPoint texs[4],
                  int level, const GPaint&) override;

private:
    bool pointInBound(const GPoint& p1);
    bool edgeInBound(const GPoint& p1, const GPoint& p2, std::vector<GPoint>& crossPoint);
    GPoint letPointInBound(const GPoint& p);
    nonstd::optional<GPoint> getIntersectionPoint(const GPoint& p0, const GPoint& p1, const GPoint& q0, const GPoint& q1);
    std::vector<std::pair<GPoint, GPoint>> createPolygon(const std::vector<std::pair<GPoint, GPoint>>& raw_polygon);

private:
    std::vector<int> getIntersections(const std::vector<std::pair<GPoint, GPoint>>& edges, int y);
    std::vector<GPoint> fillPolygon(const std::vector<std::pair<GPoint, GPoint>>& raw_polygon);

    std::vector<GPoint> bresenhamLine(const GPoint& p1, const GPoint& p2);
    std::map<float, std::pair<float, float>> getLinePolygon(const std::vector<std::pair<GPoint, GPoint>>& raw_polygon);

    float cross_product(const GPoint& p1, const GPoint& p2, const GPoint& p3) {
        return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
    }

    bool is_between(const GPoint& p, const GPoint& a, const GPoint& b) {
        return std::min(a.x, b.x) <= p.x && p.x <= std::max(a.x, b.x) &&
               std::min(a.y, b.y) <= p.y && p.y <= std::max(a.y, b.y);
    }

    bool segment_intersection(const GPoint& p1, const GPoint& p2, const GPoint& q1, const GPoint& q2, GPoint& intersection) {
        float d1 = cross_product(q1, q2, p1);
        float d2 = cross_product(q1, q2, p2);
        float d3 = cross_product(p1, p2, q1);
        float d4 = cross_product(p1, p2, q2);

        if (d1 * d2 < 0 && d3 * d4 < 0) {
            float t = d1 / (d1 - d2);
            intersection.x = p1.x + t * (p2.x - p1.x);
            intersection.y = p1.y + t * (p2.y - p1.y);
            return true;
        }

        if (d1 == 0 && is_between(p1, q1, q2)) { intersection = p1; return true; }
        if (d2 == 0 && is_between(p2, q1, q2)) { intersection = p2; return true; }
        if (d3 == 0 && is_between(q1, p1, p2)) { intersection = q1; return true; }
        if (d4 == 0 && is_between(q2, p1, p2)) { intersection = q2; return true; }

        return false;
    }

    bool is_point_in_polygon(const GPoint& point, const std::vector<std::pair<GPoint, GPoint>>& polygon_edges) {
        int count = 0;
        for (const auto& edge : polygon_edges) {
            GPoint a = edge.first;
            GPoint b = edge.second;
            if ((a.y > point.y) != (b.y > point.y) &&
                (point.x < (b.x - a.x) * (point.y - a.y) / (b.y - a.y) + a.x)) {
                count++;
            }
        }
        return (count % 2) == 1;
    }

    bool is_polygon_contained(const std::vector<std::pair<GPoint, GPoint>>& poly1,
                              const std::vector<std::pair<GPoint, GPoint>>& poly2) {
        for (const auto& edge : poly1) {
            if (!is_point_in_polygon(edge.first, poly2)) {
                return false;
            }
        }
        return true;
    }

    std::vector<std::pair<GPoint, GPoint>> compute_intersection(
            const std::vector<std::pair<GPoint, GPoint>>& poly1_edges,
            const std::vector<std::pair<GPoint, GPoint>>& poly2_edges) {

        if (is_polygon_contained(poly1_edges, poly2_edges)) {
            return poly1_edges;
        } else if (is_polygon_contained(poly2_edges, poly1_edges)) {
            return poly2_edges;
        }

        std::vector<std::pair<GPoint, GPoint>> intersection_edges;
        std::vector<GPoint> intersection_points;

        for (const auto& edge1 : poly1_edges) {
            for (const auto& edge2 : poly2_edges) {
                GPoint intersection;
                if (segment_intersection(edge1.first, edge1.second, edge2.first, edge2.second, intersection)) {
                    intersection_points.push_back(intersection);
                }
            }
        }

        auto remove_duplicates = [](std::vector<GPoint>& points) {
            std::sort(points.begin(), points.end(), [](const GPoint& a, const GPoint& b) {
                return a.x < b.x || (a.x == b.x && a.y < b.y);
            });
            points.erase(std::unique(points.begin(), points.end(), [](const GPoint& a, const GPoint& b) {
                return std::fabs(a.x - b.x) < 1e-6 && std::fabs(a.y - b.y) < 1e-6;
            }), points.end());
        };
        remove_duplicates(intersection_points);

        for (size_t i = 0; i < intersection_points.size(); ++i) {
            size_t next = (i + 1) % intersection_points.size();
            intersection_edges.emplace_back(intersection_points[i], intersection_points[next]);
        }

        return intersection_edges;
    }

    bool arePointsCollinear(const GPoint *points, int count) {
        if (count < 3) return true;
        const GPoint& A = points[0];
        const GPoint& B = points[1];

        for (size_t i = 2; i < count; ++i) {
            if (std::fabs(cross_product(A, B, points[i])) > 1e-6) {
                return false;
            }
        }

        return true;
    }

    template <typename T>
    T interpolateQuad(const T corners[4], float u, float v) {
        T top = lerp(corners[0], corners[1], u);
        T bottom = lerp(corners[3], corners[2], u);
        return lerp(top, bottom, v);
    }

    template <typename T>
    T lerp(const T& a, const T& b, float t) {
        return a * (1 - t) + b * t;
    }

private:
    GBitmap             bitmap_;
    GMatrix             currentMatrix_;
    std::stack<GMatrix> matrixStack_;
};

#endif //PA4_SHIMULANG_MAIN_MYGCANVAS_HPP
