#include <iostream>
#include "include/GPath.h"

static float find_quad_extreme(float p0, float p1, float p2) {
    if (p0 == p2) return p0;
    float t = (p0 - p1) / (p0 - 2 * p1 + p2);
    if (t > 0 && t < 1) {
        return (1 - t) * (1 - t) * p0 + 2 * (1 - t) * t * p1 + t * t * p2;
    }
    return std::min(std::min(p0, p1), p2);
}

static void find_cubic_extremes(float p0, float p1, float p2, float p3, std::vector<float>& extremes) {
    float a = -p0 + 3 * p1 - 3 * p2 + p3;
    float b = 2 * (p0 - 2 * p1 + p2);
    float c = -p0 + p1;

    if (std::fabs(a) < std::numeric_limits<float>::epsilon()) {
        if (std::fabs(b) >= std::numeric_limits<float>::epsilon()) {
            float t = -c / b;
            if (t > 0 && t < 1) {
                extremes.push_back(t);
            }
        }
    } else {
        float discriminant = b * b - 4 * a * c;
        if (discriminant >= 0) {
            float sqrt_disc = std::sqrt(discriminant);
            float t1 = (-b + sqrt_disc) / (2 * a);
            float t2 = (-b - sqrt_disc) / (2 * a);
            if (t1 > 0 && t1 < 1) extremes.push_back(t1);
            if (t2 > 0 && t2 < 1) extremes.push_back(t2);
        }
    }
}

GRect GPath::bounds() const {
    if (fPts.empty())
        return GRect::LTRB(0, 0, 0, 0);

    float left = std::numeric_limits<float>::max();
    float top = std::numeric_limits<float>::max();
    float right = std::numeric_limits<float>::lowest();
    float bottom = std::numeric_limits<float>::lowest();

    auto updateBounds = [&](float x, float y) {
        left = std::min(left, x);
        top = std::min(top, y);
        right = std::max(right, x);
        bottom = std::max(bottom, y);
    };

    bool curve_detect = false;
    Iter iter(*this);
    GPoint *pts = new GPoint[GPath::kMaxNextPoints];
    std::vector<GPoint> all_points;
    while (auto verb = iter.next(pts)) {
        switch (verb.value()) {
            case GPathVerb::kLine:{
                for (int i = 0; i < 2; ++i) {
                    all_points.push_back(pts[i]);
                }
            };break;
            case GPathVerb::kMove:{
            };break;
            case GPathVerb::kQuad:{
                curve_detect = true;
                float x_extreme = find_quad_extreme(pts[0].x, pts[1].x, pts[2].x);
                float y_extreme = find_quad_extreme(pts[0].y, pts[1].y, pts[2].y);
                updateBounds(pts[0].x, pts[0].y);
                updateBounds(pts[2].x, pts[2].y);
                updateBounds(x_extreme, y_extreme);
            };break;
            case GPathVerb::kCubic:{
                curve_detect = true;
                std::vector<float> x_extremes, y_extremes;
                find_cubic_extremes(pts[0].x, pts[1].x, pts[2].x, pts[3].x, x_extremes);
                find_cubic_extremes(pts[0].y, pts[1].y, pts[2].y, pts[3].y, y_extremes);

                updateBounds(pts[0].x, pts[0].y);
                updateBounds(pts[3].x, pts[3].y);

                for (float t : x_extremes) {
                    float x = (1 - t) * (1 - t) * (1 - t) * pts[0].x +
                              3 * (1 - t) * (1 - t) * t * pts[1].x +
                              3 * (1 - t) * t * t * pts[2].x +
                              t * t * t * pts[3].x;
                    updateBounds(x, pts[0].y);
                }

                for (float t : y_extremes) {
                    float y = (1 - t) * (1 - t) * (1 - t) * pts[0].y +
                              3 * (1 - t) * (1 - t) * t * pts[1].y +
                              3 * (1 - t) * t * t * pts[2].y +
                              t * t * t * pts[3].y;
                    updateBounds(pts[0].x, y);
                }
            };break;
        }
    }
    delete[] pts;

    if(!curve_detect){
        for (const GPoint &pt: fPts) {
            left = std::min(left, pt.x);
            top = std::min(top, pt.y);
            right = std::max(right, pt.x);
            bottom = std::max(bottom, pt.y);
        }
        return GRect::LTRB(left, top, right, bottom);
    }

    return GRect::LTRB(left, top, right, bottom);
}

static GPoint lerp(const GPoint& P0, const GPoint& P1, float t) {
    return GPoint{
            P0.x + (P1.x - P0.x) * t,
            P0.y + (P1.y - P0.y) * t
    };
}

void GPath::ChopQuadAt(const GPoint *src, GPoint *dst, float t) {
//    const GPoint& P0 = src[0];
//    const GPoint& P1 = src[1];
//    const GPoint& P2 = src[2];
//
//    dst[0] = BezierQuad(P0, P1, P2, 0);
//    dst[1] = BezierQuad(P0, P1, P2, 0.5f * t);
//    dst[2] = BezierQuad(P0, P1, P2, t);
//    dst[3] = BezierQuad(P0, P1, P2, t + (1.f-t) * 0.5f);
//    dst[4] = BezierQuad(P0, P1, P2, 1);
    const GPoint& P0 = src[0];
    const GPoint& P1 = src[1];
    const GPoint& P2 = src[2];

    GPoint Q0 = lerp(P0, P1, t);
    GPoint Q1 = lerp(P1, P2, t);
    GPoint R = lerp(Q0, Q1, t);

    dst[0] = P0;
    dst[1] = Q0;
    dst[2] = R;
    dst[3] = Q1;
    dst[4] = P2;
}

void GPath::ChopCubicAt(const GPoint *src, GPoint *dst, float t) {
//    const GPoint& P0 = src[0];
//    const GPoint& P1 = src[1];
//    const GPoint& P2 = src[2];
//    const GPoint& P3 = src[3];
//
//    dst[0] = BezierCubic(P0, P1, P2, P3, 0);
//    dst[1] = BezierCubic(P0, P1, P2, P3, 0.33f * t);
//    dst[2] = BezierCubic(P0, P1, P2, P3, 0.66f * t);
//    dst[3] = BezierCubic(P0, P1, P2, P3, t);
//    dst[4] = BezierCubic(P0, P1, P2, P3, t + 0.33f * (1.f-t));
//    dst[5] = BezierCubic(P0, P1, P2, P3, t + 0.66f * (1.f-t));
//    dst[6] = BezierCubic(P0, P1, P2, P3, 1);

    const GPoint& P0 = src[0];
    const GPoint& P1 = src[1];
    const GPoint& P2 = src[2];
    const GPoint& P3 = src[3];

    GPoint Q0 = lerp(P0, P1, t);
    GPoint Q1 = lerp(P1, P2, t);
    GPoint Q2 = lerp(P2, P3, t);

    GPoint R0 = lerp(Q0, Q1, t);
    GPoint R1 = lerp(Q1, Q2, t);

    GPoint S = lerp(R0, R1, t);

    dst[0] = P0;
    dst[1] = Q0;
    dst[2] = R0;
    dst[3] = S;

    dst[4] = R1;
    dst[5] = Q2;
    dst[6] = P3;
}