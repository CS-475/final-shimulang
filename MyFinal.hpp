#ifndef MYFINAL_HPP
#define MYFINAL_HPP

#include "include/GFinal.h"
#include "include/GColor.h"
#include "include/GBitmap.h"
#include "include/GMatrix.h"
#include "include/GMath.h"
#include "include/GShader.h"
#include "include/GPath.h"
#include "include/GPathBuilder.h"
#include "include/GPoint.h"

#include "MyGShader.hpp"

static GPoint evalQuad(const GPoint src[3], float t) {
    float mt = 1 - t;
    float a = mt * mt;
    float b = 2 * mt * t;
    float c = t * t;
    return {
        a*src[0].x + b*src[1].x + c*src[2].x,
        a*src[0].y + b*src[1].y + c*src[2].y
    };
}

static GPoint lerp(const GPoint& A, const GPoint& B, float t) {
    return {A.x + (B.x - A.x)*t, A.y + (B.y - A.y)*t};
}

static GPoint bilinear(const GPoint corners[4], float u, float v){
    const GPoint& P0 = corners[0]; // top-left
    const GPoint& P1 = corners[1]; // top-right
    const GPoint& P2 = corners[2]; // bottom-right
    const GPoint& P3 = corners[3]; // bottom-left

    float mt_u = 1 - u;
    float mt_v = 1 - v;

    return {
        P0.x * mt_u * mt_v + P1.x * u * mt_v + P2.x * u * v + P3.x * mt_u * v,
        P0.y * mt_u * mt_v + P1.y * u * mt_v + P2.y * u * v + P3.y * mt_u * v
    };
}

class MyFinal : public GFinal {
public:
    std::shared_ptr<GShader> createVoronoiShader(const GPoint points[], const GColor colors[], int count) override{
        GMatrix localMatrix = {1, 0, 0,
                               0, 1, 0};
        return std::make_shared<MyVoronoiShader>(points, colors, count, localMatrix);
    }

    std::shared_ptr<GShader> createSweepGradient(GPoint center, float startRadians, const GColor colors[], int count) override {
        return std::make_shared<MySweepGradientShader>(center, startRadians, colors, count);
    }

    std::shared_ptr<GShader> createLinearPosGradient(GPoint p0, GPoint p1,const GColor colors[],const float pos[],int count) override{
        return std::make_shared<MyLinearPosGradientShader>(p0, p1, colors, pos, count);
    }

    std::shared_ptr<GShader> createColorMatrixShader(const GColorMatrix& cm, GShader* realShader) override {
        if (realShader == nullptr) {
            return nullptr;
        }
        return std::make_shared<MyColorMatrixShader>(cm, realShader);
    }
    

    std::pair<float, float> getPerpendicularUnitVector(float dx, float dy) {
        float length = std::sqrt(dx * dx + dy * dy);
        if (length == 0) {
            throw std::invalid_argument("Points must not be the same");
        }
        return {-dy / length, dx / length};
    }

    void expandLineSegment(GPathBuilder& path_builder, const GPoint& p1, const GPoint& p2, float width) {
        float dx = p2.x - p1.x;
        float dy = p2.y - p1.y;

        auto [nx, ny] = getPerpendicularUnitVector(dx, dy);

        float offset_x = nx * (width / 2.0f);
        float offset_y = ny * (width / 2.0f);

        GPoint out1 = {p1.x + offset_x, p1.y + offset_y};
        GPoint out2 = {p1.x - offset_x, p1.y - offset_y};
        GPoint out3 = {p2.x - offset_x, p2.y - offset_y};
        GPoint out4 = {p2.x + offset_x, p2.y + offset_y};

        path_builder.moveTo(out1);
        path_builder.lineTo(out2);
        path_builder.lineTo(out3);
        path_builder.lineTo(out4);
        path_builder.lineTo(out1);
    }

    void addCap(GPathBuilder& pb, const GPoint& p, float width) {
        pb.addCircle(p, width/2.0f);
    }

    void addCorner(GPathBuilder& pb, const GPoint& p, float width) {
        pb.addCircle(p, width/2.0f);
    }
    
    std::shared_ptr<GPath> strokePolygon(const GPoint points[], int count, float width, bool isClosed) override {
        GPathBuilder path_builder;
        if (count < 2) {
            return path_builder.detach();
        }

        for (int i = 0; i < count - 1; ++i) {
            expandLineSegment(path_builder, points[i], points[i+1], width);
            if (isClosed) {
                addCorner(path_builder, points[i+1], width);
            }
        }

        if (isClosed) {
            expandLineSegment(path_builder, points[count - 1], points[0], width);
            addCorner(path_builder, points[0], width);
        } else {
            addCap(path_builder, points[0], width);
            addCap(path_builder, points[count - 1], width);
        }

        return path_builder.detach();
    }

    void drawQuadraticCoons(GCanvas *canvas, const GPoint pts[8], const GPoint tex[4], int level,
                            const GPaint &paint) override {
        // Quadratic Coons patch
        GPoint cornerPts[4] = {pts[0], pts[2], pts[4], pts[6]};

        int rows = level + 1;
        int cols = level + 1;
        int gridCount = rows * cols;

        std::vector<GPoint> gridPoints(gridCount);
        std::vector<GPoint> gridTex(gridCount);

        for (int j = 0; j <= level; j++) {
            float v = (level == 0) ? 0 : (float) j / level;

            GPoint left_src[3] = {pts[0], pts[7], pts[6]};
            GPoint leftPt = evalQuad(left_src, v);

            GPoint right_src[3] = {pts[2], pts[3], pts[4]};
            GPoint rightPt = evalQuad(right_src, v);

            for (int i = 0; i <= level; i++) {
                float u = (level == 0) ? 0 : (float) i / level;

                GPoint top_src[3] = {pts[0], pts[1], pts[2]};
                GPoint topPt = evalQuad(top_src, u);

                GPoint bottom_src[3] = {pts[6], pts[5], pts[4]};
                GPoint bottomPt = evalQuad(bottom_src, u);

                GPoint TB = lerp(topPt, bottomPt, v);
                GPoint LR = lerp(leftPt, rightPt, u);
                GPoint C = bilinear(cornerPts, u, v);

                GPoint finalPt = {TB.x + LR.x - C.x, TB.y + LR.y - C.y};
                gridPoints[j * cols + i] = finalPt;

                GPoint t0 = tex[0];
                GPoint t1 = tex[1];
                GPoint t2 = tex[2];
                GPoint t3 = tex[3];

                float mt_u = 1 - u;
                float mt_v = 1 - v;
                GPoint texPt = {
                    (mt_u * mt_v) * t0.x + (u * mt_v)* t1.x + (u * v) * t2.x + (mt_u * v)* t3.x,
                    (mt_u * mt_v) * t0.y + (u * mt_v)* t1.y + (u * v) * t2.y + (mt_u * v)* t3.y
                };
                gridTex[j * cols + i] = texPt;
            }
        }

        int quadCount = level * level;
        int triangleCount = quadCount * 2;
        int indexCount = triangleCount * 3;
        std::vector<int> indices(indexCount);

        int idx = 0;
        for (int j = 0; j < level; j++) {
            for (int i = 0; i < level; i++) {
                int base = j * (cols) + i;
                indices[idx++] = base;
                indices[idx++] = base + 1;
                indices[idx++] = base + cols;

                indices[idx++] = base + 1;
                indices[idx++] = base + cols + 1;
                indices[idx++] = base + cols;
            }
        }
        canvas->drawMesh(gridPoints.data(), nullptr, gridTex.data(), triangleCount, indices.data(), paint);
    }
};

#endif //MYFINAL_HPP