#ifndef PA4_MYGSHADER_HPP
#define PA4_MYGSHADER_HPP

#include <iostream>

#include "include/GShader.h"
#include "include/GMatrix.h"
#include "include/GMath.h"
#include "include/GBitmap.h"
#include "include/GFinal.h"

__attribute__((unused))static GPixel ColorToPixel(const GColor &color) {
    unsigned r = std::lround(color.r * color.a * 255);
    unsigned g = std::lround(color.g * color.a * 255);
    unsigned b = std::lround(color.b * color.a * 255);
    unsigned a = std::lround(color.a * 255);
    return GPixel_PackARGB(a, r, g, b);
}

static GColor PixelToUnpremulColor(GPixel p) {
    float a = GPixel_GetA(p)/255.0f;
    if (a <= 0) {
        return {0,0,0,0};
    }
    float r = GPixel_GetR(p)/255.0f / a;
    float g = GPixel_GetG(p)/255.0f / a;
    float b = GPixel_GetB(p)/255.0f / a;
    return {r,g,b,a};
}


class BitmapShader : public GShader, public std::enable_shared_from_this<BitmapShader>{
public:
    BitmapShader(const GBitmap& bm, const GMatrix local_matrix, GTileMode tileMode_)
        : bitmap(bm), localMatrix(local_matrix), tileMode(tileMode_){

    }
    bool isOpaque() override;
    bool setContext(const GMatrix& ctm) override;
    void shadeRow(int x, int y, int count, GPixel row[]) override;

private:
    int applyTileMode(float coord, int max) const;

private:
    GBitmap bitmap;
    GMatrix localMatrix;
    GMatrix invTransform;
    GTileMode tileMode;
};

class LinearGradientShader : public GShader, public std::enable_shared_from_this<LinearGradientShader>{
public:
    LinearGradientShader(const GColor* colors, GMatrix m, int count, GTileMode mode) :
    numColor(count), M(m), shadeColors(new GColor[count]), tileMode(mode){
        for (int i = 0; i < count; i++){
            shadeColors[i] = colors[i];
        }
    }

    bool isOpaque() override;

    bool setContext(const GMatrix& ctm) override;

    void shadeRow(int x, int y, int count, GPixel row[]) override;

private:
    float applyTileMode(float px) const;

private:
    int numColor;
    GMatrix M;
    GMatrix invTransform;
    GColor* shadeColors;
    GTileMode tileMode;

};

class MySweepGradientShader : public GShader, public std::enable_shared_from_this<MySweepGradientShader> {
public:
    MySweepGradientShader(GPoint center, float startRadians, const GColor colors[], int count)
        : fDeviceCenter(center), fStart(startRadians), fCount(count) {
        fColors.assign(colors, colors + count);
    }

    bool isOpaque() override {
        for (int i = 0; i < fCount; i++) {
            if (fColors[i].a < 1.0f) return false;
        }
        return true;
    }

    bool setContext(const GMatrix& ctm) override {
        auto inv = ctm.invert();
        if (!inv.has_value()) return false;
        fInvCTM = inv.value();
        fLocalCenter = fInvCTM * fDeviceCenter;

        return true;
    }

    void shadeRow(int x, int y, int count, GPixel row[]) override {

        float twoPi = 2.0f * (float)M_PI;
        float sectorCount = fCount - 1;

        for (int i = 0; i < count; i++) {
            GPoint src = fInvCTM * GPoint{x + i + 0.5f, y + 0.5f};
            float dx = src.x - fLocalCenter.x;
            float dy = src.y - fLocalCenter.y;

            float angle = (dx == 0 && dy == 0) ? fStart : std::atan2(dy, dx);

            float normAngle = angle - fStart;
            normAngle = fmod(normAngle, twoPi);
            if (normAngle < 0) {
                normAngle += twoPi;
            }

            float t = normAngle / twoPi;

            if (fCount == 1) {
                row[i] = ColorToPixel(fColors[0]);
                continue;
            }

            float x_t = t * sectorCount;
            int idx = (int)std::floor(x_t);
            if (idx >= fCount - 1) idx = fCount - 2;
            float localT = x_t - idx;

            const GColor& c0 = fColors[idx];
            const GColor& c1 = fColors[idx+1];
            GColor c = {c0.r*(1 - localT)+c1.r*localT,
                        c0.g*(1 - localT)+c1.g*localT,
                        c0.b*(1 - localT)+c1.b*localT,
                        c0.a*(1 - localT)+c1.a*localT};

            row[i] = ColorToPixel(c);
        }
    }

private:
    GPoint fDeviceCenter;
    GPoint fLocalCenter;
    float  fStart;
    int    fCount;
    std::vector<GColor> fColors;
    GMatrix fInvCTM;
};


class MyVoronoiShader : public GShader, public std::enable_shared_from_this<MyVoronoiShader>{
public:
    MyVoronoiShader(const GPoint* points, const GColor* colors, int count, const GMatrix& localMatrix)
        : numColor(count), localM(localMatrix) {
        shadeColors.resize(count);
        voronoiPoints.resize(count);
        for (int i = 0; i < count; i++){
            shadeColors[i] = colors[i];
            voronoiPoints[i] = points[i];
        }
    }

    bool isOpaque() override {
        for (int i = 0; i < numColor; i++) {
            if (shadeColors[i].a < 1) return false;
        }
        return true;
    }

    bool setContext(const GMatrix& ctm) override {
        auto inv = (ctm * localM).invert();
        if (inv.has_value()) {
            invTransform = inv.value();
            return true;
        }
        return false;
    }

    void shadeRow(int x, int y, int count, GPixel row[]) override {
        for (int i = 0; i < count; i++) {
            GPoint src = invTransform * GPoint{x + i + 0.5f, y + 0.5f};
            int closestIndex = 0;
            float minDist = std::numeric_limits<float>::max();
            for (int k = 0; k < numColor; k++){
                float dx = src.x - voronoiPoints[k].x;
                float dy = src.y - voronoiPoints[k].y;
                float dist = dx*dx + dy*dy;
                if (dist < minDist) {
                    minDist = dist;
                    closestIndex = k;
                }
            }
            row[i] = ColorToPixel(shadeColors[closestIndex]);
        }
    }

private:
    int numColor;
    GMatrix localM;
    GMatrix invTransform;
    std::vector<GColor> shadeColors;
    std::vector<GPoint> voronoiPoints;
};

class MyLinearPosGradientShader : public GShader {
public:
    MyLinearPosGradientShader(GPoint p0, GPoint p1, const GColor colors[], const float pos[], int count)
        : fCount(count), fP0(p0), fP1(p1) {
        fColors.assign(colors, colors+count);
        fPos.assign(pos, pos+count);

        fLinePts.resize(count);
        for (int i = 0; i < count; i++) {
            float t = fPos[i];
            float mt = 1 - t;
            fLinePts[i] = { mt*p0.x + t*p1.x, mt*p0.y + t*p1.y };
        }

        float dx = p1.x - p0.x;
        float dy = p1.y - p0.y;
        fLen = std::sqrt(dx*dx + dy*dy);
        if (fLen > 0) {
            fDX = dx / (fLen * fLen);
            fDY = dy / (fLen * fLen);
        } else {
            fDX = fDY = 0;
        }
    }

    bool isOpaque() override {
        for (auto& c : fColors) {
            if (c.a < 1.0f) {
                return false;
            }
        }
        return true;
    }

    bool setContext(const GMatrix& ctm) override {
        auto inv = ctm.invert();
        if (!inv.has_value()) return false;
        fInvCTM = inv.value();
        return true;
    }

    void shadeRow(int x, int y, int count, GPixel row[]) override {

        for (int i = 0; i < count; i++) {
            GPoint srcPt = fInvCTM * GPoint{x + i + 0.5f, y + 0.5f};

            float tx = srcPt.x - fP0.x;
            float ty = srcPt.y - fP0.y;
            float t = tx * (fP1.x - fP0.x) * (1/(fLen*fLen)) + ty * (fP1.y - fP0.y) * (1/(fLen*fLen));

            t = std::clamp(t, 0.0f, 1.0f);

            int idx = findInterval(t);
            float pStart = fPos[idx];
            float pEnd   = fPos[idx+1];

            float localT = (t - pStart) / (pEnd - pStart);
            GColor c = lerpColor(fColors[idx], fColors[idx+1], localT);
            row[i] = ColorToPixel(c);
        }
    }

private:
    int findInterval(float t) const {
        int low = 0;
        int high = fCount - 2;
        while (low <= high) {
            int mid = (low + high)/2;
            if (fPos[mid] <= t && t < fPos[mid+1]) {
                return mid;
            } else if (t < fPos[mid]) {
                high = mid - 1;
            } else {
                low = mid + 1;
            }
        }
        return fCount - 2;
    }

    GColor lerpColor(const GColor& c0, const GColor& c1, float t) const {
        return {c0.r*(1-t)+c1.r*t,
                c0.g*(1-t)+c1.g*t,
                c0.b*(1-t)+c1.b*t,
                c0.a*(1-t)+c1.a*t};
    }

private:
    int fCount;
    GPoint fP0, fP1;
    float fLen;
    float fDX, fDY;
    GMatrix fInvCTM;
    std::vector<GColor> fColors;
    std::vector<float> fPos;
    std::vector<GPoint> fLinePts;
};

class MyColorMatrixShader : public GShader, public std::enable_shared_from_this<MyColorMatrixShader> {
public:
    MyColorMatrixShader(const GColorMatrix& cm, GShader* real)
        : fCM(cm), fRealShader(real) {}

    bool isOpaque() override {
        return false;
    }

    bool setContext(const GMatrix& ctm) override {
        return fRealShader->setContext(ctm);
    }

    void shadeRow(int x, int y, int count, GPixel row[]) override {
        std::vector<GPixel> tmp(count);
        fRealShader->shadeRow(x, y, count, tmp.data());

        for (int i = 0; i < count; i++) {
            GColor c = PixelToUnpremulColor(tmp[i]);

            float r = fCM[0]*c.r + fCM[4]*c.g + fCM[8]*c.b  + fCM[12]*c.a + fCM[16];
            float g = fCM[1]*c.r + fCM[5]*c.g + fCM[9]*c.b  + fCM[13]*c.a + fCM[17];
            float b = fCM[2]*c.r + fCM[6]*c.g + fCM[10]*c.b + fCM[14]*c.a + fCM[18];
            float a = fCM[3]*c.r + fCM[7]*c.g + fCM[11]*c.b + fCM[15]*c.a + fCM[19];

            r = std::clamp(r, 0.0f, 1.0f);
            g = std::clamp(g, 0.0f, 1.0f);
            b = std::clamp(b, 0.0f, 1.0f);
            a = std::clamp(a, 0.0f, 1.0f);

            GColor newC = {r,g,b,a};
            row[i] = ColorToPixel(newC);
        }
    }

private:
    GColorMatrix fCM;
    GShader* fRealShader;
};

#endif //PA4_MYGSHADER_HPP
