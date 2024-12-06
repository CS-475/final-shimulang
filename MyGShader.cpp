#include "include/GShader.h"
#include "MyGShader.hpp"

bool BitmapShader::isOpaque() {
    for (int y = 0; y < bitmap.height(); ++y) {
        for (int x = 0; x < bitmap.width(); ++x) {
            GPixel pixel = *bitmap.getAddr(x, y);
            if (GPixel_GetA(pixel) != 255) {
                return false;
            }
        }
    }
    return true;
}

bool BitmapShader::setContext(const GMatrix &ctm) {
    auto inv = (ctm * localMatrix).invert();
    if (inv.has_value()){
        invTransform = inv.value();
        return true;
    }
    else
        return false;
}

void BitmapShader::shadeRow(int x, int y, int count, GPixel *row) {
    GPoint loc = invTransform * GPoint{x + 0.5f, y + 0.5f};
    GVector step{
            invTransform[0],
            invTransform[1]
    };

    // Bitmap dimensions
    int width = bitmap.width();
    int height = bitmap.height();
    for (int i = 0; i < count; ++i) {
        // Apply tiling mode for x and y coordinates
        int u = applyTileMode(loc.x, width);
        int v = applyTileMode(loc.y, height);

        // Check if the coordinates are valid
        if (u >= 0 && u < width && v >= 0 && v < height) {
            // Map the pixel using the bitmap
            row[i] = *bitmap.getAddr(u, v);
        } else {
            // Assign transparent pixel if out of bounds
            row[i] = GPixel_PackARGB(0, 0, 0, 0);
        }

        // Increment the location by the step vector
        loc.x += step.x;
        loc.y += step.y;
    }
}

int BitmapShader::applyTileMode(float coord, int max) const {
    switch (tileMode) {
        case GTileMode::kClamp:
            return std::max(0, std::min(static_cast<int>(coord), max - 1));
        case GTileMode::kRepeat: {
            float modCoord = std::fmod(coord, (float)max);
            if (modCoord < 0) {
                modCoord += max;
            }
            return static_cast<int>(modCoord);
        }
        case GTileMode::kMirror: {
            int coordInt = static_cast<int>(std::floor(coord));
            int mirrorIndex = coordInt % (2 * max);
            if (mirrorIndex < 0) mirrorIndex += 2 * max;
            return (mirrorIndex >= max) ? (2 * max - 1 - mirrorIndex) : mirrorIndex;
        }
        default:
            return 0;
    }
}

std::shared_ptr<GShader> GCreateBitmapShader(const GBitmap& bm, const GMatrix& localMatrix, GTileMode gTileMode){
    return std::make_shared<BitmapShader>(bm, localMatrix, gTileMode);
}
//!----------------------------------------------------------------------------------------------------

bool LinearGradientShader::isOpaque() {
    for (int i = 0; i < numColor; i++) {
        if (shadeColors[i].a < 1) {
            return false;
        }
    }
    return true;
}

bool LinearGradientShader::setContext(const GMatrix &ctm) {
    auto inv = (ctm * M).invert();
    if (inv.has_value()){
        invTransform = inv.value();
        return true;
    }
    else
        return false;
}

void LinearGradientShader::shadeRow(int x, int y, int count, GPixel *row) {
    GPoint p_t = invTransform * GPoint{x + 0.5f, y + 0.5f};
    float px = p_t.x;
    float dpx = invTransform.e0().x;

    float x_t;
    int k;
    float t;
    GColor c_t;

    if (numColor == 1) {
        c_t = shadeColors[0];
        GPixel src = ColorToPixel(c_t);
        std::fill_n(row, count, src);
        return;
    }

    if (numColor == 2) {
        GColor c0 = shadeColors[0];
        GColor c1 = shadeColors[1];
        for (int i = 0; i < count; i++, px += dpx) {
            x_t = applyTileMode(px);
            c_t = (1 - x_t) * c0 + x_t * c1;
            row[i] = ColorToPixel(c_t);
        }
        return;
    }

    for (int i = 0; i < count; i++, px += dpx) {
        x_t = applyTileMode(px) * (numColor - 1);
        k = static_cast<int>(floor(x_t));
        t = x_t - k;
        c_t = (1 - t) * shadeColors[k] + t * shadeColors[k + 1];
        row[i] = ColorToPixel(c_t);
    }
}

float LinearGradientShader::applyTileMode(float px) const {
    switch (tileMode) {
        case GTileMode::kClamp:
            return std::clamp(px, 0.0f, 1.0f);
        case GTileMode::kRepeat:
            return px - floor(px);
        case GTileMode::kMirror:
            px = fabs(px);
            int integerPart = static_cast<int>(floor(px));
            return (integerPart % 2 == 0) ? (px - integerPart) : (1.0f - (px - integerPart));
    }
    return 0.0f;
}




std::shared_ptr<GShader> GCreateLinearGradient(GPoint p0, GPoint p1, const GColor colors[], int count,
                                               GTileMode tileMode){
    if (count >= 1) {
        int dx = GRoundToInt(p1.x - p0.x);
        int dy = GRoundToInt(p1.y - p0.y);
        GMatrix m = GMatrix(dx, -dy, p0.x, dy, dx, p0.y);
        return std::make_shared<LinearGradientShader>(colors, m, count,tileMode);

    }
    else {
        std::cout << "nullptr" <<  std::endl;
        return nullptr;
    }
}


