#include "include/GMatrix.h"
#include "include/nonstd/optional.hpp"

GMatrix::GMatrix() : GMatrix(1, 0, 0, 0, 1, 0) {}

GMatrix GMatrix::Translate(float tx, float ty) {
    return GMatrix(1, 0, tx, 0, 1, ty);
}

GMatrix GMatrix::Scale(float sx, float sy) {
    GMatrix scale = GMatrix(sx, 0, 0, 0, sy, 0);
    return scale;
}

GMatrix GMatrix::Rotate(float radians){
    return GMatrix(cos(radians), -sin(radians), 0, sin(radians), cos(radians), 0);
}

/*
 * Matrix operation
 */
GMatrix GMatrix::Concat(const GMatrix& a, const GMatrix& b) {
    float result_a = a[0]*b[0] + a[2]*b[1];
    float result_b = a[1]*b[0] + a[3]*b[1];
    float result_c = a[0]*b[2] + a[2]*b[3];
    float result_d = a[1]*b[2] + a[3]*b[3];
    float result_e = a[0]*b[4] + a[2]*b[5] + a[4];
    float result_f = a[1]*b[4] + a[3]*b[5] + a[5];

    return GMatrix(result_a, result_c, result_e, result_b, result_d, result_f);
}

nonstd::optional<GMatrix> GMatrix::invert() const {
    // Calculate the determinant of the 2x2 top-left submatrix
    float det = fMat[0] * fMat[3] - fMat[1] * fMat[2];
    if (det == 0) {
        // Matrix is singular
        return nonstd::nullopt;
    }
    float invDet = 1.0f / det;

    float a = fMat[3] * invDet;
    float b = -fMat[1] * invDet;
    float c = -fMat[2] * invDet;
    float d = fMat[0] * invDet;

    float e = (fMat[2] * fMat[5] - fMat[3] * fMat[4]) * invDet;
    float f = (fMat[1] * fMat[4] - fMat[0] * fMat[5]) * invDet;

    return GMatrix(a, c, e, b, d, f);
}

void GMatrix::mapPoints(GPoint dst[], const GPoint src[], int count) const {
    for(int i = 0; i < count; i++){
        float src_x = src[i].x;
        float src_y = src[i].y;
        dst[i].x = fMat[0] * src_x + fMat[2] * src_y + fMat[4];
        dst[i].y = fMat[1] * src_x + fMat[3] * src_y + fMat[5];
    }
}