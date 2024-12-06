#ifndef BLEND_MODES_H
#define BLEND_MODES_H

#include "include/GCanvas.h"
#include "include/GRect.h"
#include "include/GColor.h"
#include "include/GBitmap.h"
#include "include/GPixel.h"
#include <cmath>

inline int div255(int prod) {return (prod + 128) * 257 >> 16;}

// Clear Mode: 0
static GPixel clear_mode(const GPixel& src, const GPixel& dst){
    return GPixel_PackARGB(0, 0, 0, 0);
}

// Src Mode: S
static GPixel src_mode(const GPixel& src, const GPixel& dst){
    return src;
}

// Dst Mode: D
static GPixel dst_mode(const GPixel& src, const GPixel& dst){
    return dst;
}

// SrcOver Mode: S + (1 - Sa)*D
static GPixel srcover_mode(const GPixel& src, const GPixel& dst){
    unsigned sa = GPixel_GetA(src);
    unsigned r = GPixel_GetR(src) + div255((255 - sa) * GPixel_GetR(dst));
    unsigned g = GPixel_GetG(src) + div255((255 - sa) * GPixel_GetG(dst));
    unsigned b = GPixel_GetB(src) + div255((255 - sa) * GPixel_GetB(dst));
    unsigned a = GPixel_GetA(src) + div255((255 - sa) * GPixel_GetA(dst));
    return (GPixel_PackARGB(a, r, g, b));
}

// DstOver Mode: D + (1 - Da)*S
static GPixel dstover_mode(const GPixel& src, const GPixel& dst){
    unsigned da = GPixel_GetA(dst);
    unsigned r = GPixel_GetR(dst) + div255((255 - da) * GPixel_GetR(src));
    unsigned g = GPixel_GetG(dst) + div255((255 - da) * GPixel_GetG(src));
    unsigned b = GPixel_GetB(dst) + div255((255 - da) * GPixel_GetB(src));
    unsigned a = GPixel_GetA(dst) + div255((255 - da) * GPixel_GetA(src));
    return GPixel_PackARGB(a, r, g, b);
}

// SrcIn Mode: Da * S
static GPixel srcin_mode(const GPixel& src, const GPixel& dst){
    unsigned da = GPixel_GetA(dst);
    unsigned r = div255(da * GPixel_GetR(src));
    unsigned g = div255(da * GPixel_GetG(src));
    unsigned b = div255(da * GPixel_GetB(src));
    unsigned a = div255(da * GPixel_GetA(src));
    return GPixel_PackARGB(a, r, g, b);
}

// DstIn Mode: Sa * D
static GPixel dstin_mode(const GPixel& src, const GPixel& dst){
    unsigned sa = GPixel_GetA(src);
    unsigned r = div255(sa * GPixel_GetR(dst));
    unsigned g = div255(sa * GPixel_GetG(dst));
    unsigned b = div255(sa * GPixel_GetB(dst));
    unsigned a = div255(sa * GPixel_GetA(dst));
    return GPixel_PackARGB(a, r, g, b);
}

// SrcOut Mode: (1 - Da)*S
static GPixel srcout_mode(const GPixel& src, const GPixel& dst){
    unsigned da = GPixel_GetA(dst);
    unsigned r = div255((255 - da) * GPixel_GetR(src));
    unsigned g = div255((255 - da) * GPixel_GetG(src));
    unsigned b = div255((255 - da) * GPixel_GetB(src));
    unsigned a = div255((255 - da) * GPixel_GetA(src));
    return GPixel_PackARGB(a, r, g, b);
}

// DstOut Mode: (1 - Sa)*D
static GPixel dstout_mode(const GPixel& src, const GPixel& dst){
    unsigned sa = GPixel_GetA(src);
    unsigned r = div255((255 - sa) * GPixel_GetR(dst));
    unsigned g = div255((255 - sa) * GPixel_GetG(dst));
    unsigned b = div255((255 - sa) * GPixel_GetB(dst));
    unsigned a = div255((255 - sa) * GPixel_GetA(dst));
    return GPixel_PackARGB(a, r, g, b);
}

// SrcATop Mode: Da*S + (1 - Sa)*D
static GPixel srcatop_mode(const GPixel& src, const GPixel& dst){
    unsigned sa = GPixel_GetA(src), da = GPixel_GetA(dst);
    unsigned r = div255(da * GPixel_GetR(src)) + div255((255 - sa) * GPixel_GetR(dst));
    unsigned g = div255(da * GPixel_GetG(src)) + div255((255 - sa) * GPixel_GetG(dst));
    unsigned b = div255(da * GPixel_GetB(src)) + div255((255 - sa) * GPixel_GetB(dst));
    unsigned a = GPixel_GetA(dst); // Keep destination alpha
    return GPixel_PackARGB(a, r, g, b);
}

// DstATop Mode: Sa*D + (1 - Da)*S
static GPixel dstatop_mode(const GPixel& src, const GPixel& dst){
    unsigned sa = GPixel_GetA(src), da = GPixel_GetA(dst);
    unsigned r = div255(sa * GPixel_GetR(dst)) + div255((255 - da) * GPixel_GetR(src));
    unsigned g = div255(sa * GPixel_GetG(dst)) + div255((255 - da) * GPixel_GetG(src));
    unsigned b = div255(sa * GPixel_GetB(dst)) + div255((255 - da) * GPixel_GetB(src));
    unsigned a = div255(sa * GPixel_GetA(dst)) + div255((255 - da) * GPixel_GetA(src));
    return GPixel_PackARGB(a, r, g, b);
}

// Xor Mode: (1 - Sa)*D + (1 - Da)*S
static GPixel xor_mode(const GPixel& src, const GPixel& dst){
    unsigned sa = GPixel_GetA(src), da = GPixel_GetA(dst);
    unsigned r = div255((255 - sa) * GPixel_GetR(dst)) + div255((255 - da) * GPixel_GetR(src));
    unsigned g = div255((255 - sa) * GPixel_GetG(dst)) + div255((255 - da) * GPixel_GetG(src));
    unsigned b = div255((255 - sa) * GPixel_GetB(dst)) + div255((255 - da) * GPixel_GetB(src));
    unsigned a = div255((255 - sa) * GPixel_GetA(dst)) + div255((255 - da) * GPixel_GetA(src));
    return GPixel_PackARGB(a, r, g, b);
}

static GBlendMode blendTable(const GBlendMode mode, int sa){
    if (sa == 0){ // mode 0
        switch (mode) {
            case GBlendMode::kClear:    return GBlendMode::kClear;
            case GBlendMode::kSrc:      return GBlendMode::kClear;
            case GBlendMode::kDst:      return GBlendMode::kDst;
            case GBlendMode::kSrcOver:  return GBlendMode::kDst;
            case GBlendMode::kDstOver:  return GBlendMode::kDst;
            case GBlendMode::kSrcIn:    return GBlendMode::kClear;
            case GBlendMode::kDstIn:    return GBlendMode::kClear;
            case GBlendMode::kSrcOut:   return GBlendMode::kClear;
            case GBlendMode::kDstOut:   return GBlendMode::kDst;
            case GBlendMode::kSrcATop:  return GBlendMode::kDst;
            case GBlendMode::kDstATop:  return GBlendMode::kClear;
            case GBlendMode::kXor:      return GBlendMode::kDst;
        }
    }
    if (sa == 255){ // mode 1
        switch (mode) {
            case GBlendMode::kClear:    return GBlendMode::kClear;
            case GBlendMode::kSrc:      return GBlendMode::kSrc;
            case GBlendMode::kDst:      return GBlendMode::kDst;
            case GBlendMode::kSrcOver:  return GBlendMode::kSrc;
            case GBlendMode::kDstOver:  return GBlendMode::kDstOver;
            case GBlendMode::kSrcIn:    return GBlendMode::kSrcIn;
            case GBlendMode::kDstIn:    return GBlendMode::kDst;
            case GBlendMode::kSrcOut:   return GBlendMode::kSrcOut;
            case GBlendMode::kDstOut:   return GBlendMode::kClear;
            case GBlendMode::kSrcATop:  return GBlendMode::kSrcIn;
            case GBlendMode::kDstATop:  return GBlendMode::kDstATop;
            case GBlendMode::kXor:      return GBlendMode::kXor;
        }
    }

    return mode;
}

static GPixel blend(const GPixel& src, const GPixel& dst, const GBlendMode mode) {
    switch (mode) {
        case GBlendMode::kClear:    return clear_mode(src, dst);
        case GBlendMode::kSrc:      return src_mode(src, dst);
        case GBlendMode::kDst:      return dst_mode(src, dst);
        case GBlendMode::kSrcOver:  return srcover_mode(src, dst);
        case GBlendMode::kDstOver:  return dstover_mode(src, dst);
        case GBlendMode::kSrcIn:    return srcin_mode(src, dst);
        case GBlendMode::kDstIn:    return dstin_mode(src, dst);
        case GBlendMode::kSrcOut:   return srcout_mode(src, dst);
        case GBlendMode::kDstOut:   return dstout_mode(src, dst);
        case GBlendMode::kSrcATop:  return srcatop_mode(src, dst);
        case GBlendMode::kDstATop:  return dstatop_mode(src, dst);
        case GBlendMode::kXor:      return xor_mode(src, dst);

        default:
            return src_mode(src, dst);
    }
}

#endif