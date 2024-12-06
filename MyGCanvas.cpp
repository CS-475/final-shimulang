#include "MyGCanvas.hpp"
#include <iostream>

MyGCanvas::MyGCanvas(const GBitmap &bitmap) : bitmap_(bitmap), currentMatrix_(GMatrix()) {
    // init
    matrixStack_.push(currentMatrix_);
}

void MyGCanvas::save() {
    // Update Stack
    if(!matrixStack_.empty()){
        matrixStack_.pop();
        matrixStack_.push(currentMatrix_);
    }
    else{
        matrixStack_.push(currentMatrix_);
    }
}

void MyGCanvas::restore() {
    // recover save matrix
    // Update currentMatrix_
    if(!matrixStack_.empty()){
        currentMatrix_ = matrixStack_.top();
        matrixStack_.pop();
    }
    else{
        //! Identity Matrix
        currentMatrix_ = GMatrix();
    }
}

void MyGCanvas::concat(const GMatrix &matrix) {
    // a = a * b
//     std::cout << "concat" << std::endl;
    currentMatrix_ = GMatrix::Concat(currentMatrix_, matrix);

//    std::cout << currentMatrix_.e0().x << ", " << currentMatrix_.e0().y << ", " << currentMatrix_.origin().x << std::endl;
//    std::cout << currentMatrix_.e1().x << ", " << currentMatrix_.e1().y << ", " << currentMatrix_.origin().y << std::endl;
}

void MyGCanvas::clear(const GColor &color) {
    // Fill the entire canvas with the specified color, using kSrc porter-duff mode.
    GPixel pixel = ColorToPixel(color);

    int width = bitmap_.width();
    int height = bitmap_.height();

    if( width == 0 || height == 0) {
        printf("Error Canvas\n");
        return;
    }

    for(int w = 0; w < width; w++) {
        for(int h = 0; h < height; h++) {
            *bitmap_.getAddr(w,h) = pixel;
        }
    }
}

//----------------------------------------------------------
void MyGCanvas::drawConvexPolygon(const GPoint *points, int count, const GPaint &paint) {

    if(count < 3)
        //! point values not important given the count(pa2)
        return;
    else {
        //! no area(pa2)
        if(arePointsCollinear(points, count))
            return;
        //! clipped out(pa2)
    }

    // return;
    GPoint *transformedPoints = new GPoint[count];
    currentMatrix_.mapPoints(transformedPoints, points, count);
    GBlendMode blendMode = paint.getBlendMode();
    bool isOpaque;

    auto sh = paint.shareShader();
//    int sh_weight = 1;
    if (sh != nullptr) {
        sh->setContext(currentMatrix_);
        // reset the blendMode if is opaque
        isOpaque = sh->isOpaque();
        if (isOpaque) {
            blendMode = blendTable(blendMode, 255);
        }

//        if (auto myShaderPtr = std::dynamic_pointer_cast<MyShader>(sh)) {
//            sh_weight = 1;
//        } else if (auto myGShaderPtr = std::dynamic_pointer_cast<MyGShader>(sh)) {
//            sh_weight = 2;
//        }


    }

    GPathBuilder gPathBuilder;
    gPathBuilder.addPolygon(transformedPoints, count);
    auto path = gPathBuilder.detach();
    GPath::Edger path_edger(*path);
    delete[] transformedPoints;

    GPoint *pts = new GPoint[GPath::kMaxNextPoints];
    std::vector<std::pair<GPoint, GPoint>> raw_polygon;
    while (auto verb = path_edger.next(pts)) {
        switch (verb.value()) {
            case GPathVerb::kLine:{
                raw_polygon.push_back(std::make_pair(pts[0], pts[1]));
            };break;
            case GPathVerb::kMove:{
                printf("Error");
            };break;
        }
    }
    delete[] pts;
    if(raw_polygon.empty()) {
        std::cout << "count : " << count << ", raw_polygon.size() " << raw_polygon.size() << std::endl;
        return;
    }

    auto edges = createPolygon(raw_polygon);

    bool do_flag = false;
    if(edges.size() < 3) {
//        std::cout << "edges.size() : " << edges.size() << std::endl;
//        return;
        if(raw_polygon.size() < 3) {
            // std::cout << "count : " << count << ", edges.size() : " <<edges.size() << ", raw_polygon : "  << raw_polygon.size() << std::endl;
            return;
        }
        edges = raw_polygon;
        do_flag = true;
    }
    GPixel src = ColorToPixel(paint.getColor());
    std::vector<GPoint> all_pixels = fillPolygon(edges);
    if(do_flag){
        //172755
//        std::cout << all_pixels.size() << std::endl;
//        std::cout << "w : " << bitmap_.width() << std::endl;
//        std::cout << "h : " << bitmap_.height() << std::endl;
//        std::cout << bitmap_.width()  * bitmap_.height() << std::endl;
//        if( all_pixels.size() == bitmap_.width()  * bitmap_.height())
//            return;
    }

    for (int i = 0; i < all_pixels.size(); i++) {
        int x = static_cast<int>(all_pixels[i].x);
        int y = static_cast<int>(all_pixels[i].y);
        if(sh == nullptr) {
            drawRow(x, y, 1, src, paint, bitmap_);
        }
        else {
            GPixel tnp[1];
            sh->shadeRow(x, y, 1, tnp);
            GPixel *p = bitmap_.getAddr(x, y);
            *p = blend(*tnp, *p, blendMode);
        }
    }
}

void MyGCanvas::drawRect(const GRect &rect, const GPaint &paint) {

    auto sh = paint.shareShader();
    GBlendMode blendMode = paint.getBlendMode();
    bool isOpaque;
    if (sh != nullptr) {
        sh->setContext(currentMatrix_);
        isOpaque = sh->isOpaque();
        if (isOpaque) {
            blendMode = blendTable(blendMode, 255);
        }
    }

    if (currentMatrix_ != GMatrix()) {
        // std::cout << "currentMatrix_" << std::endl;
        GPoint points[4] = {{rect.left, rect.top},
                            {rect.left, rect.bottom},
                            {rect.right, rect.bottom},
                            {rect.right, rect.top}};
        drawConvexPolygon(points, 4, paint);
        return;
    }

    GIRect ra = intersect(rect, GIRect::XYWH(0, 0, bitmap_.width(), bitmap_.height()));
    int width = ra.width();

    if (sh == nullptr) {
        GPixel src = ColorToPixel(paint.getColor());
        for (int y = ra.top; y < ra.bottom; y++) {
            drawRow(ra.left, y, width, src, paint, bitmap_);
        }
    } else {
        for (int y = ra.top; y < ra.bottom; y++) {
            GPixel tnp[width];
            sh->shadeRow(ra.left, y, width, tnp);
            for (int x = ra.left; x < ra.left + width; x++) {
                GPixel *p = bitmap_.getAddr(x, y);
                *p = blend(tnp[x - ra.left], *p, blendMode);
            }
        }
    }
}

void MyGCanvas::drawPath(const GPath &path, const GPaint &paint) {

    auto new_path = path.shared_from_this();
    new_path->transform(currentMatrix_);
    GRect bound = new_path->bounds();
    GBlendMode blendMode = paint.getBlendMode();
    bool isOpaque;
    auto sh = paint.shareShader();
    if (sh != nullptr){
        sh->setContext(currentMatrix_);
        // reset the blendMode if is opaque
        isOpaque = sh->isOpaque();
        if (isOpaque){
            blendMode = blendTable(blendMode, 255);
        }
    }

    GPoint *pts = new GPoint[GPath::kMaxNextPoints];
    std::vector<std::pair<GPoint, GPoint>> raw_polygon;
    GPath::Edger path_edger(*new_path);
//    GPath::Edger path_edger(*transformedPath);
    while (auto verb = path_edger.next(pts)) {
        switch (verb.value()) {
            case GPathVerb::kLine:{
                if(pts[0] != pts[1]) {
                    raw_polygon.push_back(std::make_pair(pts[0], pts[1]));
                }
            };break;
            case GPathVerb::kMove:{
                printf("Error\n");
            };break;
            case GPathVerb::kQuad:{
//                GPoint *rst = new GPoint[5];
//                GPath::ChopQuadAt(pts, rst, 0.5);
//                for (int i = 0; i < 4; ++i) {
//                    raw_polygon.push_back(std::make_pair(rst[i], rst[i+1]));
//                }
//                delete[] rst;
                const GPoint &P0 = pts[0];
                const GPoint &P1 = pts[1];
                const GPoint &P2 = pts[2];

                for(int i = 0 ; i < 20; i++){
                    auto p1 = BezierQuad(P0, P1, P2, float(0.05 * i));
                    auto p2 = BezierQuad(P0, P1, P2, float(0.05 * (i + 1)));
                    raw_polygon.push_back(std::make_pair(p1, p2));
                }
            };break;
            case GPathVerb::kCubic:{
                const GPoint &P0 = pts[0];
                const GPoint &P1 = pts[1];
                const GPoint &P2 = pts[2];
                const GPoint &P3 = pts[3];

                for(int i = 0 ; i < 20; i++){
                    auto p1 = BezierCubic(P0, P1, P2, P3, float(0.05 * i));
                    auto p2 = BezierCubic(P0, P1, P2, P3, float(0.05 * (i + 1)));
                    raw_polygon.push_back(std::make_pair(p1, p2));
                }

//                std::cout << "kCubic" << std::endl;
//                GPoint *rst = new GPoint[7];
//                GPath::ChopCubicAt(pts, rst, 0.5);
//                for (int i = 0; i < 6; ++i) {
//                    raw_polygon.push_back(std::make_pair(rst[i], rst[i+1]));
//                }
//                delete[] rst;
            };break;
        }
    }
    delete[] pts;

//    bool flg = false;
//    std::vector<std::pair<GPoint, GPoint>> new_raw_polygon;
//    for(int i = 0; i < raw_polygon.size() -1; i++){
//        new_raw_polygon.push_back(raw_polygon[i]);
//        if(raw_polygon[i].second != raw_polygon[i+1].first){
//            new_raw_polygon.push_back(std::make_pair(raw_polygon[i].second, raw_polygon[i+1].first));
////            std::cout << "!= " << i << std::endl;
//            flg = true;
//        }
//    }
//    if(flg){
////        new_raw_polygon.push_back(raw_polygon.back());
////        if(new_raw_polygon.back().second != new_raw_polygon.front().first){
////            new_raw_polygon.push_back(std::make_pair(new_raw_polygon.back().second, new_raw_polygon.front().first));
////        }
////        raw_polygon = new_raw_polygon;
//        std::cout << "----flg---------" << std::endl;
//         for(const auto p : raw_polygon){
//            std::cout << p.first.x << ", " <<  p.first.y << " <-> " << p.second.x << ", " <<  p.second.y << std::endl;
//         }
//    }

//    for(const auto p : raw_polygon){
//        std::cout << p.first.x << ", " <<  p.first.y << " <-> " << p.second.x << ", " <<  p.second.y << std::endl;
//    }

    if(raw_polygon.size() <= 2)
        return;

    for(auto& edge : raw_polygon){
        GPoint p1 = edge.first;
        GPoint p2 = edge.second;
        //! same point
        if(p1 == p2) continue;
        currentMatrix_.mapPoints(&edge.first, &p1, 1);
        currentMatrix_.mapPoints(&edge.second, &p2, 1);
        // std::cout << p1.x << ", " << p1.y << "<->" << p2.x << ", " << p2.y << std::endl;
    }

    GPixel src = ColorToPixel(paint.getColor());

    std::vector<std::vector<std::pair<GPoint, GPoint>>> mult_polygon;
    std::vector<std::pair<GPoint, GPoint>> one_polygon;
    for(auto& edge : raw_polygon){
        GPoint p1 = edge.first;
        GPoint p2 = edge.second;
        if(p1 == p2) {
            // std::cout << "continue" << std::endl;
            continue;
        }
        if(p1.x == p2.x && p1.y == p2.y) continue;
        // if(p1.x == 0 && p1.y == 0) continue;
        // std::cout << p1.x << ", " << p1.y << "<->" << p2.x << ", " << p2.y << std::endl;
        one_polygon.push_back(std::make_pair(p1, p2));
        if(!one_polygon.empty()) {
            if(p2 == one_polygon.begin()->first) {
                mult_polygon.push_back(one_polygon);
                one_polygon.clear();
                // std::cout << "---------" << std::endl;
            }
        }
    }

    std::vector<GPoint> all_pixels;
    for(auto single_polygon : mult_polygon) {
        auto sp = fillPolygon(single_polygon);
        // all_pixels.push_back()
        all_pixels.insert(all_pixels.end(), sp.begin(), sp.end());

    }

    // std::vector<GPoint> all_pixels = fillPolygon(raw_polygon);



    for (int i = 0; i < all_pixels.size(); i++) {
        int x = static_cast<int>(all_pixels[i].x);
        int y = static_cast<int>(all_pixels[i].y);
        if(sh == nullptr) {
            drawRow(x, y, 1, src, paint, bitmap_);
        }
        else {
            GPixel tnp[1];
            sh->shadeRow(x, y, 1, tnp);
            GPixel *p = bitmap_.getAddr(x, y);
            *p = blend(*tnp, *p, blendMode);
        }
    }
}

void MyGCanvas::drawMesh(const GPoint *verts, const GColor *colors, const GPoint *texs, int count, const int *indices,
                         const GPaint &paint) {
//    std::cout << "count : " << count << std::endl;
    if(count <= 0)
        return;

    GPoint *transformedPoints = new GPoint[count*3];
    currentMatrix_.mapPoints(transformedPoints, verts, count*3);

    GPoint *transformedPointsTex = new GPoint[count * 3];
    if(texs != nullptr) {
//        GPoint *transformedPointsTex = new GPoint[count * 3];
        currentMatrix_.mapPoints(transformedPointsTex, texs, count * 3);
//        *transformedPointsTex = *texs;
    }

    GBlendMode blendMode = paint.getBlendMode();
    bool isOpaque = false;
    auto shader = paint.shareShader();

    if (shader != nullptr) {
//        std::cout << currentMatrix_.e0().x << ", " << currentMatrix_.e0().y << ", " << currentMatrix_.origin().x << std::endl;
//        std::cout << currentMatrix_.e1().x << ", " << currentMatrix_.e1().y << ", " << currentMatrix_.origin().y << std::endl;
        shader->setContext(currentMatrix_);
        isOpaque = shader->isOpaque();
        if (isOpaque) {
            blendMode = blendTable(blendMode, 255);
        }
    }
    else{
//        std::cout << "null : " << currentMatrix_.e0().x << ", " << currentMatrix_.e0().y << ", " << currentMatrix_.origin().x << std::endl;
//        std::cout << "null : " << currentMatrix_.e1().x << ", " << currentMatrix_.e1().y << ", " << currentMatrix_.origin().y << std::endl;
    }

    auto PointInTriangles = [this](const GPoint& p0, const GPoint& p1, const GPoint& p2, const GPoint& p) -> bool {
        double cross1 = cross_product(p0, p1, p);
        double cross2 = cross_product(p1, p2, p);
        double cross3 = cross_product(p2, p0, p);
        return (cross1 >= 0 && cross2 >= 0 && cross3 >= 0) || (cross1 <= 0 && cross2 <= 0 && cross3 <= 0);
    };
    // edge function to calculate if points are inside the triangle
    auto edgeFunction = [](const GPoint& a, const GPoint& b, const GPoint& c) {
        return (c.x - a.x) * (b.y - a.y) - (c.y - a.y) * (b.x - a.x);
    };

    for (int i = 0; i < count * 3; i += 3) {
        //! Vertex
        GPoint p0 = transformedPoints[indices[i]];
        GPoint p1 = transformedPoints[indices[i + 1]];
        GPoint p2 = transformedPoints[indices[i + 2]];

        //! Color
        GColor c0 = colors ? colors[indices[i]] : GColor{1, 1, 1, 1};
        GColor c1 = colors ? colors[indices[i + 1]] : GColor{1, 1, 1, 1};
        GColor c2 = colors ? colors[indices[i + 2]] : GColor{1, 1, 1, 1};

        //! Texture
        GPoint t0 = texs ? transformedPointsTex[indices[i]] : GPoint{0, 0};
        GPoint t1 = texs ? transformedPointsTex[indices[i + 1]] : GPoint{0, 0};
        GPoint t2 = texs ? transformedPointsTex[indices[i + 2]] : GPoint{0, 0};

        float minX = std::min({p0.x, p1.x, p2.x});
        float minY = std::min({p0.y, p1.y, p2.y});
        float maxX = std::max({p0.x, p1.x, p2.x});
        float maxY = std::max({p0.y, p1.y, p2.y});

        int xStart = static_cast<int>(std::floor(minX));
        int yStart = static_cast<int>(std::floor(minY));
        int xEnd = static_cast<int>(std::ceil(maxX));
        int yEnd = static_cast<int>(std::ceil(maxY));

        //! area is zero
        if (p0 == p1 || p0 == p2 || p1 == p2) {
            std::cout << "zero area" << std::endl;
            continue;
        }

        bool same_color = false;
        if(c0 == c1 || c0 == c2 || c1 == c2) {
            same_color = true;
        }

        for (int y = yStart; y <= yEnd; ++y) {
            for (int x = xStart; x <= xEnd; ++x) {
                GPoint pixel = {x + 0.5f, y + 0.5f};

                if (PointInTriangles(p0, p1, p2, pixel)) {
                    // calculate cross product position
                    float area = cross_product(p0, p1, p2);
                    float w0 = cross_product(p1, p2, pixel) / area;
                    float w1 = cross_product(p2, p0, pixel) / area;
                    float w2 = cross_product(p0, p1, pixel) / area;

                    // interpolate color
                    GColor interpolatedColor = {
                            w0 * c0.r + w1 * c1.r + w2 * c2.r,
                            w0 * c0.g + w1 * c1.g + w2 * c2.g,
                            w0 * c0.b + w1 * c1.b + w2 * c2.b,
                            w0 * c0.a + w1 * c1.a + w2 * c2.a
                    };

                    // interpolate texture
                    GPoint interpolatedTex = {
                            w0 * t0.x + w1 * t1.x + w2 * t2.x,
                            w0 * t0.y + w1 * t1.y + w2 * t2.y
                    };

                    GPixel src = ColorToPixel(interpolatedColor);

                    if(!colors && texs){
                        GPixel ppx[1];
                        shader->shadeRow(interpolatedTex.x, interpolatedTex.y, 1, ppx);
                        GPixel *p = bitmap_.getAddr(x, y);
                        *p = blend(ppx[0], *p, blendMode);
                        continue;
                    }

                    if(shader == nullptr) {
                        if(x < 0 || x >= bitmap_.width() || y < 0 || y >= bitmap_.height() )
                            continue;

                        if(texs != nullptr)
                            drawRow(interpolatedTex.x, interpolatedTex.y, 1, src, paint, bitmap_);
                        else
                            drawRow(x, y, 1, src, paint, bitmap_);
                    }
                    else {
                        GPixel tnp[1];
                        if(texs) {
                            if (interpolatedTex.x < 0 || interpolatedTex.x >= bitmap_.width() ||
                                interpolatedTex.y < 0 || interpolatedTex.y >= bitmap_.height())
                                continue;
                        }

                        if(x < 0 || x >= bitmap_.width() || y < 0 || y >= bitmap_.height())
                            continue;

                        if(texs != nullptr)
                            shader->shadeRow(interpolatedTex.x, interpolatedTex.y, 1, tnp);
                        else
                            shader->shadeRow(x, y, 1, tnp);

                        GPixel *p = bitmap_.getAddr(x, y);

                        if(same_color) {
                            src = blend(tnp[0], src, blendMode);
                            *p = blend(src, *p, blendMode);
                        }
                        else {
                            if (texs) {
                                auto tnp_r = GPixel_GetR(tnp[0]) / 255.f;
                                auto tnp_g = GPixel_GetG(tnp[0]) / 255.f;
                                auto tnp_b = GPixel_GetB(tnp[0]) / 255.f;
                                auto tnp_a = GPixel_GetA(tnp[0]) / 255.f;
                                GColor finalColor = {
                                        interpolatedColor.r * tnp_r,
                                        interpolatedColor.g * tnp_g,
                                        interpolatedColor.b * tnp_b,
                                        interpolatedColor.a * tnp_a
                                };
//
                                GPixel finalPixel = ColorToPixel(finalColor);
                                *p = blend(finalPixel, *p, blendMode);
                            }
                            else{
                                *p = blend(src, *p, blendMode);
                            }
                        }
                    }
                }
            }
        }


    }
}

void MyGCanvas::drawQuad(const GPoint *verts, const GColor *colors, const GPoint *texs, int level, const GPaint &paint) {

//    return;

    if(level < 0)
        return;



//    return;
    std::vector<GPoint> subdivVerts;
    std::vector<GColor> subdivColors;
    std::vector<GPoint> subdivTexs;
    std::vector<int> indices;

    int gridSize = level + 1;
//    int vertexCount = (gridSize + 1) * (gridSize + 1);

    // Generate subdivided vertices, colors, and tex coords
    for (int y = 0; y <= gridSize; ++y) {
        float v = float(y) / gridSize;
        for (int x = 0; x <= gridSize; ++x) {
            float u = float(x) / gridSize;

            GPoint p = interpolateQuad(verts, u, v);
            subdivVerts.push_back(p);

            if (colors) {
                GColor c = interpolateQuad(colors, u, v);
                subdivColors.push_back(c);
            }

            if (texs) {
                GPoint t = interpolateQuad(texs, u, v);
                subdivTexs.push_back(t);
            }
        }
    }

    // Generate indices for triangles
    for (int y = 0; y < gridSize; ++y) {
        for (int x = 0; x < gridSize; ++x) {
            int topLeft = y * (gridSize + 1) + x;
            int topRight = topLeft + 1;
            int bottomLeft = topLeft + (gridSize + 1);
            int bottomRight = bottomLeft + 1;

            indices.push_back(topLeft);
            indices.push_back(bottomLeft);
            indices.push_back(topRight);

            indices.push_back(topRight);
            indices.push_back(bottomLeft);
            indices.push_back(bottomRight);
        }
    }

    // Draw the mesh
    drawMesh(subdivVerts.data(),
             colors ? subdivColors.data() : nullptr,
             texs ? subdivTexs.data() : nullptr,
//             nullptr,
             indices.size() / 3, indices.data(), paint);
}

//----------------------------------------------------------
bool MyGCanvas::pointInBound(const GPoint &p1) {
    int width = bitmap_.width();
    int height = bitmap_.height();
    if ((p1.x < 0) || (p1.x > width) ||
        (p1.y < 0) || (p1.y > height)) {
        return false;
    }
    return true;
}

bool MyGCanvas::edgeInBound(const GPoint &p1, const GPoint &p2,std::vector<GPoint>& crossPoint) {
    crossPoint.clear();
    int width = bitmap_.width();
    int height = bitmap_.height();
    if ((p1.x < 0 && p2.x < 0) || (p1.x > width && p2.x > width) ||
        (p1.y < 0 && p2.y < 0) || (p1.y > height && p2.y > height)) {
        return false;
    }
    std::vector<GPoint> bitmap_points;
    bitmap_points.push_back(GPoint{0.f * width, 0.f * height});
    bitmap_points.push_back(GPoint{1.f * width, 0.f * height});
    bitmap_points.push_back(GPoint{1.f * width, 1.f * height});
    bitmap_points.push_back(GPoint{0.f * width, 1.f * height});
    bitmap_points.push_back(GPoint{0.f * width, 0.f * height});
    for(int i = 0 ; i < 4; i++){
        auto q1 = bitmap_points[i];
        auto q2 = bitmap_points[i+1];
        // one IntersectionPoint
        auto rst = getIntersectionPoint(p1, p2, q1, q2);
        if(rst != nonstd::nullopt){
            crossPoint.push_back(rst.value());
        }
    }

    if(crossPoint.empty())
        return false;

    if(crossPoint.size() > 2){
        printf("Error crossPoint\n");
        return false;
    }

    return true;
}

GPoint MyGCanvas::letPointInBound(const GPoint &p) {
    int width = bitmap_.width();
    int height = bitmap_.height();
    GPoint adjustedPoint = p;
    adjustedPoint.x = std::max(0.0f, std::min(p.x, static_cast<float>(width)));
    adjustedPoint.y = std::max(0.0f, std::min(p.y, static_cast<float>(height)));
    return adjustedPoint;
}

nonstd::optional<GPoint> MyGCanvas::getIntersectionPoint(const GPoint &p0, const GPoint &p1, const GPoint &q0, const GPoint &q1)  {
    float A1 = p1.y - p0.y;
    float B1 = p0.x - p1.x;
    float C1 = A1 * p0.x + B1 * p0.y;

    float A2 = q1.y - q0.y;
    float B2 = q0.x - q1.x;
    float C2 = A2 * q0.x + B2 * q0.y;

    float det = A1 * B2 - A2 * B1;

    if (det == 0) {
        return nonstd::nullopt;
    } else {
        float x = (B2 * C1 - B1 * C2) / det;
        float y = (A1 * C2 - A2 * C1) / det;

        if (x < std::min(p0.x, p1.x) || x > std::max(p0.x, p1.x) ||
            x < std::min(q0.x, q1.x) || x > std::max(q0.x, q1.x) ||
            y < std::min(p0.y, p1.y) || y > std::max(p0.y, p1.y) ||
            y < std::min(q0.y, q1.y) || y > std::max(q0.y, q1.y)) {
            return nonstd::nullopt;
        }
        return GPoint{x, y};
    }
}

std::vector<std::pair<GPoint, GPoint>> MyGCanvas::createPolygon(const std::vector<std::pair<GPoint, GPoint>> &raw_polygon) {
    std::vector<std::pair<GPoint, GPoint>> new_polygon;
    for(const auto& edge : raw_polygon){
        int rst1 = pointInBound(edge.first) ? 1 : 0;
        int rst2 = pointInBound(edge.second) ? 1 : 0;
        int rst = rst1+2*rst2;
        switch (rst) {
            case 0: {
                /*std::pair<GPoint, GPoint> new_edge = edge;
                std::vector<GPoint> crossPoints;
                if(edgeInBound(edge.first, edge.second, crossPoints)){
                    if(crossPoints.size() == 2) {
                        new_edge.first = crossPoints[0];
                        new_edge.second = crossPoints[1];
                        GPoint v1{edge.second.x - edge.first.x,edge.second.y - edge.first.y };
                        GPoint v2{new_edge.second.x - new_edge.first.x,new_edge.second.y - new_edge.first.y };
                        float dir = v1.x * v2.x + v1.y * v2.y;
                        if(dir > 0){
                            new_polygon.push_back(new_edge);
                        }
                        else{
                            std::pair<GPoint, GPoint> new2_edge = new_edge;
                            new2_edge.second = new_edge.first;
                            new2_edge.first = new_edge.second;
                            new_polygon.push_back(new2_edge);
                        }
                    }
                    else if(crossPoints.size() == 0){
                        if((edge.first.x > 0 && edge.first.x < bitmap_.width()) &&
                           (edge.second.x > 0 && edge.second.x < bitmap_.width()) ||
                           (edge.first.y > 0 && edge.first.y < bitmap_.height()) &&
                           (edge.second.y > 0 && edge.second.y < bitmap_.height())){
                            break;
                        }
                        else{

                        }

                    }
                    else{
                        printf("Error\n");
                    }
                }*/
            }; break;
            case 1: {
                std::pair<GPoint, GPoint> new_edge = edge;
                std::vector<GPoint> crossPoints;
                if(edgeInBound(edge.first, edge.second, crossPoints)){
                    if(crossPoints.size() == 1) {
                        new_edge.second = crossPoints[0];
                        // --- second line
                        new_polygon.push_back(new_edge);
                        new_polygon.push_back(std::make_pair(new_edge.second,
                                                             letPointInBound(edge.second)));
                    }
                }
            }; break;
            case 2: {
                std::pair<GPoint, GPoint> new_edge = edge;
                std::vector<GPoint> crossPoints;
                if(edgeInBound(edge.first, edge.second, crossPoints)){
                    if(crossPoints.size() == 1) {
                        new_edge.first = crossPoints[0];
                        // --- first line
                        new_polygon.push_back(std::make_pair(letPointInBound(edge.first),
                                                             new_edge.first));
                        new_polygon.push_back(new_edge);
                    }
                }
            };break;
            case 3: {
                new_polygon.push_back(edge);
            };break;
        }
    }

    if(new_polygon.size() < 3){
        new_polygon.clear();
        return new_polygon;
    }

    std::vector<std::pair<GPoint, GPoint>> polygon;
    polygon.push_back(new_polygon.front());
    for(int i = 1; i < new_polygon.size(); i++){
        auto p0 = polygon.back();
        auto p1 = new_polygon[i];
        // [p0, p1] [p2, p3]    -> p1 = p2
        if(std::abs(p0.second.x - p1.first.x) + std::abs(p0.second.y - p1.first.y) < 1e-2){
            polygon.push_back(p1);
        }
        else{
            // printf("add\n");
            polygon.push_back(std::make_pair(p0.second, p1.first));
        }
    }
    return polygon;
}

std::vector<int> MyGCanvas::getIntersections(const std::vector<std::pair<GPoint, GPoint>> &edges, int y)  {
    std::vector<int> intersections;
    for (const auto& edge : edges) {
        GPoint p1 = edge.first;
        GPoint p2 = edge.second;
        if (p1.y > p2.y) std::swap(p1, p2);
        if (y < p1.y || y >= p2.y) continue;
        if (p1.y != p2.y) {
            int x = static_cast<int>(p1.x + (y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y));
            intersections.push_back(x);
        }
    }
    return intersections;
}

std::vector<GPoint> MyGCanvas::fillPolygon(const std::vector<std::pair<GPoint, GPoint>> &raw_polygon)  {
    std::vector<GPoint> filled_points;

    int minY = static_cast<int>(raw_polygon[0].first.y);
    int maxY = static_cast<int>(raw_polygon[0].first.y);
    for (const auto& edge : raw_polygon) {
        minY = std::min(GRoundToInt(minY), std::min(GRoundToInt(edge.first.y), GRoundToInt(edge.second.y)));
        maxY = std::max(GRoundToInt(maxY), std::max(GRoundToInt(edge.first.y), GRoundToInt(edge.second.y)));
    }

    for (int y = minY; y <= maxY; ++y) {
        std::vector<int> intersections = getIntersections(raw_polygon, y);

        std::sort(intersections.begin(), intersections.end());

        for (size_t i = 0; i < intersections.size(); i += 2) {
            if (i + 1 < intersections.size()) {
                int x1 = intersections[i];
                int x2 = intersections[i + 1];

                for (int x = x1; x <= x2; ++x) {
                    if(x < 0 || y < 0 || x >= bitmap_.width() || y >= bitmap_.height()) {
//                        std::cout << "small" << std::endl;
                        continue;
                    }
                    filled_points.emplace_back(GPoint{float(x), float(y)});
                }
            }
        }
    }

    return filled_points;
}

std::vector<GPoint> MyGCanvas::bresenhamLine(const GPoint &p1, const GPoint &p2) {
    std::vector<GPoint> line_points;
    int x0 = static_cast<int>(std::round(p1.x));
    int y0 = static_cast<int>(std::round(p1.y));
    int x1 = static_cast<int>(std::round(p2.x));
    int y1 = static_cast<int>(std::round(p2.y));

    int dx = std::abs(x1 - x0);
    int dy = std::abs(y1 - y0);
    int sx = (x0 < x1) ? 1 : -1;
    int sy = (y0 < y1) ? 1 : -1;
    int err = dx - dy;

    while (true) {
        line_points.push_back(GPoint{static_cast<float>(x0), static_cast<float>(y0)});

        if (x0 == x1 && y0 == y1) break;
        int e2 = 2 * err;
        if (e2 > -dy) {
            err -= dy;
            x0 += sx;
        }
        if (e2 < dx) {
            err += dx;
            y0 += sy;
        }
    }

    return line_points;
}

std::map<float, std::pair<float, float>> MyGCanvas::getLinePolygon(const std::vector<std::pair<GPoint, GPoint>> &raw_polygon) {
    std::vector<GPoint> all_points;

    for (const auto& edge : raw_polygon) {
        std::vector<GPoint> edge_points = bresenhamLine(edge.first, edge.second);
        all_points.insert(all_points.end(), edge_points.begin(), edge_points.end());
    }

    int minY = static_cast<int>(raw_polygon[0].first.y);
    int maxY = static_cast<int>(raw_polygon[0].first.y);
    for (const auto& edge : raw_polygon) {
        minY = std::min(GRoundToInt(minY), std::min(GRoundToInt(edge.first.y), GRoundToInt(edge.second.y)));
        maxY = std::max(GRoundToInt(maxY), std::max(GRoundToInt(edge.first.y), GRoundToInt(edge.second.y)));
    }

    std::map<float, std::pair<float, float>> key_map;
    for(int key = minY; key < maxY; key++) {
        std::pair<float, float> minmax_value = std::make_pair(1e6, -1);
        key_map.insert(std::make_pair(key, minmax_value));
    }

    for (const auto& edgePoint : all_points) {
        float key = GRoundToInt(edgePoint.y);
        float value = GRoundToInt(edgePoint.x);
        key_map[key].first = std::fmin(key_map[key].first, value);
        key_map[key].second = std::fmax(key_map[key].second, value);
    }
    for (auto& edge : key_map) {
        edge.second.second = edge.second.second - edge.second.first;
    }

    return key_map;
}

//----------------------------------------------------------


std::unique_ptr<GCanvas> GCreateCanvas(const GBitmap& bitmap){
    return std::unique_ptr<MyGCanvas>(new MyGCanvas(bitmap));
}

std::string GDrawSomething(GCanvas* canvas, GISize dim){
    GColor whiteColor{1.0f, 1.0f, 1.0f, 1.0f};
    GColor redColor{1.0f, 0.0f, 0.0f, 1.0f};

    // 1.White Canvas
    canvas->clear(whiteColor);

    GPoint TriPoints[3] = {
            {10, 10},
            {400, 100},
            {250, 400}
    };

    const GColor clr[] = {
            { 1, 0, 0, 1 }, { 0, 1, 0, 1 }, { 0, 0, 1, 1 },
    };
    const int indices[] = {
            0, 1, 2,
    };

    // 2.Red Paint
    GPaint paint;
    paint.setColor(redColor);
    paint.setBlendMode(GBlendMode::kSrc);


    GRect r = GRect::XYWH(10, 15, 100, 75);
    auto sh = GCreateLinearGradient({r.left, r.top}, {r.right, r.bottom},
                                    &redColor, 1);
    // paint.setShader(GCreateLinearGradient({1,0},{0,1}, &redColor, 1));
    paint.setShader(std::move(sh));
    if(paint.shareShader() == nullptr){
        printf("null\n");
    }

    // 3.Use Paint to Draw on a Canva
    canvas->scale(0.5, 0.5);
//     canvas->rotate(0.5);
    //    canvas->translate(0,-50);
    canvas->translate(250*0.25, 250*0.25);
    canvas->rotate(gFloatPI/2);
//    canvas->translate(250*0.25, 250*0.25);
//    canvas->translate(-250*0.25*0.5, -250*0.25*0.5);
    canvas->translate(-75, -450);

    //    canvas->drawConvexPolygon(HeartPoints.data(), HeartPoints.size(), paint);
//    canvas->drawConvexPolygon(TriPoints, 3, paint);

    canvas->drawMesh(TriPoints, clr, nullptr, 1, indices, GPaint());
//    GPathBuilder gPathBuilder;
//    gPathBuilder.addCircle(GPoint{400,400}, 100);
//    auto circle = gPathBuilder.detach();
//    canvas->drawPath(*circle, paint);
//
//    GPathBuilder bu;
//    bu.moveTo(10, 0);
//    bu.cubicTo({100, 100}, {100, -120}, {200, 0});
//    auto path = bu.detach();
//    canvas->drawPath(*path, paint);
//    int x;
//    std::cin >> x;
    canvas->save();
    return "Tri";
}