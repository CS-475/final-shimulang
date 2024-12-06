#include <iostream>
#include "include/GPathBuilder.h"

void GPathBuilder::addRect(const GRect &rect, GPathDirection direction) {
    //! Polygon Exist
    if(rect.isEmpty())
        return;

    //! step1 : move the brush to First Point
    moveTo(rect.left, rect.top);

    //! step2 : draw Polygon
    // CW
    if (direction == GPathDirection::kCW) {
        lineTo(rect.right, rect.top);
        lineTo(rect.right, rect.bottom);
        lineTo(rect.left, rect.bottom);
        lineTo(rect.left, rect.top);
    }
    // CCW
    else {
        lineTo(rect.left, rect.bottom);
        lineTo(rect.right, rect.bottom);
        lineTo(rect.right, rect.top);
        lineTo(rect.left, rect.top);
    }
}

void GPathBuilder::addPolygon(const GPoint *pts, int count) {
    //! >=3 point draw a Polygon
    if(count < 3)
        return;
    //! step1 : move the brush to First Point
    moveTo(pts[0]);

    //! step2 : draw Polygon
    for(int i = 1; i < count; i++) {
        lineTo(pts[i]);
    }
}

void GPathBuilder::addCircle(GPoint center, float radius, GPathDirection direction) {
    //! Circle Exist
    if(radius == 0)
        return;

    std::vector<GPoint> pts;
    const int count = 32;
    float delta_angles = 2.f * (float)M_PI / (float)count;
    for (int i = 0; i < count; i++){
        float angles = float(i) * delta_angles;
        GPoint point = {center.x + radius * cos(angles), center.y + radius * sin(angles)};
        point.x = (int)point.x;
        point.y = (int)point.y;
        pts.push_back(point);
    }
    pts.push_back(pts.front());

    if (direction == GPathDirection::kCCW)
        std::reverse(std::begin(pts), std::end(pts));

    moveTo(pts[0]);
    for(int i = 1; i <= count; i+=2){
        quadTo(pts[i], pts[i+1]);
    }
//    for(auto p : fPts){
//        std::cout << p.x << ", " << p.y << std::endl;
//    }
//    std::cout << std::endl;
}
