#ifndef Plane_h_
#define Plane_h_

#include <cmath>
#include "Object.h"
#include <algorithm>

struct Plane : public Object {
  Vector3 center; // Center of the Plane
  Vector3 v1, v2, v3, v4; // 3 vertices
  Vector3 vmin, vmax;
  int index;
  Plane(const Vector3& v1_, const Vector3& v2_, const Vector3& v3_, const Vector3& v4_, int index_)
    : v1(v1_), v2(v2_), v3(v3_), v4(v4_), index(index_) {
      double x = (v1.x + v2.x + v3.x + v4.x) / 4.0;
      double y = (v1.y + v2.y + v3.y + v4.y) / 4.0;
      double z = (v1.z + v2.z + v3.z + v4.z) / 4.0;
      center = Vector3(x, y, z);

      vmin = ::min(v1, ::min(v2, ::min(v3, v4)));
      vmax = ::max(v1, ::max(v2, ::max(v3, v4)));
     }

  bool getIntersection(const Ray& ray, IntersectionInfo* I) const {
    //Möller–Trumbore intersection algorithm
    Vector3 e1 = v2 - v1;
    Vector3 e2 = v3 - v1;
    Vector3 P = ray.d ^ e2;
    double det = e1 * P;
    const double EPSILON = 0.000001;

    if (det > -EPSILON && det < EPSILON)
      return false;
    double inv_det = 1.f / det;
    Vector3 T = ray.o - v1;
    double u = T * P * inv_det;
    if (u < 0.f || u > 1.f) return false;

    Vector3 Q = T ^ e1;
    double v = ray.d * Q * inv_det;
    if (v < 0.f || v > 1.f) return false;

    double t = e2 * Q * inv_det;
    if (t <= EPSILON) return false;

    I->object = this;
    I->t = t;
    I->hit = ray.o + ray.d * t;

    // calculate texture mapping
    Vector3 hitP = I->hit;
    I->u = std::fabs(hitP * e1) / ::length(e1);
    I->v = std::fabs(hitP * e2) / ::length(e2);
    return true;
  }

  Vector3 getNormal(const IntersectionInfo& I) const {
    Vector3 e1 = v1 - v2;
    Vector3 e2 = v3 - v1;
    return normalize(e1 ^ e2);
  }

  BBox getBBox() const {
    return BBox(vmin, vmax);
  }

  Vector3 getCentroid() const {
    return center;
  }
  int getIndex() const {
    return index;
  }
};

#endif
