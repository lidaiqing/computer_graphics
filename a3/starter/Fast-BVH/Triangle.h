#ifndef Triangle_h_
#define Triangle_h_

#include <cmath>
#include "Object.h"

//! For the purposes of demonstrating the BVH, a simple sphere
struct Triangle : public Object {
  Vector3 center; // Center of the triangle
  Vector3 v1, v2, v3; // 3 vertices
  Vector3 vmin, vmax;

  Triangle(const Vector3& v1_, const Vector3& v2_, const Vector3& v3_)
    : v1(v1_), v2(v2_), v3(v3_) {
      float x = (v1.x + v2.x + v3.x) / 3.0;
      float y = (v1.y + v2.y + v3.y) / 3.0;
      float z = (v1.z + v2.z + v3.z) / 3.0;
      center = Vector3(x, y, z);

      float min_x = min(v1.x, min(v2.x, v3.x));
      float min_y = min(v1.y, min(v2.y, v3.y));
      float min_z = min(v1.z, min(v2.z, v3.z));
      float max_x = max(v1.x, max(v2.x, v3.x));
      float max_y = max(v1.y, max(v2.y, v3.y));
      float max_z = max(v1.z, max(v2.z, v3.z));

      vmin = Vector3(min_x, min_y, min_z);
      vmax = Vector3(max_x, max_y, max_z);
     }

  bool getIntersection(const Ray& ray, IntersectionInfo* I) const {
    //Möller–Trumbore intersection algorithm
    Vector3 e1 = v2 - v1;
    Vector3 e2 = v3 - v1;
    Vector3 P = ray.d ^ e2;
    float det = e1 * P;
    const float EPSILON = 0.000001;

    if (det > -EPSILON && det < EPSILON)
      return false;
    float inv_det = 1.f / det;
    Vector3 T = ray.o - v1;
    float u = T * P * inv_det;
    if (u < 0.f || u > 1.f) return false;

    Vector3 Q = T ^ e1;
    float v = ray.d * Q * inv_det;
    if (v < 0.f || u + v > 1.f) return false;

    float t = e2 * Q * inv_det;
    if (t <= EPSILON) return false;

    I->object = this;
    I->t = t;
    I->hit = ray.o + ray.d * t;
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

};

#endif
