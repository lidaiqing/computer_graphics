#ifndef Triangle_h_
#define Triangle_h_

#include <cmath>
#include "Object.h"
#include <algorithm>

struct Triangle : public Object {
  Vector3 center; // Center of the triangle
  Vector3 v1, v2, v3; // 3 vertices
  Vector3 vmin, vmax;
  int index;
  Triangle(const Vector3& v1_, const Vector3& v2_, const Vector3& v3_, int index_)
    : v1(v1_), v2(v2_), v3(v3_), index(index_) {
      float x = (v1.x + v2.x + v3.x) / 3.0;
      float y = (v1.y + v2.y + v3.y) / 3.0;
      float z = (v1.z + v2.z + v3.z) / 3.0;
      center = Vector3(x, y, z);

      vmin = ::min(v1, ::min(v2, v3));
      vmax = ::max(v1, ::max(v2, v3));
     }

  bool getIntersection(const Ray& ray, IntersectionInfo* I) const {
    //Möller–Trumbore intersection algorithm
    Vector3 e1 = v2 - v1;
    Vector3 e2 = v3 - v1;
    Vector3 D = normalize(ray.d);
    Vector3 P = D ^ e2;
    float det = e1 * P;
    float ray_d_len = length(ray.d);
    const float EPSILON = 0.000001;

    if (det > -EPSILON && det < EPSILON)
      return false;
    float inv_det = 1.f / det;
    Vector3 T = ray.o - v1;
    float u = (T * P) * inv_det;
    if (u < 0.f || u > 1.f) return false;

    Vector3 Q = T ^ e1;
    float v = D * Q * inv_det;
    if (v < 0.f || u + v > 1.f) return false;

    float t = e2 * Q * inv_det;
    if (t <= EPSILON) return false;

    I->object = this;
    // t is for original length
    I->t = t / ray_d_len;
    I->hit = ray.o + D * t;

    // calculate texture mapping
    Vector3 hitP = I->hit;
    double BaryV1 = ((v2.y - v3.y) * (hitP.x - v3.x) + (v3.x - v2.x) * (hitP.y - v3.y)) / ((v2.y - v3.y) * (v1.x - v3.x) + (v3.x - v2.x) * (v1.y - v3.y));
    double BaryV2 = ((v3.y - v1.y) * (hitP.x - v3.x) + (v1.x - v3.x) * (hitP.y - v3.y)) / ((v2.y - v3.y) * (v1.x - v3.x) + (v3.x - v2.x) * (v1.y - v3.y));
    double BaryV3 = 1.0 - BaryV1 - BaryV2;
    I->u = fabs(BaryV3) > 1 ? 0,5 : fabs(BaryV3);
    I->v = fabs(BaryV2 + BaryV3) > 1 ? 0.5 : fabs(BaryV2 + BaryV3);  
    return true;
  }

  Vector3 getNormal(const IntersectionInfo& I) const {
    Vector3 e1 = v2 - v1;
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
