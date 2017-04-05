#ifndef Sphere_h_
#define Sphere_h_

#include <cmath>
#include "Object.h"

//! For the purposes of demonstrating the BVH, a simple sphere
struct Sphere : public Object {
  Vector3 center; // Center of the sphere
  float r, r2; // Radius, Radius^2
  int index;
  Sphere(const Vector3& center, float radius, int index_=0)
    : center(center), r(radius), r2(radius*radius), index(index_){ }

  bool getIntersection(const Ray& ray, IntersectionInfo* I) const {
    Vector3 s = center - ray.o;
    float sd = s * ray.d;
    float ss = s * s;

    // Compute discriminant
    float disc = sd*sd - ss + r2;

    // Complex values: No intersection
    if( disc < 0.f ) return false;

    // Assume we are not in a sphere... The first hit is the lesser valued
    I->object = this;
    I->t = sd - sqrt(disc);
    I->hit = ray.o + I->t * ray.d;
    return true;
  }

  Vector3 getNormal(const IntersectionInfo& I) const {
    return normalize(I.hit - center);
  }

  BBox getBBox() const {
    return BBox(center-Vector3(r,r,r), center+Vector3(r,r,r));
  }

  Vector3 getCentroid() const {
    return center;
  }
  int getIndex() const {
    return index;
  }

};

#endif
