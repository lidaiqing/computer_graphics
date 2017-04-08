#ifndef Sphere_h_
#define Sphere_h_

#include <cmath>
#include "Object.h"
#include "../RayTracer.h"
#include "../utils.h"

//! For the purposes of demonstrating the BVH, a simple sphere
struct Sphere : public Object {
  Vector3 center; // Center of the sphere
  float r, r2; // Radius, Radius^2
  int index;
  struct object3D* old_obj;
  Vector3 vmin, vmax;
  Sphere(const Vector3& center, float radius, struct object3D* old_obj_, int index_)
    : center(center), r(radius), r2(radius*radius), old_obj(old_obj_), index(index_){
      vmin = center - Vector3(r,r,r);
      vmax = center + Vector3(r,r,r);
    }

  bool getIntersection(const Ray& ray, IntersectionInfo* I) const {
    struct point3D ray_o;
    ray_o.px = ray.o.x, ray_o.py = ray.o.y, ray_o.pz = ray.o.z, ray_o.pw = 1;
    struct point3D ray_d;
    ray_d.px = ray.d.x, ray_d.py = ray.d.y, ray_d.pz = ray.d.z, ray_d.pw = 0;
    struct ray3D* old_ray = newRay(&ray_o, &ray_d);
    double lambda;
    struct point3D p;
    struct point3D n;
    double a, b;
    old_obj->intersect(old_obj, old_ray, &lambda, &p, &n, &a, &b);
    free(old_ray);
    if (lambda < 0) return false;
//    Vector3 s = center - ray.o;
//    Vector3 D(ray.d.x, ray.d.y, ray.d.z);
//    normalize(D);
//    float sd = s * D;
//    float ss = s * s;
//
//    // Compute discriminant
//    float disc = sd*sd - ss + r2;
//
//    // Complex values: No intersection
//    if( disc < 0.f ) return false;
//
//    // Assume we are not in a sphere... The first hit is the lesser valued
//    I->object = this;
//    I->t = (sd - sqrt(disc));
//    I->hit = ray.o + I->t * D;
//    I->normal = normalize(I->hit - center);
//    return true;
    I->object = this;
    I->t = lambda;
    Vector3 hitP(p.px, p.py, p.pz);
    I->hit = hitP;
    I->u = a;
    I->v = b;
    //std::cout << "u: " << I->u << "v: " << I->v << std::endl;
    Vector3 normal((float)n.px, (float)n.py, (float)n.pz);
    I->normal = normal;
    return true;
  }

  Vector3 getNormal(const IntersectionInfo& I) const {
     return I.normal;
  }

  BBox getBBox() const {
    return BBox(vmin, vmax);
    //return BBox(center-Vector3(r,r,r), center+Vector3(r,r,r));
  }

  Vector3 getCentroid() const {
    return center;
  }
  int getIndex() const {
    return index;
  }

};

#endif
