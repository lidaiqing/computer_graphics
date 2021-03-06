#ifndef Plane_h_
#define Plane_h_

#include <cmath>
#include "Object.h"
#include "../RayTracer.h"
#include "../utils.h"
#include <algorithm>

struct Plane : public Object {
  Vector3 center; // Center of the Plane
  Vector3 v1, v2, v3, v4; // 3 vertices
  Vector3 vmin, vmax;
  struct object3D* old_obj;
  int index;
  Plane(const Vector3& v1_, const Vector3& v2_, const Vector3& v3_, const Vector3& v4_, struct object3D* old_obj_, int index_)
    : v1(v1_), v2(v2_), v3(v3_), v4(v4_), old_obj(old_obj_), index(index_) {
      double x = (v1.x + v2.x + v3.x + v4.x) / 4.0;
      double y = (v1.y + v2.y + v3.y + v4.y) / 4.0;
      double z = (v1.z + v2.z + v3.z + v4.z) / 4.0;
      center = Vector3(x, y, z);

      vmin = ::min(v1, ::min(v2, ::min(v3, v4)));
      vmax = ::max(v1, ::max(v2, ::max(v3, v4)));
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
    I->object = this;
    I->t = lambda;
    Vector3 hitP(p.px, p.py, p.pz);
    I->hit = hitP;
    I->u = fabs(a) > 1 ? 0.5 : fabs(a);
    I->v = fabs(b) > 1 ? 0.5 : fabs(b);
    Vector3 normal(n.px, n.py, n.pz);
    I->normal = normal;
    return true;
  }

  Vector3 getNormal(const IntersectionInfo& I) const {
    return I.normal;
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
