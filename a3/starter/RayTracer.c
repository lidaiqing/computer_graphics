/*
  CSC418 - RayTracer code - Winter 2017 - Assignment 3&4

  Written Dec. 9 2010 - Jan 20, 2011 by F. J. Estrada
  Freely distributable for adacemic purposes only.

  Uses Tom F. El-Maraghi's code for computing inverse
  matrices. You will need to compile together with
  svdDynamic.c

  You need to understand the code provided in
  this file, the corresponding header file, and the
  utils.c and utils.h files. Do not worry about
  svdDynamic.c, we need it only to compute
  inverse matrices.

  You only need to modify or add code in sections
  clearly marked "TO DO"
*/

#include "utils.h"
#include "Fast-BVH/BVH.h"
#include "Fast-BVH/Triangle.h"
#include "Fast-BVH/Sphere.h"
using std::vector;
// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;
struct pointLS *light_list;
struct image *env_list[5];


std::unordered_map<int, object3D*> index_to_obj;
BVH *bvh;
vector<Object*> objects;
int MAX_DEPTH;


void buildScene(void)
{
 // Sets up all objects in the scene. This involves creating each object,
 // defining the transformations needed to shape and position it as
 // desired, specifying the reflectance properties (albedos and colours)
 // and setting up textures where needed.
 // Light sources must be defined, positioned, and their colour defined.
 // All objects must be inserted in the object_list. All light sources
 // must be inserted in the light_list.
 //
 // To create hierarchical objects:
 //   Copy the transform matrix from the parent node to the child, and
 //   apply any required transformations afterwards.
 //
 // NOTE: After setting up the transformations for each object, don't
 //       forget to set up the inverse transform matrix!

 struct object3D *o;
 struct pointLS *l;
 struct point3D p;


 ///////////////////////////////////////
 // TO DO: For Assignment 3 you have to use
 //        the simple scene provided
 //        here, but for Assignment 4 you
 //        *MUST* define your own scene.
 //        Part of your mark will depend
 //        on how nice a scene you
 //        create. Use the simple scene
 //        provided as a sample of how to
 //        define and position objects.
 ///////////////////////////////////////

 // Simple scene for Assignment 3:
 // Insert a couple of objects. A plane and two spheres
 // with some transformations.

 env_list[0]=readPPMimage(POS_X_PATH);
 env_list[1]=readPPMimage(NEG_X_PATH);
 env_list[2]=readPPMimage(POS_Y_PATH);
 env_list[3]=readPPMimage(NEG_Y_PATH);
 env_list[4]=readPPMimage(POS_Z_PATH);
 env_list[5]=readPPMimage(NEG_Z_PATH);



  /* Ground plane */
  //o=newPlane(.05,.35,.35,.05,.55,.8,.75,1,1,1);	// Note the plane is highly-reflective (rs=rg=.75) so we
//						// should see some reflections if all is done properly.
//						// Colour is close to cyan, and currently the plane is
//						// completely opaque (alpha=1). The refraction index is
//						// meaningless since alpha=1
// Scale(o,6,6,1.1);
// RotateZ(o,PI);
// Translate(o,0,-6,0);
// RotateX(o,PI/2.0);
// Translate(o,0,-6,0);
// //RotateX(o,PI/2.25);
// Translate(o,0,0,15);
// invert(&o->T[0][0],&o->Tinv[0][0]);		// Very important! compute
//						// and store the inverse
//						// transform for this object!
// insertObject(o,&object_list);			// Insert into object list
//
//  /* Up plane */
//  o=newPlane(.05,.75,.75,.05,.55,.8,.75,1,1,1);	// Note the plane is highly-reflective (rs=rg=.75) so we
//						// should see some reflections if all is done properly.
//						// Colour is close to cyan, and currently the plane is
//						// completely opaque (alpha=1). The refraction index is
//						// meaningless since alpha=1
// Scale(o,6,6,1.1);
// //RotateZ(o,PI);
// Translate(o,0,6,0);
// RotateX(o,-PI/2.0);
// Translate(o,0,6,0);
// //RotateX(o,PI/2.25);
// Translate(o,0,0,15);
// invert(&o->T[0][0],&o->Tinv[0][0]);		// Very important! compute
//						// and store the inverse
//						// transform for this object!
// insertObject(o,&object_list);			// Insert into object list
// loadTexture(o,"posy.PPM");
//
//

 // Let's add a couple spheres
 o=newSphere(.05,.95,.55,0.5,1,.25,.25,0.7,0.6,6);
 Scale(o,.75,.5,1.5);
 RotateY(o,PI/2);
 Translate(o,-1.45,1.1,3.5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 Vector3 pos(0, 2, 0);
 objects.push_back(new Sphere(pos, .5f, 0));
 index_to_obj[0] = o;

// o=newSphere(.05,.95,.55,0.05,0,0,0,0.6,1,6);
// Scale(o,.5,2.0,1.0);
// RotateZ(o,PI/1.5);
// Translate(o,1.75,1.25,5.0);
// invert(&o->T[0][0],&o->Tinv[0][0]);
// insertObject(o,&object_list);
// pos = Vector3(0, -2, 0);
// objects.push_back(new Sphere(pos, .5f, 1));
// index_to_obj[1] = o;

 vector<uint32_t> faces;
 vector<float> verts;
 read_ply_file(MESH_PATH, verts, faces);
 // Let's add a couple meshes
 int index = 1;
 Vector3 scale_factor(6,6,6);
 Vector3 translate(0,-1,-6);
 for (int i = 0; i < faces.size(); i += 3)
 {
    o = newTriangle(.05,.95,.55,.55,1,.25,.25,1,0.6,6);
    insertObject(o, &object_list);
    Vector3 v1(verts[faces[i]*3], verts[faces[i]*3+1], verts[faces[i]*3+2]);
    Vector3 v2(verts[faces[i+1]*3], verts[faces[i+1]*3+1], verts[faces[i+1]*3+2]);
    Vector3 v3(verts[faces[i+2]*3], verts[faces[i+2]*3+1], verts[faces[i+2]*3+2]);
    v1 = scale_factor.cmul(v1) + translate;
    v2 = scale_factor.cmul(v2) + translate;
    v3 = scale_factor.cmul(v3) + translate;
    objects.push_back(new Triangle(v1, v2, v3, index));
    index_to_obj[index] = o;
    index++;
 }

 // Insert a single point light source.
 p.px=-5;
 p.py=30;
 p.pz=1;
 p.pw=1;
 l=newPLS(&p,.65,.65,.65);
 insertPLS(l,&light_list);

 p.px=5;
 p.py=30;
 p.pz=1;
 p.pw=1;
 l=newPLS(&p,.65,.65,.65);
 insertPLS(l,&light_list);

 bvh = new BVH(&objects);

}
 void areaLighting(struct object3D* obj, struct pointLS* centre_light, struct point3D *p, struct point3D *n,struct ray3D *ray, int depth, double R, double G, double B, struct colourRGB* col, int sample_num )
 {

  struct point3D direction;
  // calculate p to light direction L
  direction.px = centre_light->p0.px, direction.py = centre_light->p0.py, direction.pz = centre_light->p0.pz, direction.pw = 1;
  subVectors(p, &direction);
  direction.pw = 0;
  /* reverse the direction (light to source direction) */
  direction.px=-direction.px;
  direction.py=-direction.py;
  direction.pz=-direction.pz;
  /* find the 'distance' between the light source and the point currently investigated */
  double distance=dot(&direction,&direction);
  direction.pw=0;
  /* Set the light source area to be 5% of the light source distance to the object */
  double areaBoundary = distance*0.002;
  double sampleBoundary = areaBoundary / (double)sample_num;
  double lambda = 0;
  double dummy_value;
  struct object3D *dummy_obj;
  struct point3D dummy_point;
  /* variable to calculate and sotre colours */
  colourRGB accumulated_colour;
  accumulated_colour.R = 0;
  accumulated_colour.G = 0;
  accumulated_colour.B = 0;
  /* iterate over an uniform surface */
  srand(time(NULL));
  int i = 0;
  int j = 0;
  for(i = 0; i < sample_num; i++)
  {
      for(j = 0; j < sample_num; j++)
      {

        double x_r=((float) rand() / (float)(RAND_MAX));
        if(i % 2 != 0)
        {
          x_r=-x_r;
        }
        double y_r=((float) rand() / (float)(RAND_MAX));
        if(j % 2 != 0)
        {
          y_r=-y_r;
        }
        /* get randomly generated x and y offset */
        double x_offset=sampleBoundary*x_r*0.1;
        double y_offset=sampleBoundary*y_r*0.1;

        double change_in_x=(-0.5*areaBoundary+(double)i*sampleBoundary+0.5*sampleBoundary+x_offset);
        double change_in_y=(0.5*areaBoundary-(double)j*sampleBoundary-0.5*sampleBoundary+y_offset);
        double change_in_z=(-(double)direction.px*change_in_x-(double)direction.py*change_in_y)/(-direction.pz);

        /* randomly offset the light source position */
        struct point3D origin;
        struct point3D changed_direction;
        origin.px = centre_light->p0.px + change_in_x;
        origin.py = centre_light->p0.py + change_in_y;
        origin.pz = centre_light->p0.pz + change_in_z;
        origin.pw = 1;

        changed_direction.px=origin.px - p->px;
        changed_direction.py=origin.py - p->py;
        changed_direction.pz=origin.pz - p->pz;
        changed_direction.pw=0;

        normalize(&changed_direction);
        ray3D *test_ray = newRay(p, &changed_direction);
      /* test if the altered ray can reach this point */
        findFirstHit_BVH(test_ray, false, &lambda, obj, &dummy_obj, &dummy_point, &dummy_point, &dummy_value, &dummy_value);
        //findFirstHit(test_ray, &lambda, obj, &dummy_obj, &dummy_point, &dummy_point, &dummy_value, &dummy_value);
        free(test_ray);

        if (lambda > 0) {
          // do not add contribute to color
        } else {
	          struct pointLS *sample_light=(struct pointLS *)malloc(sizeof(struct pointLS));
	          sample_light->col.R = centre_light->col.R;
	          sample_light->col.G = centre_light->col.G;
	          sample_light->col.B = centre_light->col.B;
	          sample_light->p0.px = origin.px;
	          sample_light->p0.py = origin.py;
	          sample_light->p0.pz = origin.pz;
	          sample_light->next = NULL;
              phongModel(obj, sample_light, p, n, ray, depth, R, G, B, &accumulated_colour);
	          free(sample_light);
        }
      }
  }

  col->R += (accumulated_colour.R / (sample_num * sample_num));
  col->G += (accumulated_colour.G / (sample_num * sample_num));
  col->B += (accumulated_colour.B / (sample_num * sample_num));

 };

void phongModel(struct object3D* obj, struct pointLS* light, struct point3D *p, struct point3D *n, struct ray3D *ray, int depth, double CR, double CG, double CB, struct colourRGB* col)
{
    struct point3D L;
    struct point3D* R;
    struct point3D N;
    struct point3D V;
    struct point3D neg_L;
    N.px = n->px, N.py = n->py, N.pz = n->pz, N.pw = 0;
    normalize(&N);

    // viewpoint V vector
    V.px = -ray->d.px, V.py = -ray->d.py, V.pz = -ray->d.pz, V.pw = 0;
    normalize(&V);

    // calculate light direction L
    L.px = p->px, L.py = p->py, L.pz = p->pz, L.pw = 1;
    subVectors(&light->p0, &L);
    // L is a direction
    L.pw = 0;
    normalize(&L);

    // calculate p0 to light
    neg_L.px = -L.px, neg_L.py = -L.py, neg_L.pz = -L.pz, neg_L.pw = 0;
    // calculate reflection direction
    R = getReflectionDirection(&L, p, n);

    // avoid redundant computation
    double c1 = max(0, dot(&N, &neg_L));
    //std::cout<<"dot1: " << c1 << std::endl;
    double c2 = pow(max(0, dot(R, &V)), obj->shinyness);

    //std::cout<<"dot2: " << c2 << std::endl;
    // multiply ambient and difuse terms by its color
    double cal_R = (obj->alb.ra * light->col.R + obj->alb.rd * c1 * light->col.R) * CR + obj->alb.rs * c2 * light->col.R;
    double cal_G = (obj->alb.ra * light->col.G + obj->alb.rd * c1 * light->col.G) * CG + obj->alb.rs * c2 * light->col.G;
    double cal_B = (obj->alb.ra * light->col.B + obj->alb.rd * c1 * light->col.B) * CB + obj->alb.rs * c2 * light->col.B;
    col->R += cal_R;
    col->G += cal_G;
    col->B += cal_B;
    free(R);
}
void rtShade(struct object3D *obj, struct point3D *p, struct point3D *n, struct ray3D *ray, int depth, double a, double b, struct colourRGB *col)
{
 // This function implements the shading model as described in lecture. It takes
 // - A pointer to the first object intersected by the ray (to get the colour properties)
 // - The coordinates of the intersection point (in world coordinates)
 // - The normal at the point
 // - The ray (needed to determine the reflection direction to use for the global component, as well as for
 //   the Phong specular component)
 // - The current racursion depth
 // - The (a,b) texture coordinates (meaningless unless texture is enabled)
 //
 // Returns:
 // - The colour for this ray (using the col pointer)
 //

 struct colourRGB tmp_col;	// Accumulator for colour components
 double R,G,B;			// Colour for the object in R G and B

 // This will hold the colour as we process all the components of
 // the Phong illumination model
 tmp_col.R=0;
 tmp_col.G=0;
 tmp_col.B=0;

 if (obj!= NULL && obj->texImg == NULL)		// Not textured, use object colour
 {
   R=obj->col.R;
   G=obj->col.G;
   B=obj->col.B;
   //std::cout<<"gotcha"<<std::endl;
 }
 else if (obj != NULL && obj->texImg != NULL)
 {
  // Get object colour from the texture given the texture coordinates (a,b), and the texturing function
  // for the object. Note that we will use textures also for Photon Mapping.
  obj->textureMap(obj->texImg,a,b,&R,&G,&B);
  //derict return --Only need texture//
  //printf("here\n");
 }

 //////////////////////////////////////////////////////////////
 // TO DO: Implement this function. Refer to the notes for
 // details about the shading model.
 //////////////////////////////////////////////////////////////

 // Be sure to update 'col' with the final colour computed here!
 // base case when depth > MAX_DEPTH
 if (depth > MAX_DEPTH) {
    return;
  }
 // Loop through each light source
 struct pointLS* lightPtr = light_list;
 while (lightPtr) {
    areaLighting(obj, lightPtr, p, n, ray, depth, R, G, B, &tmp_col, 1);
    lightPtr = lightPtr->next;
 }
     // reflection ray
    struct ray3D* reflectedRay = getReflectionRay(ray, p, n);
    struct colourRGB reflectedCol;
    reflectedCol.R = reflectedCol.G = reflectedCol.B = 0;
    rayTrace(reflectedRay, depth + 1, &reflectedCol, obj);
    free(reflectedRay);

    reflectedCol.R *= obj->alb.rg;
    reflectedCol.G *= obj->alb.rg;
    reflectedCol.B *= obj->alb.rg;


    // refraction ray
   /*struct ray3D* refractedRay = getRefractionRay(ray, obj, p, n);
   struct colourRGB refractedCol;
   refractedCol.R = refractedCol.G = refractedCol.B = 0;
   if (obj->alpha < 1) rayTrace(refractedRay, depth + 1, &refractedCol, obj);

   free(refractedRay);*/
   //refractedCol.R *= obj->r_index;
   //refractedCol.G *= obj->r_index;
   //refractedCol.B *= obj->r_index;

    col->R += (tmp_col.R + reflectedCol.R);// + refractedCol.R);
    col->G += (tmp_col.G + reflectedCol.G);// + refractedCol.G);
    col->B += (tmp_col.B + reflectedCol.B);// + refractedCol.B);
 return;

}
void isBlock(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p, struct point3D *n, double *a, double *b)
{
 // Find the first intersection between the ray and any objects in the scene.
 // It returns:
 //   - The lambda at the intersection (or < 0 if no intersection)
 //   - The pointer to the object at the intersection (so we can evaluate the colour in the shading function)
 //   - The location of the intersection point (in p)
 //   - The normal at the intersection point (in n)
 //
 // Os is the 'source' object for the ray we are processing, can be NULL, and is used to ensure we don't
 // return a self-intersection due to numerical errors for recursive raytrace calls.
 //
 /* Set lambda to equal to -1 to indicate the ray does not intersect any object */
  *lambda = -1;
  *obj = NULL;
 /* set to -1 so the default is invalid */
  struct object3D *objPtr = object_list;
  double cur_lambda;
  while (objPtr) {
    if (objPtr == Os) {
      objPtr = objPtr->next;
      continue;
    }
    cur_lambda = -1;
    objPtr->intersect(objPtr, ray, &cur_lambda, p, n, a, b);
    if (cur_lambda > 0) {
      *obj = objPtr;
      *lambda = cur_lambda;
      return;
    }
    objPtr = objPtr->next;
  }

}
void findFirstHit_BVH(struct ray3D *ray, bool occlusion, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p, struct point3D *n, double *a, double *b)
{
  //wrap BVH method
  Vector3 ray_o(ray->p0.px, ray->p0.py, ray->p0.pz);
  Vector3 ray_d(ray->d.px, ray->d.py, ray->d.pz);
  Ray rayBVH(ray_o, normalize(ray_d));
  //std::cout<< "x: " << rayBVH.d.x << " y: " << rayBVH.d.y << " z:" << rayBVH.d.z << std::endl;
  IntersectionInfo I;
  bool hit = bvh->getIntersection(rayBVH, &I, occlusion);
  if (!hit) {
    *lambda = -1;
    *obj = NULL;
  } else {
    *lambda = I.t;
    Vector3 normal = I.object->getNormal(I);
    n->px = normal.x, n->py = normal.y, n->pz = normal.z, n->pw = 0;
    p->px = I.hit.x, p->py = I.hit.y, p->pz = I.hit.z, p->pw = 1;
    *obj = index_to_obj[I.object->getIndex()];
    // TODO add a, b
  }
}
void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p, struct point3D *n, double *a, double *b)
{
 // Find the closest intersection between the ray and any objects in the scene.
 // It returns:
 //   - The lambda at the intersection (or < 0 if no intersection)
 //   - The pointer to the object at the intersection (so we can evaluate the colour in the shading function)
 //   - The location of the intersection point (in p)
 //   - The normal at the intersection point (in n)
 //
 // Os is the 'source' object for the ray we are processing, can be NULL, and is used to ensure we don't
 // return a self-intersection due to numerical errors for recursive raytrace calls.
 //
 /////////////////////////////////////////////////////////////
 // TO DO: Implement this function. See the notes for
 // reference of what to do in here
 /////////////////////////////////////////////////////////////
 // Inserts an object into the object list.
 /* Set lambda to equal to -1 to indicate the ray does not intersect any object */
 *lambda = -1;
 *obj = NULL;
 /* set to -1 so the default is invalid */
 struct object3D *objPtr = object_list;
 struct object3D *best_obj = NULL;
 double min_lambda = 10000;
 double cur_lambda;
 struct point3D cur_p, cur_n;
 double cur_a, cur_b;
 while (objPtr) {
   if (objPtr == Os) {
     objPtr = objPtr->next;
     continue;
   }
   cur_lambda = -1;
   objPtr->intersect(objPtr, ray, &cur_lambda, &cur_p, &cur_n, &cur_a, &cur_b);
   if (cur_lambda > 0 && cur_lambda < min_lambda) {
     min_lambda = cur_lambda;
     best_obj = objPtr;
     *p = cur_p;
     *n = cur_n;
     *a = cur_a;
     *b = cur_b;
   }
   objPtr = objPtr->next;
 }
if (min_lambda != 10000) {
  *obj = best_obj;
  *lambda = min_lambda;
}

}

void rayTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os)
{
 // Ray-Tracing function. It finds the closest intersection between
 // the ray and any scene objects, calls the shading function to
 // determine the colour at this intersection, and returns the
 // colour.
 //
 // Os is needed for recursive calls to ensure that findFirstHit will
 // not simply return a self-intersection due to numerical
 // errors. For the top level call, Os should be NULL. And thereafter
 // it will correspond to the object from which the recursive
 // ray originates.
 //

 double lambda;		// Lambda at intersection
 double a,b;		// Texture coordinates
 struct object3D *obj;	// Pointer to object at intersection
 struct point3D p;	// Intersection point
 struct point3D n;	// Normal at intersection
 struct colourRGB I;	// Colour returned by shading function

 ///////////////////////////////////////////////////////
 // TO DO: Complete this function. Refer to the notes
 // if you are unsure what to do here.
 ///////////////////////////////////////////////////////
 /* obj is null because it is the first recursion so not from any object */
 /* By the end of this function call, obj will point to the object this ray firstly intersects */
 //findFirstHit_BVH(ray, false, &lambda, Os, &obj, &p, &n, &a, &b);
 findFirstHit_BVH(ray, false, &lambda, Os, &obj, &p, &n, &a, &b);
    if(lambda > 0)
    {
      rtShade(obj, &p, &n, ray, depth, a, b, col);
    }
    else
    {
      if (depth == 1) {
        double R, G, B;
        int index;
        convert_xyz_to_cube_uv(ray->p0.px + 2 * ray->d.px, ray->p0.py + 2 * ray->d.py, ray->p0.pz + 2 * ray->d.pz, &index, &a, &b);
        texMap(env_list[index], a, b, &R, &G, &B);
        col->R = R, col->G = G, col->B = B;
        return;
      } else {
        double R, G, B;
        int index;
        convert_xyz_to_cube_uv(ray->p0.px + 2 * ray->d.px, ray->p0.py + 2 * ray->d.py, ray->p0.pz + 2 * ray->d.pz, &index, &a, &b);
        texMap(env_list[index], a, b, &R, &G, &B);
        col->R += R, col->G += G, col->B += B;
      }
    }

}

 void add_antialiasing(point3D eye,double x,double y,double z,int multiplier, double pixel_boundary, colourRGB *col )
 {
  srand(time(NULL));
  colourRGB accumulated_colour;
  accumulated_colour.R=0;
  accumulated_colour.G=0;
  accumulated_colour.B=0;
  colourRGB colour;
  colour.R=0;
  colour.G=0;
  colour.B=0;
  double boundary=pixel_boundary/(double)multiplier;
  int i=0;
  int j=0;
  for(i=0;i<multiplier;i++)
  {
      for(j=0;j<multiplier;j++)
      {

      double x_r=((float) rand() / (float)(RAND_MAX));
      if(i%2!=0)
      {
      x_r=-x_r;
      }
      double y_r=((float) rand() / (float)(RAND_MAX));
      if(j%2!=0)
      {
      y_r=-y_r;
      }

      double x_offset=boundary*x_r*0.5;
      double y_offset=boundary*y_r*0.5;



      point3D direction;
      direction.px=x-0.5*pixel_boundary+(double)i*boundary+0.5*boundary+x_offset-eye.px;
      direction.py=y+0.5*pixel_boundary-(double)j*boundary-0.5*boundary-y_offset-eye.py;
      direction.pz=z;
      direction.pw=0;
      //normalize(&direction);
      ray3D *ray=newRay(&eye, &direction);

      rayTrace(ray,1,&colour,NULL);
      accumulated_colour.R+=colour.R;
      accumulated_colour.G+=colour.G;
      accumulated_colour.B+=colour.B;
      free(ray);
      colour.R=0;
      colour.G=0;
      colour.B=0;

      }
  }

  col->R += (accumulated_colour.R / (multiplier*multiplier));
  col->G += (accumulated_colour.G / (multiplier*multiplier));
  col->B += (accumulated_colour.B / (multiplier*multiplier));

 };


int main(int argc, char *argv[])
{
 // Main function for the raytracer. Parses input parameters,
 // sets up the initial blank image, and calls the functions
 // that set up the scene and do the raytracing.
 struct image *im;	// Will hold the raytraced image
 struct view *cam;	// Camera and view for this scene
 int sx;		// Size of the raytraced image
 int antialiasing;	// Flag to determine whether antialiaing is enabled or disabled
 char output_name[1024];	// Name of the output file for the raytraced .ppm image
 struct point3D e;		// Camera view parameters 'e', 'g', and 'up'
 struct point3D g;
 struct point3D up;
 double du, dv;			// Increase along u and v directions for pixel coordinates
 struct point3D pc,d;		// Point structures to keep the coordinates of a pixel and
				// the direction or a ray
 struct ray3D *ray;		// Structure to keep the ray from e to a pixel
 struct colourRGB col;		// Return colour for raytraced pixels
 struct colourRGB background;   // Background colour
 int i,j;			// Counters for pixel coordinates
 unsigned char *rgbIm;

 if (argc<5)
 {
  fprintf(stderr,"RayTracer: Can not parse input parameters\n");
  fprintf(stderr,"USAGE: RayTracer size rec_depth antialias output_name\n");
  fprintf(stderr,"   size = Image size (both along x and y)\n");
  fprintf(stderr,"   rec_depth = Recursion depth\n");
  fprintf(stderr,"   antialias = A single digit, 0 disables antialiasing. Anything else enables antialiasing\n");
  fprintf(stderr,"   output_name = Name of the output file, e.g. MyRender.ppm\n");
  exit(0);
 }
 sx=atoi(argv[1]);
 MAX_DEPTH=atoi(argv[2]);
 if (atoi(argv[3])==0) antialiasing=0; else antialiasing=1;
 strcpy(&output_name[0],argv[4]);

 fprintf(stderr,"Rendering image at %d x %d\n",sx,sx);
 fprintf(stderr,"Recursion depth = %d\n",MAX_DEPTH);
 if (!antialiasing) fprintf(stderr,"Antialising is off\n");
 else fprintf(stderr,"Antialising is on\n");
 fprintf(stderr,"Output file name: %s\n",output_name);

 object_list=NULL;
 light_list=NULL;

 // Allocate memory for the new image
 im=newImage(sx, sx);
 if (!im)
 {
  fprintf(stderr,"Unable to allocate memory for raytraced image\n");
  exit(0);
 }
 else rgbIm=(unsigned char *)im->rgbdata;

 ///////////////////////////////////////////////////
 // TO DO: You will need to implement several of the
 //        functions below. For Assignment 3, you can use
 //        the simple scene already provided. But
 //        for Assignment 4 you need to create your own
 //        *interesting* scene.
 ///////////////////////////////////////////////////
 buildScene();		// Create a scene. This defines all the
		        	// objects in the world of the raytracer

 //////////////////////////////////////////
 // TO DO: For Assignment 3 you can use the setup
 //        already provided here. For Assignment 4
 //        you may want to move the camera
 //        and change the view parameters
 //        to suit your scene.
 //////////////////////////////////////////

 // Mind the homogeneous coordinate w of all vectors below. DO NOT
 // forget to set it to 1, or you'll get junk out of the
 // geometric transformations later on.

 // Camera center is at (0,0,-1)
 e.px=0;
 e.py=0;
 e.pz=-8;
 e.pw=1;

 // To define the gaze vector, we choose a point 'pc' in the scene that
 // the camera is looking at, and do the vector subtraction pc-e.
 // Here we set up the camera to be looking at the origin, so g=(0,0,0)-(0,0,-1)
 g.px=0;
 g.py=0;
 g.pz=1;
 g.pw=0;

 // Define the 'up' vector to be the Y axis
 up.px=0;
 up.py=1;
 up.pz=0;
 up.pw=0;
//
// struct point3D z_direction;
// z_direction.px=0;
// z_direction.py=0;
// z_direction.pz=1;
// z_direction.pw=0;
//
// struct point3D x_direction;
// x_direction.px=1;
// x_direction.py=0;
// x_direction.pz=0;
// x_direction.pw=0;
//
// struct point3D z_direction;
// y_direction.px=0;
// y_direction.py=1;
// y_direction.pz=0;
// y_direction.pw=0;
//
//
// object3D *camera_transform=newPlane(.05,.35,.35,.05,.55,.8,.75,1,1,1);
// Translate(o,-e.px,-e.py,-e.pz);
// double cos=((g.pz)*(g.pz))/((g.pz)*(sqrt(g.px*g.px+g.py*g.py+g.py*g.py)));
// double y_rotation=acos(cos);
//
// dot(&z_direction,&g);
// RotateZ(o,PI);
// RotateX(o,PI/2.0);
// Translate(o,0,-6,0);
 // Set up view with given the above vectors, a 4x4 window,
 // and a focal length of -1 (why? where is the image plane?)
 // Note that the top-left corner of the window is at (-2, 2)
 // in camera coordinates.
 cam=setupView(&e, &g, &up, -3, -2, 2, 4);

 if (cam==NULL)
 {
  fprintf(stderr,"Unable to set up the view and camera parameters. Our of memory!\n");
  cleanup(object_list,light_list);
  deleteImage(im);
  exit(0);
 }

 // Set up background colour here
 background.R=0;
 background.G=0;
 background.B=0;

 // Do the raytracing
 //////////////////////////////////////////////////////
 // TO DO: You will need code here to do the raytracing
 //        for each pixel in the image. Refer to the
 //        lecture notes, in particular, to the
 //        raytracing pseudocode, for details on what
 //        to do here. Make sure you undersand the
 //        overall procedure of raytracing for a single
 //        pixel.
 //////////////////////////////////////////////////////
 du=cam->wsize/(sx-1);		// du and dv. In the notes in terms of wl and wr, wt and wb,
 dv=-cam->wsize/(sx-1);		// here we use wl, wt, and wsize. du=dv since the image is
				// and dv is negative since y increases downward in pixel
				// coordinates and upward in camera coordinates.

 fprintf(stderr,"View parameters:\n");
 fprintf(stderr,"Left=%f, Top=%f, Width=%f, f=%f\n",cam->wl,cam->wt,cam->wsize,cam->f);
 fprintf(stderr,"Camera to world conversion matrix (make sure it makes sense!):\n");
 printmatrix(cam->C2W);
 fprintf(stderr,"World to camera conversion matrix\n");
 printmatrix(cam->W2C);
 fprintf(stderr,"\n");





 fprintf(stderr,"Rendering row: ");
 #pragma omp parallel for private(i)
 for (j=0;j<sx;j++)		// For each of the pixels in the image
 {
  //fprintf(stderr,"%d/%d, ",j,sx);
  for (i=0;i<sx;i++)
  {
    struct point3D ray_direction;
    ray_direction.px=(-cam->wsize/2)+i*(du)+0.5*(du)-cam->e.px;
    ray_direction.py=(cam->wsize/2)+j*(dv)+0.5*(dv)-cam->e.py;
    ray_direction.pz=(-cam->f);
    ray_direction.pw=0;
    struct ray3D* ray_thread = newRay(&cam->e, &ray_direction);
    ///////////////////////////////////////////////////////////////////
    // TO DO - complete the code that should be in this loop to do the
    //         raytracing!
    ///////////////////////////////////////////////////////////////////
    struct colourRGB col_thread;
    col_thread.R = col_thread.G = col_thread.B = 0;

    if (antialiasing) {
      add_antialiasing(cam->e,(-cam->wsize/2)+i*(du)+0.5*(du),(cam->wsize/2)+j*(dv)+0.5*(dv),(-cam->f), 3, du, &col_thread);
    }
    else rayTrace(ray_thread, 1, &col_thread, NULL);
    *(rgbIm + 3 * (j * sx  + i)) = col_thread.R * 255 > 255 ? 255 : col_thread.R * 255;
    *(rgbIm + 3 * (j * sx  + i) + 1) = col_thread.G * 255 > 255 ? 255 : col_thread.G * 255;
    *(rgbIm + 3 * (j * sx  + i) + 2) = col_thread.B * 255 > 255 ? 255 : col_thread.B * 255;
    free(ray_thread);

  } // end for i
 } // end for j

 fprintf(stderr,"\nDone!\n");

 // Output rendered image
 imageOutput(im,output_name);

 // Exit section. Clean up and return.
 cleanup(object_list, light_list);		// Object and light lists
 deleteImage(im);				// Rendered image
 free(cam);					// camera view
 exit(0);
}
