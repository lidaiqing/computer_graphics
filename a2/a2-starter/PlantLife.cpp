/***********************************************************
           CSC418/2504, Winter 2017, St. George
                 PlantLife.cpp

  Learning Goals:

   After completing this assignment you should be able to

   * Design and implement objects as a hierarchy or parts
     related by transformations
   * Do simple animation of hierarchical objects
   * Understand how plant structures can be generated by
     using a small set of rules working on strings of
     symbols
   * Understand and use illumination in OpenGL
   * Use simple functions to generate a ground surface map.

  What to do:

   * Read the handout! it contains important details about
     the work you are supposed to do here and explains
     how L-systems (used to generate the plants) work.
   * Read the comments in this starter file, which detail
     what the existing code does and indicate what needs
     doing.
   * Complete all the parts marked // TO DO:
   * Add any // CRUNCHY:  extensions you want
   * Test thoroughly and make sure it works on CDF
   * Complete the CHECKLIST and REPORT files!

Program Code V3.0: F. Estrada, Sep 2012.
***********************************************************/

// OpenGL libraries
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <imgui.h>
#include "imgui_impl_glut.h"
#include <iostream>

// Standard UNIX/C libraries
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#define MAX_PLANTS 25		// Maximum number of plants for plant forests
#define GRID_RESOLVE 64		// Size of the surface grid
/******************************************************************************
  Data structures section
*******************************************************************************/
// Tree structure used to hold the plant model for rendering (remember the
// hierarchical rendering we discussed in tutorial, and you used for A1!)
struct PlantNode{
  char type;         // Node type, a=stem, b=stem, c=leaf, d=flower    - You may add your own types as needed
  GLfloat z_ang;     // Rotation of this node w.r.t. parent's z-axis   - Rotates around the current stem direction
  GLfloat x_ang;     // Rotation of this node w.r.t. parent's x-axis   - Rotates away from current stem direction
  GLfloat scl;       // Local scale for this node

  GLfloat f_c_R;     // Colour for this node's component - You can use it if needed, for example, to give
  GLfloat f_c_G;     // flowers different colours.
  GLfloat f_c_B;
  GLint misc_A;	     // Misc variables A to C are UNDEFINED. You can use them any way you want to create
  GLint misc_B;      // different types of structures for the stems, leafs, and flowers.
  GLint misc_C;

  struct PlantNode *left;     // Left child for this node (NULL means the node is a terminal)
  struct PlantNode *right;    // Right child for this node (NULL means the node is a terminal)

  // QUESTION: We don't have a translation component... WHY?

  // NOTE: We don't have pointers to the parents. Make sure your drawing function knows what
  // it's doing!
};

/******************************************************************************
  Global data
******************************************************************************/
const float PI = 3.14159;
struct PlantNode *PlantForest[MAX_PLANTS];	// Array of pointers for a plant forest, use when drawing multiple plants
GLfloat ForestXYZ[MAX_PLANTS][3];		// Location of plants in the forest as (x, y, z)
GLfloat GroundXYZ[GRID_RESOLVE][GRID_RESOLVE][3];	// Array to store ground surface points as (x,y,z) points
GLfloat GroundNormals[GRID_RESOLVE][GRID_RESOLVE][3];   // Array to hold normal vectors at vertices (will use normal to triangle
                                                        // defined by the vertices at [i][j], [i+1][j], and [i][j+1]
ImVec4 clear_color = ImColor(0.1f,0.1f,0.1f,1.0f);

// Texture data
int textures_on;				// Flag to indicate texturing should be enabled for leafs/flowers
GLuint l_texture,p_texture;			// Identifiers for OpenGL texture data for leaf and petal
int l_sx,l_sy,p_sx,p_sy;
unsigned char *leaf_texture;			// Pointer to leaf texture data (RGBA)
unsigned char *petal_texture;			// Pointer to petal texture data (RGBA)

// Window settings and GLUI variables
int windowID;               // Glut window ID (for display)
int Win[2];                 // window (x,y) size

// GLUI interface variables
GLfloat global_Z;	    // User controlled global rotation around Z
GLfloat global_scale;       // User controlled global scale factor

// Command-line parameters - These affect the shape of the plants, and control the number of plants in the forest
GLfloat Z_angle;      // Max angle around local Z-axis for each level
GLfloat X_angle;      // Max angle around local X-axis at each branching point
GLfloat scale_mult;   // Scale multiplier for children. Children nodes have scale=(scale_mult*parent_scale);
GLint n_levels;       // Number of levels of branching in the plant
GLint n_plants;	      // Number of plants in a plant forest in [1,MAX_PLANTS]

// Transition probabilities for the L-system specification
GLfloat Paab;
GLfloat Paac;
GLfloat Paad;
GLfloat Pacd;
GLfloat Pba;
GLfloat Pbc;
GLfloat Pbd;

/******************************************************************************
  Function Prototypes
*******************************************************************************/
// Initialization functions
void initGlut(char* winName);
unsigned char *readPPM(const char *name, int *sx, int *sy);

// Callbacks for handling events in glut
void WindowReshape(int w, int h);
void WindowDisplay(void);
void MouseClick(int button, int state, int x, int y);
void MotionFunc(int x, int y) ;
void PassiveMotionFunc(int x, int y) ;
void KeyboardPress(unsigned char key, int x, int y);
void KeyboardPressUp(unsigned char key, int x, int y);

// L-system generation, rendering, and animation
struct PlantNode *MakePlant(void);
void GenerateRecursivePlant(struct PlantNode *p, int level);
void FreePlant(struct PlantNode *p);
void PrintPlant(struct PlantNode *p);
void printNodeRecursive(struct PlantNode *p, int lev, int tgt);
void RenderPlant(struct PlantNode *p);
void StemSection();
void LeafSection();
void FlowerSection();
void AnimatedRenderPlant(void);

// Surface generation
void MakeSurfaceGrid(void);
void RenderSurfaceGrid(void);
void computeNormal(double *vx, double *vy, double *vz, double wx, double wy, double wz);

/**************************************************************************
 Program Code starts
**************************************************************************/
void computeNormal(double *vx, double *vy, double *vz, double wx, double wy, double wz)
{
 //
 // I'm giving you a handy small function to compute the normal to a surface
 // given two vectors known to be on that surface (e.g. two sides of a GL_QUAD
 // or GL_TRIANGLE).
 //
 // Notice that it returns the normal in vx, vy, vz, hence these are expected
 // to be pointers
 //
 // Example call:
 //     vx = .25;
 //     vy = .12;
 //     vz = -.1;
 //     wx = .5;
 //     wy = -.3;
 //     wz = .25;
 //     computeNormal(&vx, &vy, &vz, wx, wy, wz);
 //
 // The returned normal has unit length.
 //

 double len;
 double nx,ny,nz;

 len=sqrt(((*vx)*(*vx))+((*vy)*(*vy))+((*vz)*(*vz)));
 *(vx)=(*vx)/len;
 *(vy)=(*vy)/len;
 *(vz)=(*vz)/len;

 len=sqrt((wx*wx)+(wy*wy)+(wz*wz));
 wx/=len;
 wy/=len;
 wz/=len;
et
 nx=((*vy)*wz)-(wy*(*vz));
 ny=(wx*(*vz))-((*vx)*wz);
 nz=((*vx)*wy)-(wx*(*vy));

 *(vx)=nx;
 *(vy)=ny;
 *(vz)=nz;
}

void RenderSurfaceGrid(void)
{
 // Render the surface grid defined by the vertices in GroundXYZ

 /////////////////////////////////////////////////////////////////////////
 // TODO: Write code to draw the surface map you generated.
 //       Remember that you have vertices on a square grid. You can
 //       easily determine (or you can look at the notes!) how to
 //       create GL_QUADS, or GL_TRIANGLES from the vertices in the
 //       grid (but we know that only one of these types of shapes
 //       can be reliably assumed to be flat given its vertices)
 //
 //       Don't forget to specify the normal at each vertex. Otherwise
 //       your surface won't be properly illuminated
 /////////////////////////////////////////////////////////////////////////
}

void MakeSurfaceGrid(void)
{
 // Generate an interesting surface to place the plants on

 /////////////////////////////////////////////////////////////////////////
 // TODO: Write code to generate a surface map. I am already giving you
 //       a skeleton that fills an array of vertices corresponding to
 //       locations on a square grid. The square grid is defined on the
 //       XY plane, and your joy is to determine the value of Z at each
 //       location (i.e. determine the height of the surface at that
 //       point).
 //
 //       You can use any method you like, but I expect to see an
 //       interesting surface (not random, and not a simple sphere, plane,
 //       or saddle. Use your imagination and knowledge of parametric
 //       surfaces, or do a bit of research into terrain generation.
 //
 //       Do not forget to set the normals to the surface!
 //
 //       The surface coordinates are stored in GroundXYZ[][][]
 //       The normals are stored in GroundNormals[][][]
 //
 //       NOTE: when you assign locations to the plants in your forest,
 //             make sure the pant's root location agrees with the
 //             surface height at the point where it is rooted.
 //             We don't want plants sinking into the ground!
 /////////////////////////////////////////////////////////////////////////
 double side;
 double vx,vy,vz,wx,wy,wz;

 // Assign surface heights
 side=15;				// Width of the surface - X and Y coordinates
					// will have values in [-side/2, side/2]
 for (int i=0; i<GRID_RESOLVE; i++)
  for (int j=0; j<GRID_RESOLVE; j++)
  {
   GroundXYZ[i][j][0]=(-side*.5)+(i*(side/GRID_RESOLVE));
   GroundXYZ[i][j][1]=(-side*.5)+(j*(side/GRID_RESOLVE));
   GroundXYZ[i][j][2]=0;	// <----- HERE you must define surface height in some smart way!
  }

 // Compute normals at each vertex
 // Remember we talked about how to compute the normal for a triangle in lecture. You
 // can do the same thing here.
 //
 // NOTE: Be careful with indexing along the borders of the surface grid! you will
 //       run into all sorts of problems if you don't think carefully what you're
 //       doing.
 for (int i=0; i<GRID_RESOLVE; i++)
  for (int j=0; j<GRID_RESOLVE; j++)
  {
   // Obtain two vectors on the surface the point at GroundXYZ[i][j][] is located

   // Then compute the normal
   computeNormal(&vx,&vy,&vz,wx,wy,wz);

   // And store it...
   GroundNormals[i][j][0]=0;    // <----- HEY!
   GroundNormals[i][j][1]=0;    // <----- REPLACE THESE COMPONENTS with the correct
   GroundNormals[i][j][2]=1;    // <----- normal for your surface!
  }
}

void AnimatedRenderPlant(void)
{
 // This function animates the growth of all plants in the plant
 // forest.
 // It retains control as long as the animation is incomplete,
 // drawing each frame and swapping buffers as needed. Once
 // the animation is completed, control returns to the
 // display function.

 // More information is provided under 'CRUNCHY' in the display
 // function...
}

void RenderPlant(struct PlantNode *p)
{
 // Recursive rendering function for the plant. Renders the
 // section specified by p after performing the transformations
 // required to give this part the correct orientation and
 // position.

 ////////////////////////////////////////////////////////////
 // TO DO: Complete this function to draw the plant as a
 //        hierarchical structure. Hint: This will involve
 //        performing a tree-traversal in the proper order
 //
 // NOTE: Since the plant has a pre-defined number of levels,
 //       many 'branches' will end in a symbol corresponding
 //       to a stem section (a or b). You should draw
 //       something at the end of these last-level stems,
 //       else your plant will look 'dried up'.
 ////////////////////////////////////////////////////////////

 if (p==NULL) return;		// Avoid crash if called with empty node
}

void StemSection(void)
{
  // Draws a single stem section, along the current local Z axis
  // I'm giving you this function already so you can at least see
  // the 'skeleton' of your plant to help debug the L-system.

  // Create a quadrics object to make the stem
  GLUquadric *quadObject;
  quadObject=gluNewQuadric();

  gluCylinder(quadObject,.05,.04,1,10,10);

  // Destroy our quadrics object
  gluDeleteQuadric(quadObject);
}

void LeafSection(void)
{
 // Draws a single leaf, along the current local Z axis
 // Note that we draw a little stem before the actual leaf.
 glColor3f(.25,1,.1);
 StemSection();
 // Perhaps you should translate now? :)

 ////////////////////////////////////////////////////////////
 // TO DO: Draw your own leaf design.
 //        It should be aligned with the current Z
 //        axis, and all transformations for positioning
 //        and orienting the leaf in the plant should be
 //        done outside. However, *any* transformations
 //        required to actually draw the leaf should be
 //        done here.
 //
 //        You must draw a leaf using polygons (or quads, or
 //        triangles), and must design it yourself. You
 //        are not allowed to use GLUT objects to draw the
 //        leaf. Also, rotate away from the stem so that the
 //        leaf will point *away* from the growing plant
 //
 //        Note you must provide proper normal vectors for
 //        vertices in your leaf so that it can be correctly
 //        illuminated by OpenGL.
 //
 //        How to obtain the leaf's vertext coordinates?
 //        I use quadriculated paper...
 ////////////////////////////////////////////////////////////

 ////////////////////////////////////////////////////////////
 // CRUNCHY: Use texture mapping to create nicer leaves!
 //          I am setting up the OpenGL texture mapping
 //          configuration for you, your work is in creating
 //          the polygon shape, normals, and texture coordinates
 //          for your leaf. This may in fact end up being
 //          simpler than defining a very complex, non-textured
 //          leaf.
 ////////////////////////////////////////////////////////////

 // Enable texture mapping if needed (see main() to enable texturing)
 if (textures_on)
 {
  glEnable(GL_TEXTURE_2D);
  // Enable Alpha-blending
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_CULL_FACE);
  glBindTexture(GL_TEXTURE_2D,l_texture);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
 }

 ///////////////////////////////////////////////////////////
 // DO YOUR DRAWING WORK HERE!!!!
 ///////////////////////////////////////////////////////////

 // Disable texture mapping
 if (textures_on)
 {
  glDisable(GL_CULL_FACE);
  glDisable(GL_TEXTURE_2D);
  glDisable (GL_BLEND);
 }

}

void FlowerSection()
{
 // Draws a flower perpendicular to the current local Z axis

 /////////////////////////////////////////////////////////////
 // TO DO: Add code to draw a flower,
 //        It should be perpendicular to the Z axis, i.e.
 //        it should have the same orientation as the local
 //        x-y plane.
 //        Same conditions about transformations apply as for
 //        leaves.
 //
 //        DO NOT use GLU and GLUT shapes for this, design
 //        your own leaf shape using standard GL polygons
 //        and lines. Mind the fact you will have to define
 //        surface normals for your shapes otherwise
 //        illumination will not work
 /////////////////////////////////////////////////////////////

 /////////////////////////////////////////////////////////////
 // CRUNCHY: Use texture mapping to create nicer flowers!
 //          Should be easy if you already texture mapped
 //          the leaves.
 /////////////////////////////////////////////////////////////

 // Enable texture mapping (you must also enable it in main()! )
 if (textures_on)
 {
  glEnable(GL_TEXTURE_2D);
  // Enable Alpha-blending
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_CULL_FACE);
  glBindTexture(GL_TEXTURE_2D,p_texture);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
 }

 /////////////////////////////////////////////////////////////
 // DO YOUR DRAWING WORK HERE!!
 /////////////////////////////////////////////////////////////

 // Disable texture mapping
 if (textures_on)
 {
  glDisable(GL_CULL_FACE);
  glDisable(GL_TEXTURE_2D);
  glDisable (GL_BLEND);
 }
}

void FreePlant(struct PlantNode *p)
{
 // Release all memory allocated to the plant structure
 if (p==NULL) return;
 if (p->left!=NULL) FreePlant(p->left);
 if (p->right!=NULL) FreePlant(p->right);
 free(p);
}

void PrintPlant(struct PlantNode *p)
{
 // Tree traversal to print the plant structure
 // To enable viewing the structure of the plant,
 // Nodes are printed by level, note that this is
 // different from what you may need to do to
 // render the plant (which would be a depth-first
 // hierarchical structure)
 int l;

 for (l=0;l<n_levels;l++)
 {
  fprintf(stderr,"Level %02d\n",l);
  printNodeRecursive(p,0,l);
  fprintf(stderr,"\n");
 }

}

void printNodeRecursive(struct PlantNode *p, int lev, int tgt)
{
 if (p==NULL) return;
 if (lev==tgt)
 {
  fprintf(stderr,"%c ",p->type);
  return;
 }

 if (p->left!=NULL) printNodeRecursive(p->left,lev+1,tgt);
 if (p->right!=NULL) printNodeRecursive(p->right,lev+1,tgt);

 fprintf(stderr," ");
 return;
}

struct PlantNode *MakePlant(void)
{
  // This function creates a plant with the specified
  // characteristics given by the global plant parameters
  // and a set of internal parameters that govern the
  // probability of specific parts being generated
  // at each level.
  //
  // The rules for generating the plant are as follows:
  //
  // Four types of components:  a - Plant stem
  //                            b - Plant stem
  //                            c - Leaf
  //                            d - flower
  //
  // You will use the node types to render each plant
  // component appropriately.
  //
  // The generation rules for the plant structure are
  //
  // a -> aa with p = Paaa
  // a -> ab with p = Paab
  // a -> ac with p = Paac
  // a -> ad with p = Paad
  // a -> cd with p = Pacd
  // b -> a  with p = Pba
  // b -> c  with p = Pbc
  // b -> d  with p = Pbd
  // c and d are terminal nodes.
  //
  // The probabilities for transitions are parameters to
  // the function so that you can generate different
  // types of plants. 
  //
  // NOTE: The transition probabilities from a *MUST*
  // add up to 1. The transition probabilities from
  // b *MUST* also add up to 1.
  // Also, because of the probabilistic nature of the
  // generation process, you can not guarantee the plant
  // will have the specified numer of levels!
  //
  // If you add your own node types, you will have to 
  // update the generation rules, and add appropriate
  // parameters to determine transition probabilities
  // for any additional components.
  //
  // Finally, the random number seed is set in main(), 
  // but could be made an input parameter to obtain more
  // random behaviour.

  struct PlantNode *p_root=NULL;

  p_root=(struct PlantNode *)calloc(1,sizeof(struct PlantNode)); // Allocate new tree root
  p_root->type='a';                                              // Always start with a stem
  p_root->z_ang=0;                                               // which is vertical w.r.t
  p_root->x_ang=0;                                               // global coordinate frame
  p_root->scl=.8+(.2*drand48());                                 // Initial scale (defines the size
                                                                 // of the largest component).
  p_root->left=NULL;
  p_root->left=NULL;

  // Generative part, use the laws above to generate up to two children for each node in the tree
  GenerateRecursivePlant(p_root,1);  // Level 0 is the root, next level is 1
  return(p_root);
}

void GenerateRecursivePlant(struct PlantNode *p, int level)
{
  // Generates child nodes for parent 'p' with p being at level 'level'
  // if level>n_levels, stop.
  float dice;
  struct PlantNode *q, *r;
  
  q=r=NULL;

  if (p==NULL) return;                     // Reached a terminal
  if (level>=n_levels) return;             // Reached maximum plant height
  if (p->type=='c'||p->type=='d') return;  // c and d type nodes are terminal nodes as well

  dice=drand48();           // Roll the dice...
  if (p->type=='a')
  {
   /////////////////////////////////////////////////////////////
   // TO DO: Complete this part, select a replacement rule for
   //        this node based on the probabilities stated above
   //        (also found in the handbook), and add the appropriate
   //        nodes to the plant tree. This is implemented below
   //        for 'b' type nodes, and you can look at that code
   //        to give you an idea how the process works.
   ///////////////////////////////////////////////////////////// 
  }
  else if (p->type=='b')
  {
    // Generate a single node for either left or right
    q=NULL;
    r=(struct PlantNode *)calloc(1,sizeof(struct PlantNode));
    r->x_ang=drand48()*X_angle;
    r->z_ang=drand48()*Z_angle;
    r->scl=scale_mult;
    r->left=NULL;
    r->right=NULL;

    if (dice<=Pba)
      {
	// Selected rule b -> a
        r->type='a';
      }
    else if (dice<=(Pba+Pbc))
      {
        // Selected rule b -> c
        r->type='c';
      }
    else
      {
        // Selected rule b -> d
        r->type='d';
      }
  }
  else {fprintf(stderr,"Bad node type!\n"); return;}

  // Decide which node goes to left and which goes right randomly
  if (drand48()<=.5)
    {
      p->left=q;
      p->right=r;
    }
  else
    {
      p->left=r;
      p->right=q;
    }

  // Recursive call for children
  GenerateRecursivePlant(p->left,level+1);
  GenerateRecursivePlant(p->right,level+1);
}

int main(int argc, char** argv)
{
 /*
   Parse input line parameters, enforce reasonable bounds on global variables
   and parameters for the L-system, and set up plant structures.
 */

    // Process program arguments
    if(argc != 15) {
        printf("Usage: PlantLife n_plants n_levels X_angle Z_angle scale_mult Paab Paac Paad Pacd Pba Pbc Pbd width height\n");
        exit(1);
    } else {
        n_plants=atoi(argv[1]);
        n_levels=atoi(argv[2]);
        X_angle=atof(argv[3]);
        Z_angle=atof(argv[4]);
        scale_mult=atof(argv[5]);
        Paab=atof(argv[6]);
        Paac=atof(argv[7]);
        Paad=atof(argv[8]);
        Pacd=atof(argv[9]);
        Pba=atof(argv[10]);
        Pbc=atof(argv[11]);
        Pbd=atof(argv[12]);
        Win[0] = atoi(argv[13]);
        Win[1] = atoi(argv[14]);

        // Enforce bounds on input variables
        if (n_plants>=MAX_PLANTS) n_plants=MAX_PLANTS;
        if (n_plants<=0) n_plants=1;
        if (n_levels<3) n_levels=3;
        if (n_levels>12) n_levels=12;
        if (X_angle<10) X_angle=50;
        if (X_angle>90) X_angle=90;
        if (Z_angle<10) Z_angle=10;
        if (Z_angle>360) Z_angle=360;
        if (scale_mult<.75) scale_mult=.75;
        if (scale_mult>.99) scale_mult=.99;
	Paab=Paab/(Paab+Paac+Paad+Pacd);
	Paac=Paac/(Paab+Paac+Paad+Pacd);
	Paad=Paad/(Paab+Paac+Paad+Pacd);
	Pacd=Pacd/(Paab+Paac+Paad+Pacd);
        Pba=Pba/(Pba+Pbc+Pbd);
        Pbc=Pbc/(Pba+Pbc+Pbd);
        Pbd=Pbd/(Pba+Pbc+Pbd);
        if (Win[0]<250) Win[0]=250;
        if (Win[0]>1024) Win[0]=1024;
        if (Win[1]<250) Win[1]=250;
        if (Win[1]>1024) Win[1]=1024;

	////////////////////////////////////////////////
        // CRUNCHY - If you are going to use textures
        //           for your leafs and flowers, update
        //           the code below.
        //           Make sure the input images are
        //           available and have a square size
        //           which is a power of 2. A reasonable
        //           size would be 256x256, don't use
        //           huge textures or you'll pay in
        //           rendering performance.
	//
        //           You MUST
        //           submit your texture images along
        //           with your completed code.
        ////////////////////////////////////////////////
	textures_on=0;		// Set to 1 to enable texturing
        if (textures_on)
        {
	 leaf_texture=readPPM("leaf_texture_image.ppm",&l_sx,&l_sy);	// Evidently, you must change this to be
									// your leaf texture image in .ppm format!
	 petal_texture=readPPM("petal_texture_image.ppm",&p_sx,&p_sy);	// Similarly, set this to be your petal
									// texture image.
         if (!leaf_texture||!petal_texture)
         {
          fprintf(stderr,"main(): Unable to load textures. Texture mapping disabled\n");
          textures_on=0;
         }
         else fprintf(stderr,"Textures read and stored\n");
        }
    }

    // Initialize OpenGL - Take a moment to read through these functions!
    glutInit(&argc, argv);
    initGlut(argv[0]);

    // Initialize all data arrays
    memset(&GroundXYZ[0][0][0],0,GRID_RESOLVE*GRID_RESOLVE*3*sizeof(GL_FLOAT));
    memset(&GroundNormals[0][0][0],0,GRID_RESOLVE*GRID_RESOLVE*3*sizeof(GL_FLOAT));
    memset(&ForestXYZ[0][0],0,n_plants*3*sizeof(GL_FLOAT));

    // Generate surface map
    MakeSurfaceGrid();

    // Make a plant forest!
    for (int i=0;i<n_plants;i++)
     PlantForest[i]=MakePlant();

    //////////////////////////////////////////////////////////////
    // TO DO: Set the locations of the plants in the plant forest
    //        randomly in X,Y, but at the correct height for
    //        the corresponding location in the surface grid.
    //////////////////////////////////////////////////////////////

    // Intialize global transformation variables and GLUI    
    global_Z=0;
    global_scale=15;
    ImGui_ImplGlut_Init(false);

    // Invoke the standard GLUT main event loop
    glutMainLoop();
    ImGui_ImplGlut_Shutdown();
    return 0;         // never reached
}

void timerFunc(int) {
    glutTimerFunc(13,timerFunc,13);
    glutPostRedisplay();
}

// Initialize glut and create a window with the specified caption
void initGlut(char* winName)
{
    // Set video mode: double-buffered, color, depth-buffered
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);

    // Create window
    glutInitWindowPosition (0, 0);
    glutInitWindowSize(Win[0],Win[1]);
    windowID = glutCreateWindow(winName);

    // Setup callback functions to handle events
    glutReshapeFunc(WindowReshape);   
    glutDisplayFunc(WindowDisplay);
    glutMouseFunc(MouseClick);
    glutMotionFunc(MotionFunc);
    glutKeyboardFunc(KeyboardPress);
    glutKeyboardUpFunc(KeyboardPressUp);
    glutPassiveMotionFunc(PassiveMotionFunc);

    glutTimerFunc(13,timerFunc,13);

    // Texturing stuff - load textures to the graphics-card memory
    // once!
    if (textures_on)
    {
     glGenTextures( 1, &l_texture);
     glBindTexture( GL_TEXTURE_2D, l_texture);

     glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
     glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
     glTexImage2D (GL_TEXTURE_2D, 0, GL_RGBA, l_sx, l_sy, 0, GL_RGBA, GL_UNSIGNED_BYTE, leaf_texture);

     glGenTextures( 1, &p_texture);
     glBindTexture( GL_TEXTURE_2D, p_texture);

     glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
     glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
     glTexImage2D (GL_TEXTURE_2D, 0, GL_RGBA, p_sx, p_sy, 0, GL_RGBA, GL_UNSIGNED_BYTE, petal_texture);
    }
}

// Handles the window being resized by updating the viewport
// and projection matrices
void WindowReshape(int w, int h)
{
    // Setup projection matrix for new window
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // We will use perspective projection. The camera is at (100,100,100)
    // looking at (0,0,0) - the origin -, with the Z axis pointing upward
    gluPerspective(60,1,15,500);
    gluLookAt(150,150,150,0,0,50,0,0,1);

    // Update OpenGL viewport and internal variables
    glViewport(0,0, w,h);
    Win[0] = w;
    Win[1] = h;
}

// Quit button handler.  Called when the "quit" button is pressed.
void quitButton(int)
{
  for (int i=0; i<n_plants; i++) FreePlant(PlantForest[i]);
  exit(0);
}

// Initialize GLUI and the user interface
void setupUI()
{
    glDisable(GL_LIGHTING);
    ImGui_ImplGlut_NewFrame();
    ImGui::Begin("PlantLife Window");

    ///////////////////////////////////////////////////////////
    // TO DO: Add the controls for global rotation and scale
    //        as specified. Variables are already provided:
    //        global_Z      <--- global rotation around Z
    //        global_scale  <--- global scaling
    //
    //        global_Z must be in [-180, 180]
    //        global_scale must be in [0, 20]
    ///////////////////////////////////////////////////////////

    ImGui::SetWindowFocus();
        ImGui::ColorEdit3("clear color", (float*)&clear_color);

    // Add "Quit" button
    if(ImGui::Button("Quit")) {
        quitButton(0);
    }

    //Some extra info
    ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);


    //End window
    ImGui::End();



    ImGui::Render();
    glEnable(GL_LIGHTING);
}


void drawAxisLines(void)
{
  // Generate a set of axis at current origin so we can visualize the effect
  // of the transformations... this may be useful to you when trying to
  // figure out what is going on at different points in the plant rendering
  // process.

  // Axes are color coded X=R, Y=G, Z=B
  glColor3f(1.0,0.0,0.0);
  glBegin(GL_LINES);
  glVertex3f(0.0f, 0.0f, 0.0f);
  glVertex3f(500.0f, 0.0f, 0.0f);
  glEnd( );

  glColor3f(0.0,1.0,0.0);
  glBegin(GL_LINES);
  glVertex3f(0.0f, 0.0f, 0.0f);
  glVertex3f(0.0f, 500.0f, 0.0f);
  glEnd( );

  glColor3f(0.0,0.0,1.0);
  glBegin(GL_LINES);
  glVertex3f(0.0f, 0.0f, 0.0f);
  glVertex3f(0.0f, 0.0f, 500.0f);
  glEnd( );
}

// Main drawing function. Callback for scene display
void WindowDisplay(void)
{
    static int Opening_animation=0;    
//    static int Opening_animation=1;	// Comment the line above and uncomment this line
                                        // if you implemented the plant growing animation.    

    // Clear colour buffer and Z-buffer
    glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // Setup the model-view transformation matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glClearDepth(1);
    glEnable(GL_DEPTH_TEST);    // Enable depth testing
    glEnable(GL_LIGHTING);      // Enable lighting
    glEnable(GL_LIGHT0);        // Set up 1 light sources, 1 for diffuse,
    glEnable(GL_LIGHT1);        // 1 for ambient

    glEnable(GL_NORMALIZE);	// Make sure normals stay normalized...

    // Set up light source colour, type, and position
    GLfloat light0_colour[]={.95,.95,.95};
    GLfloat light1_colour[]={.05,.05,.05};
    glLightfv(GL_LIGHT0,GL_DIFFUSE,light0_colour);
    glLightfv(GL_LIGHT1,GL_AMBIENT,light1_colour);
    GLfloat light0_pos[]={2,2,5,0};
    glLightfv(GL_LIGHT0,GL_POSITION,light0_pos);
    glShadeModel(GL_SMOOTH);

    // Enable material colour properties
    glEnable(GL_COLOR_MATERIAL);

    // First time through this function, animate the plants growing (if implemented!)

    ////////////////////////////////////////////////////////////////
    // CRUNCHY: Write a function to animate the plants as they grow.
    //          this function should be called 
    //
    //          AnimatedRenderPlant()
    //
    //          and should produce the illusion of plants growing.
    //          Note that a global scaling of the plant WILL NOT
    //          be acceptable. You need to grow the plant level
    //          by level.
    //
    //          The PrintPlant() function may give you ideas
    //          about how to do this.
    //
    //          If you implement this part. Go back to the top
    //          of this function, and uncomment the correct
    //          line that sets the opening animation flat to 1.
    ///////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////
    // SUPREME CRUNCHYNESS OF DOOM: (a.k.a. if you do this Paco will
    //                               be VERY impressed)
    //
    // a)  Add a button to your UI called 'scatter'. When pressed,
    //     all leafs and flowers become Boids and fly away in the
    //     typical Boid patterns from A1. Of course, you can use
    //     any of your A1 code here, and the Boid drawing function
    //     would simply be the same function that draws a leaf
    //     or a flower.
    //
    //    Notice that once the Boids are flying, the plants should
    //     have no leaves or flowers! You may want to scale your
    //     leafs/flowers so as to avoid having a bunch of huge
    //     object/Boids flying around.
    //
    //  or
    //
    //  b) Simulate the changing seasons. Leaves/petals should
    //     change colour accordingly, and fall to the ground
    //     in the fall. There should be a UI button to activate
    //     or trigger this part.
    //
    //     If you are really feeling crunchy, add some falling
    //     snow!
    //
    //     a) and b) are not mutually exclusive. More crunchy
    //     means more bonus!
    //
    //    If you are implementing this and get stuck, come and 
    //     talk to me.
    ///////////////////////////////////////////////////////////////

    if (Opening_animation) {AnimatedRenderPlant(); Opening_animation=0;}
    else
    {
     // synchronize variables that GLUT uses
     glPushMatrix();
      glScalef(global_scale,global_scale,global_scale);
      glRotatef(global_Z,0,0,1);
      RenderSurfaceGrid();
      for (int i=0; i<n_plants; i++)
      {
       glPushMatrix();
        glTranslatef(ForestXYZ[i][0],ForestXYZ[i][1],ForestXYZ[i][2]);
        RenderPlant(PlantForest[i]);
       glPopMatrix();
      }
     glPopMatrix();
    }

    setupUI();
    glFlush();
    glutSwapBuffers();
}

// Utility to read a .ppm file from disk for texture mapping
unsigned char *readPPM(const char *name, int *sx, int *sy)
{
 // Reads an image from a .ppm file. A .ppm file is a very simple image representation
 // format with a text header followed by the binary RGB data at 24bits per pixel.
 // The header has the following form:
 //
 // P6
 // # One or more comment lines preceded by '#'
 // 340 200
 // 255
 //
 // The first line 'P6' is the .ppm format identifier, this is followed by one or more
 // lines with comments, typically used to inidicate which program generated the
 // .ppm file.
 //
 // After the comments, a line with two integer values specifies the image resolution
 // as number of pixels in x and number of pixels in y.
 //
 // The final line of the header stores the maximum value for pixels in the image,
 // usually 255 but does not actually matter for this code as long as it's present.
 //
 // After this last header line, binary data stores the RGB values for each pixel
 // in row-major order at 24bpp.
 //
 // readPPMdata converts the image pixel data to RGBA format by treating all
 // white pixels ([RGB]=[255,255,255]) as fully transparent. It returns an
 // array of sx * sy *4 bytes, where each pixel now has RGBA values.
 //
 // The size of the input image is returned in (sx,sy)

 FILE *f;
 unsigned char *im;
 char line[1024];
 int sizx,sizy;
 int i,j,ly;
 unsigned char *tmp;

 f=fopen(name,"r");
 if (!f){fprintf(stderr,"readPPM(): Unable to open specified image file %s\n",name);return(NULL);}

 fgets(&line[0],1000,f);
 if (strcmp(&line[0],"P6\n")!=0)
 {
  fprintf(stderr,"readPPM(): Wrong file format, not a .ppm file or header data missing\n");
  fclose(f);
  return(NULL);
 }
 // Skip over comments
 fgets(&line[0],511,f);
 while (line[0]=='#')
  fgets(&line[0],511,f);
 sscanf(&line[0],"%d %d\n",&sizx,&sizy);           // Read file size

 *(sx)=sizx;
 *(sy)=sizy;

 im=(unsigned char *)calloc(sizx*sizy*4,sizeof(unsigned char));
 tmp=(unsigned char *)calloc(sizx*sizy*3,sizeof(unsigned char));
 if (!im||!tmp)
 {
  fprintf(stderr,"readPPM(): Unable to allocate memory for image data!\n");
  free(im);
  fclose(f);
  return(NULL);
 }
 fgets(&line[0],9,f);                          		// Read the remaining header line
 fread(tmp,sizx*sizy*3*sizeof(unsigned char),1,f);	// Read image data
 fclose(f);

 // Convert to RGBA
 for (ly=0;ly<3;ly++)
  for (j=0;j<sizy;j++)
   for (i=0;i<sizx;i++)
   {
    *(im+((i+(j*sizx))*4)+0)=*(tmp+((i+(j*sizx))*3)+0);
    *(im+((i+(j*sizx))*4)+1)=*(tmp+((i+(j*sizx))*3)+1);
    *(im+((i+(j*sizx))*4)+2)=*(tmp+((i+(j*sizx))*3)+2);
    if (*(tmp+((i+(j*sizx))*3)+0)==255 && *(tmp+((i+(j*sizx))*3)+1)==255 && *(tmp+((i+(j*sizx))*3)+2)==255)
     *(im+((i+(j*sizx))*4)+3)=0;
    else
     *(im+((i+(j*sizx))*4)+3)=192;
   }

 free(tmp);
 return(im);
}

void MouseClick(int button, int state, int x, int y) {
    ImGui_ImplGlut_MouseButtonCallback(button, state, x, y);
}
void MotionFunc(int x, int y) {
    ImGui_ImplGlut_MotionCallback(x, y);
}
void PassiveMotionFunc(int x, int y) {
    ImGui_ImplGlut_PassiveMotionCallback(x, y);
}
void KeyboardPress(unsigned char key, int x, int y) {
    ImGui_ImplGlut_KeyCallback(key,x,y);
}
void KeyboardPressUp(unsigned char key, int x, int y) {
    ImGui_ImplGlut_KeyUpCallback(key,x,y);
}
