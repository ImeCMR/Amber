/*
   This is an adaptation of Paul Beroza's "molsurf" program;
   modified by dac to exist as an NAB subroutine that gets passed
   a molecule, and atom expression (indicating which atoms are to
   be used), and a probe radius.  Returns the total molecular surface 
   area, using the coordinates and radii in from the protion of the
   input molecule that matches aex.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "nabc.h"
#include "molsurf.h"
#include "./molutil.h"

static int natm_sel;
static ATOM *atom;
static RES *res;

/* neighbor arrays:  these are big so amount of data stored must be small
   upper_neighbors is of the NEIGHBOR_TORUS type, which contains 2
   small ints, one for the index of the neighbor the other for the
   torus.  This last index is key.  The torus associated with
   the neighbor is unique (the index of the upper_neighbor atom is always 
   larger than the atom that points to it).  If the index for this torus
   is -1, then no torus has been instantiated yet, so we don't have to
   allocated the memory for it.  It the index is -2 the torus is buried.
   Any other index corresponds to an index in the torus array that has
   already been instantiated.
 */

static NEIGHBOR_TORUS *upper_neighbors; /* contains atoms and torus indices */
static NEIGHBOR *neighbors;     /* contains atom indices for all neighbors */
static TORUS *toruslist;
static PROBE *probelist;

static CONCAVE_FACE *concave_face;
static SADDLE_FACE *saddle_face;
static CONVEX_FACE *convex_face;
static CONE_FACE *cone_face;
static BROKEN_CONCAVE_FACE *broken_concave_face;
static CONCAVE_CYCLE *concave_cycle;

static VERTEX *vertexlist;
static EDGE *concave_edge_list;
static EDGE *convex_edge_list;
static CIRCLE *convex_circle_list;
static CIRCLE *concave_circle_list;

static CYCLE *cyclelist;
static LOW_TORUS *low_torus;
static CUSP_EDGE *cusp_edge;
static CUSP_PAIR *cusp_pair;

/*  size limitations: these are multiplied by the number of atoms */
#define NUM_NEIGHBOR 100
#define NUM_PROBE 100
#define NUM_TORUS 10
#define NUM_CIRCLE 30
#define NUM_FACE 20
#define NUM_CYCLE 20
#define NUM_VERTEX 20
#define NUM_EDGE 20
#define NUM_CUSP 20

static void check_broken_faces(int, BROKEN_CONCAVE_FACE[]);
static void add_edge(int *, EDGE[], int, int, int, VERTEX[], CIRCLE[]);
static void add_free_edge(int *, EDGE[], int);
static void add_probe(int *, PROBE[], POINT, REAL_T,
                      int, int, int, int, int, int, NEIGHBOR_TORUS[],
                      REAL_T, ATOM[]);
static void add_saddle_face(SADDLE_FACE[], int *, int, int, int, int, int);
static void cross(POINT, POINT, POINT);
#ifdef DEBUG
static void memory_usage(void);
#endif
static void vnorm(REAL_T[], int);
static void atom_vertex_match(VERTEX[], int, int);
static int getneighbors(int, ATOM[], NEIGHBOR[], NEIGHBOR_TORUS[], REAL_T);
static int get_probes(ATOM[], int *, NEIGHBOR[], NEIGHBOR_TORUS[],
                      PROBE[], REAL_T);
static int probe_pos(ATOM *, ATOM *, ATOM *,
                     NEIGHBOR_TORUS *, NEIGHBOR_TORUS *, NEIGHBOR_TORUS *,
                     REAL_T, REAL_T, POINT, POINT, REAL_T *);
static int no_bump(ATOM[], int, int, int,
                   NEIGHBOR[], POINT, POINT, REAL_T, REAL_T, int[]);
static int t_buried(ATOM[], int *, NEIGHBOR[], NEIGHBOR_TORUS[], REAL_T);
static int bury_check(int, int, ATOM[], REAL_T, REAL_T, NEIGHBOR[]);
static int is_buried(ATOM *, ATOM *, ATOM *, REAL_T, REAL_T);
static int get_torus(ATOM[], int *, NEIGHBOR[], NEIGHBOR_TORUS[],
                     REAL_T, TORUS[]);
static void torus_data(TORUS[], int, ATOM[], int, int, REAL_T);
static void sort_neighbors(int, REAL_T[], int, NEIGHBOR[]);
static int convex_circles(ATOM[], int, TORUS[], int, CIRCLE[], REAL_T);
static int concave_circles(ATOM[], int, PROBE[], TORUS[], CIRCLE[],
                           REAL_T);
static void concave_edges(REAL_T, ATOM[], int, PROBE[],
                          int *, VERTEX[], int *, EDGE[], int *,
                          CONCAVE_FACE[], int, TORUS[]);
static void addvert(REAL_T, int, ATOM[], int, PROBE[], int *, VERTEX[]);
static void convex_edges(REAL_T, int, ATOM[], int, PROBE[],
                         int, VERTEX[], int, EDGE[], int *, EDGE[], int,
                         TORUS[], int *, SADDLE_FACE[], CIRCLE[]);
static void check_data(int, TORUS[], int, ATOM[], VERTEX[], EDGE[], int);
static void add_convex_edge_to_torus(int, TORUS[], int, int);
static void add_convex_edge_to_atom(int ia, ATOM atom[], int ne);
static void convex_faces(int, ATOM[], int *, CONVEX_FACE[],
                         int, CYCLE[], EDGE[], CIRCLE[], int, VERTEX[]);
static void sort_edges(int, EDGE[], int, TORUS[], int, VERTEX[], CIRCLE[]);
static void cycles(int, ATOM[], VERTEX[], int, EDGE[],
                   int *, CYCLE[], CIRCLE[], TORUS[]);
static int new_edge(int, ATOM[], int, int[], int);
static int next_edge(int, int, ATOM[], int, EDGE[], int[], int);
static void dump_atom_edges(int, int, ATOM[], int, EDGE[], int[]);
static void draw_circle(int *, int *, int, FILE *, FILE *, FILE *,
                        POINT, REAL_T, POINT, int);
static int is_cycle_inside(int, int, CYCLE[], EDGE[], CIRCLE[],
                           VERTEX[], ATOM[]);
static REAL_T get_angle(POINT, POINT, POINT);
static REAL_T convex_area(ATOM[], RES[], CYCLE[], int, CONVEX_FACE[],
                          EDGE[], CIRCLE[], VERTEX[], TORUS[]);
static REAL_T cycle_piece(int, CYCLE[], CIRCLE[], EDGE[], VERTEX[],
                          ATOM[], TORUS[]);
static int concave_area(REAL_T, int, VERTEX[], CONCAVE_FACE[],
                        EDGE[], CIRCLE[], REAL_T *);
static int saddle_area(int, SADDLE_FACE[], EDGE[], EDGE[],
                       CIRCLE[], TORUS[], ATOM[], RES[], VERTEX[], PROBE[],
                       REAL_T, CIRCLE[], REAL_T *);
static int id_torus(int, int, ATOM[], TORUS[]);
static void make_cones(int, ATOM[], int, TORUS[],
                       int, PROBE probe[], int, CONCAVE_FACE[], int,
                       SADDLE_FACE[], int, CONVEX_FACE[], int *, VERTEX[],
                       int *, EDGE[], int, EDGE[], int, CIRCLE[], int,
                       CIRCLE[], int, CYCLE[], REAL_T, int *, LOW_TORUS[],
                       CONE_FACE[], int *);
static int one_sided_torus(int, TORUS[], ATOM[]);
static void add_2_verts(REAL_T, int, TORUS[], int *, VERTEX[]);
static int get_low_torus_index(int, LOW_TORUS[], int);
static void kill_saddle_face(int, TORUS[], int, SADDLE_FACE[], EDGE[], int,
                             int);
static void cone_init(int, CONE_FACE[], int, int, int);
static void make_broken_faces(int, ATOM[], int, TORUS[], int, PROBE[],
                              int, CONCAVE_FACE[], int, SADDLE_FACE[], int,
                              CONVEX_FACE[], int *, VERTEX[], int, EDGE[],
                              int, EDGE[], int, CIRCLE[], int, CIRCLE[],
                              int, CYCLE[], REAL_T, int *, LOW_TORUS[],
                              int *, BROKEN_CONCAVE_FACE[],
                              CONCAVE_CYCLE[], int *);
static void add_concave_cycle(int, CONCAVE_CYCLE[], int, int, int, int,
                              int);
static REAL_T get_cone_area(int, CONE_FACE[], EDGE[], EDGE[], CIRCLE[],
                            CIRCLE[], TORUS[], ATOM[], VERTEX[], PROBE[],
                            REAL_T);
static void axial_trim(int, ATOM[], RES[], int, TORUS[], int, PROBE[], int,
                       CONCAVE_FACE[], int, VERTEX[], int *, EDGE[], int *,
                       CIRCLE[], REAL_T, int, LOW_TORUS[], int,
                       BROKEN_CONCAVE_FACE[], CONCAVE_CYCLE[], int,
                       CONE_FACE[], CUSP_EDGE[], int *);
static int get_cycle_edge(int, CONCAVE_CYCLE[], int, EDGE[], CIRCLE[]);
static void sort_faces(int, LOW_TORUS[], REAL_T[], int[], int *, EDGE[],
                       VERTEX[], int, TORUS[]);
static void add_circle(PROBE[], int, int, CIRCLE[], int *, REAL_T);
static int cone_edge(int, int, LOW_TORUS[], int, EDGE[], CONE_FACE[],
                     VERTEX[], ATOM[]);
static void add_edges_2_cycle(int *, CUSP_EDGE[], CONCAVE_CYCLE[], int,
                              int, int, int, int, int);
#ifdef DEBUG
static void check_cycle(int, BROKEN_CONCAVE_FACE[], int, CONCAVE_CYCLE[],
                        EDGE[], VERTEX[]);
#endif
static void concentric_axial_cusps(int *, EDGE[], int *,
                                   BROKEN_CONCAVE_FACE[], CONCAVE_CYCLE[],
                                   int *, CUSP_EDGE[], int *);
static void reroute(int, int, int, int, int *, CONCAVE_CYCLE[],
                    CUSP_EDGE[], EDGE[], int *, BROKEN_CONCAVE_FACE[]);
static void split_cycle(int *, BROKEN_CONCAVE_FACE[], int *,
                        CONCAVE_CYCLE[], int, EDGE[], CUSP_EDGE[]);
static int next_cycle_edge(CONCAVE_CYCLE, EDGE[], int, int[]);
static void non_axial_trim(int, ATOM[], RES[], int, TORUS[], int, PROBE[],
                           int, CONCAVE_FACE[], int, VERTEX[], int *,
                           EDGE[], int *, CIRCLE[], REAL_T, int,
                           LOW_TORUS[], int, BROKEN_CONCAVE_FACE[],
                           CONCAVE_CYCLE[], int, CONE_FACE[], CUSP_EDGE[],
                           int *, CUSP_PAIR[], int *);
static int new_cusp_in_group(int, int, CUSP_PAIR[], CUSP_GROUP[]);
static void cusp_intersect(CUSP_EDGE[], int, int, PROBE[], CONCAVE_CYCLE[],
                           EDGE[], VERTEX[], CIRCLE[], REAL_T, CUSP_PAIR[],
                           int *);
static int get_cycle_id(CUSP_EDGE[], int, int, int);
static void get_probe_id(CUSP_EDGE[], int, int, CONCAVE_CYCLE[], int *,
                         int *);
static void add_new_cusp(CUSP_EDGE[], int, int, PROBE[], CONCAVE_CYCLE[],
                         EDGE[], VERTEX[], CIRCLE[], REAL_T, CUSP_PAIR[],
                         int *, REAL_T, POINT, int, POINT, POINT, POINT,
                         int);
static int center_cycle(CUSP_EDGE[], int, int);
static int number_of_cusps(CUSP_GROUP[], int, CUSP_PAIR[]);
static int new_cusp(int, int[], int);
static void trim_2_cusps(PROBE[], int *, VERTEX[], int *, EDGE[],
                         int *, CIRCLE[], REAL_T, int *,
                         BROKEN_CONCAVE_FACE[], CONCAVE_CYCLE[], int *,
                         CUSP_EDGE[], int *, CUSP_PAIR[], int *,
                         CUSP_GROUP[], int);
static void unique_cycles(CUSP_EDGE *, CUSP_EDGE *, int *, int *);
static void split_old_cusps(int[], int, CUSP_EDGE[], int *,
                            EXTREME_VERTEX[], int, CUSP_PAIR[], VERTEX[],
                            int *, EDGE[], CIRCLE[]);
static void make_new_cusp(int *, CUSP_EDGE[], int, int);
static void get_faces(int[], int *, int[], int, CUSP_EDGE[],
                      CONCAVE_CYCLE[], BROKEN_CONCAVE_FACE[]);
static int is_new_face(int, int[], int);
static void split_face(int, int *, BROKEN_CONCAVE_FACE[], int *,
                       CONCAVE_CYCLE[], CUSP_EDGE[], int[], int, int *,
                       EDGE[], VERTEX[], CIRCLE[]);
static void get_starters(int[], int *, CONCAVE_CYCLE[], int, EDGE[]);
static void copy_cycle(CONCAVE_CYCLE *, CONCAVE_CYCLE *);
static int cusp_match(int, int, int[], int, int[], CUSP_EDGE[], EDGE[],
                      int *);
static int normal_match(CONCAVE_CYCLE, EDGE[], int, int[]);
static void add_non_axial_cusp(int *, CUSP_EDGE[], int, int, int, int,
                               int);
static void add_1_vert(int *, VERTEX[], POINT);
static void add_cusp_verts(int *, VERTEX[], CUSP_PAIR[], EXTREME_VERTEX[],
                           int, int, CIRCLE[], int);
static void add_cusp_circle(int *, CIRCLE[], PROBE[], int, int, REAL_T);
static int broken_concave_area(REAL_T, int, BROKEN_CONCAVE_FACE[],
                               CONCAVE_CYCLE[], EDGE[], CIRCLE[], VERTEX[],
                               REAL_T *, PROBE[]);
static REAL_T conc_cycle_piece(int, CONCAVE_CYCLE[], CIRCLE[], EDGE[],
                               VERTEX[], int, PROBE[], REAL_T);
static REAL_T interior_angle(int, int, int, int, EDGE[], CIRCLE[],
                             VERTEX[]);
static void allocate_memory(void);
static void free_memory(void);

REAL_T molsurf(MOLECULE_T ** m, char **aex, REAL_T * probe_rad_in)
{

   int nat = 0, nres = 0, n_torus = 0, n_probes = 0, n_vertex = 0;

   int n_concave_faces = 0, n_saddle_faces = 0, n_cycles = 0;
   int n_convex_faces = 0, n_cone_faces = 0, n_broken_concave_faces = 0;
   int n_concave_cycles = 0, n_concave_edges = 0, n_convex_edges = 0;
   int n_convex_circles = 0, n_concave_circles = 0, n_low_torus = 0;
   int n_cusps = 0;
   int n_cusp_pairs = 0;

   REAL_T conv_area, conc_area, sad_area, cone_area;
   REAL_T broken_conc_area;

   int i;

   int r, a, n;
   REAL_T probe_rad;
   STRAND_T *sp;
   RESIDUE_T *resp;
   ATOM_T *ap;

/*  Take info from NAB molecule and put into Paul's data structure,
   (only for atoms matching aex).                                    */

   natm_sel = countmolatoms(*m, aex ? *aex : NULL);
   allocate_memory();
#ifdef DEBUG
   memory_usage();
#endif

/*  Set up local and compacted arrays of the atoms and their properties
    that will be used for this calculation:                              */

   for (n = 0, sp = (*m)->m_strands; sp; sp = sp->s_next) {
      for (r = 0; r < sp->s_nresidues; r++) {
         nres++;
         resp = sp->s_residues[r];
#ifdef DEBUG
         /* strcpy(res[r].nam, resp->r_name); */
         res[r].num = resp->r_tresnum;
#endif
         for (a = 0; a < resp->r_natoms; a++) {
            ap = &resp->r_atoms[a];
            if (ap->a_attr & AT_SELECT) {
#ifdef DEBUG
               atom[n].anum = n;
               strcpy(atom[n].anam, ap->a_atomname);
               atom[n].res = r + 1;
#endif
               atom[n].pos[0] = ap->a_pos[0];
               atom[n].pos[1] = ap->a_pos[1];
               atom[n].pos[2] = ap->a_pos[2];
               atom[n].q = ap->a_charge;
               atom[n].rad = ap->a_radius;
               atom[n].area = 0.0;
               ap->a_float1 = 0.0;      /* will hold convex surface areas */
               n++;
            }
         }
      }
   }
   probe_rad = *probe_rad_in;
   nat = n;
   if ((res = (RES *) malloc(10 * natm_sel * sizeof(RES))) == NULL) {
      fprintf(stderr, "Unable to allocate space for res\n");
      exit(1);
   }
   i = getneighbors(nat, atom, neighbors, upper_neighbors, probe_rad);

/* determine valid probe positions */
   n_probes = get_probes(atom, &nat, neighbors, upper_neighbors,
                         probelist, probe_rad);

/* check if any tori are completely buried by a single atom */
   i = t_buried(atom, &nat, neighbors, upper_neighbors, probe_rad);

/* identify accessible tori                                                */
   n_torus = get_torus(atom, &nat, neighbors, upper_neighbors,
                       probe_rad, toruslist);

/* create concave circles associated with probes and convex circles
   associated with the atoms */
   n_convex_circles = convex_circles(atom, nat, toruslist, n_torus,
                                     convex_circle_list, probe_rad);
   n_concave_circles =
       concave_circles(atom, n_probes, probelist, toruslist,
                       concave_circle_list, probe_rad);


/* 1. fill up concave edge list, orient edges 
   2. add edge pointers to torus array   
   3. fill up vertex list
   4. create concave face                                                       */
   concave_edges(probe_rad, atom, n_probes, probelist, &n_vertex,
                 vertexlist, &n_concave_edges, concave_edge_list,
                 &n_concave_faces, concave_face, n_torus,
                 toruslist /*, concave_circles */ );
   sort_edges(n_concave_edges, concave_edge_list, n_torus, toruslist,
              n_vertex, vertexlist, convex_circle_list);
   /*
      write_verts(n_vertex, vertexlist, atom);
    */

/* 1. use torus array concave edge pointers to fill up convex edge list 
   2. fill up atom array convex edge pointers
   3. create saddle face                                                        */
   convex_edges(probe_rad, nat, atom, n_probes, probelist,
                n_vertex, vertexlist, n_concave_edges, concave_edge_list,
                &n_convex_edges, convex_edge_list, n_torus, toruslist,
                &n_saddle_faces, saddle_face, convex_circle_list);

/* create cycles                                                         */
   cycles(nat, atom, vertexlist, n_convex_edges, convex_edge_list,
          &n_cycles, cyclelist, convex_circle_list, toruslist);

/* create convex faces                                                   */
   convex_faces(nat, atom, &n_convex_faces, convex_face, n_cycles,
                cyclelist, convex_edge_list, convex_circle_list, n_vertex,
                vertexlist);

/* handle cusps:

   1. create cone faces: whenever there is a cone face, mark
   the two concave edges as dead
   2. if a concave face has a dead edge, kill the face, and
   create a broken_concave_face
   3. initialize the cycles (1 with 3 edges) for the broken conc. faces
   3. trim b.c. faces from low probes in same torus
   4. trim faces from low probes in different tori                      */

   make_cones(nat, atom, n_torus, toruslist, n_probes, probelist,
              n_concave_faces, concave_face, n_saddle_faces, saddle_face,
              n_convex_faces, convex_face, &n_vertex, vertexlist,
              &n_concave_edges, concave_edge_list, n_convex_edges,
              convex_edge_list, n_convex_circles, convex_circle_list,
              n_concave_circles, concave_circle_list, n_cycles, cyclelist,
              probe_rad, &n_low_torus, low_torus, cone_face,
              &n_cone_faces);

   make_broken_faces(nat, atom, n_torus, toruslist, n_probes, probelist,
                     n_concave_faces, concave_face, n_saddle_faces,
                     saddle_face, n_convex_faces, convex_face, &n_vertex,
                     vertexlist, n_concave_edges, concave_edge_list,
                     n_convex_edges, convex_edge_list, n_convex_circles,
                     convex_circle_list, n_concave_circles,
                     concave_circle_list, n_cycles, cyclelist, probe_rad,
                     &n_low_torus, low_torus, &n_broken_concave_faces,
                     broken_concave_face, concave_cycle,
                     &n_concave_cycles);

   check_broken_faces(n_broken_concave_faces, broken_concave_face);

   axial_trim(nat, atom, res, n_torus, toruslist, n_probes, probelist,
              n_concave_faces, concave_face, n_vertex, vertexlist,
              &n_concave_edges, concave_edge_list,
              &n_concave_circles, concave_circle_list,
              probe_rad,
              n_low_torus, low_torus, n_broken_concave_faces,
              broken_concave_face, concave_cycle, n_concave_cycles,
              cone_face, cusp_edge, &n_cusps);

   check_broken_faces(n_broken_concave_faces, broken_concave_face);

   /*
      for (i = 0; i < n_concave_cycles; ++i) {
      printf("cycle %d edges: ", i);
      for (j = 0; j < concave_cycle[i].nedges; ++j)
      printf(" %d", concave_cycle[i].edge[j]);
      printf("\n");
      }
      printf(" %d cusp edges\n", n_cusps); 
      for (i = 0; i < n_cusps; ++i) {
      printf("cusp %d  = edge: %d   cyles: %d %d\n", 
      i, cusp_edge[i].edge, 
      cusp_edge[i].cycle1, cusp_edge[i].cycle2);
      }
    */

   concentric_axial_cusps(&n_concave_edges, concave_edge_list,
                          &n_broken_concave_faces, broken_concave_face,
                          concave_cycle, &n_concave_cycles, cusp_edge,
                          &n_cusps);

   check_broken_faces(n_broken_concave_faces, broken_concave_face);

   non_axial_trim(nat, atom, res, n_torus, toruslist, n_probes, probelist,
                  n_concave_faces, concave_face, n_vertex, vertexlist,
                  &n_concave_edges, concave_edge_list,
                  &n_concave_circles, concave_circle_list,
                  probe_rad,
                  n_low_torus, low_torus, n_broken_concave_faces,
                  broken_concave_face, concave_cycle, n_concave_cycles,
                  cone_face, cusp_edge, &n_cusps, cusp_pair,
                  &n_cusp_pairs);

   check_broken_faces(n_broken_concave_faces, broken_concave_face);
   /*
      draw_edges( n_cycles, atom, n_concave_faces, concave_face, concave_edge_list, 
      cyclelist, convex_edge_list, 
      convex_circle_list,  concave_circle_list,
      vertexlist, probelist, toruslist,
      n_cone_faces, cone_face, 
      n_broken_concave_faces, broken_concave_face,
      n_concave_cycles, concave_cycle,
      n_cusps, cusp_edge);
      exit(ERROR);
    */

   /* calculate areas                                                       */

   concave_area(probe_rad, n_concave_faces, vertexlist,
                concave_face, concave_edge_list,
                concave_circle_list, &conc_area);

   broken_concave_area(probe_rad,
                       n_broken_concave_faces, broken_concave_face,
                       concave_cycle, concave_edge_list,
                       concave_circle_list, vertexlist, &broken_conc_area,
                       probelist);

   conv_area =
       convex_area(atom, res, cyclelist, n_convex_faces, convex_face,
                   convex_edge_list, convex_circle_list, vertexlist,
                   toruslist);

   saddle_area(n_saddle_faces, saddle_face, convex_edge_list,
               concave_edge_list, convex_circle_list, toruslist,
               atom, res, vertexlist, probelist, probe_rad,
               concave_circle_list, &sad_area);

   cone_area = get_cone_area(n_cone_faces, cone_face, convex_edge_list,
                             concave_edge_list, convex_circle_list,
                             concave_circle_list, toruslist, atom,
                             vertexlist, probelist, probe_rad);

#ifdef DEBUG
   printf("need %d probes\n", n_probes);
   printf("need %d torus\n", n_torus);
   printf("need %d convex circles\n", n_convex_circles);
   printf("need %d concave circles\n", n_concave_circles);
   printf("need %d concave edges\n", n_concave_edges);
   printf("need %d vertex\n", n_vertex);
   printf("need %d concave faces\n", n_concave_faces);
   printf("need %d cycles\n", n_cycles);
   printf("need %d convex faces\n", n_convex_faces);

   printf("Contact: %10.3f   Reentrant : %10.3f  Total: %10.3f\n",
          conv_area, conc_area + sad_area + cone_area,
          conv_area + conc_area + sad_area + cone_area);

   printf
       ("\n reentrant: %8.3f  toric: %8.3f  contact: %8.3f Total: %8.3f\n",
        broken_conc_area + conc_area, sad_area + cone_area, conv_area,
        broken_conc_area + conc_area + sad_area + cone_area + conv_area);
#endif

/*   place the convex area parts (which can be associated uniquely with an
     individual atom) back into the molecule in the "float1" attribute:     */

   for (n = 0, sp = (*m)->m_strands; sp; sp = sp->s_next) {
      for (r = 0; r < sp->s_nresidues; r++) {
         resp = sp->s_residues[r];
         for (a = 0; a < resp->r_natoms; a++) {
            ap = &resp->r_atoms[a];
            if (ap->a_attr & AT_SELECT) {
               ap->a_float1 = atom[n].area;
               n++;
            }
         }
      }
   }

   free_memory();
   return (broken_conc_area + conv_area + conc_area + sad_area +
           cone_area);

}

static void check_broken_faces(int n_broken_concave_faces,
                               BROKEN_CONCAVE_FACE broken_concave_face[])
{
   int i;
   for (i = 0; i < n_broken_concave_faces; ++i) {
      if (broken_concave_face[i].n_cycles > 1) {
         printf("FACE CHECK: face %d has %d cycles\n", i,
                broken_concave_face[i].n_cycles);
      }
   }
}

static void add_edge(int *n_edges, EDGE edge[],
                     int v1, int v2, int icircle, VERTEX vertex[],
                     CIRCLE circle[])
{
   int ii;
   REAL_T d1, d2, r2;

   d1 = 0.0;
   d2 = 0.0;
   r2 = circle[icircle].rad * circle[icircle].rad;
   for (ii = 0; ii < 3; ++ii) {
      d1 = d1 + (vertex[v1].pos[ii] - circle[icircle].center[ii]) *
          (vertex[v1].pos[ii] - circle[icircle].center[ii]);
      d2 = d2 + (vertex[v2].pos[ii] - circle[icircle].center[ii]) *
          (vertex[v2].pos[ii] - circle[icircle].center[ii]);
   }
   if ((fabs(d1 - r2) > 0.1) || (fabs(d2 - r2) > 0.1)) {
      printf("add_edge(): edge vertex not on circle\n");
      exit(ERROR);
   }
   /*
      printf("adding edge: v1 %8.3f%8.3f%8.3f\n", vertex[v1].pos[0], vertex[v1].pos[1], vertex[v1].pos[2]);
      printf("adding edge: v2 %8.3f%8.3f%8.3f\n", vertex[v2].pos[0], vertex[v2].pos[1], vertex[v2].pos[2]);
    */

   edge[*n_edges].vert1 = v1;
   edge[*n_edges].vert2 = v2;
   edge[*n_edges].circle = icircle;
   edge[*n_edges].alive = 1;

   (*n_edges)++;
   if (*n_edges >= NUM_EDGE * natm_sel) {
      printf("MAX_EDGE exceeded\n");
      exit(ERROR);
   }
}
static void add_free_edge(int *n_edges, EDGE edge[], int icircle)
{
   edge[*n_edges].vert1 = -1;
   edge[*n_edges].vert2 = -1;
   edge[*n_edges].circle = icircle;
   edge[*n_edges].alive = 1;
   (*n_edges)++;
   if (*n_edges >= NUM_EDGE * natm_sel) {
      printf("MAX_EDGE exceeded\n");
      exit(ERROR);
   }
}
static void add_probe(int *np, PROBE probelist[], POINT px, REAL_T height,
                      int ia, int ja, int ka, int ij, int jk, int ik,
                      NEIGHBOR_TORUS upper_neighbors[], REAL_T probe_rad,
                      ATOM atom[])
{
   POINT r_12, r_13, r_1p, r_cross;
   int ii, counter_clockwise;

/* want to atoms to be arranged counter-clockwise */

   /* determine orientation of concave surface */
   counter_clockwise = 1;
   for (ii = 0; ii < 3; ++ii) {
      r_12[ii] = atom[ja].pos[ii] - atom[ia].pos[ii];
      r_13[ii] = atom[ka].pos[ii] - atom[ia].pos[ii];
      r_1p[ii] = px[ii] - atom[ia].pos[ii];
   }
   cross(r_12, r_13, r_cross);
   if (DOT(r_cross, r_1p) < 0)
      counter_clockwise = 0;

   if (counter_clockwise) {
      probelist[*np].a1 = ia;
      probelist[*np].a2 = ja;
      probelist[*np].a3 = ka;
   } else {
      probelist[*np].a1 = ia;
      probelist[*np].a2 = ka;
      probelist[*np].a3 = ja;
   }
   probelist[*np].pos[0] = px[0];
   probelist[*np].pos[1] = px[1];
   probelist[*np].pos[2] = px[2];
   probelist[*np].height = height;
   if (probelist[*np].height < probe_rad) {
      probelist[*np].low = 1;
   } else {
      probelist[*np].low = 0;
   }
   if (*np > NUM_PROBE * natm_sel) {
      fprintf(stderr, "MAXPROBE exceeded: %d %d %d\n", ia, ja, ka);
      exit(ERROR);
   }
   ++upper_neighbors[ij].nprobes;
   ++upper_neighbors[jk].nprobes;
   ++upper_neighbors[ik].nprobes;
   ++(*np);
}
static void add_saddle_face(SADDLE_FACE saddle_face[], int *nface,
                            int itorus, int e1_concave, int e3_concave,
                            int e2_convex, int e4_convex)
{
   saddle_face[*nface].torus = itorus;
   saddle_face[*nface].e1_concave = e1_concave; /* no concave edges */
   saddle_face[*nface].e3_concave = e3_concave;
   saddle_face[*nface].e2_convex = e2_convex;
   saddle_face[*nface].e4_convex = e4_convex;
   ++(*nface);
   if (*nface >= NUM_FACE * natm_sel) {
      printf("convex_edges() MAX_FACE exceeded\n");
      exit(ERROR);
   }
}

/**********************************************************************/
/* c = a x b */
static void cross(POINT a, POINT b, POINT c)
{
   c[0] = a[1] * b[2] - a[2] * b[1];
   c[1] = a[2] * b[0] - a[0] * b[2];
   c[2] = a[0] * b[1] - a[1] * b[0];
}

/***********************************************************************/

#ifdef DEBUG
static void memory_usage(void)
{

   int total = 0;

   printf("size requirements in bytes\n");

   printf("one atom              %ld\n", sizeof(ATOM));
   printf("atom            total %ld\n", natm_sel * sizeof(ATOM));
   total += natm_sel * sizeof(ATOM);

   printf("one NEIGHBOR_TORUS    %ld\n", sizeof(NEIGHBOR_TORUS));
   printf("NEIGHBOR_TORUS  total %ld\n",
          60 * natm_sel * sizeof(NEIGHBOR_TORUS));
   total += 60 * natm_sel * sizeof(NEIGHBOR_TORUS);

   printf("one NEIGHBOR          %ld\n", sizeof(NEIGHBOR));
   printf("NEIGHBOR        total %ld\n", 60 * natm_sel * sizeof(NEIGHBOR));
   total += 60 * natm_sel * sizeof(NEIGHBOR);

   printf("one TORUS             %ld\n", sizeof(TORUS));
   printf("TORUS           total %ld\n", 5 * natm_sel * sizeof(TORUS));
   total += 5 * natm_sel * sizeof(TORUS);

   printf("one PROBE             %ld\n", sizeof(PROBE));
   printf("PROBE           total %ld\n", 5 * natm_sel * sizeof(PROBE));
   total += 5 * natm_sel * sizeof(PROBE);

   printf("one VERTEX            %ld\n", sizeof(VERTEX));
   printf("VERTEX          total %ld\n", 10 * natm_sel * sizeof(VERTEX));
   total += 10 * natm_sel * sizeof(VERTEX);

   printf("one EDGE      %ld\n", sizeof(EDGE));
   printf("EDGE    total %ld\n", 10 * natm_sel * sizeof(EDGE));
   total += 10 * natm_sel * sizeof(EDGE);

   printf("one CIRCLE            %ld\n", sizeof(CIRCLE));
   printf("CIRCLE          total %ld\n", 5 * natm_sel * sizeof(CIRCLE));
   total += 5 * natm_sel * sizeof(CIRCLE);

   printf("one CONCAVE_FACE      %ld\n", sizeof(CONCAVE_FACE));
   printf("CONCAVE_FACE    total %ld\n",
          10 * natm_sel * sizeof(CONCAVE_FACE));
   total += 10 * natm_sel * sizeof(CONCAVE_FACE);

   printf("one SADDLE_FACE       %ld\n", sizeof(SADDLE_FACE));
   printf("SADDLE_FACE     total %ld\n",
          10 * natm_sel * sizeof(SADDLE_FACE));
   total += 10 * natm_sel * sizeof(SADDLE_FACE);

   printf("one CYCLE             %ld\n", sizeof(CYCLE));
   printf("CYCLE           total %ld\n", 10 * natm_sel * sizeof(CYCLE));
   total += 10 * natm_sel * sizeof(CYCLE);

   printf("Sum total             %d\n\n", total);

   return;
}
#endif

static void vnorm(REAL_T v[], int n)
{
   int i;
   REAL_T vn;

   for (vn = 0.0, i = 0; i < n; i++)
      vn += v[i] * v[i];
   if (vn) {
      vn = sqrt(vn);
      for (i = 0; i < n; i++)
         v[i] /= vn;
   }
}


static void atom_vertex_match(VERTEX vertex[], int iv, int iatom)
{
   if (vertex[iv].iatom != iatom) {
      fprintf(stderr, "vertex atom mismatch\n");
      fprintf(stderr, "       atom: %d\n", iatom);
      fprintf(stderr, "vertex atom: %d\n", vertex[iv].iatom);
      exit(ERROR);
   }
}


#define MAX_NSORT 500

static int getneighbors(int nat, ATOM a[], NEIGHBOR neighbors[],
                        NEIGHBOR_TORUS upper_neighbors[], REAL_T probe_rad)
{
   REAL_T probe_diam, maxrad = 0.0, d, cutoff, d_ext;
   int n_tot = 0, n_upper_tot = 0;
   int i, j;
   int nsort, istart;
   REAL_T dsort[MAX_NSORT];


   probe_diam = 2 * probe_rad;

   for (i = 0; i < nat; ++i)
      if (a[i].rad > maxrad)
         maxrad = a[i].rad;
   cutoff = (2 * maxrad + probe_diam);

   for (i = 0; i < nat; ++i) {
      a[i].buried = 0;
      d_ext = probe_diam + a[i].rad;

      a[i].neighbor_start = n_tot;
      a[i].n_neighbors = 0;
      a[i].upper_start = n_upper_tot;
      a[i].n_upper = 0;

      nsort = 0;
      istart = n_tot;           /* starting index for sort */

      for (j = 0; j < nat; ++j) {

         if (fabs(a[i].pos[0] - a[j].pos[0]) > cutoff)
            continue;
         if (fabs(a[i].pos[1] - a[j].pos[1]) > cutoff)
            continue;
         if (fabs(a[i].pos[2] - a[j].pos[2]) > cutoff)
            continue;
         if (i == j)
            continue;

         /* d = sqrt(dist2(a[i].pos,a[j].pos)); */
         d = DIST(a[i].pos, a[j].pos);

         if (d < (d_ext + a[j].rad)) {
            if (n_tot >= NUM_NEIGHBOR * natm_sel) {
               fprintf(stderr, "MAX_NEIGHBOR exceeded: %d %d %d\n", n_tot,
                       i, j);
               exit(ERROR);
            }
            if (a[i].rad + d < a[j].rad)
               a[i].buried = 1;

            neighbors[n_tot].iatom = j;
            ++a[i].n_neighbors;
            ++n_tot;

            if (j > i) {

               /* 2nd atom in torus */
               upper_neighbors[n_upper_tot].iatom = j;

               /* initial probe count for torus containing 
                  iatom and j to -1 (free) */

               upper_neighbors[n_upper_tot].nprobes = FREE_TORUS;
               ++a[i].n_upper;
               ++n_upper_tot;
               if (n_upper_tot > NUM_NEIGHBOR * natm_sel) {
                  fprintf(stderr, "MAX_NEIGHBOR exceeded\n");
                  exit(ERROR);
               }
            }
            /* set up sort */
            if (nsort >= MAX_NSORT) {
               fprintf(stderr, "MAX_NSORT exceeded: %d %d\n", nsort,
                       MAX_NSORT);
               exit(ERROR);
            }
            dsort[nsort] = d;
            ++nsort;
         }
      }
/*    printf("atom %d has %d neighbors\n", i, a[i].n_neighbors);   */
      sort_neighbors(nsort, dsort, istart, neighbors);
   }

#ifdef DEBUG
   printf("number of temporary tori %d\n", n_upper_tot);
#endif

   return (n_tot);
}


static int get_probes(ATOM atom[], int *nat,
                      NEIGHBOR neighbors[],
                      NEIGHBOR_TORUS upper_neighbors[], PROBE probelist[],
                      REAL_T probe_rad)
{
   int ia, ja, ka, ij, jk, ik, k1, bump[2];
   REAL_T probe_diam;
   POINT px1, px2;
   REAL_T height;
   int nprobes = 0;

   probe_diam = 2 * probe_rad;


   for (ia = 0; ia < *nat; ++ia) {

      if (atom[ia].buried)
         continue;

      for (ij = atom[ia].upper_start;
           ij < atom[ia].upper_start + atom[ia].n_upper; ++ij) {

         if (upper_neighbors[ij].nprobes == BURIED_TORUS)
            continue;           /* buried torus */
         ja = upper_neighbors[ij].iatom;
         if (atom[ja].buried)
            continue;

         for (ik = ij + 1;
              ik < atom[ia].upper_start + atom[ia].n_upper; ++ik) {

            if (upper_neighbors[ik].nprobes == BURIED_TORUS)
               continue;        /* buried torus */
            ka = upper_neighbors[ik].iatom;
            if (atom[ka].buried)
               continue;

            for (jk = atom[ja].upper_start;
                 jk < atom[ja].upper_start + atom[ja].n_upper; ++jk) {

               if (upper_neighbors[jk].nprobes == BURIED_TORUS)
                  continue;     /* buried torus */
               k1 = upper_neighbors[jk].iatom;

               if (k1 == ka) {
                  /* 3 intersecting tori */

                  /*
                     printf("intersecting tori atoms : %d %d %d\n",  ia, ja, ka);
                   */
                  if (probe_pos(&atom[ia], &atom[ja], &atom[ka],
                                &upper_neighbors[ij], &upper_neighbors[jk],
                                &upper_neighbors[ik], probe_rad,
                                probe_diam, px1, px2, &height)) {
                     if (no_bump
                         (atom, ia, ja, ka, neighbors, px1, px2, probe_rad,
                          probe_diam, bump)) {
                        if (!bump[0])
                           add_probe(&nprobes, probelist, px1, height,
                                     ia, ja, ka, ij, jk, ik,
                                     upper_neighbors, probe_rad, atom);
                        if (!bump[1])
                           add_probe(&nprobes, probelist, px2, height,
                                     ia, ja, ka, ij, jk, ik,
                                     upper_neighbors, probe_rad, atom);
                     }
                  }
               }
            }
         }
      }
   }
   return nprobes;
}

static int probe_pos(ATOM * ai, ATOM * aj, ATOM * ak,
                     NEIGHBOR_TORUS * upper_neighbor_ij,
                     NEIGHBOR_TORUS * upper_neighbor_jk,
                     NEIGHBOR_TORUS * upper_neighbor_ik, REAL_T probe_rad,
                     REAL_T probe_diam, POINT p1_ijk, POINT p2_ijk,
                     REAL_T * height)
{
   POINT u_ij, u_ik, u_ijk, u_tb, t_ij, t_ik, t_diff, b_ijk, r_ib;
   REAL_T d_ij, d_ik, t_fact1, t_fact2, w_ijk, w_sin;
   REAL_T b_factor, h_ijk;
   int i;

   for (i = 0; i < 3; ++i) {
      u_ij[i] = aj->pos[i] - ai->pos[i];
      u_ik[i] = ak->pos[i] - ai->pos[i];
   }

   d_ij = sqrt(u_ij[0] * u_ij[0] + u_ij[1] * u_ij[1] + u_ij[2] * u_ij[2]);
   d_ik = sqrt(u_ik[0] * u_ik[0] + u_ik[1] * u_ik[1] + u_ik[2] * u_ik[2]);

   t_fact1 = ((ai->rad + probe_rad) * (ai->rad + probe_rad) -
              (aj->rad + probe_rad) * (aj->rad +
                                       probe_rad)) / (d_ij * d_ij);

   t_fact2 = ((ai->rad + probe_rad) * (ai->rad + probe_rad) -
              (ak->rad + probe_rad) * (ak->rad +
                                       probe_rad)) / (d_ik * d_ik);


   for (i = 0; i < 3; ++i) {

      u_ij[i] = u_ij[i] / d_ij;
      u_ik[i] = u_ik[i] / d_ik;

      t_ij[i] = 0.5 * (ai->pos[i] + aj->pos[i])
          + 0.5 * (aj->pos[i] - ai->pos[i]) * t_fact1;
      t_ik[i] = 0.5 * (ai->pos[i] + ak->pos[i])
          + 0.5 * (ak->pos[i] - ai->pos[i]) * t_fact2;
      t_diff[i] = t_ik[i] - t_ij[i];
   }

   /* if w_sin is ~0 then the atoms are co-linear.  In this case the distance 
      between the probe sphere in any one of the three tori to the  3rd atom
      is constant.  If | r_pk | < r_p + probe_rad then the entire torus t_ij 
      is buried.  If not the torus burial status remains unchanged.  In any case, 
      no (3-contact) probe position can be found, so 0 is returned. Could check 
      for torus burial here, but you need to check this for all k not just those 
      where i<j<k this is done in a separate routine */

   w_ijk = acos(DOT(u_ij, u_ik));
   w_sin = sin(w_ijk);
   if (fabs(w_sin) < 1.0e-10)
      return 0;                 /* no probe possible */

   cross(u_ij, u_ik, u_ijk);
   for (i = 0; i < 3; ++i)
      u_ijk[i] = u_ijk[i] / w_sin;
   cross(u_ijk, u_ij, u_tb);
   b_factor = DOT(u_ik, t_diff) / w_sin;

   for (i = 0; i < 3; ++i) {
      b_ijk[i] = t_ij[i] + b_factor * u_tb[i];
      r_ib[i] = b_ijk[i] - ai->pos[i];
   }

   h_ijk = (ai->rad + probe_rad) * (ai->rad + probe_rad) -
       (r_ib[0] * r_ib[0] + r_ib[1] * r_ib[1] + r_ib[2] * r_ib[2]);

   if (h_ijk <= 0.0) {

      /* h_ijk < 0.0 means either k is too far away to matter (=free torus between
         i and j) or k is so close that torus is buried by it.  if xtest < 0 torus
         is completely buried.  See Connolly, J. Appl. Cryst. 16, 548-558 (1983) 
         appendix I.  Need to check all k not just those with k>j>i Done in separate 
         routine */

      return 0;
   } else {                     /* pfound = 1; */
      h_ijk = sqrt(h_ijk);
      for (i = 0; i < 3; ++i) {
         p1_ijk[i] = b_ijk[i] + h_ijk * u_ijk[i];
         p2_ijk[i] = b_ijk[i] - h_ijk * u_ijk[i];
      }

      /* we know each of the tori is not free, so we can change the number of probes 
         from -1 (free) to 0.  -2 is completely buried, and those should have been 
         screened out in the get_probes calling routine... double check to be safe.
         btw: don't increment a non-negative nprobe here, because the probe positions 
         have not been checked for overlap with 4th atom */

      if (upper_neighbor_ij->nprobes == FREE_TORUS)
         upper_neighbor_ij->nprobes = 0;
      if (upper_neighbor_jk->nprobes == FREE_TORUS)
         upper_neighbor_jk->nprobes = 0;
      if (upper_neighbor_ik->nprobes == FREE_TORUS)
         upper_neighbor_ik->nprobes = 0;

      *height = h_ijk;
      return 1;
   }
}

/* returns 1 if one or both of the probe positions
   is not bumped by a 4th atom */

static int no_bump(ATOM atom[], int i0, int j0, int k0,
                   NEIGHBOR neighbors[], POINT px1, POINT px2,
                   REAL_T probe_rad, REAL_T probe_diam, int bump[])
{
   int ij, ia;
   REAL_T mindist2;

   /* 4th atom must be a neighbor of i0 (and j0 and k0)
      so just search 1 neighbor list */

   bump[0] = 0;
   bump[1] = 0;


   /* look over neighbor atoms */
   for (ij = atom[i0].neighbor_start;
        ij < atom[i0].neighbor_start + atom[i0].n_neighbors; ++ij) {
      ia = neighbors[ij].iatom; /* simplify notation */

      if (ia == j0 || ia == k0)
         continue;              /* must be unique */

      mindist2 = (probe_rad + atom[ia].rad) * (probe_rad + atom[ia].rad);
      if (!bump[0]
          && (atom[ia].pos[0] - px1[0]) * (atom[ia].pos[0] - px1[0]) +
          (atom[ia].pos[1] - px1[1]) * (atom[ia].pos[1] - px1[1]) +
          (atom[ia].pos[2] - px1[2]) * (atom[ia].pos[2] - px1[2])
          < mindist2)
         bump[0] = 1;

      if (!bump[1]
          && (atom[ia].pos[0] - px2[0]) * (atom[ia].pos[0] - px2[0]) +
          (atom[ia].pos[1] - px2[1]) * (atom[ia].pos[1] - px2[1]) +
          (atom[ia].pos[2] - px2[2]) * (atom[ia].pos[2] - px2[2])
          < mindist2)
         bump[1] = 1;

      if (bump[0] && bump[1])
         return 0;
   }





   return 1;
}

/* only check tori that are still marked "free"
   (i.e., nprobes = -1). all others have nprobes > 0
   nprobes = 0 (buried) or nprobes = -2 (buried).
   Here is where we determine whether those that are
   currently -1 (free) should be -2 (buried)

   This routine only checks those tori that are
   buried by a single atom.  The more common case
   is for a torus to be buried by several atoms.
   These tori are readily identified because
   they have t->nprobes = 0, but have no probes
   associated with them.
 */

static int t_buried(ATOM atom[], int *nat, NEIGHBOR neighbors[],
                    NEIGHBOR_TORUS upper_neighbors[], REAL_T probe_rad)
{
   REAL_T probe_diam;
   int ia, ij, nb = 0;

   probe_diam = 2 * probe_rad;

   for (ia = 0; ia < *nat; ++ia) {
      for (ij = atom[ia].upper_start;
           ij < atom[ia].upper_start + atom[ia].n_upper; ++ij) {

         if (upper_neighbors[ij].nprobes != -1)
            continue;           /* only check free */

         if (atom[ia].buried || atom[upper_neighbors[ij].iatom].buried) {
            /* torus containing these two atoms is buried 
               because one of the two atoms is buried  */
            upper_neighbors[ij].nprobes = -2;
            continue;
         }
         if (bury_check(ia, upper_neighbors[ij].iatom, atom,
                        probe_rad, probe_diam, neighbors)) {
            upper_neighbors[ij].nprobes = -2;
            ++nb;
         }
      }
   }

   return nb;
}

static int bury_check(int ia, int ja, ATOM atom[],
                      REAL_T probe_rad, REAL_T probe_diam,
                      NEIGHBOR neighbors[])
{
   int i1, i2, k1, k2;

   for (i1 = atom[ia].neighbor_start;
        i1 < atom[ia].neighbor_start + atom[ia].n_neighbors; ++i1) {
      k1 = neighbors[i1].iatom;
      for (i2 = atom[ja].neighbor_start;
           i2 < atom[ja].neighbor_start + atom[ja].n_neighbors; ++i2) {
         k2 = neighbors[i2].iatom;
         if (k1 == k2 && is_buried(&atom[ia], &atom[ja], &atom[k1],
                                   probe_rad, probe_diam))
            return 1;
      }
   }
   return 0;
}

static int is_buried(ATOM * ai, ATOM * aj, ATOM * ak,
                     REAL_T probe_rad, REAL_T probe_diam)
{
   POINT u_ij, u_ik, u_ijk, u_tb;
   REAL_T d_ij, d_ik;
   REAL_T t_rad, t_fact1, t_fact2;
   POINT t_ij, t_ik, t_diff;
   REAL_T b_factor, h_ijk, xtest;
   REAL_T w_ijk, w_sin, r_pk2;
   REAL_T sarg1, sarg2, check1, check2, check3, check4;
   POINT b_ijk, r_ib;
   int i;


   for (i = 0; i < 3; ++i) {
      u_ij[i] = aj->pos[i] - ai->pos[i];
      u_ik[i] = ak->pos[i] - ai->pos[i];
   }

   d_ij = sqrt(u_ij[0] * u_ij[0] + u_ij[1] * u_ij[1] + u_ij[2] * u_ij[2]);
   d_ik = sqrt(u_ik[0] * u_ik[0] + u_ik[1] * u_ik[1] + u_ik[2] * u_ik[2]);

   for (i = 0; i < 3; ++i) {
      u_ij[i] = u_ij[i] / d_ij;
      u_ik[i] = u_ik[i] / d_ik;
   }

   sarg1 = (ai->rad + aj->rad + probe_diam) * (ai->rad + aj->rad +
                                               probe_diam) - d_ij * d_ij;
   /* For rounding to 6 decimal points in case sarg1 < 0 */
   if(sarg1 < 0){
       check1 = sarg1;
       check2 = 0;
       check3 = fabs(modf( check1 * 1000000, &check2));
       check4 = (floor(check1 * 1000000) + (check3 < 0.5 ? 1 : 0 ))/1000000;
       sarg1 = check4;
   }
   sarg2 = d_ij * d_ij - (ai->rad - aj->rad) * (ai->rad - aj->rad);
   assert(sarg1 >= 0.0);
   assert(sarg2 >= 0.0);
   t_rad = 0.5 * sqrt(sarg1) * sqrt(sarg2) / d_ij;

   t_fact1 = ((ai->rad + probe_rad) * (ai->rad + probe_rad) -
              (aj->rad + probe_rad) * (aj->rad +
                                       probe_rad)) / (d_ij * d_ij);

   t_fact2 = ((ai->rad + probe_rad) * (ai->rad + probe_rad) -
              (ak->rad + probe_rad) * (ak->rad +
                                       probe_rad)) / (d_ik * d_ik);

   for (i = 0; i < 3; ++i) {
      t_ij[i] = 0.5 * (ai->pos[i] + aj->pos[i])
          + 0.5 * (aj->pos[i] - ai->pos[i]) * t_fact1;
      t_ik[i] = 0.5 * (ai->pos[i] + ak->pos[i])
          + 0.5 * (ak->pos[i] - ai->pos[i]) * t_fact2;
      t_diff[i] = t_ik[i] - t_ij[i];
   }

   w_ijk = acos(DOT(u_ij, u_ik));
   w_sin = sin(w_ijk);

   /* if w_sin is ~0 then the atoms are co-linear.
      In this case the distance between the probe sphere
      in any one of the three tori to the  3rd atom is
      constant.  If | r_pk | < r_p + probe_rad then
      the entire torus t_ij is buried.  If not the torus
      burial status remains unchanged.
      In any case, no (3-contact) probe position can be found,
      so 0 is returned.  */

   if (fabs(w_sin) < 1.0e-10) {
      r_pk2 = (ak->pos[0] - t_ij[0]) * (ak->pos[0] - t_ij[0]) +
          (ak->pos[1] - t_ij[1]) * (ak->pos[1] - t_ij[1]) +
          (ak->pos[2] - t_ij[2]) * (ak->pos[2] - t_ij[2]);

      if (r_pk2 + t_rad * t_rad -
          (ak->rad + probe_rad) * (ak->rad + probe_rad) < 0.0) {
         /*
            printf("sin test: torus between atoms %d and %d is buried by atom %d\n",ai->anum+1, aj->anum+1, ak->anum+1);
          */
         return 1;
      } else
         return 0;
   }
   cross(u_ij, u_ik, u_ijk);
   for (i = 0; i < 3; ++i)
      u_ijk[i] = u_ijk[i] / w_sin;
   cross(u_ijk, u_ij, u_tb);
   b_factor = DOT(u_ik, t_diff) / w_sin;

   for (i = 0; i < 3; ++i) {
      b_ijk[i] = t_ij[i] + b_factor * u_tb[i];
      r_ib[i] = b_ijk[i] - ai->pos[i];
   }

   h_ijk = (ai->rad + probe_rad) * (ai->rad + probe_rad) -
       (r_ib[0] * r_ib[0] + r_ib[1] * r_ib[1] + r_ib[2] * r_ib[2]);

   /* h_ijk < 0.0 means either k is too far away to matter (=free torus between i and j)
      or k is so close that torus is buried by it.
      if xtest < 0 torus is completely buried.
      See Connolly, J. Appl. Cryst. 16, 548-558 (1983) appendix I */
   if (h_ijk < 0.0) {
      xtest = (t_ij[0] - ak->pos[0]) * (t_ij[0] - ak->pos[0]) +
          (t_ij[1] - ak->pos[1]) * (t_ij[1] - ak->pos[1]) +
          (t_ij[2] - ak->pos[2]) * (t_ij[2] - ak->pos[2]) -
          ((ak->rad + probe_rad) * (ak->rad + probe_rad) - t_rad * t_rad);
      if (xtest < 0.0) {
         /*
            printf("x test: torus between atoms %d and %d is buried by atom %d\n",ai->anum+1, aj->anum+1, ak->anum+1);
          */
         return 1;
      }
   }
   return 0;
}

static int get_torus(ATOM atom[], int *nat, NEIGHBOR neighbors[],
                     NEIGHBOR_TORUS upper_neighbors[], REAL_T probe_rad,
                     TORUS toruslist[])
{
   int ia, ja, in, nt = 0;

   for (ia = 0; ia < *nat; ++ia) {
      atom[ia].torus_start = nt;        /* points to start of torus list for that atom */
      atom[ia].ntorus = 0;      /* points to start of torus list for that atom */
      for (in = atom[ia].upper_start;
           in < atom[ia].upper_start + atom[ia].n_upper; ++in) {
         if (upper_neighbors[in].nprobes == BURIED_TORUS)
            continue;
         if (upper_neighbors[in].nprobes == 0)
            continue;           /* no probes -> buried */

         ja = upper_neighbors[in].iatom;

         /* we have a new accessible torus */

         toruslist[nt].a1 = ia;
         ++atom[ia].ntorus;

         toruslist[nt].a2 = ja;

         torus_data(toruslist, nt, atom, ia, ja, probe_rad);

         if (upper_neighbors[in].nprobes == FREE_TORUS) {

            toruslist[nt].n_concave_edges = 0;
            toruslist[nt].n_convex_edges = -2;  /* taken care of in convex_edges() */
            /* printf("free torus %d\n", nt); */

         } else {               /* partial torus */

            toruslist[nt].n_concave_edges = 0;
            toruslist[nt].n_convex_edges = 0;
            /* these are taken care of in concave_edge() */

            /*
               toruslist[nt].n_concave_edges = upper_neighbors[in].nprobes;
               toruslist[nt].n_convex_edges  =  upper_neighbors[in].nprobes;

               if (toruslist[nt].n_convex_edges%2 != 0) {
               fprintf(stderr,"odd number of probe positions on torus!\n");
               exit (ERROR);
               }
             */

         }
         ++nt;
         if (nt >= NUM_TORUS * natm_sel) {
            fprintf(stderr, "MAX_TORUS exceeded\n");
            exit(ERROR);
         }
      }
   }
   return nt;
}


static void torus_data(TORUS tor[], int nt, ATOM atom[],
                       int ia, int ja, REAL_T prad)
{
   REAL_T d_ij, c_fact;
   REAL_T sarg1, sarg2;
   int ii;
   REAL_T ri, rj;

   ri = atom[ia].rad;
   rj = atom[ja].rad;

   d_ij = DIST(atom[ia].pos, atom[ja].pos);

   c_fact = ((ri + prad) * (ri + prad) -
             (rj + prad) * (rj + prad)) / (d_ij * d_ij);


   sarg1 = (ri + rj + 2 * prad) * (ri + rj + 2 * prad) - d_ij * d_ij;
   sarg2 = d_ij * d_ij - (ri - rj) * (ri - rj);
   assert(sarg1 >= 0.0);
   assert(sarg2 >= 0.0);
   tor[nt].rad = (0.5 / d_ij) * sqrt(sarg1) * sqrt(sarg2);

   if (tor[nt].rad < prad) {
      tor[nt].low = 1;
   } else {
      tor[nt].low = 0;
   }

   for (ii = 0; ii < 3; ++ii) {

      tor[nt].uv[ii] = (atom[ja].pos[ii] - atom[ia].pos[ii]) / d_ij;

      tor[nt].center[ii] = 0.5 * (atom[ia].pos[ii] + atom[ja].pos[ii]) +
          0.5 * c_fact * (atom[ja].pos[ii] - atom[ia].pos[ii]);
   }

   return;

}



/* bubble sort, based on example in
   Software Tools, Kernighan and Plauger */

static void sort_neighbors(int n, REAL_T x[], int istart,
                           NEIGHBOR neighbors[])
{
   int i, j, itmp;
   REAL_T xtmp;

   i = n - 1;                   /* index of last array element */

   while (i > 0) {
      j = 0;
      while (j < i) {
         if (x[j] > x[j + 1]) {

            xtmp = x[j];
            itmp = neighbors[istart + j].iatom;

            x[j] = x[j + 1];
            neighbors[istart + j].iatom = neighbors[istart + j + 1].iatom;

            x[j + 1] = xtmp;
            neighbors[istart + j + 1].iatom = itmp;
         }
         ++j;
      }
      i = i - 1;
   }
   return;
}

static int convex_circles(ATOM atom[], int nat, TORUS toruslist[],
                          int n_torus, CIRCLE circlelist[],
                          REAL_T probe_rad)
{


   int it, nc, ii, ia1, ia2;
   REAL_T trad;

   nc = 0;

   for (it = 0; it < n_torus; ++it) {
      ia1 = toruslist[it].a1;
      ia2 = toruslist[it].a2;
      trad = toruslist[it].rad;

      toruslist[it].circle1 = nc;
      circlelist[nc].torus = it;
      circlelist[nc].atom_or_probe_num = ia1;
      circlelist[nc].rad =
          (trad * atom[ia1].rad) / (atom[ia1].rad + probe_rad);

      for (ii = 0; ii < 3; ++ii) {
         circlelist[nc].axis[ii] = -toruslist[it].uv[ii];
         circlelist[nc].center[ii] =
             (atom[ia1].rad * toruslist[it].center[ii] +
              probe_rad * atom[ia1].pos[ii]) / (atom[ia1].rad + probe_rad);
      }

      ++nc;

      toruslist[it].circle2 = nc;
      circlelist[nc].torus = it;
      circlelist[nc].atom_or_probe_num = ia2;
      circlelist[nc].rad =
          (trad * atom[ia2].rad) / (atom[ia2].rad + probe_rad);

      for (ii = 0; ii < 3; ++ii) {
         circlelist[nc].axis[ii] = toruslist[it].uv[ii];
         circlelist[nc].center[ii] =
             (atom[ia2].rad * toruslist[it].center[ii] +
              probe_rad * atom[ia2].pos[ii]) / (atom[ia2].rad + probe_rad);
      }

      ++nc;

      if (nc >= NUM_CIRCLE * natm_sel) {
         fprintf(stderr, "MAX_CIRCLE exceeded\n");
         exit(ERROR);
      }
   }
   return nc;
}

/******************************************************************/
/* this routine assumes that the orientation of the atoms in the  */
/* probe triplet is counter-clockwise.  This should be assured by */
/* the add_probe() routine                                        */

static int concave_circles(ATOM atom[], int n_probes,
                           PROBE probelist[], TORUS toruslist[],
                           CIRCLE concave_circle_list[], REAL_T probe_rad)
{

   int ic, ip, ii, jj, it;
   POINT r_pt, vec;
   int p_torus[3], p_atom[3];
   REAL_T change_sign[3];

   ic = 0;
   for (ip = 0; ip < n_probes; ++ip) {

      p_atom[0] = probelist[ip].a1;
      p_atom[1] = probelist[ip].a2;
      p_atom[2] = probelist[ip].a3;

      p_torus[0] = id_torus(p_atom[0], p_atom[1], atom, toruslist);
      p_torus[1] = id_torus(p_atom[1], p_atom[2], atom, toruslist);
      p_torus[2] = id_torus(p_atom[0], p_atom[2], atom, toruslist);

      /* if the first atom in the torus is not the 1st atom in
         the concave surface edge, then the orientation of the
         concave circle is reversed */

      for (ii = 0; ii < 3; ++ii) {
         if (toruslist[p_torus[ii]].a1 != p_atom[ii])
            change_sign[ii] = -1;
         else
            change_sign[ii] = 1;
      }

      probelist[ip].c1 = ic;
      probelist[ip].c2 = ic + 1;
      probelist[ip].c3 = ic + 2;

      for (jj = 0; jj < 3; ++jj) {
         it = p_torus[jj];
         concave_circle_list[ic].atom_or_probe_num = ip;
         concave_circle_list[ic].rad = probe_rad;
         concave_circle_list[ic].torus = it;
         for (ii = 0; ii < 3; ++ii) {
            concave_circle_list[ic].center[ii] = probelist[ip].pos[ii];
            r_pt[ii] = toruslist[it].center[ii] - probelist[ip].pos[ii];
         }
         cross(r_pt, toruslist[it].uv, vec);
         vnorm(vec, 3);

         for (ii = 0; ii < 3; ++ii)
            concave_circle_list[ic].axis[ii] = change_sign[jj] * vec[ii];

         ++ic;
         if (ic >= NUM_CIRCLE * natm_sel) {
            printf("concave_circles() MAX_CIRCLE exceeded\n");
            exit(ERROR);
         }
      }
   }
   return ic;
}

/******************************************************************************/
  /* 1. fill up concave edge list, orient edges                                 */
  /* 2. add edge pointers to torus array                                        */
  /* 3. fill up vertex list                                                     */
  /* 4. create concave face                                                     */
  /*                                                                            */
  /* this routine assumes that the atoms associated with the probes are already */
  /* in counter-clockwise orientation.  This should be taken care of by         */
  /* add_probe()                                                                */
/******************************************************************************/

static void concave_edges(REAL_T probe_rad, ATOM atom[],
                          int nprobes, PROBE probe[], int *nverts,
                          VERTEX vert[], int *nedges, EDGE edge[],
                          int *nfaces, CONCAVE_FACE face[], int ntorus,
                          TORUS torus[])
{
   int ic1, ic2, ic3;           /* 3 circles on probe */
   int iv1, iv2;                /* vertices for the 3 atoms touching the probe */
   int it1, it2, it3;           /* 3 tori touching the probe */
   int iface = 0, i, ii, it, ie;
   int probe_atom[3];           /* 3 atoms touching probe */


   *nedges = 0;
   *nfaces = 0;
   *nverts = 0;

#ifdef DEBUG
   printf("concave edges\n");
#endif
   for (i = 0; i < nprobes; ++i) {
      probe_atom[0] = probe[i].a1;
      probe_atom[1] = probe[i].a2;
      probe_atom[2] = probe[i].a3;
      ic1 = probe[i].c1;
      ic2 = probe[i].c2;
      ic3 = probe[i].c3;

      it1 = id_torus(probe_atom[0], probe_atom[1], atom, torus);
      it2 = id_torus(probe_atom[1], probe_atom[2], atom, torus);
      it3 = id_torus(probe_atom[0], probe_atom[2], atom, torus);

      face[iface].probe = i;
      face[iface].alive = 1;

      /* define 3 new edges */

      edge[*nedges].vert1 = *nverts;
      edge[*nedges].vert2 = *nverts + 1;
      edge[*nedges].alive = 1;
      edge[*nedges].circle = ic1;
      face[iface].e1 = *nedges;
      torus[it1].concave_edges[torus[it1].n_concave_edges] = *nedges;
      ++(*nedges);

      edge[*nedges].vert1 = *nverts + 1;
      edge[*nedges].vert2 = *nverts + 2;
      edge[*nedges].circle = ic2;
      edge[*nedges].alive = 1;
      face[iface].e2 = *nedges;
      torus[it2].concave_edges[torus[it2].n_concave_edges] = *nedges;
      ++(*nedges);

      edge[*nedges].vert1 = *nverts + 2;
      edge[*nedges].vert2 = *nverts;
      edge[*nedges].alive = 1;
      edge[*nedges].circle = ic3;
      face[iface].e3 = *nedges;
      torus[it3].concave_edges[torus[it3].n_concave_edges] = *nedges;
      ++(*nedges);

      ++(torus[it1].n_concave_edges);
      ++(torus[it2].n_concave_edges);
      ++(torus[it3].n_concave_edges);

      if (torus[it1].n_concave_edges >= NUM_EDGE * natm_sel ||
          torus[it2].n_concave_edges >= NUM_EDGE * natm_sel ||
          torus[it3].n_concave_edges >= NUM_EDGE * natm_sel) {
         printf("MAXTOR_EDGE exceeded\n");
         exit(ERROR);
      }
      if (*nedges >= NUM_EDGE * natm_sel) {
         printf("MAX_EDGE exceeded\n");
         exit(ERROR);
      }
      /* define 3 new vertices */

      for (ii = 0; ii < 3; ++ii)
         addvert(probe_rad, probe_atom[ii], atom, i, probe, nverts, vert);
      ++iface;
      if (iface >= NUM_FACE * natm_sel) {
         printf("MAX_FACE exceeded\n");
         exit(ERROR);
      }
   }

   *nfaces = iface;

   for (it = 0; it < ntorus; ++it) {
      if (torus[it].n_concave_edges % 2 != 0) {
         fprintf(stderr, "odd number of probe positions on torus!\n");
         exit(ERROR);
      }
      for (i = 0; i < torus[it].n_concave_edges; ++i) {
         ie = torus[it].concave_edges[i];
         iv1 = edge[ie].vert1;
         iv2 = edge[ie].vert2;
         if ((vert[iv1].iatom != torus[it].a1 &&
              vert[iv1].iatom != torus[it].a2) ||
             (vert[iv2].iatom != torus[it].a1 &&
              vert[iv2].iatom != torus[it].a2)) {
            printf("concave edge on torus has mismatched atoms\n");
            printf("torus %d atoms %d %d\n", it, torus[it].a1,
                   torus[it].a2);
            printf("edge %d vert1.atom %d vert2.atom %d\n", ie,
                   vert[iv1].iatom, vert[iv2].iatom);
            exit(ERROR);
         }
      }
   }

   return;
}


static void addvert(REAL_T probe_rad, int ia, ATOM atom[],
                    int ip, PROBE probe[], int *nverts, VERTEX vert[])
{
   REAL_T invradsum;
   int ii;

   invradsum = 1.0 / (probe_rad + atom[ia].rad);
   vert[*nverts].iatom = ia;
   vert[*nverts].iprobe = ip;
   for (ii = 0; ii < 3; ++ii) {
      vert[*nverts].pos[ii] = invradsum *
          (atom[ia].rad * probe[ip].pos[ii] +
           probe_rad * atom[ia].pos[ii]);
   }
   ++(*nverts);
   if (*nverts >= NUM_VERTEX * natm_sel) {
      printf("MAX_VERTS exceeded\n");
      exit(ERROR);
   }
   return;
}

/******************************************************************************/
/* 1. use torus array concave edge pointers to fill up convex edge list       */
/* 2. fill up atom array convex edge pointers                                 */
/* 3. create saddle face                                                      */
/* orientation of saddle faces is as follows:                                 */
/*                                                                            */
/*     looking down on torus with torus axis ^                                */
/*                                                                            */
/*                    torus.atom2                                             */
/*                                                                            */
/*                   convex edge 2                                            */
/*                   ---->---                                                 */
/*                   |       |                                                */
/*           concave ^       v  concave edge 3                                */
/*           edge 1  |       |                                                */
/*                   ----<----                                                */
/*                  convex edge 4                                             */
/*                                                                            */
/*                    torus.atom1                                             */
/*                                                                            */
/******************************************************************************/


static void convex_edges(REAL_T probe_rad, int nat, ATOM atom[],
                         int n_probes, PROBE probelist[], int n_vertex,
                         VERTEX vertexlist[], int n_concave_edges,
                         EDGE concave_edge[], int *n_convex_edges,
                         EDGE convex_edge[], int n_torus,
                         TORUS toruslist[], int *n_saddle_faces,
                         SADDLE_FACE saddle_face[], CIRCLE circlelist[])
{
   int i, icc1, icc2, ia, ja, c1, c2, it, n_free_edges = 0;
   int cc1_v1, cc1_v2, cc2_v1, cc2_v2;
   int ne = 0;
   int iface = 0;

   *n_convex_edges = 0;
   *n_saddle_faces = 0;

   for (ia = 0; ia < nat; ++ia) {
      atom[ia].n_convex_edges = 0;
   }

   for (it = 0; it < n_torus; ++it) {
      ia = toruslist[it].a1;
      ja = toruslist[it].a2;
      c1 = toruslist[it].circle1;
      c2 = toruslist[it].circle2;
      if (toruslist[it].n_concave_edges == 0 && toruslist[it].n_convex_edges == -2) {   /* free torus */
         n_free_edges = n_free_edges + 2;
         toruslist[it].n_convex_edges = 0;

         add_saddle_face(saddle_face, &iface, it, -1, -1, ne + 1, ne);

         add_convex_edge_to_torus(it, toruslist, 0, ne);
         add_convex_edge_to_atom(ia, atom, ne);
         add_free_edge(&ne, convex_edge, c1);

         add_convex_edge_to_torus(it, toruslist, 1, ne);
         add_convex_edge_to_atom(ja, atom, ne);
         add_free_edge(&ne, convex_edge, c2);

         continue;
      }
      for (i = 0; i < toruslist[it].n_concave_edges; i = i + 2) {

         icc1 = toruslist[it].concave_edges[i];
         icc2 = toruslist[it].concave_edges[i + 1];

         cc1_v1 = concave_edge[icc1].vert1;
         cc1_v2 = concave_edge[icc1].vert2;
         cc2_v1 = concave_edge[icc2].vert1;
         cc2_v2 = concave_edge[icc2].vert2;

         atom_vertex_match(vertexlist, cc1_v1, ia);
         atom_vertex_match(vertexlist, cc1_v2, ja);
         atom_vertex_match(vertexlist, cc2_v1, ja);
         atom_vertex_match(vertexlist, cc2_v2, ia);

         add_saddle_face(saddle_face, &iface, it, icc1, icc2, ne + 1, ne);

         add_convex_edge_to_torus(it, toruslist, i, ne);
         add_convex_edge_to_atom(ia, atom, ne);
         add_edge(&ne, convex_edge,
                  concave_edge[icc2].vert2, concave_edge[icc1].vert1,
                  c1, vertexlist, circlelist);

         add_convex_edge_to_torus(it, toruslist, i + 1, ne);
         add_convex_edge_to_atom(ja, atom, ne);
         add_edge(&ne, convex_edge,
                  convex_edge[ne].vert1 = concave_edge[icc1].vert2,
                  convex_edge[ne].vert2 = concave_edge[icc2].vert1,
                  c2, vertexlist, circlelist);
      }
   }

   *n_convex_edges = ne;
   *n_saddle_faces = iface;

   if (ne - n_free_edges != n_concave_edges) {
      printf("convex_edges() convex and concave edges not paired\n");
      exit(ERROR);
   }
   check_data(n_torus, toruslist, nat, atom, vertexlist, convex_edge, ne);

   /*
      printf("number of convex edges %d\n", ne);
      printf("number of free   edges %d\n", n_free_edges);
      for (iface = 0; iface < *n_saddle_faces; ++iface)
      face_info(iface, saddle_face, concave_edge, convex_edge,vertexlist,toruslist);
      check_convex_edges(ne, convex_edge, circlelist);
    */

   return;
}

/********************************************************************************/
static void check_data(int n_torus, TORUS toruslist[],
                       int nat, ATOM atom[], VERTEX vertexlist[],
                       EDGE convex_edge[], int ne)
{
   int ntmp, it, i, ia, ie;

   ntmp = 0;
   for (it = 0; it < n_torus; ++it) {
      for (i = 0; i < toruslist[it].n_convex_edges; ++i) {
         ntmp = ntmp + 1;
      }
   }
   if (ntmp != ne) {
      printf("convex_edges(): number of convex edges from toruslist %d\n",
             ntmp);
      printf(" not equal to total number of convex edges %d\n", ne);
      exit(ERROR);
   }
   ntmp = 0;
   for (ia = 0; ia < nat; ++ia) {
      for (i = 0; i < atom[ia].n_convex_edges; ++i) {
         ntmp = ntmp + 1;
         ie = atom[ia].convex_edges[i];
         if (vertexlist[convex_edge[ie].vert1].iatom != ia &&
             convex_edge[ie].vert1 >= 0) {
            printf("convex vertex atom mismatch 2 \n");
            printf("iatom: %d\n", ia);
            printf("number of convex edges %d\n", atom[ia].n_convex_edges);
            printf("ie: %d\n", ie);
            printf("edge vertices: v1 %d v2 %d\n",
                   convex_edge[ie].vert1, convex_edge[ie].vert2);
            printf("edge v1.atom %d\n",
                   vertexlist[convex_edge[ie].vert1].iatom);
            printf("edge v2.atom %d\n",
                   vertexlist[convex_edge[ie].vert2].iatom);
            exit(ERROR);
         }
         if (vertexlist[convex_edge[ie].vert2].iatom != ia &&
             convex_edge[ie].vert1 >= 0) {
            printf("convex vertex atom mismatch 3\n");
            printf("iatom: %d\n", ia);
            printf("number of convex edges %d\n", atom[ia].n_convex_edges);
            printf("ie: %d\n", ie);
            printf("edge vertices: v1 %d v2 %d\n",
                   convex_edge[ie].vert1, convex_edge[ie].vert2);
            printf("edge v1.atom %d\n",
                   vertexlist[convex_edge[ie].vert1].iatom);
            printf("edge v2.atom %d\n",
                   vertexlist[convex_edge[ie].vert2].iatom);
            printf("convex vertex atom mismatch\n");
            exit(ERROR);
         }
      }
   }
   if (ntmp != ne) {
      printf("convex_edges(): number of convex edges from atomlist %d\n",
             ntmp);
      printf(" not equal to total number of convex edges %d\n", ne);
      exit(ERROR);
   }
   return;
}

/********************************************************************************/


static void add_convex_edge_to_torus(int it, TORUS toruslist[],
                                     int iedge, int ne)
{
   toruslist[it].convex_edges[iedge] = ne;
   ++toruslist[it].n_convex_edges;
   if (toruslist[it].n_convex_edges > NUM_EDGE * natm_sel) {
      printf("convex_edges() MAXTOR_EDGE exceeded\n");
      exit(ERROR);
   }
   return;
}

static void add_convex_edge_to_atom(int ia, ATOM atom[], int ne)
{
   atom[ia].convex_edges[atom[ia].n_convex_edges] = ne;
   ++atom[ia].n_convex_edges;
   if (atom[ia].n_convex_edges >= NUM_EDGE * natm_sel) {
      printf("convex_edges() MAXAT_EDGE exceeded\n");
      exit(ERROR);
   }
   return;
}

static void convex_faces(int nat, ATOM atom[],
                         int *nface, CONVEX_FACE face[], int n_cycles,
                         CYCLE cycle[], EDGE edge[], CIRCLE circle[],
                         int nverts, VERTEX vert[])
{

   int ia, ic, i, j, jc, k, kc, add_to_face, nfc;
   int iface = 0;
   int inside[MAXAT_CYCLES][MAXAT_CYCLES];
   int in_face[MAXAT_CYCLES];   /* has face number for cycle, -1 if none set */

   /*  inside[i][j] = 1 if cycle i is inside cycle j */

   for (ia = 0; ia < nat; ++ia) {

      if (atom[ia].n_cycles == 0 && atom[ia].n_neighbors == 0) {
         face[iface].n_cycles = 0;
         face[iface].atom = ia;
         /* printf("free atom\n"); */
         ++iface;
      } else if (atom[ia].n_cycles == 1) {
         /* printf("atom %d has %d cycle\n", ia, atom[ia].n_cycles); */
         face[iface].n_cycles = 1;
         face[iface].cycle[0] = atom[ia].cycle_start;
         face[iface].atom = ia;
         ++iface;
      } else if (atom[ia].n_cycles > 1) {
         /* printf("atom %d has %d cycles\n", ia, atom[ia].n_cycles); */

         for (i = 0; i < atom[ia].n_cycles; ++i) {
            ic = atom[ia].cycle_start + i;
            in_face[i] = -1;
            for (j = 0; j < atom[ia].n_cycles; ++j) {
               jc = atom[ia].cycle_start + j;
               inside[i][j] = is_cycle_inside(ic, jc,
                                              cycle, edge, circle, vert,
                                              atom);
               /* printf("inside[%d][%d] = %d\n", i, j, inside[i][j]); */
            }
         }

         for (i = 0; i < atom[ia].n_cycles; ++i) {
            ic = atom[ia].cycle_start + i;
            if (in_face[i] == -1) {     /* not part of a face yet */
               nfc = 0;         /* initialize cycle count */
               in_face[i] = iface;      /* mark as part of a face */
               face[iface].atom = ia;
               face[iface].cycle[nfc] = ic;
               ++nfc;

               /* find a cycle jc that is inside ic and that ic is inside of */

               for (j = i + 1; j < atom[ia].n_cycles; ++j) {
                  jc = atom[ia].cycle_start + j;

                  /* jc must be in no face */
                  if (in_face[j] != -1)
                     continue;

                  /* ic and jc must be inside each other */
                  if (!inside[i][j] || !inside[j][i])
                     continue;

                  k = 0;
                  add_to_face = 1;      /* assume jc can be added */

                  while (add_to_face && k < atom[ia].n_cycles) {
                     kc = atom[ia].cycle_start + k;
                     if ((kc != ic && kc != jc) &&      /* must be unique */
                         (inside[k][i] && inside[k][j]) &&      /* must be inside both ic and jc */
                         (!inside[i][k] || !inside[j][k])       /* ic and jc must be inside kc */
                         )
                        add_to_face = 0;
                     ++k;
                  }
                  if (add_to_face) {
                     in_face[j] = iface;
                     face[iface].cycle[nfc] = jc;
                     ++nfc;
                  }
               }
               face[iface].n_cycles = nfc;
               ++iface;
            }
         }

      } else {
         /*
            printf("buried atom %d\n", ia);
          */
      }
   }
   *nface = iface;

   /*
      printf("number of faces %d\n", *nface);
      for (iface = 0; iface < *nface; ++iface) {
      printf("---------------\n");
      printf("face %d on atom %d has %d cycles\n", iface, face[iface].atom, face[iface].n_cycles);
      for (i = 0; i < face[iface].n_cycles; ++i) {
      ic = face[iface].cycle[i];
      printf("  cycle %d (atom %d)\n", ic, cycle[ic].atom);
      for (j = 0; j < cycle[ic].nedges; ++j) {
      ie = cycle[ic].edge[j];
      printf("     edge verts %d %d  (atom %d)\n", 
      edge[ie].vert1, edge[ie].vert2, edge[ie].atom);
      }
      }
      }
    */

   return;

}

/************************************************************/
/* sort_edges() sorts the concave edges associated with     */
/* a torus such that the edges are paired by the saddle     */
/* surface that joins them.  I'm using the vertex on the    */
/* 1st atom of the torus to do the sorting.                 */
/************************************************************/
static void sort_edges(int n_concave_edges, EDGE concave_edge_list[],
                       int n_torus, TORUS toruslist[], int n_vertex,
                       VERTEX vertexlist[], CIRCLE circlelist[])
{
   int it, ie, iv;

   /* vsort[]:  for each edge there is a vertex associated with atom 1
      in the torus.  Vectors are constructed from the atom 1 circle center
      to each of these edges, and the angle between the vector to the 1st vertex
      and the subsequent vertex is the value we sort on.  The 1st vertex is
      arbitrary, so the pairing may be off by one (i.e., the sorting is correct,
      but the 1st vertex should be pair with the last vertex, not the 2nd).
      clear as mud, right? */

   REAL_T vsort[MAXTOR_EDGE];
   int ntmp;
   int iatom, jatom, iedge, itmp, icircle;
   /* u[3] is the vector from the circle center to the first vertex    */
   /* v[3] is the vector from the circle center each successive vertex */
   POINT u, v;
   int i, j, n;
   REAL_T vtmp, twopi;

   twopi = 2.0 * acos(-1.0);

   for (it = 0; it < n_torus; ++it) {

      ntmp = toruslist[it].n_concave_edges;
      if (ntmp < 1)
         continue;
      iatom = toruslist[it].a1;
      jatom = toruslist[it].a2;
      icircle = toruslist[it].circle1;

      if (ntmp >= MAXTOR_EDGE) {
         fprintf(stderr, "sort_edges() MAXTOR_EDGE exceeded\n");
         exit(ERROR);
      }
/****************************************************************/
      /* identify the vertices that we are using to sort the edges    */
      /* and calculate angles between the vectors from the circle     */
      /* center to the vertices                                       */
/****************************************************************/

      /* separate 1st one */

      iedge = toruslist[it].concave_edges[0];

      iv = concave_edge_list[iedge].vert1;
      if (vertexlist[iv].iatom != iatom)
         iv = concave_edge_list[iedge].vert2;
      if (vertexlist[iv].iatom != iatom) {
         fprintf(stderr, "sort_edges(): iatom not found\n");
         exit(ERROR);
      }
      u[0] = vertexlist[iv].pos[0] - circlelist[icircle].center[0];
      u[1] = vertexlist[iv].pos[1] - circlelist[icircle].center[1];
      u[2] = vertexlist[iv].pos[2] - circlelist[icircle].center[2];

      vsort[0] = 0.0;

      /* calculate the angles between the 1st vertex and the others     */

      for (itmp = 1; itmp < ntmp; ++itmp) {
         iedge = toruslist[it].concave_edges[itmp];
         iv = concave_edge_list[iedge].vert1;
         if (vertexlist[iv].iatom != iatom)
            iv = concave_edge_list[iedge].vert2;
         if (vertexlist[iv].iatom != iatom) {
            fprintf(stderr, "sort_edges(): iatom not found\n");
            exit(ERROR);
         }
         v[0] = vertexlist[iv].pos[0] - circlelist[icircle].center[0];
         v[1] = vertexlist[iv].pos[1] - circlelist[icircle].center[1];
         v[2] = vertexlist[iv].pos[2] - circlelist[icircle].center[2];
         /* sort by clockwise angle from u */
         vsort[itmp] = get_angle(u, v, circlelist[icircle].axis);
         if (vsort[itmp] < 0.0)
            vsort[itmp] = twopi + vsort[itmp];
      }

/*****************************************************************/
      /* sort the edges based on the angles                            */
/*****************************************************************/

      n = toruslist[it].n_concave_edges;
      i = n - 1;

      while (i > 0) {
         j = 0;
         while (j < i) {
            if (vsort[j] > vsort[j + 1]) {
               vtmp = vsort[j];
               itmp = toruslist[it].concave_edges[j];

               vsort[j] = vsort[j + 1];
               toruslist[it].concave_edges[j] =
                   toruslist[it].concave_edges[j + 1];

               vsort[j + 1] = vtmp;
               toruslist[it].concave_edges[j + 1] = itmp;

            }
            ++j;
         }
         i = i - 1;
      }

      /* check for the "off by one error" (i.e., 1st edge may not be paired with
         2nd edge, but rather the last edge */

      ie = toruslist[it].concave_edges[0];
      iv = concave_edge_list[ie].vert1;
      if (vertexlist[iv].iatom != iatom) {
         itmp = toruslist[it].concave_edges[0];
         n = toruslist[it].n_concave_edges;
         for (iedge = 0; iedge < n - 1; ++iedge) {
            toruslist[it].concave_edges[iedge] =
                toruslist[it].concave_edges[iedge + 1];
         }
         toruslist[it].concave_edges[n - 1] = itmp;
      }
      /*
         counter = 0;
         printf("edges associated with torus %d between atoms %d and %d\n", it, toruslist[it].a1,
         toruslist[it].a2);
         for (i = 0; i < toruslist[it].n_concave_edges; ++i) {
         ie = toruslist[it].concave_edges[i];
         printf("edge: %d vertex atoms: %d %d\n",ie, 
         vertexlist[concave_edge_list[ie].vert1].iatom,  
         vertexlist[concave_edge_list[ie].vert2].iatom); 

         iv1 = concave_edge_list[ie].vert1;
         iv2 = concave_edge_list[ie].vert2;
         ++counter;
         printf("ATOM  %5d O    PRO  1001    %8.3f%8.3f%8.3f\n", counter,
         vertexlist[iv1].pos[0], vertexlist[iv1].pos[1], vertexlist[iv1].pos[2]);
         ++counter;
         printf("ATOM  %5d C    PRO  1001    %8.3f%8.3f%8.3f\n", counter,
         vertexlist[iv2].pos[0], vertexlist[iv2].pos[1], vertexlist[iv2].pos[2]);
         printf("BND %d %d\n", counter-1, counter);
         }
       */
   }
   return;
}

static void cycles(int nat, ATOM atom[], VERTEX vertexlist[],
                   int n_convex_edges, EDGE convex_edge_list[],
                   int *n_cycles, CYCLE clist[], CIRCLE circlelist[],
                   TORUS toruslist[])
{
   int atom_edge_cycle[MAXAT_EDGE];
   /* atom_edge_cycle[] is a temporary array
      to keep track of which cycles the edges of
      an atom belong to.  Initialized to -1 for
      each element, as each edge is added to a face
      it is change to the face index.  When
      all elements are non-negative, you're finished
      making cycles for that atom */
   int ice;                     /* edge list index for a cycle */
   int nc;                      /* index for cyclelist */

   int ia, i, ie, je, last_vert, second_vert;

   nc = 0;

   for (ia = 0; ia < nat; ++ia) {
      atom[ia].n_cycles = 0;
      atom[ia].cycle_start = nc;
      if (atom[ia].n_convex_edges == 0)
         continue;

/************************************************/
      /* initialize atom_edge_cycle elements to -1    */
/************************************************/
      for (i = 0; i < atom[ia].n_convex_edges; ++i)
         atom_edge_cycle[i] = -1;

      /*
         printf("atom %d has %d convex edges\n", ia, atom[ia].n_convex_edges);
         for (i = 0; i < atom[ia].n_convex_edges; ++i) {
         ie = atom[ia].convex_edges[i];
         printf("edge %d v1: %d v2: %d\n", ie, convex_edge_list[ie].vert1, 
         convex_edge_list[ie].vert2);
         }
       */

      while ((ie = new_edge(ia, atom,
                            atom[ia].n_convex_edges,
                            atom_edge_cycle, nc)) != -1) {
         /* we have found an edge not on a cycle */
         ice = 0;               /* initialize cycle index               */

         /* printf("new_edge: %d\n", ie); */

         clist[nc].edge[ice] = ie;      /* new edge is first edge in new cycle  */
         clist[nc].nedges = 1;  /* initialize edge count to 1           */
         clist[nc].atom = ia;   /* associate cycle with atom ia         */
         ++ice;
         if (convex_edge_list[ie].vert1 == -1) {        /* a free edge */
            ++nc;               /* finsihed with the cycle already */
            if (nc >= NUM_CYCLE * natm_sel) {
               fprintf(stderr, "MAX_CYCLES exceeded\n");
               exit(ERROR);
            }
            atom[ia].n_cycles++;        /* increment atom cycle count */
            continue;           /* move on to next new_edge   */
         }
         /* once the 2nd vertex in an edge is equal to the 1st vertex in the 1st
            edge, you have finished a cycle.  first need to id the vertex */

         last_vert = convex_edge_list[ie].vert1;
         second_vert = convex_edge_list[ie].vert2;

         while (second_vert != last_vert) {
            je = next_edge(second_vert, ia, atom,
                           atom[ia].n_convex_edges, convex_edge_list,
                           atom_edge_cycle, nc);
            clist[nc].edge[ice] = je;
            clist[nc].nedges++;
            ++ice;
            second_vert = convex_edge_list[je].vert2;
         }
         ++nc;                  /* when second_vert = last_vert you've finished the cycle */
         if (nc >= NUM_CYCLE * natm_sel) {
            fprintf(stderr, "MAX_CYCLES exceeded\n");
            exit(ERROR);
         }
         atom[ia].n_cycles++;   /* increment atom cycle count */


      }
   }
   *n_cycles = nc;

   return;
}


static int new_edge(int ia, ATOM atom[], int nae, int atom_edge_cycle[],
                    int ncycle)
{
   int iae, ie;

   for (iae = 0; iae < nae; ++iae) {
      if (atom_edge_cycle[iae] == -1) {
         atom_edge_cycle[iae] = ncycle;
         ie = atom[ia].convex_edges[iae];
         /* printf("new_edge() found edge %d\n", ie); */
         return ie;
      }
   }
   return -1;                   /* no new edges found */
}

static int next_edge(int ivert, int ia, ATOM atom[],
                     int nae, EDGE convex_edges[], int atom_edge_cycle[],
                     int ncycle)
{

   int iae, ie;

   for (iae = 0; iae < nae; ++iae) {
      ie = atom[ia].convex_edges[iae];
      if (atom_edge_cycle[iae] == -1 && convex_edges[ie].vert1 == ivert) {
         atom_edge_cycle[iae] = ncycle;
         return ie;
      }
   }
   printf("next_edge(): No edge found\n");
   dump_atom_edges(ivert, ia, atom, nae, convex_edges, atom_edge_cycle);

   exit(ERROR);
   return -1;
}

static void dump_atom_edges(int ivert, int ia, ATOM atom[],
                            int nae, EDGE convex_edges[],
                            int atom_edge_cycle[])
{
   int iae, ie;


   printf("looking for 1st vertex of : %d\n", ivert);
   printf("amoung edges for atom %d\n", ia);

   for (iae = 0; iae < nae; ++iae) {
      ie = atom[ia].convex_edges[iae];
      printf("edge: %d  v1: %d v2: %d\n",
             ie, convex_edges[ie].vert1, convex_edges[ie].vert2);
   }
   return;
}

#define MKAXIS(ax,x,y,z)        ((ax)[0]=(x),(ax)[1]=(y),(ax)[2]=(z))

static void draw_circle(int *narc, int *npt, int n,
                        FILE * atomfile, FILE * bondfile, FILE * propfile,
                        POINT center, REAL_T rad, POINT zaxis, int resnum)
{
   REAL_T aij[3][3];
   REAL_T dx, dy, dz;
   REAL_T xax[3], yax[3], zax[3];
   REAL_T theta, dtheta;
   REAL_T xc[3], xpoint[3];
   int i, ii;


   (*narc)++;
   dtheta = (2.0 * acos(-1.0)) / n;

   for (ii = 0; ii < 3; ++ii)
      zax[ii] = zaxis[ii];

   vnorm(zax, 3);

/* generate arbitrary x and y ax: borrowed from
   Tom Macke's ax code in NAB */
   dx = zax[0];
   dy = zax[1];
   dz = zax[2];

   if (dz) {                    /* z value non-zero */

      MKAXIS(yax, 0.0, 1.0, -dy / dz);
      vnorm(yax, 3);
      cross(yax, zax, xax);
      vnorm(xax, 3);

   } else if (dy) {

      MKAXIS(yax, 1.0, -dx / dy, 0.0);
      vnorm(yax, 3);
      cross(yax, zax, xax);
      vnorm(xax, 3);

   } else if (dx) {

      MKAXIS(yax, 0.0, 1.0, 0.0);
      MKAXIS(xax, 0.0, 0.0, -1.0);

   } else {
      MKAXIS(xax, 1.0, 0.0, 0.0);
      MKAXIS(yax, 0.0, 1.0, 0.0);
      MKAXIS(zax, 0.0, 0.0, 1.0);
   }

   aij[0][0] = xax[0];
   aij[0][1] = xax[1];
   aij[0][2] = xax[2];

   aij[1][0] = yax[0];
   aij[1][1] = yax[1];
   aij[1][2] = yax[2];

   aij[2][0] = zax[0];
   aij[2][1] = zax[1];
   aij[2][2] = zax[2];


   xc[2] = 0.0;

   for (i = 0; i <= n; ++i) {
      theta = i * dtheta;
      xc[0] = rad * cos(theta);
      xc[1] = rad * sin(theta);

      for (ii = 0; ii < 3; ++ii) {
         xpoint[ii] =
             center[ii] + aij[0][ii] * xc[0] + aij[1][ii] * xc[1] +
             aij[2][ii] * xc[2];
      }
      fprintf(atomfile, "ATOM      1 N    VRT  %4d    %8.3f%8.3f%8.3f\n",
              resnum, xpoint[0], xpoint[1], xpoint[2]);
      fprintf(propfile, "0.5\n");
      (*npt)++;
      if (i > 0)
         fprintf(bondfile, "%d %d\n", *npt - 1, *npt);
   }
}


#define NPOINTS 20

/**************************************************************/
/* is_cycle_inside()  determines wheter cycle ic is           */
/* inside cycle jc.  It does this by taking a point on        */
/* cycle ic  and using it to project cycle jc onto the plane  */
/* that is tangent to the antipode of ic.  If the orientation */
/* of the polygon is counter-clockwise, the point is inside   */
/* and cycle ic is contained by jc.  Note if cycle ic is a    */
/* circle use the center of the circle projected out to the   */
/* sphere surface along the radius containing that point      */
/*                                                            */
/* cycles appear counterclockwise when viewed from above      */
/* their interior                                             */
/*                                                            */
/* P is the North Pole, Q is the South Pole u_pq is the unit  */
/* vector pointing from North to South                        */
/**************************************************************/

static int is_cycle_inside(int icycle, int jcycle, CYCLE cycle[],
                           EDGE edge[], CIRCLE circle[], VERTEX vert[],
                           ATOM atom[])
{
   int ic, ia, iv, i, j, ii, ie, je, ne;
   POINT p, u_qp, u_pq, u_pr, uvec, vvec;
   POINT pt[MAXAT_EDGE];
   REAL_T dist, sdist, theta, theta_tot, dpq;

   if (icycle == jcycle)
      return 0;

/****************************************************/
   /*  special cases                                   */
/****************************************************/
   if (cycle[jcycle].nedges < 3) {
      return 1;
   }
   for (i = 0; i < cycle[icycle].nedges; ++i) {
      ie = cycle[icycle].edge[i];
      for (j = 0; j < cycle[jcycle].nedges; ++j) {
         je = cycle[jcycle].edge[j];
         if (edge[ie].circle == edge[je].circle) {
            return 0;
         }
      }
   }


/****************************************************/
   /* first determine point P, the north pole of the   */
   /*   stereographic projection                       */
/****************************************************/

   ie = cycle[icycle].edge[0];
   ia = cycle[icycle].atom;

   if (cycle[icycle].nedges == 1) {
      ic = edge[ie].circle;
      /* you have a circle: set p to center of exterior of circle */

      p[0] = atom[ia].pos[0] - atom[ia].rad * circle[ic].axis[0];
      p[1] = atom[ia].pos[1] - atom[ia].rad * circle[ic].axis[1];
      p[2] = atom[ia].pos[2] - atom[ia].rad * circle[ic].axis[2];

   } else {
      /* you have a cycle: set p to 1st vertex */
      iv = edge[ie].vert1;
      p[0] = vert[iv].pos[0];
      p[1] = vert[iv].pos[1];
      p[2] = vert[iv].pos[2];
   }

   if (cycle[jcycle].nedges > MAXAT_EDGE) {
      printf("is_cycle_inside(): MAXAT_EDGE exceeded\n");
      exit(ERROR);
   }
   for (ii = 0; ii < 3; ++ii)
      u_pq[ii] = atom[ia].pos[ii] - p[ii];
   for (ii = 0; ii < 3; ++ii)
      u_qp[ii] = p[ii] - atom[ia].pos[ii];
   vnorm(u_qp, 3);
   vnorm(u_pq, 3);


   ne = cycle[jcycle].nedges;


   for (i = 0; i < ne; ++i) {
      ie = cycle[jcycle].edge[i];
      iv = edge[ie].vert1;

      u_pr[0] = vert[iv].pos[0] - p[0];
      u_pr[1] = vert[iv].pos[1] - p[1];
      u_pr[2] = vert[iv].pos[2] - p[2];

      dist = sqrt(u_pr[0] * u_pr[0] +
                  u_pr[1] * u_pr[1] + u_pr[2] * u_pr[2]);

      u_pr[0] = u_pr[0] / dist;
      u_pr[1] = u_pr[1] / dist;
      u_pr[2] = u_pr[2] / dist;


      dpq = 2.0 * atom[ia].rad;
      sdist = dpq / (DOT(u_pq, u_pr));

      if (sdist < 0) {
         printf("is_cycle_inside(): sdist < 0\n");
         exit(ERROR);
      }
      for (ii = 0; ii < 3; ++ii)
         pt[i][ii] = p[ii] + sdist * u_pr[ii];
   }

   theta_tot = 0.0;
   for (i = 1; i < ne - 1; ++i) {
      for (ii = 0; ii < 3; ++ii) {
         uvec[ii] = pt[i][ii] - pt[i - 1][ii];
         vvec[ii] = pt[i + 1][ii] - pt[i][ii];
      }
      theta = get_angle(uvec, vvec, u_qp);
      theta_tot += theta;
   }

   for (ii = 0; ii < 3; ++ii) {
      uvec[ii] = pt[ne - 1][ii] - pt[ne - 2][ii];
      vvec[ii] = pt[0][ii] - pt[ne - 1][ii];
   }
   theta = get_angle(uvec, vvec, u_qp);
   theta_tot += theta;

   for (ii = 0; ii < 3; ++ii) {
      uvec[ii] = pt[0][ii] - pt[ne - 1][ii];
      vvec[ii] = pt[1][1] - pt[0][ii];
   }
   theta = get_angle(uvec, vvec, u_qp);
   theta_tot += theta;

   if (theta_tot < 0) {
      return 1;                 /* inside */
   } else {
      return 0;
   }


}

/************************************************************************/
/* get_angle() returns the angle between two vectors u and v            */
/* as measured looking down from the +z axis                            */
/* from -pi to pi  measured clockwise                                   */
/*                                                                      */
/* relying on calling routine to be suze z in perp to u v plane         */
/*                                                                      */
/*                                                                      */
/*             u(Y)                                                     */
/*         v2                                                           */
/*         ^    ^    ^ v1                                               */
/*          \   |   /                                                   */
/*           \  |  /                                                    */
/*            \ | /                                                     */
/*              o ------> (X)                                           */
/*           / (Z)\                                                     */
/*          /      \                                                    */
/*         /        \                                                   */
/*     v3 v          v v4                                               */
/*                                                                      */
/*                                                                      */
/*  set up u to lie along +y axis  atan2 returns arctan between -pi     */
/*  ang pi with the sign determined by the sign of the first arguemnt   */
/* in this case v_y  so call atan2(v_x,v_y)                             */
/*   so       atan(u,v4) > atan(u,v1) > 0                               */
/*   so       atan(u,v3) < atan(u,v2) < 0                               */
/*                                                                      */
/************************************************************************/

static REAL_T get_angle(POINT u, POINT v, POINT zaxis)
{
   REAL_T x, y;
   POINT xaxis, yaxis;

   yaxis[0] = u[0];
   yaxis[1] = u[1];
   yaxis[2] = u[2];
   vnorm(yaxis, 3);

   cross(yaxis, zaxis, xaxis);
   vnorm(xaxis, 3);

   x = DOT(v, xaxis);
   y = DOT(v, yaxis);

   return atan2(x, y);
}

static REAL_T convex_area(ATOM atom[], RES res[], CYCLE cycle[],
                          int n_faces, CONVEX_FACE face[],
                          EDGE convex_edge[], CIRCLE circle[],
                          VERTEX vert[], TORUS torus[])
{
   int iface, ic, icycle, chi, ia;
   REAL_T arad, area, total_area = 0;

   for (iface = 0; iface < n_faces; ++iface) {

      chi = 2 - face[iface].n_cycles;

      area = 0;
      ia = face[iface].atom;
      arad = atom[ia].rad;

      /*
         printf("max area of spherical atom %f\n", 4.0*PI*arad*arad);
       */

      for (ic = 0; ic < face[iface].n_cycles; ++ic) {
         icycle = face[iface].cycle[ic];
         area =
             area + cycle_piece(icycle, cycle, circle, convex_edge, vert,
                                atom, torus);
      }

      face[iface].area = arad * arad * (2 * PI * chi + area);

      total_area += face[iface].area;
      atom[ia].area += face[iface].area;
#ifdef DEBUG
      printf("convex face atom %-4s %-4s %4d area: %10.3f\n",
             atom[ia].anam,
             res[atom[ia].res].nam,
             res[atom[ia].res].num, face[iface].area);

      printf("area of face %d = %f\n\n\n", iface, face[iface].area);
#endif

   }
   return total_area;
}



static REAL_T cycle_piece(int ic, CYCLE cycle[], CIRCLE circle[],
                          EDGE convex_edge[], VERTEX vert[], ATOM atom[],
                          TORUS torus[])
{
   REAL_T sum, wrap_angle;
   int i, ii, ie, icircle, iatom, itorus, iv1, iv2;
   POINT v1, v2, d_c;
   REAL_T cos_theta;
   sum = 0;

   /*
      printf("   cycle_piece()\n");
      printf("   number of edges %d\n", cycle[ic].nedges);
    */

   for (i = 0; i < cycle[ic].nedges; ++i) {
      ie = cycle[ic].edge[i];
      icircle = convex_edge[ie].circle;
      iatom = circle[icircle].atom_or_probe_num;
      itorus = circle[icircle].torus;
      iv1 = convex_edge[ie].vert1;
      iv2 = convex_edge[ie].vert2;
      /* write_info(ic,ie,icircle,iatom,itorus,iv1,iv2); */

/********************************************/
      /* 1st take care of the interior angle part */

      if (iv1 > -1)
         sum -= vert[iv1].beta;
      /* printf("beta: %f\n", Rad2Deg*vert[iv1].beta); */

      /* next do the saddle wrap stuff. This part is a little confusing and 
         Connolly's paper (J. Appl. Cryst.  16, 548-558 (1983)) is unclear about 
         the sign of theta.  I think he assumes the torus center always lies between 
         atoms, but this is not necessarily so.  You basically want the projection 
         of the vector that points from the atom center to the circle (or vertex, 
         either way) onto the interatomic axis), but the sign of that term is given 
         by the opposite dot product of the d_c vector (atom to circle center) and 
         the circle unit vector.  sign = -dot(d_c, c_axis) if the d_c and c_axis are 
         in opposite directions, then the torus center is indeed between the two atom 
         and this term should be added.  If not, the term should be subtracted.  
         the MSEED paper (Perrot, et al.  J Comp Chem 13, 1-11 (1992)) has it right.
       */

      if (iv1 == -1) {
         if (cycle[ic].nedges != 1) {
            printf("cycle_area(): vert = -1 but n_edges > 1\n");
            exit(ERROR);
         }
         wrap_angle = 2.0 * PI;
      } else {
         for (ii = 0; ii < 3; ++ii) {
            v1[ii] = vert[iv2].pos[ii] - circle[icircle].center[ii];
            v2[ii] = vert[iv1].pos[ii] - circle[icircle].center[ii];
         }
         wrap_angle = get_angle(v1, v2, circle[icircle].axis);
         if (wrap_angle < 0.0)
            wrap_angle += 2.0 * PI;
      }

      /* printf("wrap angle %f\n",Rad2Deg*wrap_angle); 
         printf("circle center %8.3f%8.3f%8.3f\n", circle[icircle].center[0], 
         circle[icircle].center[1], circle[icircle].center[2]);
         printf("circle axis %8.3f%8.3f%8.3f\n", circle[icircle].axis[0],
         circle[icircle].axis[1], circle[icircle].axis[2]); */

      d_c[0] = circle[icircle].center[0] - atom[iatom].pos[0];
      d_c[1] = circle[icircle].center[1] - atom[iatom].pos[1];
      d_c[2] = circle[icircle].center[2] - atom[iatom].pos[2];
      cos_theta = -(d_c[0] * circle[icircle].axis[0] +
                    d_c[1] * circle[icircle].axis[1] +
                    d_c[2] * circle[icircle].axis[2]) / atom[iatom].rad;
      sum += wrap_angle * cos_theta;
      /* printf("d_c %8.3f%8.3f%8.3f\n", d_c[0], d_c[1], d_c[2]); 
         printf("cos theta %f\n", cos_theta);
         printf("    theta %f\n", Rad2Deg*acos(cos_theta));
         printf("  contribution from wrap angle term %f\n", wrap_angle*cos_theta); */
   }
   /* printf("contribution from edge %f\n", sum); */
   return sum;
}



static int concave_area(REAL_T probe_rad, int n_c_faces, VERTEX vertex[],
                        CONCAVE_FACE c_face[], EDGE c_edge[],
                        CIRCLE c_circle[], REAL_T * c_area)
{
   int iface, ie1, ie2, ie3, ii, iv, jv, kv;
   POINT n_ij, n_jk, n_ik, n_ji, n_kj, n_ki;
   REAL_T total_area = 0;
   int ic1, ic2, ic3, n_surfaced = 0;

   *c_area = 0.0;
   for (iface = 0; iface < n_c_faces; ++iface) {
      ++n_surfaced;

      ie1 = c_face[iface].e1;
      ie2 = c_face[iface].e2;
      ie3 = c_face[iface].e3;

      iv = c_edge[ie1].vert1;
      jv = c_edge[ie2].vert1;
      kv = c_edge[ie3].vert1;

      ic1 = c_edge[ie1].circle;
      ic2 = c_edge[ie2].circle;
      ic3 = c_edge[ie3].circle;

      for (ii = 0; ii < 3; ++ii) {
         n_ij[ii] = c_circle[ic1].axis[ii];
         n_jk[ii] = c_circle[ic2].axis[ii];
         n_ik[ii] = c_circle[ic3].axis[ii];
         n_ji[ii] = -n_ij[ii];
         n_ki[ii] = -n_ik[ii];
         n_kj[ii] = -n_jk[ii];
      }
      vertex[iv].beta = acos(DOT(n_ik, n_ji));
      vertex[jv].beta = acos(DOT(n_ij, n_kj));
      vertex[kv].beta = acos(DOT(n_jk, n_ki));

      if (!c_edge[ie1].alive || !c_edge[ie2].alive || !c_edge[ie3].alive)
         continue;              /* only surface intact surfaces */
      if (!c_face[iface].alive)
         continue;              /* only surface intact pieces */

      c_face[iface].area = probe_rad * probe_rad *
          (vertex[iv].beta + vertex[jv].beta + vertex[kv].beta - PI);
      total_area += c_face[iface].area;
   }
   *c_area = total_area;
   return n_surfaced;
}

static int saddle_area(int n_faces, SADDLE_FACE saddle_face[],
                       EDGE convex_edge[], EDGE concave_edge[],
                       CIRCLE convex_circle_list[], TORUS torus[],
                       ATOM atom[], RES res[], VERTEX vertex[],
                       PROBE probe[], REAL_T probe_rad,
                       CIRCLE concave_circle_list[], REAL_T * sad_area)
{
   int iface, itorus, icv1, icv2, icircle, jcircle;
   int ia, ja, ii, iv1, circle1, circle2;
   int icc1, icc2;
   POINT d_c1, d_c2, zaxis, uvec, vvec;
   REAL_T theta1, theta2, wrap_angle;
   REAL_T sin_theta1, sin_theta2;
   REAL_T area1, area2, total_area = 0.0;
   int n_surfaced = 0;

   for (iface = 0; iface < n_faces; ++iface) {
      saddle_face[iface].area = 0.0;

      /* identify edges */
      icv1 = saddle_face[iface].e2_convex;
      icv2 = saddle_face[iface].e4_convex;
      icc1 = saddle_face[iface].e1_concave;
      icc2 = saddle_face[iface].e3_concave;
      itorus = saddle_face[iface].torus;

      if (!saddle_face[iface].alive) {
         continue;              /* only surface intact saddles */
      }
      ++n_surfaced;

      ia = torus[itorus].a1;
      ja = torus[itorus].a2;

      icircle = torus[itorus].circle1;
      jcircle = torus[itorus].circle2;

      if (convex_circle_list[icircle].atom_or_probe_num != ia ||
          convex_circle_list[jcircle].atom_or_probe_num != ja) {
         printf("saddle_area() circle mismatch\n");
         exit(ERROR);
      }
      for (ii = 0; ii < 3; ++ii) {

         d_c1[0] = convex_circle_list[icircle].center[0] - atom[ia].pos[0];
         d_c1[1] = convex_circle_list[icircle].center[1] - atom[ia].pos[1];
         d_c1[2] = convex_circle_list[icircle].center[2] - atom[ia].pos[2];

         d_c2[0] = convex_circle_list[jcircle].center[0] - atom[ja].pos[0];
         d_c2[1] = convex_circle_list[jcircle].center[1] - atom[ja].pos[1];
         d_c2[2] = convex_circle_list[jcircle].center[2] - atom[ja].pos[2];

      }

      sin_theta1 = -(d_c1[0] * convex_circle_list[icircle].axis[0] +
                     d_c1[1] * convex_circle_list[icircle].axis[1] +
                     d_c1[2] * convex_circle_list[icircle].axis[2]) /
          atom[ia].rad;

      sin_theta2 = -(d_c2[0] * convex_circle_list[jcircle].axis[0] +
                     d_c2[1] * convex_circle_list[jcircle].axis[1] +
                     d_c2[2] * convex_circle_list[jcircle].axis[2]) /
          atom[ja].rad;

      theta1 = asin(sin_theta1);
      theta2 = asin(sin_theta2);

      /*
         if (theta1 < 0.0 || theta2 < 0.0) {
         printf("negative torus angle\n");
         printf("theta1 %f    sin(theta1) %f\n", theta1, sin_theta1);
         printf("theta2 %f    sin(theta2) %f\n", theta2, sin_theta2);
         printf("one sided torus = %d\n", one_sided_torus(itorus, torus, atom));
         }
       */

      /* identify a probe position to construct vectors */

      iv1 = convex_edge[icv1].vert1;

      if (iv1 == -1) {          /* free torus */
         wrap_angle = 2.0 * PI;
      } else {

         circle1 = concave_edge[icc1].circle;
         circle2 = concave_edge[icc2].circle;

         for (ii = 0; ii < 3; ++ii) {
            /* sign difference between uvec and vvec is correct */
            uvec[ii] = -concave_circle_list[circle1].axis[ii];
            vvec[ii] = concave_circle_list[circle2].axis[ii];
            zaxis[ii] = -torus[itorus].uv[ii];
         }

         wrap_angle = get_angle(uvec, vvec, zaxis);
         if (wrap_angle < 0.0)
            wrap_angle += 2.0 * PI;
      }


      area1 = wrap_angle * (torus[itorus].rad * probe_rad * theta1 -
                            probe_rad * probe_rad * sin_theta1);

      area2 = wrap_angle * (torus[itorus].rad * probe_rad * theta2 -
                            probe_rad * probe_rad * sin_theta2);

      /*
         printf("torus rad %f\n", torus[itorus].rad);
         printf("wrap angle %f\n", Rad2Deg*wrap_angle);
         printf("theta1 %f\n", Rad2Deg*theta1);
         printf("sin theta1 %f\n", sin_theta1);
         printf("area1 %f\n", area1);
         printf("theta2 %f\n", Rad2Deg*theta1);
         printf("sin theta2 %f\n", sin_theta1);
         printf("area2 %f\n", area2);
       */

      saddle_face[iface].area = area1 + area2;

      if (saddle_face[iface].area < 0.0) {
         printf("negative saddle face area!\n");
         printf("area1 %f area 2 %f\n", area1, area2);
         printf("theta1 %f theta2 %f\n", Rad2Deg * theta1,
                Rad2Deg * theta2);
         printf("sin(theta1) %f sin(theta2) %f\n", sin_theta1, sin_theta2);
         printf("wrap angle %f \n", Rad2Deg * wrap_angle);
         printf("torus rad %f\n", torus[itorus].rad);
         printf("probe rad %f\n", probe_rad);
         exit(ERROR);
      }
      total_area += saddle_face[iface].area;
#ifdef DEBUG
      printf
          ("convex face atom %-4s %-4s %4d : %-4s %-4s %4d area: %10.3f\n",
           atom[ia].anam, res[atom[ia].res].nam, res[atom[ia].res].num,
           atom[ja].anam, res[atom[ja].res].nam, res[atom[ja].res].num,
           saddle_face[iface].area);
#endif


   }
   *sad_area = total_area;
   return n_surfaced;
}

static int id_torus(int ia1, int ia2, ATOM atom[], TORUS torus[])
{
   int it, itmp;

   if (ia1 >= ia2) {
      itmp = ia2;
      ia2 = ia1;
      ia1 = itmp;
      /*
         fprintf(stderr, "id_torus(): atom indices reversed\n");
         exit(ERROR);
       */
   }
   for (it = atom[ia1].torus_start; it < atom[ia1].torus_start +
        atom[ia1].ntorus; ++it)
      if (torus[it].a2 == ia2)
         return it;

   fprintf(stderr,
           "id_torus(): Could not find torus for atoms %d and %d\n", ia1,
           ia2);
   exit(ERROR);
   return -1;
}


static void make_cones(int nat, ATOM atom[], int n_torus, TORUS torus[],
                       int n_probes, PROBE probe[],
                       int n_concave_faces, CONCAVE_FACE concave_face[],
                       int n_saddle_faces, SADDLE_FACE saddle_face[],
                       int n_convex_faces, CONVEX_FACE convex_face[],
                       int *n_vertex, VERTEX vertex[],
                       int *n_concave_edges, EDGE concave_edge[],
                       int n_convex_edges, EDGE convex_edge[],
                       int n_convex_circles, CIRCLE convex_circle[],
                       int n_concave_circles, CIRCLE concave_circle[],
                       int n_cycles, CYCLE cycle[], REAL_T probe_rad,
                       int *n_low_torus, LOW_TORUS low_torus[],
                       CONE_FACE cone_face[], int *n_cone_faces)
{
   int i, it, ilow, iface;
   int nlow = 0, e1, e3;
   REAL_T probe_diam;
   void add_edge();
   int ncones;
   int vertex1, vertex2;

   probe_diam = 2 * probe_rad;

   for (i = 0; i < n_torus; ++i) {
      if (torus[i].low && !one_sided_torus(i, torus, atom)) {
         low_torus[nlow].itorus = i;
         low_torus[nlow].vert1 = *n_vertex;
         low_torus[nlow].vert2 = (*n_vertex) + 1;
         add_2_verts(probe_rad, i, torus, n_vertex, vertex);
         /*
            if (torus[i].n_concave_edges ==0) {
            printf("low torus %d is free!\n", i);
            printf("atoms %d %d\n", torus[i].a1, torus[i].a2);
            }
          */
         low_torus[nlow].ncones = 0;
         ++nlow;
         if (nlow >= NUM_TORUS * natm_sel) {
            printf("MAX_LOW_TORUS exceeded %d\n", nlow);
            exit(ERROR);
         }
      }
   }
   *n_low_torus = nlow;
#ifdef DEBUG
   printf("number of low tori %d\n", *n_low_torus);
#endif
       ncones = 0;
#ifdef DEBUG
   printf("processing %d saddle faces\n", n_saddle_faces);
#endif
   for (iface = 0; iface < n_saddle_faces; ++iface) {

      it = saddle_face[iface].torus;
      saddle_face[iface].alive = 1;

      if (!torus[it].low || one_sided_torus(it, torus, atom))
         continue;

      ilow = get_low_torus_index(nlow, low_torus, it);
      vertex1 = low_torus[ilow].vert1;
      vertex2 = low_torus[ilow].vert2;
      e1 = saddle_face[iface].e1_concave;
      e3 = saddle_face[iface].e3_concave;

      kill_saddle_face(it, torus, iface, saddle_face, concave_edge, e1,
                       e3);

      low_torus[ilow].cone[low_torus[ilow].ncones] = ncones;
      low_torus[ilow].ncones = low_torus[ilow].ncones + 1;
      cone_init(ncones, cone_face, it, saddle_face[iface].e4_convex,
                vertex1);

      low_torus[ilow].cone[low_torus[ilow].ncones] = ncones + 1;
      low_torus[ilow].ncones = low_torus[ilow].ncones + 1;
      cone_init(ncones + 1, cone_face, it, saddle_face[iface].e2_convex,
                vertex2);

      if (low_torus[ilow].ncones >= MAXTOR_PROBE) {
         fprintf(stderr, "make_cones() MAXTOR_PROBE exceeded\n");
         exit(ERROR);
      }
      if (torus[it].n_concave_edges == 0) {
         ncones = ncones + 2;
      } else {

         cone_face[ncones].e2_concave = *n_concave_edges;
         add_edge(n_concave_edges, concave_edge,
                  concave_edge[e1].vert1, vertex1, concave_edge[e1].circle,
                  vertex, concave_circle);

         cone_face[ncones].e3_concave = *n_concave_edges;
         add_edge(n_concave_edges, concave_edge,
                  vertex1, concave_edge[e3].vert2, concave_edge[e3].circle,
                  vertex, concave_circle);
         ++ncones;

         cone_face[ncones].e2_concave = *n_concave_edges;
         add_edge(n_concave_edges, concave_edge,
                  concave_edge[e3].vert1, low_torus[ilow].vert2,
                  concave_edge[e3].circle, vertex, concave_circle);

         cone_face[ncones].e3_concave = *n_concave_edges;
         add_edge(n_concave_edges, concave_edge,
                  low_torus[ilow].vert2, concave_edge[e1].vert2,
                  concave_edge[e1].circle, vertex, concave_circle);

         ++ncones;
      }
   }

   /*
      for (ilow = 0; ilow < *n_low_torus; ++ilow) {
      printf("low torus %d has %d cones: ", ilow, low_torus[ilow].ncones);
      for (i = 0; i < low_torus[ilow].ncones; ++i)
      printf(" %d", low_torus[ilow].cone[i]);
      printf("\n");
      }
    */

#ifdef DEBUG
   printf("%d cone faces\n", ncones);
#endif
   /* write_cone_info(ncones, cone_face, concave_edge, convex_edge); */
   *n_cone_faces = ncones;

   return;
}

/************************************************************************************/
static int one_sided_torus(int i, TORUS torus[], ATOM atom[])
{
   int ia, ja;

   ia = torus[i].a1;
   ja = torus[i].a2;

   if (DIST(atom[ia].pos, torus[i].center) >
       DIST(atom[ia].pos, atom[ja].pos) ||
       DIST(atom[ja].pos, torus[i].center) >
       DIST(atom[ia].pos, atom[ja].pos)
       ) {
      return 1;
   }
   return 0;
}

/************************************************************************************/
static void add_2_verts(REAL_T prad, int i, TORUS torus[], int *n_vertex,
                        VERTEX vertex[])
{
   int ii, nv;
   REAL_T d_tq;

   nv = *n_vertex;

   d_tq = sqrt(prad * prad - torus[i].rad * torus[i].rad);

   for (ii = 0; ii < 3; ++ii) {
      vertex[nv].pos[ii] = torus[i].center[ii] - d_tq * torus[i].uv[ii];
      vertex[nv + 1].pos[ii] =
          torus[i].center[ii] + d_tq * torus[i].uv[ii];
   }

   vertex[nv].iatom = torus[i].a1;
   vertex[nv + 1].iatom = torus[i].a2;
   vertex[nv].iprobe = -1;      /* no probe associated with cusp point */
   vertex[nv + 1].iprobe = -1;

   nv = nv + 2;
   if (nv >= NUM_VERTEX * natm_sel) {
      printf("MAX_VERTS exceeded %d\n", nv);
      exit(ERROR);
   }
   *n_vertex = nv;

   return;
}

/************************************************************************************/
static int get_low_torus_index(int nlow, LOW_TORUS low_torus[], int it)
{
   int i;

   for (i = 0; i < nlow; ++i) {
      if (low_torus[i].itorus == it)
         return i;
   }

   printf("get_low_torus_index() low torus not found!\n");
   exit(ERROR);
   return -1;
}

/************************************************************************************/
static void kill_saddle_face(int it, TORUS torus[], int iface,
                             SADDLE_FACE saddle_face[],
                             EDGE concave_edge[], int e1, int e3)
{
   saddle_face[iface].alive = 0;
   if (torus[it].n_concave_edges > 0) { /* mark concave edges as dead */
      concave_edge[e1].alive = 0;
      concave_edge[e3].alive = 0;
   }
   return;
}

/************************************************************************************/
static void cone_init(int ncones, CONE_FACE cone_face[], int it,
                      int iedge_convex, int ivertex)
{

   if (ncones >= NUM_FACE * natm_sel) {
      printf("MAX_CONE_FACE exceeded\n");
      exit(ERROR);
   }
   cone_face[ncones].itorus = it;
   cone_face[ncones].e1_convex = iedge_convex;
   cone_face[ncones].e2_concave = -1;
   cone_face[ncones].e3_concave = -1;
   cone_face[ncones].cusp_vertex = ivertex;
   return;
}


/*********************************************************************/
/*  make_broken_faces()                                              */
/* for each concave face with a dead edge (because of a low torus)   */
/* initialize a broken_concave_face, and a 3-edged concave_cycle     */
/* associated with it.  No trimming of a broken_concave_face         */
/* occurs in this routine, only initialization of the cycle,         */
/* so this is just recasting the concave edge using the cycle        */
/* data structures.                                                  */
/*                                                                   */
/* Also, when broken_concave_face keep track of the low tori         */
/* and associate that particular probe with the probe list for       */
/* the low torus.  This will facilitate the axial cusp trimming      */
/* in later routines.                                                */
/* as low_torus[i].probe[count] = broken_concave_face[iface].probe;  */
/*********************************************************************/
static void make_broken_faces(int nat, ATOM atom[], int n_torus,
                              TORUS torus[], int n_probes, PROBE probe[],
                              int n_concave_faces,
                              CONCAVE_FACE concave_face[],
                              int n_saddle_faces,
                              SADDLE_FACE saddle_face[],
                              int n_convex_faces,
                              CONVEX_FACE convex_face[], int *n_vertex,
                              VERTEX vertex[], int n_concave_edges,
                              EDGE concave_edge[], int n_convex_edges,
                              EDGE convex_edge[], int n_convex_circles,
                              CIRCLE convex_circle[],
                              int n_concave_circles,
                              CIRCLE concave_circle[], int n_cycles,
                              CYCLE cycle[], REAL_T probe_rad,
                              int *n_low_torus, LOW_TORUS low_torus[],
                              int *n_broken_concave_faces,
                              BROKEN_CONCAVE_FACE broken_concave_face[],
                              CONCAVE_CYCLE concave_cycle[],
                              int *n_concave_cycles)
{
   int n_low_probes = 0, n_dead;
   int i, iface;
   int e1, e2, e3;
   int c1, c2, c3;
   int n, count, icycle;

   for (i = 0; i < n_probes; ++i)
      if (probe[i].low)
         ++n_low_probes;

#ifdef DEBUG
   printf("number of low probes %d\n", n_low_probes);
#endif


   n = 0;
   icycle = 0;
   for (iface = 0; iface < n_concave_faces; ++iface) {
      n_dead = 0;
      e1 = concave_face[iface].e1;
      e2 = concave_face[iface].e2;
      e3 = concave_face[iface].e3;

      if (!concave_edge[e1].alive)
         ++n_dead;
      if (!concave_edge[e2].alive)
         ++n_dead;
      if (!concave_edge[e3].alive)
         ++n_dead;
      if (n_dead) {
         concave_face[iface].alive = 0;
         c1 = concave_edge[e1].circle;
         c2 = concave_edge[e2].circle;
         c3 = concave_edge[e3].circle;
         broken_concave_face[n].alive = 1;
         broken_concave_face[n].itorus[0] = concave_circle[c1].torus;
         broken_concave_face[n].itorus[1] = concave_circle[c2].torus;
         broken_concave_face[n].itorus[2] = concave_circle[c3].torus;
         broken_concave_face[n].probe = concave_face[iface].probe;
         broken_concave_face[n].n_cycles = 1;
         broken_concave_face[n].concave_cycle[0] = icycle;
         add_concave_cycle(icycle, concave_cycle, e1, e2, e3,
                           concave_face[iface].probe, n);
         /* broken_face_info(n, broken_concave_face, concave_cycle); */
         ++n;
         ++icycle;
      }
   }
   *n_broken_concave_faces = n;
   *n_concave_cycles = icycle;

#ifdef DEBUG
   printf("%d concave faces with low tori\n", *n_broken_concave_faces);
#endif
   for (iface = 0; iface < *n_broken_concave_faces; ++iface) {
      count = 0;
      for (i = 0; i < 3; ++i) {
         if (torus[broken_concave_face[iface].itorus[i]].low)
            ++count;
      }
      /*
         if (count > 1) printf("face %d cycle %d has %d low tori\n", iface, 
         broken_concave_face[iface].concave_cycle[0], count);
       */
   }

   for (i = 0; i < *n_low_torus; ++i) {
      count = 0;
      for (iface = 0; iface < *n_broken_concave_faces; ++iface) {
         if (broken_concave_face[iface].itorus[0] == low_torus[i].itorus ||
             broken_concave_face[iface].itorus[1] == low_torus[i].itorus ||
             broken_concave_face[iface].itorus[2] == low_torus[i].itorus) {
            low_torus[i].face[count] = iface;
            ++count;
            if (count >= MAXTOR_PROBE) {
               printf("trim_probes(): MAXTOR_PROBE exceeded\n");
               exit(ERROR);
            }
         }
      }
      low_torus[i].nfaces = count;
      if (low_torus[i].nfaces % 2 != 0) {
         printf("trim_probes(): n broken faces on low torus is odd\n");
         exit(ERROR);
      }
      if (low_torus[i].nfaces == 0 &&
          torus[low_torus[i].itorus].n_concave_edges > 0) {
         printf("trim_probes() low torus has no faces, but is not free\n");
         exit(ERROR);
      }
      /*
         printf("%d probes associated with low torus %d\n", 
         low_torus[i].nfaces, low_torus[i].itorus);
       */
   }
   return;
}

static void add_concave_cycle(int i, CONCAVE_CYCLE concave_cycle[],
                              int e1, int e2, int e3, int iprobe,
                              int iface)
{
   if (i > NUM_CYCLE * natm_sel) {
      fprintf(stderr, "add_concave_cycle() MAX_CYCLES exceeded\n");
      exit(ERROR);
   }
   concave_cycle[i].nedges = 3;
   concave_cycle[i].edge[0] = e1;
   concave_cycle[i].edge_direction[0] = 1;
   concave_cycle[i].cusp_edge[0] = -1;
   concave_cycle[i].edge[1] = e2;
   concave_cycle[i].edge_direction[1] = 1;
   concave_cycle[i].cusp_edge[1] = -1;
   concave_cycle[i].edge[2] = e3;
   concave_cycle[i].edge_direction[2] = 1;
   concave_cycle[i].cusp_edge[2] = -1;
   concave_cycle[i].iprobe = iprobe;
   concave_cycle[i].area = 0.0;
   concave_cycle[i].iface = iface;
   return;
}


static REAL_T get_cone_area(int n_faces, CONE_FACE cone_face[],
                            EDGE convex_edge[], EDGE concave_edge[],
                            CIRCLE convex_circle[],
                            CIRCLE concave_circle[], TORUS torus[],
                            ATOM atom[], VERTEX vertex[], PROBE probe[],
                            REAL_T probe_rad)
{
   int iface, itorus, icv1, icircle;
   int ii, iv1, iv2, ia;
   POINT zaxis, uvec, vvec;
   REAL_T theta1, theta2, wrap_angle;
   REAL_T total_area = 0.0;
   int cusp_vert;


   for (iface = 0; iface < n_faces; ++iface) {

      /* identify edges */
      icv1 = cone_face[iface].e1_convex;
      itorus = cone_face[iface].itorus;
      cusp_vert = cone_face[iface].cusp_vertex;

      icircle = convex_edge[icv1].circle;

      if (convex_edge[icv1].vert1 == -1) {      /* free torus */
         wrap_angle = 2.0 * PI;
      } else {
         iv1 = convex_edge[icv1].vert2;
         iv2 = convex_edge[icv1].vert1;
         for (ii = 0; ii < 3; ++ii) {
            uvec[ii] =
                vertex[iv1].pos[ii] - convex_circle[icircle].center[ii];
            vvec[ii] =
                vertex[iv2].pos[ii] - convex_circle[icircle].center[ii];
            zaxis[ii] = -convex_circle[icircle].axis[ii];
         }
         vnorm(uvec, 3);
         vnorm(vvec, 3);
         wrap_angle = get_angle(uvec, vvec, convex_circle[icircle].axis);
         if (wrap_angle < 0.0)
            wrap_angle += 2.0 * PI;
      }

      ia = convex_circle[icircle].atom_or_probe_num;

      theta1 = acos(torus[itorus].rad / probe_rad);
      theta2 = acos(torus[itorus].rad / (atom[ia].rad + probe_rad));

      if (theta1 < 0.0 || theta2 < 0.0 || theta2 < theta1) {
         printf("theta negative for cone face\n");
         exit(ERROR);
      }
      cone_face[iface].area = wrap_angle * probe_rad * probe_rad *
          (cos(theta1) * (theta2 - theta1) - (sin(theta2) - sin(theta1)));

      /*
         printf("cone face area: face %d theta1 %f theta2 %f wrap_angle %f area %f\n",
         iface, Rad2Deg*theta1, Rad2Deg*theta2, 
         Rad2Deg*wrap_angle, cone_face[iface].area);

         printf("  torus center %8.3f%8.3f%8.3f\n", 
         torus[itorus].center[0],  
         torus[itorus].center[1], 
         torus[itorus].center[2]);

         printf("  vertex center %8.3f%8.3f%8.3f\n", 
         vertex[cusp_vert].pos[0],  
         vertex[cusp_vert].pos[1], 
         vertex[cusp_vert].pos[2]);
       */

      total_area += cone_face[iface].area;

   }
   return total_area;
}

/* axial_trim(): creates cusp edges for intersecting
   concave faces associated with low tori.  These
   cusp edges can intersect on broken concave faces that
   have two or more low tori associated with them.
   If they do intersect, the borken concave face must be
   broken up into more than one cycle (two or more separate faces),
   but this is done in a later routine.


   need to properly pair broken concave faces:
   faces must be sorted by counterclockwise angle
   from the starting probe, ala saddle surfaces.
   The starting probe is determined as the probe
   (associated with with a broken face) that has the
   solvent accessible volume to the right of the
   probe when viewed down the torus axis.  This
   is determined by looking at cross products of
   the vectors that point from the probe to the verts
   associated with the edge that is part of the torus
   in question... simple right?  Here's a picture:
   /~~~~*
   #

   *
   |
   |        (x) torus axis (into the page)
   |
   #

   * and # are probe positions.  The # probes can be starting
   probes, the * probes cannot.  Once probes are sorted this
   way, then probes are pair 1,2  3,4  5,6 etc. for cusp
   trimming between the broken faces.
 */

static void axial_trim(int nat, ATOM atom[], RES res[], int n_torus,
                       TORUS toruslist[], int n_probes, PROBE probelist[],
                       int n_concave_faces, CONCAVE_FACE concave_face[],
                       int n_vertex, VERTEX vertexlist[],
                       int *n_concave_edges, EDGE concave_edge[],
                       int *n_concave_circles, CIRCLE concave_circle[],
                       REAL_T probe_rad, int n_low_torus,
                       LOW_TORUS low_torus[], int n_broken_concave_faces,
                       BROKEN_CONCAVE_FACE broken_concave_face[],
                       CONCAVE_CYCLE concave_cycle[], int n_concave_cycles,
                       CONE_FACE cone_face[], CUSP_EDGE cusp_edge[],
                       int *n_cusps)
{

   int face_edge[MAXTOR_PROBE];
   REAL_T face_angle[MAXTOR_PROBE];     /* angle from starting probe */
   POINT face_vector[MAXTOR_PROBE];     /* vector from torus center to probe pos */
   int it, iface, jface, iprobe, itorus, iedge, jedge;
   int a1, a2;
   int icycle, jcycle, i, ii;
   int j;
   int f1_edge1, f1_edge2, f1_edge3;
   int f2_edge1, f2_edge2, f2_edge3;
   int jprobe;

#ifdef DEBUG
   printf("axial trim\n");
   printf("number of concave_circles %d\n", *n_concave_circles);
#endif


/* first do some some  checks */
   for (it = 0; it < n_low_torus; ++it) {
      itorus = low_torus[it].itorus;
      if (low_torus[it].nfaces >= MAXTOR_PROBE) {
         fprintf(stderr, "axial_trim(): MAXTOR_PROBE exceeded\n");
         exit(ERROR);
      }
      for (i = 0; i < low_torus[it].nfaces; ++i) {
         iface = low_torus[it].face[i];
         if (broken_concave_face[iface].n_cycles != 1) {
            fprintf(stderr, "axial_trim(): n_cycles != 1\n");
            exit(ERROR);
         }
         icycle = broken_concave_face[iface].concave_cycle[0];
         if (concave_cycle[icycle].nedges != 3) {
            fprintf(stderr, "axial_trim(): n_edges != 3\n");
            exit(ERROR);
         }
         for (j = 0; j < concave_cycle[icycle].nedges; ++j) {
            if (concave_cycle[icycle].edge_direction[j] != 1) {
               fprintf(stderr,
                       "axial_trim(): bad edge direction on cycle\n");
               exit(ERROR);
            }
         }
      }
   }

   *n_cusps = 0;
   for (it = 0; it < n_low_torus; ++it) {
      itorus = low_torus[it].itorus;
      if (toruslist[itorus].n_concave_edges == 0)
         continue;              /* free torus */
      a1 = toruslist[itorus].a1;
      a2 = toruslist[itorus].a2;
      for (i = 0; i < low_torus[it].nfaces; ++i) {
         iface = low_torus[it].face[i];
         iprobe = broken_concave_face[iface].probe;
         icycle = broken_concave_face[iface].concave_cycle[0];
         face_edge[i] = -1;

         /* identify the edge on the cycle associated with torus itorus */
         face_edge[i] = get_cycle_edge(itorus, concave_cycle, icycle,
                                       concave_edge, concave_circle);

         for (ii = 0; ii < 3; ++ii)
            face_vector[i][ii] =
                probelist[iprobe].pos[ii] - toruslist[itorus].center[ii];
      }
      face_angle[0] = 0.0;

      for (i = 1; i < low_torus[it].nfaces; ++i) {
         face_angle[i] =
             get_angle(face_vector[i], face_vector[0],
                       toruslist[itorus].uv);
         if (face_angle[i] < 0.0)
            face_angle[i] = face_angle[i] + TWOPI;
      }

      sort_faces(it, low_torus, face_angle, face_edge,
                 n_concave_edges, concave_edge, vertexlist, itorus,
                 toruslist);

      if (low_torus[it].nfaces % 2 != 0) {
         fprintf(stderr, "odd number of faces on torus\n");
         exit(ERROR);
      }
      /* now for each pair of faces, you need to create a new edge on
         the circle of intersection of the probe positions for the pair
         of faces.  You also need to hunt down the cone edges that were
         already created in order to complete the cycle */

      for (i = 0; i < low_torus[it].nfaces; i = i + 2) {
         iface = low_torus[it].face[i];
         iedge = face_edge[i];
         icycle = broken_concave_face[iface].concave_cycle[0];
         jface = low_torus[it].face[i + 1];
         jedge = face_edge[i + 1];
         jcycle = broken_concave_face[jface].concave_cycle[0];

         if (concave_edge[iedge].alive != 0 ||
             concave_edge[jedge].alive != 0) {
            printf("concave edge should already be dead\n");
            exit(ERROR);
         }
         /* 1. using previously found vertices create new circle and edge */
         iprobe = broken_concave_face[iface].probe;
         jprobe = broken_concave_face[jface].probe;

         add_circle(probelist, iprobe, jprobe, concave_circle,
                    n_concave_circles, probe_rad);

         f1_edge2 = *n_concave_edges;
         f2_edge2 = *n_concave_edges;
         cusp_edge[*n_cusps].edge = *n_concave_edges;
         cusp_edge[*n_cusps].probe1 = iprobe;
         cusp_edge[*n_cusps].probe2 = jprobe;
         cusp_edge[*n_cusps].alive = 1;

         /*
            printf("adding edge between vertices %d %d on circle %d \n",
            low_torus[it].vert1, low_torus[it].vert2, *n_concave_circles - 1);
            printf("atom 1 %8.3f%8.3f%8.3f\n", atom[a1].pos[0], atom[a1].pos[1], atom[a1].pos[2]);
            printf("atom 2 %8.3f%8.3f%8.3f\n", atom[a2].pos[0], atom[a2].pos[1], atom[a2].pos[2]);
            printf("circle center %8.3f%8.3f%8.3f\n", 
            concave_circle[*n_concave_circles - 1].center[0],
            concave_circle[*n_concave_circles - 1].center[1],
            concave_circle[*n_concave_circles - 1].center[2]);
            printf("torus center %8.3f%8.3f%8.3f\n",
            toruslist[itorus].center[0],
            toruslist[itorus].center[1],
            toruslist[itorus].center[2]);
            printf("probe 1 %8.3f%8.3f%8.3f\n", 
            probelist[iprobe].pos[0],
            probelist[iprobe].pos[1],
            probelist[iprobe].pos[2]);
            printf("probe 2 %8.3f%8.3f%8.3f\n", 
            probelist[jprobe].pos[0],
            probelist[jprobe].pos[1],
            probelist[jprobe].pos[2]);
          */
         add_edge(n_concave_edges, concave_edge,
                  low_torus[it].vert1, low_torus[it].vert2,
                  *n_concave_circles - 1, vertexlist, concave_circle);

         /* this picture might help here:
            /O\ . /0\
            /   \ /.  \
            /   . ^ .   \                a2
            Face i+1  O    . ^ .    O  Face i       ^
            (aka j)    \    .^ .   /                |torus axis
            \   /.\   /                 a1
            \O/. .\O/
            the "O"'s are the original concave vertices on the soon-to-be-broken
            concave faces that are intersecting. The ^ is the new cusp edge
            where the probes intersect.  The orientation is along "up" or from
            atom 1 in the torus towards atom 2.  An old edge (dotted lines)
            is replaced by 3 edges, two from the cone face, and the new cusp edge.
            The orientation of the middle edges is reversed for the cycle
            associated with face i+1, in order to keep the entire cycle counter
            clockwise.  It may occur that two cusp edges on one face intersect
            each other, but that will be handled later.
          */
         /* 2. add this edge and the 2 cone edges into the 2 broken face cycles */
         f1_edge1 =
             cone_edge(concave_edge[iedge].vert1, low_torus[it].vert2,
                       low_torus, it, concave_edge, cone_face, vertexlist,
                       atom);
         f1_edge3 =
             cone_edge(low_torus[it].vert1, concave_edge[iedge].vert2,
                       low_torus, it, concave_edge, cone_face, vertexlist,
                       atom);
         f2_edge1 =
             cone_edge(concave_edge[jedge].vert1, low_torus[it].vert1,
                       low_torus, it, concave_edge, cone_face, vertexlist,
                       atom);
         f2_edge3 =
             cone_edge(low_torus[it].vert2, concave_edge[jedge].vert2,
                       low_torus, it, concave_edge, cone_face, vertexlist,
                       atom);

         /* insert the 3 edges where there was once one.  */
         add_edges_2_cycle(n_cusps, cusp_edge,
                           concave_cycle, icycle, iedge, f1_edge1,
                           f1_edge2, f1_edge3, -1);
         add_edges_2_cycle(n_cusps, cusp_edge, concave_cycle, jcycle,
                           jedge, f2_edge1, f2_edge2, f2_edge3, 1);
         ++(*n_cusps);
         if (*n_cusps >= NUM_CUSP * natm_sel) {
            fprintf(stderr, "axial_trim(): MAX_CUSPS exceeded\n");
            exit(ERROR);
         }
      }
   }

#ifdef DEBUG
   for (iface = 0; iface < n_broken_concave_faces; ++iface) {
      for (i = 0; i < broken_concave_face[iface].n_cycles; ++i) {
         icycle = broken_concave_face[iface].concave_cycle[i];
         check_cycle(iface, broken_concave_face, icycle, concave_cycle,
                     concave_edge, vertexlist);
      }
   }
#endif

   return;
}

/********************************************************************/
/* identify the edge on the cycle associated with torus itorus */
static int get_cycle_edge(int itorus, CONCAVE_CYCLE concave_cycle[],
                          int icycle, EDGE concave_edge[],
                          CIRCLE concave_circle[])
{
   int i, iedge, icircle;
   for (i = 0; i < concave_cycle[icycle].nedges; ++i) {
      iedge = concave_cycle[icycle].edge[i];
      icircle = concave_edge[iedge].circle;
      if (concave_circle[icircle].torus == itorus)
         return iedge;
   }
   fprintf(stderr, "get_cycle_edge(): could not find edge\n");
   fprintf(stderr, "face edges: ");
   for (i = 0; i < concave_cycle[icycle].nedges; ++i)
      fprintf(stderr, " %d", concave_cycle[icycle].edge[i]);
   fprintf(stderr, "\n");
   exit(ERROR);
   return -1;
}

/********************************************************************/
static void sort_faces(int it, LOW_TORUS low_torus[], REAL_T face_angle[],
                       int face_edge[], int *n_concave_edges,
                       EDGE concave_edge[], VERTEX vertexlist[],
                       int itorus, TORUS toruslist[])
{
   int i, j, n, iv1;
   REAL_T vtmp;
   int iface_tmp, iedge_tmp, itmp, iedge, jtmp, ie;

   n = low_torus[it].nfaces;
   i = n - 1;

   while (i > 0) {
      j = 0;
      while (j < i) {
         if (face_angle[j] > face_angle[j + 1]) {
            vtmp = face_angle[j];
            iface_tmp = low_torus[it].face[j];
            iedge_tmp = face_edge[j];

            face_angle[j] = face_angle[j + 1];
            face_edge[j] = face_edge[j + 1];
            low_torus[it].face[j] = low_torus[it].face[j + 1];

            face_angle[j + 1] = vtmp;
            face_edge[j + 1] = iedge_tmp;
            low_torus[it].face[j + 1] = iface_tmp;
         }
         ++j;
      }
      i = i - 1;
   }

   /* check for off by one error */

   iedge = face_edge[0];
   iv1 = concave_edge[iedge].vert2;
   if (vertexlist[iv1].iatom != toruslist[itorus].a1) { /* we're off by one */
      if (vertexlist[iv1].iatom != toruslist[itorus].a2) {
         fprintf(stderr, "bad vertex\n");
         fprintf(stderr,
                 "concave edge %d ( vert %d atom %d ) (vert %d atom %d ) \n",
                 iedge, concave_edge[iedge].vert1,
                 vertexlist[concave_edge[iedge].vert1].iatom,
                 concave_edge[iedge].vert2,
                 vertexlist[concave_edge[iedge].vert2].iatom);
         fprintf(stderr, "torus atoms %d %d\n", toruslist[itorus].a1,
                 toruslist[itorus].a2);
         fprintf(stderr, "here are all the %d edges vertices and atoms:\n",
                 *n_concave_edges);
         for (ie = 0; ie < *n_concave_edges; ++ie) {
            fprintf(stderr,
                    "edge: %10d  ( v1 %10d a1 %10d )  (v2 %10d a2 %d )\n",
                    ie, concave_edge[ie].vert1,
                    vertexlist[concave_edge[ie].vert1].iatom,
                    concave_edge[ie].vert2,
                    vertexlist[concave_edge[ie].vert2].iatom);
         }
         exit(ERROR);
      }
      itmp = low_torus[it].face[0];
      jtmp = face_edge[0];
      vtmp = face_angle[0];
      n = low_torus[it].nfaces;
      for (i = 0; i < n - 1; ++i) {
         low_torus[it].face[i] = low_torus[it].face[i + 1];
         face_angle[i] = face_angle[i + 1];
         face_edge[i] = face_edge[i + 1];
      }
      low_torus[it].face[n - 1] = itmp;
      face_angle[n - 1] = vtmp;
      face_edge[n - 1] = jtmp;
   }
   return;
}

/************************************************************************/
static void add_circle(PROBE probelist[], int iprobe, int jprobe,
                       CIRCLE concave_circle[], int *n_concave_circles,
                       REAL_T probe_rad)
{
   REAL_T d_pp;
   int ii;
   POINT axis;

   /*  create circle */
   d_pp = 0.0;
   for (ii = 0; ii < 3; ++ii) {
      axis[ii] = probelist[iprobe].pos[ii] - probelist[jprobe].pos[ii];
      d_pp = d_pp + axis[ii] * axis[ii];
   }
   d_pp = sqrt(d_pp);
   vnorm(axis, 3);
   concave_circle[*n_concave_circles].torus = -1;
   concave_circle[*n_concave_circles].atom_or_probe_num = -1;
   concave_circle[*n_concave_circles].rad =
       sqrt(probe_rad * probe_rad - (d_pp * d_pp) / 4.0);
   for (ii = 0; ii < 3; ++ii) {
      concave_circle[*n_concave_circles].center[ii] =
          0.5 * (probelist[iprobe].pos[ii] + probelist[jprobe].pos[ii]);
      concave_circle[*n_concave_circles].axis[ii] = axis[ii];
   }
   ++(*n_concave_circles);
   if (*n_concave_circles >= NUM_CIRCLE * natm_sel) {
      fprintf(stderr, "axial_trim(): MAX_CIRCLE exceeded\n");
      exit(ERROR);
   }
}

/************************************************************************/
/* cone_edge() look over the edges on cone faces associated with a      */
/*              low torus to find which edge has v1 and v2 for vertices */
/*              1 and 2, respectively                                   */
/************************************************************************/
static int cone_edge(int v1, int v2, LOW_TORUS low_torus[], int it,
                     EDGE concave_edge[], CONE_FACE cone_face[],
                     VERTEX vertexlist[], ATOM atom[])
{
   int ic, icone, e2, e3, iv1, iv2;

   for (ic = 0; ic < low_torus[it].ncones; ++ic) {
      icone = low_torus[it].cone[ic];
      e2 = cone_face[icone].e2_concave;
      e3 = cone_face[icone].e3_concave;
      if (concave_edge[e2].vert1 == v1 && concave_edge[e2].vert2 == v2)
         return e2;
      if (concave_edge[e3].vert1 == v1 && concave_edge[e3].vert2 == v2)
         return e3;
   }
   fprintf(stderr, "cone_edge(): could not fine cone edges\n");
   fprintf(stderr, "low torus: %d = torus %d\n", it, low_torus[it].itorus);
   fprintf(stderr, "  looking for edge with verts %d %d\n", v1, v2);
   fprintf(stderr, "  and found:\n");

   for (ic = 0; ic < low_torus[it].ncones; ++ic) {
      icone = low_torus[it].cone[ic];
      e2 = cone_face[icone].e2_concave;
      e3 = cone_face[icone].e3_concave;
      iv1 = concave_edge[e2].vert1;
      iv2 = concave_edge[e2].vert2;
      fprintf(stderr, "ic %d cone %d edge %d verts %d %d\n", ic, icone, e2,
              iv1, iv2);
      fprintf(stderr, "iv1: %8.3f%8.3f%8.3f atom %d\n",
              vertexlist[iv1].pos[0], vertexlist[iv1].pos[1],
              vertexlist[iv1].pos[2], vertexlist[iv1].iatom);
      fprintf(stderr, "iv2: %8.3f%8.3f%8.3f atom %d\n",
              vertexlist[iv2].pos[0], vertexlist[iv2].pos[1],
              vertexlist[iv2].pos[2], vertexlist[iv2].iatom);
      iv1 = concave_edge[e3].vert1;
      iv2 = concave_edge[e3].vert2;
      fprintf(stderr, "ic %d cone %d edge %d verts %d %d\n", ic, icone, e3,
              iv1, iv2);
      fprintf(stderr, "iv1: %8.3f%8.3f%8.3f atom %d\n",
              vertexlist[iv1].pos[0], vertexlist[iv1].pos[1],
              vertexlist[iv1].pos[2], vertexlist[iv1].iatom);
      fprintf(stderr, "iv2: %8.3f%8.3f%8.3f atom %d\n",
              vertexlist[iv2].pos[0], vertexlist[iv2].pos[1],
              vertexlist[iv2].pos[2], vertexlist[iv2].iatom);
   }


   exit(ERROR);
   return -1;
}

/************************************************************************/
static void add_edges_2_cycle(int *n_cusps, CUSP_EDGE cusp_edge[],
                              CONCAVE_CYCLE concave_cycle[], int icycle,
                              int iedge, int new_edge1, int new_edge2,
                              int new_edge3, int direction_of_cusp_edge)
{
   int i, j, itmp, n;

   itmp = -1;
   for (i = 0; i < concave_cycle[icycle].nedges; ++i) {
      if (concave_cycle[icycle].edge[i] == iedge) {
         itmp = i;
         continue;
      }
   }
   if (itmp < 0)
      fprintf(stderr,
              "add_edges_2_cycle(): could not find edge to replace\n");
   n = concave_cycle[icycle].nedges + 2;
   if (n >= NUM_FACE * natm_sel) {
      fprintf(stderr, "add_edges_2_cycle(): MAX_FACE_EDGE exceeded\n");
      exit(ERROR);
   }
   /* shift edges that are to the right of itmp by 2 places */
   for (j = 0; j < 2; ++j) {
      for (i = n - 1; i > itmp + 1; i = i - 1) {
         concave_cycle[icycle].edge[i] = concave_cycle[icycle].edge[i - 1];
         concave_cycle[icycle].edge_direction[i] =
             concave_cycle[icycle].edge_direction[i - 1];
         concave_cycle[icycle].cusp_edge[i] =
             concave_cycle[icycle].cusp_edge[i - 1];
      }
   }
   concave_cycle[icycle].edge[itmp] = new_edge1;
   concave_cycle[icycle].edge_direction[itmp] = 1;
   concave_cycle[icycle].cusp_edge[itmp] = -1;

   concave_cycle[icycle].edge[itmp + 1] = new_edge2;
   concave_cycle[icycle].edge_direction[itmp + 1] = direction_of_cusp_edge;
   if (direction_of_cusp_edge == -1) {
      cusp_edge[*n_cusps].cycle1 = icycle;
   } else {
      cusp_edge[*n_cusps].cycle2 = icycle;
   }
   concave_cycle[icycle].cusp_edge[itmp + 1] = *n_cusps;

   concave_cycle[icycle].edge[itmp + 2] = new_edge3;
   concave_cycle[icycle].edge_direction[itmp + 2] = 1;
   concave_cycle[icycle].cusp_edge[itmp + 2] = -1;
   concave_cycle[icycle].nedges = n;
   /* printf("cycle %d has %d edges:", icycle, concave_cycle[icycle].nedges);
      for (i = 0; i < concave_cycle[icycle].nedges; ++i)
      printf(" %d", concave_cycle[icycle].edge[i]);
      printf("\n"); printf("--------------------------\n"); */

   return;
}

#ifdef DEBUG
/************************************************************************/
static void check_cycle(int iface,
                        BROKEN_CONCAVE_FACE broken_concave_face[],
                        int icycle, CONCAVE_CYCLE concave_cycle[],
                        EDGE concave_edge[], VERTEX vertexlist[])
{
   int i, iprobe, iedge, v1, v2;

   printf("face %d cycle %d\n", iface, icycle);

   iprobe = broken_concave_face[iface].probe;
   if (iprobe != concave_cycle[icycle].iprobe) {
      fprintf(stderr,
              "check_cycle(): face and cycle have different probes\n");
      exit(ERROR);
   }
   printf("vertices: ");
   for (i = 0; i < concave_cycle[icycle].nedges; ++i) {
      iedge = concave_cycle[icycle].edge[i];
      if (concave_cycle[icycle].edge_direction[i] == 1) {
         v1 = concave_edge[iedge].vert1;
         v2 = concave_edge[iedge].vert2;
      } else {
         v2 = concave_edge[iedge].vert1;
         v1 = concave_edge[iedge].vert2;
      }
      printf(" %d %d   ", v1, v2);
   }
   printf("\n");
}
#endif

/* concentric_axial_cusps(): at this stage if a face has 2 or 3 
   axial cusps from the same probe, then it must be split up, and
   so must the face adjoining those cusps. */

static void concentric_axial_cusps(int *n_concave_edges,
                                   EDGE concave_edge[],
                                   int *n_broken_concave_faces,
                                   BROKEN_CONCAVE_FACE
                                   broken_concave_face[],
                                   CONCAVE_CYCLE concave_cycle[],
                                   int *n_concave_cycles,
                                   CUSP_EDGE cusp_edge[], int *n_cusps)
{
   int icycle, jcycle, iprobe, jprobe, icusp, jcusp;

#ifdef DEBUG
   printf("trimming concentric axial cusps\n");
#endif
   for (icusp = 0; icusp < *n_cusps; ++icusp) {
      iprobe = cusp_edge[icusp].probe1;
      jprobe = cusp_edge[icusp].probe2;
      icycle = cusp_edge[icusp].cycle1;
      jcycle = cusp_edge[icusp].cycle2;
      for (jcusp = icusp + 1; jcusp < *n_cusps; ++jcusp) {
         if ((cusp_edge[jcusp].probe1 == iprobe ||
              cusp_edge[jcusp].probe1 == jprobe) &&
             (cusp_edge[jcusp].probe2 == iprobe ||
              cusp_edge[jcusp].probe2 == jprobe)) {

            if ((cusp_edge[jcusp].cycle1 != icycle &&
                 cusp_edge[jcusp].cycle1 != jcycle) ||
                (cusp_edge[jcusp].cycle2 != icycle &&
                 cusp_edge[jcusp].cycle2 != jcycle)) {

               printf("concentric_axial_cusps(): cycle mismatch\n");
               exit(ERROR);
            }
#ifdef DEBUG
            printf("rerouting cycles %d and %d\n", icycle, jcycle);
#endif
            reroute(icycle, jcycle, icusp, jcusp, n_concave_cycles,
                    concave_cycle, cusp_edge, concave_edge,
                    n_broken_concave_faces, broken_concave_face);
            cusp_edge[icusp].concentric_pair = 1;
            cusp_edge[jcusp].concentric_pair = 1;
         }
      }
   }
}

/***************************************************************/

static void reroute(int icycle, int jcycle, int cusp1, int cusp2,
                    int *n_concave_cycles, CONCAVE_CYCLE concave_cycle[],
                    CUSP_EDGE cusp_edge[], EDGE concave_edge[],
                    int *n_broken_concave_faces,
                    BROKEN_CONCAVE_FACE broken_concave_face[])
{
   int ii, jj, iedge1, iedge2, ie1, ie2, je;
   int direction1, endpoint1, direction2, endpoint2;

   for (ii = 0; ii < concave_cycle[icycle].nedges; ++ii) {
      if (concave_cycle[icycle].cusp_edge[ii] == cusp1)
         iedge1 = ii;
      if (concave_cycle[icycle].cusp_edge[ii] == cusp2)
         iedge2 = ii;
   }
   ie1 = cusp_edge[concave_cycle[icycle].cusp_edge[iedge1]].edge;
   ie2 = cusp_edge[concave_cycle[icycle].cusp_edge[iedge2]].edge;
   direction1 = concave_cycle[icycle].edge_direction[iedge1];
   direction2 = concave_cycle[icycle].edge_direction[iedge2];

   if (direction1 == 1)
      endpoint1 = concave_edge[ie1].vert2;
   else
      endpoint1 = concave_edge[ie1].vert1;

   if (direction2 == 1)
      endpoint2 = concave_edge[ie2].vert2;
   else
      endpoint2 = concave_edge[ie2].vert1;

   if (direction1 == 1)
      concave_edge[ie1].vert2 = endpoint2;
   else
      concave_edge[ie1].vert1 = endpoint2;

   if (direction2 == 1)
      concave_edge[ie2].vert2 = endpoint1;
   else
      concave_edge[ie2].vert1 = endpoint1;

   /* need to correct jcycle now */
   for (jj = 0; jj < concave_cycle[jcycle].nedges; ++jj) {
      if (concave_cycle[jcycle].cusp_edge[jj] != -1) {
         je = cusp_edge[concave_cycle[jcycle].cusp_edge[jj]].edge;
         if (je == ie1) {
            concave_cycle[jcycle].edge_direction[jj] =
                -concave_cycle[icycle].edge_direction[iedge1];
            concave_cycle[jcycle].edge[jj] =
                concave_cycle[icycle].edge[iedge1];
         } else if (je == ie2) {
            concave_cycle[jcycle].edge_direction[jj] =
                -concave_cycle[icycle].edge_direction[iedge2];
            concave_cycle[jcycle].edge[jj] =
                concave_cycle[icycle].edge[iedge2];
         } else {
            printf("reroute(): cusp edge mismatch\n");
            exit(ERROR);
         }
      }
   }
   /* dump_cycle(&concave_cycle[icycle], concave_edge);
      dump_cycle(&concave_cycle[jcycle], concave_edge); */
   split_cycle(n_broken_concave_faces, broken_concave_face,
               n_concave_cycles, concave_cycle, icycle, concave_edge,
               cusp_edge);
   split_cycle(n_broken_concave_faces, broken_concave_face,
               n_concave_cycles, concave_cycle, jcycle, concave_edge,
               cusp_edge);
   return;
}

/*****************************************************************/

static void split_cycle(int *n_broken_concave_faces,
                        BROKEN_CONCAVE_FACE broken_concave_face[],
                        int *n_concave_cycles,
                        CONCAVE_CYCLE concave_cycle[], int icycle,
                        EDGE concave_edge[], CUSP_EDGE cusp_edge[])
{
   CONCAVE_CYCLE tmp_cycle;
   int *edge_used;
   int vertex2, iface, ntot, i, first_vert, next_vert, ne, iedge;
   int ncycle, istart, no_start, icusp;

   if ((edge_used =
        (int *) malloc(NUM_EDGE * natm_sel * sizeof(int))) == NULL) {
      fprintf(stderr, "Unable to allocate space for edge_used\n");
      exit(1);
   }
   vertex2 = -1;
   iface = concave_cycle[icycle].iface;
   ntot = concave_cycle[icycle].nedges;

/********** shorten icycle **************************************/
   tmp_cycle.nedges = concave_cycle[icycle].nedges;
   tmp_cycle.iprobe = concave_cycle[icycle].iprobe;
   tmp_cycle.iface = concave_cycle[icycle].iface;
   for (i = 0; i < concave_cycle[icycle].nedges; ++i) {
      tmp_cycle.edge[i] = concave_cycle[icycle].edge[i];
      tmp_cycle.cusp_edge[i] = concave_cycle[icycle].cusp_edge[i];
      tmp_cycle.edge_direction[i] =
          concave_cycle[icycle].edge_direction[i];
      edge_used[i] = 0;
   }

   concave_cycle[icycle].edge[0] = tmp_cycle.edge[0];
   concave_cycle[icycle].cusp_edge[0] = tmp_cycle.cusp_edge[0];
   concave_cycle[icycle].edge_direction[0] = tmp_cycle.edge_direction[0];
   concave_cycle[icycle].iprobe = tmp_cycle.iprobe;
   concave_cycle[icycle].iface = concave_cycle[icycle].iface;

   if (concave_cycle[icycle].edge_direction[0] == 1) {
      first_vert = concave_edge[concave_cycle[icycle].edge[0]].vert1;
      next_vert = concave_edge[concave_cycle[icycle].edge[0]].vert2;
   } else {
      next_vert = concave_edge[concave_cycle[icycle].edge[0]].vert1;
      first_vert = concave_edge[concave_cycle[icycle].edge[0]].vert2;
   }
   edge_used[0] = 1;
   ne = 1;
   while (next_vert != first_vert) {
      iedge =
          next_cycle_edge(tmp_cycle, concave_edge, next_vert, edge_used);
      concave_cycle[icycle].edge[ne] = tmp_cycle.edge[iedge];
      concave_cycle[icycle].cusp_edge[ne] = tmp_cycle.cusp_edge[iedge];
      concave_cycle[icycle].edge_direction[ne] =
          tmp_cycle.edge_direction[iedge];
      if (concave_cycle[icycle].edge_direction[ne] == 1)
         next_vert = concave_edge[concave_cycle[icycle].edge[ne]].vert2;
      else
         next_vert = concave_edge[concave_cycle[icycle].edge[ne]].vert1;
      ++ne;
   }
   concave_cycle[icycle].nedges = ne;

/********** create new cycle ************************************/
   ncycle = *n_concave_cycles;
   ne = 0;
   no_start = 1;

   while (no_start) {
      if (!edge_used[ne]) {
         istart = ne;
         no_start = 0;
      }
      ++ne;
      if (ne > ntot) {
         printf
             ("split_cycle(): could not find a starting edge for 2nd cycle\n");
         exit(ERROR);
      }
   }

   concave_cycle[ncycle].edge[0] = tmp_cycle.edge[istart];
   icusp = tmp_cycle.cusp_edge[istart];
   concave_cycle[ncycle].cusp_edge[0] = icusp;
   if (icusp != -1) {
      if (cusp_edge[icusp].cycle1 == icycle)
         cusp_edge[icusp].cycle1 = ncycle;
      if (cusp_edge[icusp].cycle2 == icycle)
         cusp_edge[icusp].cycle2 = ncycle;
   }
   concave_cycle[ncycle].edge_direction[0] =
       tmp_cycle.edge_direction[istart];
   concave_cycle[ncycle].iprobe = tmp_cycle.iprobe;
   concave_cycle[ncycle].iface = tmp_cycle.iface;
   edge_used[istart] = 1;

   if (tmp_cycle.edge_direction[istart] == 1) {
      first_vert = concave_edge[tmp_cycle.edge[istart]].vert1;
      next_vert = concave_edge[tmp_cycle.edge[istart]].vert2;
   } else {
      next_vert = concave_edge[tmp_cycle.edge[istart]].vert1;
      first_vert = concave_edge[tmp_cycle.edge[istart]].vert2;
   }

   ne = 1;
   while (next_vert != first_vert) {
      iedge =
          next_cycle_edge(tmp_cycle, concave_edge, next_vert, edge_used);
      concave_cycle[ncycle].edge[ne] = tmp_cycle.edge[iedge];
      icusp = tmp_cycle.cusp_edge[iedge];
      if (icusp != -1) {
         if (cusp_edge[icusp].cycle1 == icycle)
            cusp_edge[icusp].cycle1 = ncycle;
         if (cusp_edge[icusp].cycle2 == icycle)
            cusp_edge[icusp].cycle2 = ncycle;
      }
      concave_cycle[ncycle].cusp_edge[ne] = icusp;
      concave_cycle[ncycle].edge_direction[ne] =
          tmp_cycle.edge_direction[iedge];
      if (concave_cycle[ncycle].edge_direction[ne] == 1)
         next_vert = concave_edge[concave_cycle[ncycle].edge[ne]].vert2;
      else
         next_vert = concave_edge[concave_cycle[ncycle].edge[ne]].vert1;
      ++ne;
   }
   concave_cycle[ncycle].iprobe = concave_cycle[icycle].iprobe;
   concave_cycle[ncycle].iface = concave_cycle[icycle].iface;
   concave_cycle[ncycle].nedges = ne;
   *n_concave_cycles = ncycle + 1;

   /* check to be sure all edges were used */
   for (i = 0; i < ntot; ++i) {
      if (!edge_used[i]) {
         printf("edge %d not used\n", i);
         exit(ERROR);
      }
   }
   if (broken_concave_face[iface].n_cycles != 1) {
      printf("concentric_axial_cusps(): n_cycles != 1\n");
      exit(ERROR);
   }
   i = broken_concave_face[iface].n_cycles;
   i = *n_broken_concave_faces;
   /* create a new face */
   broken_concave_face[i].itorus[0] = broken_concave_face[iface].itorus[0];
   broken_concave_face[i].itorus[1] = broken_concave_face[iface].itorus[1];
   broken_concave_face[i].itorus[2] = broken_concave_face[iface].itorus[2];
   broken_concave_face[i].probe = broken_concave_face[iface].probe;
   broken_concave_face[i].n_cycles = 1;
   broken_concave_face[i].alive = 1;
   broken_concave_face[i].area = 0.0;
   broken_concave_face[i].concave_cycle[0] = ncycle;
   /*
      dump_cycle(&concave_cycle[icycle], concave_edge);
      dump_cycle(&concave_cycle[ncycle], concave_edge);
    */

   free(edge_used);
   return;
}

/******************************************************************/

static int next_cycle_edge(CONCAVE_CYCLE cycle, EDGE concave_edge[],
                           int next_vert, int edge_used[])
{
   int i;
   for (i = 0; i < cycle.nedges; ++i) {
      if (!edge_used[i] &&
          (concave_edge[cycle.edge[i]].vert1 == next_vert ||
           concave_edge[cycle.edge[i]].vert2 == next_vert)) {
         edge_used[i] = 1;
         return i;
      }
   }
   printf("next_cycle_edge(): could not find next edge with vertex %d\n",
          next_vert);
   exit(ERROR);
   return -1;
}

/*****************************************************************/

static void non_axial_trim(int nat, ATOM atom[], RES res[], int n_torus,
                           TORUS toruslist[], int n_probes,
                           PROBE probelist[], int n_concave_faces,
                           CONCAVE_FACE concave_face[], int n_vertex,
                           VERTEX vertexlist[], int *n_concave_edges,
                           EDGE concave_edge[], int *n_concave_circles,
                           CIRCLE concave_circle_list[], REAL_T probe_rad,
                           int n_low_torus, LOW_TORUS low_torus[],
                           int n_broken_concave_faces,
                           BROKEN_CONCAVE_FACE broken_concave_face[],
                           CONCAVE_CYCLE concave_cycle[],
                           int n_concave_cycles, CONE_FACE cone_face[],
                           CUSP_EDGE cusp_edge[], int *n_cusps,
                           CUSP_PAIR cusp_pair[], int *n_cusp_pairs)
{
   CUSP_GROUP *group;
   int n_groups;
   int icusp, jcusp, ig, i, j, ng;
   int iedge, jedge, icircle, jcircle;
   int cycle_intersect();
   void trim_1_pair(), trim_2_pair(), trim_3_pair();
   int nc;
   int cusp_added;
   int n;

#ifdef DEBUG
   printf("non axial trim\n");
   printf("There are %d cusp edges\n", *n_cusps);
#endif
   if ((group =
        (CUSP_GROUP *) malloc(NUM_CUSP * natm_sel * sizeof(CUSP_GROUP))) ==
       NULL) {
      fprintf(stderr, "Unable to allocate space for group\n");
      exit(1);
   }
   ng = 0;
   *n_cusp_pairs = 0;
   for (icusp = 0; icusp < *n_cusps; ++icusp) {
      iedge = cusp_edge[icusp].edge;
      icircle = concave_edge[iedge].circle;
      if (cusp_edge[icusp].concentric_pair == 1)
         continue;

      /*
         for (jcusp = icusp+1; jcusp < *n_cusps; ++jcusp) {
         if (cusp_edge[jcusp].concentric_pair == 1) continue;
         cusp_intersect(cusp_edge, icusp, jcusp, probelist, concave_cycle,
         concave_edge, vertexlist, concave_circle_list, 
         probe_rad, cusp_pair, n_cusp_pairs);

         }
       */

      for (jcusp = 0; jcusp < *n_cusps; ++jcusp) {
         if (icusp == jcusp)
            continue;
         if (cusp_edge[jcusp].concentric_pair == 1)
            continue;
         jedge = cusp_edge[jcusp].edge;
         jcircle = concave_edge[jedge].circle;
         if (concave_circle_list[icircle].rad <
             concave_circle_list[jcircle].rad) {
            cusp_intersect(cusp_edge, icusp, jcusp, probelist,
                           concave_cycle, concave_edge, vertexlist,
                           concave_circle_list, probe_rad, cusp_pair,
                           n_cusp_pairs);
         }
      }

   }
   for (i = 0; i < *n_cusp_pairs; ++i) {
      cusp_pair[i].group = -1;
   }

#ifdef DEBUG
   printf("number of cusp pairs %d\n", *n_cusp_pairs);
#endif

   ig = 0;
   for (i = 0; i < *n_cusp_pairs; ++i) {
      if (cusp_pair[i].group >= 0)
         continue;
      cusp_pair[i].group = ig;
      nc = 0;
      group[ig].cusp_pair[nc] = i;
      ++nc;
      group[ig].n_pairs = nc;

      cusp_added = 1;

      while (cusp_added) {
         cusp_added = 0;
         for (j = i + 1; j < *n_cusp_pairs; ++j) {
            if (cusp_pair[j].group >= 0)
               continue;
            if (new_cusp_in_group(j, ig, cusp_pair, group)) {
               cusp_pair[j].group = ig;
               group[ig].cusp_pair[nc] = j;
               ++nc;
               group[ig].n_pairs = nc;
               cusp_added = 1;
            }
         }
      }
      ++ig;
   }

   n_groups = ig;

#ifdef DEBUG
   printf("cusp groups\n");
   i = 0;
   for (ig = 0; ig < n_groups; ++ig) {
      printf("group %d has %d cusp_pairs:", ig, group[ig].n_pairs);
      for (nc = 0; nc < group[ig].n_pairs; ++nc) {
         int inew = group[ig].cusp_pair[nc];
         printf("  %d (%d %d)", inew, cusp_pair[inew].cusp1,
                cusp_pair[inew].cusp2);
      }
      printf("\n");
   }
#endif

   for (ig = 0; ig < n_groups; ++ig) {
/******* TEMP ************/
      /*
         if (ig == 0) continue;
         itmp = group[ig].cusp_pair[0];
         group[ig].cusp_pair[0] = group[ig].cusp_pair[1]; 
         group[ig].cusp_pair[1] = itmp;
         group[ig].n_pairs = 1;
         n_groups = 0;
         printf("ig = %d\n", ig);
         cusp1 = cusp_pair[group[ig].cusp_pair[0]].cusp1;
         cusp2 = cusp_pair[group[ig].cusp_pair[0]].cusp2;
         cusp3 = cusp_pair[group[ig].cusp_pair[1]].cusp1;
         cusp4 = cusp_pair[group[ig].cusp_pair[1]].cusp2;
         printf("cusp: %d probes: %d %d\n",  cusp1,
         cusp_edge[cusp1].probe1,
         cusp_edge[cusp1].probe2);
         printf("cusp: %d probes: %d %d\n",  cusp2,
         cusp_edge[cusp2].probe1,
         cusp_edge[cusp2].probe2);
         printf("cusp: %d probes: %d %d\n",  cusp3,
         cusp_edge[cusp3].probe1,
         cusp_edge[cusp3].probe2);
         printf("cusp: %d probes: %d %d\n",  cusp4,
         cusp_edge[cusp4].probe1,
         cusp_edge[cusp4].probe2);

         make_circle(probelist,probe_rad,
         n_concave_circles, concave_circle_list);
         exit(ERROR);
       */

/******* END TEMP ************/
#ifdef DEBUG
      printf("ig: %d\n", ig);
      printf("number of cusps %d\n", *n_cusps);
#endif
      n = number_of_cusps(group, ig, cusp_pair);
      if (n == 1) {
         continue;
      } else if (n == 2) {
         trim_2_cusps(probelist, &n_vertex, vertexlist, n_concave_edges,
                      concave_edge, n_concave_circles, concave_circle_list,
                      probe_rad, &n_broken_concave_faces,
                      broken_concave_face, concave_cycle,
                      &n_concave_cycles, cusp_edge, n_cusps, cusp_pair,
                      n_cusp_pairs, group, ig);
      } else if (n == 3) {
         printf("WARNING: 3 cusps intersection: not trimmed yet\n");
         /*
            trim_3_cusps( probelist, &n_vertex, vertexlist, n_concave_edges,
            concave_edge, n_concave_circles, concave_circle_list,
            probe_rad, &n_broken_concave_faces, broken_concave_face,
            concave_cycle, &n_concave_cycles, cusp_edge, n_cusps,
            cusp_pair, n_cusp_pairs, group, ig);
          */
      } else if (n == 4) {
         printf("WARNING: 4 cusps intersection: not trimmed yet\n");
         /*
            printf("trim 4 cusps\n");
            trim_4_cusps( probelist, &n_vertex, vertexlist, n_concave_edges,
            concave_edge, n_concave_circles, concave_circle_list,
            probe_rad, &n_broken_concave_faces, broken_concave_face,
            concave_cycle, &n_concave_cycles, cusp_edge, n_cusps,
            cusp_pair, n_cusp_pairs, group, ig);
          */
      } else {
         printf("two many interesecting cusps %d\n", n);
         exit(ERROR);
      }
   }

   /*
      for (ig = 0; ig < n_groups; ++ig) {
      for (nc = 0; nc < group[ig].n_pairs; ++ nc) {
      inew = group[ig].cusp_pair[nc]; 
      printf("ATOM      1 O    VAN  %4d    %8.3f%8.3f%8.3f\n", ig,
      cusp_pair[inew].circle_center[0], 
      cusp_pair[inew].circle_center[1], 
      cusp_pair[inew].circle_center[2]);
      printf("ATOM      1 C    VAN  %4d    %8.3f%8.3f%8.3f\n", ig,
      cusp_pair[inew].vert1[0], 
      cusp_pair[inew].vert1[1], 
      cusp_pair[inew].vert1[2]);
      printf("ATOM      1 N    VAN  %4d    %8.3f%8.3f%8.3f\n", ig,
      cusp_pair[inew].vert2[0], 
      cusp_pair[inew].vert2[1], 
      cusp_pair[inew].vert2[2]);
      ++i;
      }
      }
    */
   /*
      i = 0;
      for (ig = 0; ig < n_groups; ++ig) {
      for (nc = 0; nc < group[ig].n_pairs; ++ nc) {
      inew = group[ig].cusp_pair[nc]; 
      printf("%d %d\n", i+1, i+2);
      printf("%d %d\n", i+1, i+3);
      i = i + 3;
      }
      }
    */

   free(group);
   return;
}

/*****************************************************************/
static int new_cusp_in_group(int icusp, int igroup, CUSP_PAIR cusp_pair[],
                             CUSP_GROUP group[])
{
   int nc;
   int jcusp;

   for (nc = 0; nc < group[igroup].n_pairs; ++nc) {
      jcusp = group[igroup].cusp_pair[nc];
      if (cusp_pair[icusp].cusp1 == cusp_pair[jcusp].cusp1 ||
          cusp_pair[icusp].cusp2 == cusp_pair[jcusp].cusp1 ||
          cusp_pair[icusp].cusp1 == cusp_pair[jcusp].cusp2 ||
          cusp_pair[icusp].cusp2 == cusp_pair[jcusp].cusp2)
         return 1;
   }
   return 0;
}

static void cusp_intersect(CUSP_EDGE cusp_edge[], int icusp, int jcusp,
                           PROBE probe[], CONCAVE_CYCLE concave_cycle[],
                           EDGE concave_edge[], VERTEX vertex[],
                           CIRCLE concave_circle[], REAL_T probe_rad,
                           CUSP_PAIR cusp_pair[], int *n_cusp_pairs)
{
   int icycle, iprobe, ie1, ie2, c1, c2, ii;
   POINT uv1, uv2;
   REAL_T h1, h2;
   POINT xaxis, yaxis, zaxis, midpoint;
   REAL_T alpha, offset;
   REAL_T v1_angle, v2_angle, angle_range;
   int ivert1, ivert2;
   POINT v1_vector, v2_vector, midpoint_axis;
   int flag;
   int left_flag;
   FILE *atomfile, *bondfile, *propfile;
   int npt = 0, narc = 0;

   left_flag = 0;

   flag = 0;

   /*
      if ( (icusp ==  4 || icusp == 7) && (jcusp == 4 || jcusp == 7) ) flag = 1;
    */

   if (flag)
      printf("flag icusp %d jcusp %d\n", icusp, jcusp);

   icycle = get_cycle_id(cusp_edge, icusp, jcusp, flag);
   if (icycle == -1)
      return;
   iprobe = concave_cycle[icycle].iprobe;
   ie1 = cusp_edge[icusp].edge;
   ie2 = cusp_edge[jcusp].edge;
   c1 = concave_edge[ie1].circle;
   c2 = concave_edge[ie2].circle;
   h1 = 0.0;
   h2 = 0.0;
   for (ii = 0; ii < 3; ++ii) {
      uv1[ii] = concave_circle[c1].center[ii] - probe[iprobe].pos[ii];
      uv2[ii] = concave_circle[c2].center[ii] - probe[iprobe].pos[ii];
      h1 = h1 + uv1[ii] * uv1[ii];
      h2 = h2 + uv2[ii] * uv2[ii];
   }
   h1 = sqrt(h1);
   h2 = sqrt(h2);
   for (ii = 0; ii < 3; ++ii) {
      uv1[ii] = uv1[ii] / h1;
      uv2[ii] = uv2[ii] / h2;
      yaxis[ii] = uv1[ii];
   }
   vnorm(yaxis, 3);
   cross(uv2, uv1, zaxis);
   vnorm(zaxis, 3);
   cross(yaxis, zaxis, xaxis);
   vnorm(xaxis, 3);

   alpha = get_angle(uv1, uv2, zaxis);

   if (flag)
      printf("alpha %8.1f\n", Rad2Deg * alpha);

   if (alpha < 0.0) {
      printf("alpha < 0 should not happen\n");
      exit(ERROR);
   } else if (alpha > asin(concave_circle[c1].rad / probe_rad) +
              asin(concave_circle[c2].rad / probe_rad)) {
      if (flag)
         printf("circles don't intersect\n");
      return;                   /* circles don't intersect */
   } else if (alpha > 0.5 * PI) {
      if (flag)
         printf("alpha %8.1f > %8.1f\n", Rad2Deg * alpha,
                Rad2Deg * 0.5 * PI);
      alpha = PI - alpha;
      offset = h2 * sin(alpha) + (h1 + h2 * cos(alpha)) / tan(alpha);
      if (flag) {
         printf("alpha > pi/2\n");
         printf("alpha %8.1f\n", Rad2Deg * alpha);
         printf("h1:%8.3f h2:%8.3f\n", h1, h2);
         printf("sin(alpha) %8.3f\n", sin(alpha));
         printf("cos(alpha) %8.3f\n", cos(alpha));
         printf("tan(alpha) %8.3f\n", tan(alpha));
      }
   } else if (h1 * cos(alpha) < h2) {
      offset = h2 * sin(alpha) - (h1 - h2 * cos(alpha)) / tan(alpha);
      if (flag) {
         printf("h1*cos(alpha) < h2\n");
         printf("h1:%8.3f h2:%8.3f\n", h1, h2);
         printf("sin(alpha) %8.3f\n", sin(alpha));
         printf("cos(alpha) %8.3f\n", cos(alpha));
         printf("tan(alpha) %8.3f\n", tan(alpha));
      }
   } else {
      left_flag = 1;
      offset = (h1 - h2 / cos(alpha)) / tan(alpha);
   }

   if (flag) {

      printf("\n------------------\n");
      printf("cusps %d %d\n", icusp, jcusp);
      printf("icycle %d\n", icycle);
      printf("circles %d %d\n", c1, c2);
      printf("probe %d\n", iprobe);
      printf("probe center %8.3f%8.3f%8.3f\n", probe[iprobe].pos[0],
             probe[iprobe].pos[1], probe[iprobe].pos[2]);
      printf("offset %8.3f\n", offset);
      printf("icusp %d jcusp %d\n", icusp, jcusp);
      printf("c1 radius %8.3f c2 radius %8.3f\n", concave_circle[c1].rad,
             concave_circle[c2].rad);
      printf("c1 height %8.3f c2 height %8.3f\n", h1, h2);
      printf("uv1: %8.3f%8.3f%8.3f\n", uv1[0], uv1[1], uv1[2]);
      printf("uv2: %8.3f%8.3f%8.3f\n", uv2[0], uv2[1], uv2[2]);
      printf("alpha %8.1f\n", Rad2Deg * alpha);
      printf("cos(alpha) %8.3f\n", cos(alpha));
      printf("tan(alpha) %8.3f\n", tan(alpha));
      printf("\n");
      printf("drawing circles\n");
      atomfile = fopen("circle.pdb", "w");
      bondfile = fopen("circle.bnd", "w");
      propfile = fopen("circle.prop", "w");

      printf("axis: %8.3f%8.3f%8.3f\n",
             concave_circle[c1].center[0],
             concave_circle[c1].center[1], concave_circle[c1].center[2]);
      printf("axis: %8.3f%8.3f%8.3f\n",
             concave_circle[c1].center[0] + concave_circle[c1].axis[0],
             concave_circle[c1].center[1] + concave_circle[c1].axis[1],
             concave_circle[c1].center[2] + concave_circle[c1].axis[2]);
      draw_circle(&npt, &narc, 100, atomfile, bondfile, propfile,
                  concave_circle[c1].center, concave_circle[c1].rad,
                  concave_circle[c1].axis, 1);

      printf("axis: %8.3f%8.3f%8.3f\n",
             concave_circle[c2].center[0],
             concave_circle[c2].center[1], concave_circle[c2].center[2]);
      printf("axis: %8.3f%8.3f%8.3f\n",
             concave_circle[c2].center[0] + concave_circle[c2].axis[0],
             concave_circle[c2].center[1] + concave_circle[c2].axis[1],
             concave_circle[c2].center[2] + concave_circle[c2].axis[2]);
      draw_circle(&npt, &narc, 100, atomfile, bondfile, propfile,
                  concave_circle[c2].center, concave_circle[c2].rad,
                  concave_circle[c2].axis, 1);



   }
   if (left_flag) {
      for (ii = 0; ii < 3; ++ii)
         midpoint[ii] = concave_circle[c1].center[ii] - offset * xaxis[ii];
   } else {
      for (ii = 0; ii < 3; ++ii)
         midpoint[ii] = concave_circle[c1].center[ii] + offset * xaxis[ii];
   }
   angle_range = acos(offset / concave_circle[c1].rad);
   if (flag) {
      printf("angle range %8.1f\n", Rad2Deg * angle_range);
      printf("offset %8.3f\n", offset);
      printf("circle rad %8.3f\n", concave_circle[c1].rad);
      printf("offset/rad: %8.3f\n", offset / concave_circle[c1].rad);
      printf("xaxis %8.3f%8.3f%8.3f\n", xaxis[0], xaxis[1], xaxis[2]);
   }
   if (flag) {
      printf("midpoint %8.3f%8.3f%8.3f\n", midpoint[0], midpoint[1],
             midpoint[2]);
      printf("angle range %8.1f\n", Rad2Deg * angle_range);
   }
   if (DOT(concave_circle[c1].axis, yaxis) > 0.0) {
      ivert1 = concave_edge[ie1].vert1;
      ivert2 = concave_edge[ie1].vert2;
   } else {
      ivert2 = concave_edge[ie1].vert1;
      ivert1 = concave_edge[ie1].vert2;
      if (flag)
         printf("reversing vertices\n");
   }
   for (ii = 0; ii < 3; ++ii) {
      v1_vector[ii] =
          vertex[ivert1].pos[ii] - concave_circle[c1].center[ii];
      v2_vector[ii] =
          vertex[ivert2].pos[ii] - concave_circle[c1].center[ii];
   }
   if (flag) {
      printf("v1: %8.3f%8.3f%8.3f\n", vertex[ivert1].pos[0],
             vertex[ivert1].pos[1], vertex[ivert1].pos[2]);
      printf("v2: %8.3f%8.3f%8.3f\n", vertex[ivert2].pos[0],
             vertex[ivert2].pos[1], vertex[ivert2].pos[2]);
   }
   for (ii = 0; ii < 3; ++ii)
      midpoint_axis[ii] = midpoint[ii] - concave_circle[c1].center[ii];
   vnorm(v1_vector, 3);
   vnorm(v2_vector, 3);
   vnorm(midpoint_axis, 3);

   v1_angle = get_angle(midpoint_axis, v1_vector, yaxis);
   v2_angle = get_angle(midpoint_axis, v2_vector, yaxis);

   /* first make sure you are both on one side or the other */
   if ((fabs(v1_angle) < angle_range && fabs(v2_angle) > angle_range) ||
       (fabs(v1_angle) > angle_range && fabs(v2_angle) < angle_range)) {
      printf("cusp straddles vertex, no bueno\n");
      exit(ERROR);
   }
   if (flag)
      printf("angle1: %8.1f\n", Rad2Deg * v1_angle);
   if (flag)
      printf("angle2: %8.1f\n", Rad2Deg * v2_angle);
   if (fabs(v1_angle) < angle_range && fabs(v2_angle) < angle_range) {
      if (v2_angle > v1_angle) {
         if (flag)
            printf(" add new cusp 1\n");
         add_new_cusp(cusp_edge, icusp, jcusp, probe, concave_cycle,
                      concave_edge, vertex, concave_circle,
                      probe_rad, cusp_pair, n_cusp_pairs, angle_range,
                      midpoint, c1, xaxis, yaxis, zaxis, flag);
         return;
      } else {
         return;
      }
   }
   if (v1_angle < 0.0) {
      if (flag) {
         printf("v1 < 0, changing to positive angle:");
         printf(" before %8.1f", Rad2Deg * v1_angle);
         printf(" + %8.1f", Rad2Deg * TWOPI);
      }
      v1_angle = TWOPI + v1_angle;
      if (flag)
         printf(" after %8.3f\n", Rad2Deg * v1_angle);
   }
   if (v2_angle < 0.0) {
      if (flag) {
         printf("v2 < 0, changing to positive angle:");
         printf(" before %8.1f", Rad2Deg * v2_angle);
         printf(" + %8.1f", Rad2Deg * TWOPI);
      }
      v2_angle = TWOPI + v2_angle;
      if (flag)
         printf(" after %8.3f\n", Rad2Deg * v2_angle);
   }
   if (v2_angle > v1_angle) {
      if (flag) {
         printf(" add new cusp 2\n");
         printf("v2_angle > v1_angle\n");
         printf("v1 angle %8.1f\n", Rad2Deg * v1_angle);
         printf("v2 angle %8.1f\n", Rad2Deg * v2_angle);
      }
      add_new_cusp(cusp_edge, icusp, jcusp, probe, concave_cycle,
                   concave_edge, vertex, concave_circle,
                   probe_rad, cusp_pair, n_cusp_pairs, angle_range,
                   midpoint, c1, xaxis, yaxis, zaxis, flag);
      return;
   }
   if (flag) {
      printf("no cusp \n");
      printf("v1_angle %8.1f\n", Rad2Deg * v1_angle);
      printf("v2_angle %8.1f\n", Rad2Deg * v2_angle);
      printf("angle range %8.1f\n", Rad2Deg * angle_range);
   }
   return;
}

/*****************************************************************/
static int get_cycle_id(CUSP_EDGE cusp_edge[], int icusp, int jcusp,
                        int flag)
{

   if (cusp_edge[icusp].cycle1 != cusp_edge[jcusp].cycle1 &&
       cusp_edge[icusp].cycle2 != cusp_edge[jcusp].cycle1 &&
       cusp_edge[icusp].cycle1 != cusp_edge[jcusp].cycle2 &&
       cusp_edge[icusp].cycle2 != cusp_edge[jcusp].cycle2)
      return -1;

   if ((cusp_edge[icusp].cycle1 == cusp_edge[jcusp].cycle1 &&
        cusp_edge[icusp].cycle2 == cusp_edge[jcusp].cycle2) ||
       (cusp_edge[icusp].cycle2 == cusp_edge[jcusp].cycle1 &&
        cusp_edge[icusp].cycle1 == cusp_edge[jcusp].cycle2)) {
      printf("cusps have identical cycles\n");
      printf("icusp %d cycles %d %d\n", icusp,
             cusp_edge[icusp].cycle1, cusp_edge[icusp].cycle2);
      printf("jcusp %d cycles %d %d\n", jcusp,
             cusp_edge[jcusp].cycle1, cusp_edge[jcusp].cycle2);
      exit(ERROR);
   }
   if (cusp_edge[icusp].cycle1 == cusp_edge[jcusp].cycle1 ||
       cusp_edge[icusp].cycle1 == cusp_edge[jcusp].cycle2) {
      return cusp_edge[icusp].cycle1;
   } else if (cusp_edge[icusp].cycle2 == cusp_edge[jcusp].cycle1 ||
              cusp_edge[icusp].cycle2 == cusp_edge[jcusp].cycle2) {
      return cusp_edge[icusp].cycle2;
   } else {
      printf("no cycles found\n");
      exit(ERROR);
      return -1;
   }
}

/*****************************************************************/
static void get_probe_id(CUSP_EDGE cusp_edge[], int icusp, int jcusp,
                         CONCAVE_CYCLE concave_cycle[], int *probe1,
                         int *probe2)
{
   if (cusp_edge[icusp].cycle1 == cusp_edge[jcusp].cycle1) {
      *probe1 = concave_cycle[cusp_edge[icusp].cycle2].iprobe;
      *probe2 = concave_cycle[cusp_edge[jcusp].cycle2].iprobe;
   } else if (cusp_edge[icusp].cycle1 == cusp_edge[jcusp].cycle2) {
      *probe1 = concave_cycle[cusp_edge[icusp].cycle2].iprobe;
      *probe2 = concave_cycle[cusp_edge[jcusp].cycle1].iprobe;
   } else if (cusp_edge[icusp].cycle2 == cusp_edge[jcusp].cycle1) {
      *probe1 = concave_cycle[cusp_edge[icusp].cycle1].iprobe;
      *probe2 = concave_cycle[cusp_edge[jcusp].cycle2].iprobe;
   } else if (cusp_edge[icusp].cycle2 == cusp_edge[jcusp].cycle2) {
      *probe1 = concave_cycle[cusp_edge[icusp].cycle1].iprobe;
      *probe2 = concave_cycle[cusp_edge[jcusp].cycle1].iprobe;
   } else {
      printf("get_probeid(): no cycles match\n");
      exit(ERROR);
   }
}

/*****************************************************************/
static void add_new_cusp(CUSP_EDGE cusp_edge[], int icusp, int jcusp,
                         PROBE probe[], CONCAVE_CYCLE concave_cycle[],
                         EDGE concave_edge[], VERTEX vertex[],
                         CIRCLE concave_circle[], REAL_T probe_rad,
                         CUSP_PAIR cusp_pair[], int *n, REAL_T angle_range,
                         POINT midpoint, int c1, POINT xaxis, POINT yaxis,
                         POINT zaxis, int flag)
{
   REAL_T d_pp;
   int p1, p2, ii;
   POINT v1, v2;
   int ie, iv1;
   POINT v_a, v_b, v_c;
   REAL_T angle1, angle2;


   for (ii = 0; ii < 3; ++ii) {
      v1[ii] = midpoint[ii] +
          concave_circle[c1].rad * sin(angle_range) * zaxis[ii];
      v2[ii] = midpoint[ii] -
          concave_circle[c1].rad * sin(angle_range) * zaxis[ii];
   }
   if (flag) {
      printf("angle range %8.1f\n", angle_range * Rad2Deg);
      printf("sin angle %8.3f\n", sin(angle_range));
      printf("radius %8.3f\n", concave_circle[c1].rad);
      printf("zaxis: %8.3f%8.3f%8.3f\n", zaxis[0], zaxis[1], zaxis[2]);

      printf("v1: %8.3f%8.3f%8.3f\n", v1[0], v1[1], v1[2]);
      printf("v2: %8.3f%8.3f%8.3f\n", v2[0], v2[1], v2[2]);
   }
   get_probe_id(cusp_edge, icusp, jcusp, concave_cycle, &p1, &p2);
   d_pp = sqrt((probe[p1].pos[0] - probe[p2].pos[0]) *
               (probe[p1].pos[0] - probe[p2].pos[0]) +
               (probe[p1].pos[1] - probe[p2].pos[1]) *
               (probe[p1].pos[1] - probe[p2].pos[1]) +
               (probe[p1].pos[2] - probe[p2].pos[2]) *
               (probe[p1].pos[2] - probe[p2].pos[2]));

   cusp_pair[*n].cusp1 = icusp;
   cusp_pair[*n].cusp2 = jcusp;
   cusp_pair[*n].circle_rad =
       sqrt(probe_rad * probe_rad - d_pp * d_pp / 4);

   cusp_pair[*n].cycle2 = center_cycle(cusp_edge, icusp, jcusp);
   if (cusp_edge[icusp].cycle1 != cusp_pair[*n].cusp2) {
      cusp_pair[*n].cycle1 = cusp_edge[icusp].cycle1;
   } else {
      cusp_pair[*n].cycle1 = cusp_edge[icusp].cycle2;
   }
   if (cusp_edge[jcusp].cycle1 != cusp_pair[*n].cusp2) {
      cusp_pair[*n].cycle3 = cusp_edge[icusp].cycle1;
   } else {
      cusp_pair[*n].cycle3 = cusp_edge[icusp].cycle2;
   }

   for (ii = 0; ii < 3; ++ii) {
      cusp_pair[*n].circle_center[ii] =
          0.5 * (probe[p1].pos[ii] + probe[p2].pos[ii]);
      cusp_pair[*n].circle_axis[ii] =
          probe[p1].pos[ii] - probe[p2].pos[ii];
   }
   vnorm(cusp_pair[*n].circle_axis, 3);
   cusp_pair[*n].cusp1 = icusp;
   cusp_pair[*n].cusp2 = jcusp;

   /* now make sure you have the verts in the right order */
   ie = cusp_edge[icusp].edge;
   iv1 = concave_edge[ie].vert1;

   for (ii = 0; ii < 3; ++ii) {
      v_a[ii] = vertex[iv1].pos[ii] - concave_circle[c1].center[ii];
      v_b[ii] = v1[ii] - concave_circle[c1].center[ii];
      v_c[ii] = v2[ii] - concave_circle[c1].center[ii];
   }

   angle1 = get_angle(v_a, v_b, yaxis);
   angle2 = get_angle(v_a, v_c, yaxis);
   if (angle1 < 0.0)
      angle1 = angle1 + 2 * PI;
   if (angle2 < 0.0)
      angle2 = angle2 + 2 * PI;

   if (angle2 > angle1) {       /* order is correct */
      for (ii = 0; ii < 3; ++ii) {
         cusp_pair[*n].vert1[ii] = v1[ii];
         cusp_pair[*n].vert2[ii] = v2[ii];
      }
   } else {
      for (ii = 0; ii < 3; ++ii) {
         cusp_pair[*n].vert1[ii] = v2[ii];
         cusp_pair[*n].vert2[ii] = v1[ii];
      }
   }


   ++(*n);
   if (*n >= NUM_CUSP * natm_sel) {
      printf("add_new_cusp(): MAX_CUSP_PAIRS exceeded\n");
      exit(ERROR);
   }
   return;
}

/************************************************************************/
static int center_cycle(CUSP_EDGE cusp[], int icusp, int jcusp)
{
   if (cusp[icusp].cycle1 == cusp[jcusp].cycle1 ||
       cusp[icusp].cycle1 == cusp[jcusp].cycle2) {
      return cusp[icusp].cycle1;
   } else if
       (cusp[icusp].cycle2 == cusp[jcusp].cycle1 ||
        cusp[icusp].cycle2 == cusp[jcusp].cycle2) {
      return cusp[icusp].cycle2;
   } else {
      printf("center_cycle():no cusp match\n");
      exit(ERROR);
      return -1;
   }
}

/************************************************************************/
static int number_of_cusps(CUSP_GROUP group[], int ig,
                           CUSP_PAIR cusp_pair[])
{
   int i, ic, ip;
   int cusps[10];

   i = 0;
   for (ic = 0; ic < group[ig].n_pairs; ++ic) {
      ip = group[ig].cusp_pair[ic];
      if (new_cusp(cusp_pair[ip].cusp1, cusps, i)) {
         cusps[i] = cusp_pair[ip].cusp1;
         ++i;
      }
      if (i > 10) {
         printf("number_of_cusps(): too many unique cusps\n");
      }
      if (new_cusp(cusp_pair[ip].cusp2, cusps, i)) {
         cusps[i] = cusp_pair[ip].cusp2;
         ++i;
      }
      if (i > 10) {
         printf("number_of_cusps(): too many unique cusps\n");
      }
   }
   return i;
}

/*********************************************************************/
static int new_cusp(int icusp, int cusplist[], int ncusps)
{
   int i;
   for (i = 0; i < ncusps; ++i) {
      if (cusplist[i] == icusp)
         return 0;
   }
   return 1;
}


static void trim_2_cusps(PROBE probe[],
                         int *nverts, VERTEX vertex[],
                         int *n_concave_edges, EDGE concave_edge[],
                         int *n_concave_circles, CIRCLE concave_circle[],
                         REAL_T probe_rad,
                         int *n_broken_concave_faces,
                         BROKEN_CONCAVE_FACE broken_concave_face[],
                         CONCAVE_CYCLE concave_cycle[],
                         int *n_concave_cycles, CUSP_EDGE cusp_edge[],
                         int *n_cusps, CUSP_PAIR cusp_pair[],
                         int *n_cusp_pairs, CUSP_GROUP group[], int ig)
{
   int i, ip;
   int two_cusp[2];
   int trim_face[MAXTMP], nfaces;
   EXTREME_VERTEX xvertex[2];
   int icycle, jcycle;
   int cusp_start, cusp_stop;   /*pointers to begging and end of new cusps added */
   int iprobe, jprobe;

   int iface;
   int cusp_list[MAXTMP], n_list;

   ip = group[ig].cusp_pair[0];
   two_cusp[0] = cusp_pair[ip].cusp1;
   two_cusp[1] = cusp_pair[ip].cusp2;

   /* kill the intersecting cusps and associated edges */

   for (i = 0; i < 2; ++i) {
      cusp_edge[two_cusp[i]].alive = 0;
      concave_edge[cusp_edge[two_cusp[i]].edge].alive = 0;
   }

   xvertex[0].cusp_pair = ip;
   xvertex[0].vert = 1;
   xvertex[1].cusp_pair = ip;
   xvertex[1].vert = 2;

   unique_cycles(&cusp_edge[two_cusp[0]],
                 &cusp_edge[two_cusp[1]], &icycle, &jcycle);
   iprobe = concave_cycle[icycle].iprobe;
   jprobe = concave_cycle[jcycle].iprobe;

   cusp_start = *n_cusps;

   add_cusp_circle(n_concave_circles, concave_circle, probe,
                   iprobe, jprobe, probe_rad);
   add_cusp_verts(nverts, vertex, cusp_pair, xvertex, 0, 1,
                  concave_circle, (*n_concave_circles) - 1);
   add_edge(n_concave_edges, concave_edge, (*nverts) - 2,
            (*nverts) - 1, (*n_concave_circles) - 1,
            vertex, concave_circle);
   add_non_axial_cusp(n_cusps, cusp_edge, icycle, jcycle,
                      iprobe, jprobe, (*n_concave_edges) - 1);

   split_old_cusps(two_cusp, 2, cusp_edge, n_cusps, xvertex, 2,
                   cusp_pair, vertex, n_concave_edges,
                   concave_edge, concave_circle);
   cusp_stop = *n_cusps;

   get_faces(trim_face, &nfaces, two_cusp, 2, cusp_edge,
             concave_cycle, broken_concave_face);
   n_list = cusp_stop - cusp_start;
   if (n_list > MAXTMP) {
      printf("MAXTMP exceeded\n");
      exit(ERROR);
   }
   for (i = 0; i < n_list; ++i)
      cusp_list[i] = cusp_start + i;

   for (i = 0; i < nfaces; ++i) {
      iface = trim_face[i];
      if (broken_concave_face[iface].n_cycles != 1) {
         printf("broken concave face num of cycles != 1\n");
         exit(ERROR);
      }
      split_face(iface, n_broken_concave_faces, broken_concave_face,
                 n_concave_cycles, concave_cycle,
                 cusp_edge, cusp_list, n_list,
                 n_concave_edges, concave_edge, vertex, concave_circle);
   }
   return;
}


/******************************************************/
static void unique_cycles(CUSP_EDGE * cusp1, CUSP_EDGE * cusp2,
                          int *icycle, int *jcycle)
{
   int c1, c2, c3, c4;

   c1 = cusp1->cycle1;
   c2 = cusp1->cycle2;
   c3 = cusp2->cycle1;
   c4 = cusp2->cycle2;

   if (c1 != c3 && c1 != c4) {
      *icycle = c1;
   } else if (c2 != c3 && c2 != c4) {
      *icycle = c2;
   } else {
      printf("no unique cycle found\n");
   }

   if (c3 != c1 && c3 != c2) {
      *jcycle = c3;
   } else if (c4 != c1 && c4 != c2) {
      *jcycle = c4;
   } else {
      printf("no unique cycle found\n");
   }
   return;
}


static void split_old_cusps(int cusp_list[], int nlist,
                            CUSP_EDGE cusp_edge[], int *n_cusps,
                            EXTREME_VERTEX xvertex[], int nx,
                            CUSP_PAIR cusp_pair[], VERTEX vertex[],
                            int *n_concave_edges, EDGE concave_edge[],
                            CIRCLE concave_circle[])
{

   int cut_vert[MAXTMP], icut, ncut, i, ix, icusp, ip;
   REAL_T cut_angle[2];
   int iva, ivb, iedge, icircle, ii, itmp;
   POINT va, vb;

   for (i = 0; i < nlist; ++i) {
      icusp = cusp_list[i];
      iedge = cusp_edge[icusp].edge;
      icircle = concave_edge[iedge].circle;
      ncut = 0;
      for (ix = 0; ix < nx; ++ix) {
         ip = xvertex[ix].cusp_pair;
         if (cusp_pair[ip].cusp1 == icusp || cusp_pair[ip].cusp2 == icusp) {
            cut_vert[ncut] = xvertex[ix].vert_index;
            ncut = ncut + 1;
            if (ncut >= MAXTMP) {
               printf("MAXTMP exceeded\n");
               exit(ERROR);
            }
         }
      }
      if (ncut != 2) {
         printf("split_old_cusps: not cutting with 2 verts\n");
         printf("ncut %d\n", ncut);
         exit(ERROR);
      }
      iva = concave_edge[iedge].vert1;
      for (ii = 0; ii < 3; ++ii)
         va[ii] = vertex[iva].pos[ii] - concave_circle[icircle].center[ii];
      for (icut = 0; icut < ncut; ++icut) {
         ivb = cut_vert[icut];
         for (ii = 0; ii < 3; ++ii)
            vb[ii] =
                vertex[ivb].pos[ii] - concave_circle[icircle].center[ii];
         cut_angle[icut] = get_angle(vb, va, concave_circle[icircle].axis);
         if (cut_angle[icut] < 0.0)
            cut_angle[icut] = cut_angle[icut] + TWOPI;
      }

      if (cut_angle[0] > cut_angle[1]) {
         itmp = cut_vert[0];
         cut_vert[0] = cut_vert[1];
         cut_vert[1] = itmp;
      }
      /* now we are ready to split the cusp */
      add_edge(n_concave_edges, concave_edge,
               concave_edge[iedge].vert1, cut_vert[0],
               icircle, vertex, concave_circle);
      make_new_cusp(n_cusps, cusp_edge, icusp, (*n_concave_edges) - 1);

      add_edge(n_concave_edges, concave_edge,
               cut_vert[1], concave_edge[iedge].vert2,
               icircle, vertex, concave_circle);
      make_new_cusp(n_cusps, cusp_edge, icusp, (*n_concave_edges) - 1);
   }
   return;
}

static void make_new_cusp(int *nc, CUSP_EDGE cusp_edge[], int iold,
                          int new_edge)
{
   if (cusp_edge[iold].alive) {
      printf("new_cusp(): old cusp not dead\n");
      exit(ERROR);
   }
   cusp_edge[*nc].cycle1 = cusp_edge[iold].cycle1;
   cusp_edge[*nc].cycle2 = cusp_edge[iold].cycle2;
   cusp_edge[*nc].edge = new_edge;
   cusp_edge[*nc].probe1 = cusp_edge[iold].probe1;
   cusp_edge[*nc].probe2 = cusp_edge[iold].probe2;
   cusp_edge[*nc].alive = 1;
   cusp_edge[*nc].concentric_pair = 0;

   ++(*nc);

   if (*nc > NUM_CUSP * natm_sel) {
      printf("MAX_CUSPS exceeded\n");
      exit(ERROR);
   }
   return;
}
static void get_faces(int face_list[], int *nfaces, int cusp_list[],
                      int nlist, CUSP_EDGE cusp_edge[],
                      CONCAVE_CYCLE concave_cycle[],
                      BROKEN_CONCAVE_FACE broken_concave_face[])
{
   int i, icusp, cycle[2], face[2], j;
   int nf;

   nf = 0;
   for (i = 0; i < nlist; ++i) {
      icusp = cusp_list[i];
      cycle[0] = cusp_edge[icusp].cycle1;
      face[0] = concave_cycle[cycle[0]].iface;
      cycle[1] = cusp_edge[icusp].cycle2;
      face[1] = concave_cycle[cycle[1]].iface;

      for (j = 0; j < 2; ++j) {
         if (broken_concave_face[face[j]].n_cycles != 1) {
            printf("broken concave face has more than one cycle\n");
            printf("face %d number of cycles %d\n", face[j],
                   broken_concave_face[face[j]].n_cycles);
            exit(ERROR);
         }
         if (broken_concave_face[face[j]].concave_cycle[0] != cycle[j]) {
            printf("face cycle mismatch\n");
            printf("face %d face.cycle = %d  cycle = %d\n",
                   face[j],
                   broken_concave_face[face[j]].concave_cycle[0],
                   cycle[j]);
            exit(ERROR);
         }
         if (is_new_face(face[j], face_list, nf)) {
            face_list[nf] = face[j];
            ++nf;
            if (nf > MAXTMP) {
               printf("MAXTMP exceeded\n");
               exit(ERROR);
            }
         }
      }
   }
   *nfaces = nf;
   return;
}

/*********************************************/

static int is_new_face(int iface, int four_face[], int nf)
{
   int i;
   for (i = 0; i < nf; ++i)
      if (iface == four_face[i])
         return 0;
   return 1;
}

/************************************************************/
static void split_face(int iface,
                       int *n_broken_faces,
                       BROKEN_CONCAVE_FACE broken_concave_face[],
                       int *n_concave_cycles,
                       CONCAVE_CYCLE concave_cycle[],
                       CUSP_EDGE cusp_edge[], int cusp_list[], int n_list,
                       int *n_concave_edges, EDGE concave_edge[],
                       VERTEX vertex[], CIRCLE concave_circle[])
{
   int icycle, ie;
   int start[MAXTMP], ns;       /* starting edges for the 2 subcycles */
   CONCAVE_CYCLE tmp_cycle;
   int cusp_used[MAXTMP];
   int edge_used[MAXTMP];
   int i, direction;
   int firstvert, nextvert;
   int icusp;
   int two_cycle[2];
   int jj;
   int nnc;
   int narc, npt;
   int ncirc;

   if (n_list >= MAXTMP) {
      printf("split_face(): MAXTMP exceeded\n");
      exit(ERROR);
   }
   narc = 0;
   npt = 0;
   ncirc = 0;


   for (i = 0; i < n_list; ++i)
      cusp_used[i] = 0;

   icycle = broken_concave_face[iface].concave_cycle[0];
   two_cycle[0] = icycle;
   two_cycle[1] = *n_concave_cycles;


   /* copy cycle to tmp_cycle */

   copy_cycle(&concave_cycle[icycle], &tmp_cycle);

   /*
      dump_cycle(&tmp_cycle, concave_edge);
    */
   if (tmp_cycle.nedges > MAXTMP) {
      printf("MAXTMP exceeded\n");
      exit(ERROR);
   }
   for (i = 0; i < tmp_cycle.nedges; ++i)
      edge_used[i] = 0;

   /* find the starting edges for the subcycles */
   get_starters(start, &ns, concave_cycle, icycle, concave_edge);

   /* the first cycle will go right back in the original
      cycle spot. the 2nd one will be new */

   for (jj = 0; jj < ns; ++jj) {
      ie = start[jj];
      i = 0;
      concave_cycle[two_cycle[jj]].edge[i] = tmp_cycle.edge[ie];
      concave_cycle[two_cycle[jj]].edge_direction[i] =
          tmp_cycle.edge_direction[ie];
      concave_cycle[two_cycle[jj]].cusp_edge[i] = tmp_cycle.cusp_edge[ie];
      ++i;
      firstvert = concave_edge[tmp_cycle.edge[ie]].vert1;
      nextvert = concave_edge[tmp_cycle.edge[ie]].vert2;

      while (nextvert != firstvert) {
         /* 1st look through "normal" edges for next vertex 
            also include cusp edges that are not dead */
         ie = normal_match(tmp_cycle, concave_edge, nextvert, edge_used);
         if (ie != -1) {        /* normal cusp was found */
            concave_cycle[two_cycle[jj]].edge[i] = tmp_cycle.edge[ie];
            concave_cycle[two_cycle[jj]].edge_direction[i] =
                tmp_cycle.edge_direction[ie];
            concave_cycle[two_cycle[jj]].cusp_edge[i] =
                tmp_cycle.cusp_edge[ie];
         } else {               /* look through new cusps */
            icusp = cusp_match(nextvert, n_list, cusp_list, icycle,
                               cusp_used, cusp_edge, concave_edge,
                               &direction);
            concave_cycle[two_cycle[jj]].edge[i] = cusp_edge[icusp].edge;
            concave_cycle[two_cycle[jj]].edge_direction[i] = direction;
            concave_cycle[two_cycle[jj]].cusp_edge[i] = icusp;
         }

         if (concave_cycle[two_cycle[jj]].edge_direction[i] == 1)
            nextvert =
                concave_edge[concave_cycle[two_cycle[jj]].edge[i]].vert2;
         else
            nextvert =
                concave_edge[concave_cycle[two_cycle[jj]].edge[i]].vert1;
         ++i;
      }
      concave_cycle[two_cycle[jj]].nedges = i;
   }
   if (ns == 2) {
      nnc = broken_concave_face[iface].n_cycles;

      broken_concave_face[*n_broken_faces].itorus[0] =
          broken_concave_face[iface].itorus[0];
      broken_concave_face[*n_broken_faces].itorus[1] =
          broken_concave_face[iface].itorus[1];
      broken_concave_face[*n_broken_faces].itorus[2] =
          broken_concave_face[iface].itorus[2];
      broken_concave_face[*n_broken_faces].probe =
          broken_concave_face[iface].probe;
      broken_concave_face[*n_broken_faces].n_cycles = 1;
      broken_concave_face[*n_broken_faces].concave_cycle[0] =
          *n_concave_cycles;
      broken_concave_face[*n_broken_faces].alive = 1;
      broken_concave_face[*n_broken_faces].area = 0.0;

      /* reset n_cycle count to 1 for iface */
      broken_concave_face[iface].n_cycles = 1;

      ++(*n_concave_cycles);
      if (*n_concave_cycles > NUM_CYCLE * natm_sel) {
         printf("MAX_CYCLES exceeded\n");
         exit(ERROR);
      }
      ++(*n_broken_faces);
      if (*n_broken_faces > NUM_FACE * natm_sel) {
         printf("MAX_FACE exceeded\n");
         exit(ERROR);
      }
   } else if (ns > 2) {
      printf("split_face(): face has more than 2 cycles\n");
      exit(ERROR);
   }
   /*
      printf(" remade icycle %d\n", icycle);
      dump_cycle(&concave_cycle[icycle], concave_edge);
      printf(" new cycle %d\n", *n_concave_cycles - 1);
      dump_cycle(&concave_cycle[*n_concave_cycles - 1], concave_edge);
    */

   return;
}

/************************************************************/
static void get_starters(int start[], int *ns,
                         CONCAVE_CYCLE concave_cycle[], int icycle,
                         EDGE concave_edge[])
{
   int j;

   *ns = 0;
   for (j = 0; j < concave_cycle[icycle].nedges; ++j) {
      if (concave_cycle[icycle].cusp_edge[j] != -1 &&
          concave_edge[concave_cycle[icycle].edge[j]].alive == 0) {
         /* next edge is starting edge */
         if (j == concave_cycle[icycle].nedges - 1)
            start[*ns] = 0;
         else
            start[*ns] = j + 1;
         ++(*ns);
         if (*ns > MAXTMP) {
            printf("too many starting edges\n");
            exit(ERROR);
         }
      }
   }
   if (*ns > 2) {
      printf("split_face(): num starting edges  > 2  (%d)\n", *ns);
      exit(ERROR);
   }
}

/************************************************************/
static void copy_cycle(CONCAVE_CYCLE * c1, CONCAVE_CYCLE * c2)
{
   int i;

   c2->nedges = c1->nedges;
   c2->iprobe = c1->iprobe;
   c2->iface = c1->iface;
   c2->intersects_self = c1->intersects_self;

   for (i = 0; i < c2->nedges; ++i) {
      c2->edge[i] = c1->edge[i];
      c2->edge_direction[i] = c1->edge_direction[i];
      c2->cusp_edge[i] = c1->cusp_edge[i];
   }
   return;
}

/************************************************************/
static int cusp_match(int lastvert, int ncusps, int cusp_index[],
                      int icycle, int cusp_used[], CUSP_EDGE cusp_edge[],
                      EDGE concave_edge[], int *direction)
{
   int i, icusp, ie;

   for (i = 0; i < ncusps; ++i) {
      icusp = cusp_index[i];
      if (cusp_used[i])
         continue;
      if (cusp_edge[icusp].cycle1 != icycle &&
          cusp_edge[icusp].cycle2 != icycle)
         continue;
      ie = cusp_edge[icusp].edge;
      if (concave_edge[ie].vert1 == lastvert) {
         *direction = 1;
         cusp_used[i] = 1;
         return icusp;
      }
      if (concave_edge[ie].vert2 == lastvert) {
         *direction = -1;
         cusp_used[i] = 1;
         return icusp;
      }
   }

   printf("cusp_match(): could not find match for vertex %d\n", lastvert);
   exit(ERROR);
   return -1;
}

/**************************************************/
static int normal_match(CONCAVE_CYCLE tmp_cycle, EDGE concave_edge[],
                        int lastvert, int edge_used[])
{
   int i, ie;

   for (i = 0; i < tmp_cycle.nedges; ++i) {
      if (edge_used[i])
         continue;
      ie = tmp_cycle.edge[i];
      if (concave_edge[ie].alive == 0)
         continue;
      if (concave_edge[ie].vert1 == lastvert ||
          concave_edge[ie].vert2 == lastvert) {
         edge_used[i] = 1;
         return i;
      }
   }
   return -1;
}

/***************************************************************/
static void add_non_axial_cusp(int *nc, CUSP_EDGE cusp_edge[],
                               int icycle, int jcycle, int iprobe,
                               int jprobe, int iedge)
{
   cusp_edge[*nc].cycle1 = icycle;
   cusp_edge[*nc].cycle2 = jcycle;
   cusp_edge[*nc].probe1 = iprobe;
   cusp_edge[*nc].probe2 = jprobe;
   cusp_edge[*nc].edge = iedge;
   cusp_edge[*nc].alive = 1;
   cusp_edge[*nc].concentric_pair = 0;
   ++(*nc);
   if (*nc > NUM_CUSP * natm_sel) {
      printf("MAX_CUSPS exceeded\n");
      exit(ERROR);
   }
   return;
}

/***************************************************************/
static void add_1_vert(int *nv, VERTEX vertex[], POINT v)
{
   int ii;

   for (ii = 0; ii < 3; ++ii)
      vertex[*nv].pos[ii] = v[ii];
   ++(*nv);
   if (*nv > NUM_VERTEX * natm_sel) {
      printf("MAX_VERTS exceeded\n");
      exit(ERROR);
   }
   return;
}
static void add_cusp_verts(int *nverts, VERTEX vertex[],
                           CUSP_PAIR cusp_pair[], EXTREME_VERTEX xvertex[],
                           int ix, int jx, CIRCLE concave_circle[],
                           int icircle)
{
   int ip1, ip2, ii;
   POINT vert1, vert2;

   ip1 = xvertex[ix].cusp_pair;
   if (xvertex[ix].vert == 1) {
      for (ii = 0; ii < 3; ++ii)
         vert1[ii] = cusp_pair[ip1].vert1[ii];
   } else {
      for (ii = 0; ii < 3; ++ii)
         vert1[ii] = cusp_pair[ip1].vert2[ii];
   }

   ip2 = xvertex[jx].cusp_pair;
   if (xvertex[jx].vert == 1) {
      for (ii = 0; ii < 3; ++ii)
         vert2[ii] = cusp_pair[ip2].vert1[ii];
   } else {
      for (ii = 0; ii < 3; ++ii)
         vert2[ii] = cusp_pair[ip2].vert2[ii];
   }

   /* START HERE:  I THINK COMMENTED  CRITERION IS BOGUS! 
      and the one just below is correct */

   add_1_vert(nverts, vertex, vert1);
   xvertex[ix].vert_index = (*nverts) - 1;
   add_1_vert(nverts, vertex, vert2);
   xvertex[jx].vert_index = (*nverts) - 1;
   /*
      for (ii = 0; ii < 3; ++ii) {
      va[ii] = vert1[ii] - concave_circle[icircle].center[ii];
      vb[ii] = vert2[ii] - concave_circle[icircle].center[ii];
      }

      vnorm(va,3);
      vnorm(vb,3);

      cross (va, vb, vc);

      if (DOT(vc, concave_circle[icircle].axis) > 0.0 ) {
      add_1_vert(nverts, vertex, vert1);
      xvertex[ix].vert_index = (*nverts) - 1;
      add_1_vert(nverts, vertex, vert2);
      xvertex[jx].vert_index = (*nverts) - 1;
      } else {
      add_1_vert(nverts, vertex, vert2);
      xvertex[jx].vert_index = (*nverts) - 1;
      add_1_vert(nverts, vertex, vert1);
      xvertex[ix].vert_index = (*nverts) - 1;
      }
    */
   return;
}

/***************************************************************/

/***************************************************************************/
static void add_cusp_circle(int *nc, CIRCLE circle[], PROBE probe[],
                            int ip, int jp, REAL_T probe_rad)
{
   int ii;
   REAL_T d_pp2;

   d_pp2 = 0.0;
   for (ii = 0; ii < 3; ++ii) {
      circle[*nc].center[ii] =
          0.5 * (probe[ip].pos[ii] + probe[jp].pos[ii]);
      /*
         circle[*nc].axis[ii] = probe[ip].pos[ii] - probe[jp].pos[ii];
       */
      circle[*nc].axis[ii] = probe[jp].pos[ii] - probe[ip].pos[ii];
      d_pp2 = d_pp2 + (probe[ip].pos[ii] - probe[jp].pos[ii]) *
          (probe[ip].pos[ii] - probe[jp].pos[ii]);
   }
   vnorm(circle[*nc].axis, 3);
   circle[*nc].rad = sqrt(probe_rad * probe_rad - d_pp2 / 4.0);
   circle[*nc].torus = -1;
   circle[*nc].atom_or_probe_num = -1;
   ++(*nc);
   if (*nc >= NUM_CIRCLE * natm_sel) {
      printf("MAX_CIRCLE exceeded\n");
      exit(ERROR);
   }
   return;
}


static int broken_concave_area(REAL_T probe_rad,
                               int n_broken_concave_faces,
                               BROKEN_CONCAVE_FACE broken_concave_face[],
                               CONCAVE_CYCLE concave_cycle[],
                               EDGE concave_edge[], CIRCLE circle[],
                               VERTEX vertex[], REAL_T * broken_conc_area,
                               PROBE probe[])
{
   int iface, ic, icycle, chi, iprobe;
   REAL_T area, total_area = 0;
   int n_surfaced = 0;

   *broken_conc_area = 0;

   for (iface = 0; iface < n_broken_concave_faces; ++iface) {
      ++n_surfaced;
      chi = 2 - broken_concave_face[iface].n_cycles;
      area = 0;
      iprobe = broken_concave_face[iface].probe;
      for (ic = 0; ic < broken_concave_face[iface].n_cycles; ++ic) {
         icycle = broken_concave_face[iface].concave_cycle[ic];
         area = area + conc_cycle_piece(icycle, concave_cycle,
                                        circle, concave_edge, vertex,
                                        iprobe, probe, probe_rad);
      }
      broken_concave_face[iface].area =
          probe_rad * probe_rad * (2 * PI * chi + area);

      total_area += broken_concave_face[iface].area;
      *broken_conc_area += broken_concave_face[iface].area;
#ifdef DEBUG
      printf
          ("-----\nbroken concave face probe %8.3f%8.3f%8.3f n_cycles %d\n",
           probe[iprobe].pos[0], probe[iprobe].pos[1],
           probe[iprobe].pos[2], broken_concave_face[iface].n_cycles);
      for (ic = 0; ic < broken_concave_face[iface].n_cycles; ++ic) {
         printf("   cycle %d ( ", ic);
         icycle = broken_concave_face[iface].concave_cycle[ic];
         for (int ie = 0; ie < concave_cycle[icycle].nedges; ++ie) {
            int iedge = concave_cycle[icycle].edge[ie];
            if (concave_cycle[icycle].edge_direction[ie] == -1) {
               printf(" !", iedge);
            }
            printf("%d:%d ", concave_edge[iedge].vert1,
                   concave_edge[iedge].vert2);
         }
         printf(" )\n");
      }
      printf(" area: %10.3f\n", broken_concave_face[iface].area);
#endif
   }
   return n_surfaced;;
}



static REAL_T conc_cycle_piece(int ic, CONCAVE_CYCLE concave_cycle[],
                               CIRCLE circle[], EDGE concave_edge[],
                               VERTEX vertex[], int iprobe, PROBE probe[],
                               REAL_T probe_rad)
{
   REAL_T sum = 0, wrap_angle;
   int i, ii, ie, icircle, itorus, iv1, iv2;
   POINT v1, v2;
   REAL_T cos_theta;
   int edge1, edge2, direction1, direction2;
   REAL_T int_angle;
   REAL_T cdist;

   for (i = 0; i < concave_cycle[ic].nedges; ++i) {
      ie = concave_cycle[ic].edge[i];
      icircle = concave_edge[ie].circle;
      itorus = circle[icircle].torus;
      iv1 = concave_edge[ie].vert1;
      iv2 = concave_edge[ie].vert2;

      /* 1st take care of the interior angle part */
      edge1 = ie;
      direction1 = concave_cycle[ic].edge_direction[i];

      if (i < concave_cycle[ic].nedges - 1) {
         edge2 = concave_cycle[ic].edge[i + 1];
         direction2 = concave_cycle[ic].edge_direction[i + 1];
      } else {
         edge2 = concave_cycle[ic].edge[0];
         direction2 = concave_cycle[ic].edge_direction[0];
      }
      int_angle =
          interior_angle(edge1, direction1, edge2, direction2,
                         concave_edge, circle, vertex);

      sum = sum - (PI - int_angle);

      /* next do the saddle wrap stuff. This part is a little confusing and 
         Connolly's paper (J. Appl. Cryst.  16, 548-558 (1983)) is unclear about 
         the sign of theta.  I think he assumes the torus center always lies 
         between atoms, but this is not necessarily so.  You basically want the 
         projection of the vector that points from the atom center to the circle 
         (or vertex, either way) onto the interatomic axis), but the sign of 
         that term is given by the opposite dot product of the d_c vector 
         (atom to circle center) and the circle unit vector.  
         sign = -dot(d_c, c_axis) if the d_c and c_axis are 
         in opposite directions, then the torus center is indeed between 
         the two atom and this term should be added.  If not, the term 
         should be subtracted.  the MSEED paper 
         (Perrot, et al.  J Comp Chem 13, 1-11 (1992)) has it right.
       */

      if (iv1 == -1) {
         if (concave_cycle[ic].nedges != 1) {
            printf("concave_cycle(): vert = -1 but n_edges > 1\n");
            exit(ERROR);
         }
         wrap_angle = 2.0 * PI;
      } else {
         for (ii = 0; ii < 3; ++ii) {
            v1[ii] = vertex[iv1].pos[ii] - circle[icircle].center[ii];
            v2[ii] = vertex[iv2].pos[ii] - circle[icircle].center[ii];
         }
         wrap_angle = get_angle(v2, v1, circle[icircle].axis);
         if (wrap_angle < 0.0)
            wrap_angle += 2.0 * PI;
      }

      /*
         printf("wrap angle %8.2f\n",Rad2Deg*wrap_angle); 
         printf(" circle center %8.3f%8.3f%8.3f ", circle[icircle].center[0], 
         circle[icircle].center[1], circle[icircle].center[2]);
         printf(" circle axis %8.3f%8.3f%8.3f\n", circle[icircle].axis[0],
         circle[icircle].axis[1], circle[icircle].axis[2]);
       */

      cdist = sqrt((circle[icircle].center[0] - probe[iprobe].pos[0]) *
                   (circle[icircle].center[0] - probe[iprobe].pos[0]) +
                   (circle[icircle].center[1] - probe[iprobe].pos[1]) *
                   (circle[icircle].center[1] - probe[iprobe].pos[1]) +
                   (circle[icircle].center[2] - probe[iprobe].pos[2]) *
                   (circle[icircle].center[2] - probe[iprobe].pos[2])
          );

      cos_theta = cdist / probe_rad;
      sum += wrap_angle * cos_theta;
   }
   /* printf("contribution from edge %f\n", sum); */
   return sum;
}

/**********************************************************************/
static REAL_T interior_angle(int e1, int dir1, int e2, int dir2,
                             EDGE concave_edge[], CIRCLE circle[],
                             VERTEX vertex[])
{
   int c1, c2, ii;
   POINT v1, v2, zaxis, vtmp1, vtmp2;

   c1 = concave_edge[e1].circle;
   c2 = concave_edge[e2].circle;

   /* e1 and e2 share a vertex: need to find the two tangent
      vectors at this vertex. the common vertex should 
      be the 2nd vertex of e1 and the 1st vertex of e2,
      unless dir1 or dir2 are -1, which means the
      directions are reversed */

   if (dir1 > 0) {              /* normal direction edge */
      for (ii = 0; ii < 3; ++ii)
         vtmp1[ii] =
             vertex[concave_edge[e1].vert2].pos[ii] -
             circle[c1].center[ii];
      vnorm(vtmp1, 3);
      cross(circle[c1].axis, vtmp1, v1);
   } else {
      for (ii = 0; ii < 3; ++ii)
         vtmp1[ii] =
             vertex[concave_edge[e1].vert1].pos[ii] -
             circle[c1].center[ii];
      vnorm(vtmp1, 3);
      cross(vtmp1, circle[c1].axis, v1);
   }

   if (dir2 > 0) {              /* normal direction edge */
      for (ii = 0; ii < 3; ++ii)
         vtmp2[ii] =
             vertex[concave_edge[e2].vert1].pos[ii] -
             circle[c2].center[ii];
      vnorm(vtmp2, 3);
      cross(vtmp2, circle[c2].axis, v2);
   } else {
      for (ii = 0; ii < 3; ++ii)
         vtmp2[ii] =
             vertex[concave_edge[e2].vert2].pos[ii] -
             circle[c2].center[ii];
      vnorm(vtmp2, 3);
      cross(circle[c2].axis, vtmp2, v2);
   }

   vnorm(v1, 3);
   vnorm(v2, 3);

   /* angle between v1 and v2 is always < 180,
      so no ambiguity in cross product */
   cross(v1, v2, zaxis);
   vnorm(zaxis, 3);

   return get_angle(v2, v1, zaxis);

}

static void allocate_memory(void)
{
   if ((atom = (ATOM *) malloc(natm_sel * sizeof(ATOM))) == NULL) {
      fprintf(stderr, "Unable to allocate space for atom: %d\n", natm_sel);
      exit(1);
   }
   if ((upper_neighbors =
        (NEIGHBOR_TORUS *) malloc(NUM_NEIGHBOR * natm_sel *
                                  sizeof(NEIGHBOR_TORUS))) == NULL) {
      fprintf(stderr, "Unable to allocate space for upper_neighbors\n");
      exit(1);
   }
   if ((neighbors =
        (NEIGHBOR *) malloc(NUM_NEIGHBOR * natm_sel * sizeof(NEIGHBOR))) ==
       NULL) {
      fprintf(stderr, "Unable to allocate space for neighbors\n");
      exit(1);
   }
   if ((probelist =
        (PROBE *) malloc(NUM_PROBE * natm_sel * sizeof(PROBE))) == NULL) {
      fprintf(stderr, "Unable to allocate space for probelist\n");
      exit(1);
   }
   if ((toruslist =
        (TORUS *) malloc(NUM_TORUS * natm_sel * sizeof(TORUS))) == NULL) {
      fprintf(stderr, "Unable to allocate space for toruslist\n");
      exit(1);
   }
   if ((convex_circle_list =
        (CIRCLE *) malloc(NUM_CIRCLE * natm_sel * sizeof(CIRCLE))) ==
       NULL) {
      fprintf(stderr, "Unable to allocate space for convex_circle_list\n");
      exit(1);
   }
   if ((concave_circle_list =
        (CIRCLE *) malloc(NUM_CIRCLE * natm_sel * sizeof(CIRCLE))) ==
       NULL) {
      fprintf(stderr,
              "Unable to allocate space for concave_circle_list\n");
      exit(1);
   }
   if ((concave_face =
        (CONCAVE_FACE *) malloc(NUM_FACE * natm_sel *
                                sizeof(CONCAVE_FACE))) == NULL) {
      fprintf(stderr, "Unable to allocate space for concave_face\n");
      exit(1);
   }
   if ((convex_face =
        (CONVEX_FACE *) malloc(NUM_FACE * natm_sel *
                               sizeof(CONVEX_FACE))) == NULL) {
      fprintf(stderr, "Unable to allocate space for convex_face\n");
      exit(1);
   }
   if ((saddle_face =
        (SADDLE_FACE *) malloc(NUM_FACE * natm_sel *
                               sizeof(SADDLE_FACE))) == NULL) {
      fprintf(stderr, "Unable to allocate space for saddle_face\n");
      exit(1);
   }
   if ((cone_face =
        (CONE_FACE *) malloc(NUM_FACE * natm_sel * sizeof(CONE_FACE))) ==
       NULL) {
      fprintf(stderr, "Unable to allocate space for cone_face\n");
      exit(1);
   }
   if ((broken_concave_face =
        (BROKEN_CONCAVE_FACE *) malloc(NUM_FACE * natm_sel *
                                       sizeof(BROKEN_CONCAVE_FACE))) ==
       NULL) {
      fprintf(stderr,
              "Unable to allocate space for broken_concave_face\n");
      exit(1);
   }
   if ((concave_cycle =
        (CONCAVE_CYCLE *) malloc(NUM_CYCLE * natm_sel *
                                 sizeof(CONCAVE_CYCLE))) == NULL) {
      fprintf(stderr, "Unable to allocate space for concave_cycle\n");
      exit(1);
   }
   if ((cyclelist =
        (CYCLE *) malloc(NUM_CYCLE * natm_sel * sizeof(CYCLE))) == NULL) {
      fprintf(stderr, "Unable to allocate space for cyclelist\n");
      exit(1);
   }
   if ((vertexlist =
        (VERTEX *) malloc(NUM_VERTEX * natm_sel * sizeof(VERTEX))) ==
       NULL) {
      fprintf(stderr, "Unable to allocate space for vertexlist\n");
      exit(1);
   }
   if ((concave_edge_list =
        (EDGE *) malloc(NUM_EDGE * natm_sel * sizeof(EDGE))) == NULL) {
      fprintf(stderr, "Unable to allocate space for concave_edge_list\n");
      exit(1);
   }
   if ((convex_edge_list =
        (EDGE *) malloc(NUM_EDGE * natm_sel * sizeof(EDGE))) == NULL) {
      fprintf(stderr, "Unable to allocate space for convex_edge_list\n");
      exit(1);
   }
   if ((low_torus =
        (LOW_TORUS *) malloc(NUM_TORUS * natm_sel * sizeof(LOW_TORUS))) ==
       NULL) {
      fprintf(stderr, "Unable to allocate space for low_torus\n");
      exit(1);
   }
   if ((cusp_edge =
        (CUSP_EDGE *) malloc(NUM_EDGE * natm_sel * sizeof(CUSP_EDGE))) ==
       NULL) {
      fprintf(stderr, "Unable to allocate space for cusp_edge\n");
      exit(1);
   }
   if ((cusp_pair =
        (CUSP_PAIR *) malloc(NUM_CUSP * natm_sel * sizeof(CUSP_PAIR))) ==
       NULL) {
      fprintf(stderr, "Unable to allocate space for cusp_pair\n");
      exit(1);
   }
}

static void free_memory(void)
{
   free(atom);
   free(upper_neighbors);
   free(neighbors);
   free(probelist);
   free(toruslist);
   free(convex_circle_list);
   free(concave_circle_list);
   free(concave_face);
   free(convex_face);
   free(saddle_face);
   free(cone_face);
   free(broken_concave_face);
   free(concave_cycle);
   free(cyclelist);
   free(vertexlist);
   free(concave_edge_list);
   free(convex_edge_list);
   free(low_torus);
   free(cusp_edge);
   free(cusp_pair);
   free(res);
}
