#ifndef DEM_h
#define DEM_h

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct dem_collision{
    int idc;
    double dem_overlapn, dem_overlapt;
    struct dem_collision *next;
} dem_collision_t;

typedef struct dem_particle{
    int idp;
    double u, v, w;
    double x, y, z;
    double mp, Ip;
    double diameter, rho;
    double omega_x, omega_y, omega_z;
    double xmax, xmin, ymax, ymin, zmax, zmin;
    double Ep, poisson, coef_rest, coef_fric;
    double Fcontact_x, Fcontact_y, Fcontact_z, Mcontact_x, Mcontact_y, Mcontact_z;
    int contact_dim;
    dem_collision_t *LIST;
    struct dem_particle *next;
} dem_particle_t;


void dem_init_part(dem_particle_t *root);
void dem_init_coll(dem_collision_t *LIST);
void dem_add_particles(dem_particle_t *root, dem_collision_t *LIST, int id_part, double x0, double y0, double z0, double u0, double v0, double w0, double omega_x0, double omega_y0, double omega_z0);
dem_particle_t *search_particle(dem_particle_t *root, int id_part);
void add_contact(dem_collision_t *LIST, int id_part);
void delete_contact(dem_collision_t *LIST, int aux);
void delete_all_contacts(dem_collision_t *LIST);
int search_contact(dem_collision_t *LIST, int id_part);
void dem_bounding_box(dem_particle_t *root);
void dem_overlap_save(dem_collision_t *LIST, int type, int aux, double overlap);
double dem_overlap_value(dem_collision_t *LIST, int type, int aux);
int dem_testing_contacts(dem_particle_t *root,int coll);
int dem_collision(dem_particle_t *root, int coll);
void dem_contact_results(dem_particle_t *root);
void dem_adv(dem_particle_t *root, int coll);
void dem_print_vtk(dem_particle_t *root);

#endif /* DEM_h */
