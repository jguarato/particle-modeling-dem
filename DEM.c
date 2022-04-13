//
//  Computational code for modeling spherical particles
//  considering interparticle collision based on 
//  the Discrete Element Method (DEM).
//
//  The initial interest for the development of this
//  code was to apply it in multiphase flows.
//
//  The setup of the present case corresponds to
//  two particles, one falling on top of the other.
//
//  Created by: Jessica Guarato
//  Last modified on: September, 2017
//

#include "DEM.h"

// ==========================================================================================================
// SETUP
// ==========================================================================================================

// Gravity
double gx = 0;
double gy = -9.81;
double gz = 0;

// Particles properties
int n_particles = 0;                    // Initial number of particles
int id_iter = 0;                        // Initial particle id
int np_iter = 2;                        // Number of particles per iteration
double d_particle = 0.068;              // Diameter
double rho = 1768;                      // Denisty
double E_modulus = 5.4e9;               // Young Modulus
double poisson = 0.34;                  // Poisson
double coef_rest = 1;                   // Restitution coefficient
double coef_fric = 0.40;                // Dynamic friction coefficient between particles

// Domain
double xd1 = 0, xd2 = 1, xd3 = 0, xd4 = 1;
double yd1 = 0, yd2 = 2, yd3 = 0, yd4 = 2;
double zd1 = 0, zd2 = 0, zd3 = 0, zd4 = 0;

// Time parameters
double t = 0;
double tfinal = 2;
double dt = 0.00001;
int iter = 1;

// Output files
int print_step   = 100;                 // Saving frequency of results
char output[20]  = "output"; 


// ==========================================================================================================
// PROGRAM
// ==========================================================================================================
int main() {
    
    dem_particle_t *ptr_dem;
    dem_particle_t *dem_particles = (dem_particle_t *) malloc(sizeof(dem_particle_t));
    
    int coll, id;
    double x0, y0, z0, u0, v0, w0, omega_x0, omega_y0, omega_z0;
    char command1[50], command2[50];


    // Creating output directory
    sprintf(command1, "rm -rf %s",output);
    system(command1);
    sprintf(command2, "mkdir -p %s",output);
    system(command2);
    
    dem_init_part(dem_particles);

    // Hard-coded initial conditions setup
    n_particles = n_particles+np_iter;
    x0 = -0.03; y0 = 2;
    v0 = -10;

    for (id=0;id<n_particles;id++) {
        dem_collision_t *ptr_collision = (dem_collision_t *) malloc(sizeof(dem_collision_t));
        dem_init_coll(ptr_collision);
        
        // Hard-coded initial conditions setup
        x0 = x0+0.03; y0 = y0-0.5; z0 = 0;
        u0 = 0; v0 = v0+5; w0 = 0;
        omega_x0 = 0; omega_y0 = 0; omega_z0 = 0;
        
        dem_add_particles(dem_particles,ptr_collision,id,x0,y0,z0,u0,v0,w0,omega_x0,omega_y0,omega_z0);
    }
    id_iter = id_iter+np_iter;
    
    
    while (t <= tfinal) {
        if(iter%print_step == 0) {
            
            printf("------------------------------\n");
            printf("time iter        = %d\n",iter);
            printf("time step        = %.5e\n",dt);
            printf("elapsed time     = %.5e\n",t);
            printf("------------------------------\n\n");
        }
        
        coll = 0;

        
        for (id=0;id<n_particles;id++) {
            ptr_dem = search_particle(dem_particles,id);
            dem_bounding_box(ptr_dem);
            
        }
        
        coll = dem_testing_contacts(dem_particles,coll);

        if (coll > 0) {
            coll = dem_collision(dem_particles,coll);
        }
        
        if (coll > 0) {
            dem_contact_results(dem_particles);
        }
        
        for (id=0;id<n_particles;id++) {
            ptr_dem = search_particle(dem_particles,id);
            dem_adv(ptr_dem,coll);

        }
        
        if(iter%print_step == 0) {
            dem_print_vtk(dem_particles);
        }
        
        
        t = t+dt;
        iter = iter+1;
        
    }
    
    return 0;
}


void dem_init_part(dem_particle_t *root) {
    root->next = NULL;
}

void dem_init_coll(dem_collision_t *LIST) {
    LIST->next = NULL;
}

void dem_add_particles(dem_particle_t *root, dem_collision_t *LIST, int id_part, 
    double x0, double y0, double z0, double u0, double v0, double w0, 
    double omega_x0, double omega_y0, double omega_z0) {
    
    double mp, Ip;
    mp = rho*((1./6.)*M_PI*pow(d_particle,3));
    Ip = (1./10.)*mp*pow(d_particle,2);
    
    dem_particle_t *new = (dem_particle_t *) malloc(sizeof(dem_particle_t));
    new->idp = id_part;
    new->x = x0;
    new->y = y0;
    new->z = z0;
    new->u = u0;
    new->v = v0;
    new->w = w0;
    new->omega_x = omega_x0;
    new->omega_y = omega_y0;
    new->omega_z = omega_z0;
    new->diameter = d_particle;
    new->rho = rho;
    new->Ep = E_modulus;
    new->poisson = poisson;
    new->coef_rest = coef_rest;
    new->coef_fric = coef_fric;
    new->Fcontact_x = 0;
    new->Fcontact_y = 0;
    new->Fcontact_z = 0;
    new->Mcontact_x = 0;
    new->Mcontact_y = 0;
    new->Mcontact_z = 0;
    new->contact_dim = 0;
    new->mp = mp;
    new->Ip = Ip;
    new->LIST = LIST;
    
    new->next = NULL;
    
    if(root->next == NULL) {
        root->next = new;
    }
    
    else{
        dem_particle_t *tmp = root->next;
        
        while(tmp->next != NULL) {
            tmp = tmp->next;
        }
        
        tmp->next = new;
    }

}


dem_particle_t *search_particle(dem_particle_t *root, int id_part) {
    
    dem_particle_t *tmp = root->next;
    
    while(tmp != NULL) {
        
        if(tmp->idp == id_part) {
            break;
        }
        
        tmp = tmp->next;
    }
    
    return tmp;
}


void add_contact(dem_collision_t *LIST, int id_part) {
    
    dem_collision_t *new = (dem_collision_t *) malloc(sizeof(dem_collision_t));
    
    new->idc = id_part;
    new->dem_overlapn = 0.0;
    new->dem_overlapt = 0.0;
    
    new->next = NULL;
    
    if(LIST->next == NULL) {
        LIST->next = new;
    }
    
    else{
        dem_collision_t *tmp = LIST->next;
        
        while(tmp->next != NULL) {
            tmp = tmp->next;
        }
        
        tmp->next = new;
    }
    
}


void delete_contact(dem_collision_t *LIST, int aux) {
    
    int count;
    
    if(LIST->next != NULL) {
        
        if(aux==1) {
            dem_collision_t *tmp = LIST->next;
            LIST->next = tmp->next;
        }
        
        else{
            dem_collision_t *current = LIST->next;
            dem_collision_t *previous = LIST;
            
            for(count=1; count<aux; count++) {
                previous = current;
                current = current->next;
            }
            
            previous->next = current->next;
            
        }
    }
}


void delete_all_contacts(dem_collision_t *LIST) {
    
    if(LIST->next != NULL) {
        dem_collision_t *current, *nextContact;
        
        current = LIST->next;
        
        while(current != NULL) {
            nextContact = current->next;
            current = nextContact;
        }
    }
}


int search_contact(dem_collision_t *LIST, int id_part) {
    
    int aux = 1;
    dem_collision_t *tmp = LIST->next;
    
    while(tmp != NULL) {
        
        if(tmp->idc == id_part) {
            return aux;
            break;
        }
        
        tmp = tmp->next;
        aux = aux+1;
    }
    
    return 0;
}


void dem_bounding_box(dem_particle_t *root) {
    
    dem_particle_t *ptr_dem = root;

        ptr_dem->xmax = ptr_dem->x + ptr_dem->diameter/2.0;
        ptr_dem->xmin = ptr_dem->x - ptr_dem->diameter/2.0;
        ptr_dem->ymax = ptr_dem->y + ptr_dem->diameter/2.0;
        ptr_dem->ymin = ptr_dem->y - ptr_dem->diameter/2.0;
        ptr_dem->zmax = ptr_dem->z + ptr_dem->diameter/2.0;
        ptr_dem->zmin = ptr_dem->z - ptr_dem->diameter/2.0;

}


void dem_overlap_save(dem_collision_t *LIST, int type, int aux, double overlap) {
    
    int count;
    
    if(aux==1) {
        if(type==1) LIST->dem_overlapn = overlap;
        else LIST->dem_overlapt = overlap;
    }
    
    else{
        dem_collision_t *current = LIST->next;
        dem_collision_t *previous = LIST;
        
        for(count=1; count<aux; count++) {
            previous = current;
            current = current->next;
        }
        
        if(type==1) current->dem_overlapn = overlap;
        else current->dem_overlapt = overlap;
    }
}


double dem_overlap_value(dem_collision_t *LIST, int type, int aux) {
    
    int count;
    double overlap;
    
    if(aux==1) {
        if(type==1) overlap = LIST->dem_overlapn;
        else overlap = LIST->dem_overlapt;
    }
    
    else{
        dem_collision_t *current = LIST->next;
        dem_collision_t *previous = LIST;
        
        for(count=1; count<aux; count++) {
            previous = current;
            current = current->next;
        }
        
        if(type==1) overlap = current->dem_overlapn;
        else overlap = current->dem_overlapt;
    }
    
    return overlap;
}


int dem_testing_contacts(dem_particle_t *root,int coll) {
    
    dem_particle_t *dem_particles = root;
    dem_particle_t *ptr_dem;
    
    int i, j, aux1, id_i, id_j, contact;
    double xmax_i, xmin_i, ymax_i, ymin_i, zmax_i, zmin_i;
    double xmax_j, xmin_j, ymax_j, ymin_j, zmax_j, zmin_j;
    
    
    for (i=1;i<=n_particles;i++) {
        aux1 = 1;
        
        id_i = i-1;
        
        ptr_dem = search_particle(dem_particles,id_i);
        
        xmin_i = ptr_dem->xmin;
        xmax_i = ptr_dem->xmax;
        ymin_i = ptr_dem->ymin;
        ymax_i = ptr_dem->ymax;
        zmin_i = ptr_dem->zmin;
        zmax_i = ptr_dem->zmax;
        
        
        for (j=1;j<=n_particles;j++) {
            
            id_j = j-1;
            
            ptr_dem = search_particle(dem_particles,id_j);
            
            xmin_j = ptr_dem->xmin;
            xmax_j = ptr_dem->xmax;
            ymin_j = ptr_dem->ymin;
            ymax_j = ptr_dem->ymax;
            zmin_j = ptr_dem->zmin;
            zmax_j = ptr_dem->zmax;
            
            if (id_i==id_j) continue;
            
            else if (((xmin_j>=xmin_i && xmin_j<=xmax_i) ||
                (xmax_j>=xmin_i && xmax_j<=xmax_i)) && ((ymin_j>=ymin_i && ymin_j<=ymax_i) || 
                (ymax_j>=ymin_i && ymax_j<=ymax_i)) && ((zmin_j>=zmin_i && zmin_j<=zmax_i) || 
                (zmax_j>=zmin_i && zmax_j<=zmax_i))) {
                
                ptr_dem = search_particle(dem_particles,id_i);
                contact = search_contact(ptr_dem->LIST,id_j);
                
                if(contact == 0) {
                    add_contact(ptr_dem->LIST,id_j);
                    ptr_dem->contact_dim = ptr_dem->contact_dim+1;
                    
                }
                
                aux1 = aux1+1;
                coll = coll+1;
                
                continue;
            }
            
        }
    }
    
    if(coll <= 0) {
        ptr_dem = search_particle(dem_particles,id_i);
        delete_all_contacts(ptr_dem->LIST);
    }
    
    return(coll);
}


int dem_collision(dem_particle_t *root, int coll) {
    
    dem_particle_t *dem_particles = root;
    dem_particle_t *ptr_dem;
    
    int i, j, aux1, contact, id_i, id_j;
    double xi, yi, zi, di, xj, yj, zj, dj;
    double posr, overlap;
    
    for (i=1;i<=n_particles;i++) {
        
        id_i = i-1;
        
        ptr_dem = search_particle(dem_particles,id_i);
        
        xi = ptr_dem->x;
        yi = ptr_dem->y;
        zi = ptr_dem->z;
        di = ptr_dem->diameter;
        
        for (j=1; j<=n_particles;j++) {
            aux1 = 1;
            
            id_j = j-1;
            
            ptr_dem = search_particle(dem_particles,id_j);
            
            xj = ptr_dem->x;
            yj = ptr_dem->y;
            zj = ptr_dem->z;
            dj = ptr_dem->diameter;
            
            if (id_i == id_j) continue;
            
            while (aux1>0) {
                
                ptr_dem = search_particle(dem_particles,id_i);
                contact = search_contact(ptr_dem->LIST,id_j);
                
                if(contact == 0) {
                    aux1 = -1;
                    continue;
                }
                
                else if (contact != 0) {
                    posr = sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj));
                    overlap = 0.5*(di+dj)-posr;
                    
                    if (overlap>0) {
                        dem_overlap_save(ptr_dem->LIST,1,contact,overlap);
                    }
                    
                    else {
                        delete_contact(ptr_dem->LIST,contact);
                        ptr_dem->contact_dim = ptr_dem->contact_dim-1;
                        coll = coll-1;
                    }
                    
                    aux1 = -1;
                }
                
                aux1 = aux1+1;
            }
        }
    }

    return (coll);
}


int sign(double num) {
    return (num > 0) - (num < 0);
}


void dem_contact_results(dem_particle_t *root) {
    
    dem_particle_t *dem_particles = root;
    dem_particle_t *ptr_dem;
    
    int aux1, contact;

    int i, j, id_i, id_j;
    double xi, yi, zi, di, ui, vi, wi, omega_xi, omega_yi, omega_zi, mi, Ipi, Ei, poisson_i, restc, fric;
    double xj, yj, zj, dj, uj, vj, wj, omega_xj, omega_yj, omega_zj, mj, Ipj;
    double Fcx, Fcy, Fcz, Mcx, Mcy, Mcz, overlapn, overlapt, posr;
    double nx, ny, nz, tx, ty, tz;
    double vrx, vry, vrz, vrn, vrnx, vrny, vrnz, vrt, vrtx, vrty, vrtz;
    double meff, kn, dampcn, beta, kt, dampct;
    double feln_x, feln_y, feln_z, felt_x, felt_y, felt_z;
    double fdissn_x, fdissn_y, fdissn_z, fdisst_x, fdisst_y, fdisst_z;
    double fn, fnx, fny, fnz, ft, ftx, fty, ftz;
    double fcx, fcy, fcz, Ttx, Tty, Ttz;
    
    
    for (i=1;i<=n_particles;i++) {
        
        id_i = i-1;
        
        ptr_dem = search_particle(dem_particles,id_i);
        
        xi = ptr_dem->x;
        yi = ptr_dem->y;
        zi = ptr_dem->z;
        di = ptr_dem->diameter;
        ui = ptr_dem->u;
        vi = ptr_dem->v;
        wi = ptr_dem->w;
        omega_xi = ptr_dem->omega_x;
        omega_yi = ptr_dem->omega_y;
        omega_zi = ptr_dem->omega_z;
        mi = ptr_dem->mp;
        Ipi = ptr_dem->Ip;
        Ei = ptr_dem->Ep;
        poisson_i = ptr_dem->poisson;
        restc = ptr_dem->coef_rest;
        fric = ptr_dem->coef_fric;
        
        Fcx = 0;
        Fcy = 0;
        Fcz = 0;
        
        Mcx = 0;
        Mcy = 0;
        Mcz = 0;
        
        for (j=1;j<=n_particles;j++) {
            
            aux1 = 1;
            
            id_j = j-1;
            
            ptr_dem = search_particle(dem_particles,id_j);
            
            xj = ptr_dem->x;
            yj = ptr_dem->y;
            zj = ptr_dem->z;
            dj = ptr_dem->diameter;
            uj = ptr_dem->u;
            vj = ptr_dem->v;
            zj = ptr_dem->w;
            omega_xj = ptr_dem->omega_x;
            omega_yj = ptr_dem->omega_y;
            omega_zj = ptr_dem->omega_z;
            mj = ptr_dem->mp;
            Ipj = ptr_dem->Ip;
            
            if (id_i==id_j) continue;
            
            while (aux1>0) {
                
                ptr_dem = search_particle(dem_particles,id_i);
                contact = search_contact(ptr_dem->LIST, id_j);
                
                if(contact == 0) {
                    aux1 = -1;
                    continue;
                }
                
                else if (contact != 0) {
                    
                    // Normal forces and velocities
                    overlapn = dem_overlap_value(ptr_dem->LIST,1,contact);
                    
                    posr =  sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj));
                    nx = (xj-xi)/posr;
                    ny = (yj-yi)/posr;
                    nz = (zj-zi)/posr;
                    
                    vrx = (ui-uj) + (0.5*((di*omega_yi+dj*omega_yj)*nz - (di*omega_zi+dj*omega_zj)*ny));
                    vry = (vi-vj) + (0.5*((di*omega_zi+dj*omega_zj)*nx - (di*omega_xi+dj*omega_xj)*nz));
                    vrz = (wi-wj) + (0.5*((di*omega_xi+dj*omega_xj)*ny - (di*omega_yi+dj*omega_yj)*nx));
                    
                    vrn = vrx*nx+vry*ny+vrz*nz;
                    vrnx = vrn*nx;
                    vrny = vrn*ny;
                    vrnz = vrn*nz;
                    
                    meff = (mi*mj)/(mi+mj);
                    
                    kn = Ei*di; // Non-cohesive dry friction model
                    
                    dampcn = -2*log(restc)*sqrt(meff*kn)/sqrt((log(restc)*log(restc))+M_PI*M_PI);
                    beta = dampcn/kn;
                    
                    feln_x = -kn*overlapn*nx;
                    feln_y = -kn*overlapn*ny;
                    feln_z = -kn*overlapn*nz;
                    
                    fdissn_x = -dampcn*vrnx;
                    fdissn_y = -dampcn*vrny;
                    fdissn_z = -dampcn*vrnz;
                    
                    fnx = feln_x+fdissn_x;
                    fny = feln_y+fdissn_y;
                    fnz = feln_z+fdissn_z;
                    fn = sqrt(fnx*fnx+fny*fny+fnz*fnz);
                    
                    // Tangencial forces and velocities
                    vrtx = vrx-vrnx;
                    vrty = vry-vrny;
                    vrtz = vrz-vrnz;
                    vrt = sqrt(vrtx*vrtx+vrty*vrty+vrtz*vrtz);
                    
                    overlapt = dem_overlap_value(ptr_dem->LIST, 2, contact);
                    
                    if (vrt == 0) {
                        tx = 0;
                        ty = 0;
                        tz = 0;
                    }
                    
                    else {
                        tx = vrtx/vrt;
                        ty = vrty/vrt;
                        tz = vrtz/vrt;
                    }
                    
                    kt = kn/poisson_i;
                    dampct = beta*kt;
                    
                    // Need to check
                    
                    felt_x = -kt*overlapt*tx;
                    felt_y = -kt*overlapt*ty;
                    felt_z = -kt*overlapt*tz;
                    
                    fdisst_x = -dampct*vrtx;
                    fdisst_y = -dampct*vrty;
                    fdisst_z = -dampct*vrtz;
                    
                    ftx = felt_x+fdisst_x;
                    fty = felt_y+fdisst_y;
                    ftz = felt_z+fdisst_z;
                    ft = sqrt(ftx*ftx+fty*fty+ftz*ftz);
                    
                    if (ft >= fric*fn) {
                        
                        ftx = -fric*fn*sign(overlapt)*tx;
                        fty = -fric*fn*sign(overlapt)*ty;
                        ftz = -fric*fn*sign(overlapt)*tz;
                    }
                    
                    fcx = fnx+ftx;
                    fcy = fny+fty;
                    fcz = fnz+ftz;

                    Fcx = Fcx+fcx;
                    Fcy = Fcy+fcy;
                    Fcz = Fcz+fcz;
                    
                    // Torque
                    Ttx = 0.5*di*(ny*fcz - nz*fcy);
                    Tty = 0.5*di*(nz*fcx - nx*fcz);
                    Ttz = 0.5*di*(nx*fcy - ny*fcx);
                    
                    Mcx = Mcx+Ttx;
                    Mcy = Mcy+Tty;
                    Mcz = Mcz+Ttz;
                    
                    overlapt = overlapt+vrt*dt;
                    dem_overlap_save(ptr_dem->LIST,2,contact,overlapt);
                    
                    aux1 = -1;
                } 
                
                aux1 = aux1+1;
            }
        }

        ptr_dem = search_particle(dem_particles,id_i);
        ptr_dem->Fcontact_x = Fcx;
        ptr_dem->Fcontact_y = Fcy;
        ptr_dem->Fcontact_z = Fcz;
        ptr_dem->Mcontact_x = Mcx;
        ptr_dem->Mcontact_y = Mcy;
        ptr_dem->Mcontact_z = Mcz;
    }
    
}


void dem_adv(dem_particle_t *root, int coll) {
    
    dem_particle_t *ptr_dem = root;
    
    double Fcontact_x = 0;
    double Fcontact_y = 0;
    double Fcontact_z = 0;
    
    double Mcontact_x = 0;
    double Mcontact_y = 0;
    double Mcontact_z = 0;
    
    ptr_dem->mp = ptr_dem->rho*((1./6.)*M_PI*pow(ptr_dem->diameter,3));
    ptr_dem->Ip = (1./10.)*ptr_dem->mp*pow(ptr_dem->diameter,2);
    
    if (coll>0) {
        
        Fcontact_x = ptr_dem->Fcontact_x;
        Fcontact_y = ptr_dem->Fcontact_y;
        Fcontact_z = ptr_dem->Fcontact_z;
        
        Mcontact_x = ptr_dem->Mcontact_x;
        Mcontact_y = ptr_dem->Mcontact_y;
        Mcontact_z = ptr_dem->Mcontact_z;
        
    }
    
    ptr_dem->u = ptr_dem->u + (gx + Fcontact_x/ptr_dem->mp)*dt;
    ptr_dem->v = ptr_dem->v + (gy + Fcontact_y/ptr_dem->mp)*dt;
    ptr_dem->w = ptr_dem->w + (gz + Fcontact_z/ptr_dem->mp)*dt;
    
    ptr_dem->x = ptr_dem->x + ptr_dem->u*dt;
    ptr_dem->y = ptr_dem->y + ptr_dem->v*dt;
    ptr_dem->z = ptr_dem->z + ptr_dem->w*dt;
    
    ptr_dem->omega_x = ptr_dem->omega_x + (Mcontact_x/ptr_dem->Ip)*dt;
    ptr_dem->omega_y = ptr_dem->omega_y + (Mcontact_y/ptr_dem->Ip)*dt;
    ptr_dem->omega_x = ptr_dem->omega_z + (Mcontact_z/ptr_dem->Ip)*dt;

    
    // Walls
    if (ptr_dem->xmin < xd1) {
        ptr_dem->x = xd1+ptr_dem->diameter/2.0;
        ptr_dem->u = -ptr_dem->u;
    }
    
    else if (ptr_dem->xmax > xd2) {
        ptr_dem->x = xd2-ptr_dem->diameter/2.0;
        ptr_dem->u = -ptr_dem->u;
    }
    
    if (ptr_dem->ymin < yd1) {
        ptr_dem->y = yd1+ptr_dem->diameter/2.0;
        ptr_dem->v = -ptr_dem->v;
    }
    
    else if (ptr_dem->ymax > yd2) {
        ptr_dem->y = yd2-ptr_dem->diameter/2.0;
        ptr_dem->v = -ptr_dem->v;
    }
    
    if (ptr_dem->zmin < zd1) {
        ptr_dem->z = zd1+ptr_dem->diameter/2.0;
        ptr_dem->w = -ptr_dem->w;
    }
    
    else if (ptr_dem->zmax > zd2) {
        ptr_dem->z = zd2-ptr_dem->diameter/2.0;
        ptr_dem->w = -ptr_dem->w;
    }
    
}


// Print vtk file to run using paraview
void dem_print_vtk(dem_particle_t *root) {
    
    dem_particle_t *dem_particles = root;
    dem_particle_t *ptr_dem;
    
    char string[500];
    int id;
    
    sprintf(string, "%s/dem_output_%i.vtk",output,iter);
    
    FILE *pvtk;
    pvtk = fopen(string, "w");
    
    fprintf(pvtk,"# vtk DataFile Version 4.1\nPoints\nASCII\nDATASET POLYDATA\n");
    fprintf(pvtk,"POINTS	%i	double\n",n_particles);
    
    for (id=0;id<n_particles;id++) {
        ptr_dem = search_particle(dem_particles,id);
        fprintf(pvtk,"%6.16lf	%6.16lf	%6.16lf\n",ptr_dem->x, ptr_dem->y, ptr_dem->z);
    }
    
    fprintf(pvtk,"POINT_DATA	%i\n", n_particles);
    
    fprintf(pvtk,"SCALARS	id float\n");
    fprintf(pvtk,"LOOKUP_TABLE default\n");
    
    for (id=0;id<n_particles;id++) {
        ptr_dem = search_particle(dem_particles,id);
        fprintf(pvtk,"%i\n",ptr_dem->idp);
    }
    
    fprintf(pvtk,"SCALARS	diameter float\n");
    fprintf(pvtk,"LOOKUP_TABLE default\n");
    
    for (id=0;id<n_particles;id++) {
        ptr_dem = search_particle(dem_particles,id);
        fprintf(pvtk,"%6.16f\n",ptr_dem->diameter);
    }
    
    fprintf(pvtk,"SCALARS	u float\n");
    fprintf(pvtk,"LOOKUP_TABLE default\n");
    
    for (id=0;id<n_particles;id++) {
        ptr_dem = search_particle(dem_particles,id);
        fprintf(pvtk,"%6.16f\n",ptr_dem->u);
    }
    
    fprintf(pvtk,"SCALARS	v float\n");
    fprintf(pvtk,"LOOKUP_TABLE default\n");
    
    for (id=0;id<n_particles;id++) {
        ptr_dem = search_particle(dem_particles,id);
        fprintf(pvtk,"%6.16f\n",ptr_dem->v);
    }
    
    fprintf(pvtk,"SCALARS	w float\n");
    fprintf(pvtk,"LOOKUP_TABLE default\n");
    
    for (id=0;id<n_particles;id++) {
        ptr_dem = search_particle(dem_particles,id);
        fprintf(pvtk,"%6.16f\n",ptr_dem->w);
    }
    
    fprintf(pvtk,"SCALARS	omega_x float\n");
    fprintf(pvtk,"LOOKUP_TABLE default\n");
    
    for (id=0;id<n_particles;id++) {
        ptr_dem = search_particle(dem_particles,id);
        fprintf(pvtk,"%6.16f\n",ptr_dem->omega_x);
    }
    
    fprintf(pvtk,"SCALARS	omega_y float\n");
    fprintf(pvtk,"LOOKUP_TABLE default\n");
    
    for (id=0;id<n_particles;id++) {
        ptr_dem = search_particle(dem_particles,id);
        fprintf(pvtk,"%6.16f\n",ptr_dem->omega_y);
    }
    
    fprintf(pvtk,"SCALARS	omega_z float\n");
    fprintf(pvtk,"LOOKUP_TABLE default\n");
    
    for (id=0;id<n_particles;id++) {
        ptr_dem = search_particle(dem_particles,id);
        fprintf(pvtk,"%6.16f\n",ptr_dem->omega_z);
    }
    
    fprintf(pvtk,"SCALARS	Fcontact_x float\n");
    fprintf(pvtk,"LOOKUP_TABLE default\n");
    
    for (id=0;id<n_particles;id++) {
        ptr_dem = search_particle(dem_particles,id);
        fprintf(pvtk,"%6.16f\n",ptr_dem->Fcontact_x);
    }
    
    fprintf(pvtk,"SCALARS	Fcontact_y float\n");
    fprintf(pvtk,"LOOKUP_TABLE default\n");
    
    for (id=0;id<n_particles;id++) {
        ptr_dem = search_particle(dem_particles,id);
        fprintf(pvtk,"%6.16f\n",ptr_dem->Fcontact_y);
    }
    
    fprintf(pvtk,"SCALARS	Fcontact_z float\n");
    fprintf(pvtk,"LOOKUP_TABLE default\n");
    
    for (id=0;id<n_particles;id++) {
        ptr_dem = search_particle(dem_particles,id);
        fprintf(pvtk,"%6.16f\n",ptr_dem->Fcontact_z);
    }
    
    fprintf(pvtk,"SCALARS	Mcontact_x float\n");
    fprintf(pvtk,"LOOKUP_TABLE default\n");
    
    for (id=0;id<n_particles;id++) {
        ptr_dem = search_particle(dem_particles,id);
        fprintf(pvtk,"%6.16f\n",ptr_dem->Mcontact_x);
    }
    
    fprintf(pvtk,"SCALARS	Mcontact_y float\n");
    fprintf(pvtk,"LOOKUP_TABLE default\n");
    
    for (id=0;id<n_particles;id++) {
        ptr_dem = search_particle(dem_particles,id);
        fprintf(pvtk,"%6.16f\n",ptr_dem->Mcontact_y);
    }
    
    fprintf(pvtk,"SCALARS	Mcontact_z float\n");
    fprintf(pvtk,"LOOKUP_TABLE default\n");
    
    for (id=0;id<n_particles;id++) {
        ptr_dem = search_particle(dem_particles,id);
        fprintf(pvtk,"%6.16f\n",ptr_dem->Mcontact_z);
    }
    
    fclose(pvtk);
}