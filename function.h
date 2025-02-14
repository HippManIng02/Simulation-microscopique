#pragma once 
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib> 
#include <ctime> 
#include <iomanip>

#define N_particules_total 1000
#define r_etoile 3.0
#define my_epsilon 0.2 
#define r_etoile_2  r_etoile * r_etoile // r étoile au carré
#define GAMMA 0.01 //Facteur gamma du thermostat Berendsen
#define dt 1.0 // Pas de temps en femtosecondes (fs)
#define CONVERSION_FORCE 0.0001 * 4.186  // Facteur de conversion
#define CONSTANTE_R  0.00199 // Constante R
#define MASSE_PARTICULE 18.0 // Masse des particules
#define T0 300.0 // Température initiale
#define L 42.0 // Dimension de la boîte

#define r_etoile_6  r_e_carree * r_e_carree * r_e_carree // (r^*)^6
#define r_etoile_12 r_etoile_6 * r_etoile_6 // (r^*)^12


typedef double f64;


struct Particules{
 std::vector<f64> x, y, z; //Tableaux contenant l'ensemble des coordonnées des particules x, y et z 
 std::vector<f64> fx, fy,fz; //Tableaux contenant l'ensemble des forces des particules dans chaque composant x, y et z
 std::vector<f64> Mx, My,Mz; // Tableaux contenant l'ensemble des moments cinétiques des particules dans chaque composant x, y et z
};

// Définition des valeurs du vecteur de translation de façon statique
struct vecteurs_translation {
    static constexpr int N = 27;

    static constexpr double x[N] = {
          0.0, -L, -L, -L,   0.0,   0.0,   0.0,  L,  L,
        -L, -L, -L,   0.0,  L,  L, -L,   0.0,  L,
        -L, -L, -L,   0.0,   0.0,   0.0,  L,  L,  L
    };

    static constexpr double y[N] = {
          0.0, -L,   0.0,  L, -L,   0.0,  L, -L,   0.0,
        -L,   0.0,  L, -L, -L,   0.0,   0.0,   0.0,   0.0,
         L,   0.0, -L, -L,   0.0,  L, -L,   0.0,  L
    };

    static constexpr double z[N] = {
          0.0, -L, -L, -L,   0.0,   0.0,   0.0, -L, -L,
        -L,   0.0,  L,  L, -L, -L,   0.0,   0.0,   0.0,
         L,  L,  L,   0.0,   0.0,   0.0,  L,  L,  L
    };
};




unsigned int get_nombre_particules(const std::string& nom_fichier);

//calcul de distance au carrée
f64 calcul_distance_carree(f64 dx, f64 dy, f64 dz);


//Calcul de l'énergie et de la force en utilisant le potentiel de Lennard Jones
void calcul_energie_force_LJ(Particules& p, uint taille, bool mode_periodique, f64 *U,f64 r_cut, uint N_sym = 27,f64 r_e_carree = r_etoile_2, f64 n_epsilon = my_epsilon);

// Algorithme de velocity-verlet
void velocity_verlet(Particules &p, uint taille,f64 r_cut, f64 *U,bool mode_periodique, uint N_sym = 27,f64 r_e_carree = r_etoile_2, f64 n_epsilon = my_epsilon);

//Calcul de l'énergie cinétique
void calcul_energie_cinetique_temperature(const Particules &p,uint taille , f64 *EC, f64 *TC);

//Génération des moments cinétiques et recabrage pour T0 = 300 K
void calcul_moments_cinetiques_init(Particules &p, uint taille);

//calcul de l'énergie cinétique
void calcul_energie_cinetique(const Particules &p, uint &taille, f64 *EC);

//correction rapport
void correction_rapport(Particules &p, uint &taille);

//Correction du moment cinétique par rapport au centre de masse
void correction_moments_cinetiques(Particules &p, uint taille);

//Fonction signe
f64 fonction_signe(f64 valeur, f64 s);

//Correction moment cinétique à l'aide du thermostat de Berendsen
void correction_moment_cinetique_thermostat_berendsen(Particules &p, uint taille, f64 t, f64 gamma = GAMMA, f64 t0 = T0);

//Fonction de sauvegarde de donnée dans un fichier PDB
void sauvegarder_trajectoire_PDB(const Particules &p, uint taille, const std::string& nom_fichier, uint iteration, uint XDIM, uint YDIM, uint ZDIM);

//Initialisation des tableau à 0
void initialisation(Particules& p, uint taille);
//Initialisation des forces à 0
void init_force(Particules& p, uint taille);

//Verification de la nullitée des forces 
void verifier_valeur_force(const std::vector<f64>& fx, const std::vector<f64>& fy, const std::vector<f64>& fz, uint taille );

