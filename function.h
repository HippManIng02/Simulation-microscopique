#pragma once 
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib> 
#include <ctime> 

#define N_particules_total 1000
#define r_etoile 3.0
#define my_epsilon 0.2
#define r_etoile_2 r_etoile * r_etoile
#define GAMMA 0.01


typedef double f64;

const f64 dt = 1.0; // Pas de temps en femtosecondes (fs)
const f64 CONVERSION_FORCE = 0.0001 * 4.186; // Facteur de conversion
const f64 CONSTANTE_R = 0.00199; // Constante R
const f64 MASSE_PARTICULE = 18.0; // Masse des particules
const f64 T0 = 300.0; // Température initiale

struct Particules{
 std::vector<f64> x, y, z; //Tableaux contenant l'ensemble des coordonnées des particules x, y et z 
 std::vector<f64> fx, fy,fz; //Tableaux contenant l'ensemble des forces des particules dans chaque composant x, y et z
 std::vector<f64> vx, vy,vz; // Tableaux contenant l'ensemble des vitesses des particules dans chaque composant x, y et z
 std::vector<f64> Mx, My,Mz; // Tableaux contenant l'ensemble des moments cinétiques des particules dans chaque composant x, y et z
};

// Définition des valeurs du vecteur de translation de façon statique
struct vecteurs_translation {
    static constexpr int N = 27;

    static constexpr double x[N] = {
          0.0, -42.0, -42.0, -42.0,   0.0,   0.0,   0.0,  42.0,  42.0,
        -42.0, -42.0, -42.0,   0.0,  42.0,  42.0, -42.0,   0.0,  42.0,
        -42.0, -42.0, -42.0,   0.0,   0.0,   0.0,  42.0,  42.0,  42.0
    };

    static constexpr double y[N] = {
          0.0, -42.0,   0.0,  42.0, -42.0,   0.0,  42.0, -42.0,   0.0,
        -42.0,   0.0,  42.0, -42.0, -42.0,   0.0,   0.0,   0.0,   0.0,
         42.0,   0.0, -42.0, -42.0,   0.0,  42.0, -42.0,   0.0,  42.0
    };

    static constexpr double z[N] = {
          0.0, -42.0, -42.0, -42.0,   0.0,   0.0,   0.0, -42.0, -42.0,
        -42.0,   0.0,  42.0,  42.0, -42.0, -42.0,   0.0,   0.0,   0.0,
         42.0,  42.0,  42.0,   0.0,   0.0,   0.0,  42.0,  42.0,  42.0
    };
};




unsigned int get_nombre_particules(const std::string& nom_fichier);

//calcul de distance au carrée
f64 calcul_distance_carree(f64 x_1, f64 y_1, f64 z_1, f64 x_2, f64 y_2, f64 z_2);

//calcul du terme (r*)²/(r)²
f64 calcul_de_r_etoile_div_par_r__2(f64 r_carree, f64 r_e_carree = r_etoile_2);

//Calcul du potentiel de Lennard jones
f64 potentiel_lennard_jones(f64 r_carree, f64 r_e_carree = r_etoile_2, f64 n_epsilon = my_epsilon);


//Calcul de la force de Lennard jones, (dérivée du potentiel)
f64 force_lennard_jones(f64 r_carree, f64 r_e_carree = r_etoile_2, f64 n_epsilon = my_epsilon);


//Calcul de l'énergie et de la force en utilisant le potentiel de Lennard Jones
void calcul_energie_force_LJ(Particules& p, uint taille, bool mode_periodique, f64 *U,f64 r_cut, uint N_sym = 27,f64 r_e_carree = r_etoile_2, f64 n_epsilon = my_epsilon);

//Calcul de l'énergie en utilisant le potentiel de Lennard Jones
void calcul_energie_LJ(Particules& p, uint taille, bool mode_periodique, f64 *U,f64 r_cut, uint N_sym = 27,f64 r_e_carree = r_etoile_2, f64 n_epsilon = my_epsilon);


// Algorithme de velocity-verlet
void velocity_verlet(Particules &p, uint taille,f64 r_cut, uint N_sym = 27,f64 r_e_carree = r_etoile_2, f64 n_epsilon = my_epsilon);

//Calcul de l'énergie cinétique
void calcul_energie_cinetique(const Particules &p,uint taille , f64 *EC, f64 *TC);

//Génération des moments cinétiques et recabrage pour T0 = 300 K
void calcul_moments_cinetiques(Particules &p, uint taille);

//Correction du moment cinétique par rapport au centre de masse
void correction_moments_cinetiques(Particules &p, uint taille);

//Fonction signe
f64 fonction_signe(f64 valeur, f64 s);

//Correction moment cinétique à l'aide du thermostat de Berendsen
void correction_moment_cinetique_thermostat_berendsen(Particules &p, uint taille, f64 t, f64 gamma = GAMMA, f64 t0 = T0);

//Initialisation des tableau à 0
void initialisation(Particules& p, uint taille);
//Initialisation des forces à 0
void init_force(Particules& p, uint taille);

//Verification de la nullitée des forces 
void verifier_valeur_force(const std::vector<f64>& fx, const std::vector<f64>& fy, const std::vector<f64>& fz, uint taille );

