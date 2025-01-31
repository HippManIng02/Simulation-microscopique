#pragma once 

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#define N_particules_total 1000
#define r_etoile 3.0
#define my_epsilon 0.2
#define r_etoile_2 r_etoile * r_etoile

typedef double f64;

struct Particules{
 std::vector<f64> x, y, z; //Tableaux contenant l'ensemble des coordonnées des particules x, y et z 
 std::vector<f64> fx, fy,fz; //Tableaux contenant l'ensemble des forces des particules dans chaque composant fx, fy et fz
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
void calcul_energie_force_LJ(Particules& p, unsigned int taille, bool mode_periodique, f64 *U,f64 r_cut, unsigned int N_sym = 27,f64 r_e_carree = r_etoile_2, f64 n_epsilon = my_epsilon);

//Initialisation des forces à 0
void init_force(Particules& p, unsigned int taille);

//Verification de la nullitée des forces 
void verifier_valeur_force(const std::vector<f64>& fx, const std::vector<f64>& fy, const std::vector<f64>& fz, unsigned int taille );

