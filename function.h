#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <iomanip>

// Définition du type flottant 64 bits
using f64 = double;

// Nombre total de particules dans la simulation
constexpr uint N_PARTICULES_TOTAL = 1000;

// Déclarations avec extern
extern f64 R_ETOILE;
extern f64 MY_EPSILON;
extern f64 GAMMA;
extern f64 DT;
extern f64 CONVERSION_FORCE;
extern f64 CONSTANTE_R;
extern f64 MASSE_PARTICULE;
extern f64 T0;
extern f64 L;
extern uint N_particules_local;
extern uint N_sym;
extern uint n_iterations;
extern uint m_step;
extern f64 r_cut;
extern bool mode_periodique;
extern std::string nom_fichier;

// Déclarations des puissances de R_ETOILE
extern f64 R_ETOILE_2;
extern f64 R_ETOILE_6;
extern f64 R_ETOILE_7;
extern f64 R_ETOILE_8;
extern f64 R_ETOILE_12;
extern f64 R_ETOILE_13;
extern f64 R_ETOILE_14;

// Autres constantes spécifiques
extern f64 FACTEUR_CONVERSION;
extern f64 FACTEUR_TEMPERATURE;

// // Paramètres du potentiel de Lennard-Jones
// f64 R_ETOILE = 3.0;   // Distance caractéristique (r*)
// f64 MY_EPSILON = 0.2; // Profondeur du puits de potentiel (ε)

// // Facteurs liés à la simulation
// f64 GAMMA = 0.01; // Facteur gamma du thermostat Berendsen
// f64 DT = 1.0;     // Pas de temps en femtosecondes (fs)

// // Facteurs de conversion
// f64 CONVERSION_FORCE = 0.0001 * 4.186; // Facteur de conversion des unités de force
// f64 CONSTANTE_R = 0.00199;             // Constante R en kcal/(mol·K)

// // Caractéristiques des particules
// f64 MASSE_PARTICULE = 18.0; // Masse d'une particule en unités de simulation
// f64 T0 = 300.0;             // Température initiale en Kelvin
// f64 L = 42.0;               // Longueur de la boîte de simulation (cube de côté L)

// uint N_particules_local = 1000;
// uint N_sym = 27;
// uint n_iterations = 1000;
// uint m_step = 10;
// bool mode_periodique = false;
// f64 r_cut = 10.0;

// // Calculs des puissances de r*
// f64 R_ETOILE_2 = R_ETOILE * R_ETOILE;                  // r_etoile^2
// f64 R_ETOILE_6 = R_ETOILE_2 * R_ETOILE_2 * R_ETOILE_2; // r_etoile^6
// f64 R_ETOILE_7 = R_ETOILE_6 * R_ETOILE;                // r_etoile^7
// f64 R_ETOILE_8 = R_ETOILE_6 * R_ETOILE_2;              // r_etoile^8
// f64 R_ETOILE_12 = R_ETOILE_6 * R_ETOILE_6;             // r_etoile^12
// f64 R_ETOILE_13 = R_ETOILE_12 * R_ETOILE;              // r_etoile^13
// f64 R_ETOILE_14 = R_ETOILE_12 * R_ETOILE_2;            // r_etoile^14

// // Constantes pour les calculs
// f64 FACTEUR_TEMPERATURE = 1.0 / CONSTANTE_R;             // Facteur pour le calcul de la température
// f64 FACTEUR_CONVERSION = 1.0 / (CONVERSION_FORCE * 2.0); // Facteur pour normaliser l'énergie cinétique

constexpr uint max_voisins = 100; // Définir le nombre max de voisins estimé

// Structure contenant les informations des particules (positions, forces et moments)
struct Particules
{
  std::vector<f64> x, y, z;    // Coordonnées des particules
  std::vector<f64> fx, fy, fz; // Forces appliquées aux particules
  std::vector<f64> Mx, My, Mz; // Moments cinétiques des particules
};

// Définition des vecteurs de translation utilisés pour les conditions périodiques
struct vecteurs_translation
{
  static constexpr uint N = 27;
  static f64 x[N];
  static f64 y[N];
  static f64 z[N];
};

// Récupère le nombre de particules à partir d'un fichier
uint get_nombre_particules(const std::string &nom_fichier);

// Calcule la distance au carré entre deux points
f64 calcul_distance_carree(f64 dx, f64 dy, f64 dz);

// Calcule l'énergie et la force avec le potentiel de Lennard-Jones
void calcul_energie_force_LJ(Particules &p, uint taille, bool mode_periodique, f64 *U, f64 r_cut, uint N_sym = 27, f64 n_epsilon = MY_EPSILON);

// Vérifie si les distances sont bien corrigées dans les conditions périodiques
void verifier_distances_corrigees(f64 dx, f64 dy, f64 dz);

// Algorithme de velocity-verlet pour intégrer les équations du mouvement
void velocity_verlet(Particules &p, uint taille, f64 r_cut, f64 *U, bool mode_periodique, uint N_sym = 27);

// Applique les conditions périodiques aux particules
void appliquer_conditions_periodiques(Particules &p, uint &taille, f64 Lx, f64 Ly, f64 Lz);

// Calcule l'énergie cinétique et la température du système
void calcul_energie_cinetique_temperature(const Particules &p, uint taille, f64 *EC, f64 *TC);

// Initialise les moments cinétiques des particules
void calcul_moments_cinetiques_init(Particules &p, uint taille);

// Calcule l'énergie cinétique du système
void calcul_energie_cinetique(const Particules &p, uint &taille, f64 *EC);

// Correction du rapport d'énergie cinétique
void correction_rapport(Particules &p, uint &taille);

// Corrige les moments cinétiques par rapport au centre de masse
void correction_moments_cinetiques(Particules &p, uint taille);

// Fonction mathématique pour déterminer le signe d'une valeur
f64 fonction_signe(f64 valeur, f64 s);

// Correction du moment cinétique avec le thermostat de Berendsen
void correction_moment_cinetique_thermostat_berendsen(Particules &p, uint taille, f64 t, f64 gamma = GAMMA, f64 t0 = T0);

// Sauvegarde des trajectoires des particules dans un fichier PDB
void sauvegarder_trajectoire_PDB(const Particules &p, uint taille, const std::string &nom_fichier, uint iteration, uint XDIM, uint YDIM, uint ZDIM);

// Initialise les moments cinétiques en fonction des vitesses des particules
void init_moment_avec_vitesse(Particules &p, const std::string &nom_fichier);

// Initialise tous les tableaux à zéro
void initialisation(Particules &p, uint taille);

// Initialise les forces des particules à zéro
void init_force(Particules &p, uint taille);

// Vérifie que les forces appliquées sur les particules sont nulles ou cohérentes
void verifier_valeur_force(const std::vector<f64> &fx, const std::vector<f64> &fy, const std::vector<f64> &fz, uint taille);

// Fonction de construction de la liste des voisins
void construire_liste_voisins(Particules &p, uint taille, f64 r_cut, uint N_sym = 27);

// Calcul des forces et de l’énergie potentielle avec la liste de Verlet
void calcul_energie_force_LJ_voisins(Particules &p, uint taille, f64 *U, f64 n_epsilon, uint N_sym = 27);

/*Fonction utilitaires*/

// Fonction permettant parser les arguments
void parse_arguments(int argc, char* argv[]);
// Fonction permettant d'afficher l'aide
void afficher_aide();
