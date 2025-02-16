#include "function.h"
#include <limits>
#include <cmath>
#include <unistd.h>

// Tableau global pour stocker les voisins
int liste_voisins[N_PARTICULES_TOTAL][2 * max_voisins];

// Définition des variables globales
f64 R_ETOILE = 3.0;
f64 MY_EPSILON = 0.2;
f64 GAMMA = 0.01;
f64 DT = 1.0;
f64 CONVERSION_FORCE = 0.0001 * 4.186;
f64 CONSTANTE_R = 0.00199;
f64 MASSE_PARTICULE = 18.0;
f64 T0 = 300.0;
f64 L = 42.0;
uint N_particules_local = 1000;
uint N_sym = 27;
uint n_iterations = 1000;
uint m_step = 10;
f64 r_cut = 10.0;
bool mode_periodique = false;
std::string nom_fichier = "particule.xyz";

// Définition des puissances de R_ETOILE
f64 R_ETOILE_2 = R_ETOILE * R_ETOILE;
f64 R_ETOILE_6 = R_ETOILE_2 * R_ETOILE_2 * R_ETOILE_2;
f64 R_ETOILE_7 = R_ETOILE_6 * R_ETOILE;
f64 R_ETOILE_8 = R_ETOILE_6 * R_ETOILE_2;
f64 R_ETOILE_12 = R_ETOILE_6 * R_ETOILE_6;
f64 R_ETOILE_13 = R_ETOILE_12 * R_ETOILE;
f64 R_ETOILE_14 = R_ETOILE_12 * R_ETOILE_2;

// Autres constantes spécifiques
f64 FACTEUR_TEMPERATURE = 1.0 / CONSTANTE_R;  // Facteur pour le calcul de la température
f64 FACTEUR_CONVERSION = 1.0 / (CONVERSION_FORCE * 2.0);

f64 vecteurs_translation::x[N] = {
    -L, -L, -L, 0.0, 0.0, 0.0, L, L, L, -L, -L, -L, 0.0, 0.0, 0.0, L, L, L, -L, -L, -L, 0.0, 0.0, 0.0, L, L, L
};

f64 vecteurs_translation::y[N] = {
    -L, 0.0, L, -L, 0.0, L, -L, 0.0, L, -L, 0.0, L, -L, 0.0, L, -L, 0.0, L, -L, 0.0, L, -L, 0.0, L, -L, 0.0, L
};

f64 vecteurs_translation::z[N] = {
    -L, -L, -L, -L, -L, -L, -L, -L, -L, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, L, L, L, L, L, L, L, L, L
};


// Fonction pour compter le nombre de particules
uint get_nombre_particules(const std::string &nom_fichier)
{
    std::ifstream fichier(nom_fichier);
    uint compteur = 0;
    if (!fichier.is_open())
    {
        std::cerr << "Une erreur est survenue lors de l'ouverture du fichier pour le comptage du nombre de particule" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Compter le nombre de ligne dans le fichier
    std::string tmp;
    while (std::getline(fichier, tmp))
        compteur++;

    // Retrait de la première ligne du fichier qui sert d'indication
    compteur > 0 && compteur--;

    fichier.close();
    return compteur;
}

// calcul de distance au carrée
f64 calcul_distance_carree(f64 dx, f64 dy, f64 dz)
{
    f64 r_carree = (dx * dx) + (dy * dy) + (dz * dz);
    if (r_carree == 0)
    {
        std::cerr << "Attention, la valeur de r² = 0" << std::endl;
        exit(EXIT_FAILURE);
    }
    return r_carree;
}

// Calcul de l'énergie et de la force en utilisant le potentiel de Lennard Jones
void calcul_energie_force_LJ(Particules &p, uint taille, bool mode_periodique, f64 *U, f64 r_cut, uint N_sym, f64 n_epsilon)
{
    init_force(p, taille);

    vecteurs_translation vecteurs_t;
    const uint n_sym_effectif = mode_periodique == 1 ? N_sym : 1;
    const f64 r_cut_effectif = mode_periodique == 1 ? r_cut * r_cut : std::numeric_limits<f64>::max();
    const f64 multiplicateur_potentiel = mode_periodique == 1 ? 2.0 : 4.0;
    const f64 facteur = mode_periodique == 1 ? -12.0 : -24.0;

    for (uint i_sym = 0; i_sym < n_sym_effectif; i_sym++)
    {
        for (uint i = 0; i < taille; i++)
        {
            f64 x = p.x[i];
            f64 y = p.y[i];
            f64 z = p.z[i];

            for (uint j = 0; j < taille; j++)
            {
                if (i == j)
                    continue;
                f64 x_j = p.x[j] + (mode_periodique == 1 ? vecteurs_t.x[i_sym] : 0.0);
                f64 y_j = p.y[j] + (mode_periodique == 1 ? vecteurs_t.y[i_sym] : 0.0);
                f64 z_j = p.z[j] + (mode_periodique == 1 ? vecteurs_t.z[i_sym] : 0.0);

                f64 dx = x - x_j;
                f64 dy = y - y_j;
                f64 dz = z - z_j;

                f64 r_carree = calcul_distance_carree(dx, dy, dz);

                if (r_carree <= r_cut_effectif)
                {
                    f64 inv_r2 = 1.0 / r_carree;
                    f64 inv_r6 = inv_r2 * inv_r2 * inv_r2;  // (1 / r^6)
                    f64 inv_r14 = inv_r6 * inv_r6 * inv_r2; // (1 / r^14)

                    f64 factor = R_ETOILE_6 * inv_r6;
                    // Mise à jour de l'énergie potentielle
                    *U += n_epsilon * (factor * factor - factor);

                    f64 factor1 = R_ETOILE_13 * inv_r14;
                    f64 factor2 = R_ETOILE_7 * inv_r6 * inv_r2;
                    // Calcul de la force
                    f64 force = facteur * n_epsilon * (2 * factor1 - factor2);

                    p.fx[i] += (force * dx);
                    p.fy[i] += (force * dy);
                    p.fz[i] += (force * dz);

                    p.fx[j] -= (force * dx);
                    p.fy[j] -= (force * dy);
                    p.fz[j] -= (force * dz);
                }
            }
        }
    }
    *U = (*U) * multiplicateur_potentiel;
}

// Algorithme de velocity-verlet
void velocity_verlet(Particules &p, uint taille, f64 r_cut, f64 *U, bool mode_periodique, uint N_sym)
{
    // Mise à jour des vitesses intermédiaires
    for (uint i = 0; i < taille; i++)
    {
        p.Mx[i] -= (p.fx[i] * DT * CONVERSION_FORCE * 0.5);
        p.My[i] -= (p.fy[i] * DT * CONVERSION_FORCE * 0.5);
        p.Mz[i] -= (p.fz[i] * DT * CONVERSION_FORCE * 0.5);
    }

    // Mise à jour des positions
    for (uint i = 0; i < taille; i++)
    {
        p.x[i] += (p.Mx[i] / MASSE_PARTICULE) * DT;
        p.y[i] += (p.My[i] / MASSE_PARTICULE) * DT;
        p.z[i] += (p.Mz[i] / MASSE_PARTICULE) * DT;
    }

    if (mode_periodique)
    {
        appliquer_conditions_periodiques(p, taille, L, L, L);
    }

    // Mise à jour des forces
    *U = 0.0;
    calcul_energie_force_LJ(p, taille, mode_periodique, U, r_cut, N_sym);
    // Correction des vitesses finales et calcul du moment cinétique

    for (uint i = 0; i < taille; i++)
    {
        // Moment cinétique Moment = Masse * Vitesse
        p.Mx[i] -= (p.fx[i] * DT * CONVERSION_FORCE * 0.5);
        p.My[i] -= (p.fy[i] * DT * CONVERSION_FORCE * 0.5);
        p.Mz[i] -= (p.fz[i] * DT * CONVERSION_FORCE * 0.5);
    }
}

// Mise des positions en tenant compte des conditions périodiques
void appliquer_conditions_periodiques(Particules &p, uint &taille, f64 Lx, f64 Ly, f64 Lz)
{
    for (uint i = 0; i < taille; i++)
    {
        p.x[i] = p.x[i] - std::floor(p.x[i] / Lx) * Lx;
        p.y[i] = p.y[i] - std::floor(p.y[i] / Ly) * Ly;
        p.z[i] = p.z[i] - std::floor(p.z[i] / Lz) * Lz;
    }
}
// Calcul de l'énergie cinétique
void calcul_energie_cinetique_temperature(const Particules &p, uint taille, f64 *EC, f64 *TC)
{
    // Initialisation de l'énergie cinétique et la température
    *TC = 0.0f;

    // Nombre de degrés de liberté
    uint Ndl = 3 * taille - 3;

    // Calcul de l'energie cinématique
    calcul_energie_cinetique(p, taille, EC);

    // Calcul de la température
    *TC = (FACTEUR_TEMPERATURE / Ndl) * (*EC);
}

// Fonction signe
f64 fonction_signe(f64 valeur, f64 s)
{
    return (s < 0.0) ? -valeur : valeur;
}

// Génération des moments cinétiques et recabrage pour T0 = 300 K
void calcul_moments_cinetiques_init(Particules &p, uint taille)
{
    std::srand(42);
    f64 c, s;
    // Initialisation de comments cinématiques
    for (uint i = 0; i < taille; i++)
    {
        c = (f64)rand() / RAND_MAX;
        s = (f64)rand() / RAND_MAX;
        p.Mx[i] = fonction_signe(1.0, 0.5 - s) * c;

        c = (f64)rand() / RAND_MAX;
        s = (f64)rand() / RAND_MAX;
        p.My[i] = fonction_signe(1.0, 0.5 - s) * c;

        c = (f64)rand() / RAND_MAX;
        s = (f64)rand() / RAND_MAX;
        p.Mz[i] = fonction_signe(1.0, 0.5 - s) * c;
    }
}

// calcul de l'énergie cinétique
void calcul_energie_cinetique(const Particules &p, uint &taille, f64 *EC)
{
    *EC = 0.0;
    for (uint i = 0; i < taille; i++)
    {
        *EC = *EC + ((p.Mx[i] * p.Mx[i] + p.My[i] * p.My[i] + p.Mz[i] * p.Mz[i]) / MASSE_PARTICULE);
    }
    *EC *= FACTEUR_CONVERSION;
}

// correction rapport
void correction_rapport(Particules &p, uint &taille)
{
    uint Ndl = 3 * taille - 3;
    f64 energie_cinetique_initial = 0.0;
    calcul_energie_cinetique(p, taille, &energie_cinetique_initial);

    f64 RAPPORT = (Ndl * CONSTANTE_R * T0) / energie_cinetique_initial;

    RAPPORT = sqrt(RAPPORT);

    for (uint i = 0; i < taille; i++)
    {
        p.Mx[i] *= RAPPORT;
        p.My[i] *= RAPPORT;
        p.Mz[i] *= RAPPORT;
    }
}

// Correction du moment cinétique par rapport au centre de masse
void correction_moments_cinetiques(Particules &p, uint taille)
{
    // Calcul du moment cinétique total du centre de masse
    f64 M_moy_x = 0.0, M_moy_y = 0.0, M_moy_z = 0.0;
    for (uint i = 0; i < taille; i++)
    {
        M_moy_x += p.Mx[i];
        M_moy_y += p.My[i];
        M_moy_z += p.Mz[i];
    }

    M_moy_x /= taille;
    M_moy_y /= taille;
    M_moy_z /= taille;

    // Correction des moments cinétiques des particules
    for (uint i = 0; i < taille; i++)
    {
        p.Mx[i] -= M_moy_x;
        p.My[i] -= M_moy_y;
        p.Mz[i] -= M_moy_z;
    }
}

// Correction moment cinétique à l'aide du thermostat de Berendsen
void correction_moment_cinetique_thermostat_berendsen(Particules &p, uint taille, f64 t, f64 gamma, f64 t0)
{
    f64 val_thermostat_berendsen = 0.0;

    if (t0 == 0)
    {
        std::cerr << "Thermostat Berendsen : Attention, la température initiale = 0" << std::endl;
        exit(EXIT_FAILURE);
    }

    val_thermostat_berendsen = -gamma * ((t / t0) - 1);

    for (uint i = 0; i < taille; i++)
    {
        p.Mx[i] += val_thermostat_berendsen * p.Mx[i];
        p.My[i] += val_thermostat_berendsen * p.My[i];
        p.Mz[i] += val_thermostat_berendsen * p.Mz[i];
    }
}

// Fonction de construction de la liste des voisins
void construire_liste_voisins(Particules &p, uint taille, f64 r_cut, uint N_sym)
{
    f64 r_cut_carre = r_cut * r_cut;
    vecteurs_translation vecteurs_t;

    for (uint i = 0; i < taille; i++)
    {
        uint ni = 0; // Indice de voisin pour la particule i
        for (uint j = i + 1; j < taille; j++)
        {
            for (uint i_sym = 0; i_sym < N_sym; i_sym++)
            {
                f64 x_j_loc = p.x[j] + vecteurs_t.x[i_sym];
                f64 y_j_loc = p.y[j] + vecteurs_t.y[i_sym];
                f64 z_j_loc = p.z[j] + vecteurs_t.z[i_sym];

                f64 dx = p.x[i] - x_j_loc;
                f64 dy = p.y[i] - y_j_loc;
                f64 dz = p.z[i] - z_j_loc;

                f64 r_carre = dx * dx + dy * dy + dz * dz;

                if (r_carre <= r_cut_carre)
                {
                    liste_voisins[i][ni * 2] = j;
                    liste_voisins[i][ni * 2 + 1] = i_sym;
                    ni++;

                    if (ni >= max_voisins)
                    {
                        std::cerr << "Erreur : Dépassement de la capacité de la liste des voisins pour la particule " << i << std::endl;
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
    }
}

// Calcul des forces et de l’énergie potentielle avec la liste de Verlet
void calcul_energie_force_LJ_voisins(Particules &p, uint taille, f64 *U, f64 n_epsilon)
{
    init_force(p, taille);
    *U = 0.0;

    vecteurs_translation vecteurs_t;

    for (uint i = 0; i < taille; i++)
    {
        for (uint ni = 0; ni < max_voisins && liste_voisins[i][ni * 2] != -1; ni++) // Utilisation d'une boucle for
        {
            uint j = liste_voisins[i][ni * 2];         // Indice du voisin
            uint i_sym = liste_voisins[i][ni * 2 + 1]; // Indice de translation périodique

            f64 x_j_loc = p.x[j] + vecteurs_t.x[i_sym];
            f64 y_j_loc = p.y[j] + vecteurs_t.y[i_sym];
            f64 z_j_loc = p.z[j] + vecteurs_t.z[i_sym];

            f64 dx = p.x[i] - x_j_loc;
            f64 dy = p.y[i] - y_j_loc;
            f64 dz = p.z[i] - z_j_loc;

            f64 r_carree = dx * dx + dy * dy + dz * dz;
            f64 inv_r2 = 1.0 / r_carree;
            f64 inv_r6 = inv_r2 * inv_r2 * inv_r2;
            f64 inv_r14 = inv_r6 * inv_r6 * inv_r2;

            f64 factor = R_ETOILE_6 * inv_r6;
            *U += n_epsilon * (factor * factor - factor);

            f64 factor1 = R_ETOILE_13 * inv_r14;
            f64 factor2 = R_ETOILE_7 * inv_r6 * inv_r2;
            // Calcul de la force
            f64 force = -24 * n_epsilon * (2 * factor1 - factor2);

            p.fx[i] += force * dx;
            p.fy[i] += force * dy;
            p.fz[i] += force * dz;

            p.fx[j] -= force * dx;
            p.fy[j] -= force * dy;
            p.fz[j] -= force * dz;
        }
    }
    *U *= 2;
}

void sauvegarder_trajectoire_PDB(const Particules &p, uint taille, const std::string &nom_fichier, uint iteration, uint XDIM, uint YDIM, uint ZDIM)
{
    std::ofstream fichier(nom_fichier, std::ios::app);

    if (!fichier.is_open())
    {
        std::cerr << "Erreur lors de l'ouverture du fichier " << nom_fichier << std::endl;
        return;
    }

    // Écrire les lignes de tête CRYST1 et MODEL
    fichier << "CRYST1" << std::setw(9) << XDIM << std::setw(9) << YDIM << std::setw(9) << ZDIM
            << "  90.00  90.00  90.00 P\n";
    fichier << "MODEL" << std::setw(10) << iteration << "\n";

    // Écrire les coordonnées des particules
    for (uint i = 0; i < taille; ++i)
    {
        fichier << "ATOM  " << std::setw(6) << i << "  C           0   "
                << std::fixed << std::setprecision(3)
                << std::setw(8) << p.x[i]
                << std::setw(8) << p.y[i]
                << std::setw(8) << p.z[i]
                << "                  MRES\n";
    }

    // Écrire les lignes de fin TER et ENDMDL
    fichier << "TER\n";
    fichier << "ENDMDL\n";

    fichier.close();
}

// Initialisation du moment avec les vitesses
void init_moment_avec_vitesse(Particules &p, const std::string &nom_fichier)
{
    std::ifstream fichier(nom_fichier);
    if (!fichier)
    {
        std::cerr << "Erreur : Impossible d'ouvrir le fichier " << nom_fichier << std::endl;
        return;
    }
    f64 Vx, Vy, Vz;
    size_t i = 0;
    while (fichier >> Vx >> Vy >> Vz)
    {
        if (i >= p.x.size())
        {
            std::cerr << "Erreur : Plus de vitesses que de particules enregistrées." << std::endl;
            break;
        }
        // Mise à jour du moment linéaire P = m * v
        p.Mx[i] = MASSE_PARTICULE * Vx;
        p.My[i] = MASSE_PARTICULE * Vy;
        p.Mz[i] = MASSE_PARTICULE * Vz;
        i++;
    }

    fichier.close();
}

// Sauvegarde des données de la simulation
void sauvegarder_donnee_simulation(const std::string &nom_fichier, const f64 &temperature, const f64 &energie_cinetique, const f64 &U, uint iteration)
{
    std::ofstream fichier(nom_fichier, std::ios::app);

    if (!fichier)
    {
        std::cerr << "Erreur lors de l'ouverture du fichier " << nom_fichier << std::endl;
        exit(EXIT_FAILURE);
    }

    fichier << "# Iteration :       " << iteration << "\n";
    fichier << "#   Temp :          " << temperature << "\n";
    fichier << "#   Etot :          " << std::fixed << std::setprecision(4)
            << (energie_cinetique + U)
            << "     Epot :          " << U
            << "     Ecin :          " << energie_cinetique << "\n";

    fichier.close();
}

// Initialisation des force à 0
void init_force(Particules &p, uint taille)
{
    for (uint i = 0; i < taille; i++)
    {
        p.fx[i] = 0.0;
        p.fy[i] = 0.0;
        p.fz[i] = 0.0;
    }
}

// Initialisation des forces à 0
void initialisation(Particules &p, uint taille)
{
    p.fx.resize(taille);
    p.fy.resize(taille);
    p.fz.resize(taille);
    p.Mx.resize(taille);
    p.My.resize(taille);
    p.Mz.resize(taille);
    for (uint i = 0; i < taille; i++)
    {
        p.Mx[i] = 0.0;
        p.My[i] = 0.0;
        p.Mz[i] = 0.0;
    }
}

// Verification de la nullitée des forces
void verifier_valeur_force(const std::vector<f64> &fx, const std::vector<f64> &fy, const std::vector<f64> &fz, uint taille)
{
    f64 total_x = 0.0;
    f64 total_y = 0.0;
    f64 total_z = 0.0;
    for (uint i = 0; i < taille; i++)
    {
        total_x += fx[i];
        total_y += fy[i];
        total_z += fz[i];
    }

    std::cout << "AFFICHAGE DE LA VALEUR TOTALE DES FORCES!!" << std::endl;
    std::cout << " Fx : " << total_x << std::endl;
    std::cout << " Fy : " << total_y << std::endl;
    std::cout << " Fz : " << total_z << std::endl;
}




/* Fonction utilitaire */

//Fonction permettant parser les arguments
void parse_arguments(int argc, char* argv[]) {
    int opt;
    while ((opt = getopt(argc, argv, "h:f:p:i:s:o:r:d:e:t:l:g:m:c:R:")) != -1) {
        switch (opt) {
            case 'h': afficher_aide(); break;
            case 'f': nom_fichier = optarg; break;
            case 'p': N_particules_local = static_cast<uint>(std::stoul(optarg)); break;
            case 'i': n_iterations = static_cast<uint>(std::stoul(optarg)); break;
            case 's': m_step = static_cast<uint>(std::stoul(optarg)); break; 
            case 'o': mode_periodique = std::stoi(optarg) != 0; break;
            case 'r': r_cut = std::stod(optarg); break; 
            case 'd': DT = std::stod(optarg); break;
            case 'e': MY_EPSILON = std::stod(optarg); break;
            case 't': T0 = std::stod(optarg); break;
            case 'l': L = std::stod(optarg); break;
            case 'g': GAMMA = std::stod(optarg); break;
            case 'm': MASSE_PARTICULE = std::stod(optarg); break;
            case 'c': CONVERSION_FORCE = std::stod(optarg); break;
            case 'R': CONSTANTE_R = std::stod(optarg); break;
            default:
                std::cerr << "Utilisation incorrecte. Tapez ./main -h pour l'aide.\n";
                exit(EXIT_FAILURE);
        }
    }
}

//Fonction permettant d'afficher l'aide
void afficher_aide() {
    std::cout << "Usage: ./main [options]\n";
    std::cout << "Options disponibles :\n";
    std::cout << "  -h           Afficher cette aide\n";
    std::cout << "  -f <nom>     Nom du fichier contenant les particules (par défaut: particule.xyz)\n";
    std::cout << "  -p <val>     Nombre de particules locales (par défaut: 1000)\n";
    std::cout << "  -i <val>     Nombre d'itérations (par défaut: 1000)\n";
    std::cout << "  -s <val>     Nombre de pas pour la correction (par défaut: 10)\n";
    std::cout << "  -o <0|1>    Mode périodique (1 = activé, 0 = désactivé) (par défaut: 0)\n";
    std::cout << "  -r <val>    Rayon de coupure (par défaut: 10.0)\n";
    std::cout << "  -d <val>     Pas de temps en femtosecondes (par défaut: 1.0)\n";
    std::cout << "  -e <val>     Profondeur du puits de potentiel epsilon (par défaut: 0.2)\n";
    std::cout << "  -t <val>     Température initiale en Kelvin (par défaut: 300.0)\n";
    std::cout << "  -l <val>     Longueur de la boîte de simulation (par défaut: 42.0)\n";
    std::cout << "  -g <val>     Facteur gamma du thermostat Berendsen (par défaut: 0.01)\n";
    std::cout << "  -m <val>     Masse des particules (par défaut: 18.0)\n";
    std::cout << "  -c <val>    Facteur de conversion de la force (par défaut: 0.0001 * 4.186)\n";
    std::cout << "  -R <val>     Constante R en kcal/(mol·K) (par défaut: 0.00199)\n";
    exit(0);
}
