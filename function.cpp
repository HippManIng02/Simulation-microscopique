#include "function.h"
#include <limits>
#include <cmath>

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

//Calcul de l'énergie et de la force en utilisant le potentiel de Lennard Jones
void calcul_energie_force_LJ(Particules &p, uint taille, bool mode_periodique, f64 *U, f64 r_cut, uint N_sym, f64 r_e_carree, f64 n_epsilon)
{
    init_force(p, taille);

    vecteurs_translation vecteurs_t;
    const uint n_sym_effectif = mode_periodique == 1 ? N_sym : 1;
    const f64 r_cut_effectif = mode_periodique == 1 ? r_cut : std::numeric_limits<f64>::max();
    const f64 multiplicateur_potentiel = mode_periodique == 1 ? 2.0 : 4.0;

    for (uint i_sym = 0; i_sym < n_sym_effectif; i_sym++)
    {
        for (uint i = 0; i < taille; i++)
        {
            f64 x = p.x[i];
            f64 y = p.y[i];
            f64 z = p.z[i];
            for (uint j = 0; j < taille; j++)
            {
                if (i != j)
                {
                   // std::cout<< "x : " << vecteurs_t.x[i_sym] << "y : " <<  vecteurs_t.y[i_sym] << "z :" << vecteurs_t.z[i_sym] << std::endl;
                    f64 x_j = p.x[j] +  vecteurs_t.x[i_sym];
                    f64 y_j = p.y[j] +  vecteurs_t.y[i_sym];
                    f64 z_j = p.z[j] +  vecteurs_t.z[i_sym]; 

                   // std::cout<< "x : " << x_j << "y : " <<  y_j << "z :" << z_j << std::endl;
                   f64 dx = x - x_j;
                   f64 dy = y - y_j;
                   f64 dz = z - z_j;

                    f64 r_carree = calcul_distance_carree(dx, dy, dz);

                    if (r_carree <= r_cut_effectif * r_cut_effectif)
                    {

                        f64 r_3 = r_carree * r_carree * r_carree; // (r)^3
                        f64 r_6 = r_3 * r_3; // (r)^6
                        std::cout << std::fixed << std::setprecision(15);

                        *U += n_epsilon * ((r_etoile_12 / r_6) - (r_etoile_6 / r_3));

                        f64 force = -48.0 * n_epsilon * (((r_etoile_12 * r_e_carree)/(r_6 * r_carree)) - ((r_etoile_6 * r_e_carree)/(r_3 * r_carree))); 

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
    }
    *U = (*U) * multiplicateur_potentiel;
}


// Algorithme de velocity-verlet
void velocity_verlet(Particules &p, uint taille, f64 r_cut, f64 *U, bool mode_periodique,uint N_sym, f64 r_e_carree, f64 n_epsilon)
{
    // Mise à jour des vitesses intermédiaires
    for (uint i = 0; i < taille; i++)
    {
        p.Mx[i] -= (p.fx[i] * dt * CONVERSION_FORCE * 0.5);
        p.My[i] -= (p.fy[i] * dt * CONVERSION_FORCE * 0.5);
        p.Mz[i] -= (p.fz[i] * dt * CONVERSION_FORCE * 0.5);

        // Mise à jour des positions
        p.x[i] += (p.Mx[i] / MASSE_PARTICULE) * dt;
        p.y[i] += (p.My[i] / MASSE_PARTICULE) * dt;
        p.z[i] += (p.Mz[i] / MASSE_PARTICULE) * dt;
    }
    // Mise à jour des forces
    *U = 0.0;
    calcul_energie_force_LJ(p, taille, mode_periodique, U, r_cut, N_sym);
    // Correction des vitesses finales et calcul du moment cinétique

    for (uint i = 0; i < taille; i++)
    {
        // Moment cinétique Moment = Masse * Vitesse
        //std::cout<< "Mx : " << p.Mx[i] << " My : " <<  p.My[i] << " Mz :" << p.Mz[i] << std::endl;
        p.Mx[i] -= (p.fx[i] * dt * CONVERSION_FORCE * 0.5);
        p.My[i] -= (p.fy[i] * dt * CONVERSION_FORCE * 0.5);
        p.Mz[i] -= (p.fz[i] * dt * CONVERSION_FORCE * 0.5);
        //std::cout<< "Mx2 : " << p.Mx[i] << " My2 : " <<  p.My[i] << " Mz2 :" << p.Mz[i] << std::endl;
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
    *TC = (1.0 / (Ndl * CONSTANTE_R)) * (*EC);
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
        c = (f64) rand()/RAND_MAX;
        s = (f64) rand()/RAND_MAX;
        p.Mx[i] = fonction_signe(1.0, 0.5 - s) * c;

        c = (f64) rand()/RAND_MAX;
        s = (f64) rand()/RAND_MAX;
        p.My[i] = fonction_signe(1.0, 0.5 - s) * c;

        c = (f64) rand()/RAND_MAX;
        s = (f64) rand()/RAND_MAX;
        p.Mz[i] = fonction_signe(1.0, 0.5 - s) * c;
    }
}


//calcul de l'énergie cinétique
void calcul_energie_cinetique(const Particules &p, uint &taille, f64 *EC){
    *EC = 0.0;
    for (uint i = 0; i < taille; i++)
    {
        *EC = *EC + ((p.Mx[i] * p.Mx[i] + p.My[i] * p.My[i] + p.Mz[i] * p.Mz[i]) / MASSE_PARTICULE);
    }
    *EC = *EC * (1.0 / (CONVERSION_FORCE * 2.0));
}

//correction rapport
void correction_rapport(Particules &p, uint &taille){
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

//Correction moment cinétique à l'aide du thermostat de Berendsen
void correction_moment_cinetique_thermostat_berendsen(Particules &p, uint taille, f64 t, f64 gamma, f64 t0){
    f64 val_thermostat_berendsen = 0.0;

    if (t0 == 0)
    {
        std::cerr << "Thermostat Berendsen : Attention, la température initiale = 0" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    val_thermostat_berendsen = -gamma * ((t / t0) - 1);
    
    for(uint i = 0; i < taille ; i++){
        p.Mx[i] += val_thermostat_berendsen * p.Mx[i];
        p.My[i] += val_thermostat_berendsen * p.My[i];
        p.Mz[i] += val_thermostat_berendsen * p.Mz[i];
    }
}

void sauvegarder_trajectoire_PDB(const Particules &p, uint taille, const std::string& nom_fichier, uint iteration, uint XDIM, uint YDIM, uint ZDIM) {
    std::ofstream fichier(nom_fichier, std::ios::app); 

    if (!fichier.is_open()) {
        std::cerr << "Erreur lors de l'ouverture du fichier " << nom_fichier << std::endl;
        return;
    }

    // Écrire les lignes de tête CRYST1 et MODEL
    fichier << "CRYST1 " << XDIM << " " << YDIM << " " << ZDIM << " 90.00 90.00 90.00 P 1\n";
    fichier << "MODEL " << iteration << "\n";

    // Écrire les coordonnées des particules
    for (uint i = 0; i < taille; ++i) {
        fichier << "ATOM  " << std::setw(5) << i + 1 << "  C   0     "
                << std::fixed << std::setprecision(3)
                << std::setw(8) << p.x[i]
                << std::setw(8) << p.y[i]
                << std::setw(8) << p.z[i]
                << "  1.00  0.00           C\n";
    }

    // Écrire les lignes de fin TER et ENDMDL
    fichier << "TER\n";
    fichier << "ENDMDL\n";

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