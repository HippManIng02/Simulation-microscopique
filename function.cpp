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
f64 calcul_distance_carree(f64 x_1, f64 y_1, f64 z_1, f64 x_2, f64 y_2, f64 z_2)
{
    f64 x = x_1 - x_2;
    f64 y = y_1 - y_2;
    f64 z = z_1 - z_2;
    //std::cout << " r.x :" << x << " | r.y :" << y << " | r.z :" << z << std::endl;
    return (x * x) + (y * y) + (z * z);
}

// calcul du terme (r*)²/(r)²
f64 calcul_de_r_etoile_div_par_r__2(f64 r_carree, f64 r_e_carree)
{
    if (r_carree == 0)
    {
        std::cerr << "Attention, la valeur de r² = 0" << std::endl;
        exit(EXIT_FAILURE);
    }
    return r_e_carree / r_carree;
}

// Calcul du potentiel de Lennard jones
f64 potentiel_lennard_jones(f64 r_carree, f64 r_e_carree, f64 n_epsilon)
{
    f64 tmp_2 = calcul_de_r_etoile_div_par_r__2(r_carree, r_e_carree);
    f64 tmp_6 = tmp_2 * tmp_2 * tmp_2;   // ((r*)/(r))⁶
    f64 tmp_12 = tmp_6 * tmp_6;          //((r*)²/(r)²)⁶ * ((r*)²/(r)²)⁶ = ((r*)²/(r)²)¹²
    return n_epsilon * (tmp_12 - tmp_6); // 4 * epsilon * [(r*)²/(r)²)¹² - 2* ((r*)²/(r)²)⁶]
}

// Calcul de la force de Lennard jones, (dérivée du potentiel)
f64 force_lennard_jones(f64 r_carree, f64 r_e_carree, f64 n_epsilon)
{
    f64 tmp_2 = calcul_de_r_etoile_div_par_r__2(r_carree, r_e_carree); // (r*)²/(r)²
    f64 tmp_6 = tmp_2 * tmp_2 * tmp_2;                                 // ((r*)/(r))⁶
    f64 tmp_8 = tmp_6 * tmp_2;                                         // ((r*)²/(r)²)⁴
    f64 tmp_14 = tmp_6 * tmp_6 * tmp_2;                                //(r*)²/(r)²)¹⁴
    return -48 * n_epsilon * (tmp_14 - tmp_8);                         // -48 * epsilon * [(r*)²/(r)²)¹⁴ - ((r*)²/(r)²)⁸]
}

// Calcul de l'énergie et de la force en utilisant le potentiel de Lennard Jones
void calcul_energie_force_LJ(Particules &p, uint taille, bool mode_periodique, f64 *U, f64 r_cut, uint N_sym, f64 r_e_carree, f64 n_epsilon)
{
    vecteurs_translation vecteurs_t;
    const uint n_sym_effectif = mode_periodique ? N_sym : 1;
    const f64 r_cut_effectif = mode_periodique ? r_cut : std::numeric_limits<f64>::max();
    const f64 multiplicateur_potentiel = mode_periodique ? 2.0 : 4.0;

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
                    f64 x_j = p.x[j] + (mode_periodique ? vecteurs_t.x[i_sym] : 0.0);
                    f64 y_j = p.y[j] + (mode_periodique ? vecteurs_t.y[i_sym] : 0.0);
                    f64 z_j = p.z[j] + (mode_periodique ? vecteurs_t.z[i_sym] : 0.0);

                    f64 r_carree = calcul_distance_carree(x, y, z, x_j, y_j, z_j);

                    if (r_carree < r_cut_effectif * r_cut_effectif)
                    {
                        *U += potentiel_lennard_jones(r_carree, r_e_carree, n_epsilon);

                        f64 force = force_lennard_jones(r_carree, r_e_carree, n_epsilon);

                        f64 dx = x - x_j;
                        f64 dy = y - y_j;
                        f64 dz = z - z_j;

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
    *U *= multiplicateur_potentiel;
}


//Calcul de l'énergie en utilisant le potentiel de Lennard Jones
void calcul_energie_LJ(Particules& p, uint taille, bool mode_periodique, f64 *U,f64 r_cut, uint N_sym,f64 r_e_carree, f64 n_epsilon){
    
    vecteurs_translation vecteurs_t;
    const uint n_sym_effectif = mode_periodique ? N_sym : 1;
    const f64 r_cut_effectif = mode_periodique ? r_cut : std::numeric_limits<f64>::max();
    const f64 multiplicateur_potentiel = mode_periodique ? 2.0 : 4.0;

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
                    f64 x_j = p.x[j] + (mode_periodique ? vecteurs_t.x[i_sym] : 0.0);
                    f64 y_j = p.y[j] + (mode_periodique ? vecteurs_t.y[i_sym] : 0.0);
                    f64 z_j = p.z[j] + (mode_periodique ? vecteurs_t.z[i_sym] : 0.0);

                    f64 r_carree = calcul_distance_carree(x, y, z, x_j, y_j, z_j);

                    if (r_carree < r_cut_effectif * r_cut_effectif)
                    {
                        *U += potentiel_lennard_jones(r_carree, r_e_carree, n_epsilon);
                    }
                }
            }
        }
    }
    *U *= multiplicateur_potentiel;
}

// Algorithme de velocity-verlet
void velocity_verlet(Particules &p, uint taille, f64 r_cut, uint N_sym, f64 r_e_carree, f64 n_epsilon)
{
    // Mise à jour des vitesses intermédiaires
    // f64 inv_masse = 1.0 / MASSE_PARTICULE;

    // f64 dt_par_2 = dt * 0.5;
    // f64 tmp =  inv_masse * dt_par_2;

    for (uint i = 0; i < taille; i++)
    {
        p.vx[i] += ((p.fx[i] * dt * CONVERSION_FORCE) / (MASSE_PARTICULE * 2.0f));
        p.vy[i] += ((p.fy[i] * dt * CONVERSION_FORCE) / (MASSE_PARTICULE * 2.0f));
        p.vz[i] += ((p.fz[i] * dt * CONVERSION_FORCE) / (MASSE_PARTICULE * 2.0f));
    }

    // Mise à jour des positions
    for (uint i = 0; i < taille; i++)
    {
        p.x[i] += p.vx[i] * dt;
        p.y[i] += p.vy[i] * dt;
        p.z[i] += p.vz[i] * dt;
    }

    // Mise à jour des forces
    f64 U_temp = 0.0f;

    init_force(p, taille);
    calcul_energie_force_LJ(p, taille, 1, &U_temp, r_cut, N_sym, r_e_carree, n_epsilon);

    // Correction des vitesses finales et calcul du moment cinétique
    // f64 masse_mul_conv_force = MASSE_PARTICULE * CONVERSION_FORCE;
    for (uint i = 0; i < taille; i++)
    {
        p.vx[i] += ((p.fx[i] * dt * CONVERSION_FORCE) / (MASSE_PARTICULE * 2.0f));
        p.vy[i] += ((p.fy[i] * dt * CONVERSION_FORCE) / (MASSE_PARTICULE * 2.0f));
        p.vz[i] += ((p.fz[i] * dt * CONVERSION_FORCE) / (MASSE_PARTICULE * 2.0f));

        // Moment cinétique Moment = Masse * Vitesse
        p.Mx[i] = p.vx[i] * MASSE_PARTICULE;
        p.My[i] = p.vy[i] * MASSE_PARTICULE;
        p.Mz[i] = p.vz[i] * MASSE_PARTICULE;
    }
}

// Calcul de l'énergie cinétique
void calcul_energie_cinetique(const Particules &p, uint taille, f64 *EC, f64 *TC)
{
    // Initialisation de l'énergie cinétique et la température
    *EC = 0.0f;
    *TC = 0.0f;

    // Nombre de degrés de liberté
    int Ndl = 3 * taille - 3;

    // Calcul de l'energie cinématique
    for (uint i = 0; i < taille; i++)
    {
        *EC += (p.Mx[i] * p.Mx[i] + p.My[i] * p.My[i] + p.Mz[i] * p.Mz[i]);
    }

    // Application du facteur de conversion
    *EC *= (1.0f / (CONVERSION_FORCE * 2.0f * MASSE_PARTICULE));

    // Calcul de la température
    *TC = (1.0f / (Ndl * CONSTANTE_R)) * (*EC);
}

// Fonction signe
f64 fonction_signe(f64 valeur, f64 s)
{
    return (s < 0.5) ? valeur : -valeur;
}

// Génération des moments cinétiques et recabrage pour T0 = 300 K
void calcul_moments_cinetiques(Particules &p, uint taille)
{
    std::srand(42);
    double c, s;
    // Initialisation de comments cinématiques
    for (uint i = 0; i < taille; i++)
    {
        c = (f64) rand()/RAND_MAX;
        s = (f64) rand()/RAND_MAX;
        p.Mx[i] = fonction_signe(1.0f, 0.5f - s) * c;

        c = (f64) rand()/RAND_MAX;
        s = (f64) rand()/RAND_MAX;
        p.My[i] = fonction_signe(1.0f, 0.5f - s) * c;

        c = (f64) rand()/RAND_MAX;
        s = (f64) rand()/RAND_MAX;
        p.Mz[i] = fonction_signe(1.0f, 0.5f - s) * c;
 
    }

    // calcul de l'énergie cinétique initiale
    f64 energie_cinetique_initial = 0.0;
    int Ndl = 3 * taille - 3;

    for (uint i = 0; i < taille; i++)
    {
        energie_cinetique_initial += (p.Mx[i] * p.Mx[i] + p.My[i] * p.My[i] + p.Mz[i] * p.Mz[i]);
    }

    // Application du facteur de conversion
    energie_cinetique_initial *= 1.0 / (2.0 * CONVERSION_FORCE * MASSE_PARTICULE);

    // Calcul du facteur de correction
    f64 RAPPORT = (Ndl * CONSTANTE_R * T0) / energie_cinetique_initial;

    RAPPORT = sqrt(RAPPORT);

    for (uint i = 0; i < taille; i++)
    {
        p.Mx[i] *= RAPPORT;
        p.My[i] *= RAPPORT;
        p.Mz[i] *= RAPPORT;
    }

    // for (uint i = 0; i < taille; i++)
    // {
    //     std::cout<<" My[] = " << p.My[i] << std::endl;
    //     std::cout<<" Mz[] = " << p.Mz[i] << std::endl;
    // }
    
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
    f64 val_thermostat_berendsen = gamma * ((t0 / t) - 1);
    for(uint i = 0; i < taille ; i++){
        p.Mx[i] += val_thermostat_berendsen * p.Mx[i];
        p.My[i] += val_thermostat_berendsen * p.My[i];
        p.Mz[i] += val_thermostat_berendsen * p.Mz[i];
    }
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
void initialisation(Particules &p, unsigned int taille)
{
    p.fx.resize(taille);
    p.fy.resize(taille);
    p.fz.resize(taille);
    p.Mx.resize(taille);
    p.My.resize(taille);
    p.Mz.resize(taille);
    p.vx.resize(taille);
    p.vy.resize(taille);
    p.vz.resize(taille);
    for (uint i = 0; i < taille; i++)
    {
        p.fx[i] = 0.0;
        p.fy[i] = 0.0;
        p.fz[i] = 0.0;
        p.Mx[i] = 0.0;
        p.My[i] = 0.0;
        p.Mz[i] = 0.0;
        p.vx[i] = 0.0;
        p.vy[i] = 0.0;
        p.vz[i] = 0.0;
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