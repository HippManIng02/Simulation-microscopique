#include <vector>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <set>
#include <tuple>
#include "function.h"

int main(int argc, char *argv[])
{
    // std::cout << std::setprecision(16);
    //std::string nom_fichier;

    // if (argc != 6)
    // {
    //     fprintf(stderr, "Utilisez: %s particule.xyz N_particules_local nbr_iter m_step mode_periodique(1|0)\n", argv[0]);
    //     return EXIT_FAILURE;
    // }

    parse_arguments(argc, argv);

    std::cout << "Paramètres après parsing :\n";
    std::cout << "  Fichier : " << nom_fichier << "\n";
    std::cout << "  N_particules_local = " << N_particules_local << "\n";
    std::cout << "  n_iterations = " << n_iterations << "\n";
    std::cout << "  m_step = " << m_step << "\n";
    std::cout << "  Mode périodique = " << (mode_periodique ? "Activé" : "Désactivé") << "\n";
    std::cout << "  r_cut = " << r_cut << "\n";
    std::cout << "\nOptions facultatives :\n";
    std::cout << "  dt = " << DT << "\n";
    std::cout << "  epsilon = " << MY_EPSILON << "\n";
    std::cout << "  T0 = " << T0 << " K\n";
    std::cout << "  L = " << L << " nm\n";
    std::cout << "  gamma = " << GAMMA << "\n";
    std::cout << "  Masse particule = " << MASSE_PARTICULE << "\n";
    std::cout << "  Conversion force = " << CONVERSION_FORCE << " (kcal/mol)\n";
    std::cout << "  Constante R = " << CONSTANTE_R << " (kcal/(mol·K))\n";

    // Lecture des arguments
    // nom_fichier = argv[1];
    // N_particules_local = std::stoi(argv[2]);

    // n_iterations = std::stoi(argv[3]) > 1 ? std::stoi(argv[3]) : 1000;
    // m_step = std::stoi(argv[4]);
    // mode_periodique = (std::stoi(argv[5]) != 0); // 1 = périodique, 0 = non périodique

    // Affichage des paramètres d'entrée
    std::cout << "Mode périodique : " << (mode_periodique ? "Oui" : "Non") << std::endl;
    std::cout << "m step : " << m_step << std::endl;

    // Vérifier le nombre de particules dans le fichier
    uint nbr_particules = get_nombre_particules(nom_fichier);
    std::cout << "Nombre de particules : " << nbr_particules << std::endl;

    if (nbr_particules > N_PARTICULES_TOTAL)
    {
        std::cerr << "Attention, le nombre de particules dépasse la limite de 1000." << std::endl;
        return EXIT_FAILURE;
    }

    if (N_particules_local > nbr_particules)
    {
        std::cerr << "Le nombre de particules locales doit être inférieur ou égal à " << nbr_particules << std::endl;
        return EXIT_FAILURE;
    }

    // Lecture des particules depuis le fichier
    Particules particules;
    std::ifstream fichier(nom_fichier);

    if (!fichier.is_open())
    {
        std::cerr << "Erreur lors de l'ouverture du fichier : " << nom_fichier << std::endl;
        return EXIT_FAILURE;
    }

    std::string ligne;
    std::getline(fichier, ligne); // Ignorer la première ligne (nombre de particules)

    while (std::getline(fichier, ligne))
    {
        std::stringstream ss(ligne);
        ss >> std::setprecision(15);
        f64 x, y, z;
        int tmp;
        ss >> tmp >> x >> y >> z;
        particules.x.push_back(x);
        particules.y.push_back(y);
        particules.z.push_back(z);
    }

    fichier.close();

    // Initialisation des tableaux
    initialisation(particules, N_particules_local);

    // Calcul de l'énergie et des forces selon le mode

    if (mode_periodique)
    {
        std::cout << "Calcul avec conditions périodiques..." << std::endl;
    }
    else
    {
        std::cout << "Calcul sans conditions périodiques..." << std::endl;
    }

    // Calcul de moment cinétique initiale
    calcul_moments_cinetiques_init(particules, N_particules_local);
    // correction rapport
    correction_rapport(particules, N_particules_local);
    // Correction du moment cinétique pour la conservation du centre de masse
    correction_moments_cinetiques(particules, N_particules_local);
    // correction rapport
    correction_rapport(particules, N_particules_local);

    // init_moment_avec_vitesse(particules, "velocities_init.txt");
    // correction_moments_cinetiques(particules, N_particules_local);
    // correction_rapport(particules, N_particules_local);

    // f64 Px = 0.0, Py = 0.0, Pz = 0.0;
    // for (uint i = 0; i < N_particules_local; i++)
    // {
    //     Px += particules.Mx[i];
    //     Py += particules.My[i];
    //     Pz += particules.Mz[i];
    // }
    // std::cout << "Px = " << Px << " Py = " << Py << " Pz = " << Pz << std::endl;

    f64 U = 0.0;
    // Boucle de simulation
    f64 EC, TC;

    calcul_energie_cinetique_temperature(particules, N_particules_local, &EC, &TC);

    // Mise à jour des forces et de l'énergie potentielle
    calcul_energie_force_LJ(particules, N_particules_local, mode_periodique, &U, r_cut, N_sym);

    std::cout << " | E_cinétique = " << EC
              << " | E_total = " << U + EC
              << " | Température = " << TC << " K"
              << std::endl;

    for (uint iter = 1; iter <= n_iterations; iter++)
    {
        // correction après m_step itération
        if (iter % m_step == 0)
        {
            correction_moment_cinetique_thermostat_berendsen(particules, N_particules_local, TC, GAMMA, T0);
        }
        // Intégration avec Velocity-Verlet
        velocity_verlet(particules, N_particules_local, r_cut, &U, mode_periodique, N_sym);

        // Calcul de l'énergie cinétique
        calcul_energie_cinetique_temperature(particules, N_particules_local, &EC, &TC);

        f64 E_totale = U + EC;

        std::cout << "Iteration " << iter
                  << " | E_potentiel = " << U
                  << " | E_cinétique = " << EC
                  << " | Température = " << TC << " K"
                  << " | E_totale = " << E_totale << std::endl;

        // Vérifier que la somme des forces est proche de zéro
        verifier_valeur_force(particules.fx, particules.fy, particules.fz, N_particules_local);
        sauvegarder_trajectoire_PDB(particules, N_particules_local, "trajectoires.pdb", iter, L, L, L);
    }

    return EXIT_SUCCESS;
}
