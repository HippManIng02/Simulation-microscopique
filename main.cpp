#include <vector>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include "function.h"

int main(int argc, char *argv[])
{
    //std::cout << std::setprecision(16);
    f64 r_cut = 10.0;
    bool mode_periodique = false;
    std::string nom_fichier;
    unsigned int N_particules_local;
    unsigned int N_sym = 27;

    if (argc < 3) 
    {
        fprintf(stderr, "Utilisez: %s particule.xyz N_particules_local mode_periodique(1|0)\n", argv[0]);
        return EXIT_FAILURE;
    }

    // Lecture des arguments
    nom_fichier = argv[1];
    N_particules_local = atoi(argv[2]);
    mode_periodique = atoi(argv[3]); // 1 = périodique, 0 = non périodique

    // Affichage des paramètres d'entrée
    std::cout << "Mode périodique : " << (mode_periodique ? "Oui" : "Non") << std::endl;

    // Vérifier le nombre de particules dans le fichier
    unsigned int nbr_particules = get_nombre_particules(nom_fichier);
    std::cout << "Nombre de particules : " << nbr_particules << std::endl;

    if (nbr_particules > N_particules_total)
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
        f64 x, y, z;
        int tmp;
        ss >> tmp >> x >> y >> z;
        particules.x.push_back(x);
        particules.y.push_back(y);
        particules.z.push_back(z);
    }

    fichier.close();

    // Initialisation des forces
    f64 U = 0.0;
    init_force(particules, N_particules_local);

    // Calcul de l'énergie et des forces selon le mode

    if (mode_periodique)
    {
        std::cout << "Calcul avec conditions périodiques..." << std::endl;
    }
    else
    {
        std::cout << "Calcul sans conditions périodiques..." << std::endl;
    }

    calcul_energie_force_LJ(particules, N_particules_local, mode_periodique, &U, r_cut, N_sym);

    // Afficher les résultats
    std::cout << "Énergie totale du système : " << U << std::endl;

    // Vérifier que la somme des forces est proche de zéro
    verifier_valeur_force(particules.fx, particules.fy, particules.fz, N_particules_local);

    return EXIT_SUCCESS;
}
