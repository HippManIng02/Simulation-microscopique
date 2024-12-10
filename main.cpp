#include <vector>
#include <sstream>
#include <iomanip>
#include "function.h"

int main(int argc, char *argv[])
{
    //std::cout << std::setprecision(16);
    if (argc < 2)
    {
        fprintf(stderr, "Utilisez: %s particule.xyz N_particules_local\n", argv[0]);
        return EXIT_FAILURE;
    }

    // Déclaration des variables
    std::string nom_fichier = argv[1];
    unsigned int N_particules_local = atoi(argv[2]);

    unsigned int nbr_particules = getNbrParticules(nom_fichier);
    std::cout << "nbr particule :" << nbr_particules << std::endl;
    if (nbr_particules > N_particules_total)
    {
        std::cerr << "Attention le nombre de particule est supérieur à 1000" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (N_particules_local > nbr_particules)
    {
        std::cout << "Le nombre de particule local doit être inférieur à " << nbr_particules << std::endl;
        return EXIT_FAILURE;
    }

    Particules particules;
    std::string ligne;
    // Ouverture et lecture des articules dans le fichier d'entrée
    std::ifstream fichier(nom_fichier);
    // Activation de l'exception
    if (!fichier.is_open())
    {
        std::cerr << "Une erreur est survenue lors de l'ouverture du fichier pour le comptage du nombre de particule" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::getline(fichier, ligne);
    
    //variable temporaire contenant les types de particules, sans importance dans notre cas.
    while (std::getline(fichier, ligne))
    {
        std::stringstream ss(ligne);
        f64 x, y, z;
        int tmp;
        ss >> tmp >> x >> y >> z;
        particules.x.push_back(x);
        particules.y.push_back(y);
        particules.z.push_back(z);
        //std::cout<<"x: "<<x<< " y: "<<y<<" z: "<<z<<std::endl;
    }
    
    std::cout<<"x : "<<particules.x.back()<<"y :"<<particules.y.back()<<"z :"<<particules.z.back()<<std::endl;
    
    fichier.close();

    f64 U = 0.0;
    init_force(particules, N_particules_local);
    calcul_energie_force_LJ(particules, N_particules_local, &U);
    std::cout<<"Energie totale du système: " << U <<std::endl;
    //std::cout<<"Energie totale du système (LJ) : " << calcul_potentiel_LJ(particules, N_particules_local)<<std::endl;
    //calcul_force_LJ(particules, N_particules_local);
    verifier_valeur_force(particules.fx, particules.fy, particules.fz, N_particules_local);

    return EXIT_SUCCESS;
}
