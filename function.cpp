#include "function.h"

//Fonction pour compter le nombre de particules
unsigned int getNbrParticules(const std::string& nom_fichier){
    std::ifstream fichier(nom_fichier);
    unsigned int compteur = 0;
    if(!fichier.is_open()){
        std::cerr<<"Une erreur est survenue lors de l'ouverture du fichier pour le comptage du nombre de particule"<<std::endl;
        exit(EXIT_FAILURE);
    }
    
    //Compter le nombre de ligne dans le fichier
    std::string tmp;
    while(std::getline(fichier, tmp))
        compteur ++;
    
    //Retrait de la première ligne du fichier qui sert d'indication
    compteur > 0 && compteur --;

    fichier.close();
    return compteur;
}


//calcul de distance au carrée
f64 calcul_distance_carree(f64 x_1, f64 y_1, f64 z_1, f64 x_2, f64 y_2, f64 z_2){
    f64 x = x_1 - x_2;
    f64 y = y_1 - y_2;
    f64 z = z_1 - z_2;
    return (x*x) + (y*y) + (z*z);
}


//calcul du terme (r*)²/(r)²
f64 calcul_de_r_etoile_div_par_r__2(f64 r_carree, f64 r_e_carree){
    if( r_carree == 0){
        std::cerr<<"Attention, la valeur de r² = 0"<<std::endl;
        exit(EXIT_FAILURE); 
    }
    return r_e_carree / r_carree;
}

// //calcul du terme ((r*)²/(r)²)³
// f64 calcul_de_r_etoile_div_par_r__6(f64 r_carree , f64 r_e_carree){
//     f64 tmp = calcul_de_r_etoile_div_par_r__2(r_carree, r_e_carree); //(r*)²/(r)²
//     return tmp * tmp *tmp;
// }


//Calcul du potentiel de Lennard jones
f64 potentiel_lennard_jones(f64 r_carree, f64 r_e_carree, f64 n_espsilon){
    f64 tmp_2 = calcul_de_r_etoile_div_par_r__2(r_carree, r_e_carree); 
    f64 tmp_6 = tmp_2 * tmp_2 * tmp_2;// ((r*)/(r))⁶
    f64 tmp_12 = tmp_6 * tmp_6; //((r*)²/(r)²)⁶ * ((r*)²/(r)²)⁶ = ((r*)²/(r)²)¹²
    return 4 * n_espsilon * (tmp_12 - (2 * tmp_6)); // 4 * epsilon * [(r*)²/(r)²)¹² - 2* ((r*)²/(r)²)⁶]
}


//Calcul de la force de Lennard jones, (dérivée du potentiel)
f64 force_lennard_jones(f64 r_carree, f64 r_e_carree, f64 n_espsilon){
    f64 tmp_2 = calcul_de_r_etoile_div_par_r__2(r_carree, r_e_carree); // (r*)²/(r)²
    f64 tmp_6 = tmp_2 * tmp_2 * tmp_2; // ((r*)/(r))⁶
    f64 tmp_8 = tmp_6 * tmp_2; // ((r*)²/(r)²)⁴
    f64 tmp_14 = tmp_6 * tmp_6 * tmp_2; //(r*)²/(r)²)¹⁴
    return  - 48 * n_espsilon * (tmp_14 - tmp_8); // -48 * epsilon * [(r*)²/(r)²)¹⁴ - ((r*)²/(r)²)⁸]
}

// Calcul de l'énergie total du système 
f64 calcul_potentiel_LJ(const Particules& p, unsigned int taille, f64 r_e_carree, f64 n_espsilon){
    f64 U = 0.0;
    for(unsigned int i=0; i < taille; i++){
        f64 x = p.x[i];
        f64 y = p.y[i];
        f64 z = p.z[i];
        for(unsigned int j = i + 1; j< taille ; j++){
            f64 r_carree = calcul_distance_carree(x, y, z, p.x[j], p.y[j], p.z[j]);
            U += potentiel_lennard_jones(r_carree, r_e_carree, n_espsilon);
        }
    }
    return U;
}

// Calcul de la force pour chaque composant
void calcul_force_LJ(Particules& p, unsigned int taille, f64 r_e_carree, f64 n_espsilon){
    p.fx.resize(taille);
    p.fy.resize(taille);
    p.fz.resize(taille);
    for(unsigned int i=0; i < taille; i++){
        //Initiaisation des forces
        p.fx[i] = 0.0;
        p.fy[i] = 0.0;
        p.fz[i] = 0.0;

        f64 x = p.x[i];
        f64 y = p.y[i];
        f64 z = p.z[i];

        for(unsigned int j = i + 1; j< taille ; j++){
            f64 x_j = p.x[j];
            f64 y_j = p.y[j];
            f64 z_j = p.z[j];

            f64 r_carree = calcul_distance_carree(x, y, z, x_j, y_j, z_j);

            f64 force = force_lennard_jones(r_carree, r_e_carree, n_espsilon);

            f64 dx = x - x_j;
            f64 dy = y - y_j;
            f64 dz = z - z_j;

            p.fx[i] += (force * dx);
            p.fy[i] += (force * dy);
            p.fz[i] += (force * dz);
        }
    }
}

//Calcul de l'énergie et de la force en utilisant le potentiel de Lennard Jones
void calcul_energie_force_LJ(Particules& p, unsigned int taille, f64 *U,f64 r_e_carree, f64 n_espsilon){
    for(unsigned int i=0; i < taille; i++){
        
        f64 x = p.x[i];
        f64 y = p.y[i];
        f64 z = p.z[i];

        for(unsigned int j = i + 1; j< taille ; j++){   

            f64 x_j = p.x[j];
            f64 y_j = p.y[j];
            f64 z_j = p.z[j];

            f64 r_carree = calcul_distance_carree(x, y, z, x_j, y_j, z_j);

            //Calcul de l'énergie avec le potentiel de Lennard Jones
            *U += potentiel_lennard_jones(r_carree, r_e_carree, n_espsilon);

            f64 force = force_lennard_jones(r_carree, r_e_carree, n_espsilon);

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

//Initialisation des forces à 0
void init_force(Particules& p, unsigned int taille){
    p.fx.resize(taille);
    p.fy.resize(taille);
    p.fz.resize(taille);
    for (unsigned int i = 0; i < taille; i++)
    {
        p.fx[i] = 0.0;
        p.fy[i] = 0.0;
        p.fz[i] = 0.0;
    }
}

//Verification de la nullitée des forces 
void verifier_valeur_force(const std::vector<f64>& fx, const std::vector<f64>& fy, const std::vector<f64>& fz, unsigned int taille ){
    f64 total_x = 0.0;
    f64 total_y = 0.0;
    f64 total_z = 0.0;
    for (unsigned int i = 0; i < taille; i++)
    {
        total_x += fx[i];
        total_y += fy[i];
        total_z += fz[i];
    }

    std::cout << "AFFICHAGE DE LA VALEUR TOTALE DES FORCES!!"<<std::endl;
    std::cout << " Fx : " <<total_x<< std::endl;
    std::cout << " Fy : " <<total_y<< std::endl;
    std::cout << " Fz : " <<total_z<< std::endl;
}