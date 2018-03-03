#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gmp.h"

// La fonction euclide prend en paramètre a et p et renvoie u et v de façon à ce que a*u + p*v = 1
void euclide(mpz_t a, mpz_t p, mpz_t u, mpz_t v) {
    // déclaration de variables 
    mpz_t r0,r1,q,r,u0,u1,v0,v1,tmp0,tmp1;
    // initialisation
    mpz_init(r0);mpz_init(r1);
    mpz_init(q);mpz_init(r);
    mpz_init(u0);mpz_init(u1);
    mpz_init(v0);mpz_init(v1);
    mpz_init(tmp0);mpz_init(tmp1);
    
    /* Première étape :
     * r0 prend la valeur de a
     * r1 prend la valeur de p
     * u0 prend la valeur 1
     * u1 prend la valeur 0
     * v0 prend la valeur 0
     * v1 prend la valeur 1
     */
    mpz_set(r0,a);
    mpz_set(r1,p);
    mpz_set_d(u0,1);
    mpz_set_d(u1,0);
    mpz_set_d(v0,0);
    mpz_set_d(v1,1);
    // Première itération r0 = q*r1 + r
    /* q est le quotient
     * r est le reste
     */
    
    mpz_fdiv_qr( q, r, r0, r1 );
    /* u = u0 - q*u1
     * v = v0 - q*v1
     */
    
    /* Le résultat de q*u1 est stocké dans tmp0
     * Le résultat de q*v1 est stocké dans tmp1
     */
    mpz_mul(tmp0,q,u1);
    mpz_mul(tmp1,q,v1);
    
    /* Le résultat de u0-tmp0 est stocké dans u
     * Le résultat de v0-tmp1 est stocké dans v
     */
    mpz_sub(u,u0,tmp0);
    mpz_sub(v,v0,tmp1);

    
    // mpz_cmp_d(a,b) renvoie 0 si a=b, >0 si a>b, <0 si a<b
    // On boucle tant que le reste est différent de 0
    
    /* A chaque boucle :
     * r0 prend la valeur de r1
     * r1 prend la valeur de r
     * u0 prend la valeur de u1
     * u1 prend la valeur de u
     * vO prend la valeur de v1
     * v1 prend la valeur de v
     * q et r sont recalculés à partir des nouvelles valeurs de r0 et r1
     * u et v sont recalculés à partir des nouvelles valeurs de de q,u0,v0,u1,v1
     */

    while (mpz_cmp_d(r,1)>0) {
        mpz_set(r0,r1);
        mpz_set(r1,r);
        mpz_set(u0,u1);
        mpz_set(u1,u);
        mpz_set(v0,v1);
        mpz_set(v1,v);
    
        mpz_fdiv_qr( q, r, r0, r1 );
    
        mpz_mul(tmp0,q,u1);
        mpz_mul(tmp1,q,v1);
        

        mpz_sub(u,u0,tmp0);
        mpz_sub(v,v0,tmp1);
    }

    
//     gmp_printf("q : %Zd, r : %Zd\n",q,r);

// Libère la mémoire
    mpz_clear(r0);mpz_clear(r1);
    mpz_clear(q);mpz_clear(r);
    mpz_clear(u0);mpz_clear(u1);
    mpz_clear(v0);mpz_clear(v1);
    mpz_clear(tmp0);mpz_clear(tmp1);
}


// A = g^a mod p
void expMod(mpz_t res,mpz_t p,mpz_t g,mpz_t a) {
//         On applique le modulo après chaque calcul
    mpz_t mod2,gg,aa,un,deux;
//     initialisation des variables 
    mpz_init(gg);
    mpz_init(aa);
    mpz_init(deux);
    mpz_init(un);
    mpz_init(mod2);
    mpz_set_d(mod2,2);
    mpz_set_d(deux,2);    
    mpz_set_d(un,1);
    //     Si a = 1 le résultat est g
    if (mpz_cmp_d(a,1)==0) {
        mpz_mod(gg,g,p);
        mpz_set(res,gg);
    }
    else {
//     Si a est pair le résultat est expMod(res,p,g²,a/2)
        mpz_mod(mod2,a,mod2);
        if (mpz_cmp_d(mod2,0)==0) {
//            g²
            mpz_mul(gg,g,g);
//         a/2
            mpz_fdiv_q(aa,a,deux);
//         Appel récursif
            expMod(res,p,gg,aa);
            mpz_mod(res,res,p);
        } 
//     Si a > 2 est impair le résultat est g*expMod(res,p,g²,(a-1)/2)
        else { 
            if (mpz_cmp_d(a,2)>0) {
//         g²
                mpz_mul(gg,g,g);
                mpz_mod(gg,gg,p);
//         a-1
                mpz_sub(aa,a,un);
//         (a-1)/2
                mpz_fdiv_q(aa,aa,deux);
//         Appel récursif
                expMod(res,p,gg,aa);
                mpz_mod(res,res,p);
//         g*expMod(...)
                mpz_mul(res,g,res);
                mpz_mod(res,res,p);
            }
        }
    
    }
//     Libère la mémoire
    mpz_clear(mod2);
    mpz_clear(gg);
    mpz_clear(aa);
    mpz_clear(deux);
    mpz_clear(un);
}


void keyGen(mpz_t p, mpz_t g, mpz_t x, mpz_t X, gmp_randstate_t state) {
    mpz_t un,k;
    mpz_init(un);
    mpz_init(k);
    mpz_set_d(un,1);
// //     initialisation pour l'aléatore
//     gmp_randstate_t state;
//     gmp_randinit_default (state);
//     gmp_randseed_ui (state, (unsigned) time(NULL));
//     Tire au hasard un x entre 0 et p-2
    mpz_sub(k,p,un);
    mpz_urandomm(x,state,k);
//     Calcul X = g^x mod p
    expMod(X, p, g, x);
}

void encrypt(mpz_t C, mpz_t B, mpz_t p, mpz_t g, mpz_t X, mpz_t m,mpz_t r,gmp_randstate_t state) {
    
    mpz_t un,k,y;
    mpz_init(un);
    mpz_init(k);
    mpz_init(y);
    mpz_set_d(un,1);
    
//     initialisation pour l'aléatore
//     gmp_randstate_t state;
//     gmp_randinit_default (state);
//     gmp_randseed_ui (state, (unsigned) time(NULL));
    
//     Tire au hasard un nombre r entre 0 et p-2
    mpz_sub(k,p,un);
    mpz_urandomm(r,state,k);
    mpz_urandomm(r,state,k);
//     gmp_printf("    r = %Zd\n",r);
//     Calcul y = X^r mod p
    expMod(y, p, X, r);
//     gmp_printf("    y = %Zd\n",y);
    
//     C = m * y mod p
    mpz_mul(C,m,y);
    mpz_mod(C,C,p);
    
//     B = g^r mod p
    expMod(B,p,g,r);
}

void decrypt(mpz_t C, mpz_t B,mpz_t x,mpz_t m, mpz_t p) {
    mpz_t D,u,v,t;
    mpz_init(D);
    mpz_init(u);
    mpz_init(v);
//     D = B^x mod p
    expMod(D,p,B,x);
//     gmp_printf("    D = %Zd\n",D);
//     u = (D)^-1
    euclide(D, p, u, v);
//     gmp_printf("    u = %Zd\n   v = %Zd\n",u,v);
//     C * (D)^-1 mod p
    mpz_mul(m,C,u);
    mpz_mod(m,m,p);
    
    mpz_clear(D);mpz_clear(u);mpz_clear(v);
}



int main( int argc, char ** argv ) {
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui (state, (unsigned) time(NULL));
    FILE * fp;
    fp = fopen ("test.txt", "w+");
       
//     initialisation de p = 2^1024 − 2^960 − 1 + 2^64 ∗ ([2^894 π] + 129093) et g = 2
    mpz_t k,b,c,d,p,q,un,g,deux,truc,cent;
    mpz_init(k);mpz_init(b);mpz_init(c);mpz_init(d);mpz_init(q);mpz_init(un);mpz_init(deux);mpz_init(truc);mpz_init(cent);
    mpz_init(g);mpz_init(p);
    mpz_set_d(un,1);
    mpz_set_d(deux,2);
    mpz_set_d(truc,129093); 
    mpz_set_d(cent,100);
    
//     2^1024
    mpz_pow_ui(k,deux,1024);
//     2^960
    mpz_pow_ui(b,deux,960);
//     2^64
    mpz_pow_ui(c,deux,64);
//     2^894
    mpz_pow_ui(d,deux,894);
//     2^1024 − 2^960
    mpz_sub(p,k,b);
//    (2^1024 − 2^960) - 1
    mpz_sub(p,p,un);
//     (2^1024 − 2^960 − 1) + 2^64
    mpz_add(p,p,c);
//     2^894 * 314
    mpz_mul_ui(q,d,314);
//     (2^894 * 314) / 100
    mpz_fdiv_q(q,q,cent);
//     (2^894 * π) + 129093
    mpz_add(q,q,truc);
//     (2^1024 − 2^960 − 1 + 2^64) ∗ ([2^894 π] + 129093)
    mpz_mul(p,p,q);
    
    
    mpz_clear(k);mpz_clear(b);mpz_clear(c);mpz_clear(d);mpz_clear(q);mpz_clear(truc);mpz_clear(un);mpz_clear(deux);mpz_clear(cent);
//     
// //    p =  
// // 74552348966919475999252352065578124041644418824828005444807371572971918118807129822843933514921654543589158889382724163333540486848544797300790583472480953829458644890488139850347078442895190178137605518160621347885987509079573008588779202218426541234334203658740434002771140502224523168013812274864773843243435524417268618512771558980788120310982900822695893919136238191653865219655852574475844837862711973004021255940937431491801716342376266668110205626780626493974387752884772478512132126093497864917801347061932974000514178185105720025718058922503179907990813359665487929030
//     
//         
    mpz_set_d(g,2);
    
//     p pour les tests
//     mpz_set_d(p,1344567754356789876);
    mpz_nextprime(p,p);
    
//     initialisation pour l'aléatore
    
//     Question 3 
//     Test sur 5 valeur de a différentes pour Euclide()
    mpz_t a,u,v;
    mpz_init(a);mpz_init(u);mpz_init(v);
    int i = 0;
    gmp_printf("Partie tests sur Euclide en cours\n");
    gmp_fprintf(fp,"Euclide\n\n");
    for (i=0;i<5;i++) {
        mpz_urandomm(a,state,p);
        euclide(a,p,u,v);
        gmp_fprintf(fp,"essaie %d\n    a = %Zd\n    u = %Zd\n    v = %Zd\n\n",i,a,u,v);
    }
    mpz_clear(u);mpz_clear(v);
    
//     Question 4
//     Test sur 5 valeur de a différentes pour expMod()
    gmp_printf("Partie tests sur Exponentiation modulaire en cours\n");
    mpz_t A;
    mpz_init(A);
    gmp_fprintf(fp,"Exponentiation modulaire\n\n");
    for (i=0;i<5;i++) {
        mpz_urandomm(a,state,p);
        expMod(A,p,g,a);
        gmp_fprintf(fp,"essaie %d\n    a = %Zd\n    A = %Zd\n\n",i,a,A);
    }
    mpz_clear(A);
    
//     Question 5
//     Test sur 5 valeur de m pour le chiffrement + dechriffrement
    gmp_printf("Partie tests sur chiffrement + déchiffrement en cours\n");
    mpz_t x,X,m,C,B,r;
    mpz_init(x);mpz_init(X);mpz_init(m);mpz_init(C);mpz_init(B);mpz_init(r);
    gmp_fprintf(fp,"chiffrement + déchiffrement\n");
    for (i=0;i<5;i++) {
        mpz_urandomb(m,state,50);
        gmp_fprintf(fp,"essaie %d\n    m = %Zd\n\n",i,m);
        
        keyGen(p, g, x, X, state) ;
//     gmp_printf("Clé secrète :\n    x = %Zd\nClé publique :\n   p = %Zd\n   g = %Zd\n   X = %Zd\n",x,p,g,X); 
        encrypt(C, B, p, g, X, m, r,state) ;
        gmp_fprintf(fp,"    r = %Zd\n\n",r);
//     gmp_printf("   m = %Zd\n   C = %Zd\n   B = %Zd\n",m,C,B);
        decrypt(C, B, x, m, p) ;
//     gmp_printf("    m = %Zd\n",m);
        gmp_fprintf(fp,"    m déchiffré = %Zd\n\n",m);
    }
    
    
// Question 6
// Test sur 5 valeurs différentes pour les deux messages à chiffrés
    gmp_printf("Partie tests sur Propriété homomorphique en cours\n");
    gmp_fprintf(fp,"propriété homomorphique du chiffrement El Gamal\n");
    mpz_t mtmp,m1,m2;
    
    mpz_init(mtmp);
    mpz_init(m1);
    mpz_init(m2);
    
    mpz_t C1,C2,B1,B2;
    mpz_init(C1);mpz_init(B1);mpz_init(C2);mpz_init(B2);
    
    for (i=0;i<5;i++) {
        mpz_urandomb(m1,state,50);
        mpz_urandomb(m2,state,50);
        
        mpz_mod(m1,m1,p);
        mpz_mod(m2,m2,p);
            
    //     Génération de la clé
        keyGen(p, g, x, X, state) ;
    //     chiffrement des deux messages
        encrypt(C1, B1, p, g, X, m1, r, state);
        encrypt(C2, B2, p, g, X, m2, r, state);
        
    //     multiplication des chiffrés
        mpz_mul(C,C1,C2);
        mpz_mod(C,C,p);
        mpz_mul(B,B1,B2);
        mpz_mod(B,B,p);
    //     déchiffrement du couple (C,B)
        decrypt(C, B, x, m, p);
        
    //     Vérification m = m1*m2 ?
        mpz_mul(mtmp,m1,m2);
        mpz_mod(mtmp,mtmp,p);
        
        gmp_fprintf(fp,"Essaie %d\n    m1*m2 = %Zd\n    m déchiffré = %Zd\n\n",i,mtmp,m);
        
    //     Si m = mtmp ( Vérif de la propriété homomorphique d'ElGamal
        if (mpz_cmp(m,mtmp)==0)
            gmp_fprintf(fp,"Propriété montrée pour\n    m1 = %Zd\n    m2 = %Zd\n\n",m1,m2);
    }
    printf("Création d'un fichier test.txt avec les résultats des tests\n");
    fclose(fp);
    mpz_clear(r);mpz_clear(m);mpz_clear(mtmp);mpz_clear(m1);mpz_clear(m2);mpz_clear(C);mpz_clear(B);mpz_clear(C1);mpz_clear(B1);mpz_clear(C2);mpz_clear(B2);mpz_clear(x);mpz_clear(X);
}
