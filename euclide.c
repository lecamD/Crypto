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
    while (mpz_cmp_d(r,1)!=0) {
        
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
//     printf("Nouvelle boucle\n");
//     gmp_printf("a = %Zd\n",a);
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
//     Si a est pair le résultat est expMod(res,p,g²,a/2)
    mpz_mod(mod2,a,mod2);
    if (mpz_cmp_d(mod2,0)==0) {
//         g²
        mpz_mul(gg,g,g);
//         a/2
        mpz_fdiv_q(aa,a,deux);
//         Appel récursif
        expMod(res,p,gg,aa);
        mpz_mod(res,res,p);
    } 
//     Si a > 2 est impair le résultat est g*expMod(res,p,g²,(a-1)/2)
    else if (mpz_cmp_d(a,2)>0) {
//         g²
        mpz_mul(gg,g,g);
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
    
//     Libère la mémoire
    mpz_clear(mod2);
    mpz_clear(gg);
    mpz_clear(aa);
    mpz_clear(deux);
    mpz_clear(un);
}

int main( int argc, char ** argv ) {
    
//     Test de la fonction euclide()
    
    mpz_t a,p,u,v;
    mpz_init(a);
    mpz_init(p);
    mpz_init(u);
    mpz_init(v);
    
    mpz_set_d(a,81);
    mpz_set_d(p,11);
    
    euclide(a,p,u,v);
    
    gmp_printf("%Zd*(%Zd) + %Zd*(%Zd) = 1\n",a,u,p,v);
    
//     Test de la fonction expMod()
    mpz_t res,pp,gg,aa;
    mpz_init(res);
    mpz_init(pp);
    mpz_init(gg);
    mpz_init(aa);
    
    mpz_set_d(pp,12349);
    mpz_set_d(gg,2);
    mpz_set_d(aa,34567);
    
    expMod(res,pp,gg,aa);
    gmp_printf("%Zd^%Zd = %Zd mod %Zd\n",gg,aa,res,pp);
    
    
//     Libère la mémoire
    mpz_clear(a);
    mpz_clear(p);
    mpz_clear(u);
    mpz_clear(v);
    mpz_clear(res);
    mpz_clear(pp);
    mpz_clear(gg);
    mpz_clear(aa);
}
