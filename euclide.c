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
    
    gmp_printf("q : %Zd, r : %Zd\n",q,r);
}

// A = g^a mod p
void expMod(mpz_t p,mpz_t g,mpz_t a) {

    
    
    
    
}


int main( int argc, char ** argv ) {
    
    mpz_t a,p,u,v;
    mpz_init(a);
    mpz_init(p);
    mpz_init(u);
    mpz_init(v);
    
    mpz_set_d(a,81);
    mpz_set_d(p,11);
    
    euclide(a,p,u,v);
    
    gmp_printf("%Zd*%Zd + %Zd*%Zd = 1\n",a,u,p,v);


}
