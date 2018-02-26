#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gmp.h"


void Euclide(mpz_t a, mpz_t p, mpz_t u, mpz_t v) {
    mpz_t r0,r1,q,r,u0,u1,v0,v1,u,v;
    
    mpz_init(r0);
    mpz_init(r1);
    mpz_init(q);
    mpz_init(r);
    mpz_init(u0);
    mpz_init(u1);
    mpz_init(v0);
    mpz_init(v1);
    mpz_init(u);
    mpz_init(v);
    
    mpz_set(r0,a);
    mpz_set(r1,p);
    mpz_set_d(u0,1);
    mpz_set_d(u1,0);
    mpz_set_d(v0,0);
    mpz_set_d(v1,1);
}