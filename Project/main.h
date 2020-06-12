#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <time.h>

void verification_primalite(int p);
void Puissance(mpz_t resultat,mpz_t n, mpz_t p);
void Decomposition (mpz_t t,mpz_t s,mpz_t n);
int TestFermat(mpz_t n,int k);
int TestMillerRabin(mpz_t n,int k);
void squareAndMultiply( mpz_t res, mpz_t a,mpz_t n,mpz_t H);
