#ifndef _MOD_MATH_H_
#define _MOD_MATH_H_

#include <gmp.h>
#include <stdbool.h>

int root_mod(mpz_t, const mpz_t, const mpz_t);
void mod_random_prime(mpz_t, unsigned);
bool is_quadratic_residue(mpz_t, mpz_t);
int root_mod_pq(mpz_t, const mpz_t, const mpz_t, const mpz_t, const mpz_t);

#endif
