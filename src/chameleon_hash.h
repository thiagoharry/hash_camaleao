#ifndef _CHAMELEON_HASH_H_
#define _CHAMELEON_HASH_H_

// Please, define the following elements before including this file:
// PK           The public key
// SK           The secret key
// PAIR_OF_KEYS (should countain PK and SK)
// MSG          A encoded message
// RND          The random parameter
// DIGEST       The hash digest

// A chameleon hash scheme is defined by 3 algorithms:
struct chameleon_hash_scheme{
  PAIR_OF_KEYS (*KeyGen)(unsigned);
  DIGEST *(*Hash)(PK, MSG, unsigned, RND);
  void (*Collision)(SK, MSG, unsigned, RND, MSG, RND *);
  // Auxiliary functions:
  void (*FreePairOfKeys)(PAIR_OF_KEYS);
  void (*RandomR)(PK, RND *);
  // The following functions are optional and you should check if they
  // are null before using:
  void (*IForge)(PK, MSG, RND, MSG, RND, MSG *, RND *);
  void (*FirstPreImage)(SK, MSG, unsigned, DIGEST, RND *);
};

#endif
