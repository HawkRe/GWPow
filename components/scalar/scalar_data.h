#ifndef COSMO_SCALAR_DATA
#define COSMO_SCALAR_DATA

namespace cosmo
{
typedef struct {

  // field values
  real_t phi, Pi;
  real_t h11, h12, h13, h22, h23, h33;
  real_t w11, w12, w13, w22, w23, w33;
  // derivatives of fields
  real_t d1phi, d2phi, d3phi;
  real_t d1Pi, d2Pi, d3Pi;
  real_t d1psi1, d2psi1, d3psi1;
  real_t d1psi2, d2psi2, d3psi2;
  real_t d1psi3, d2psi3, d3psi3;
} ScalarData;

}
#endif
