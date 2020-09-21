#ifndef SCALAR_MACROS
#define SCALAR_MACROS


#define SCALAR_APPLY_TO_FIELDS_ARGS(function, ...)        \
  function(phi, __VA_ARGS__);                           \
  function(Pi, __VA_ARGS__);                            \
  function(psi1, __VA_ARGS__);                          \
  function(psi2, __VA_ARGS__);                          \  
  function(psi3, __VA_ARGS__);                          \
  function(h11, __VA_ARGS__);                   \
  function(h12, __VA_ARGS__);                   \
  function(h13, __VA_ARGS__);                   \
  function(h22, __VA_ARGS__);                   \
  function(h23, __VA_ARGS__);                   \
  function(h33, __VA_ARGS__);                   \
  function(w11, __VA_ARGS__);                   \
  function(w12, __VA_ARGS__);                   \
  function(w13, __VA_ARGS__);                   \
  function(w22, __VA_ARGS__);                   \     
  function(w23, __VA_ARGS__);                   \
  function(w33, __VA_ARGS__);                                                  

#define SCALAR_APPLY_TO_FIELDS(function)          \
  function(phi);                                \
  function(Pi);                                 \
  function(psi1);                               \
  function(psi2);                               \
  function(psi3);                               \
  function(h11);                       \ 
  function(h12);                       \
  function(h13);                       \
  function(h22);                       \
  function(h23);                       \
  function(h33);                       \
  function(w11);                       \
  function(w12);                       \
  function(w13);                       \                          
  function(w22);                       \
  function(w23);                       \
  function(w33); 

#define SCALAR_RK_EVOLVE_PT_FIELD(field)               \
  field##_s(i,j,k) = ev_##field(&bd, &sd, dx) * dt;

#define SCALAR_RK_EVOLVE_BD_FIELD(field)               \
  field##_s(i,j,k) = ev_##field##_bd(&bd, &sd, dx, l_idx, codim) * dt;


#define SCALAR_RK_EVOLVE_PT \
  SCALAR_APPLY_TO_FIELDS(SCALAR_RK_EVOLVE_PT_FIELD)

#define SCALAR_RK_EVOLVE_BD \
  SCALAR_APPLY_TO_FIELDS(SCALAR_RK_EVOLVE_BD_FIELD)



#endif
