fp t91;
fp t95;
fp t104;
fp t94;
fp t86;
fp t103;
fp t93;
fp t88;
fp t102;
fp t89;
fp t101;
fp t92;
fp t100;
fp t99;
fp t87;
fp t98;
fp t97;
fp t96;
fp t90;
      t91 = RATIONAL(1.0,9.0)*x;
      t95 = RATIONAL(1.0,18.0);
      t104 = t91+t95;
      t94 = RATIONAL(1.0,12.0);
      t86 = t94*y;
      t103 = t91+t86;
      t93 = RATIONAL(-1.0,12.0);
      t88 = t93*y;
      t102 = t91+t88;
      t89 = t94*z;
      t101 = t91+t89;
      t92 = RATIONAL(-1.0,18.0);
      t100 = t92+t91;
      t99 = t86+t101;
      t87 = t93*z;
      t98 = t87+t103;
      t97 = t87+t102;
      t96 = t88+t101;
      t90 = RATIONAL(-2.0,9.0)*x;
      coeffs_dx->coeff_m1_m1_m1 = t92+t99;
      coeffs_dx->coeff_0_m1_m1 = t90;
      coeffs_dx->coeff_p1_m1_m1 = t95+t97;
      coeffs_dx->coeff_m1_0_m1 = t89+t100;
      coeffs_dx->coeff_0_0_m1 = t90;
      coeffs_dx->coeff_p1_0_m1 = t87+t104;
      coeffs_dx->coeff_m1_p1_m1 = t92+t96;
      coeffs_dx->coeff_0_p1_m1 = t90;
      coeffs_dx->coeff_p1_p1_m1 = t95+t98;
      coeffs_dx->coeff_m1_m1_0 = t86+t100;
      coeffs_dx->coeff_0_m1_0 = t90;
      coeffs_dx->coeff_p1_m1_0 = t95+t102;
      coeffs_dx->coeff_m1_0_0 = t100;
      coeffs_dx->coeff_0_0_0 = t90;
      coeffs_dx->coeff_p1_0_0 = t104;
      coeffs_dx->coeff_m1_p1_0 = t88+t100;
      coeffs_dx->coeff_0_p1_0 = t90;
      coeffs_dx->coeff_p1_p1_0 = t95+t103;
      coeffs_dx->coeff_m1_m1_p1 = t92+t98;
      coeffs_dx->coeff_0_m1_p1 = t90;
      coeffs_dx->coeff_p1_m1_p1 = t95+t96;
      coeffs_dx->coeff_m1_0_p1 = t87+t100;
      coeffs_dx->coeff_0_0_p1 = t90;
      coeffs_dx->coeff_p1_0_p1 = t95+t101;
      coeffs_dx->coeff_m1_p1_p1 = t92+t97;
      coeffs_dx->coeff_0_p1_p1 = t90;
      coeffs_dx->coeff_p1_p1_p1 = t95+t99;
