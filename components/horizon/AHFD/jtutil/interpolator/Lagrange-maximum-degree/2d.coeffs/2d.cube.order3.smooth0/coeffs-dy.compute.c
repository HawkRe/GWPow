fp t240;
fp t227;
fp t219;
fp t226;
fp t222;
fp t211;
fp t239;
fp t201;
fp t238;
fp t223;
fp t213;
fp t237;
fp t200;
fp t236;
fp t203;
fp t235;
fp t202;
fp t234;
fp t224;
fp t216;
fp t221;
fp t233;
fp t225;
fp t217;
fp t232;
fp t218;
fp t231;
fp t230;
fp t220;
fp t229;
fp t228;
fp t215;
fp t210;
fp t207;
fp t206;
fp t205;
fp t204;
      t240 = x*y;
      t227 = y*y;
      t219 = RATIONAL(-3.0,8.0)*t227;
      t226 = x*x;
      t222 = RATIONAL(1.0,40.0);
      t211 = t222*t226;
      t239 = t219+t211;
      t201 = RATIONAL(1.0,20.0)*t240;
      t238 = t201+RATIONAL(9.0,40.0)*y;
      t223 = RATIONAL(-1.0,40.0);
      t213 = t223*t226;
      t237 = t219+t213;
      t200 = RATIONAL(3.0,20.0)*t240;
      t236 = t200+RATIONAL(7.0,40.0)*y;
      t203 = RATIONAL(-3.0,20.0)*t240;
      t235 = t203+RATIONAL(13.0,40.0)*y;
      t202 = RATIONAL(-1.0,20.0)*t240;
      t234 = t202+RATIONAL(11.0,40.0)*y;
      t224 = RATIONAL(3.0,40.0);
      t216 = t224*t226;
      t221 = RATIONAL(-1.0,8.0)*t227;
      t233 = t216+t221;
      t225 = RATIONAL(-3.0,40.0);
      t217 = t225*t226;
      t232 = t217+t221;
      t218 = RATIONAL(1.0,8.0)*t227;
      t231 = t218+t217;
      t230 = t216+t218;
      t220 = RATIONAL(3.0,8.0)*t227;
      t229 = t220+t213;
      t228 = t220+t211;
      t215 = RATIONAL(-1.0,50.0)*x;
      t210 = RATIONAL(2.0,25.0)*x;
      t207 = RATIONAL(-9.0,100.0)*x;
      t206 = RATIONAL(-1.0,100.0)*x;
      t205 = RATIONAL(7.0,100.0)*x;
      t204 = RATIONAL(-13.0,100.0)*x;
      coeffs_dy->coeff_m1_m1 = RATIONAL(6.0,25.0)*x+RATIONAL(-109.0,1200.0)+
t232+t235;
      coeffs_dy->coeff_0_m1 = t215+RATIONAL(-223.0,1200.0)+t233+t234;
      coeffs_dy->coeff_p1_m1 = RATIONAL(-157.0,1200.0)+t204+t233+t238;
      coeffs_dy->coeff_p2_m1 = t207+RATIONAL(89.0,1200.0)+t232+t236;
      coeffs_dy->coeff_m1_0 = RATIONAL(-23.0,40.0)*y+RATIONAL(-31.0,400.0)+t200
+t215+t229;
      coeffs_dy->coeff_0_0 = RATIONAL(-21.0,40.0)*y+RATIONAL(-57.0,400.0)+
RATIONAL(-1.0,25.0)*x+t201+t228;
      coeffs_dy->coeff_p1_0 = t202+RATIONAL(-19.0,40.0)*y+RATIONAL(-63.0,400.0)
+t206+t228;
      coeffs_dy->coeff_p2_0 = t203+RATIONAL(-49.0,400.0)+RATIONAL(-17.0,40.0)*y
+t205+t229;
      coeffs_dy->coeff_m1_p1 = t204+RATIONAL(111.0,400.0)+t236+t239;
      coeffs_dy->coeff_0_p1 = t206+RATIONAL(117.0,400.0)+t237+t238;
      coeffs_dy->coeff_p1_p1 = RATIONAL(3.0,50.0)*x+RATIONAL(103.0,400.0)+t234+
t237;
      coeffs_dy->coeff_p2_p1 = t210+RATIONAL(69.0,400.0)+t235+t239;
      coeffs_dy->coeff_m1_p2 = t207+t203+RATIONAL(-131.0,1200.0)+t224*y+t230;
      coeffs_dy->coeff_0_p2 = t205+t202+t222*y+RATIONAL(43.0,1200.0)+t231;
      coeffs_dy->coeff_p1_p2 = t201+RATIONAL(37.0,1200.0)+t223*y+t210+t231;
      coeffs_dy->coeff_p2_p2 = t225*y+t200+RATIONAL(-3.0,50.0)*x+RATIONAL(
-149.0,1200.0)+t230;
