 b=BarrierOption(100,100,260,200,0.01,1,0.3,0.2,0.5,3,100,500);
 price = b.Finite_diff_price(100,100);
 [delta1,delta2] = b.Finite_diff_delta(100,100);
 [gamma1,gamma2] = b.Finite_diff_gamma(100,100);
 [vega1,vega2] = b.Finite_diff_vega();
 theta = b.Finite_diff_theta();