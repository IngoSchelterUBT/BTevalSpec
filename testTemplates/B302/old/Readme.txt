
 - eval.yaml.0:     Add initial excitations
 - eval.yaml.0.res: Some have low significances. Don't get discouraged: This means that there are still excitations missing
 - eval.yaml.1:     Add one batch of additional excitations
 - eval.yaml.1.res: Suddenly, almost all excitations have high significance values. There is one with an angle close to 90° which shows a large, negative signifErr and some excitations are moved to the edge of their energy-range fit interval, e.g., S40.
 - eval.yaml.2:     Reset the energy-fit constraint to a +/- pi/T interval of the current energies and re-fit
 - eval.yaml.2.res:
 	-> S39's energy still moves to the edge of the energy interval
 	-> S35 hat too large strength and negative signifErr
	-> S9 and S13 have negative signifExc
 	-> S0 (Qy) is removed somehow (0.0 strength) and the fit error is very large
 	-> roll back to eval.yaml.1.res, fit all excitations apart from those that are close to the boundary
 - eval.yaml.2b:
 	-> Fix all excitations apart from
		o S9
		o S12
		o S40
		o S50
	-> Force set energy-fit constraint interval of the 0.1:0.35 (which is the total fit range)
	-> re-fit
 - eval.yaml.2b.res:
 	-> S9 and S35 have very low signifAng, signif Err
	-> S50 has very low signifRng (Keep this excitation as unsignificant excitation at the high end of the fit interval)
 - eval.yaml.3:
 	-> Remove S9 and S35
	-> Refit
 - eval.yaml.3.res:
  -> Almost perfect; overall fit error ~0.025; Almost all significances high, especially signifErr
	-> Last excitation significance is low (since it's the boundary of the interval, this is fine)
	-> By looking at the difference between fit and ft, there may be additional excitations at
		o 0.1681 Ry
		o 0.2972 Ry
		o 0.3154 Ry
		Some of them are, however, smaller than the error of some larger lines
		Before adding additional lines, I try to refine the fit setting firstsingle: true
 - eval.yaml.4:
 	-> re-fit with firstsingle: true
 - eval.yaml.5:
 	-> re-fit with firstsingle: false
 - eval.yaml.5.res:
	-> Almost perfect; S37 has >85° (signifAng~0.2) and signifRng=~0
 - eval.yaml.6:
	-> Remove S37 and re-fit
 - eval.yaml.6.res:
 	-> Final fit (error ~ 0.025)
