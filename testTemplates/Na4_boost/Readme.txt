
Note: The original data are errorneous, i.e., the lines especially of the two largest excitations (one Z and one Y) are no proper sine cardinals, which is the difficulty here. 
	A possible reason may be a too large excitation energy.

eval.yaml.0 - eval.yaml.1: Add excitations
eval.yaml.1.res:
	-> S0(Y) and S1(Z) have low signifFit. The reason ist that S0 is a small Y-excitation next to a large Z-excitation S1. They are slightly mixed due to the bad line-shape of th elarge Z-excitation which gives S0 a large Z-component. The fit is not good due to the line-shape error, however, both lines are significant.
	-> S3(X) has low signifRng, i.e., the fitted energy gets very close the the energy-range boundary.
	-> S6(Y) and S7(Y) have a low signifFit. This is the large Y-Line with the large line-shape error. There should only one excitation

Refine the result as far as possible:

eval.yaml.2:
	-> Remove S7
	-> Reset the energy fit range
	-> refit remaining excitations

eval.yaml.2.res:
	-> S0 and S1: as before; cannot get better with error-prone input data
	-> S3(Z): Small signifExc and Very close to S2(Z) - Probably erroneous line
	-> S6(Y): Rather small signifFit (~0.6). This is the large but error-prone Y-Line

eval.yaml.3:
	-> Remove S3 and re-fit

eval.yaml.3.res:
	-> Very good
	-> Now, S4 and S5 (large Y line) show not perfect signifFit. Since S4 also points in Y direction, it should be an artifact from the bad S5 line shape

eval.yaml.4:
	-> Remove S4 and re-fit

eval.yaml.4.res:
	-> Final result, all significances are either ok or the reason for a low significance is known
