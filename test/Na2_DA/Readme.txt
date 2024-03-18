eval.yaml.0:
 - Initial fit 16 excitations
 
eval.yaml.0.res:
 - Overall fit-error is very small with about 0.0012 (also the fit represents the data perfectly)
 - Some excitations have low significance:
 	-> S0: at the lower boundary of the fit interval

eval.yaml.1:
 - Remove S0
 - Remove S3 maybe
 - Remove S8 (extremly bad significances)
 - Reset erange
 - refit

 eval.yaml.1.res:
 - Overall fit-error is now "much" larger with 0.013
 - Low signifExc of S4,S8,S3, and S5 (even negative)
 - Low signifErr of S5 and S12
 - Low signifRng of S2
 - The reason is probably, that the above excitations are really small compared to the dominant excitations.
 	However, there really seem to be excitations by looking very closely at the FT

eval.yaml.2:
 - Remove S4, S8, S3, S5, and S12
 - Reset erange
 - refit
