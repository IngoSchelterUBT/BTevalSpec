DESCRIPTION:                        # Description of evaluation 
DIPOLE:                             # List of dipole moment files; different calculations are only allowed to differ in the boost direction or laser polarization
- [dipole_global1_01.dat]
OPT:
  FT:
    calc: false                     # Turn fourier transformation on/off
    minpw: 17                       # Min. length of fourier transform vector (2^n). Can increase sampling rate.
    smooth: 0                       # Artificial decay rate
    window: 0                       # Kaiser-Bessel windowing parameter (>0)
  Pade:
    calc: false                     # Turn Pade Approximation on/off
    wmax: 0.6                       # Maximum energy for Pade Approximation in Ry
    dw: 1.0e-05                     # Step of Pade Approximation
    smooth: 0.004455749593509885    # Smooth for Pade Approximation
    thin: 0                         # Only keep every 2^n data point
  Fit:
    calc: true                      # Turn fit on/off
    skipfirst: true                 # Skip first fit of existing excitations
    guess: false                    # Turn guess for fit via Pade Approximation on/off
    plot_result: false              # If True: The fit results are plotted without fitting again
    gnuplot_spectrum: false         # If True: A gnuplot script for plotting the resulting spectrum is created.
    dat_spectrum: false             # If True: A .dat file for the spectrum is created
    guess_thres: 0.1                # Relative height of line in Pade Approximation compared to highest line 
                                    #   which should be identified as a line for fitting (only relevent, if 
                                    #   fit_guess == True).
    relerr_crit: 0.01               # Criterium for relative error between fit and data (only relevant, if 
                                    #   fit_guess == True).
    max_excit: 24                   # Maximum numer of iterations used to reach relative error between fit and 
                                    #   data (only relevant, if fit_guess == True)
    relspacing_lines: 0.01          # Threshold for relative error between two lines which should be identified 
                                    #   as one in fit of trace (only relevant, if number of dipole files == 3).
    range:                          # Range of spectrum which should be fitted in Ry
    - 0.10
    - 0.40
    significances: true
EXT: {}
SPEC:
- name: S0
  fix: false
  energy: 0.12844309500292983
  erange:
  - 0.12590098276091613
  - 0.1319803021449414
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - -0.0027844174744710987
    - 0.1334526746986033
    - 1.0756340883878452
  dipole:
  - -0.0027844174744710987
  - 0.1334526746986033
  - 1.0756340883878452
  strength: 0.025149287760500395
  strengths:
  - 0.025149287760500395
  signif: 1.0002989781340892
- name: S1
  fix: false
  energy: 0.1318404690282606
  erange:
  - 0.12866034030798737
  - 0.13473965969201265
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - -0.0008383546378367826
    - 0.03452614676710996
    - 6.9928827771917685
  dipole:
  - -0.0008383546378367826
  - 0.03452614676710996
  - 6.9928827771917685
  strength: 1.0745350304280348
  strengths:
  - 1.0745350304280348
  signif: 0.986478406644503
- name: S2
  fix: false
  energy: 0.1545365372451253
  erange:
  - 0.15272034030798737
  - 0.15879965969201265
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - -0.021714822155749183
    - -0.03948811736688099
    - 1.8678862859394187
  dipole:
  - -0.021714822155749183
  - -0.03948811736688099
  - 1.8678862859394187
  strength: 0.08991528183940402
  strengths:
  - 0.08991528183940402
  signif: 0.4392474422034526
- name: S3
  fix: false
  energy: 0.15602034030800316
  erange:
  - 0.15602034030798737
  - 0.16209965969201265
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.08183029952592068
    - 0.034477663571591724
    - -0.9215360248554766
  dipole:
  - 0.08183029952592068
  - 0.034477663571591724
  - -0.9215360248554766
  strength: 0.02228785801906086
  strengths:
  - 0.02228785801906086
  signif: 0.038452700501163825
- name: S4
  fix: false
  energy: 0.16475999696302554
  erange:
  - 0.1618466654591247
  - 0.16792598484315
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 2.083163297327073
    - 0.04255778999499489
    - -0.006693246690311343
  dipole:
  - 2.083163297327073
  - 0.04255778999499489
  - -0.006693246690311343
  strength: 0.11921553622120462
  strengths:
  - 0.11921553622120462
  signif: 1.0
- name: S5
  fix: false
  energy: 0.1934061186692168
  erange:
  - 0.19024034030798737
  - 0.19631965969201265
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.000108771179631521
    - 5.316937250664893
    - -0.005341871699232213
  dipole:
  - 0.000108771179631521
  - 5.316937250664893
  - -0.005341871699232213
  strength: 0.9112603361724267
  strengths:
  - 0.9112603361724267
  signif: 1.0
- name: S6
  fix: false
  energy: 0.20499291450588975
  erange:
  - 0.2021985931332427
  - 0.208277912517268
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.013414412563233382
    - 0.09651423044057322
    - 0.7889009310895325
  dipole:
  - 0.013414412563233382
  - 0.09651423044057322
  - 0.7889009310895325
  strength: 0.02158779091943742
  strengths:
  - 0.02158779091943742
  signif: 1.0
- name: S7
  fix: false
  energy: 0.21822672218456315
  erange:
  - 0.2148503403079874
  - 0.22092965969201267
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 4.098654948475863
    - -0.026676295608382253
    - -0.015583599369787694
  dipole:
  - 4.098654948475863
  - -0.026676295608382253
  - -0.015583599369787694
  strength: 0.6110321619323927
  strengths:
  - 0.6110321619323927
  signif: 1.0
- name: S8
  fix: false
  energy: 0.22431063317927724
  erange:
  - 0.2214469264490576
  - 0.2275262458330829
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.02825974109456004
    - -0.07435594591871546
    - 1.0731591643842058
  dipole:
  - 0.02825974109456004
  - -0.07435594591871546
  - 1.0731591643842058
  strength: 0.04329187787311326
  strengths:
  - 0.04329187787311326
  signif: 1.0
- name: S9
  fix: false
  energy: 0.24264198914891516
  erange:
  - 0.2411903403079874
  - 0.24726965969201267
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.017600922121602136
    - 2.2083219605161504
    - -0.020815982822565197
  dipole:
  - 0.017600922121602136
  - 2.2083219605161504
  - -0.020815982822565197
  strength: 0.19724484489540736
  strengths:
  - 0.19724484489540736
  signif: 1.0
- name: S10
  fix: false
  energy: 0.2550580564178699
  erange:
  - 0.2526703403079874
  - 0.2587496596920127
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 2.7044754714390096
    - -0.015690414252963067
    - 0.015228244655607485
  dipole:
  - 2.7044754714390096
  - -0.015690414252963067
  - 0.015228244655607485
  strength: 0.31094406794076984
  strengths:
  - 0.31094406794076984
  signif: 1.0
- name: S11
  fix: false
  energy: 0.2650688342955755
  erange:
  - 0.2617988541231756
  - 0.2678781735072009
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.004261159308202119
    - 0.16233931655560013
    - 0.27691770977513963
  dipole:
  - 0.004261159308202119
  - 0.16233931655560013
  - 0.27691770977513963
  strength: 0.004552805915082597
  strengths:
  - 0.004552805915082597
  signif: 1.0
- name: S12
  fix: false
  energy: 0.3012915244312743
  erange:
  - 0.29844025971231725
  - 0.30451957909634253
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 1.0655836379235637
    - -0.03029865969090581
    - -0.007044222782484175
  dipole:
  - 1.0655836379235637
  - -0.03029865969090581
  - -0.007044222782484175
  strength: 0.057066428465306594
  strengths:
  - 0.057066428465306594
  signif: 0.23046228522324866
- name: S13
  fix: false
  energy: 0.3028455296289227
  erange:
  - 0.30122315127604954
  - 0.3073024706600748
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.11536294426458153
    - 0.6617200822880821
    - 0.010596419467891032
  dipole:
  - 0.11536294426458153
  - 0.6617200822880821
  - 0.010596419467891032
  strength: 0.022778746594855823
  strengths:
  - 0.022778746594855823
  signif: 0.5316358467078013
- name: S14
  fix: false
  energy: 0.31166602987782877
  erange:
  - 0.30864419544600225
  - 0.31472351483002753
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - -0.06721420165572545
    - -0.07793311381803075
    - 0.5665854084643507
  dipole:
  - -0.06721420165572545
  - -0.07793311381803075
  - 0.5665854084643507
  strength: 0.017225279908156598
  strengths:
  - 0.017225279908156598
  signif: 0.4873069805332686
- name: S15
  fix: false
  energy: 0.31440800506340344
  erange:
  - 0.30864419544600225
  - 0.31472351483002753
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.5235098311875785
    - 0.04453170752407276
    - 0.05602612827028034
  dipole:
  - 0.5235098311875785
  - 0.04453170752407276
  - 0.05602612827028034
  strength: 0.014629642554460184
  strengths:
  - 0.014629642554460184
  signif: 0.128300873269343
- name: S16
  fix: false
  energy: 0.31613643667800745
  erange:
  - 0.3130504404219117
  - 0.319129759805937
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.023740975318627275
    - -0.06344813350500618
    - -0.21510184421914352
  dipole:
  - 0.023740975318627275
  - -0.06344813350500618
  - -0.21510184421914352
  strength: 0.0026796832399859746
  strengths:
  - 0.0026796832399859746
  signif: 0.8303339886842622
- name: S17
  fix: false
  energy: 0.318196342300208
  erange:
  - 0.3149057014643999
  - 0.3209850208484252
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.006650084425273872
    - 0.7604645545346816
    - 0.03526939236745864
  dipole:
  - 0.006650084425273872
  - 0.7604645545346816
  - 0.03526939236745864
  strength: 0.030737474581286423
  strengths:
  - 0.030737474581286423
  signif: 0.25240362518248394
- name: S18
  fix: false
  energy: 0.3335484975695162
  erange:
  - 0.3304435126952384
  - 0.3365228320792637
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.954622129089065
    - -0.00916318078690435
    - -0.016567788152853567
  dipole:
  - 0.954622129089065
  - -0.00916318078690435
  - -0.016567788152853567
  strength: 0.05068057421774235
  strengths:
  - 0.05068057421774235
  signif: 0.9436491989625484
- name: S19
  fix: false
  energy: 0.3368662018752612
  erange:
  - 0.33392212714990377
  - 0.34000144653392905
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.08436510343047166
    - 0.5673435479658224
    - 0.017380198085367646
  dipole:
  - 0.08436510343047166
  - 0.5673435479658224
  - 0.017380198085367646
  strength: 0.01848824142197912
  strengths:
  - 0.01848824142197912
  signif: 0.9945831472527588
- name: S20
  fix: false
  energy: 0.3449291611793061
  erange:
  - 0.3418069865804786
  - 0.3478863059645039
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.7883298861884864
    - -0.019982928243103434
    - 0.03240768307856136
  dipole:
  - 0.7883298861884864
  - -0.019982928243103434
  - 0.03240768307856136
  strength: 0.035810176715470786
  strengths:
  - 0.035810176715470786
  signif: 1.0
- name: S21
  fix: false
  energy: 0.37220629061215826
  erange:
  - 0.37149116326028947
  - 0.37757048264431475
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 1.231385613546026
    - -0.9138626020990602
    - -0.3175267582720609
  dipole:
  - 1.231385613546026
  - -0.9138626020990602
  - -0.3175267582720609
  strength: 0.15212558854941718
  strengths:
  - 0.15212558854941718
  signif: 1157.7058400623148
- name: S22
  fix: false
  energy: 0.3744327120694623
  erange:
  - 0.3719549785209116
  - 0.37803429790493687
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.0071542858807303175
    - 0.7073787068027645
    - 0.02462279314788993
  dipole:
  - 0.0071542858807303175
  - 0.7073787068027645
  - 0.02462279314788993
  strength: 0.03126775876395046
  strengths:
  - 0.03126775876395046
  signif: 0.9971597340154543
- name: S23
  fix: false
  energy: 0.3860043885285013
  erange:
  - 0.38285463714552964
  - 0.3889339565295549
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.9329358867597572
    - -0.00317073868090815
    - 0.02287081216953694
  dipole:
  - 0.9329358867597572
  - -0.00317073868090815
  - 0.02287081216953694
  strength: 0.056028697600870314
  strengths:
  - 0.056028697600870314
  signif: 1.0
