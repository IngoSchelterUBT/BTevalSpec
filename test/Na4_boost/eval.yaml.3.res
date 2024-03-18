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
    smooth: 0.                      # Smooth for Pade Approximation
    thin: 0                         # Only keep every 2^n data point
  Fit:
    calc: true                      # Turn fit on/off
    skipfirst: true                 # Skip first fit of existing excitations
    firstsingle: false
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
    fiterr: 0.004913860277095138
    reset_erange: false
EXT: {}
SPEC:
- name: S0
  fix: false
  energy: 0.13094740206464803
  energyErr: 1.1473672094258704e-05
  erange:
  - 0.12791048560930424
  - 0.13398980499332952
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - -0.000567587094925619
    - 0.1124013011644809
    - 3.380672635121145
  dipolesErr:
  - - 0.002217631138313371
    - 0.0021919420411705966
    - 0.03216079182844161
  ampl:
  - - - -1.2795702055324906e-05
      - 0.0025339786143658213
      - 0.07621399459631684
  dipole:
  - -0.000567587094925619
  - 0.1124013011644809
  - 3.380672635121145
  dipoleErr:
  - 0.002217631138313371
  - 0.0021919420411705966
  - 0.03216079182844161
  strength: 0.24970756958759366
  strengthErr: 0.00475645691853103
  strengths:
  - 0.24970756958759366
  strengthsErr:
  - 0.00475645691853103
  signifFit: 0.1320997519696956
  signifErr: 0.9890025693380489
  signifAng: 0.772087943179549
  signifExc: 0.9999507735917735
  signifRng: 0.9999999999993366
- name: S1
  fix: false
  energy: 0.1320890533676624
  energyErr: 3.4511959228733146e-06
  erange:
  - 0.1290524587637086
  - 0.1351317781477339
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - -0.0007097547454771877
    - -0.009695029558515195
    - 6.273826401787114
  dipolesErr:
  - - 0.00123516629915109
    - 0.0014027383404283597
    - 0.01726108244887432
  ampl:
  - - - -2.8695536849652324e-05
      - -0.00039197212801697715
      - 0.25365214934880653
  dipole:
  - -0.0007097547454771877
  - -0.009695029558515195
  - 6.273826401787114
  dipoleErr:
  - 0.00123516629915109
  - 0.0014027383404283597
  - 0.01726108244887432
  strength: 0.8665260335950197
  strengthErr: 0.004767471703654248
  strengths:
  - 0.8665260335950197
  strengthsErr:
  - 0.004767471703654248
  signifFit: 0.13210011942994002
  signifErr: 0.9968235217815206
  signifAng: 0.759204896390138
  signifExc: 0.9999913010138146
  signifRng: 0.9999999999989662
- name: S2
  fix: false
  energy: 0.15420744115182247
  energyErr: 3.0268042994863566e-06
  erange:
  - 0.15108195632629087
  - 0.15716127571031616
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - 0.002634761435135556
    - 0.0017297477194595308
    - 1.660690588422732
  dipolesErr:
  - - 0.002887871812985411
    - 0.0028942451134160645
    - 0.002505766386391039
  ampl:
  - - - 2.8318112940523157e-05
      - 1.859112958959308e-05
      - 0.01784891148590812
  dipole:
  - 0.002634761435135556
  - 0.0017297477194595308
  - 1.660690588422732
  dipoleErr:
  - 0.002887871812985411
  - 0.0028942451134160645
  - 0.002505766386391039
  strength: 0.07088153165626017
  strengthErr: 0.00021454983623056447
  strengths:
  - 0.07088153165626017
  strengthsErr:
  - 0.00021454983623056447
  signifFit: 1.0
  signifErr: 0.9982524304595567
  signifAng: 0.7608328184636618
  signifExc: 0.9995774712015046
  signifRng: 0.9999993644377512
- name: S3
  fix: false
  energy: 0.16475745506087427
  energyErr: 1.8293213382049059e-06
  erange:
  - 0.16171739475593006
  - 0.16779671413995534
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - 2.1091108480873917
    - -0.009151166767696256
    - -0.0036632331583691864
  dipolesErr:
  - - 0.0020039730952763077
    - 0.0023153140156028746
    - 0.002316045500555853
  ampl:
  - - - 0.028539509599867476
      - -0.00012382934356128914
      - -4.956917176005076e-05
  dipole:
  - 2.1091108480873917
  - -0.009151166767696256
  - -0.0036632331583691864
  dipoleErr:
  - 0.0020039730952763077
  - 0.0023153140156028746
  - 0.002316045500555853
  strength: 0.12215243298131159
  strengthErr: 0.00023049190819366078
  strengths:
  - 0.12215243298131159
  strengthsErr:
  - 0.00023049190819366078
  signifFit: 0.9996681913631267
  signifErr: 0.9989105860440607
  signifAng: 0.7575197519599204
  signifExc: 0.9990235906302053
  signifRng: 0.9999999999999997
- name: S4
  fix: false
  energy: 0.19244077153742412
  energyErr: 9.934994001039144e-06
  erange:
  - 0.18861323555790627
  - 0.19469255494193155
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - 1.214859635695733e-05
    - 3.1037294316807706
    - -0.004735313609335037
  dipolesErr:
  - - 0.002190785138328204
    - 0.01999544681076266
    - 0.0021959267849028245
  ampl:
  - - - 2.430205276212755e-07
      - 0.06208700510893036
      - -9.472521581752715e-05
  dipole:
  - 1.214859635695733e-05
  - 3.1037294316807706
  - -0.004735313609335037
  dipoleErr:
  - 0.002190785138328204
  - 0.01999544681076266
  - 0.0021959267849028245
  strength: 0.3089687522392748
  strengthErr: 0.00398032325941347
  strengths:
  - 0.3089687522392748
  strengthsErr:
  - 0.00398032325941347
  signifFit: 0.7859250524064791
  signifErr: 0.9925622229153257
  signifAng: 0.759256875792806
  signifExc: 0.9995595038680507
  signifRng: 0.9954863044606483
- name: S5
  fix: false
  energy: 0.1938199787454929
  energyErr: 4.577208445779438e-06
  erange:
  - 0.19045738441768245
  - 0.19653670380170774
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - -0.00012315315732474682
    - 4.553125622020686
    - 0.003642149398064886
  dipolesErr:
  - - 0.0014866752543155214
    - 0.013647053957500004
    - 0.001492150827302572
  ampl:
  - - - -3.6223057997517633e-06
      - 0.13392115724774653
      - 0.00010712659890142977
  dipole:
  - -0.00012315315732474682
  - 4.553125622020686
  - 0.003642149398064886
  dipoleErr:
  - 0.0014866752543155214
  - 0.013647053957500004
  - 0.001492150827302572
  strength: 0.6696792383762986
  strengthErr: 0.0040147883769903366
  strengths:
  - 0.6696792383762986
  strengthsErr:
  - 0.0040147883769903366
  signifFit: 0.785926584691766
  signifErr: 0.9965387322506573
  signifAng: 0.7601291361450171
  signifExc: 0.9999343775358024
  signifRng: 0.9998726031616302
- name: S6
  fix: false
  energy: 0.20540987523542362
  energyErr: 1.1182295869112925e-05
  erange:
  - 0.20234218755605862
  - 0.2084215069400839
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - 0.012028944175770901
    - -0.01577320238109348
    - 0.8631705491907088
  dipolesErr:
  - - 0.005638065818435868
    - 0.005783363444281179
    - 0.004929422022299454
  ampl:
  - - - 6.673144773049886e-05
      - -8.750299401643565e-05
      - 0.0047885017624279
  dipole:
  - 0.012028944175770901
  - -0.01577320238109348
  - 0.8631705491907088
  dipoleErr:
  - 0.005638065818435868
  - 0.005783363444281179
  - 0.004929422022299454
  strength: 0.02552070100485422
  strengthErr: 0.0002897328052544874
  strengths:
  - 0.02552070100485422
  strengthsErr:
  - 0.0002897328052544874
  signifFit: 1.0
  signifErr: 0.9934454263981576
  signifAng: 0.7580858163744567
  signifExc: 0.9977810005253906
  signifRng: 0.9999999927711737
- name: S7
  fix: false
  energy: 0.21822679297028866
  energyErr: 5.159440529007684e-07
  erange:
  - 0.21518700797897883
  - 0.2212663273630041
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - 4.076517239537305
    - 0.00012558516856691184
    - 0.0018971143130413425
  dipolesErr:
  - - 0.0010378402327208842
    - 0.0011831437748778255
    - 0.0012369817101164026
  ampl:
  - - - 0.10732191294531739
      - 3.306263591244107e-06
      - 4.994506957479417e-05
  dipole:
  - 4.076517239537305
  - 0.00012558516856691184
  - 0.0018971143130413425
  dipoleErr:
  - 0.0010378402327208842
  - 0.0011831437748778255
  - 0.0012369817101164026
  strength: 0.6044153440205048
  strengthErr: 0.0003079377218372136
  strengths:
  - 0.6044153440205048
  strengthsErr:
  - 0.0003079377218372136
  signifFit: 0.9995844804540815
  signifErr: 0.9997058514011012
  signifAng: 0.7600241298042723
  signifExc: 0.9998164235289679
  signifRng: 1.0
- name: S8
  fix: false
  energy: 0.22427955722061993
  energyErr: 7.429248999941751e-06
  erange:
  - 0.2212308883589597
  - 0.227310207742985
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - 0.027278167213184844
    - 0.0036942063047126377
    - 1.0619871354176094
  dipolesErr:
  - - 0.00452685244657784
    - 0.004408435064602583
    - 0.003831278555554696
  ampl:
  - - - 0.00019244811536694732
      - 2.6062712922113368e-05
      - 0.007492344377751178
  dipole:
  - 0.027278167213184844
  - 0.0036942063047126377
  - 1.0619871354176094
  dipoleErr:
  - 0.00452685244657784
  - 0.004408435064602583
  - 0.003831278555554696
  strength: 0.042186028585494306
  strengthErr: 0.00031462986907966736
  strengths:
  - 0.042186028585494306
  strengthsErr:
  - 0.00031462986907966736
  signifFit: 0.9997664606659726
  signifErr: 0.9956940331740326
  signifAng: 0.7707067666370301
  signifExc: 0.9981585881314589
  signifRng: 0.9999999999228316
- name: S9
  fix: false
  energy: 0.24268255625009402
  energyErr: 1.6792923861016269e-06
  erange:
  - 0.2396400849302561
  - 0.24571940431428138
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - 0.01696182100237366
    - 2.2121972925123674
    - 0.000331352793417358
  dipolesErr:
  - - 0.0021690864391435993
    - 0.001869247433682289
    - 0.0021619880907036114
  ampl:
  - - - 0.00024410265197187254
      - 0.03183639455408098
      - 4.768597404733024e-06
  dipole:
  - 0.01696182100237366
  - 2.2121972925123674
  - 0.000331352793417358
  dipoleErr:
  - 0.0021690864391435993
  - 0.001869247433682289
  - 0.0021619880907036114
  strength: 0.1979523054803337
  strengthErr: 0.0003375434602703956
  strengths:
  - 0.1979523054803337
  strengthsErr:
  - 0.0003375434602703956
  signifFit: 0.999946223921937
  signifErr: 0.9990155163529031
  signifAng: 0.7627885804236985
  signifExc: 0.9993683835520671
  signifRng: 0.9999999999992679
- name: S10
  fix: false
  energy: 0.2550586762388026
  energyErr: 1.1330757038104041e-06
  erange:
  - 0.2520189277104945
  - 0.2580982470945198
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - 2.6814035181377545
    - 0.013342351852300656
    - 0.03259660977591735
  dipolesErr:
  - - 0.001528759964339388
    - 0.0017794114068577659
    - 0.0017622385702652258
  ampl:
  - - - 0.04720589695949235
      - 0.00023489104958526402
      - 0.0005738607381925885
  dipole:
  - 2.6814035181377545
  - 0.013342351852300656
  - 0.03259660977591735
  dipoleErr:
  - 0.001528759964339388
  - 0.0017794114068577659
  - 0.0017622385702652258
  strength: 0.3056948539267035
  strengthErr: 0.00035541651734238605
  strengths:
  - 0.3056948539267035
  strengthsErr:
  - 0.00035541651734238605
  signifFit: 1.0
  signifErr: 0.9993287429627079
  signifAng: 0.7662839081747959
  signifExc: 0.9995271207680868
  signifRng: 1.0
- name: S11
  fix: false
  energy: 0.2660823133157124
  energyErr: 6.767145680483105e-05
  erange:
  - 0.26294176705944855
  - 0.26902108644347383
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - 0.008211941737951235
    - 0.13293542739423242
    - 0.29132799338111137
  dipolesErr:
  - - 0.011344770172015272
    - 0.009978540221781343
    - 0.009355794158614223
  ampl:
  - - - 2.2924592956619118e-05
      - 0.0003711047471809364
      - 0.0008132760653020438
  dipole:
  - 0.008211941737951235
  - 0.13293542739423242
  - 0.29132799338111137
  dipoleErr:
  - 0.011344770172015272
  - 0.009978540221781343
  - 0.009355794158614223
  strength: 0.004550507063872639
  strengthErr: 0.00036766103634059967
  strengths:
  - 0.004550507063872639
  strengthsErr:
  - 0.00036766103634059967
  signifFit: 0.9992319851241003
  signifErr: 0.9533526274495429
  signifAng: 0.8828798664603229
  signifExc: 0.9615720501470234
  signifRng: 0.9999987865176312
- name: S12
  fix: false
  energy: 0.2707250760465725
  energyErr: 0.00012462488702897346
  erange:
  - 0.2676274565496985
  - 0.2737067759337238
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - -0.059252159482692704
    - -0.0032533542024681325
    - 0.30506248845616885
  dipolesErr:
  - - 0.023125730815136233
    - 0.02047037549499134
    - 0.019596262450479872
  ampl:
  - - - -9.277102452434967e-05
      - -5.093772195623125e-06
      - 0.0004776359181692435
  dipole:
  - -0.059252159482692704
  - -0.0032533542024681325
  - 0.30506248845616885
  dipoleErr:
  - 0.023125730815136233
  - 0.02047037549499134
  - 0.019596262450479872
  strength: 0.00435797546030706
  strengthErr: 0.0004098093351369719
  strengths:
  - 0.00435797546030706
  strengthsErr:
  - 0.0004098093351369719
  signifFit: 0.9995075023621262
  signifErr: 0.9457079251333
  signifAng: 0.6712737791534394
  signifExc: 0.9279354193513788
  signifRng: 0.9999998678072003
- name: S13
  fix: false
  energy: 0.3013030899964411
  energyErr: 1.536402865017967e-05
  erange:
  - 0.2982635364452134
  - 0.3043428558292387
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - 1.052394321893605
    - 0.01823793538775578
    - -0.009860865505008334
  dipolesErr:
  - - 0.006123095903376236
    - 0.010161149379243936
    - 0.0052721564388912906
  ampl:
  - - - 0.007206007290219442
      - 0.00012487970775654363
      - -6.751981385561604e-05
  dipole:
  - 1.052394321893605
  - 0.01823793538775578
  - -0.009860865505008334
  dipoleErr:
  - 0.006123095903376236
  - 0.010161149379243936
  - 0.0052721564388912906
  strength: 0.055638812780782934
  strengthErr: 0.0006605816378782157
  strengths:
  - 0.055638812780782934
  strengthsErr:
  - 0.0006605816378782157
  signifFit: 0.9998465646229266
  signifErr: 0.9931453068929756
  signifAng: 0.7627798391844056
  signifExc: 0.9951636541167208
  signifRng: 1.0
- name: S14
  fix: false
  energy: 0.3032232549079305
  energyErr: 3.219417789929421e-05
  erange:
  - 0.30016129033765315
  - 0.30624060972167844
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - 0.09617201948543143
    - 0.6538918604270951
    - 0.0160692452701724
  dipolesErr:
  - - 0.0136157526647129
    - 0.008241336130371936
    - 0.00715741806512556
  ampl:
  - - - 0.0004756060586190583
      - 0.003233736092522472
      - 7.946833651639581e-05
  dipole:
  - 0.09617201948543143
  - 0.6538918604270951
  - 0.0160692452701724
  dipoleErr:
  - 0.0136157526647129
  - 0.008241336130371936
  - 0.00715741806512556
  strength: 0.02208889652140768
  strengthErr: 0.0006886619561956124
  strengths:
  - 0.02208889652140768
  strengthsErr:
  - 0.0006886619561956124
  signifFit: 1.000020336877065
  signifErr: 0.9820000439856584
  signifAng: 0.8179582084089013
  signifExc: 0.9966494219609088
  signifRng: 0.9999999971006571
- name: S15
  fix: false
  energy: 0.31153022657384316
  energyErr: 4.18904167545251e-05
  erange:
  - 0.30850657081229577
  - 0.31458589019632105
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - -0.05113054422558606
    - -0.008440959734876522
    - 0.529169048610112
  dipolesErr:
  - - 0.015424615683216163
    - 0.011113174499275884
    - 0.011108091788913263
  ampl:
  - - - -0.00015498891397110082
      - -2.558656869385802e-05
      - 0.0016040380049418858
  dipole:
  - -0.05113054422558606
  - -0.008440959734876522
  - 0.529169048610112
  dipoleErr:
  - 0.015424615683216163
  - 0.011113174499275884
  - 0.011108091788913263
  strength: 0.014678549561006683
  strengthErr: 0.0005187587854336362
  strengths:
  - 0.014678549561006683
  strengthsErr:
  - 0.0005187587854336362
  signifFit: 0.9999248122716263
  signifErr: 0.9795957002992509
  signifAng: 0.7140835759174916
  signifExc: 0.9915940851801909
  signifRng: 0.9999999992315657
- name: S16
  fix: false
  energy: 0.3144543907108879
  energyErr: 4.007420891002545e-05
  erange:
  - 0.3114102770950735
  - 0.3174895964790988
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - 0.5743129003012405
    - -0.02319285161006319
    - 0.020408244074156808
  dipolesErr:
  - - 0.008432737686506585
    - 0.009842200285008101
    - 0.01038009716201349
  ampl:
  - - - 0.002118754843472665
      - -8.556305570881155e-05
      - 7.529008308225643e-05
  dipole:
  - 0.5743129003012405
  - -0.02319285161006319
  - 0.020408244074156808
  dipoleErr:
  - 0.008432737686506585
  - 0.009842200285008101
  - 0.01038009716201349
  strength: 0.017336379569541472
  strengthErr: 0.0005059157591677484
  strengths:
  - 0.017336379569541472
  strengthsErr:
  - 0.0005059157591677484
  signifFit: 0.9999223557706813
  signifErr: 0.9831515802609336
  signifAng: 0.7574440440312158
  signifExc: 0.9914925504018859
  signifRng: 0.9999999999953904
- name: S17
  fix: false
  energy: 0.31793573304302797
  energyErr: 1.703994017976589e-05
  erange:
  - 0.31492208913769465
  - 0.32100140852171993
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - 0.006076978628339058
    - 0.77207867546858
    - -0.020386700929158074
  dipolesErr:
  - - 0.007839291028380173
    - 0.0061543969142677235
    - 0.006635063471336666
  ampl:
  - - - 2.9724798480428778e-05
      - 0.0037765285091371967
      - -9.97190567651579e-05
  dipole:
  - 0.006076978628339058
  - 0.77207867546858
  - -0.020386700929158074
  dipoleErr:
  - 0.007839291028380173
  - 0.0061543969142677235
  - 0.006635063471336666
  strength: 0.03161118567474912
  strengthErr: 0.0004942897600471261
  strengths:
  - 0.03161118567474912
  strengthsErr:
  - 0.0004942897600471261
  signifFit: 0.9999606444528724
  signifErr: 0.9909722359370647
  signifAng: 0.7526185468397287
  signifExc: 0.9990959835757605
  signifRng: 0.999999994634033
- name: S18
  fix: false
  energy: 0.3335276519248526
  energyErr: 1.1281418237328155e-05
  erange:
  - 0.3304875215987394
  - 0.3365668409827647
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - 0.925087157227307
    - 0.019307281120082072
    - 0.010083363910943685
  dipolesErr:
  - - 0.004487196994435593
    - 0.005806563675802482
    - 0.005024982853169693
  ampl:
  - - - 0.0056995804326944845
      - 0.00011895463126996172
      - 6.21248962257844e-05
  dipole:
  - 0.925087157227307
  - 0.019307281120082072
  - 0.010083363910943685
  dipoleErr:
  - 0.004487196994435593
  - 0.005806563675802482
  - 0.005024982853169693
  strength: 0.047597769773049324
  strengthErr: 0.00047959371941276986
  strengths:
  - 0.047597769773049324
  strengthsErr:
  - 0.00047959371941276986
  signifFit: 0.9996482565077203
  signifErr: 0.9941826357763132
  signifAng: 0.7717046205365898
  signifExc: 0.9982677984957616
  signifRng: 0.9999999999999994
- name: S19
  fix: false
  energy: 0.33676333851552964
  energyErr: 2.8140958448076647e-05
  erange:
  - 0.3337417998961406
  - 0.3398211192801659
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - 0.10300369472498076
    - 0.5441713572486168
    - -0.028048139697105987
  dipolesErr:
  - - 0.008450038023171076
    - 0.006886030575190377
    - 0.007985284374345163
  ampl:
  - - - 0.00041164883224850847
      - 0.002174752122752273
      - -0.00011209291068523726
  dipole:
  - 0.10300369472498076
  - 0.5441713572486168
  - -0.028048139697105987
  dipoleErr:
  - 0.008450038023171076
  - 0.006886030575190377
  - 0.007985284374345163
  strength: 0.017260183324425117
  strengthErr: 0.0004932006295673126
  strengths:
  - 0.017260183324425117
  strengthsErr:
  - 0.0004932006295673126
  signifFit: 0.9996486219287443
  signifErr: 0.9835025207500232
  signifAng: 0.8028632247253219
  signifExc: 0.9937763967870188
  signifRng: 0.9999999987369058
- name: S20
  fix: false
  energy: 0.34490794024289195
  energyErr: 1.3363446859535156e-05
  erange:
  - 0.34186675890591506
  - 0.34794607828994034
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - 0.7751605599287221
    - -0.003996837960058472
    - 0.03791930497806111
  dipolesErr:
  - - 0.0051656360894059205
    - 0.006036819243001659
    - 0.005860672206555256
  ampl:
  - - - 0.004048360331060867
      - -2.0873920944412864e-05
      - 0.00019803769436966347
  dipole:
  - 0.7751605599287221
  - -0.003996837960058472
  - 0.03791930497806111
  dipoleErr:
  - 0.0051656360894059205
  - 0.006036819243001659
  - 0.005860672206555256
  strength: 0.03462460347824161
  strengthErr: 0.00048313597334489415
  strengths:
  - 0.03462460347824161
  strengthsErr:
  - 0.00048313597334489415
  signifFit: 1.0
  signifErr: 0.9919439168612733
  signifAng: 0.7758147191678821
  signifExc: 0.9956729319947444
  signifRng: 0.9999999999999372
- name: S21
  fix: false
  energy: 0.36579385036838413
  energyErr: 6.562189444265093e-05
  erange:
  - 0.36276487657065504
  - 0.3688441959546803
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - 0.3484382230232492
    - 0.03872871397091953
    - -0.022926617794444545
  dipolesErr:
  - - 0.01139106549253601
    - 0.012785449698389007
    - 0.013687981340289382
  ampl:
  - - - 0.000819234457027183
      - 9.105745255511283e-05
      - -5.3904175946415326e-05
  dipole:
  - 0.3484382230232492
  - 0.03872871397091953
  - -0.022926617794444545
  dipoleErr:
  - 0.01139106549253601
  - 0.012785449698389007
  - 0.013687981340289382
  strength: 0.0075252780352222624
  strengthErr: 0.0005060671961174331
  strengths:
  - 0.0075252780352222624
  strengthsErr:
  - 0.0005060671961174331
  signifFit: 0.9999857273787037
  signifErr: 0.9611737891225794
  signifAng: 0.773667484932966
  signifExc: 0.9738289770128454
  signifRng: 0.9999999998472632
- name: S22
  fix: false
  energy: 0.3744895597589312
  energyErr: 1.6973365360367372e-05
  erange:
  - 0.3714424843692348
  - 0.3775218037532601
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - -0.005897018569624863
    - 0.6982251396993999
    - 0.009235926224521284
  dipolesErr:
  - - 0.00705646455254966
    - 0.005985248356468117
    - 0.006840486029510462
  ampl:
  - - - -2.6705100861455692e-05
      - 0.0031619660951589363
      - 4.182560025251792e-05
  dipole:
  - -0.005897018569624863
  - 0.6982251396993999
  - 0.009235926224521284
  dipoleErr:
  - 0.00705646455254966
  - 0.005985248356468117
  - 0.006840486029510462
  strength: 0.0304359163836553
  strengthErr: 0.0005243626571140931
  strengths:
  - 0.0304359163836553
  strengthsErr:
  - 0.0005243626571140931
  signifFit: 1.0
  signifErr: 0.9900531688475551
  signifAng: 0.7616033913725753
  signifExc: 0.9947570223079883
  signifRng: 0.9999999999645751
- name: S23
  fix: false
  energy: 0.38600778563145977
  energyErr: 9.456284135329386e-06
  erange:
  - 0.38296812949097503
  - 0.3890474488750003
  phase: 0.0
  phaseErr: 0.0
  tmod: 1.0
  dipoles:
  - - 0.9431928104063375
    - -0.004546658713196588
    - -0.0033625907518559156
  dipolesErr:
  - - 0.004494261511320379
    - 0.00521063553513648
    - 0.005168802749402652
  ampl:
  - - - 0.005694271693106399
      - -2.744922323741765e-05
      - -2.0300732917530992e-05
  dipole:
  - 0.9431928104063375
  - -0.004546658713196588
  - -0.0033625907518559156
  dipoleErr:
  - 0.004494261511320379
  - 0.00521063553513648
  - 0.005168802749402652
  strength: 0.05723496065682965
  strengthErr: 0.0005401388548773856
  strengths:
  - 0.05723496065682965
  strengthsErr:
  - 0.0005401388548773856
  signifFit: 0.9999633678665129
  signifErr: 0.9945514191032053
  signifAng: 0.7566363352027575
  signifExc: 0.9958135791393952
  signifRng: 1.0
