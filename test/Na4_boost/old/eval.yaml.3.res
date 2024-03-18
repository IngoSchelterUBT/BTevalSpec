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
    max_excit: 22                   # Maximum numer of iterations used to reach relative error between fit and 
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
  energy: 0.1309234285701168
  erange:
  - 0.12590098276091613
  - 0.1319803021449414
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - -0.0010086862320531674
    - 0.08347238765024205
    - 3.33292162301137
  dipole:
  - -0.0010086862320531674
  - 0.08347238765024205
  - 3.33292162301137
  strength: 0.24254296572717327
  strengths:
  - 0.24254296572717327
  signif: 0.07362279690666908
- name: S1
  fix: false
  energy: 0.13208263383856472
  erange:
  - 0.12866034030798737
  - 0.13473965969201265
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - -0.0004533105836031835
    - 0.001442038253450001
    - 6.3021139646013875
  dipole:
  - -0.0004533105836031835
  - 0.001442038253450001
  - 6.3021139646013875
  strength: 0.8743131293449911
  strengths:
  - 0.8743131293449911
  signif: 0.34336851195667295
- name: S2
  fix: false
  energy: 0.1542078618456928
  erange:
  - 0.15272034030798737
  - 0.15879965969201265
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.002230089618001626
    - -0.02056757955136604
    - 1.6727177429283078
  dipole:
  - 0.002230089618001626
  - -0.02056757955136604
  - 1.6727177429283078
  strength: 0.07192287178708325
  strengths:
  - 0.07192287178708325
  signif: 1.0
- name: S3
  fix: false
  energy: 0.16475639423307759
  erange:
  - 0.1618466654591247
  - 0.16792598484315
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 2.0818228063353343
    - 0.04513242400135177
    - -0.0032652446578530503
  dipole:
  - 2.0818228063353343
  - 0.04513242400135177
  - -0.0032652446578530503
  strength: 0.11906488220901144
  strengths:
  - 0.11906488220901144
  signif: 1.0
- name: S4
  fix: false
  energy: 0.19340503798353503
  erange:
  - 0.19024034030798737
  - 0.19631965969201265
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.00015773566471091301
    - 5.314824521095464
    - -0.0011795580916942334
  dipole:
  - 0.00015773566471091301
  - 5.314824521095464
  - -0.0011795580916942334
  strength: 0.9105303246156392
  strengths:
  - 0.9105303246156392
  signif: 1.0
- name: S5
  fix: false
  energy: 0.20558733001968335
  erange:
  - 0.2021985931332427
  - 0.208277912517268
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.014108179395346403
    - -0.13716282915592257
    - 0.9212389747313137
  dipole:
  - 0.014108179395346403
  - -0.13716282915592257
  - 0.9212389747313137
  strength: 0.029731146424485027
  strengths:
  - 0.029731146424485027
  signif: 1.0
- name: S6
  fix: false
  energy: 0.2182268618583869
  erange:
  - 0.2148503403079874
  - 0.22092965969201267
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 4.087252710874259
    - -0.02289929959999129
    - 0.003760210959170294
  dipole:
  - 4.087252710874259
  - -0.02289929959999129
  - 0.003760210959170294
  strength: 0.6076226266255067
  strengths:
  - 0.6076226266255067
  signif: 0.9934977675601832
- name: S7
  fix: false
  energy: 0.2242798009496503
  erange:
  - 0.2214469264490576
  - 0.2275262458330829
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.02932982705846306
    - -0.057075225293588994
    - 1.0936811136220927
  dipole:
  - 0.02932982705846306
  - -0.057075225293588994
  - 1.0936811136220927
  strength: 0.04486553687729574
  strengths:
  - 0.04486553687729574
  signif: 1.000652524602864
- name: S8
  fix: false
  energy: 0.2426431921771962
  erange:
  - 0.2411903403079874
  - 0.24726965969201267
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.018008890938405367
    - 2.200075475166826
    - 0.0009877237272523333
  dipole:
  - 0.018008890938405367
  - 2.200075475166826
  - 0.0009877237272523333
  strength: 0.19575876031251135
  strengths:
  - 0.19575876031251135
  signif: 1.0
- name: S9
  fix: false
  energy: 0.25505812991166943
  erange:
  - 0.2526703403079874
  - 0.2587496596920127
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 2.6929647913285506
    - -0.009613764519766511
    - 0.03275378136173966
  dipole:
  - 2.6929647913285506
  - -0.009613764519766511
  - 0.03275378136173966
  strength: 0.3083323171357762
  strengths:
  - 0.3083323171357762
  signif: 1.0
- name: S10
  fix: false
  energy: 0.2655249731525376
  erange:
  - 0.2617988541231756
  - 0.2678781735072009
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.011690295036932771
    - 0.16820164620898806
    - 0.26409882167411475
  dipole:
  - 0.011690295036932771
  - 0.16820164620898806
  - 0.26409882167411475
  strength: 0.0043447251445367435
  strengths:
  - 0.0043447251445367435
  signif: 1.0
- name: S11
  fix: false
  energy: 0.3012532718674523
  erange:
  - 0.30122315127604954
  - 0.3073024706600748
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 1.088534338341487
    - -0.09304454579786386
    - -0.01900982199297554
  dipole:
  - 1.088534338341487
  - -0.09304454579786386
  - -0.01900982199297554
  strength: 0.05994566890413342
  strengths:
  - 0.05994566890413342
  signif: 0.2259130070308804
- name: S12
  fix: false
  energy: 0.3027517531669455
  erange:
  - 0.293
  - 0.313
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.15360594040033668
    - 0.6733838379548424
    - 0.021293658528913142
  dipole:
  - 0.15360594040033668
  - 0.6733838379548424
  - 0.021293658528913142
  strength: 0.024093690838813932
  strengths:
  - 0.024093690838813932
  signif: 0.4740054677045553
- name: S13
  fix: false
  energy: 0.31149486574851637
  erange:
  - 0.30864419544600225
  - 0.31472351483002753
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - -0.06336470037842784
    - -0.08719505292769052
    - 0.5787285615291219
  dipole:
  - -0.06336470037842784
  - -0.08719505292769052
  - 0.5787285615291219
  strength: 0.017991154867240143
  strengths:
  - 0.017991154867240143
  signif: 1.0014478902564166
- name: S14
  fix: false
  energy: 0.31444920954904454
  erange:
  - 0.3130504404219117
  - 0.319129759805937
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.5320442501329992
    - 0.06445626059129854
    - 0.023871907951967617
  dipole:
  - 0.5320442501329992
  - 0.06445626059129854
  - 0.023871907951967617
  strength: 0.015082847810431182
  strengths:
  - 0.015082847810431182
  signif: 0.994775529619566
- name: S15
  fix: false
  energy: 0.3182162500253837
  erange:
  - 0.307
  - 0.327
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.004893792454140703
    - 0.7832648271280315
    - -0.0199126490215393
  dipole:
  - 0.004893792454140703
  - 0.7832648271280315
  - -0.0199126490215393
  strength: 0.03256011222306591
  strengths:
  - 0.03256011222306591
  signif: 0.9705074499403322
- name: S16
  fix: false
  energy: 0.3335411711382863
  erange:
  - 0.3304435126952384
  - 0.3365228320792637
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.9293406006235783
    - 0.01700325685247705
    - 0.010627654176938598
  dipole:
  - 0.9293406006235783
  - 0.01700325685247705
  - 0.010627654176938598
  strength: 0.04803415404043457
  strengths:
  - 0.04803415404043457
  signif: 0.9906764800721718
- name: S17
  fix: false
  energy: 0.3369870406600552
  erange:
  - 0.33392212714990377
  - 0.34000144653392905
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.09261365167603572
    - 0.5764103327663503
    - -0.02962019261703398
  dipole:
  - 0.09261365167603572
  - 0.5764103327663503
  - -0.02962019261703398
  strength: 0.019191609440403445
  strengths:
  - 0.019191609440403445
  signif: 0.9844337188592404
- name: S18
  fix: false
  energy: 0.34491879695898564
  erange:
  - 0.3418069865804786
  - 0.3478863059645039
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.7723188388371662
    - 0.006497568957048946
    - 0.03902771015838901
  dipole:
  - 0.7723188388371662
  - 0.006497568957048946
  - 0.03902771015838901
  strength: 0.034379307937936206
  strengths:
  - 0.034379307937936206
  signif: 1.0
- name: S19
  fix: false
  energy: 0.3735564406557027
  erange:
  - 0.37149116326028947
  - 0.37757048264431475
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.13533383848161915
    - -3.4107407344092646
    - -0.08871607533698174
  dipole:
  - 0.13533383848161915
  - -3.4107407344092646
  - -0.08871607533698174
  strength: 0.7259034764461275
  strengths:
  - 0.7259034764461275
  signif: 0.00015602906783944548
- name: S20
  fix: false
  energy: 0.3735967254379836
  erange:
  - 0.3719549785209116
  - 0.37803429790493687
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - -0.12862438950845512
    - 3.474316211082447
    - 0.08720157070154218
  dipole:
  - -0.12862438950845512
  - 3.474316211082447
  - 0.08720157070154218
  strength: 0.7531100699690754
  strengths:
  - 0.7531100699690754
  signif: 0.007159735632921456
- name: S21
  fix: false
  energy: 0.386000781939887
  erange:
  - 0.38285463714552964
  - 0.3889339565295549
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.9284082801614867
    - 0.03139156242935161
    - -0.004154070209283246
  dipole:
  - 0.9284082801614867
  - 0.03139156242935161
  - -0.004154070209283246
  strength: 0.05551621642331977
  strengths:
  - 0.05551621642331977
  signif: 1.0
