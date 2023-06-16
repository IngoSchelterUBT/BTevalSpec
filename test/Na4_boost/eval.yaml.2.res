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
    max_excit: 19                   # Maximum numer of iterations used to reach relative error between fit and 
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
  energy: 0.13092095209679538
  erange:
  - 0.12590098276091613
  - 0.1319803021449414
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - -0.0024708896619320036
    - 0.08385391861464857
    - 3.3248203015034186
  dipole:
  - -0.0024708896619320036
  - 0.08385391861464857
  - 3.3248203015034186
  strength: 0.24136297878218538
  strengths:
  - 0.24136297878218538
  signif: 0.07125983425881766
- name: S1
  fix: false
  energy: 0.13208157632422374
  erange:
  - 0.12866034030798737
  - 0.13473965969201265
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.0004830675721209304
    - 0.0014569013026724848
    - 6.306120168488658
  dipole:
  - 0.0004830675721209304
  - 0.0014569013026724848
  - 6.306120168488658
  strength: 0.8754180629511863
  strengths:
  - 0.8754180629511863
  signif: 0.3428508872673884
- name: S2
  fix: false
  energy: 0.15420933174963564
  erange:
  - 0.15272034030798737
  - 0.15879965969201265
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.0006487840899202583
    - -0.020543036759329492
    - 1.6734030873357764
  dipole:
  - 0.0006487840899202583
  - -0.020543036759329492
  - 1.6734030873357764
  strength: 0.07198235438632843
  strengths:
  - 0.07198235438632843
  signif: 1.0
- name: S3
  fix: false
  energy: 0.1647545278868802
  erange:
  - 0.1618466654591247
  - 0.16792598484315
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 2.082096596302268
    - 0.044818152152435384
    - -0.0039023223392553276
  dipole:
  - 2.082096596302268
  - 0.044818152152435384
  - -0.0039023223392553276
  strength: 0.11909418704048802
  strengths:
  - 0.11909418704048802
  signif: 1.0
- name: S4
  fix: false
  energy: 0.19340492005025509
  erange:
  - 0.19024034030798737
  - 0.19631965969201265
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.0008062179902821363
    - 5.314377590214725
    - -0.0009329040061134493
  dipole:
  - 0.0008062179902821363
  - 5.314377590214725
  - -0.0009329040061134493
  strength: 0.9103766438056086
  strengths:
  - 0.9103766438056086
  signif: 1.0
- name: S5
  fix: false
  energy: 0.20557959786616353
  erange:
  - 0.2021985931332427
  - 0.208277912517268
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.018661195800857573
    - -0.1336642518838075
    - 0.9177510305556585
  dipole:
  - 0.018661195800857573
  - -0.1336642518838075
  - 0.9177510305556585
  strength: 0.02948290084712357
  strengths:
  - 0.02948290084712357
  signif: 1.0
- name: S6
  fix: false
  energy: 0.21822671470361918
  erange:
  - 0.2148503403079874
  - 0.22092965969201267
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 4.087748610162382
    - -0.023066539075571494
    - 0.0039031868496358976
  dipole:
  - 4.087748610162382
  - -0.023066539075571494
  - 0.0039031868496358976
  strength: 0.6077699840328544
  strengths:
  - 0.6077699840328544
  signif: 0.9930995744371214
- name: S7
  fix: false
  energy: 0.22427463332428535
  erange:
  - 0.2214469264490576
  - 0.2275262458330829
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.033060014900543055
    - -0.05701619752653122
    - 1.0916945947177759
  dipole:
  - 0.033060014900543055
  - -0.05701619752653122
  - 1.0916945947177759
  strength: 0.044710676939984345
  strengths:
  - 0.044710676939984345
  signif: 1.00064299256848
- name: S8
  fix: false
  energy: 0.24264133301752686
  erange:
  - 0.2411903403079874
  - 0.24726965969201267
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.020646658805740972
    - 2.1985490007766155
    - 0.0007647747416327695
  dipole:
  - 0.020646658805740972
  - 2.1985490007766155
  - 0.0007647747416327695
  strength: 0.1954898367207505
  strengths:
  - 0.1954898367207505
  signif: 1.0
- name: S9
  fix: false
  energy: 0.25505657082657374
  erange:
  - 0.2526703403079874
  - 0.2587496596920127
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 2.694791913992428
    - -0.010790672167828806
    - 0.03240552284770469
  dipole:
  - 2.694791913992428
  - -0.010790672167828806
  - 0.03240552284770469
  strength: 0.30874895533573676
  strengths:
  - 0.30874895533573676
  signif: 1.0
- name: S10
  fix: false
  energy: 0.265495868233634
  erange:
  - 0.2617988541231756
  - 0.2678781735072009
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.017216671974404123
    - 0.17515040436771268
    - 0.2643191215710493
  dipole:
  - 0.017216671974404123
  - 0.17515040436771268
  - 0.2643191215710493
  strength: 0.0044620419714588395
  strengths:
  - 0.0044620419714588395
  signif: 1.0
- name: S11
  fix: false
  energy: 0.3015859421928086
  erange:
  - 0.30122315127604954
  - 0.3073024706600748
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.9438843205744626
    - 0.27375433120445297
    - 0.005048966949810828
  dipole:
  - 0.9438843205744626
  - 0.27375433120445297
  - 0.005048966949810828
  strength: 0.04854953300065377
  strengths:
  - 0.04854953300065377
  signif: 1.0
- name: S12
  fix: false
  energy: 0.312300579654844
  erange:
  - 0.30864419544600225
  - 0.31472351483002753
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - -0.888264274980792
    - 0.6993119844977432
    - 0.7267647121319729
  dipole:
  - -0.888264274980792
  - 0.6993119844977432
  - 0.7267647121319729
  strength: 0.09401484265491672
  strengths:
  - 0.09401484265491672
  signif: 0.16157530777399598
- name: S13
  fix: false
  energy: 0.31314230453896874
  erange:
  - 0.3130504404219117
  - 0.319129759805937
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 3.61785626918299
    - -2.427720437257384
    - -1.0070089165000937
  dipole:
  - 3.61785626918299
  - -2.427720437257384
  - -1.0070089165000937
  strength: 1.0436394963973876
  strengths:
  - 1.0436394963973876
  signif: 0.2232314422671581
- name: S14
  fix: false
  energy: 0.33353079607888325
  erange:
  - 0.3304435126952384
  - 0.3365228320792637
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.9175327817966318
    - 0.0330397016976465
    - 0.007778718675386285
  dipole:
  - 0.9175327817966318
  - 0.0330397016976465
  - 0.007778718675386285
  strength: 0.04686210723822134
  strengths:
  - 0.04686210723822134
  signif: 0.9829679078929691
- name: S15
  fix: false
  energy: 0.337071784855746
  erange:
  - 0.33392212714990377
  - 0.34000144653392905
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.1101346837192566
    - 0.5677558871331766
    - -0.030701144863128894
  dipole:
  - 0.1101346837192566
  - 0.5677558871331766
  - -0.030701144863128894
  strength: 0.018843377708383353
  strengths:
  - 0.018843377708383353
  signif: 0.9778679288034723
- name: S16
  fix: false
  energy: 0.34488073000094865
  erange:
  - 0.3418069865804786
  - 0.3478863059645039
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.7658931766660395
    - 0.015062945166829088
    - 0.03376508892044207
  dipole:
  - 0.7658931766660395
  - 0.015062945166829088
  - 0.03376508892044207
  strength: 0.0337959739161166
  strengths:
  - 0.0337959739161166
  signif: 1.0
- name: S17
  fix: false
  energy: 0.3735472025809825
  erange:
  - 0.37149116326028947
  - 0.37757048264431475
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.12621298972465117
    - -2.4714928530026685
    - -0.023149423788004685
  dipole:
  - 0.12621298972465117
  - -2.4714928530026685
  - -0.023149423788004685
  strength: 0.38131340687335397
  strengths:
  - 0.38131340687335397
  signif: 0.00035712920941162837
- name: S18
  fix: false
  energy: 0.3736245954132776
  erange:
  - 0.3719549785209116
  - 0.37803429790493687
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - -0.11319128888417114
    - 2.5572196954376696
    - 0.022442068988128318
  dipole:
  - -0.11319128888417114
  - 2.5572196954376696
  - 0.022442068988128318
  strength: 0.40804093068417263
  strengths:
  - 0.40804093068417263
  signif: 0.0041366798830327
- name: S19
  fix: false
  energy: 0.38599335482028185
  erange:
  - 0.38285463714552964
  - 0.3889339565295549
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.9315782971441362
    - 0.03100951943126902
    - -0.0056066667680175625
  dipole:
  - 0.9315782971441362
  - 0.03100951943126902
  - -0.0056066667680175625
  strength: 0.05589384164430929
  strengths:
  - 0.05589384164430929
  signif: 1.0
