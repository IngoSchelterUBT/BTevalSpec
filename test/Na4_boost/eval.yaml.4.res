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
    max_excit: 21                   # Maximum numer of iterations used to reach relative error between fit and 
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
  energy: 0.13092306419927804
  erange:
  - 0.12590098276091613
  - 0.1319803021449414
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - -0.0010526003407142497
    - 0.08391891256824972
    - 3.3317956823724666
  dipole:
  - -0.0010526003407142497
  - 0.08391891256824972
  - 3.3317956823724666
  strength: 0.24238018087372867
  strengths:
  - 0.24238018087372867
  signif: 0.07364994865724547
- name: S1
  fix: false
  energy: 0.1320825365004415
  erange:
  - 0.12866034030798737
  - 0.13473965969201265
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - -0.0004390253774127231
    - 0.0011014914860236523
    - 6.3027766689375495
  dipole:
  - -0.0004390253774127231
  - 0.0011014914860236523
  - 6.3027766689375495
  strength: 0.8744963535262137
  strengths:
  - 0.8744963535262137
  signif: 0.3435492899797113
- name: S2
  fix: false
  energy: 0.15420785044447738
  erange:
  - 0.15272034030798737
  - 0.15879965969201265
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.0022161324612049774
    - -0.019993318238723457
    - 1.6724561366133057
  dipole:
  - 0.0022161324612049774
  - -0.019993318238723457
  - 1.6724561366133057
  strength: 0.07189977453214008
  strengths:
  - 0.07189977453214008
  signif: 1.0
- name: S3
  fix: false
  energy: 0.1647563904054251
  erange:
  - 0.1618466654591247
  - 0.16792598484315
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 2.0816891645804456
    - 0.04543669004007996
    - -0.0032665813869185594
  dipole:
  - 2.0816891645804456
  - 0.04543669004007996
  - -0.0032665813869185594
  strength: 0.11905035745189914
  strengths:
  - 0.11905035745189914
  signif: 1.0
- name: S4
  fix: false
  energy: 0.19340504902111144
  erange:
  - 0.19024034030798737
  - 0.19631965969201265
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.00015323693508719518
    - 5.314732956574166
    - -0.0012169180389210378
  dipole:
  - 0.00015323693508719518
  - 5.314732956574166
  - -0.0012169180389210378
  strength: 0.910499006209223
  strengths:
  - 0.910499006209223
  signif: 1.0
- name: S5
  fix: false
  energy: 0.20559014140030235
  erange:
  - 0.2021985931332427
  - 0.208277912517268
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.014109395682623286
    - -0.13992293709261241
    - 0.9226387292289961
  dipole:
  - 0.014109395682623286
  - -0.13992293709261241
  - 0.9226387292289961
  strength: 0.029846196808989856
  strengths:
  - 0.029846196808989856
  signif: 1.0
- name: S6
  fix: false
  energy: 0.21822688233661278
  erange:
  - 0.2148503403079874
  - 0.22092965969201267
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 4.087380257015956
    - -0.02317603211204372
    - 0.003795432250568683
  dipole:
  - 4.087380257015956
  - -0.02317603211204372
  - 0.003795432250568683
  strength: 0.6076610791932807
  strengths:
  - 0.6076610791932807
  signif: 0.9934953916037415
- name: S7
  fix: false
  energy: 0.2242796212038212
  erange:
  - 0.2214469264490576
  - 0.2275262458330829
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.029441906859889014
    - -0.05828300525905339
    - 1.09426443968159
  dipole:
  - 0.029441906859889014
  - -0.05828300525905339
  - 1.09426443968159
  strength: 0.04491866267197483
  strengths:
  - 0.04491866267197483
  signif: 1.0006614100157079
- name: S8
  fix: false
  energy: 0.24264390006855244
  erange:
  - 0.2411903403079874
  - 0.24726965969201267
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.0180314612583441
    - 2.1997970940566463
    - 0.0009839587186468639
  dipole:
  - 0.0180314612583441
  - 2.1997970940566463
  - 0.0009839587186468639
  strength: 0.19570983063555947
  strengths:
  - 0.19570983063555947
  signif: 1.0
- name: S9
  fix: false
  energy: 0.2550581685114844
  erange:
  - 0.2526703403079874
  - 0.2587496596920127
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 2.6931458406417694
    - -0.009924612794117667
    - 0.03273536443202829
  dipole:
  - 2.6931458406417694
  - -0.009924612794117667
  - 0.03273536443202829
  strength: 0.30837402410619547
  strengths:
  - 0.30837402410619547
  signif: 1.0
- name: S10
  fix: false
  energy: 0.26551097181155864
  erange:
  - 0.2617988541231756
  - 0.2678781735072009
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.011210061624636345
    - 0.16564247075879732
    - 0.265255953462534
  dipole:
  - 0.011210061624636345
  - 0.16564247075879732
  - 0.265255953462534
  strength: 0.004333307856264985
  strengths:
  - 0.004333307856264985
  signif: 1.0
- name: S11
  fix: false
  energy: 0.3012553842361737
  erange:
  - 0.30122315127604954
  - 0.3073024706600748
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 1.0795265039382713
    - -0.07279766989665061
    - -0.01813721222621833
  dipole:
  - 1.0795265039382713
  - -0.07279766989665061
  - -0.01813721222621833
  strength: 0.05879530699364401
  strengths:
  - 0.05879530699364401
  signif: 0.22986392970828828
- name: S12
  fix: false
  energy: 0.302789568758025
  erange:
  - 0.293
  - 0.313
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.15093887245483145
    - 0.6646869633661225
    - 0.020979607326462528
  dipole:
  - 0.15093887245483145
  - 0.6646869633661225
  - 0.020979607326462528
  strength: 0.02346777851958688
  strengths:
  - 0.02346777851958688
  signif: 0.49035359616962737
- name: S13
  fix: false
  energy: 0.3114942162184934
  erange:
  - 0.30864419544600225
  - 0.31472351483002753
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - -0.059830741687682154
    - -0.07845309535940333
    - 0.5717997570886518
  dipole:
  - -0.059830741687682154
  - -0.07845309535940333
  - 0.5717997570886518
  strength: 0.01747947551702331
  strengths:
  - 0.01747947551702331
  signif: 1.002074599642955
- name: S14
  fix: false
  energy: 0.31445782072962686
  erange:
  - 0.3130504404219117
  - 0.319129759805937
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.5342315345778449
    - 0.05891814148573598
    - 0.023984635909252663
  dipole:
  - 0.5342315345778449
  - 0.05891814148573598
  - 0.023984635909252663
  strength: 0.015169966364938168
  strengths:
  - 0.015169966364938168
  signif: 0.9957556188314319
- name: S15
  fix: false
  energy: 0.3182066739135317
  erange:
  - 0.307
  - 0.327
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.0039949647113622245
    - 0.7848146467292302
    - -0.019958575273895204
  dipole:
  - 0.0039949647113622245
  - 0.7848146467292302
  - -0.019958575273895204
  strength: 0.032687692232496623
  strengths:
  - 0.032687692232496623
  signif: 0.9713263704214933
- name: S16
  fix: false
  energy: 0.33354282332307117
  erange:
  - 0.3304435126952384
  - 0.3365228320792637
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.9319695689039267
    - 0.011452432389605355
    - 0.010661474272278768
  dipole:
  - 0.9319695689039267
  - 0.011452432389605355
  - 0.010661474272278768
  strength: 0.04829767360664688
  strengths:
  - 0.04829767360664688
  signif: 0.9569853403061291
- name: S17
  fix: false
  energy: 0.33695470484499773
  erange:
  - 0.33392212714990377
  - 0.34000144653392905
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.09068948715727002
    - 0.5783530409785912
    - -0.029320348109462952
  dipole:
  - 0.09068948715727002
  - 0.5783530409785912
  - -0.029320348109462952
  strength: 0.019294953356386604
  strengths:
  - 0.019294953356386604
  signif: 0.9863079533782362
- name: S18
  fix: false
  energy: 0.3449252550846697
  erange:
  - 0.3418069865804786
  - 0.3478863059645039
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.7756191454439738
    - 0.0002176625789278859
    - 0.03905566999987315
  dipole:
  - 0.7756191454439738
  - 0.0002176625789278859
  - 0.03905566999987315
  strength: 0.03467133772219531
  strengths:
  - 0.03467133772219531
  signif: 1.0
- name: S19
  fix: false
  energy: 0.37443895024164087
  erange:
  - 0.3719549785209116
  - 0.37803429790493687
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.00688366015110618
    - 0.7182345848490089
    - 0.007292226014298243
  dipole:
  - 0.00688366015110618
  - 0.7182345848490089
  - 0.007292226014298243
  strength: 0.0321993458358601
  strengths:
  - 0.0321993458358601
  signif: 1.0
- name: S20
  fix: false
  energy: 0.3860032559767218
  erange:
  - 0.38285463714552964
  - 0.3889339565295549
  phase: 0.0
  tmod: 1.0
  dipoles:
  - - 0.9310502341240384
    - 0.02866380090378274
    - -0.004825505390160169
  dipole:
  - 0.9310502341240384
  - 0.02866380090378274
  - -0.004825505390160169
  strength: 0.05582246801061536
  strengths:
  - 0.05582246801061536
  signif: 1.0
