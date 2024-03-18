DESCRIPTION:                        # Description of evaluation 
DIPOLE:                             # List of dipole moment files; different calculations are only allowed to differ in the boost direction or laser polarization
- [dipole_atomgroups1_01.dat, dipole_atomgroups1_02.dat]
EXT:
  profile: laser_profile.dat      # Profile of excitation
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
    smooth: 0.002475691517674122    # Smooth for Pade Approximation (0:auto)
    thin: 0                         # Only keep every 2^n data point
  Fit:
    calc: true                      # Turn fit on/off
    skipfirst: true                 # Skip first fit of existing excitations
    firstsingle: false              # If not skipfirst: Fit input excitations all at once (False) or one after the other (by strength) (True)
    guess: false                    # Turn guess for fit via Pade Approximation on/off
    plot_result: true               # If True: The fit results are plotted without fitting again
    gnuplot_spectrum: false         # If True: A gnuplot script for plotting the resulting spectrum is created.
    dat_spectrum: false             # If True: A .dat file for the spectrum is created
    guess_thres: 0.1                # Relative height of line in Pade Approximation compared to highest line 
                                    #   which should be identified as a line for fitting (only relevent, if 
                                    #   fit_guess == True).
    relerr_crit: 0.001               # Criterium for relative error between fit and data (only relevant, if 
                                    #   fit_guess == True).
    max_excit: 16                   # Maximum numer of iterations used to reach relative error between fit and 
                                    #   data (only relevant, if fit_guess == True)
    relspacing_lines: 0.01          # Threshold for relative error between two lines which should be identified 
                                    #   as one in fit of trace (only relevant, if number of dipole files == 3).
    range:                          # Range of spectrum which should be fitted in Ry
    - 0.10
    - 0.40
    nsigma: 2.                      # New excitations are guessed at energies where the peak height is larger than
                                    # the mean value plus nsigma times the standard deviation of the distribution of
                                    # maxima of Fourier Transform minus current fit
                                    # (and the latter scaled for error compensation around existing lines)
    significances: true
    fiterr: 0.0012123154951597852   # Current fit error
    reset_erange: false             # Reset the excitation-specific energy fit range around the current energy
SPEC:
- name: S0
  fix: false
  energy: 0.10143135326089937
  energyErr: 5.4090065718404174e-05
  erange:
  - 0.1001185622234784
  - 0.1034963371895998
  phase: -10.015986580448212
  phaseErr: 0.05554630340655411
  tmod: 1.0
  dipoles:
  - - -0.00010636543046231801
    - 0.013036719606364891
    - 0.11227814932882531
  - - -0.00010535705353016296
    - -0.031918371806817676
    - -0.1367058608500365
  dipolesErr:
  - - 0.004162780355815126
    - 0.004452298983101061
    - 0.016755376033465264
  - - 0.004162787438735439
    - 0.005670670739936351
    - 0.013867922367603716
  ampl:
  - - - 1.1247909313304102e-09
      - -1.3786043006455412e-07
      - -1.187316627241637e-06
    - - 1.1141275680203486e-09
      - 3.375297311832866e-07
      - 1.4456342805693095e-06
  dipole:
  - -0.00021172248399248098
  - -0.018881652200452785
  - -0.02442771152121119
  dipoleErr:
  - 0.005887065444757164
  - 0.007209679110452352
  - 0.02174998613371295
  strength: 1.6115323916464747e-05
  strengthErr: 2.2608354341587702e-05
  strengths:
  - 0.00021598708155752403
  - 0.0003331561138288615
  strengthsErr:
  - 6.555384076258748e-05
  - 7.023321562107444e-05
  signifFit: 1.0
  signifErr: 0.19002301669506816
  signifAng: 0.8894814058471977
  signifExc: -0.21048246469222032
  signifRng: 0.9975407988780665
- name: S1
  fix: false
  energy: 0.1490962938624608
  energyErr: 2.167148396737459e-08
  erange:
  - 0.1474311125169393
  - 0.1508088874830607
  phase: 4.977784228945919
  phaseErr: 2.3172100562891603e-05
  tmod: 1.0
  dipoles:
  - - -3.9112125782732427e-07
    - 0.03776523684341318
    - 1.9892493923598278
  - - -4.753890764770753e-07
    - -0.4218937079897963
    - -1.317685520958409
  dipolesErr:
  - - 2.749693146007809e-05
    - 2.7550826155428196e-05
    - 4.301587418386984e-05
  - - 2.74969430815025e-05
    - 3.0140445386990956e-05
    - 6.11005672213552e-05
  ampl:
  - - - -5.760356459152537e-10
      - 5.561989322462819e-05
      - 0.0029297271260065646
    - - -7.001436210619517e-10
      - -0.0006213567013449751
      - -0.0019406611504444481
  dipole:
  - -8.665103343043996e-07
  - -0.38412847114638315
  - 0.6715638714014187
  dipoleErr:
  - 3.888654161207474e-05
  - 4.0834966266336146e-05
  - 7.472378969627988e-05
  strength: 0.01487366526580261
  strengthErr: 1.7144031794769844e-06
  strengths:
  - 0.09836725790067813
  - 0.04756890997302348
  strengthsErr:
  - 4.304398320825358e-06
  - 4.6332872247797815e-06
  signifFit: 1.0482207634887983
  signifErr: 0.9999333978539211
  signifAng: 0.9316827628134388
  signifExc: 0.9999757782454244
  signifRng: 0.9999999611813586
- name: S2
  fix: false
  energy: 0.1593670344291849
  energyErr: 3.8089947550088076e-08
  erange:
  - 0.1575411125169393
  - 0.1609188874830607
  phase: 3.9148106369894715
  phaseErr: 4.1915630415614765e-05
  tmod: 1.0
  dipoles:
  - - -2.169135624675573e-06
    - -0.03210521155658637
    - 0.5002150954592623
  - - -6.376283524195027e-07
    - 0.26524502955648227
    - 0.6625473512108062
  dipolesErr:
  - - 1.701385220899459e-05
    - 1.710753103729557e-05
    - 1.4009992487797834e-05
  - - 1.7013878830976126e-05
    - 1.7262864086941507e-05
    - 1.321654751895819e-05
  ampl:
  - - - -5.164649806391532e-09
      - -7.644158934261323e-05
      - 0.0011909978179922211
    - - -1.5181748478109445e-09
      - 0.0006315408197447734
      - 0.0015775062703457735
  dipole:
  - -2.8067639770950757e-06
  - 0.23313981799989592
  - 1.1627624466700686
  dipoleErr:
  - 2.406123936676263e-05
  - 2.4303787661108423e-05
  - 1.9260244490375287e-05
  strength: 0.03735485415322966
  strengthErr: 1.4906775330407123e-06
  strengths:
  - 0.006673385303066195
  - 0.01352824040392127
  strengthsErr:
  - 3.431030783072889e-07
  - 7.084114421870153e-07
  signifFit: 1.0
  signifErr: 0.9999761454682211
  signifAng: 0.9901945860165712
  signifExc: 0.9999805172116067
  signifRng: 0.9999566573049415
- name: S3
  fix: false
  energy: 0.1642358211252503
  energyErr: 2.914712363303095e-06
  erange:
  - 0.1625017147771436
  - 0.165879489743265
  phase: 3.4672340763087863
  phaseErr: 0.0030636927268759273
  tmod: 1.0
  dipoles:
  - - 4.8393070924920765e-06
    - -0.06630573464684497
    - -0.1224162289237051
  - - -2.648144770090775e-06
    - -0.0073937165881612825
    - -0.08117829671709643
  dipolesErr:
  - - 0.00011772427203178719
    - 0.00012186561119403334
    - 9.348814082436637e-05
  - - 0.00011772459916849398
    - 0.0001226932233018002
    - 0.00017036929504379126
  ampl:
  - - - -1.7379336397909554e-09
      - 2.3812286459891203e-05
      - 4.396317039539188e-05
    - - 9.51024555998412e-10
      - 2.6552951767788524e-06
      - 2.9153449034977638e-05
  dipole:
  - 2.1911623224013015e-06
  - -0.07369945123500625
  - -0.20359452564080155
  dipoleErr:
  - 0.00016648749344859754
  - 0.0001729302004737192
  - 0.00019433406589817174
  strength: 0.0012832932673682948
  strengthErr: 2.8637193912656587e-06
  strengths:
  - 0.000530541472423043
  - 0.00018187975504003443
  strengthsErr:
  - 1.068863435221626e-06
  - 8.068228713659632e-07
  signifFit: 0.9627419892636032
  signifErr: 0.9986142944637787
  signifAng: 0.9696851128579872
  signifExc: -0.15521868694413565
  signifRng: 0.999999486103378
- name: S4
  fix: false
  energy: 0.16680508950831596
  energyErr: 1.0240306457527805e-05
  erange:
  - 0.16598032923180894
  - 0.16935810419793035
  phase: 3.0321388110402947
  phaseErr: 0.010743917595650602
  tmod: 1.0
  dipoles:
  - - -5.900980010175687e-06
    - 0.017938155188275597
    - 0.03152209894997998
  - - 1.1223467561524798e-05
    - -0.015063337308925926
    - 0.08089083669812473
  dipolesErr:
  - - 0.00023782092997140508
    - 0.0002713445557561979
    - 0.00031899261958484685
  - - 0.00023782186908864302
    - 0.00024513063748264814
    - 0.00018159369715758787
  ampl:
  - - - -1.0433112701229648e-09
      - 3.171520567917577e-06
      - 5.573203270598141e-06
    - - 1.9843433085023633e-09
      - -2.663243995567823e-06
      - 1.430174673211903e-05
  dipole:
  - 5.322487551349111e-06
  - 0.002874817879349671
  - 0.11241293564810471
  dipoleErr:
  - 0.000336330248638566
  - 0.00036567321117519377
  - 0.00036705934424417606
  strength: 0.0003515398553416665
  strengthErr: 2.3528025777256947e-06
  strengths:
  - 3.6569796493854136e-05
  - 0.00018821818463472332
  strengthsErr:
  - 8.296515852908817e-07
  - 6.115911071914854e-07
  signifFit: 0.9812694379402457
  signifErr: 0.9956277577206067
  signifAng: 0.9998365623456156
  signifExc: 0.9909827946794583
  signifRng: 0.931465721121619
- name: S5
  fix: false
  energy: 0.19832267244656362
  energyErr: 4.152482208838663e-05
  erange:
  - 0.19659213643286394
  - 0.19996991139898534
  phase: 2.9989673733306743
  phaseErr: 0.04467222515881847
  tmod: 1.0
  dipoles:
  - - 8.299859915605709e-05
    - -0.10119994316423325
    - -0.010628359177855578
  - - 4.99298025868512e-05
    - 0.07143886099577394
    - -0.02161941868842341
  dipolesErr:
  - - 0.0027952528814365364
    - 0.006814331125594299
    - 0.0023827431260076896
  - - 0.0027952510452290266
    - 0.005201414350642358
    - 0.0020838392431064357
  ampl:
  - - - -1.1992002424228418e-09
      - 1.4621812610058946e-06
      - 1.5356320506900143e-07
    - - -7.214077342883319e-10
      - -1.032180064430543e-06
      - 3.123668639690076e-07
  dipole:
  - 0.0001329284017429083
  - -0.02976108216845931
  - -0.032247777866278984
  dipoleErr:
  - 0.0039530832367954535
  - 0.00857262036581065
  - 0.0031654147904575315
  strength: 6.365030185987077e-05
  strengthErr: 2.357941398356466e-05
  strengths:
  - 0.00034225196079000014
  - 0.00018413971305420235
  strengthsErr:
  - 4.72472433869593e-05
  - 2.1595442437373812e-05
  signifFit: 0.9999708194487542
  signifErr: 0.7859463090066549
  signifAng: 0.8572451936085812
  signifExc: -0.09454403991430405
  signifRng: 0.9999996301753385
- name: S6
  fix: false
  energy: 0.20562884265313605
  energyErr: 4.831321811483922e-06
  erange:
  - 0.2040131806028167
  - 0.2073909555689381
  phase: -10.30435614387754
  phaseErr: 0.005195734396578725
  tmod: 1.0
  dipoles:
  - - 0.00010273572677940876
    - -0.22758620584294653
    - 0.008687111340683562
  - - -4.345186657827569e-05
    - -0.29125368245688754
    - 0.09487962977233119
  dipolesErr:
  - - 0.0009805612313576008
    - 0.0018171619496880277
    - 0.0009414557976293227
  - - 0.0009805614103695853
    - 0.0021854059535533032
    - 0.0006970001109383345
  ampl:
  - - - 4.230106530061226e-09
      - -9.370780016529386e-06
      - 3.57688679114496e-07
    - - -1.789114948792607e-09
      - -1.1992265423990698e-05
      - 3.906634566682914e-06
  dipole:
  - 5.928386020113307e-05
  - -0.5188398882998341
  - 0.10356674111301475
  dipoleErr:
  - 0.0013867231187037754
  - 0.002842195759130644
  - 0.0011713872858870864
  strength: 0.00959330188879546
  strengthErr: 9.275545075611486e-05
  strengths:
  - 0.0017776941613196918
  - 0.003215721637750362
  strengthsErr:
  - 2.777920524663089e-05
  - 3.909820449417921e-05
  signifFit: 1.0
  signifErr: 0.9944101460662516
  signifAng: 0.4424368816601379
  signifExc: 0.9917542766399631
  signifRng: 0.9999964661759696
- name: S7
  fix: false
  energy: 0.22681655797116607
  energyErr: 1.676593459268737e-05
  erange:
  - 0.2251167749611198
  - 0.2284945499272412
  phase: 3.5797911258769712
  phaseErr: 0.017996943376191638
  tmod: 1.0
  dipoles:
  - - -0.0018327794331229256
    - -0.7187581530589874
    - -0.999472089606432
  - - -0.00037005404721541694
    - -0.21031039260778658
    - 0.9831246553881887
  dipolesErr:
  - - 0.014312616423048857
    - 0.44430389996954733
    - 0.6078675845692513
  - - 0.01426964763117734
    - 0.13072266966436236
    - 0.617966183680933
  ampl:
  - - - 5.191125053761056e-09
      - 2.0357951363421274e-06
      - 2.8308832536658416e-06
    - - 1.0481331255842312e-09
      - 5.956786334470946e-07
      - -2.7845811325260186e-06
  dipole:
  - -0.0022028334803383424
  - -0.929068545666774
  - -0.01634743421824325
  dipoleErr:
  - 0.02021073557274185
  - 0.463135371022693
  - 0.866824783069354
  strength: 0.03264043195083646
  strengthErr: 0.033606607043623296
  strengths:
  - 0.057292391999965886
  - 0.03820969450753019
  strengthsErr:
  - 0.07008026551047482
  - 0.0438542386529994
  signifFit: 0.9999974871403478
  signifErr: 0.40555275271391233
  signifAng: 0.13263761410872615
  signifExc: 0.12982016634333093
  signifRng: 0.9999999982678317
- name: S8
  fix: false
  energy: 0.27047085069982385
  energyErr: 4.154785586863437e-05
  erange:
  - 0.26871540945959216
  - 0.2720931844257135
  phase: -7.570495028382724
  phaseErr: 0.044546515300976905
  tmod: 1.0
  dipoles:
  - - 0.03390530114401514
    - 17.636657957640633
    - -1.7643644822716864
  - - -0.002016206415176394
    - -5.133799043160107
    - 1.7650576981141908
  dipolesErr:
  - - 9.18356450986085
    - 4772.327573866873
    - 465.75519025237395
  - - 0.681567783174598
    - 1389.162342008962
    - 465.5676120513668
  ampl:
  - - - 3.3511718194776023e-09
      - 1.7431926319239746e-06
      - -1.7438832078681463e-07
    - - -1.9928016837513784e-10
      - -5.074204357372614e-07
      - 1.744568376652328e-07
  dipole:
  - 0.031889094728838746
  - 12.502858914480527
  - 0.0006932158425043689
  dipoleErr:
  - 9.208821409379013
  - 4970.400636241096
  - 658.5446823399586
  strength: 7.0467798556615735
  strengthErr: 5602.800794252341
  strengths:
  - 14.162116798560172
  - 1.3285224068189945
  strengthsErr:
  - 7514.262514576361
  - 568.8838627545498
  signifFit: 1.0
  signifErr: -458.04350825818915
  signifAng: 0.007446100057637386
  signifExc: 0.6324112477739964
  signifRng: 0.9999975884990173
- name: S9
  fix: false
  energy: 0.2949434728972713
  energyErr: 9.935844771311808e-06
  erange:
  - 0.2935295259028716
  - 0.296907300868993
  phase: 2.4406448677316255
  phaseErr: 0.01076439863088217
  tmod: 1.0
  dipoles:
  - - -5.4454268186397735e-05
    - 0.023732477164642475
    - 0.49422823296323154
  - - 0.0005907413363250225
    - -0.039249719852384214
    - -0.18084716360778053
  dipolesErr:
  - - 0.002009436503557412
    - 0.002058113810566853
    - 0.002250283798035948
  - - 0.002009752757210849
    - 0.002176946456322401
    - 0.0038306442634711948
  ampl:
  - - - -1.1490893892024515e-09
      - 5.008007378968354e-07
      - 1.042916262134594e-05
    - - 1.2465774014458667e-08
      - -8.282442885314454e-07
      - -3.816221642309566e-06
  dipole:
  - 0.0005362870681386248
  - -0.015517242687741739
  - 0.313381069355451
  dipoleErr:
  - 0.00284199598995946
  - 0.0029958184742972446
  - 0.004442703303729371
  strength: 0.004839470209600302
  strengthErr: 0.0001324587680345907
  strengths:
  - 0.012034911795788064
  - 0.0016834678383343543
  strengthsErr:
  - 0.0001141322247619193
  - 7.639217465888718e-05
  signifFit: 1.000004251538465
  signifErr: 0.9839937814792561
  signifAng: 0.999387259696014
  signifExc: 0.9840885008712201
  signifRng: 0.9992976526746816
- name: S10
  fix: false
  energy: 0.29739667571137596
  energyErr: 2.9280077767420637e-05
  erange:
  - 0.2965443250969149
  - 0.29992210006303627
  phase: 2.1235017213233824
  phaseErr: 0.03133275924428927
  tmod: 1.0
  dipoles:
  - - 0.0008565568063924746
    - -0.008992978802807094
    - 0.04058516916773636
  - - -0.000696551308709328
    - 0.06323049750820645
    - 0.19716330754691522
  dipolesErr:
  - - 0.002252286895529536
    - 0.002304162449495129
    - 0.0029984155126281924
  - - 0.002252338544610786
    - 0.002328316297042037
    - 0.002227478734543008
  ampl:
  - - - 1.6254313820981845e-08
      - -1.706538300266414e-07
      - 7.701579991038632e-07
    - - -1.3217994976726967e-08
      - 1.1998834658541893e-06
      - 3.7414381053693935e-06
  dipole:
  - 0.0001600054976831466
  - 0.05423751870539935
  - 0.23774847671465157
  dipoleErr:
  - 0.0031852511956380314
  - 0.0032757016611307775
  - 0.0037352586388643434
  strength: 0.0029475020325020617
  strengthErr: 0.0001056975537314534
  strengths:
  - 8.568807709690531e-05
  - 0.002124995943597206
  strengthsErr:
  - 1.0200619883684531e-05
  - 5.797541403671461e-05
  signifFit: 0.9999993876886754
  signifErr: 0.9776120595183628
  signifAng: 0.9873964422065495
  signifExc: 0.9827689348724088
  signifRng: 0.939808179695204
- name: S11
  fix: false
  energy: 0.30540535904288285
  energyErr: 3.993448695650034e-05
  erange:
  - 0.30535681504873374
  - 0.3087345900148551
  phase: 0.7546167546645042
  phaseErr: 0.050190220558370274
  tmod: 1.0
  dipoles:
  - - -6.087381357104584e-05
    - 0.14022781338093412
    - 0.13260111923306542
  - - -0.005111386220798187
    - 0.04636451601111122
    - 0.07263031151402627
  dipolesErr:
  - - 0.0023378677554415418
    - 0.004965051137478188
    - 0.0018033301533832702
  - - 0.0023504652676070515
    - 0.0025271467720468766
    - 0.0024815086307418167
  ampl:
  - - - -1.0787375275826444e-09
      - 2.4849602124618536e-06
      - 2.3498084829074213e-06
    - - -9.057826035342149e-08
      - 8.216200108938363e-07
      - 1.2870730134026642e-06
  dipole:
  - -0.005172260034369233
  - 0.18659232939204534
  - 0.2052314307470917
  dipoleErr:
  - 0.003315164040611016
  - 0.005571194091506811
  - 0.0030675535409422695
  strength: 0.003917505283945445
  strengthErr: 0.00016817189966593406
  strengths:
  - 0.0018959004867656973
  - 0.0003792603543210331
  strengthsErr:
  - 9.520705052371936e-05
  - 2.9053067740743002e-05
  signifFit: 0.9999803774314732
  signifErr: 0.9542961100073195
  signifAng: 0.8601038456021405
  signifExc: 0.9258368696552521
  signifRng: 0.11010999609658134
- name: S12
  fix: false
  energy: 0.3072758241994766
  energyErr: 3.176569299937656e-05
  erange:
  - 0.3067482608305999
  - 0.31012603579672127
  phase: 0.6504751761334361
  phaseErr: 0.051252989588481054
  tmod: 1.0
  dipoles:
  - - -0.0028033708144496183
    - -0.2489689039883877
    - 0.10033466200703978
  - - 0.010645094839374605
    - -0.0188122533113839
    - 0.09723409910845145
  dipolesErr:
  - - 0.002891277696819195
    - 0.010805508719046853
    - 0.003323858680533693
  - - 0.0028962930215688095
    - 0.004194687158569785
    - 0.004607987452831827
  ampl:
  - - - -4.478947428381526e-08
      - -3.977777847004824e-06
      - 1.6030475670846607e-06
    - - 1.7007675156615004e-07
      - -3.005634971898187e-07
      - 1.5535098528815048e-06
  dipole:
  - 0.007841724024924987
  - -0.2677811572997716
  - 0.1975687611154912
  dipoleErr:
  - 0.004092432038154365
  - 0.011591135364392387
  - 0.0056816885600686285
  strength: 0.005674445039058067
  strengthErr: 0.0001996547758713513
  strengths:
  - 0.0036904040781115723
  - 0.0005081158581693491
  strengthsErr:
  - 0.00024221969664087824
  - 4.0967378849354587e-05
  signifFit: 0.9999893254650561
  signifErr: 0.9501791945639977
  signifAng: 0.7704111777498245
  signifExc: 0.9652834857363117
  signifRng: 0.7764313481495492
- name: S13
  fix: false
  energy: 0.30883543015093734
  energyErr: 6.92179928332086e-05
  erange:
  - 0.3088354295033991
  - 0.31221320446952044
  phase: 3.0051608124310527
  phaseErr: 0.08758410091417487
  tmod: 1.0
  dipoles:
  - - 0.010487748039009613
    - 0.021469638356913116
    - 0.0639041524868629
  - - 0.021397887149772626
    - 0.047199023463382964
    - 0.14248862269715204
  dipolesErr:
  - - 0.0029383022641794863
    - 0.005923227868048179
    - 0.004133391346744704
  - - 0.0029435372406884275
    - 0.003075492172214038
    - 0.0026113323599925697
  ampl:
  - - - 1.6156446392106714e-07
      - 3.307412228832883e-07
      - 9.84447767096159e-07
    - - 3.296358907114794e-07
      - 7.271041309435817e-07
      - 2.195046816083743e-06
  dipole:
  - 0.031885635188782235
  - 0.06866866182029607
  - 0.20639277518401494
  dipoleErr:
  - 0.004159090247037438
  - 0.006674075252660281
  - 0.004889169737254819
  strength: 0.002487672297846911
  strengthErr: 0.00016471276838743914
  strengths:
  - 0.00023958829519298708
  - 0.0011832836478171376
  strengthsErr:
  - 4.345588703722988e-05
  - 5.973200484860908e-05
  signifFit: 0.999991560910632
  signifErr: 0.9582332594358227
  signifAng: 0.9689311532157446
  signifExc: 0.918140038754178
  signifRng: 1.5336436071100579e-06
- name: S14
  fix: false
  energy: 0.3299565949687838
  energyErr: 3.286616771963219e-05
  erange:
  - 0.3283156704495251
  - 0.33169344541564644
  phase: 2.0170623335675804
  phaseErr: 0.03525936880348554
  tmod: 1.0
  dipoles:
  - - 0.0005428452766400157
    - 0.10774046485000942
    - 0.23201420625465485
  - - 0.0014211687509638805
    - -0.005929886529460539
    - -0.03395685936701998
  dipolesErr:
  - - 0.004489494965532028
    - 0.004821276706233972
    - 0.003238092754074228
  - - 0.004489700678172634
    - 0.00449247379527366
    - 0.004904462647223339
  ampl:
  - - - 4.8824153653038765e-09
      - 9.690306311673115e-07
      - 2.086763529744612e-06
    - - 1.2782161777926174e-08
      - -5.33341089106428e-08
      - -3.0541205582036366e-07
  dipole:
  - 0.001964014027603896
  - 0.10181057832054888
  - 0.19805734688763488
  dipoleErr:
  - 0.00634925013093052
  - 0.006589918806730081
  - 0.0058769889009591455
  strength: 0.0027274192268612516
  strengthErr: 0.00020318446328574242
  strengths:
  - 0.003598665793100064
  - 6.54552260273947e-05
  strengthsErr:
  - 0.00014003002295121935
  - 2.0545229886114594e-05
  signifFit: 1.0
  signifErr: 0.9558209398667361
  signifAng: 0.9430482094404874
  signifExc: 0.9313021476215373
  signifRng: 0.9999993495391244
- name: S15
  fix: false
  energy: 0.3755726450827728
  energyErr: 0.00011444590134954053
  erange:
  - 0.3730738430995525
  - 0.37645161806567384
  phase: 5.477410219555111
  phaseErr: 0.12267304459893122
  tmod: 1.0
  dipoles:
  - - 0.0008279177253497121
    - 0.06380213350030879
    - 0.19168455427607098
  - - -6.364036271260373e-05
    - 0.00020569569102758678
    - -0.0726208950328979
  dipolesErr:
  - - 0.013036941878196049
    - 0.013946906863404513
    - 0.010819016571175575
  - - 0.013036802558864018
    - 0.013037775135831563
    - 0.017480924702004866
  ampl:
  - - - 2.568315326205775e-09
      - 1.979230451241443e-07
      - 5.946320068654235e-07
    - - -1.9742120975974285e-10
      - 6.380964915052262e-10
      - -2.2528006347126984e-07
  dipole:
  - 0.0007642773626371083
  - 0.06400782919133638
  - 0.11906365924317308
  dipoleErr:
  - 0.018436921502633385
  - 0.019091877632883246
  - 0.020558060414458745
  strength: 0.001143852890725074
  strengthErr: 0.0004611830039190927
  strengths:
  - 0.00255479323485485
  - 0.0003301183870914089
  strengthsErr:
  - 0.0003723771747361828
  - 0.00015869550284284883
  signifFit: 0.9999932173271324
  signifErr: 0.7668628445947987
  signifAng: 0.9384966894768424
  signifExc: 0.8581019151649543
  signifRng: 0.9471124036661455
