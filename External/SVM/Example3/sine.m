Y = [
1e-05
0.0011892614
0.0023181608
0.0033953433
0.0044204449
0.0053927571
0.0063114912
0.0071763802
0.0079869842
0.0087427508
0.0094438593
0.010089584
0.010680256
0.011215696
0.01169618
0.012121766
0.012492828
0.012809865
0.013072854
0.013282873
0.013440063
0.013545674
0.013600136
0.013604433
0.013559587
0.013466676
0.01332655
0.013140837
0.01291068
0.012637533
0.012322326
0.011967227
0.011573481
0.011142636
0.01067633
0.010176786
0.0096452974
0.0090840592
0.0084950133
0.0078796654
0.0072408432
0.0065796546
0.0058989007
0.0052002844
0.0044860237
0.003758456
0.0030199188
0.0022719417
0.0015175893
0.00075857677
-2.65359e-06
-0.00076388395
-0.0015227379
-0.0022773416
-0.0030250405
-0.0037637632
-0.0044913309
-0.0052055916
-0.005904142
-0.0065848959
-0.0072459649
-0.0078849725
-0.0085003205
-0.0090893664
-0.0096506046
-0.010182027
-0.010681769
-0.01114785
-0.011578788
-0.011972534
-0.012327633
-0.012642682
-0.012915987
-0.013146145
-0.013331857
-0.013471984
-0.013564829
-0.01360974
-0.013605284
-0.013550981
-0.01344537
-0.01328818
-0.013078161
-0.012815014
-0.01249834
-0.012127185
-0.011701646
-0.011221003
-0.010685722
-0.010094891
-0.0094490079
-0.0087480774
-0.0079922914
-0.007181641
-0.0063167983
-0.0053979057
-0.0044257521
-0.0034006042
-0.0023233094
-0.0011947272
-1.530718e-05
];
X = 0:100;
X = .01 * X;
YY = sin(2 * pi * X);
subplot(2,2,1)
plot(X,Y,X,YY)
title('Y=sin(2\piX), n=5')
grid

Y = [
0.81442389
0.80884276
0.80273673
0.79610505
0.78894818
0.78126585
0.77305937
0.7643305
0.75508076
0.7453129
0.73503063
0.72423634
0.71293558
0.70113198
0.68883137
0.67603993
0.66276357
0.64900943
0.63478473
0.62009716
0.60495626
0.5893699
0.57334807
0.55690057
0.54003856
0.52277231
0.50511302
0.48707341
0.46866528
0.44990151
0.43079487
0.41135959
0.39160904
0.37155794
0.35122055
0.33061246
0.30974825
0.2886444
0.26731627
0.24577995
0.22405271
0.20214957
0.18008969
0.15788817
0.13556315
0.11313173
0.090612197
0.068020773
0.045376555
0.022695984
-1.9e-06
-0.022699784
-0.045379857
-0.068024684
-0.090615886
-0.11313553
-0.13556695
-0.15789197
-0.18009349
-0.20215337
-0.22405625
-0.24578375
-0.26732007
-0.2886482
-0.30975205
-0.33061622
-0.35122439
-0.37156159
-0.39161284
-0.41136339
-0.43079867
-0.44990516
-0.46866896
-0.4870771
-0.50511682
-0.52277611
-0.54004252
-0.55690437
-0.57335141
-0.58937374
-0.60496006
-0.62010096
-0.63478853
-0.64901304
-0.66276755
-0.67604426
-0.68883532
-0.70113578
-0.71293938
-0.72424029
-0.73503401
-0.74531689
-0.75508444
-0.76433414
-0.77306301
-0.78126934
-0.78895187
-0.79610916
-0.80274053
-0.80884671
-0.81442769
];
subplot(2,2,2)
plot(X,Y,X,YY)
title('Y=sin(2\piX), n=17')
grid

Y = [
0.81301598
0.81610409
0.81829966
0.81959428
0.81998407
0.81946249
0.81802725
0.81567544
0.81240639
0.80821897
0.80311464
0.79709612
0.79016703
0.78233079
0.77359504
0.76396553
0.75345199
0.74206243
0.72980829
0.71670081
0.70275331
0.68798073
0.67239766
0.65602045
0.63886645
0.62095476
0.60230467
0.58293738
0.56287472
0.54213891
0.52075266
0.49874276
0.47613292
0.45294963
0.42922072
0.40497444
0.38023781
0.35504185
0.32941598
0.30339012
0.27699669
0.25026563
0.22323181
0.19592556
0.16838176
0.14063236
0.11271247
0.084653723
0.056492799
0.028262763
-1.65e-06
-0.028266063
-0.056495894
-0.084657476
-0.11271532
-0.14063566
-0.16838506
-0.19592886
-0.22323481
-0.25026918
-0.27699974
-0.30339342
-0.32941928
-0.35504515
-0.38024111
-0.40497804
-0.42922477
-0.45295278
-0.47613622
-0.49874606
-0.52075596
-0.54214176
-0.56287762
-0.58294053
-0.60230797
-0.62095806
-0.63886945
-0.6560233
-0.67240045
-0.68798408
-0.70275727
-0.71670411
-0.72981159
-0.74206588
-0.75345521
-0.76396958
-0.77359849
-0.78233409
-0.79017041
-0.79709979
-0.80311789
-0.80822217
-0.81240961
-0.81567843
-0.81803025
-0.81946542
-0.81998684
-0.81959758
-0.81830231
-0.81610779
-0.81301928
];
subplot(2,2,3)
plot(X,Y,X,YY)
title('Y=sin(2\piX), n=65')
grid

Y = [
0.4327151
0.47412782
0.51337175
0.55040987
0.58521631
0.61776221
0.64802729
0.67599298
0.70164202
0.72496105
0.74594404
0.76458374
0.78088179
0.79483638
0.8064576
0.81575072
0.82273399
0.82742166
0.8298337
0.82999538
0.8279337
0.82368142
0.81727173
0.80874604
0.79814269
0.78550974
0.77089137
0.75434213
0.73591841
0.71567587
0.69367523
0.6699817
0.64465982
0.61778135
0.58941645
0.55963832
0.52852508
0.49615517
0.46261155
0.42797042
0.39232584
0.35575537
0.31835017
0.28020039
0.24139628
0.20202949
0.16219109
0.12197425
0.081473715
0.040783999
-1.060496e-06
-0.040784729
-0.081474709
-0.12197585
-0.16219186
-0.20203082
-0.24139803
-0.28020202
-0.31835143
-0.35575652
-0.39232638
-0.42797251
-0.46261258
-0.49615691
-0.52852698
-0.55964093
-0.58941733
-0.61778178
-0.64466112
-0.669983
-0.69367697
-0.71567773
-0.73591903
-0.75434349
-0.77089236
-0.78551062
-0.79814278
-0.80874733
-0.81727299
-0.8236829
-0.82793624
-0.82999728
-0.82983514
-0.82742295
-0.82273574
-0.8157539
-0.80645904
-0.79483726
-0.78088361
-0.76458632
-0.74594345
-0.72496268
-0.7016429
-0.67599383
-0.64802798
-0.61776343
-0.58521677
-0.55041078
-0.51337203
-0.47413039
-0.43271692
];
subplot(2,2,4)
plot(X,Y,X,YY)
title('Y=sin(2\piX), n=257')
grid

