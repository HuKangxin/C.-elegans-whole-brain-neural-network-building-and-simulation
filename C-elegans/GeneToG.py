# Created by Kangxin Hu at Beihang University, on Oct.2023.

def GeneToG(shl1, shk1, kvs1, kqt3, egl2, egl36, irk, egl19, unc2, cca1, slo1, slo2, kcnl, nca, itr1, mec10, osm9, trp4):
    gshl1 = -0.01126 * shl1 + 4.589
    gshk1 = -0.001643 * shk1 + 1.108
    gkvs1 = -0.0009522 * kvs1 + 0.8884
    gkqt3 = 0.0104 * kqt3
    gegl2 = 0.04258 * egl2
    gegl36 = 0.05598 * egl36
    girk = -0.0005724 * irk + 0.2853
    gegl19 = -0.006216 * egl19 + 1.782
    gunc2 = -0.0007945 * unc2 + 1.026
    gcca1 = -0.8816 * cca1 + 10.12  # correct
    gslo1 = 0.005892 * slo1 - 0.02507
    gslo2 = 0.001368 * slo2 + 0.06728
    gkcnl = 0.06
    gnca = -8.143e-05 * nca + 0.06038
    gitr1 = 0.005 * itr1
    gmec10 = 0.0003 * mec10
    gosm9 = 0.005 * osm9
    gtrp4 = 0.005 * trp4

    if gshl1 < 0:
        gshl1 = 0
    if gshk1 < 0:
        gshk1 = 0
    if gkvs1 < 0:
        gkvs1 = 0
    if gkqt3 < 0:
        gkqt3 = 0
    if gegl2 < 0:
        gegl2 = 0
    if gegl36 < 0:
        gegl36 = 0
    if girk < 0:
        girk = 0
    if gegl19 < 0:
        gegl19 = 0
    if gunc2 < 0:
        gunc2 = 0
    if gcca1 < 0:
        gcca1 = 0
    if gslo1 < 0:
        gslo1 = 0
    if gslo2 < 0:
        gslo2 = 0
    if gkcnl < 0:
        gkcnl = 0
    if gnca < 0:
        gnca = 0
    if gtrp4 < 0:
        gtrp4 = 0
    if gitr1 < 0:
        gitr1 = 0
    if gosm9 < 0:
        gosm9 = 0
    if gmec10 < 0:
        gmec10 = 0
    thre = 10
#    if gshl1 > thre:
#        print(gshl1)
#    if gshk1 > thre:
#        print(gshk1)
#    if gkvs1 > thre:
#        print(gkvs1)
#    if gkqt3 > thre:
#        print(gkqt3)
#    if gegl2 > thre:
#        print(gegl2)
#    if gegl36 > thre:
#        print(gegl36)
#    if girk > thre:
#        print(girk)
#    if gegl19 > thre:
#        print(gegl19)
#    if gunc2 > thre:
#        print(gunc2)
#    if gcca1 > thre:
#        print(gcca1)
#    if gslo1 > thre:
#        print(gslo1)
#    if gslo2 > thre:
#        print(gslo2)
#    if gkcnl > thre:
#        print(gkcnl)
#    if gnca > thre:
#        print(gnca)

    return gshl1, gshk1, gkvs1, gkqt3, gegl2, gegl36, girk, gegl19, gunc2, gcca1, gslo1, gslo2, gkcnl, gnca, gitr1, gmec10, gosm9, gtrp4
