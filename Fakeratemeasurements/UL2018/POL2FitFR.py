import sys, os,array
import math
lepFlav = ['Ele']#['Tau','Ele','Mu']
etas = ['barrel','endcap']
DMs = ['DM0','DM1','DM10','DM11']
era = 'UL2018'
for ilep in lepFlav :
    print '<<<<<<<<<<<<<< lepton : ', ilep
    if ilep is 'Ele' :
        wps = ['Fall17MVAv2WP90_noIso_Iso0p15']#['Fall17MVAv2WP80_noIso_Iso0p15','Fall17MVAv2WP90_noIso_Iso0p15']
    elif ilep is 'Mu' :
        wps = ['Loose_Iso0p15','Medium_Iso0p15']
    elif ilep is 'Tau' :
        wps = ['Medium_VLoose_VLoose','Medium_Tight_VLoose','Medium_VLoose_Tight']
    for iwps in wps :
        print '<<<<<<<<<<<<<< wps : ', iwps
        if ilep is 'Tau' :
            for iDM in DMs :
                print '<<<<<<<<<<<<<< DM : ', iDM
                FitFR_cmd = 'root -l -q \'FitFR.C("%s","%s", "%s","%s")\'' %(ilep,iwps,iDM,era)
                os.system(FitFR_cmd)
        else :
            for ieta in etas :
                print '<<<<<<<<<<<<<< eta : ', ieta
                FitFR_cmd = 'root -l -q \'FitFR.C("%s","%s", "%s","%s")\'' %(ilep,iwps,ieta,era)
                os.system(FitFR_cmd)

