import os

# Bounds on |Umu4|^2 vs m4
muon_bounds = {} 
PATH = os.path.dirname(os.path.abspath(__file__))+'/digitized'

muon_bounds['BEBC'] = {
         'file' : [PATH+"/umu4/BEBC.dat", PATH+"/Umu4/BEBC_top.dat"], ##file(s), if >1 then it is a closed regin
         'units' : 'GeV', # units of m 
         'year' : '', # year of th constrait
         'ref' : '', # reference to the paper used indigitizig
         'invisible' : False, # does the bound apply to fully invisibleHNL?
         'CL' : 90, # confidence level 
        }

muon_bounds['CHARMII'] = {
         'file' : [PATH+"/umu4/CHARMII.dat", PATH+"/Umu4/CHARMII_top.dat"],
         'units' : 'GeV',
         'year' : '',
         'ref' : '',
         'invisible' : False,
         'CL' : 90,
        }

muon_bounds['DELPHI_97'] = {
         'file' : PATH+"/umu4/DELPHI_1997.dat",
         'units' : 'GeV',
         'year' : '1997',
         'ref' : '',
         'invisible' : False,
         'CL' : 90,
        }

muon_bounds['FMMF'] = {
         'file' : [PATH+"/umu4/FMMF.dat", PATH+"/Umu4/FMMF_top.dat"],
         'units' : 'GeV',
         'year' : '',
         'ref' : '',
         'invisible' : False,
         'CL' : 90,
        }

muon_bounds['atre_knumu'] = {
         'file' : PATH+"/umu4/K_to_numu.dat",
         'units' : 'GeV',
         'year' : '',
         'ref' : '',
         'invisible' : True,
         'CL' : 90,
        }

muon_bounds['L3'] = {
         'file' : PATH+"/umu4/L3.dat",
         'units' : 'GeV',
         'year' : '',
         'ref' : '',
         'invisible' : False,
         'CL' : 90,
        }

muon_bounds['NA3'] = {
         'file' : [PATH+"/umu4/NA3.dat",PATH+"/Umu4/NA3_top.dat"],
         'units' : 'GeV',
         'year' : '',
         'ref' : '',
         'invisible' : False, ## doublecheck #
         'CL' : 90,
        }

muon_bounds['NuTeV'] = {
         'file' : [PATH+"/umu4/NuTeV.dat", PATH+"/Umu4/NuTeV_top.dat"],
         'units' : 'GeV',
         'year' : '',
         'ref' : '',
         'invisible' : False,
         'CL' : 90,
        }

muon_bounds['PS191'] = {
         'file' : [PATH+"/umu4/PS191_bottom.dat", PATH+"/Umu4/PS191_top.dat"],
         'units' : 'GeV',
         'year' : '',
         'ref' : '',
         'invisible' : False,
         'CL' : 90,
        }

muon_bounds['E949'] = {   ## WHERE FROM? ##
         'file' : PATH+"/peak_searches/E949.dat",
         'units' : 'MeV',
         'year' : '',
         'ref' : '',
         'invisible' : True,
         'CL' : 90,
        }

muon_bounds['KEK'] = {    ## WHERE FROM? ##
         'file' : PATH+"/peak_searches/KEK_K_nu_mu.dat",
         'units' : 'MeV',
         'year' : '',
         'ref' : '',
         'invisible' : True,
         'CL' : 90,
        }

muon_bounds['NA62_17'] = {
         'file' : PATH+"/peak_searches/NA62_(2017)_K_nu_mu.dat",
         'units' : 'MeV',
         'year' : '',
         'ref' : '',
         'invisible' : True,
         'CL' : 90,
        }

muon_bounds['NA62_18'] = {
         'file' : PATH+"/peak_searches/NA62_2018_K_nu_mu.dat",
         'units' : 'MeV',
         'year' : '',
         'ref' : '',
         'invisible' : True,
         'CL' : 90,
        }

muon_bounds['deGouvea_lowmass'] = {
         'file' : PATH+"/deGouvea/lowmassUmu4.dat",
         'units' : 'GeV',
         'year' : '',
         'ref' : '',
         'invisible' : True,
         'CL' : 99,
        }

muon_bounds['deGouvea_LNU'] = {
         'file' : PATH+"/deGouvea/lepton_non_universality.dat",
         'units' : 'GeV',
         'year' : '',
         'ref' : '',
         'invisible' : True,
         'CL' : 99,
        }

# muon_bounds['kusenko_all'] = {  # WHAT FIG???
#          'file' : PATH+"/Kusenko/Umu4.dat",
#          'units' : 'MeV',
#          'year' : '',
#          'ref' : '',
#          'invisible' : False,
#          'CL' : 90,
#         }

# muon_bounds['kusenko_all_v2'] = {  # WHAT FIG???
#          'file' : PATH+"/Kusenko/Umu4_v2.dat",
#          'units' : 'MeV',
#          'year' : '',
#          'ref' : '',
#          'invisible' : False,
#          'CL' : 90,
#         }

# muon_bounds['kusenko_all2'] = {  # WHAT FIG???
#          'file' : PATH+"/Kusenko/Umu4_2.dat",
#          'units' : 'MeV',
#          'year' : '',
#          'ref' : '',
#          'invisible' : False,
#          'CL' : 90,
#         }

# muon_bounds['kusenko_all3'] = {  # WHAT FIG???
#          'file' : PATH+"/Kusenko/Umu4_3.dat",
#          'units' : 'MeV',
#          'year' : '',
#          'ref' : '',
#          'invisible' : False,
#          'CL' : 90,
#         }

# muon_bounds['kusenko_all4'] = {  # WHAT FIG???
#          'file' : PATH+"/Kusenko/Umu4_4.dat",
#          'units' : 'MeV',
#          'year' : '',
#          'ref' : '',
#          'invisible' : False,
#          'CL' : 90,
#         }

# muon_bounds['kusenko_all5'] = {  # WHAT FIG???
#          'file' : PATH+"/Kusenko/Umu4_5.dat",
#          'units' : 'MeV',
#          'year' : '',
#          'ref' : '',
#          'invisible' : False,
#          'CL' : 90,
#         }

muon_bounds['SuperK_19'] = {
         'file' : PATH+"/umu4/SuperK_2019.dat",
         'units' : 'GeV',
         'year' : '2019',
         'ref' : '',
         'invisible' : False,
         'CL' : 90,
        }

muon_bounds['T2K_19'] = {
         'file' : PATH+"/umu4/T2K_2019.dat",
         'units' : 'MeV',
         'year' : '2019',
         'ref' : '',
         'invisible' : False,
         'CL' : 90,
        }

muon_bounds['NA62_21'] = {
         'file' : PATH+"/umu4/NA62_2021.dat",
         'units' : 'GeV',
         'year' : '2021',
         'ref' : '',
         'invisible' : True,
         'CL' : 90,
        }

muon_bounds['PIENU_19_lowT'] = {
         'file' : PATH+"/umu4/PIENU_2019_lowT.dat",
         'units' : 'MeV',
         'year' : '2019',
         'ref' : '',
         'invisible' : True,
         'CL' : 90,
        }

muon_bounds['PIENU_19_highT'] = {
         'file' : PATH+"/umu4/PIENU_2019_highT.dat",
         'units' : 'MeV',
         'year' : '2019',
         'ref' : '',
         'invisible' : True,
         'CL' : 90,
        }

muon_bounds['PSI_87'] = {
         'file' : PATH+"/umu4/PSI_1987.dat",
         'units' : 'MeV',
         'year' : '1987',
         'ref' : '',
         'invisible' : True,
         'CL' : 90,
        }

muon_bounds['DELPHI_shortlived'] = {
         'file' : PATH+"/umu4/DELPHI_short_95CL.dat",
         'units' : 'GeV',
         'year' : '',
         'ref' : '',
         'invisible' : False,
         'CL' : 95,
        }

muon_bounds['DELPHI_longlived'] = {
         'file' : PATH+"/umu4/DELPHI_long_95CL.dat",
         'units' : 'GeV',
         'year' : '',
         'ref' : '',
         'invisible' : False,
         'CL' : 95,
        }

muon_bounds['SIN_87'] = {
         'file' : PATH+"/umu4/SIN_1987.dat",
         'units' : 'MeV',
         'year' : '1987',
         'ref' : '',
         'invisible' : False,
         'CL' : 90,
        }
