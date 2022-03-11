import os

# Bounds on |Umu4|^2 vs m4
bounds = {} 
PATH = os.path.dirname(os.path.abspath(__file__))+'/../digitized'

# files_utau4 = [ 
#         PATH+"/utau4/EWPD_Bolton.dat",
#         PATH+"/utau4/Bdecays_Bolton.dat",
#         PATH+"/utau4/l_universality_Bolton.dat",

#         ]

bounds['Borexino'] = {
         'file' : PATH+"/utau4/Borexino_RPlestid.dat", ##file(s), if >1 then it is a closed regin
         'units' : 'MeV', # units of m 
         'year' : '2020', # year of th constrait
         'ref' : 'arxiv.org/abs/2010.09523', # reference to the paper used indigitizig
         'invisible' : False, # does the bound apply to fully invisibleHNL?
         'CL' : None, # confidence level 
        }

bounds['NOMAD'] = {
        
         'file' : PATH+"/utau4/NOMAD_UL_90CL_2002.dat",
         'units' : 'MeV',
         'year' : '2002',
         'ref' : 'https://arxiv.org/abs/hep-ex/0101041',
         'invisible' : False, # does the bound apply to fully invisibleHNL?
         'CL' : 90, # confidence level 
        }

bounds['CHARM'] = {
        
         'file' : PATH+"/utau4/CHARM_UL_90CL_2002.dat",
         'units' : 'MeV',
         'year' : '2002',
         'ref' : 'https://arxiv.org/abs/hep-ph/0208075',
         'invisible' : False, # does the bound apply to fully invisibleHNL?
         'CL' : 90, # confidence level 
        }

bounds['T2K_marginalized'] = {
         'file' : PATH+"/utau4/T2K_UL_90CL_2019.dat",
         'units' : 'MeV',
         'year' : '2019',
         'ref' : 'https://arxiv.org/abs/1902.07598',
         'invisible' : False,
         'CL' : 90,
        }

bounds['DELPHI_short'] = {
         'file' : PATH+"/utau4/DELPHI_short_95CL.dat",
         'units' : 'GeV',
         'year' : '1997',
         'ref' : '10.1007/s002880050370',
         'invisible' : False,
         'CL' : 95,
        }

bounds['DELPHI_long'] = {
         'file' : PATH+"/utau4/DELPHI_long_95CL.dat",
         'units' : 'GeV',
         'year' : '1997',
         'ref' : '10.1007/s002880050370',
         'invisible' : True,
         'CL' : 95,
        }

bounds['argoneut'] = {
         'file' : [PATH+"/2106.13684/argoneut_bottom.dat",PATH+"/2106.13684/argoneut_top.dat"],
         'units' : 'GeV',
         'year' : '2021',
         'ref' : 'https://arxiv.org/abs/2106.13684',
         'invisible' : False,
         'CL' : 90,
        }

bounds['charm_boiarska'] = {
         'file' : [PATH+"/2107.14685/CHARM_bottom.dat",PATH+"/2107.14685/CHARM_top.dat"],
         'units' : 'GeV',
         'year' : '2021',
         'ref' : 'arxiv.org/abs/2107.14685',
         'invisible' : False,
         'CL' : 90,
        }







