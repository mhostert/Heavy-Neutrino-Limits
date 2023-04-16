import os
import numpy as np
import pandas as pd
from . import decay_functions
from . import plot_tools

from hepunits import units as u
from hepunits.units import prefixes as _pre

from math import sqrt

operator_units_dict = {'U2'} #Transforms the bounds on the mixing U2 to bounds on the coefficient of the operator [GeV^-2]

dirac_to_majorana_dic = {'PS': 1, 'BD': 1/sqrt(2), 'C': 1/2, 'UPMNS': 1}
majorana_to_dirac_dic = {'PS': 1, 'BD': sqrt(2), 'C': 2, 'UPMNS': 1}

#function_dic = {'A', 'ReT2KCC', 'ReLHCCC', 'ReCHARMCC', 1}
function_dic = {10001,10002,10003,10004,10005,10006,10007}

#cl90_dict = {'1sigma': sqrt(4.61/2.30), '68.27': sqrt(4.61/2.30), 90: 1, 95: sqrt(4.61/5.99), '2sigma': sqrt(4.61/6.18), 95.45: sqrt(4.61/6.18), 99: sqrt(4.61/9.21)} # 2dof chi^2 relations 
cl90_dict = {'1sigma': sqrt(2.71/1), 68.27: sqrt(2.71/1), 90: 1, 95: sqrt(2.71/3.84), '2sigma': sqrt(2.71/4), 95.45: sqrt(2.71/4), 99: sqrt(2.71/6.63)} # 1dof chi^2 relations 

unit_dict = {'microeV': _pre.micro*u.eV, 'millieV': _pre.milli*u.eV, 'eV': u.eV, 'keV': u.keV, 'MeV': u.MeV, 'GeV': u.GeV, 'TeV': u.TeV, 'PeV': u.PeV}

#https://docs.google.com/spreadsheets/d/1p_fslIlThKMOThGl4leporUsogq9TmgXwILntUZOscg/edit?usp=sharing original version
#https://docs.google.com/spreadsheets/d/1TIpmkgOa63-8Sy75qh0YutI5XdRtiClU3aquUdmjqpc/edit?usp=sharing Josu final version

def load_google_sheet(sheet_id="1TIpmkgOa63-8Sy75qh0YutI5XdRtiClU3aquUdmjqpc", sheet_name = "Umu4", drop_empty_cols=True):
    url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"
    if drop_empty_cols: 
        cols = pd.read_csv(url, header=0, nrows=1, low_memory=False).columns
        cols = cols.str.strip()
        nonempty_cols = [i for i in range(len(cols)) if 'Unnamed' not in cols[i] ]
        df = pd.read_csv(url, header=0, usecols=nonempty_cols, dtype={'year': object})
        df = df[df['id'].isna() == False]
        return df.where(df.notnull(), None)
    else:
        return pd.read_csv(url, header=0)

class Limits:
    def __init__(self, flavor='e', invisible=False, nature='dirac'):
        """ Class that contains all limits on a HNL given a certain flavor structure

        Args:
            flavor (str, optional): The single flavor dominance to be used. Can be 'e', 'mu', or 'tau'. Defaults to 'e'.
            invisible (bool, optional): If True, include only limits that apply to an invisible HNL. Defaults to False.
            nature (str, optional): define the nature of the HNL and applies the corresponding re-scalings. Can be either 'dirac' or 'majorana'. Default to 'dirac'.        
     
        """
        self.flavor = flavor
        subscript = 'e' if self.flavor == 'e' else f'\{self.flavor}'
        self.latexflavor = fr'$C_{{\rm{{HN}}\ell}}^{{{subscript}}}/\Lambda^2 ~(\rm{{GeV}}^{{-2}})$'
        self.invisible = invisible
        self.nature = nature 
        if self.nature != 'dirac':
            if self.nature != 'majorana':
                raise ValueError(f"Define the right nature of the HNL.")
        _df = load_google_sheet(sheet_name=f'Bosonic_charged_{flavor}')
        self.limits = _df.set_index('id')

        self.num_of_limits = self.limits.index.size

        self.limits = self.limits.apply(self.insert_limit, axis = 1)
        self.interp_func_all = self.get_combined_limit_func()

        self.path_to_plot = os.path.join(os.getcwd(), f'plots/bosonic_currents/charged/CHN{self.flavor}_{self.nature}.pdf')

    def insert_limit(self, df):
        """ After the data in the Google Spreadsheet, this function can be used to load the limit,
        its description, raw data, and interpolating function.

        Args:
            df (pandas DataFrame): the dataframe to which the limit is being inserted

        Returns:
            pandas DataFrame: after data is inserted.
        """
        self.get_data(df)
        self.get_data(df, top = True)
        
        return df

    def get_data(self, df, top = False):

        global_path='src/HNLimits/include/data/'
        suffix = "_top" if top else ""
        limit_path = df.file_top if top else df.file_bottom
        
        if (limit_path is None):
            m4, ualpha4 = None, None 
            interp_func = lambda x: np.ones(np.size(x))
        else:
            if (self.invisible and not df.is_invisible):
                m4, ualpha4, interp_func = None, None, None
            else:
                m4, ualpha4 = np.genfromtxt(f"{global_path}{limit_path}", unpack=True)
                # m4, ualpha4 = m4[~np.isnan(m4)], ualpha4[~np.isnan(m4)]
                # m4, ualpha4 = m4[~np.isnan(ualpha4)], ualpha4[~np.isnan(ualpha4)]
                
                if df['CL'] in cl90_dict:
                    # fix the CL to 90%CL
                    ualpha4 = ualpha4 * cl90_dict[df['CL']]
#                else:
#                    raise ValueError(f"CL of experimental data {df['CL']} not defined.")

                if self.nature == 'dirac':
                    if df['hnl_type'] == 'Majorana':
                        # fix the Dirac bounds with the corresponding Majorana factor
                        ualpha4 = ualpha4 * majorana_to_dirac_dic[df['type']]
                if self.nature == 'majorana':
                    if df['hnl_type'] == 'Dirac':
                        # fix the Dirac bounds with the corresponding Majorana factor
                        ualpha4 = ualpha4 * dirac_to_majorana_dic[df['type']]

                if df['units'] in unit_dict:
                    # fix units to HEP units (MeV)
                    m4 = m4 * unit_dict[df['units']]
                    m4 = m4/1000 # set units to GeV

                #if top:
                    coeficiente = 2 if top else 1
                #else:
                #    coeficiente = 1    

                    if df['r_function'] in function_dic:
                    # Re-scale the bound to CC only if needed (if NC was included in the experimental result)
                        if self.flavor == 'e':
                            if df['r_function'] == 10001:
                                for r in range(0,len(ualpha4)):
                                    ualpha4[r] = ualpha4[r]*(decay_functions.ReCHARMCC(m4[r]))**coeficiente
                            if df['r_function'] == 10002:
                                for r in range(0,len(ualpha4)):
                                    ualpha4[r] = ualpha4[r]*(decay_functions.ReT2KCC(m4[r]))**coeficiente
                            if df['r_function'] == 10003:
                                for r in range(0,len(ualpha4)):
                                    ualpha4[r] = ualpha4[r]*(decay_functions.ReLHCCC(m4[r]))**coeficiente
                            if df['r_function'] == 10004:
                                for r in range(0,len(ualpha4)):
                                    ualpha4[r] = ualpha4[r]*(decay_functions.ReBorexinoCC(m4[r]))**coeficiente
                            if df['r_function'] == 10005:
                                for r in range(0,len(ualpha4)):
                                    ualpha4[r] = ualpha4[r]*(decay_functions.ReBEBCCC(m4[r]))**coeficiente        
                        if self.flavor == 'mu':
                            if df['r_function'] == 10005:
                                for r in range(0,len(ualpha4)):
                                    ualpha4[r] = ualpha4[r]*(decay_functions.RmuNUTEVCC(m4[r]))**coeficiente                      
                            if df['r_function'] == 10006:
                                for r in range(0,len(ualpha4)):
                                    ualpha4[r] = ualpha4[r]*(decay_functions.RmuT2KCC(m4[r]))**coeficiente
                            if df['r_function'] == 10007:
                                for r in range(0,len(ualpha4)):
                                    ualpha4[r] = ualpha4[r]*(decay_functions.RmuLHCCC(m4[r]))**coeficiente

                    if df['operator_units'] in operator_units_dict:
                    # we have to re-scale the mixing^2 ualpha4 to the coefficient of the Boson Charged Current operator C_{HN\ell}^\alpha/Lambda^2
                        ualpha4 = (sqrt(2)/(246.22)**2)*np.sqrt(ualpha4) #[GeV^-2] with m4 in GeV. See Eq. 4.2 of the draft
            
                    # order data points
                    order = np.argsort(m4)
                    _m4 = m4[order]
                    _ualpha4 = ualpha4[order]
                    # interpolation 
                    interp_func = plot_tools.log_interp1d(_m4, _ualpha4, kind='linear', bounds_error=False, fill_value=None, assume_sorted=False)    

                else:
                    raise ValueError(f"HNL mass units of {df['units']} not defined.")
        
        df[f'm4{suffix}'] = m4
        df[f'ualpha4{suffix}'] = ualpha4
        df[f'interp_func{suffix}'] = interp_func
        
        return df

    def get_combined_limit_func(self, xrange=[1, 1e5], npoints = 1000, logx=True, include_cosmo=False):
        """ Returns a function that takes m4 [GeV] and gives bound on Ualpha4^2

        Args:
            xrange (list, optional): range of HNL mass in GeV. Defaults to [1, 1e5].
            npoints (int, optional): number of mass points to return the limit for. Defaults to 1000.
            logx (bool, optional): If True, return logarithmically spaced points. Defaults to True.
            include_cosmo (bool, optional): If True, include cosmological bounds in the combination. Defaults to False.

        Returns:
            scipy.interpolate.interp1d(logx, logy, kind=kind, **kwargs): interpolating function.
        """

        y=[]
        if logx:
            x=np.geomspace(*xrange,npoints)
        else:
            x=np.linspace(*xrange,npoints)
        
        for _, limit in self.limits.iterrows():
            
            ## FIX ME -- no functionality for gaps between top of constraint and other bounds
            # check if it is a closed contour with top and bottom files
            if (not (limit['file_top'] is None) or (limit['interp_func'] is None)) or (include_cosmo and 'bbn' in limit.id):
                continue
            else:
                y.append(limit.interp_func(x))
        
        if len(y)==0:
            print("No limits were found when constructing the combination")
            return None
        else:
            y = np.array(y)
            z = np.ones(len(y[0]))
            for i in range(0,np.shape(y)[0]):
                for j in range(0, np.size(y[i])):
                    if y[i,j] < z[j]:
                        z[j] = y[i,j]

            return plot_tools.log_interp1d(x, z, kind='linear', bounds_error=False, fill_value=None, assume_sorted=False)

