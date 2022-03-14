import numpy as np
import pandas as pd
from . import plot_tools

from hepunits import units as u
from hepunits.units import prefixes as _pre

unit_dict ={'microeV': _pre.micro*u.eV, 'millieV': _pre.milli*u.eV, 'eV': u.eV, 'keV': u.keV, 'MeV': u.MeV, 'GeV': u.GeV, 'TeV': u.TeV, 'PeV': u.PeV}

def load_google_sheet(sheet_id="1p_fslIlThKMOThGl4leporUsogq9TmgXwILntUZOscg", sheet_name = "Umu4", drop_empty_cols=True):
    url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"
    if drop_empty_cols: 
        cols = pd.read_csv(url, header=0, low_memory=False, nrows=1).columns
        nonempty_cols = [i for i in range(len(cols)) if 'Unnamed' not in cols[i] ]
        df = pd.read_csv(url, header=0, usecols=nonempty_cols, dtype={'year': object})
        return df.where(df.notnull(), None)
    else:
        return pd.read_csv(url, header=0)

class limits:
    def __init__(self, flavor='e', invisible=False):
        self.flavor = flavor
        subscript = 'e' if self.flavor == 'e' else f'\{self.flavor}'
        self.latexflavor = fr'$|U_{{{subscript} N}}|^2$'
        self.invisible = invisible
        self.limits = load_google_sheet(sheet_name=f'U{flavor}4')

        self.num_of_limits = self.limits.index.size

        self.limits = self.limits.apply(self.load_limit, axis = 1)
        self.interp_func_all = self.get_combined_limit_func()

    # get the data for the limits and interpolations
    def load_limit(self, df):        
        
        self.get_data(df)
        self.get_data(df, top=True)
        
        return df

    def get_data(self, df, top = False):

        global_path='data/'
        suffix = "_top" if top else ""
        limit_path = df.file_top if top else df.file_bottom
        
        if (limit_path is None):
            m4, ualpha4 = None, None 
            interp_func = lambda x: np.ones(np.size(x))
        else:
            if self.invisible & ~df.is_invisible:
                m4, ualpha4, interp_func = None, None, None
            else:
                m4, ualpha4 = np.genfromtxt(f"{global_path}{limit_path}", unpack=True)
                
                if df['units'] in unit_dict:
                    # fix units to HEP units (MeV)
                    m4 = m4 * unit_dict[df['units']]
            
                    # order data points
                    order = np.argsort(m4)
                    m4 = m4[order]
                    ualpha4 = ualpha4[order]
                    # interpolation 
                    interp_func = plot_tools.log_interp1d(m4, ualpha4, kind='linear', bounds_error=False, fill_value=None, assume_sorted=False)    

                    

                else:
                    raise ValueError(f"HNL mass units of {df['units']} not defined.")
        
        df[f'm4{suffix}'] = m4
        df[f'ualpha4{suffix}'] = ualpha4
        df[f'interp_func{suffix}'] = interp_func
        
        return df

    ## Return a function that takes m4 [GeV] and gives bound on Ualpha4^2
    def get_combined_limit_func(self, xrange=[1, 1e5], npoints = 1000, logx=True, no_cosmo=True):
        y=[]
        if logx:
            x=np.geomspace(*xrange,npoints)
        else:
            x=np.linspace(*xrange,npoints)
        
        for _, limit in self.limits.iterrows():
            ## FIX ME -- no functionality for gaps between top of constraint and other bounds
            # check if it is a closed contour with top and bottom files
            if ((limit['file_top'] is None) | (limit['interp_func'] is None)) | (no_cosmo and 'bbn' in limit.id):
                continue
            else:
                y.append(limit.interp_func(x))
        y = np.array(y)
        z = np.ones(len(y[0]))
        for i in range(0,np.shape(y)[0]):
            for j in range(0, np.size(y[i])):
                if y[i,j] < z[j]:
                    z[j] = y[i,j]

        return plot_tools.log_interp1d(x, z, kind='linear', bounds_error=False, fill_value=None, assume_sorted=False)

