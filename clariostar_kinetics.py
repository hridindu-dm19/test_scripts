# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 14:21:46 2022

@author: hroychow
"""

from pathlib import Path
import csv
# import openpyxl
import re

import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt
import seaborn as sns
#%%


wdir = Path(r'Y:\Research\Hridindu\Clariostar')
time_course = wdir / '220422_0916_20220421 - LgHT [3-19] fluorescence timecourse_.csv'

ndir = Path(r'H:\Notebook')
notebook = ndir / '0003_20220417_HT LgHT complementation pilot assay.xlsx'
# %%
def load_clariostar_timecourse(time_course, sheet_name = None):
    if str(time_course).endswith('.csv'):
        
        reader = csv.reader(open(time_course,'r', encoding = "ISO-8859-1"))
        
        for line in reader:
            if len(line)>0:
                
                if line[0].strip().startswith('Number of cycles'):
                    n_cycles = int(' '.join(line[0].split()).split(' ')[3])
                if line[0].strip().startswith('Cycle time'):
                    cycle_time = int(' '.join(line[0].split()).split(' ')[3])
            
            else:
                pass
                
        time=list(range(0,n_cycles*cycle_time,cycle_time))    
        
        data = pd.read_csv(time_course,
                           encoding = "ISO-8859-1",
                           engine='python',
                           skiprows=(59)
                           )
        
        data = data.T.reset_index(drop = False).T.reset_index(drop = True)
        data[['Well', 0]] = data[0].str.split(':', expand = True)
        data = data.set_index('Well', drop = True)
        data.columns = time
        
    elif str(time_course).endswith('.xlsx'):
        meta_data = pd.read_excel(time_course, 
                              sheet_name,
                              nrows= 58
                              )
        
        mask = meta_data['Testname: Janelia_646_kinetic'].str.contains(pat = 'Number of cycles')
        n_cycles = meta_data[mask == True].dropna(axis = 1).iloc[0][0]
        n_cycles =int(' '.join(n_cycles.split()).split(' ')[3])
        
        mask = meta_data['Testname: Janelia_646_kinetic'].str.contains(pat = 'Cycle time')
        cycle_time = meta_data[mask == True].dropna(axis = 1).iloc[0][0]
        cycle_time =int(' '.join(cycle_time.split()).split(' ')[3])
        
        time=list(range(0,n_cycles*cycle_time,cycle_time))
        
        data = pd.read_excel(time_course, 
                             sheet_name,
                             header= 59
                             )
        data = data.T.reset_index(drop = False).T.reset_index(drop = True)
        data = data.iloc[1:,:]
        data[['Well', 0]] = data[0].str.split(':', expand = True)
        data = data.set_index('Well', drop = True)
        data.columns = time
    
    for col in data.columns:
        data[col] = pd.to_numeric(data[col], errors='coerce')
   
    return data

def load_platemap(platemap_file, platemap_sheet):
    
    """Read platemap file """
    platemap = pd.read_excel(platemap_file, 
                             sheet_name = platemap_sheet)
    platemap.columns = platemap.columns.map(str)
    
    """Read conditions into their own series"""
    row_conditions = platemap[[col for col in platemap.columns if bool(re.match('0[1-9]|1[0-2]', col)) == False]]
    row_conditions = row_conditions.iloc[:8]  
    col_conditions = platemap[~platemap['Row'].str.contains(pat = r'\b[A-H]\b')].T
    col_con_header = col_conditions.iloc[0]
    col_conditions = col_conditions.iloc[1:13]
    col_conditions.columns = col_con_header
    col_conditions['Column'] = col_conditions.index 
    
    '''Load only what's in the 96-well map'''
    platemap = platemap[platemap['Row'].str.contains(pat = r'\b[A-H]\b')]
    cols = [col for col in platemap.columns if bool(re.match('Row|0[1-9]|1[0-2]', col)) == True]
    platemap = platemap[cols]
    platemap.columns = platemap.columns[platemap.columns.str.contains(pat = 'Row|0[1-9]|1[0-2]')]
    
    
    '''Melt platemap and merge conditions'''
    platemap = pd.melt(platemap, 
                    id_vars = 'Row',
                    # value_vars= platemap.columns,
                    value_name= 'Condition',
                    var_name = 'Column'
                   ).dropna(axis = 0, how = 'any')
    platemap = platemap.merge(row_conditions, on = 'Row')    
    platemap = platemap.merge(col_conditions, on = 'Column')  
    platemap['Well'] = platemap['Row'].astype(str) + platemap['Column'].astype(str)
    platemap = platemap.drop(columns = ['Row','Column'])
    '''Rename conditions'''
    new_cols = ['Condition %s' % i for i in list(range(1,len(platemap.columns)))]
    new_cols = {i : j for i,j in zip(platemap.columns, new_cols)}
    platemap.rename(columns = new_cols, inplace=True)
    platemap = platemap.fillna('')

    
    return platemap

def label_plate(platemap, data, long = True):
    labeled_platemap_data = platemap.merge(data, left_on = 'Well', right_on = 'Well')
    if long == True:
        pattern = r'Condition|Well'
        id_vars = [col for col in labeled_platemap_data.columns if re.search(pattern, str(col))]
        labeled_platemap_data = labeled_platemap_data.melt(id_vars = id_vars,
                                                           var_name = 'Time (s)',
                                                           value_name = 'signal')
    else:
        pass
    return labeled_platemap_data

# %%
class Plate96:
        
    def combine_conditions(self): 
        """Merge the platemap with the kinetics data"""
        if self.data is not None:
            pattern = r'Condition'
            cols = [col for col in self.data.columns if re.search(pattern, str(col))]
            if len(cols)>=2:
                self.data['Condition'] = self.data.filter(regex = pattern).apply(' '.join, 1)
                self.data = self.data.drop(columns = cols)
            else:
                print("""There's only one condition column, so I refuse to delete it into oblivion.""")
                pass
        else:
            print("""self.data is None, check if u fucked up""")
            pass
        return self
    
    def format_pzfx(self): 
        """Creates a dataframe that is easily copied/pasted into graphpad prism9 x vs ny plot"""
        data = self.combine_conditions().data
        pivot = pd.pivot_table(data=data,
                                index = 'Time (s)',
                                columns = ['Condition', 'Well'],
                                values = 'signal'
                                )
        pivot = pivot.reset_index()
        # pivot.rename(columns = {'x': 'Time'}, inplace=True)
        pivot = pivot.T 
        pivot.index = [i for (i,j) in pivot.index]
        pivot = pivot.T.set_index('Time (s)')
        
        return pivot
    
    def __init__(self, platemap_file=None, platemap_sheet=None, data_file=None, data_sheet = None):
        
        self.platemap = load_platemap(platemap_file, platemap_sheet) if (platemap_file and platemap_sheet) is not None else None
        self.data = load_clariostar_timecourse(data_file, data_sheet)            if data_file is not None else None
        self.data = label_plate(self.platemap, self.data)            if (platemap_file and platemap_sheet and data_file) is not None else None
        
        return None
        

   #%%
    
if __name__ == "__main__":    
    print('poop poop')
    
    
    
expt_4_path = ndir / '0004_20220425_HT complementation' / '0004_20220425_HT LgHT complementation pilot assay.xlsx'
data_sheet = '220428_0846_20220427'

expt4 = Plate96(
                platemap_file = expt_4_path,
                platemap_sheet='PlateMap',
                data_file = expt_4_path,
                data_sheet = data_sheet
                )

#%%
pivoted = expt4.format_pzfx()








    