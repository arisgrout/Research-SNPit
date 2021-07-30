"""
In:
* missing_rsids.v   12345
* missing_genos.v   [C:T]

Out:
*

Trigger:
* display available research report in web-app
* collect research pubs for missing rsids
* collect metadata on rsid genotypes from SNPedia
"""

# TODO Find alternative to SNPedia for gathering metadata.

# Adds this file to top-level interpretter path, so runs "as if" from top.
# import sys
# from pathlib import Path
# sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# %%
# // LOAD MODULES / PACKAGES
import sys
#import os.path
from os import path
from importlib import reload
from collections import OrderedDict
import numpy as np
import pandas as pd
import pickle
import copy
import modules.functions as fn

# from modules.functions import create_connection, execute_query, check_db

# // RELOAD ALL CHANGED MODULES (JUPYTER SPECIAL FUNCTION) 
# // autoreloads on each cell execution.
# https://stackoverflow.com/questions/34377270/reload-module-from-a-folder
%load_ext autoreload
%autoreload 2

# // RELOAD ALL MODULES
# for module in sys.modules.values():
#     reload(module)
# https://stackoverflow.com/questions/45405600/how-to-reload-all-imported-modules

# // RELOAD ONE MODULE
# reload(functions)

# %%

# // Import Variables
with open("./data/temp/missing_genos.v", "rb") as f:
    df = pickle.load(f)

# add columns for genotypes table
add_cols = ["query","orient37","orient38","magnitude", "repute", "summary", "fullurl", "fulltext"]
for col in add_cols:
    df[col] = ""

# API request SNP orientations
df_orient = fn.get_orientation(df)

# save df_orient for backup
with open("./data/temp/df_orient.v", "wb") as f:
    pickle.dump(df_orient, f)

# TODO some values casted as 'missing' when not actually missing from HTML. API query just failing to find for unknown reason. E.g. uncomment print(response.text) in function, then run line below. Compare to HTML at addresses:
# fn.get_orientation(df_orient[df_orient.orient37 == 'missing'])
# 'plus' stabilizedorientation: https://snpedia.com/index.php/rs28940586
# FAILED 'minus' stabilizedorientation: https://www.snpedia.com/index.php/rs3094315
# see github issue for details.
# save missing to CSV for later eval
#df_orient[df_orient.orient37 == 'missing'].to_csv('./data/#sOrient37_missing.csv', index=False)

# flip SNP genotypes if orient37 is 'minus'
df_flip = copy.deepcopy(df_orient)
df_flip['genotype'] = df_orient.apply(lambda x: fn.snp_flipper(x.genotype) if x.orient37 == 'minus' else x.genotype, axis=1)

# construct queries for requesting genotype metadata
df_q = fn.get_queries(df_flip)

# Queries missing orient37 can't be used so they are deleted here.
df_q.loc[(df_q.orient37 == 'missing'), 'query'] = ''

#%%
# Query remaining info from SNPedia & save to file: df_snpedia.
df_g = fn.join_snpedia(df_q)

# %%
# JUST THIS TIME THEN DELETE
df_g = df_snpedia
fn.create_db()
# %%
df_g["date"] = pd.datetime.now()
df_g = df_g.replace(r"^\s*$", np.nan, regex=True) # replace empty strings with np.nan
cnx = fn.create_connection("data/SNP_db.sqlite")
df_g.to_sql('genotypes', cnx, if_exists='append', index=False)