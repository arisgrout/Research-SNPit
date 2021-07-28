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

# // Import Variables
with open("./data/temp/missing_genos.v", "rb") as f:
    missing_genos = pickle.load(f)
with open("./data/temp/missing_rsids.v", "rb") as f:
    missing_rsids = pickle.load(f)

#  ---  design SNPedia genotype queries for API calls (reorienting all "minus" genotypes to "positive" equivalents)
df_queries = design_queries(df_no_indels)
#  TODO: split into sub-functions, setup DB & checkpoints prior to collecting, fix somewhat flawed query design (notes in functions.py)

with open("./data/start_df_queries.v", "wb") as f:
    pickle.dump(df_queries, f)

# %%
"""Get SNPedia Data"""  # ---------------------------------------- BREAK --------------------------------------------


# %%
with open("./data/start_df_queries.v", "rb") as f:
    df = pickle.load(f)
# %%
df.head()
# %%
df = df.set_index("rsid")
# %%
df.tail()
# %%
""" JOIN SNPEDIA INFO TO DF """
#%%
def join_snpedia(df):
    """Use genotype queries to add variant: magnitude, repute, summary, SNPedia_url and fulltext data. Data from SNPedia API (webscrape).

    Args:
        df (pd.DataFrame): 23andMe raw data + query & orientation

    Returns:
        pd.DataFrame: ['rsid', 'chromosome', 'position', 'genotype', 'query', 'orientation', 'magnitude', 'repute', 'summary', 'fullurl', 'fulltext']
    """

    #  ---  add cols for appending genotype metadata from SNPedia
    add_cols = ["magnitude", "repute", "summary", "fullurl", "fulltext"]
    for col in add_cols:
        df[col] = ""

    df = df.set_index("rsid")
    query = df["query"].tolist()

    url = "https://bots.snpedia.com/api.php"
    counter = 0

    for rsid in query:

        # GET REQUEST
        # rsid = 'rs53576(G;G)' # for testing
        params = f"?action=ask&query=[[{rsid}]][[Category:Is a genotype]]|?Genotype|?Magnitude|?Repute|?Summary&format=jsonfm&api_version=2"
        response = re.get(url + params)
        print(response.status_code)

        # slice applicable JSON from text response
        start = response.text.find("results") - 1
        stop = response.text.find("serializer") - 11
        result = "{" + response.text[start:stop].lower() + "}"  # dict

        # convert: str -> dict(JSON) -> dataframe
        result = json.loads(result)
        result = pd.json_normalize(result).iloc[:, 1:-3]

        # insert each col / val pair
        rsid = rsid.partition("(")[0]  # result.columns[-1].split('/')[-1].partition('(')[0].partition('.')[-1]
        for col in result:
            k = col.split(".")[-1]
            v = "".join(str(i) for i in result[col][0])
            df.at[rsid, k] = v

        counter += 1
        print(df[df.index == rsid].iloc[:, 5:-1])
        print(counter / len(query))

    print("complete!")
    return df


# %%

# %%
df_done = join_snpedia(df)
# %%
with open("./data/df_returned.v", "wb") as f:
    pickle.dump((df_done), f)
with open("./data/df_original.v", "wb") as f:
    pickle.dump((df), f)
# %%
""" SORT GENOTYPES: MAGNITUDE & REPUTE ############################################################### """
# Top 10 "Good" based on Mag
# Top 10 "Bad" based on Mag
with open("./data/df_returned.v", "rb") as f:
    df = pickle.load(f)
# %%
df.magnitude.unique()
# %%
rsid_nonZero = df[(df.repute != "") & (~df.magnitude.isin(["", "0"]))].index.tolist()
# %%
rsid_nonZero
# %%
with open("data/rsid_nonZero.v", "wb") as f:
    pickle.dump(rsid_nonZero, f)