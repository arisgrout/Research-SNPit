#!/usr/bin/env python3
# %%

#  -----------------------  mvp-origins  -----------------------  Research-SNPit  -----------------------
#  TODO: 'origins' folder v 'operations' folder (with DB setup & checks).


# %%
# FOR OPTIONAL MODULE RELOAD
import sys, importlib

# %%

#  ----------------------  packages

# RELOAD MODULE OPTION IF NEEDED
importlib.reload(sys.modules["helper_functions"])

import requests as re
import pickle  # save API responses to variables
import numpy as np
import pandas as pd
import json
from api_functions import all_snpedia_rsids, all_chromosome_accessions_23me, all_23me_rsids, design_queries
from helper_functions import create_connection, execute_query, create_db

# %%

# // DATABASE SETUP
create_db()

# %%
#  ---  build list of all SNP rsids available on ** SNPedia **
wiki_rsids = all_snpedia_rsids()  # function saves list of rsids to file: "./data/SNPedia_ids.txt"
#  TODO: rename functions (to exclude "SNPedia" & otherwise (to avoid LICENSING dispute) - do so "project-wide" at END)

# %%
#  //--  create prioritized list of rsids & genotypes (for collecting metadata)


#  ---  collect list of SNP rsids ** 23andMe ** (from a 'test' 23andMe raw file)
t3m_df = pd.read_csv(
    "./data/start_23andMeTestData.txt", sep="\t", skiprows=20, names=["rsid", "chromosome", "position", "genotype"]
)


#  ---  Keep only rsids - if they exist in both SNPedia & the 23andMe raw test file
df_select = t3m_df[t3m_df.rsid.isin(wiki_rsids)]
print(f"RSIDs --> 23andMe Raw: {len(t3m_df)} | SNPedia Total: {len(wiki_rsids)} | Both: {len(df_select)}\n\n")

#  --- rm INDEL MARKERS proprietary to 23andMe.
# df_no_indels = exclude_indels(df_select)
#  TODO: design method to handle 23andMe indels and include in query (see api_functions.py for details)

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
