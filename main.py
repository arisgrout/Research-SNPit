"""
In:
* 23andMe_raw_genotype.txt

Out:
* missing_rsids.v   12345
* missing_genos.v   [C:T]
* dataframe.v       immediate report for app dashboard
* dataframe.csv     immediate csv for app dashboard

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

# %%
# // RELOAD ALL CHANGED MODULES (JUPYTER SPECIAL FUNCTION) 
# // autoreloads on each cell execution.
# https://stackoverflow.com/questions/34377270/reload-module-from-a-folder
%load_ext autoreload
%autoreload 2

# %%
# // RELOAD ALL MODULES
# for module in sys.modules.values():
#     reload(module)
# https://stackoverflow.com/questions/45405600/how-to-reload-all-imported-modules

# // RELOAD ONE MODULE
# reload(functions)

# %%
# // FIRST TIME - SETUP FUNCTIONS
if not path.exists("./data/SNP_db.sqlite"):
    # // DATABASE SETUP
    fn.create_db()

if not path.exists("./data/SNPedia_ids.txt"):
    #  ---  build list of all SNP rsids available on ** SNPedia **
    # function saves list of rsids to file: "./data/SNPedia_ids.txt"
    wiki_rsids = fn.all_snpedia_rsids()  
    #  TODO: rename functions (to exclude "SNPedia" & otherwise (to avoid LICENSING dispute) - do so "project-wide" at END)
    #  TODO: run this function separately on sleep timer to update list of SNPedia RSIDS periodically. 

# %%

# %%
# load 23andMe Data (raw test file)
df = pd.read_csv(
    f"./data/temp/{fn.find_files('data/temp', '*.txt')[0]}", sep="\t", skiprows=20, names=["rsid", "chromosome", "position", "genotype"]
    )

print(len(df))

# -------- check if RSIDs are missing from publications table
query = "SELECT DISTINCT rsid FROM rsid_pubs"
rsid_compare = fn.check_db(query=query, compare=[df.rsid], path="data/SNP_db.sqlite")

with open("./data/temp/missing_rsids.v", "wb") as f:
    pickle.dump(rsid_compare["missing"], f)  # save missing_rsids for later collection

print(len(rsid_compare), len(rsid_compare["missing"]), len(rsid_compare["available"]))


# -------- check if SNPedia data is missing from genotype table
query = "SELECT rsid, genotype FROM genotypes"
geno_compare = fn.check_db(query=query, compare=[df.rsid.tolist(), df.genotype.tolist()], path="data/SNP_db.sqlite")

with open("./data/temp/missing_genos.v", "wb") as f:
    pickle.dump(geno_compare["missing"], f)  # save missing_genos for later collection

print(len(geno_compare), len(geno_compare["missing"]), len(geno_compare["available"]))


# ------- select available genotypes from db to match on user
cnx = fn.create_connection("data/SNP_db.sqlite")
genos = f"{geno_compare['available']}"[1:-1]  # convert to str, rm []
query = f"SELECT * FROM genotypes WHERE (rsid, genotype) IN (VALUES {genos});"
geno_df = pd.read_sql_query(query, cnx)


# ------- sort data by magnitude & split by repute [good & bad]
geno_df = geno_df.sort_values("magnitude", ascending=False)
g_df = copy.deepcopy(geno_df[geno_df.repute == "good"][:50])
b_df = copy.deepcopy(geno_df[geno_df.repute == "bad"][:50])
df = g_df.append(b_df, ignore_index=True).iloc[:, 1:]


# ------- import summaries and publication metadata
rsids = f"{df.rsid.tolist()}"[1:-1]  # convert to str and remove [] brackets

# most_cited summaries
query = f"""
SELECT rsid, abs_summary, most_cited, most_recent, pmids
FROM rsid_summaries 
WHERE (rsid IN ({rsids}) AND most_cited = 1);
"""
c_df = pd.read_sql_query(query, cnx)

# most_recent summaries
query = f"""
SELECT rsid, abs_summary, most_cited, most_recent, pmids
FROM rsid_summaries 
WHERE (rsid IN ({rsids}) AND most_recent = 1);
"""
r_df = pd.read_sql_query(query, cnx)

# publications metadata
query = f"""
SELECT rsid, pmid, date, n_citedby, ncbi_url 
FROM rsid_pubs 
WHERE rsid IN ({rsids})
"""
pub_df = pd.read_sql_query(query, cnx)


# ------- construct final dataframe

# %%
c_df = c_df.set_index("rsid")
r_df = r_df.set_index("rsid")
df = df.set_index("rsid")
# %%
df = df.join(c_df[["abs_summary", "pmids"]], how="left", rsuffix="_mC")
df = df.join(c_df[["abs_summary", "pmids"]], how="left", rsuffix="_mR")
df = df.rename(columns={"abs_summary": "abs_summary_mC", "pmids": "pmids_mC"})
df["fullurl"] = df.fullurl.apply(lambda x: x.replace("bots.", ""))
df["fullurl"] = df.fullurl.apply(lambda x: x.rsplit("(")[0] + "(" + x.rsplit("(")[1].upper())
# %%
df.columns
# %%
with open("data/temp/dataframe.v", "wb") as f:
    pickle.dump(df, f)
# %%
df.to_csv("data/temp/dataframe.csv")


# ------- trigger app.py to display dataframe in browser


# ------- trigger scrape for missing information


# %%
# with open("data/temp/dataframe.v", "rb") as f:
#    df = pickle.load(f)
