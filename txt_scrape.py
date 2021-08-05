#!/usr/bin/env python3

"""
In:
(0) missing_rsids.v
    ** List: 
    [   rsids,...   ]

Out:
(1) SNP_db.sqlite/rsid_pubs
    ** SQL Table:
    [   rsid, date, title, journal, authors, abstract, body, pmid,
    pmc, doi, ncbi_url, pdf_url, citation, n_citedby, pmc_citedby     ]

Trigger:
# TODO signal completion?
"""

# TODO Find alternative to SNPedia for gathering metadata.
# TODO Refactor to non-blocking I/O
# TODO asyncio api calls

# Adds this file to top-level interpretter path, so runs "as if" from top.
# import sys
# from pathlib import Path
# sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# %%
# // LOAD MODULES / PACKAGES
import sys
#import os.path
import os
from os import path
from importlib import reload
import numpy as np
import pandas as pd
import requests as re
import pickle
import json
import pubmed_parser as pp  # TODO try using
from metapub import PubMedFetcher, FindIt
import scipdf
import modules.functions as fn

NCBI_API_KEY = os.environ["NCBI_API_KEY"]
# http://biopython.org/DIST/docs/tutorial/Tutorial.html # 9.16.3  Searching for citations
# from Bio import Entrez

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
# // GET TEXT FUNCTIONS / API CALLS


def test_from_rsid(rsids, start_rsid):
    """Continue collecting publications for rsids in list, beginning with start_rsid

    Args:
        rsids (list): list of rsids to collect publications on
        start_rsid (str): rsid identifier to resume collecting publications on

    Returns:
        runtime_rsids (list): [start_rsid, onward...]
        start_rsid (str): starting rsid
        start_idx (str): starting rsid index
        rsids (list): [original list of ALL rsids]
    """
    start_idx = rsids.index(start_rsid)  # start_rsid index
    print(f"STARTING POINT SET TO: | INDEX: {start_idx} / {len(rsids)} | RSID: {rsids[start_idx]}")
    runtime_rsids = rsids[start_idx:]  # runtime rsids
    return runtime_rsids, start_rsid, start_idx, rsids


def get_pmids(rsid):
    """Get all PMID identifiers for papers referencing RSID in question from LitVar archive.

    Args:
        rsid (str): 'rs12345' identifier

    Returns:
        rsid, pmids (tuple): ('rs12345', ['pmids',...])
    """

    url = "https://www.ncbi.nlm.nih.gov/research/bionlp/litvar/api/v1"
    endpoint = "/public/rsids2pmids"
    params = {"rsids": rsid}
    response = re.get(url + endpoint, params=params)
    try:
        pmids = response.json()[0]["pmids"]
    except Exception as e:
        msg = f"{fn.date_time()} - RSID: {rsid}. No articles for RSID. Function: get_pmids. {e}"
        print(msg)
        errors.append(msg)
        pmids = None
        pass

    return (rsid, pmids)


def missing_pmids(rsid, pmids):
    """Check DB, return list of missing pmids (article ids) for given rsid

    Args:
        rsid (str): rs12345
        pmids (list): list of all pmids for given rsid; pulled from LitVar

    Returns:
        list: list of pmids not yet in DB
    """    

    cnx = fn.create_connection("data/SNP_db.sqlite")
    in_db = pd.read_sql_query(f"SELECT pmid FROM rsid_pubs WHERE rsid='{rsid}';", cnx)
    in_db = [] if len(in_db) == 0 else in_db.pmid.apply(int)
    pmids = pd.Series(pmids, dtype=int)
    missing_pmids = list(set(pmids).difference(in_db))
    missing_pmids = None if (len(missing_pmids) == 0) else [str(i) for i in missing_pmids]

    # TEMP STOP CONDITION TO PREVENT OVERCOLLECTION.
    # missing_pmids = None if (len(in_db) >= 40) else missing_pmids

    return missing_pmids


# %%
def get_articles(pmids):
    """For each ID in pmids, parse-out metadata from Entrez.   

    Args:
        pmids (list): list of article ids (pmids) in LitVar referencing a given rsid id.

    Returns:
        [list of dicts()]: Each dict contains metadata for each pmid parsed. [{title, journal, authors, citation, date, abstract, ncbi_url, pmid, pmc, doi},...]
    """
    # GET ARTICLES
    fetch = PubMedFetcher()
    pub_dates = [] # [publication_date, pmid, {metadata}]

    for i in pmids:
        try:
            article = fetch.article_by_pmid(i)
            pub_dates.append(
                (
                    article.history["pubmed"],
                    i,
                    {
                        "title": article.title,
                        "journal": f"{article.journal} {article.year} {article.volume} {article.issue}",
                        "authors": article.authors,
                        "citation": article.citation,
                        "date": article.history["pubmed"],
                        "abstract": article.abstract,
                        "ncbi_url": article.url,
                        "pmid": article.pmid,
                        "pmc": article.pmc,
                        "doi": article.doi,
                    },
                )
            )
        except Exception as e:
            msg = f"{fn.date_time()} - Article retrieval failed for pmid: {i}. Function: get_articles. {e}"
            print(msg)
            errors.append(msg)
            continue

    if not pub_dates:
        return None  # if list is empty, returns False
    pubs = sorted(pub_dates, reverse=True)
    pubs = [i[2] for i in pubs]

    return pubs


def add_article_bodies(pubs):
    """For each article in pubs list, parse-out the PDF body.

    Args:
        pubs (list): list of articles by id (pmids) for a given rsid.

    Returns:
        [list of dicts]: each dict countains source url and body text [{pdf_url (source), body (text)},...]
    """
    # GET BODY TEXT
    for pub in pubs:
        src = FindIt(pub["pmid"])
        if src.url:
            article_dict = scipdf.parse_pdf_to_dict(src.url)
# TODO CHANGE FROM TEST-ARTICLE TO REAL!!!

            j = json.dumps(article_dict["sections"], indent=4)
            j = json.loads(j)
            body = "\n\n ".join([idx["heading"] + " \n\n" + idx["text"] for idx in j])
            pub.update({"pdf_url": src.url, "body": body})

    return pubs


"""
def add_citedby2(pubs):

    for pub in pubs:
        #Entrez.email = "aris.grout@gmail.com"  # Always tell NCBI who you are
        pmid = pub['pmid']
        results = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs", id=pmid))
        pmc_ids = [link["Id"] for link in results[0]["LinkSetDb"][0]["Link"]]

        # get PMIDs instead of PMCids
        results2 = Entrez.read(Entrez.elink(dbfrom="pmc", db="pubmed", LinkName="pmc_pubmed", id=",".join(pmc_ids)))
        pubmed_ids = [link["Id"] for link in results2[0]["LinkSetDb"][0]["Link"]]
        pub.update({'citedby_pmids': pubmed_ids, 'citedby': len(pubmed_ids)})

    return pubs
"""


def add_citedby(pubs):
    """Count number of citations per article in list of articles (pubs). 

    Args:
        pubs (list of dicts): Each dict in this list holds metadata for a single publication (pmid).

    Returns:
        [list of dicts]: add n_citations and pmc_citedby fields to each dict in list. Note, each dict contains metadata for a single publication (pmid). 
    """
    for pub in pubs:
        try:
            d = pp.parse_citation_web(pub["pmid"], id_type="PMID", api_key=NCBI_API_KEY)
            pub.update({"n_citedby": d["n_citations"], "pmc_citedby": d["pmc_cited"]})
        except Exception as e:
            msg = f"{fn.date_time()} - RSID: {rsid}. Citing PMC ids not found. PMID: {pub['pmid']}. Function: add_citedby. {e}"
            print(msg)
            errors.append(msg)
            pass

    return pubs


def dict_to_df(dictionary):
    """convert article metadata from dict() to DataFrame() for sending to SNP_db/rsid_pubs table. 

    Args:
        dictionary (dict): list of dicts

    Returns:
        [type]: [description]
    """    

    df = pd.DataFrame.from_dict(dictionary)
    df["authors"] = [",".join(authors) for authors in df.authors]
    df["pmc_citedby"] = [",".join(pmc) for pmc in df.pmc_citedby]
    df = df.drop_duplicates()

    return df


def check_pmids(rsid, missing=False):
    """check DB for missing articles (by pmid id).

    Args:
        rsid (str): 'rs12345'
        missing (bool, optional): flag indicating whether rsid was found in DB yet or not. Defaults to False.

    Returns:
        (list): list of article ids (by pmid) that are missing from the DB. 
    """
    # FOR GIVEN RSID:
    try:
        # GET ALL PMIDS
        rsid, pmids = get_pmids(rsid)

        # IF RSID EXISTS IN DB, ONLY QUERY FOR MISSING PMIDS
        if missing is False:
            pmids = missing_pmids(rsid, pmids)

        # RETURN PMIDS DEPENDING ON COUNT
        if (pmids is None) and (missing is True):
            return [{"rsid": rsid}]  # blank entry if no articles exist for rsid
        elif pmids is None:
            return None  # rsid in db, but no new articles to add
        else:
            return pmids  # return list of pmids for collection

    except Exception as e:
        msg = f"{fn.date_time()} - RSID: {rsid}. Unknown error. Function: check_pmids. {e}"
        print(msg)
        errors.append(msg)
        pass

    return pmids

def hundred_chunks(lst: list):
    """Generator: yield successive 100-sized chunks from lst"""
    for i in range(0, len(lst), 100):
        yield lst[i : i + 100]

def scrape_articles(rsid, pmids):
    """Super-function that collects all metadata for each article in list of rsid pmids.

    Args:
        pmids (list): list of all articles (by pmid) missing from SNP_db/rsid_pubs for given rsid (per LitVar).

    Returns:
        [list of dicts]: each dict countains all metadata for a single article (pmid) of all those available for the given rsid.
    """
    try:
        if type(pmids[0]) == dict:
            return pmids  # not actual pmid [{rsid: rs}]
        chunk = 1
        total_chunks = round((len(pmids) / 100) + 1)
        pubs = []
        for hun_pmids in hundred_chunks(lst=pmids):
            # Print progress
            print(f"chunk: {chunk} / {total_chunks}: {round(chunk/total_chunks * 100, 2)}%")
            # Per 100 PMIDs, get the articles for the date and n_citedby
            new_pubs = get_articles(hun_pmids)
            if new_pubs is None:
                continue
            new_pubs = add_citedby(new_pubs)
            # Add new pubs to pubs
            pubs.extend(new_pubs)
            # Get the most recent 20 pubs of the combined list
            most_recent = [pub["pmid"] for pub in sorted(pubs, key=lambda x: x["date"], reverse=True)]
            # Get the most cited 20 pubs of the combined list
            most_cited = [pub["pmid"] for pub in sorted(pubs, key=lambda x: x.get("n_citedby", 0), reverse=True)]
            # Get the top 20 most recent and the top 20 most cited and remove duplicates
            most_recent_and_most_cited = set(most_recent[:20]).union(set(most_cited[:20]))
            # Filter pubs to include the most recent 10 and most cited 10
            pubs = [pub for pub in pubs if pub["pmid"] in list(most_recent_and_most_cited)]
            # Increment counter
            chunk += 1

        # Compare with the db, and find which PMIDs are missing
        check_pmids = [pub['pmid'] for pub in pubs]
        need_to_scrape_pmids = missing_pmids(rsid=rsid, pmids=check_pmids)
        # No PMIDs are missing for RSID, return
        if need_to_scrape_pmids is None:
            return
        # Only scrape article bodies of the most_recent_and_most_cited PMIDs
        pubs = [pub for pub in pubs if pub["pmid"] in need_to_scrape_pmids]
        pubs = add_article_bodies(pubs)
    except Exception as e:
        msg = f"{fn.date_time()} - RSID: {rsid}. Unknown error. Function: scrape_articles: {e}"
        print(msg)
        errors.append(msg)
        pass

    print(f"{len(pubs)} articles added:")
    count = 1
    for title in [i["title"] for i in pubs]:
        print(count, title)
        count += 1

    for pub in pubs:
        pub.update({"rsid": rsid})
    
    return pubs

# %%

# // Import Variables
with open("./data/temp/missing_rsids.v", "rb") as f:
    rsids = pickle.load(f)
errors = []
cnx = fn.create_connection("data/SNP_db.sqlite")

# TESTING
#rsids, start_rsid, start_idx, og_rsids = fn.test_from_rsid("rs72551375")
cols = ["rsid", "date", "title", "journal", "authors", "abstract", "body", "pmid", "pmc", "doi", "ncbi_url", "pdf_url", "citation", "n_citedby", "pmc_citedby"]

# %%
# EXECUTION & PRINT CHECKPOINTS
for rsid in rsids:

    rsid_missing = pd.read_sql_query(f"SELECT EXISTS(SELECT 1 FROM rsid_pubs WHERE rsid='{rsid}');", cnx)
    rsid_missing = True if (rsid_missing.iloc[0, :][0] == 0) else False

    pmids = check_pmids(rsid, rsid_missing)
    if pmids is None:
        continue  # move to next rsid if no new articles to add

    print(f"RSID: {rsid}")
    print(f"{len(pmids)} new articles found.")

    result = scrape_articles(rsid, pmids)
    if result is None:
        continue

    for pub in result:
        missing_cols = set(cols).difference(pub)
        for col in missing_cols:
            pub[col] = "" if col in ["pmc_citedby", "authors"] else None
    df = dict_to_df(result)
    fn.send_data_to_db(df, "rsid_pubs", db_path="data/SNP_db.sqlite", if_exist="append", mode=None)

    current = rsids.index(rsid) + 1
    print(f"rsid: {current} / {len(rsids)}: {round(current/len(rsids) * 100, 2)}%")
    print("\n\n")

    # write errors to file ################################################## FWD ERROR LOG ##########
    with open("data/error_log_fwd.txt", "a") as f:
        for e in errors:
            f.write(e + "\n\n")
    errors.clear()
# %%
src = FindIt('31548494')
# %%
src.url
# %%
in_db = pd.read_sql_query(f"SELECT pmid FROM rsid_pubs WHERE rsid='{'rs11676382'}';", cnx)
# %%
in_db = pd.read_sql_query(f"SELECT pmid FROM rsid_pubs WHERE rsid='{'rs6152'}';", cnx)
# %%
in_db = [] if len(in_db) == 0 else in_db.pmid.apply(int)
pmids = pd.Series(pmids, dtype=int)
missing_pmids = list(set(pmids).difference(in_db))
missing_pmids = None if (len(missing_pmids) == 0) else [str(i) for i in missing_pmids]

# %%
set(pmids).difference([])
# %%
missing_pmids
# %%
in_db.dtypes
# %%
pmids
# %%
missing_pmids
# %%
def missing_pmids(rsid, pmids):
    """Check DB, return list of missing pmids (article ids) for given rsid

    Args:
        rsid (str): rs12345
        pmids (list): list of all pmids for given rsid; pulled from LitVar

    Returns:
        list: list of pmids not yet in DB
    """    

    cnx = fn.create_connection("data/SNP_db.sqlite")
    in_db = pd.read_sql_query(f"SELECT pmid FROM rsid_pubs WHERE rsid='{rsid}';", cnx)
    in_db = [] if (len(in_db) == 0) else in_db.pmid.apply(int)
    pmids = pd.Series(pmids, dtype=int)
    missing_pmids = list(set(pmids).difference(in_db))
    missing_pmids = None if (len(missing_pmids) == 0) else [str(i) for i in missing_pmids]

    # TEMP STOP CONDITION TO PREVENT OVERCOLLECTION.
    # missing_pmids = None if (len(in_db) >= 40) else missing_pmids

    return missing_pmids

# %%
check_pmids = ['30835568', '28342179', '22959728']
pmids = None #[22959728, 28342179, 30835568]
rsid = 'rs3011225'
# %%
need_to_scrape_pmids = missing_pmids(rsid='rs12681963', pmids=check_pmids)
# %%
pmids = None
pmids = pd.Series(pmids, dtype=int)
# %%
pmids
# %%
print(need_to_scrape_pmids)
# %%
'''
2021-08-04 21:57:18 - RSID: rs10509373. Unknown error. Function: scrape_articles: HTTPConnectionPool(host='localhost', port=8070): Max retries exceeded with url: /api/processFulltextDocument (Caused by NewConnectionError('<urllib3.connection.HTTPConnection object at 0x7fd6fc099fd0>: Failed to establish a new connection: [Errno 111] Connection refused'))
''' # TODO resolve this issue.