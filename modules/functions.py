#!/usr/bin/env python3
# %%
#  ---------------------- packages

# // RELOAD MODULE OPTION WHEN NEEDED
# importlib.reload(sys.modules['transformers'])

# %%
# General packages
import numpy as np
import pandas as pd
import pickle  # save API responses to variables
import math
import json
import os
import fnmatch
from itertools import chain

# DB handling packages
import sqlite3
from sqlite3 import Error

# Date / Time packages
from datetime import datetime, timedelta
import time
from time import localtime, strftime

# API packages
import requests as re
from mediawiki import MediaWiki  # for calling snpedia MediaWiki API

# API keys
NCBI_API_KEY = os.environ["NCBI_API_KEY"]
# http://biopython.org/DIST/docs/tutorial/Tutorial.html # 9.16.3  Searching for citations
# from Bio import Entrez

# // Adds this file to top-level interpretter path, so runs "as if" from top.
# import sys
# from pathlib import Path

# sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# %%

# // MISC FUNCTIONS ****************************************************


def date_time():
    return strftime("%Y-%m-%d %H:%M:%S", localtime())


def log_error(error, file_path="logs/", message=None):
    """log Exception to error log with optional message.

    Args:
        exception (var, str): output from except statement
        file_path (str): path to error log
        message (str, optional): custom message. Defaults to None.
    """
    msg = f"{date_time()} - {message}. {error}"
    print(msg)

    with open(f"{file_path}", "a") as f:
        f.write(msg + "\n\n")

    return


def find_files(base, pattern):
    """Return list of files matching pattern in base folder."""
    return [n for n in fnmatch.filter(os.listdir(base), pattern) if os.path.isfile(os.path.join(base, n))]


def interval_calc(start, end, intervals, roundup=False):
    """returns numpy array of n 'intervals' between 'start' and 'end'

    Args:
        start (int): the start indices of 'querys'
        end (int): the end indices of 'querys'
        intervals (list/tuple): number of intervals to split 'querys' into
        roundup (boolean): when True, round up to nearest hundred. Default = False.
    """
    result = []

    # round 'end' up by nearest hundred
    if roundup:
        end = int(math.ceil(end / 100) * 100)

    # compute intervals
    interval_size = (end - start) / intervals
    end = start + interval_size
    while True:
        result.append([int(start), int(end)])
        start = end + 1
        end = end + interval_size
        if len(result) == intervals:
            break

    return np.array(result)


# // DATABASE FUNCTIONS ****************************************************


def create_connection(path):
    """create DB connection

    Args:
        path (string): path to SQLite database

    Returns:
        connection (var): sqlite3.connect(path)
    """
    connection = None
    try:
        connection = sqlite3.connect(path)
    except Error as e:
        print(f"The error '{e}' occurred")

    return connection


def execute_query(connection, query):
    cursor = connection.cursor()
    try:
        cursor.execute(query)
        connection.commit()
        print("Query executed successfully")
    except Error as e:
        print(f"The error '{e}' occurred")


def create_db():
    """create SNPit Database"""
    cnx = create_connection("data/SNP_db.sqlite")

    execute_query(
        cnx,
        """
        CREATE TABLE IF NOT EXISTS genotypes (
            id INTEGER NOT NULL UNIQUE PRIMARY KEY AUTOINCREMENT,
            rsid TEXT,
            chromosome TEXT,
            position INTEGER,
            genotype TEXT,
            query TEXT,
            orientation TEXT,
            magnitude TEXT,
            repute TEXT,
            summary TEXT,
            fullurl TEXT,
            date TIMESTAMP
        );
        """,
    )

    execute_query(
        cnx,
        """
        CREATE TABLE IF NOT EXISTS rsid_pubs (
            id INTEGER NOT NULL UNIQUE PRIMARY KEY AUTOINCREMENT,
            rsid TEXT,
            date TIMESTAMP,
            title TEXT,
            journal TEXT,
            authors TEXT,
            abstract TEXT,
            body TEXT,
            pmid TEXT,
            pmc TEXT,
            doi TEXT,
            ncbi_url TEXT,
            pdf_url TEXT,
            citation TEXT,
            n_citedby INTEGER,
            pmc_citedby TEXT
        );
        """,
    )

    execute_query(
        cnx,
        """
        CREATE TABLE IF NOT EXISTS rsid_summaries (
            id INTEGER NOT NULL UNIQUE PRIMARY KEY AUTOINCREMENT,
            rsid TEXT,
            abs_summary TEXT,
            most_cited INTEGER,
            most_recent INTEGER,
            truncated INTEGER,
            pmids TEXT,
            model TEXT,
            run TEXT,
            date DATETIME
        );
        """,
    )


def check_db(query, compare, path):
    """compare a 'list of values' to DB query results.

    Args:
        query (str): sqlite query statement
        compare (list): field values to compare (e.g. df.col0 or [df.col1, df.col2])
        path (str): path to db

    Returns:
        list of lists: values/rows missing from table/query.
    """

    try:
        cnx = sqlite3.connect(path)
        cur = cnx.cursor()
        db = cur.execute(query).fetchall()
        if len(db) == 0:
            db = [""]
        else:
            db = [idx[0] for idx in db] if (len(db[0]) < 2) else db
        if len(compare) > 1:
            missing = set(zip(*compare)).difference(set(db))
            available = set(zip(*compare)).intersection(set(db))
        else:
            missing = set(compare[0]).difference(set(db))
            available = set(compare[0]).intersection(set(db))

        print(f"DB Check: missing: {len(missing)} | available {len(available)}\n")

        return {"missing": list(missing), "available": list(available)}

    except Error as e:
        print(f"The error '{e}' occurred")
        raise


def sample_data(db_path, table, size):
    """Generate random sample from db table: of a particular size

    Args:
        db_path (str): path to database
        table (str): name of table in DB
        size (int): number of rows to sample from DB table

    Returns:
        sample [pd.DataFrame]: DF of random rows, index same as DB table
    """
    cnx = create_connection(db_path)
    sample = pd.read_sql_query(f"SELECT * FROM flights_X ORDER BY RANDOM() LIMIT {size};", cnx, index_col="index")
    return sample


def send_data_to_db(data, table_name, db_path="data/training.sqlite", if_exist="fail", mode="interval"):
    """Send data to training.sqlite database

    Args:
        data (DataFrame): pd.DataFrame, for adding to DB as table_name
        table_name (str): name of table in database
        db_path (str): path to sqlite database. Default 'data/training.sqlite'
        if_exist (str): {‘fail’, ‘replace’, ‘append’}. Default ‘fail’
        mode (str): {'interval', 'None'}. Default 'interval'

    Returns:
        msg (str): confirmation message
    """
    cnx = create_connection(db_path)

    if if_exist == "replace":
        execute_query(cnx, f"DROP TABLE IF EXISTS {table_name}")
        if_exist = "append"

    if mode == "interval":
        intervals = interval_calc(0, len(data) + 1, 30)
        print(intervals)
        [data.iloc[inv[0] : inv[1], :].to_sql(table_name, cnx, if_exists=if_exist) for inv in intervals]

    else:
        data.to_sql(table_name, cnx, if_exists=if_exist, index=False)

    print(f"Saved DF to table: {table_name}, in database: {db_path}")

    return


# // API FUNCTIONS


def all_snpedia_rsids():
    """Get all SNP rsids that exist on SNPedia API

    Returns:
        list: rsid ids
    """
    snpedia = MediaWiki(url="https://bots.snpedia.com/api.php")

    # pull list of all available SNP markers
    all_snps = snpedia.categorymembers("Is a snp", results=None, subcategories=False)
    all_snps = [s.lower() for s in all_snps]  # convert to lowercase

    # save to file
    with open("./data/SNPedia_ids.txt", "wb") as f:
        pickle.dump((all_snps), f)

    return all_snps


def all_chromosome_accessions_23me():
    """Get all accesssion (chromosome) ids from 23andMe API. NOT YET USED.

    Returns:
        list: acession ids
    """
    # Pull accessions (chromosome IDs)
    url = "https://api.23andme.com/3/accession/"
    accessions = re.get(url)
    accessions = [obj["id"] for obj in accessions.json()["data"]]

    # storing to file
    with open("./data/23me_accessions.txt", "wb") as f:
        pickle.dump((accessions), f)

    return accessions


def all_23me_rsids():
    """Get all SNP rsids that exist on 23andME API. NOT YET USED.

    Returns:
        list: rsid ids
    """
    # load accession ids
    with open("./data/23me_accessions.txt", "rb") as f:
        accessions = pickle.load(f)

    # Pull marker ids (by accession id)
    all_markers = []
    count = 0
    end = len(accessions)

    while count < end:

        a = accessions[count]
        url = f"https://api.23andme.com/3/marker/?accession_id={a}"
        markers = re.get(url)

        if markers.status_code == 200:
            ids = [obj["id"] for obj in markers.json()["data"]]
            all_markers.extend(ids)
            count + +1
        elif markers.status_code == 429:
            dt = datetime.now() + timedelta(hours=4)
            while datetime.now() < dt:
                time.sleep(1)
        else:
            print(f"unaccounted for error code: {markers.status_code}")
            break

    with open("./data/23me_ids.txt", "wb") as f:
        pickle.dump((all_markers), f)

    return all_markers


def exclude_indels(df):
    """# INDEL MARKERS are proprietary to 23andME: http://www.enlis.com/blog/2015/10/29/reverse-engineering-23andmes-proprietary-insertions-and-deletions/. exclude for now, will include later.

    Args:
        df (pd.DataFrame): df from join_snpedia

    Returns:
        pd.DataFrame: rsid genotypes with indels removed
    """
    count = len(df[df.genotype.isin(["II", "DD", "DI", "I", "D"])])
    df = df[~df.genotype.isin(["II", "DD", "DI", "I", "D"])]
    print(f"{count} INDELs removed from query DataFrame")

    return df


#  ---  find SNPedia RSID orientations (positive or negative)
#  ---  design SNPedia genotype queries for API calls (reorienting all "minus" genotypes to "positive" equivalents)
#  !!!  REFACTOR: not accurate enough queries. Consider: https://www.snpedia.com/index.php/Orientation
#  !!!  This "quick-swap" isn't accurate to current genotype annotations (yet).
#  !!!  Split into: get_SNPedia_orientations (save DF), design_genotype_queries (save DF)


def design_queries(df):
    """get query strings for collecting genotype data.
    SNP rsid genotypes can be reported in (minus/positive) orientation.
    If minus, then query string needs to be inversed.

    Args:
        df (pd.DataFrame): 23andMe dataframe

    Returns:
        SNPedia genotype query string, orientations: [list], [list]
    """

    q_snpedia = []
    orientations = []

    for index, row in df.iterrows():
        r = row.rsid.lower()
        g = list(row.genotype)

        # GET REQUEST
        url = f"https://bots.snpedia.com/api.php?action=ask&query=[[{r}]][[Category:Is%20a%20snp]]|?Orientation&format=jsonfm&api_version=2"
        response = re.get(url)

        # slice applicable JSON from text response
        start = response.text.find("results") - 1
        stop = response.text.find("serializer") - 11
        result = "{" + response.text[start:stop].lower() + "}"

        # convert: str(JSON) -> dict
        result = json.loads(result)  # dict

        # if minus oriented, flip nucleotide values
        try:
            orient = result["results"][f"{r}"]["printouts"]["orientation"][0]
            orientations.append(orient)
            if orient == "minus":
                for i in range(len(g)):
                    if g[i] == "A":
                        g[i] = "T"
                    elif g[i] == "T":
                        g[i] = "A"
                    elif g[i] == "C":
                        g[i] = "G"
                    elif g[i] == "G":
                        g[i] = "C"
        except Exception:
            print(r)
            orientations.append("missing")
            pass

        # create queries
        if len(g) == 2:
            q_snpedia.append(f"{r}({g[0]};{g[1]})")
        else:
            q_snpedia.append(f"{r}({g[0]})")

    # append data into dataframe
    df["query"] = q_snpedia
    df["orientation"] = orientations

    return df


# //----- pull Mag & Summary per Geno from SNPedia


# // ML FUNCTIONS ****************************************************


def save_model(model, filepath="models/"):
    """Save a trained model to filepath (e.g. 'model/filename')

    Args:
        model (var): variable-held trained model (e.g. Linear_Regression)
        filepath (str): path to save model (excluding file extension)

    Returns:
        msg (str): confirmation message
    """
    pickle.dump(model, open(f"{filepath}.sav", "wb"))
    return f"model saved to: {filepath}.sav"


def load_model(filepath="models/"):
    """Load model into variable

    Args:
        filepath (string): path to model (e.g. models/LinRegression.sav)
    """
    loaded_model = pickle.load(open(filepath, "rb"))
    return loaded_model
