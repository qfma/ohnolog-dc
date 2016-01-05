#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-
#
#       corum-uniprot-to-ensembl.py
#==============================================================================
import argparse
import sys
import time
import requests
import json
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("corumdb", type=str,
                    help="The corum database file")
parser.add_argument("jsondump", type=str,
                    help="The corum database file")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================


def get_complexes(corumdb):
    ids = {}
    species = []
    with open(corumdb, "rb") as handle:
        handle.next()
        for line in handle:
            line = line.split(";")
            corumid = line[0]
            species.append(line[3])
            # Only use unambiguous complexes
            if not "(" in line[4]:
                subunits = line[4].split(",")
                ids["CORUM"+corumid] = {"uniprot": subunits,
                                        "ensembl": [],
                                        "species": line[3].lower(),
                                        "chicken_orthologs": []}
    return ids, list(set(species))


def ensembl_xref(xref, species):
    payload = {"content-type": "application/json"}
    r = requests.get('http://rest.ensembl.org/xrefs/symbol/'+species+'/'+xref, params=payload)
    return r.json()


def ensembl_homologs(gid):
    payload = {"content-type": "application/json", "format": "condensed"}
    r = requests.get('http://rest.ensembl.org/homology/id/'+gid, params=payload)
    return r.json()

def map_complex_to_ensembl(corum):
    for c in corum:
        pcomplex = corum[c]
        if pcomplex["species"] in ["human", "mouse", "rat", "pig", "rabbit", "dog", "bovine"]:
            for uniprot_id in pcomplex["uniprot"]:
                time.sleep(0.1)
                mapping = ensembl_xref(uniprot_id, pcomplex["species"])
                if len(mapping) >= 1:
                    pcomplex["ensembl"].append(mapping[0]["id"])
                else:
                    print pcomplex, mapping
        print pcomplex
        corum[c] = pcomplex
    return corum


def map_mammalia_to_chicken(corum):
    chicken = "gallus_gallus"
    for c in corum:
        pcomplex = corum[c]
        if pcomplex["ensembl"]:
            for gid in pcomplex["ensembl"]:
                time.sleep(0.1)
                result = ensembl_homologs(gid)
                for homolog in result["data"][0]["homologies"]:
                    if homolog["species"] == chicken:
                        print homolog
                        pcomplex["chicken_orthologs"].append(homolog)
        corum[c] = pcomplex
    return corum


def main():
    corum, species = get_complexes(args.corumdb)
    print len(corum), species
    corum = map_complex_to_ensembl(corum)
    corum = map_mammalia_to_chicken(corum)
    with open(args.jsondump, 'w') as outfile:
        json.dump(corum, outfile)

if __name__ == "__main__":
    main()
