#!/usr/bin/env python

# fetch_query.py -- Fetches a list of Enterez lookup IDs corresponding to a
# given query. This tool features a set of default queries corresponding to various 
# *-seq protocols.

# This script, like the others in this directory, is in gross violation of pep8 formatting rules.

import os
import sys
import urllib
import xml.etree.ElementTree as ETree


# Define the esearch URL here in case it is changed in the future:
esearch="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"


def printHelp(progName):
    print("Usage: %s <args> query" % progName)
    print("Query the GEO database for numbers usable by fetch_groseq.py")
    print("If the -o switch is not specified, it is recommended that this program's output")
    print("be piped into either another command or a file.")
    print("")
    print("<args> may be one of the following:")
    print("-h, --help\tPrints this message.")
    print("-o=<filename>\tWrites to an output file instead of stdout.")
    print("-v\tVerbose output. All logging messages are written to stderr.")
    print("")
    print("Examples:")
    print("%s gro-seq >out.txt" % progName)
    print("%s -o=out.txt pro-seq" % progName)


def main(args):
    query="gro-seq"
    out=sys.stdout
    verbose=False
    log=sys.stderr.write
    
    if len(args)==1:
        # Write to stderr so that users can pipe query output into a file or other commands:
       log("[Message] No arguments specified. Defaulting to GRO-seq.\n")
    
    else:
        for a in args[1:]:
            atoks=a.split('=')
            if a=="-h" or a=="--help":
                printHelp(args[0])
                # Ensure that nothing interesting happens after printing the help message:
                return
            elif atoks[0]=="-o":
                try:
                    out=open(atoks[1], "w")
                except:
                    log("Unable to open %s for writing. Exiting...\n")
                    return
            elif a=="-v":
                verbose=True
            else:
                query=a
    
    # Generate a query string:
    argSet={"db": "gds", "term": query, "retmax": 100000}
    params=urllib.urlencode(argSet)
    
    if(verbose):
        log("[Message] About to make query...\n")
    
    urlFile=urllib.urlopen(esearch, params)
    searchTree=ETree.parse(urlFile)
    urlFile.close()
    
    if(verbose):
        log("[Message] Query successful!\n")
    
    root=searchTree.getroot()
    
    for child in root:
        if child.tag=="Count" and verbose:
            log("[Message] Found %s elements.\n" % child.text)
        elif child.tag=="IdList":
            for idEntry in child:
                out.write("%s\n" % idEntry.text)
    
    out.close()
    
    if(verbose):
        log("[Message] Done.\n")
    


if __name__=="__main__":
    main(sys.argv)