#!/usr/bin/env python 
# fetch_groseq.py -- A python script to fetch and maintain a list of GRO-Seq datasets on GEO

import os
import sys
import httplib
import urllib
import xml.etree.ElementTree as ETree

def main(args):
    if len(args)!=2:
        print("Usage: %s outfile" % args[0])
        print("Fetch a list of GRO-Seq datasets present on GEO.")
        return
    
    systemargs={"db": "gds", "term": "gro-seq", "retmax": 100000}
    paramlist=urllib.urlencode(systemargs)
    
    print("Making search request...")
    urlFile=urllib.urlopen("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", paramlist)
    
    print("Fetching data...")
    searchTree=ETree.parse(urlFile)
    urlFile.close()
    
    out=open(args[1], "w")
    
    root=searchTree.getroot()
    for c in root:
        if c.tag=="Count":
            print("Number of elements found: %s" % c.text)
            
        elif c.tag=="IdList":
            for idElem in c:
                out.write("%s\n" % idElem.text)
    
    out.close()
    print("Done.")
            
if __name__=="__main__":
    main(sys.argv)