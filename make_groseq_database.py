#!/usr/bin/env python

# make_groseq_databse.py -- Takes a fetched groseq id file and generates a set of directories containing appropriate
# information. 

import os
import sys
import urllib
import xml.etree.ElementTree as ETree


def genDirs(idlist, dbdir):
    for i in idlist:
        dpath=os.path.join(dbdir, i)
        try:
            os.makedirs(dpath)
        except:
            pass


def genQuery(idlist):
    """This function generates a query string."""
    rStr="db=gds&id=%s" % idlist[0]
    
    for i in idlist[1:]:
        rStr+=",%s" % i
        
    return rStr


def main(args):
    if len(args)!=3:
        print("Usage: %s idfile dbdir")
        print("Fetches project summaries and adds appropriate metadata to a database of GEO GRO-Seq data.")
        print("WARNING: This script may generate several thousand directories under dbdir.")
        return
    
    infile=open(args[1], "r")
    ids=infile.read().splitlines()
    
    genDirs(ids, args[2])
    qstr=genQuery(ids)
    
    # Now let's query the database:
    urlFile=urllib.urlopen("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi", qstr)
    projects=[]
    
    print("Fetching data...")
    searchTree=ETree.parse(urlFile)
    root=searchTree.getroot()
    
    for doc in root:
        # Terrible assumption: the first element in a child should always be its id:
        curid=doc[0].text
        summaryfile=open(os.path.join(args[2], curid, "summary.txt"), "w")
        taxonfile=open(os.path.join(args[2], curid, "taxon.txt"), "w")
        datalistfile=open(os.path.join(args[2], curid, "datalist.txt"), "w")
        relationfile=open(os.path.join(args[2], curid, "relations.txt"), "w")
        
        accession=""
        postedDate=""
        title=u""
        summary=u""
        datalist=[]
        relations=[]
        taxon=u""
        matrixURL=""
        
        for c in doc[1:]:
            curname=c.attrib["Name"]
            if curname=="title":
                title=c.text
            elif curname=="summary":
                summary=c.text
            elif curname=="taxon":
                taxon=c.text
            elif curname=="PDAT":
                postedDate=c.text
            elif curname=="Accession":
                accession=c.text
            elif curname=="Samples":
                for sample in c:
                    # This should be (Title, Accession#)
                    datalist.append((sample[1].text, sample[0].text))
            elif curname=="FTPLink":
                matrixURL=c.text
            elif curname=="ExtRelations":
                for r in c:
                    # This should be SRP#, URL to all the SRA reads.
                    relations.append((r[1].text, r[2].text))
        
        summaryfile.write((u"Title: "+title+u"\n").encode('utf8'))
        summaryfile.write(u"Posted: %s\n" % postedDate)
        summaryfile.write(u"Accession nr: %s\n" % accession)
        summaryfile.write(u"Species: %s\n" % taxon)
        summaryfile.write(u"Matrix URL/FTP Link: %s\n" % matrixURL)
        summaryfile.write(u"Summary (begins on next line):\n")
        summaryfile.write(summary.encode('utf8'))
        summaryfile.close()
        
        taxonfile.write("%s" % taxon)
        taxonfile.close()
        
        for d in datalist:
            # It appears that sample titles are all one "word"
            datalistfile.write((d[0]+u" "+d[1]).encode('utf8'))
        
        datalistfile.close()
        
        for r in relations:
            relationfile.write("%s %s\n" % (r[0], r[1]))
        
        relationfile.close()

if __name__=="__main__":
    main(sys.argv)