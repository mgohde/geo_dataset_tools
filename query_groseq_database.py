#!/usr/bin/env python
# query_groseq_database.py -- A script to query the data fetched by make_groseq_database.py

import os
import sys
import urllib
from ftplib import FTP
import xml.etree.ElementTree as ETree
from multiprocessing import Pool


def isPresent(reflist, item):
    for i in range(len(reflist)):
        r=reflist[i]
        if r==item:
            return (True, i)
    return (False, 0)


def isSeries(pathToElem):
    typefile=open(os.path.join(pathToElem, "type.txt"))
    contents=typefile.read()
    typefile.close()
    
    return contents=="GSE"


def genSpeciesList(pathlist, seriesOnly):
    # This should be O(n^2). Ewwww...
    curSpeciesList=[]
    curSpeciesCount=[]
    for p in pathlist:
        if isSeries(p) or not seriesOnly:
            taxon=open(os.path.join(p, "taxon.txt"), "r")
            contents=taxon.read().split(';')
            
            for c in contents:
                c=c.strip()
                (truth, idx)=isPresent(curSpeciesList, c)
                if not truth:
                    curSpeciesList.append(c)
                    curSpeciesCount.append(1)
                else:
                    curSpeciesCount[idx]+=1
            taxon.close()
    
    print("List of available species:")
    i=0
    for c in curSpeciesList:
        print("%d %s" % (curSpeciesCount[i], c))
        i+=1
    print("\nTotal number of elements: %d" % sum(curSpeciesCount))


def printTitle(path):
    idnum=path.split('/')[-1]
    f=open(os.path.join(path, "summary.txt"), "r")
    title=f.readline()
    title=title[(len(title.split()[0])+1):]
    f.close()
    print("%s: %s" % (idnum, title.strip()))


def findSpecies(pathlist, speciesName, seriesOnly):
    curFoundList=[]
    # Making this case-sensitive would just be cruel.
    upperName=speciesName.upper()
    
    for p in pathlist:
        if isSeries(p) or not seriesOnly:
            taxon=open(os.path.join(p, "taxon.txt"), "r")
            contents=taxon.read().split(';')
            
            for c in contents:
                c=c.strip()
                if c.upper()==upperName:
                    #curFoundList.append(p.split('/')[-1])
                    curFoundList.append(p)
            taxon.close()
    
    print("Found %d elements that match \"%s\". List of paths:" % (len(curFoundList), speciesName))
    
    for c in curFoundList:
        printTitle(c)


def getSummary(basedir, idlist):
    for i in idlist:
        try:
            curdir=os.path.join(basedir, i)
            f=open(os.path.join(curdir, "summary.txt"), "r")
            contents=f.read()
            f.close()
            
            print("----Summary for %s----" % i)
            print(contents)
            
            print("")
            
            if os.path.exists(os.path.join(curdir, "matrixpath.txt")):
                print("Data matrix URL defined? YES")
            else:
                print("Data matrix URL defined? NO. This software cannot fetch data for this element.")
            
            if os.path.exists(os.path.join(curdir, "matrices")):
                print("Ready to fetch data? YES")
                # TODO: Add additional information as appropriate based on the matrix files.
            else:
                print("Ready to fetch data? NO (run fetchmatrices %s)" % i)
            
            datapath=os.path.join(curdir, "data")
            if os.path.exists(datapath):
                print("Data fetched? YES in %s" % datapath)
            else:
                print("Data fetched? NO")
            print("")
        except:
            print("ERROR: Could not look up %s\n" % i)


def fetchMatrices(basedir, idlist):
    for i in idlist:
        try:
            matUrlFile=open(os.path.join(basedir, i, "matrixpath.txt"), "r")
            matUrl=matUrlFile.read()
            matUrlFile.close()
            matDir=os.path.join(basedir, i, "matrices")
            
            if not os.path.exists(matDir):
                try:
                    os.makedirs(matDir)
                except:
                    print("ERROR: Could not create matrix storage directory for %s! This is very bad." % i)
                    return
            
            print("Fetching matrix file(s)...")
            os.system("wget -r -nH -nd -np -R .listing -P %s %sminiml/" % (matDir, matUrl.strip()))
            
            print("Unpacking matrix file(s)...")
            for f in os.listdir(matDir):
                baseName=f[:-4]
                # For some reason, all of the files have a tgz extension when they're really just gzip'd archives.
                #realName=baseName+".gz"
                #os.system("mv %s %s" % (os.path.join(matDir, f), os.path.join(matDir, realName)))
                os.system("tar -xzf %s -C %s" % (os.path.join(matDir, f), matDir))
            os.system("rm %s" % os.path.join(matDir, "*.tgz"))
                
            print("Done with %s." % i)
        except:
            print("ERROR: Could not fetch data matrices for %s\n" % i)
            print("It may be possible that the requested data element doesn't have a matrix link.")


def sraHelper(url):
    f=urllib.urlopen(url.strip())
    content=f.read()
    f.close()
    return content


def getSraList(basedir, idlist):
    for i in idlist:
        print("Finding SRAs for %s..." % i)
        matDir=os.path.join(basedir, i, "matrices")
        # if the matrix directory exists, then this element probably has the right data:
        if os.path.exists(matDir):
            sralist=[]
            sraURLlist=[]
            # Read all files in the directory:
            dircontents=os.listdir(matDir)
            # filter out all non-xml files:
            flist=[]
            for d in dircontents:
                if d.split(".")[-1]=="xml":
                    flist.append(d)
            
            for fname in flist:
                treeFile=open(os.path.join(matDir, fname), "r")
                curTree=ETree.parse(treeFile)
                treeFile.close()
                
                # This should be the <MiniML> tag. 
                root=curTree.getroot()
                
                for tag in root:
                    if tag.tag=="{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Sample":
                        for child in tag:
                            if child.tag=="{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Type":
                                if child.text!="SRA":
                                    print("Found non-SRA sample.")
                                    break;
                            # We want relation URLs:
                            elif child.tag=="{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Supplementary-Data":
                                if child.attrib["type"]=="SRA Experiment":
                                    sraURLlist.append(child.text)
            # Now that we (hopefully) have a list of urls:
            if len(sraURLlist)==0:
                print("No SRA links found in matrix file. Please check for data manually.")
                return
            
            # We now need to enumerate all SRAs defined in the FTP links:
            print("Fetching data from specified SRA index URLs (this may take a while)...")
            p=Pool(500)
            contentSet=p.map(sraHelper, sraURLlist)
            for c in contentSet:
                lines=c.splitlines()
                
                for l in lines:
                    sralist.append(l.split()[-1].strip())
            
            print("\nFound the following SRA elements for download:")
            for s in sralist:
                print(s)
            
            outPath=os.path.join(basedir, i, "%s.sralist" % i)
            print("\nWriting list of links to %s" % outPath)
            f=open(outPath, "w")
            
            for u in sraURLlist:
                # Each URL already has a newline character in it, for better or for worse.
                f.write("%s\n" % u.strip())
            
            f.close()
                
        else:
            print("ERROR: No matrices defined in %s" % matDir)


def getReadyToSra(pathlist):
    readylist=[]
    for p in pathlist:
        if os.path.exists(os.path.join(p, "matrices")):
            readylist.append(p)
    
    print("%d of %d elements are ready to be fetched. They are:" % (len(readylist), len(pathlist)))
    
    for r in readylist:
        printTitle(r)


def getReadyToDownload(pathlist):
    print("The following elements have SRA list files:")
    numFound=0
    for p in pathlist:
        # Determine if there exists a .sralist file in this directory:
        contents=os.listdir(p)
        
        for c in contents:
            ctoks=c.split('.')
            try:
                if ctoks[1]=="sralist":
                    numFound+=1
                    printTitle(p)
            except:
                # Do nothing
                pass
    
    print("\n%d of %d elements have SRA lists and can be immediately downloaded." % (numFound, len(pathlist)))


def download(basedir, elem, outdir):
    # Attempt to open the file:
    datafile=os.path.join(basedir, elem, "%s.sralist" % elem)
    
    try:
        f=open(datafile, "r")
        lines=f.read().splitlines()
        f.close()
        
        # Each line should represent one URL.
        for u in lines:
            os.system("wget -r -nH -nd -np -R index.html* -P %s %s" % (outdir, u))
    except:
        print("Error: can't open %s" % datafile)


def fetchspmats(basedir, pathlist, speciesname, seriesOnly):
    print("Attempting to fetch all matrices for '%s'..." % speciesname)
    upperName=speciesname.upper()
    
    curFoundList=[]
    # Todo: just write a function to do this lookup so we don't repeat code.
    for p in pathlist:
        if isSeries(p) or not seriesOnly:
            taxon=open(os.path.join(p, "taxon.txt"), "r")
            contents=taxon.read().split(';')
            
            for c in contents:
                c=c.strip()
                if c.upper()==upperName:
                    curFoundList.append(p.split('/')[-1])
                    #curFoundList.append(p)
            taxon.close()
            
    fetchMatrices(basedir, curFoundList)

    
def main(args):
    progName=args[0]
    seriesOnly=False
    
    # This kludge allows for switches to be specified without disrupting any other behavior.
    if len(args)!=1:
        if args[1]=='-s':
            seriesOnly=True
            args=args[1:]
    
    if len(args)<3:
        print("Usage: %s [-s] dbdir command <args>" % progName)
        print("Query a GRO-Seq metadata database fetched with make_groseq_database.py")
        print("If -s is specified, then only series IDs will be reported on")
        print("")
        print("List of commands:")
        print("  listspecies -- Print a listing of all species defined in the database.")
        print("  findspecies <species name in quotes> -- Find all projects that match a given species.")
        print("  getsummary <list of id numbers> -- Retrieves a summary for a given data element.")
        print("  fetchmatrices <id> -- Downloads data matrices necessary to fetch data for a set.")
        print("  fetchspmats <species name> -- Fetches all matrices for a given species.")
        print("  fetchallmatrices -- Fetch as many matrix files as possible. This will be slow!")
        print("  getsralist <id> -- Retrieves all SRAs for a given element given that matrices are present.")
        print("  getreadytosra -- Retrieves a list of all projects with fetched matrix files.")
        print("  getreadytodownload -- Retrieves a list of all projects that can be downloaded immediately.")
        print("  download <id> <outputdir> -- Downloads data into the specified directory")
        return
    
    curlist=os.listdir(args[1])
    pathlist=[]
    
    for c in curlist:
        pathlist.append(os.path.join(args[1], c))
    
    if args[2]=="listspecies":
        genSpeciesList(pathlist, seriesOnly)
    
    elif args[2]=="findspecies":
        try:
            findSpecies(pathlist, args[3], seriesOnly)
        except:
            print("You must specify a species name.")
            
    elif args[2]=="getsummary":
        try:
            getSummary(args[1], args[3:])
        except:
            print("You must specify an element ID to get a summary.")
        
    elif args[2]=="fetchmatrices":
        try:
            fetchMatrices(args[1], args[3:])
        except:
            print("You must specify an element ID to fetch its data matrices.")
    
    elif args[2]=="fetchspmats":
        #try:
            fetchspmats(args[1], pathlist, args[3], seriesOnly)
        #except:
        #    print("You must specify a species name to fetch matrices for that species.")
            
    elif args[2]=="fetchallmatrices":
        fetchMatrices(args[1], curlist)
        
    elif args[2]=="getsralist":
        getSraList(args[1], args[3:])
        
    elif args[2]=="getreadytosra":
        getReadyToSra(pathlist)
    
    elif args[2]=="getreadytodownload":
        getReadyToDownload(pathlist)
        
    else:
        print("Unknown command: %s" % args[2])


if __name__=="__main__":
    main(sys.argv)