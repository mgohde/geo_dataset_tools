#!/usr/bin/env python
# query_groseq_database.py -- A script to query the data fetched by make_groseq_database.py

import os
import sys
import urllib


def isPresent(reflist, item):
    for i in range(len(reflist)):
        r=reflist[i]
        if r==item:
            return (True, i)
    return (False, 0)


def genSpeciesList(pathlist):
    # This should be O(n^2). Ewwww...
    curSpeciesList=[]
    curSpeciesCount=[]
    for p in pathlist:
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


def findSpecies(pathlist, speciesName):
    curFoundList=[]
    # Making this case-sensitive would just be cruel.
    upperName=speciesName.upper()
    
    for p in pathlist:
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
        print(c)


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
                realName=baseName+".gz"
                os.system("mv %s %s" % (os.path.join(matDir, f), os.path.join(matDir, realName)))
                os.system("gzip -d %s" % os.path.join(matDir, realName))
                
            print("Done with %s." % i)
        except:
            print("ERROR: Could not fetch data matrices for %s\n" % i)
            print("It may be possible that the requested data element doesn't have a matrix link.")


def main(args):
    if len(args)<3:
        print("Usage: %s dbdir command <args>" % args[0])
        print("Query a GRO-Seq metadata database fetched with make_groseq_database.py")
        print("")
        print("List of commands:")
        print("  listspecies -- Print a listing of all species defined in the database.")
        print("  findspecies <species name in quotes> -- Find all projects that match a given species.")
        print("  getsummary <list of id numbers> -- Retrieves a summary for a given data element.")
        print("  fetchmatrices <id> -- Downloads data matrices necessary to fetch data for a set.")
        print("  getsralist <id> -- Retrieves all SRAs for a given element given that matrices are present.")
        
        return
    
    curlist=os.listdir(args[1])
    pathlist=[]
    
    for c in curlist:
        pathlist.append(os.path.join(args[1], c))
    
    if args[2]=="listspecies":
        genSpeciesList(pathlist)
    
    elif args[2]=="findspecies":
        try:
            findSpecies(pathlist, args[3])
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
    
    else:
        print("Unknown command: %s" % args[2])


if __name__=="__main__":
    main(sys.argv)