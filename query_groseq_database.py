#!/usr/bin/env python
# query_groseq_database.py -- A script to query the data fetched by make_groseq_database.py
# Note: this file is getting large enough that it may be useful to split it up into multiple modules.

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
    
    # Generate and sort a list of species:
    newSpeciesList=[(curSpeciesCount[i], val) for i, val in enumerate(curSpeciesList)]
    newSpeciesList.sort()
    newSpeciesList.reverse()
    
    print("List of available species:")
    for elem in newSpeciesList:
        print("%d %s" % (elem[0], elem[1]))
    
    print("\nTotal number of elements: %d" % sum(curSpeciesCount))


def findContribNameDate(matFile):
    firstContrib=None
    pubYear=None
    
    f=open(matFile, "r")
    eTree=curTree=ETree.parse(f)
    f.close()
    
    root=eTree.getroot()
    
    for elem in root:
        # Tag, text, attrib.
        if elem.tag=="{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Contributor" and firstContrib is None:
            if elem.attrib['iid']=="contrib1":
                # This is exactly as painful as it looks:
                for child in elem:
                    if child.tag=="{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Person":
                        for nameChunk in child:
                            if nameChunk.tag=="{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Last":
                                firstContrib=nameChunk.text
        elif elem.tag=="{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Series": # and pubYear is None:
            for child in elem:
                if child.tag=="{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Status":
                    for date in child:
                        if date.tag=="{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Release-Date":
                            # Dates are formatted as yyyy-mm-dd
                            newPubYear=date.text.split('-')[0]
                            
                            if pubYear is None:
                                pubYear=newPubYear
                            elif int(pubYear)<int(newPubYear):
                                pubYear=newPubYear
        
        # Terminate the loop early if we've found both elements:
        if firstContrib is not None and pubYear is not None:
            break
    
    if firstContrib is None or pubYear is None:
        return None
    
    else:
        return firstContrib+pubYear


def findContribByPaper(pathList, paperName):
    adjustedPaperName=paperName.strip().upper()
    for p in pathList:
        matDir=os.path.join(p, "matrices")
        ncName=os.path.join(p, "namecache.txt")
        if os.path.exists(ncName):
            with open(ncName, "r") as ncFile:
                # Be as forgiving as possible:
                if ncFile.read().strip().upper()==adjustedPaperName:
                    return p.split('/')[-1]
        elif os.path.exists(matDir):
            for matFile in os.path.listdir(matDir):
                computedName=findContribNameDate(os.path.join(matDir, matFile))
                if computedName is not None:
                    if computedName.strip().upper()==adjustedPaperName:
                        return p.split('/')[-1]
    return None
    

def printTitle(path):
    pathToks=path.split('/')
    idnum=pathToks[-1]
    protocol=pathToks[-2]
    #idnum=path.split('/')[-1]
    
    matDir=os.path.join(path, "matrices")
    contribName=""
    if os.path.exists(matDir):
        if os.path.exists(os.path.join(path, "namecache.txt")):
            with open(os.path.join(path, "namecache.txt"), "r") as nameFile:
                contribName='"%s" ' % nameFile.read()
        else:
            for mat in os.listdir(matDir):
                pName=findContribNameDate(os.path.join(matDir, mat))
                if pName is not None:
                    pName=pName.strip()
                    with open(os.path.join(path, "namecache.txt"), "w") as nameFile:
                        nameFile.write(pName)
                    
                    contribName='"%s" ' % pName
                    break
                
    f=open(os.path.join(path, "summary.txt"), "r")
    title=f.readline()
    title=title[(len(title.split()[0])+1):]
    f.close()
    print("[%s] %s%s: %s\n" % (protocol, contribName, idnum, title.strip()))


def genProtoSetStr(protocolSet):
    retStr=protocolSet[0]
    for p in protocolSet[1:]:
        retStr+=",%s" % p
    return retStr

        
def findSpecies(pathlist, speciesName, seriesOnly, protocolSet):
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
    
    print("Found %d elements that match \"%s\" given protocol(s) %s" % (len(curFoundList), speciesName, genProtoSetStr(protocolSet)))
    print("List of paths:\n")
    
    for c in curFoundList:
        printTitle(c)


def findProto(basedir, idNum, protocolSet):
    for p in protocolSet:
        if os.path.exists(os.path.join(basedir, p, idNum)):
            return p
    return None


def getSummary(basedir, idlist, pathList, protocolSet):
    for i in idlist:
        try:
            # Fortunately, while this is linear, it's linear relative to the (very small) set of protocols.
            p=findProto(basedir, i, protocolSet)
            
            # Ie. if the value given is a paper name:
            if not (i[0]>='0' and i[0]<='9'):
                paperName=i
                i=findContribByPaper(pathList, paperName)
                if i is None:
                    print("ERROR: couldn't find element matching title: %s" % paperName)
                    continue
            
            curdir=os.path.join(basedir, p, i)
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
                
            print("Protocol: %s" % p)
            print("")
        except:
            print("ERROR: Could not look up %s\n" % i)


def fetchMatrices(basedir, idlist, protocolSet):
    for i in idlist:
        try:
            p=findProto(basedir, i, protocolSet)
            matUrlFile=open(os.path.join(basedir, p, i, "matrixpath.txt"), "r")
            matUrl=matUrlFile.read()
            matUrlFile.close()
            matDir=os.path.join(basedir, p, i, "matrices")
            
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
            
            # Cache the paper's "name":
            f=open(os.path.join(basedir, p, i, "namecache.txt"), "w")
            f.write(findContribNameData(os.path.join(matDir, os.listdir(matDir)[0])).strip())
            f.close()
            
        except:
            print("ERROR: Could not fetch data matrices for %s\n" % i)
            print("It may be possible that the requested data element doesn't have a matrix link.")


def sraHelper(url):
    f=urllib.urlopen(url.strip())
    content=f.read()
    f.close()
    return content


def getSraList(basedir, idlist, pathList, protocolSet):
    for i in idlist:
        p=findProto(basedir, i, protocolSet)
        
        # Ie. if the value given is a paper name:
        if not (i[0]>='0' and i[0]<='9'):
            paperName=i
            i=findContribByPaper(pathList, paperName)
            if i is None:
                print("ERROR: couldn't find element matching title: %s" % paperName)
                continue
        
        print("Finding SRAs for %s..." % i)
        matDir=os.path.join(basedir, p, i, "matrices")
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
            
            outPath=os.path.join(basedir, p, i, "%s.sralist" % i)
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


def download(basedir, elem, outdir, protocolSet):
    # Attempt to open the file:
    p=findProto(basedir, elem, protocolSet)
    datafile=os.path.join(basedir, p, elem, "%s.sralist" % elem)
    
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


def listProtocols(protoList):
    print("Set of protocols currently defined in the database:")
    for p in protoList:
        print("   %s" % p)
    print("Specific protocols can be selected with the -pt parameter.")
    print("Example: -pt=gro-seq,pro-seq")
    

def genIdSet(idSet, basedir, protoName, seriesOnly):
    outSet=[]
    
    for i in idSet:
        if isSeries(os.path.join(basedir, protoName, i)) or not seriesOnly:
            outSet.append(i)
    
    return outSet


def queryProtocol(basedir, protoName, idsByProto, seriesOnly):
    try:
        idset=genIdSet(idsByProto[protoName], basedir, protoName, seriesOnly)
        
        print("Found %d elements matching protocol %s.\n" % (len(idset), protoName))
        
        for i in idset:
            printTitle(os.path.join(basedir, protoName, i))
    except:
        print("ERROR: Query failed for protocol %s. Is it shown by the listprotocols command?" % protoName)


def protocolOverlap(basedir, protocolSet, seriesOnly, idsByProto):
    # First step: determine which elements overlap by hammering the hell out of the idsByProto strucutre:
    # This implementation is... horribly, painfully inefficient.
    # TODO: Make this less awful.
    
    # So, what do we know about this data?
    # 1. The lists are not sorted. Python can sort them in n*log(n) time.
    # 2. There is likely some overlap.
    #    We can mark every element that's been visited. This will reduce complexity by a fair deal, but
    #    not much in an asymptotic sense.
    # 3. ID numbers are unique in that separate projects or elements never map to the same ID.
    # 4. A lot of performance could be gained by not comparing lists of strings vs integers.
    visitedList=[]
    
    for p in protocolSet:
        idsByProto[p].sort()
    
    for p in protocolSet:
        plist=idsByProto[p]
        for q in protocolSet:
            if p==q:
                continue
            qlist=idsByProto[q]
            # Determine which IDs match:
            for r in plist:
                for s in qlist:
                    if r==s:
                        print("[%s] %s <----> [%s] %s" % (p, r, q, s))
                        visitedList.append(r)
    

def main(args):
    progName=args[0]
    seriesOnly=False
    
    protocolSet=None
    
    # This kludge allows for switches to be specified without disrupting any other behavior.
    newArgs=[]
    if len(args)!=1:
        for a in args:
            aToks=a.split("=")
            if a=='-s':
                seriesOnly=True
                #args=args[1:]
            elif aToks[0]=="-pt":
                protocolSet=aToks[1].split(',')
            else:
                newArgs.append(a)
        args=newArgs
                
        
    if len(args)<3:
        print("Usage: %s [-s] dbdir command <args>" % progName)
        print("Query a GRO-Seq metadata database fetched with make_groseq_database.py")
        print("If -s is specified, then only series IDs will be reported on")
        print("If -pt=<comma separated list of protocols> is specified, then only IDs with a specific protocol will be reported on.")
        print("")
        print("List of commands:")
        print("  listprotocols -- List all protocols in the current database.")
        print("  queryprotocol -- Print all elements matching a given protocol.")
        print("  protocoloverlap -- Print out the set of elements that overlap between different protocols.")
        print("  listspecies -- Print a listing of all species defined in the database.")
        print("  findspecies <species name in quotes> -- Find all projects that match a given species.")
        print("  getsummary <list of id numbers or paper names> -- Retrieves a summary for a given data element.")
        print("  fetchmatrices <id> -- Downloads data matrices necessary to fetch data for a set.")
        print("  fetchspmats <species name> -- Fetches all matrices for a given species.")
        print("  fetchallmatrices -- Fetch as many matrix files as possible. This will be slow!")
        print("  getsralist <id or paper name> -- Retrieves all SRAs for a given element given that matrices are present.")
        print("  getreadytosra -- Retrieves a list of all projects with fetched matrix files.")
        print("  getreadytodownload -- Retrieves a list of all projects that can be downloaded immediately.")
        print("  download <id or paper name> <outputdir> -- Downloads data into the specified directory")
        return
    
    if protocolSet is None:
        protocolSet=os.listdir(args[1])
    
    pathlist=[]
    # Store a set of IDs by requested protocol to make certain operations faster and easier.
    idsByProto={'%s' % p: [] for p in protocolSet}
    
    for p in protocolSet:
        modifiedList=[]
        tmpList=os.listdir(os.path.join(args[1], p))
        
        idsByProto[p].extend(tmpList)
        
        for t in tmpList:
            modifiedList.append(os.path.join(args[1], p, t))
        
        pathlist.extend(modifiedList)
    
    if args[2]=="listspecies":
        genSpeciesList(pathlist, seriesOnly)
    
    elif args[2]=="listprotocols":
        listProtocols(protocolSet)
    
    elif args[2]=="queryprotocol":
        queryProtocol(args[1], args[3], idsByProto, seriesOnly)
    
    elif args[2]=="protocoloverlap":
        protocolOverlap(args[1], protocolSet, seriesOnly, idsByProto)
    
    elif args[2]=="findspecies":
        try:
            findSpecies(pathlist, args[3], seriesOnly, protocolSet)
        except:
            print("You must specify a species name.")
            
    elif args[2]=="getsummary":
        try:
            getSummary(args[1], args[3:], pathlist, protocolSet)
        except:
            print("You must specify an element ID or paper name to get a summary.")
        
    elif args[2]=="fetchmatrices":
        try:
            fetchMatrices(args[1], args[3:], protocolSet)
        except:
            print("You must specify an element ID to fetch its data matrices.")
    
    elif args[2]=="fetchspmats":
        try:
            fetchspmats(args[1], pathlist, args[3], seriesOnly)
        except:
            print("You must specify a species name to fetch matrices for that species.")
            
    elif args[2]=="fetchallmatrices":
        fetchMatrices(args[1], curlist)
        
    elif args[2]=="getsralist":
        getSraList(args[1], args[3:], pathlist, protocolSet)
        
    elif args[2]=="getreadytosra":
        getReadyToSra(pathlist)
    
    elif args[2]=="getreadytodownload":
        getReadyToDownload(pathlist)
        
    elif args[2]=="download":
        download(args[1], args[3], protocolSet)
        
    else:
        print("Unknown command: %s" % args[2])


if __name__=="__main__":
    main(sys.argv)