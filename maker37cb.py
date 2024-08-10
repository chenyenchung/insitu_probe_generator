def maker(name,fullseq,amplifier,pause,choose,polyAT,polyCG,polyTriplet,BlastProbes,db,dropout,show,report,maxprobe,numbr): 
    from Bio.Seq import Seq
    from Bio.Blast.Applications import NcbiblastnCommandline as bn
    import io
    import numpy as np
    import pandas as pd
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.max_rows',5000)
    pd.set_option('display.width', 80)
    from datetime import date
    BlastcDNA = 'n'

    __author__ = "Ryan W Null - ORCID_0000-0002-3830-4152"
    __copyright__ = "Copyright 2019-2021  The Ozpolat Lab,  https://bduyguozpolat.org/ "
    __credits__ = ['M. Desmond Ramirez - ORCID_0000-0003-3873-6999','Dennis Sun - ORCID_0000-0003-1000-7276','B. Duygu Ozpolat - ORCID_0000-0002-1900-965X']
    __source__ = "https://github.com/rwnull/insitu_probe_generator"
    __license__ = "GPL 3.0"
    __version__ = "2021_0.3.2"
    __DOI__ = "https://doi.org/10.5281/zenodo.3871970"

    def amp(ampl): 
        if ampl == "B1":
            upspc= "aa"
            dnspc= "ta"
            up = "GAGGAGGGCAGCAAACGG"
            dn = "GAAGAGTCTTCCTTTACG"
        elif ampl == "B2":
            upspc= "aa"
            dnspc= "aa"
            up = "CCTCGTAAATCCTCATCA"
            dn = "ATCATCCAGTAAACCGCC"
        elif ampl == "B3":
            upspc= "tt"
            dnspc= "tt"
            up = "GTCCCTGCCTCTATATCT"
            dn = "CCACTCAACTTTAACCCG"
        elif ampl == "B4":
            upspc= "aa"
            dnspc= "at"
            up = "CCTCAACCTACCTCCAAC"
            dn = "TCTCACCATATTCGCTTC"
        elif ampl == "B5":
            upspc= "aa"
            dnspc= "aa"
            up = "CTCACTCCCAATCTCTAT"
            dn = "CTACCCTACAAATCCAAT"
        elif ampl == "B7":
            upspc= "ww"
            dnspc= "ww"
            up = "CTTCAACCTCCACCTACC"
            dn = "TCCAATCCCTACCCTCAC"
        elif ampl == "B9":
            upspc= "ww"
            dnspc= "ww"
            up = "CACGTATCTACTCCACTC"
            dn = "TCAGCACACTCCCAACCC"
        elif ampl == "B10":
            upspc= "ww"
            dnspc= "ww"
            up = "CCTCAAGATACTCCTCTA"
            dn = "CCTACTCGACTACCCTAG"
        elif ampl == "B11":
            upspc= "ww"
            dnspc= "ww"
            up = "CGCTTAGATATCACTCCT"
            dn = "ACGTCGACCACACTCATC"
        elif ampl == "B13":
            upspc= "ww"
            dnspc= "ww"
            up = "AGGTAACGCCTTCCTGCT"
            dn = "TTATGCTCAACATACAAC"
        elif ampl == "B14":
            upspc= "ww"
            dnspc= "ww"
            up = "AATGTCAATAGCGAGCGA"
            dn = "CCCTATATTTCTGCACAG"
        elif ampl == "B15":
            upspc= "ww"
            dnspc= "ww"
            up = "CAGATTAACACACCACAA"
            dn = "GGTATCTCGAACACTCTC"
        elif ampl == "B17":
            upspc= "ww"
            dnspc= "ww"
            up = "CGATTGTTTGTTGTGGAC"
            dn = "GCATGCTAATCGGATGAG"
        else:
            print ("Please try again")
        return([upspc,dnspc,up,dn])


    def max33(maxprobe,seqs,numbr):
        if maxprobe == 'y':
            if numbr == 0:
                # By default, keep 33 probes
                numbr = 33
                
            if int(numbr) >=  int(len(seqs)):
                print("There was were fewer than "+str(numbr)+" pairs, no action taken.")
                return(seqs)
            if int(numbr) < int(len(seqs)):
                # If there are more than set upper limit number of probes
                # (non-overlapping and without long repeats)
                # Define a binary array to keep track of which probes to keep
                reduced = []
                entry = np.zeros(len(seqs))

                # Otherwise, keep user-provided number of probes
                keep = numbr             # this is the max number of probe pairs that ensures the cheapest opool at 50pmol

                # Define the number of probes to skip
                skip = (len(seqs))-keep
                zeroesperones = int(skip/keep)
                addtnl0s = skip-(keep*zeroesperones)
              
                # a: Count the numbers of probes kept
                a = 0
                # c: Count the probes skipped
                c = 0
                pos = 0
                # 
                # Always keep the first probe
                while (a < keep) & (pos < len(entry)):
                    entry[pos] += 1
                    a += 1
                    pos += 1
                    if c < addtnl0s:
                        # If extra probes need to be skipped
                        # skip one probe per probe kept
                        entry[pos] += 0
                        c += 1
                        pos += 1
                    b = 0
                    while b < zeroesperones:
                        # If there are more probes to skip than kept, than
                        # skip some probes (defined by int(skip/kept)) whenever
                        # there's one probe kept.
                        entry[pos] += 0
                        pos += 1
                        b += 1    
                a = 0
                while a < len(seqs):
                    if entry[a] == 1:
                        reduced.append(seqs[a])
                        a+=1
                    else:
                        a+=1
                        pass     
                return(reduced)
        else:
            return(seqs)

        
    def output(cdna,g,fullseq,count,amplifier,name,pause,seqs):

        amplifier=str((amplifier).upper())
        test=amp(amplifier)
        uspc=test[0]
        dspc=test[1]
        upinit=test[2]
        dninit=test[3]

        if int(count) > 0:
            tab = " "
            print()
            print()
            print("Figure Layout of Probe Sequences:")
            print("")
            print(str(amplifier+"_"+str(name)+"_PP"+str(count)+"_Dla"+str(pause)))
            print("Pair# ","Initiator ","Spacer ","Probe ","Probe ","Spacer ","Initiator")
            i=0
            while i < len(seqs):
                print(str(i+1),tab,upinit,tab,uspc,tab,seqs[i][1][27:52],tab,tab,seqs[i][1][0:25],tab,dspc,tab,dninit)
                i+=1
            print()
            print()
            print()
            print()
            print("Below are the hybridizing sequences and where they align to the cDNA:")
            print()
            print("Pair# ","cDNAcoord ","Probe ","cDNAcoord ","cDNAcoord ","Probe ","cDNAcoord")
            #i=0
            i=len(seqs)-1
            while i >= 0:
                pair = i+1
                coord1 = cdna - int(seqs[i][0])
                coord2 = coord1 - 25
                coord3 = coord2 - 2
                coord4 = cdna - int(seqs[i][2])
                print(pair,tab,coord1,tab,seqs[i][1][0:25],tab,coord2,tab,tab,coord3,tab,seqs[i][1][27:52],tab,coord4)
                i-=1
            print()
            print()
            print()
            print()
            print("This is the in-place localization of the probe pairs along the full-length sense cDNA.")
            print()
            print(">"+name+" Sense Strand")
            print(g)
            print()
            print()
            print()
            print()
            print("Anti-sense sequence used to create probes:")
            print()
            print(">"+name+" Anti-Sense Strand")
            print(fullseq)    
            return()


        
        
        

    #Printing out header
    
    print()
    print()
    print("HCR3.0 Probe Maker Output version "+__version__)
    print(__DOI__)
    print()
    print("Written by "+__author__) 
    print(" with "+__credits__[0])
    print(__credits__[1])
    print(' and '+__credits__[2])
    print()
    print(__copyright__)
    print(" with licensing provided under "+__license__)
    print()
    print("For more information visit: ")
    print(" "+__source__)
    print()
    print()
    print()
    print()
    print()
    print(date.today())
    print()
    print()


    name=str(name)
    
    # Convert input sequence with bio.Seq()
    # Remove spaces, newlines, and tabs
    fullseq = fullseq.replace(" ", "")
    fullseq = fullseq.replace("\n", "")
    fullseq = fullseq.replace("\t", "")
    fullseq = Seq(fullseq)

    # Get reverse_complement cDNA
    fullseq = fullseq.reverse_complement()
    fullseq = str(fullseq)
    
    # Get length
    cdna = len(fullseq)
    
    # 5' region to AVOID probes
    pause = int(pause)    
    
    
    amplifier=str((amplifier).upper())
    
    # Get amplifier sequence as a class
    # TODO: Generalize this part?
    test=amp(amplifier)
    uspc=test[0]
    dspc=test[1]
    upinit=test[2]
    dninit=test[3]

    
    # Define maximal number of same nucleotide repeat to tolerate
    # i.e., dismiss anything >= tolerance + 1
    hpA = "A"*(polyAT+1)
    hpT = "T"*(polyAT+1)
    hpC = "C"*(polyCG+1)
    hpG = "G"*(polyCG+1)
    hpCAG = "GTC"*(polyTriplet + 1)
    hpCTG = "GAC"*(polyTriplet + 1)
    hpCGG = "GCC"*(polyTriplet + 1)

    # Get effective length (exclude 5' pause)
    position = cdna-pause
    
    # Define a 52-bp sweeping window from position 0 to the end.
    # 52 = 25bp left probe + 2bp space + 25bp right probe
    start = np.arange(0,cdna-52,1)
    end = np.arange(52,cdna,1)
    table = np.vstack([start,end])

    cull = 'y'
    seqs={}

    # Define a list of good positions (list of tuples (start, end)) to design
    # probes
    pos=[]


    # Iterating pointer
    a=0
    if cull == "y":
        # Scan the sequence first for mono-base repeats
        while a < (position-52):
            inputseq = str(fullseq[table[0][a]:table[1][a]])
            if (inputseq.find(hpA) + inputseq.find(hpT) + inputseq.find(hpC) + inputseq.find(hpG) + inputseq.find(hpCAG) + inputseq.find(hpCTG) + inputseq.find(hpCGG) > -7):
                # If long repeats of the same base is found, move the pointer
                # forward and skip this window.
                a += 1
            else:
                # If no long repeat is found, add the current window to the
                # position list
                pos.append([table[0][a],table[1][a]])
                a += 1

        ## Creating the first trace through the sequence looking for max number of probe sequences 
        a = 0

        newlist = []

        # Create lists to keep non-overlapping probe regions that are spaced by
        # at least 2bp
        
        newlista = []
        newlista2 = []

        newlistb = []
        newlistb2 = []

        # Initialize the list with the first window
        strt=pos[0][0]
        stp=pos[0][1]
        newlista2.append([cdna-strt,cdna-stp])
        newlista.append([strt,stp])
        while a < len(pos):
            if pos[a][0] > (stp + 2):
                # Only keep a position if it is spaced with the previous
                # region by at least 2bp    
                strt = pos[a][0]
                stp  = pos[a][1]
                newlista2.append([cdna-strt,cdna-stp])
                newlista.append([strt,stp])
                a+=1
            else :
                # Otherwise skip
                a+=1
        lists = {}
        listz = {}
        listz[0] = newlista2
        lists[0] = newlista
                
        if choose == 'y':
            # If the user choose to select from multiple sets of probes
            # create a recursive search by using each window as the start (
            # as opposed to using the first window only) for overlap removal
            # and look for a path that results in the greatest number of probe
            # sequences given the cull
            b = 1
            for x in np.arange(1,len(pos),1):   
                c=0
                strt=pos[x][0]
                stp=pos[x][1]
                newlistb = []
                newlistb2 = []
                newlistb2.append([cdna-strt,cdna-stp])
                newlistb.append([strt,stp])
                while c < len(pos):
                    if pos[c][0] > (stp + 2):
                        strt = pos[c][0]
                        stp  = pos[c][1]
                        newlistb.append([strt,stp])
                        newlistb2.append([cdna-strt,cdna-stp])
                        c+=1
                    else :
                        c+=1
                if len(newlistb2) >= len(newlista2):
                    lists[b] = newlistb
                    listz[b] = newlistb2
                    b+=1
                else:
                    b+=1

            print()
            print("The numbered groups below are your choices for probe sets.")
            print("  The numbers in the columns represent the probe start and stop positions along the parent cDNA.")
            print()
            print((pd.DataFrame(listz)))    ## Returns a comprehensive matrix for all of the longest possibilities, If user wants to define the probeset to use, this can be prompted to show
            print()
            print()
            
            choice = int(input('Enter the number of the probe group with which you wish to proceed. '))
            newlist = np.array(lists[choice])
        else:
            choice = "Default"
            newlist = np.array(lists[0])
    
        # Generate a probe mask that are all N's that has the same length
        # as the cDNA provided
        graphic = ['n']*cdna
        
        count = str(len(newlist))
        
        
        
        pd.set_option('display.width', 1000) 
        
        
        
        if BlastProbes == 'n':
            newlist = max33(maxprobe,newlist,numbr)
            count = str(len(newlist))
            print()
            print("From the given parameters, we were able to make "+count+" probe pairs.")
            print()
            print()
            print("Below is in IDT oPool submission_format.")
            print("Copy and Paste the lines below into an XLSX file for submission to IDT starting from 'Pool name'.")
            print()
            print("Pool name, Sequence")
            a=0
            while a < len(newlist):
                seqs[a] = [newlist[a][0],str(fullseq[newlist[a][0]:(newlist[a][0]+25)]+"nn"+fullseq[(newlist[a][0]+27):newlist[a][1]]),newlist[a][1]]
                graphic[newlist[a][0]:newlist[a][1]] = str(fullseq[newlist[a][0]:(newlist[a][0]+25)]+"nn"+fullseq[(newlist[a][0]+27):newlist[a][1]])
                print(str(amplifier+'_'+name+'_'+count+'_Dla'+str(pause)+','+upinit+uspc+str(seqs[a][1][27:52])))
                print(str(amplifier+'_'+name+'_'+count+'_Dla'+str(pause)+','+str(seqs[a][1][0:25])+dspc+dninit))
                a+=1
            g = ''
            g = g.join(graphic)
            g = Seq(g)
            g = g.reverse_complement()
        else:
            graphic = ['n']*cdna





    ## THE FOLLOWING SECTION CREATES A FASTA FILE FROM THE POTENTIAL PROBE SEQUENCES (BOTH 25BP PROBES COUPLED AS A SINGLE 52BP SEQUENCE INCLUDING A 2BP "nn" SPACER)        
        ## THE RESULTANT FASTA FILE IS BLASTED AGAINST THE USER SPECIFIED TRANSCRIPTOME FASTA 
        ## PROBES THAT MATCH A SEQUENCE IN BLAST WITH A LENGTH MATCH, 60BP > X > 40BP, AND AN E-VALUE < 1E-15 ARE KEPT, OTHERS ARE DISCARDED
    
    
        if BlastProbes == "y":
            print()
            print("BLASTn of probes in progress, this may take a few minutes.")
            print()
            seqs={} 
            remove = pd.DataFrame(columns = ["pos1","seq","pos2","fasta","num"])
            a=0
            tmpFA = open((str(name)+"PrelimProbes.fa"), "w")
            while a < len(newlist):
                nm = str('>'+str(a))
                seqs[a] = [newlist[a][0],str(fullseq[newlist[a][0]:(newlist[a][0]+25)]+"nn"+fullseq[(newlist[a][0]+27):newlist[a][1]]),newlist[a][1],nm,a]
                # remove = remove.append({'pos1' : newlist[a][0], 'seq' : str(fullseq[newlist[a][0]:(newlist[a][0]+25)]+"nn"+fullseq[(newlist[a][0]+27):newlist[a][1]]), 'pos2' : newlist[a][1], 'fasta':nm, 'num':a},  
                # ignore_index = True) 
                new_row = {'pos1' : newlist[a][0], 'seq' : str(fullseq[newlist[a][0]:(newlist[a][0]+25)]+"nn"+fullseq[(newlist[a][0]+27):newlist[a][1]]), 'pos2' : newlist[a][1], 'fasta':nm, 'num':a}
                remove = pd.concat(
                    [remove, pd.DataFrame([new_row])],  
                    ignore_index = True
                ) 
                tmpFA.write(nm)
                tmpFA.write('\n')
                tmpFA.write(seqs[a][1])
                tmpFA.write('\n')
                a+=1
            tmpFA.close()
            
            
            
            
        ## Probe BLAST setup and execution from FASTA file prepared in previous step

            cline = bn(query = str(name)+"PrelimProbes.fa", subject = db, outfmt = 6, task = 'blastn-short') #this uses biopython's blastn formatting function and creates a commandline compatible command 
            stdout, stderr = cline() #cline() calls the string as a command and passes it to the command line, outputting the blast results to one variable and errors to the other

            ## From results of blast creating a numpy array (and Pandas database)
            dt = [(np.unicode_,8),(np.unicode_,40),(np.int32),(np.int32),(np.int32),(np.int32),(np.int32),(np.int32),(np.int32),(np.int32),(np.float64),(np.float64)]
            blastresult = (np.genfromtxt(io.StringIO(stdout),delimiter = '\t',dtype = dt))# "qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore")
           


            ## This loop takes the data from the blast result and filters out probe pairs that do not meet criteria
                ## by setting a length match requirement this eliminates off-target pairs and half-pairs
                ## the e-value threshold ensures that the probe is a good match to the target

            i=0
            filterblast = []
            filterblastbad = []
            uniques = []
            uniquesbad = []
            filterblastexo = []
            uniquesexo = []
            while i < len(blastresult):
                if (blastresult[i][11]>=75.0 and blastresult[i][10]<=float(1e-13)):  #abs(blastresult[i][9]-blastresult[i][8])>40 and abs(blastresult[i][9]-blastresult[i][8])<=60
                    ## Good matches
                    filterblast.append(blastresult[i])
                    uniques.append((blastresult[i][0])) #str
                    i+=1
                elif (blastresult[i][11]<30.0 and blastresult[i][10]>float(0)):
                    ## Likely not in the genome (GFP & etc...)
                    filterblastexo.append(blastresult[i])
                    uniquesexo.append((blastresult[i][0])) #str
                    i+=1
                elif (blastresult[i][11]>=60.0 and blastresult[i][10]>float(1e-12)):
                    ## Bad matches
                    filterblastbad.append(blastresult[i])
                    uniquesbad.append((blastresult[i][0])) #str
                    i+=1
                else:
                    i+=1
            

        
            if len(filterblast) + len(filterblastexo) != 0:
                filterblast = np.array(filterblast)
                uniques = np.unique((uniques))     ##
                count = str((len(uniques)))
                filterblastbad = np.array(filterblastbad)
                uniquesbad = np.unique((uniquesbad))
                filterblastexo = np.array(filterblastexo) 
                uniquesexo = np.unique((uniquesexo))
                if len(uniquesbad) > 0:
                    if dropout == 'y':
                        ind = 0
                        for a in uniquesbad:
                            remove = remove.drop(remove.index[int(a)-ind])
                            ind += 1
                        remove = remove.reset_index()
                        
                print()
                print()
                print("Probe pairs that had possible off-target matches to the provided database (lower e-value but with high site coverage). ")
                print("   Pairs ")
                print(uniquesbad)  
                print()
                print("Probe pairs that had good matches to the provided database were determined to be the following.")
                print("   Pairs ")
                print(uniques)
                print()
                print()
                if show == 'y':
                    print()
                    print()
                    print("The list below shows potentially problematic binding between a probe and some sequence in the provided database.  ")
                    print("    Consider checking if the resulting subject ID is a cDNA you intend to amplify. ")
                    print()
                    print("index qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue qcovs")
                    print()
                    print(pd.DataFrame(filterblastbad))
                    print()
                    print()
                    print()
                    print()
                    print()
                    print("This is a detailed look at the probes with good matches.")
                    print()    
                    print("index qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue qcovs")
                    print()
                    print(pd.DataFrame(filterblast))
                    print()

                if BlastcDNA == 'n':   
                    if len(uniquesbad) > 0:
                        if dropout == 'y':
                            remove = remove.to_dict()
                            a=0
                            seqs={}
                            while a <len(remove["pos1"]):
                                seqs[a] = (remove["pos1"][a],remove["seq"][a],remove["pos2"][a],remove["fasta"][a],remove["num"][a])
                                a += 1
                    
                    seqs = max33(maxprobe,seqs,numbr)
                    count = str(len(seqs))
                    print()
                    print()
                    print()
                    print()
                    print("The probes are provided below in IDT oPool submission format.")
                    print("Copy and Paste the lines below into an XLSX file for submission to IDT starting from 'Pool name'.")
                    print()
                    print()
                    print("Pool name, Sequence")
                    a=0
                    b=0
                    seqs1={}
                    while a < len(seqs):
                    #while a < len(uniques):
                        tmp = (seqs[a])
                        print(str(amplifier+'_'+name+'_'+count+'_Dla'+str(pause)+','+upinit+uspc+str(tmp[1][27:52])))
                        print(str(amplifier+'_'+name+'_'+count+'_Dla'+str(pause)+','+str(tmp[1][0:25])+dspc+dninit))
                        graphic[tmp[0]:tmp[2]] = str(tmp[1])
                        seqs1[b]=tmp
                        b+=1
                        a+=1
                    seqs=seqs1
                    
                    
                    g = ''
                    g = g.join(graphic) 
                    g = Seq(g)
                    g = g.reverse_complement()
            else:
                print()
                print()
                print("Hmm.... There were no probes that fit the parameters specified within the FASTA file.   ")
                print()
                print()
                print("  Try increasing the length of homopolymers tolerated, or BLAST against a different FASTA file.")
                print()
                print("  If BLASTing a heterologous sequence, i.e. GFP, this error could be because the RNA doesn't exist in your species. ")
                print()
                print()


    
    
    
    

    


    
    output(cdna,g,fullseq,count,amplifier,name,pause,seqs)
    
    
        
    print()
    print()
    if report == 'y':
        print("Run "+str(date.today())+"\n   with settings: \n\t5'Pause:\t"+str(pause)+" \n\tChoice of probe set:\t"+str((choose))+"\tPair used: "+str(choice)+" \n\tLength of acceptable polyA/polyT runs:\t"+str(polyAT)+" \n\tLength of acceptable polyC/polyG runs:\t"+str(polyCG)+" \n\tNumber of acceptable triplets (CAG/CTG/CGG) runs:\t"+str(polyTriplet)+" \n\tBLASTn of Probes:\t"+str((BlastProbes))+" \n\tRemoval of probes with low quality BLAST hits:\t"+str((dropout)) )

    print()
    print()
    print("References: ")
    print()
    print(" See Choi et al. 2018 Development for HCR3.0 methodology details ")
    print(" https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6031405/ ")
    print()
    print(" See Wang et al. 2020 BioRxiv for expanded amplifier information ")
    print(' https://www.biorxiv.org/content/10.1101/274456v3 ')