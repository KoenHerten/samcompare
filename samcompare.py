#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 11:35:05 2017

This is samcompare v1. A tool to compare the output of multiple read mappers

Copyright 2017, Koen Herten, All rights reserved

This file is part of aftermerge.

samcompare is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

samcompare is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with samcompare.  If not, see <http://www.gnu.org/licenses/>.
"""

from argparse import ArgumentParser
from argparse import FileType
import sys
import samRead
import samFile
   
   
def containsNone(mydict):
    for name in mydict:
        if (mydict[name] is None):
            return True
    return False
    


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

if __name__ == '__main__':
    #add all arguments
    parser = ArgumentParser(description='samcompare v1.0 comparing samfiles of the sam sample on mapping parameters')
    parser.add_argument('-v', '--verbose', action="store_true")
    parser.add_argument('-skip', '--skip', action="store_true", help="Skip reads that are not found in all samples (if false, an error is occuring)")
    parser.add_argument('-ccs', '--ccs', action="store_true", help="These are PacBio CCS reads (fastq description should end with ccs (some mappers add extra numbers, these are removed to get the original query name))")
    parser.add_argument('samfiles', type=FileType('r'), nargs='+')
    
    
    args = parser.parse_args()
    #if verbose, print the parameters
    if args.verbose:
        eprint("Verbose: {}".format(args.verbose))
        eprint("Is CCS reads: {}".format(args.ccs))
        eprint("Skip reads: {}".format(args.skip))
        eprint("Files to parse:")
        for f in args.samfiles:
            eprint("\t{}".format(f.name))
        
            
    
    samfiles = {}
    for f in args.samfiles:
        name = str(f.name).rsplit(".",1)[0]
        samfile = samFile.SamFile(f.name, args.ccs)
        samfiles[name] = samfile
        if args.verbose:
            eprint("Mapping name: {}".format(name))

    count = 0
    samreads = {}
    prevname = ""
    for sfile in samfiles:
        samreads[sfile] = samfiles[sfile].nextPrimaryRead()
        prevname = samreads[sfile].qname

    #initiate the position dict
    position_dict = {}
    for map1 in samfiles:
        p_dict = {}
        for map2 in samfiles:
            t_dict = {}
            t_dict["same"] = 0
            t_dict["diff"] = 0
            t_dict["this_unmap"] = 0
            t_dict["other_unmap"] = 0
            t_dict["both_unmap"] = 0
            p_dict[map2] = t_dict
        position_dict[map1] = p_dict

    #initiate the mapq dict
    mapq_dict = {}
    mapq_range = {}
    for map1 in samfiles:
        p_dict = {}
        for map2 in samfiles:
            t_dict = {}
            t_dict["same"] = 0
            t_dict["higher"] = 0
            t_dict["lower"] = 0
            p_dict[map2] = t_dict
        mapq_dict[map1] = p_dict
        mapq_range[map1] = 0

    #initiate the cigar dict
    cigar_dict = {}
    for map1 in samfiles:
        p_dict = {}
        for map2 in samfiles:
            t_dict = {}
            t_dict["same"] = 0
            t_dict["softclip"] = 0
            t_dict["forcemap"] = 0
            t_dict["notcomp"] = 0
            p_dict[map2] = t_dict
        cigar_dict[map1] = p_dict
        

    
    count = 0
    while (not containsNone(samreads)):
        count = count + 1
        if (args.verbose and count%10==0):
            #print("processed {}".format(count))
            eprint("processed {}".format(count))
        #no none in the samreads, so continue
        
    
        #check if the same name and read
        name = None
        isfirst = True
        same = True
        check = True
        while (check):
            for sam in samreads:
                if (name is None):
                    name = samreads[sam].qname
                    isfirst = samreads[sam].isfirst()
                same = True
                if (name != samreads[sam].qname or isfirst != samreads[sam].isfirst()):
                    same = False
                    if args.skip:
                        if args.verbose:
                            eprint("SKIP")
                        #skip different reads
                        name = max(name, samreads[sam].qname)
                        while (name != samreads[sam].qname and max(name, samreads[sam].qname) == name):
                            #samread before other
                            samreads[sam] = samfiles[sam].nextPrimaryRead()
                    else:
                        eprint("FIRST: {}".format(name))
                        for s  in samreads:
                            eprint("{}\t{}".format(s, samreads[s].qname))
                        eprint("ERROR: not sorted or same readset")
                        exit(1)
            if (same):
                check = False
        #all same read
        
        #check mapping position
        for map1 in samreads:
            read1 = samreads[map1]
            for map2 in samreads:
                read2 = samreads[map2]
                p_dict = position_dict[map1]
                t_dict = p_dict[map2]
                if (read1.ismapped() and read2.ismapped()):
                    #both mapped
                    if (read1.rname == read2.rname and read1.posOfFirstBaseOfRead() == read2.posOfFirstBaseOfRead()):
                        #same position
                        t_dict["same"] = t_dict["same"] + 1
                    else:
                        #different position
                        t_dict["diff"] = t_dict["diff"] + 1
                else:
                    #at least one unmapped
                    if (not read1.ismapped() and not read2.ismapped()):
                        #both unmapped
                        t_dict["both_unmap"] = t_dict["both_unmap"] + 1
                    elif(not read1.ismapped()):
                        #map1 unmapped
                        t_dict["this_unmap"] = t_dict["this_unmap"] + 1
                    else:
                        #map2 unmapped
                        t_dict["other_unmap"] = t_dict["other_unmap"] + 1
        
        #check mapping quality
        for map1 in samreads:
            read1 = samreads[map1]
            if (mapq_range[map1] < read1.mapq):
                mapq_range[map1] = read1.mapq
            for map2 in samreads:
                read2 = samreads[map2]
                p_dict = mapq_dict[map1]
                t_dict = p_dict[map2]
                if (read1.mapq == read2.mapq):
                    #same mapping quality
                    t_dict["same"] = t_dict["same"] + 1
                elif(int(read1.mapq) > int(read2.mapq)):
                    #read1 has higher quality
                    t_dict["higher"] = t_dict["higher"] + 1
                else:
                    #read2 has higher quality
                    t_dict["lower"] = t_dict["lower"] + 1
                
        
        #check cigar
        for map1 in samreads:
            read1 = samreads[map1]
            for map2 in samreads:
                read2 = samreads[map2]
                p_dict = cigar_dict[map1]
                t_dict = p_dict[map2]
                if (read1.ismapped() and read2.ismapped()) and (read1.rname == read2.rname and read1.posOfFirstBaseOfRead() == read2.posOfFirstBaseOfRead()):
                    #both mapped and same location
                    if (read1.cigar == read2.cigar):
                        #same cigar
                        t_dict["same"] = t_dict["same"] + 1
                    else:
                        #different cigar
                        #check the number of matches and indels
                        bases1 = read1.longcigar().count("M") + read1.longcigar().count("D") + read1.longcigar().count("I")
                        bases2 = read2.longcigar().count("M") + read2.longcigar().count("D") + read2.longcigar().count("I")
                        if (bases1 > bases2):
                            #read1 has longer mapping string
                            t_dict["forcemap"] = t_dict["forcemap"] + 1
                        else:
                            #read2 has longer mapping string
                            t_dict["softclip"] = t_dict["softclip"] + 1
                else:
                    t_dict["notcomp"] = t_dict["notcomp"] + 1
                
        for sfile in samfiles:
            samreads[sfile] = samfiles[sfile].nextPrimaryRead()
            
    
            
    print("Statistics:")
    print("Total number of mappers: {}".format(len(samfiles)))
    print("Total number of assessed reads: {}".format(count))
    print("")
   
    for f in samfiles:
        samfiles[f].close()
     
    #print mapping positions
    print("Position statistics")
    print ("Mapper1\tMapper2\t{c1}\t{c2}\t{c3}\t{c4}\t{c5}".format(c1="Same_location", c2="Different_location", c3="Both_unmapped", c4="Mapper1_unmapped", c5="Mapper2_unmapped"))
    for map1 in samfiles:
        for map2 in samfiles:
            p_dict = position_dict[map1]
            t_dict = p_dict[map2]
            print("{map1}\t{map2}\t{c1} ({p1}%)\t{c2} ({p2}%)\t{c3} ({p3}%)\t{c4} ({p4}%)\t{c5} ({p5}%)".format(map1=map1, map2=map2, 
                  c1=t_dict["same"], c2=t_dict["diff"], c3=t_dict["both_unmap"], c4=t_dict["this_unmap"], c5=t_dict["other_unmap"],
                    p1=(t_dict["same"]/count)*100, p2=(t_dict["diff"]/count)*100, p3=(t_dict["both_unmap"]/count)*100, 
                    p4=(t_dict["this_unmap"]/count)*100, p5=(t_dict["other_unmap"]/count)*100))
    
    print()
    print()
    #print mapping quality
    print("Mappinq quality statistics")
    print("Mapper1\tMapper2\t{c1}\t{c2}\t{c3}".format(c1="Same_quality", c2="Mapper1_higher_quality", c3="Mapper1_lower_quality"))
    for map1 in samfiles:
        for map2 in samfiles:
            p_dict = mapq_dict[map1]
            t_dict = p_dict[map2]
            print("{map1}\t{map2}\t{c1} ({p1}%)\t{c2} ({p2}%)\t {c3} ({p3}%)".format(map1=map1, map2=map2,
                  c1=t_dict["same"], c2=t_dict["higher"], c3=t_dict["lower"],
                    p1=(t_dict["same"]/count)*100, p2=(t_dict["higher"]/count)*100, p3=(t_dict["lower"]/count)*100))
        
    print()
    print("Mapping ranges:")
    for map1 in samfiles:
        print("{}\t0-{}".format(map1, mapq_range[map1]))
        
    print()
    print()
    print("Alignment statistics")
    print("Mapper1\tMapper2\t{c1}\t{c2}\t{c3}\t{c4}".format(c1="Same_alignment", c2="Mapper1_clipping", c3="Mapper1_alignment", c4="Not_same_location"))
    for map1 in samfiles:
        for map2 in samfiles:
            p_dict = cigar_dict[map1]
            t_dict = p_dict[map2]
            print("{map1}\t{map2}\t{c1} ({p1}%)\t{c2} ({p2}%)\t{c3} ({p3}%)\t{c4} ({p4}%)".format(map1=map1, map2=map2,
                  c1=t_dict["same"], c2=t_dict["softclip"], c3=t_dict["forcemap"], c4=t_dict["notcomp"],
                    p1=(t_dict["same"]/count)*100, p2=(t_dict["softclip"]/count)*100, p3=(t_dict["forcemap"]/count)*100, p4=(t_dict["notcomp"]/count)*100))
            
    
