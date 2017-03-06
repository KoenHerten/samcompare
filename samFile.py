#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 11:47:58 2017

@author: Koen Herten

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

import samRead

class samFile:
    
    def __init__(self, file_name):
        '''
        initiate all variables
        '''
        self._file = open(file_name)
        self._samread = None
        
    def _nextSamRead(self):
        line = self._file.readline()
        while (line.startswith("@")):
            #header
            line = self._file.readline()
        if (line is None or line is ""):
            return None
        line = line.rstrip()
        samread = samRead.samRead(line)
        return samread
        
    def nextPrimaryRead(self):
        samread = self._nextSamRead()
        if (self._samread is None):
            self._samread = samread
            return samread
        while (samread is not None and (samread.isSecondaryAlignment() or (samread.qname == self._samread.qname and not samread.issecond()))):
        #while (samread is not None and (samread.isSecondaryAlignment())):
            samread = self._nextSamRead()
        self._samread = samread
        return samread
        
    def getCurrentSamRead(self):
        return self._samread
        
    def close(self):
        self._file.close()
        