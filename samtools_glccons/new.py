#!/usr/bin/env python3

##%
with open('fa.fa', 'r') as contamix:
    with open('samtools.fa', 'r') as glcont:
        f1 = contamix.read()[5:].replace('\n', '')
        f2 = glcont.read()[5:].replace('\n', '')
        
##%
f1[1718]
f2[1718]

pos = [73, 263, 750, 1719, 6221, 12705, 14543, 14544, 14545, 14991, 14993, 14994, 14995, 14996, 14997, 16189]

##%
for p in range(len(f2)):
    if f1[p].upper() != f2[p].upper():
        print(p+1, f1[p], f2[p])