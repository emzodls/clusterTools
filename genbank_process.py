import os
from glob import glob

genbank_files  = glob('/Users/emzodls/Documents/unmapped_burk/*.final.gbk')

for gbk in genbank_files:
    path,name = os.path.split(gbk)
    name_base = name.split('_unmapped')[0]
    ctg_ctr = 0
    for line in open(gbk):
        if 'LOCUS' in line:
            line_parse  = line.split('..')
            scafNameChange = line_parse[0][12:]
            ctg_ctr += 1
            lineToWrite = line.replace('%s..' % scafNameChange,'%s_ctg%.3i' %(name_base,ctg_ctr))
            #print lineToWrite
        elif scafNameChange in line:
            lineToWrite = line.replace('%s' % scafNameChange, '%s_ctg%.3i' % (name_base, ctg_ctr))
            #print lineToWrite
        else:
            lineToWrite = line
        with open(os.path.join(path,'%s_edited.gbk' % name_base),'ab') as outfile:
            outfile.write('%s' % lineToWrite)


