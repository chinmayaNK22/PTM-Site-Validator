import re
from itertools import islice
ptm_sites = {} #PTM list
aa_list = {} # Amino acid_list
query_list = {}
query_list_full_info = {}
dicts_final = {}
mgf_file = {}
import time

def read_mgf_file(infile):
    for i in open(infile):
        if i.startswith('TITLE'):
            scan = re.search('scan=(.*?)"', i).group(1)
            mgf_file[scan] = []
        elif re.search('^\s*[0-9]',i):
            split_items = i.split(' ')
            if float(split_items[1].rstrip())>0.0:
                mgf_file[scan].append(i.rstrip())

def read_ptm_list(infile):
    for i in open(infile):
        split_i = i.split('\t')
        ptm_sites[split_i[0]] = split_i[1] + '@' + split_i[2].rstrip()

def read_aa_list(infile):
    for i in open(infile):
        split_i = i.split('\t')
        aa_list[split_i[0]] = split_i[1] + '@' + split_i[2] + '@' + split_i[3].rstrip()

def get_header_idx(infile):
    for i in open(infile):
        split_i = i.rstrip().split('\t')
        raw_file = split_i.index("Spectrum File")
        scan = split_i.index("First Scan")
        seq = split_i.index("Annotated Sequence")
        modification = split_i.index("ptmRS: Best Site Probabilities")
        #modification = split_i.index("Modifications")
        return raw_file, scan, seq, modification
        break

infile="PXD001188_Multi-PTM_PSMs_PTM.txt" # Cahnge the input file name here
a = get_header_idx(infile)
with open(infile) as f:
    for i in islice(f, 1, None):
        split_i = i.split('\t')
        if split_i[a[0]] not in query_list:
            query_list[split_i[a[0]]] = [split_i[a[1]] + '@' + split_i[a[2]] + '@' + split_i[a[3]]] #Raw file and correspoinding scan_num, seq , modification
            query_list_full_info[split_i[a[0]]] = [i.rstrip()]
        else:
             query_list[split_i[a[0]]].append(split_i[a[1]] + '@' + split_i[a[2]] + '@' + split_i[a[3]])
             query_list_full_info[split_i[a[0]]].append([i.rstrip()])


read_ptm_list('PTM_Modification.txt')
read_aa_list('Immonium_Ions-Amino_acids.txt')
for k, v in query_list.items():
    call_mgf_file = read_mgf_file(k.split('.')[0]  + ".mgf")
    for mod_iter in range(len(v)):
        if len(v[mod_iter].split('@')) > 2 and len(v[mod_iter].split('@')[2]) > 0 :
            #print len(mod_iter.split('@')[2].split(';')), mod_iter.split('@')[2]
            get_mz = mgf_file[v[mod_iter].split('@')[0]] #dictionary of scan number
            if len(v[mod_iter].split('@')[2].split(';')) <= 1:
                pass
                scan = v[mod_iter].split('@')[0]
                get_ptm_aa = re.findall(r'\d+', v[mod_iter].split('@')[2].split(':')[0])[0] #Fetch the ptm aminoacid
                mod_aa = v[mod_iter].split('@')[2].split(';')[0][0]#v[mod_iter].split('@')[1].split('.')[1][int(get_ptm_aa[0])  - 1].upper() #Get the modified amino acid
                #print (v[mod_iter].split('@')[2].split('(')[1].replace(')', '').split(':')[0])
                get_da = aa_list[mod_aa].split('@')[0]
                get_mod_da = ptm_sites[v[mod_iter].split('@')[2].split('(')[1].replace(')', '').split(':')[0]].split('@')[0]
                #get_mod_da = ptm_sites[mod_iter.split('@')[2].split('(')[1].replace(')', '')].split('@')[0]
                immonium_ion = float(get_da) +  float(get_mod_da)
                for iter_items in range(len(get_mz)):
                    #print get_mz[iter_items].split(' ')[0], get_mz[iter_items].split(' ')[1]
                    if abs(float(get_mz[iter_items].split(' ')[0]) - float(immonium_ion)) <= 0.02:
                        if k.split('.')[0]  + ".mgf" + "_" + scan + "_" + v[mod_iter].split('@')[2] not in dicts_final:
                            dicts_final[k.split('.')[0]  + ".mgf" + "_" + scan + "_" + v[mod_iter].split('@')[2]] = ''.join(query_list_full_info[k][mod_iter]) + "\t" + str(get_mz[iter_items].split(' ')[0]) + '\t' + str(get_mz[iter_items].split(' ')[1])
                            #print (''.join(query_list_full_info[k][mod_iter]) + '\t' + k + '\t' + str(scan) + '\t' + str(get_mz[iter_items].split(' ')[0]) + '\t' + str(get_mz[iter_items].split(' ')[1]) + '\t' + " Diff: " + '\t' + str(float(get_mz[iter_items].split(' ')[0]) - float(immonium_ion)))
            else:
                #immonium_ion = 0.0
                for iters in range(len(v[mod_iter].split('@')[2].split(';'))):
                    #print (mod_iter.split('@')[2].split(';')[iters], mod_iter.split('@')[2])
                    scan = v[mod_iter].split('@')[0]
                    get_ptm_aa = re.findall(r'\d+', v[mod_iter].split('@')[2].split(';')[iters].split(':')[0])[0] #Fetch the ptm aminoacid
                    mod_aa = v[mod_iter].split('@')[2].split(';')[iters].strip()[0] #Get the modified amino acid
                    get_da = aa_list[mod_aa].split('@')[0]
                    get_mod_da = ptm_sites[re.search('\((.*?)\)', v[mod_iter].split('@')[2].split(';')[iters]).group(1)].split('@')[0]
                    #get_mod_da = ptm_sites[mod_iter.split('@')[2].split('(')[1].replace(')', '')].split('@')[0]
                    #print (scan, get_ptm_aa, mod_aa, get_da, get_mod_da, mod_iter)
                    immonium_ion = float(get_da) +  float(get_mod_da)
                    
                    for iter_items in range(len(get_mz)):
                        if abs(float(get_mz[iter_items].split(' ')[0]) - float(immonium_ion)) <= 0.02:
                            #file + scan + modification(R52(Demadation) not in dictionary
                            if k.split('.')[0]  + ".mgf" + "_" + scan + "_" + v[mod_iter].split('@')[2].split(';')[iters] not in dicts_final:
                                dicts_final[k.split('.')[0]  + ".mgf" + "_" + scan + "_" + v[mod_iter].split('@')[2].split(';')[iters]] = ''.join(query_list_full_info[k][mod_iter]) + '\t' + str(get_mz[iter_items].split(' ')[0]) + '\t' + str(get_mz[iter_items].split(' ')[1])
                            #   print (''.join(query_list_full_info[k][mod_iter]) + '\t' + k + '\t' + str(scan) + '\t'  + str(get_mz[iter_items].split(' ')[0]) + '\t' + str(get_mz[iter_items].split(' ')[1]) + '\t' +  " Diff: " + '\t' + str(float(get_mz[iter_items].split(' ')[0]) - float(immonium_ion)))
                

# Reading the input file to map the results
write_file = open('PXD001188_PTM_Immonium_Ion_Validated_130620.txt', 'w')
read_header = get_header_idx(infile)
with open(infile) as f:
    for i in islice(f, 1, None):
        split_i = i.split('\t')
        if len(split_i[read_header[3]]) > 0:
            if ";" not in split_i[read_header[3]]:
                if split_i[read_header[0]].split('.')[0] + ".mgf" + "_" + split_i[read_header[1]] + "_" + split_i[read_header[3]] in dicts_final:
                    write_file.write(dicts_final[split_i[read_header[0]].split('.')[0] + ".mgf" + "_" + split_i[read_header[1]] + "_" + split_i[read_header[3]]] + '\n')
                else:
                    write_file.write(i.rstrip() + "\t" + "-" + "\t" + "-" + '\n')
            else:
                items = ""
                for ik in range(len(split_i[read_header[3]].split(';'))):
                    #print (split_i[read_header[0]].split('.')[0] + ".mgf" + "_" + split_i[read_header[1]] + "_" + split_i[read_header[3]].split(';')[ik])
                    if split_i[read_header[0]].split('.')[0] + ".mgf" + "_" + split_i[read_header[1]] + "_" + split_i[read_header[3]].split(';')[ik] in dicts_final:
                        write_file.write(dicts_final[split_i[read_header[0]].split('.')[0] + ".mgf" + "_" + split_i[read_header[1]] + "_" + split_i[read_header[3]].split(';')[ik]] + '\n')
                    else:
                        write_file.write(i.rstrip() + '\t' + "-" + "\t" + "-" + '\n')
        else:
            write_file.write(i.rstrip() + '\t' + "-" + "\t" + "-" + '\n')
write_file.close()
