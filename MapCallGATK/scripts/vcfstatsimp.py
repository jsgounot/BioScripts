"""
Improve VCF stats results by adding information
from the reference
"""

vcfstat   = snakemake.input[0]
reference = snakemake.input[1]
outfile   = snakemake.output[0]

elements =  ["SNPs", "MNPs", "Insertions", "Deletions", "Indels"]
elements += ["Haploid " + element for element in elements]


def reference_size(fname) :
    size = 0
    with open(fname) as f :
        for line in f :
            if line.startswith(">") :
                continue
            size += len(line.strip())
    return size

def load_vcfstat(fname) : 
    with open(vcfstat) as f :
        data = {}
        for line in f :
            key, value = line.strip().split(":")
            key, value = key.strip(), value.strip()
            try : value = int(value)
            except ValueError : value = value
            data[key] = value

    return data

def make_new_stat(data, refsize) :
    return {element : data.get(element, 0) * 1000 / refsize
        for element in elements if element in elements}

def write_data(subdata, outfile) :
    if subdata :
        maxlen = max(len(element) for element in subdata)
    
    with open(outfile, "w") as f :
        for element in elements :
            welement = element + " (variants/kb)" + " " *( maxlen - len(element))
            f.write("%s : %f\n" %(welement, subdata[element]))

refsize = reference_size(reference)
vcfstat = load_vcfstat(vcfstat)
subdata = make_new_stat(vcfstat, refsize)
write_data(subdata, outfile)


