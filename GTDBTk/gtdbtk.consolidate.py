import click
import pandas as pd

from collections import defaultdict

@click.command()
@click.argument('bac', type=str)
@click.argument('arc', type=str)
@click.argument('field', type=str)
@click.argument('outfile', type=str)
def run(bac, arc, field, outfile) :

	bac = pd.read_csv(bac, sep="\t")
	arc = pd.read_csv(arc, sep="\t")
	df = pd.concat([bac, arc])

	df = make_tax_dataframe(df, field)
	df.to_csv(outfile, sep="\t", index=False)
	
def make_tax_dataframe(df, column) :

    data = df.set_index("user_genome")[column].to_dict()

    def extract(taxonomy) :
        values = {}
        taxonomy = str(taxonomy)
        for element in taxonomy.split(";") :
            if not "__" in element : continue
            k, v = element.split("__")
            if v : values[k] = v
        return values

    data = {k : extract(taxonomy) for k, taxonomy in data.items()}

    reverse = defaultdict(lambda : defaultdict(list))
    for mag, taxonomy in data.items() :
        for rank, value in taxonomy.items() :
            reverse[rank][value].append(mag)

    cnames = {"d" : "division", 
              "p" : "phylum",
              "c" : "class",
              "o" : "order",
              "f" : "family",
              "g" : "genus",
              "s"  : "species"
            }

    ssdf = pd.DataFrame([{"name" : name, ** values} for name, values in data.items()])
    ssdf.columns = [cnames.get(column, column) for column in ssdf.columns]

    print ("Count: # MAGs with rank assignment (not NaN)")
    print (ssdf[ssdf.columns[1:]].agg(['nunique','count','size']))
       
    return ssdf

if __name__ == '__main__':
    run()