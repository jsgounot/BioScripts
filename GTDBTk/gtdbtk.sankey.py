import click
import pandas as pd

@click.command()
@click.argument('table', type=str)
def run(table) :

	print ("Produce input for sankeymatic: http://sankeymatic.com/build/")
	print ("Note that you might want to dereplicate results before (TODO)")

	df = pd.read_csv(table, sep='\t')

	colnames = ["division", "phylum", "class", "order"]
	for idx in range(len(colnames) - 1) :
	    pair = colnames[idx:idx+2]
	    sdf = df.groupby(pair).size().rename("count").reset_index()
	    left, right = pair
	    for _, row in sdf.iterrows() :
	        print ("%s [%i] %s" %(row[left], row["count"], row[right]))

if __name__ == '__main__':
    run()