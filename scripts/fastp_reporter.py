import json
#	https://stackoverflow.com/questions/19483351/converting-json-string-to-dictionary-not-list
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("json_in", help="JSON file to be condensed")
parser.add_argument("flat_out", help="flatfile summary")
parser.add_argument("-t", "--tag", help="line-name for the processed JSON", default=None)
args = parser.parse_args()

def prune_jason(jsn):
	sparse_dict = {}
	sparse_dict = dict([(k,{}) for k in ['prefiltered', 'postfiltered', 'filtration', 'duplication']])

	for i in ['gc_content','q20_rate', 'q30_rate', 'read1_mean_length','total_reads']:
		try:
			sparse_dict['prefiltered'][i] =  jsn['summary']['before_filtering'][i] 
		except KeyError:
			sparse_dict['prefiltered'][i] = "NA"

		try:
			sparse_dict['postfiltered'][i] =  jsn['summary']['after_filtering'][i] 
		except KeyError:
			sparse_dict['postfiltered'][i] = "NA"

	for i in ['adaptor_trimmed_reads']:
		try: 
			sparse_dict['filtration'][i] = jsn['adapter_cutting']['adapter_trimmed_reads']
		except KeyError:
			sparse_dict['filtration'][i] = "NA"

	for i in ['corrected_reads','low_quality_reads','passed_filter_reads']:
		try:
			sparse_dict['filtration'][i] = jsn['filtering_result'][i]
		except KeyError:
			sparse_dict['filtration'][i] = "NA"

	for i in ['rate']:
		sparse_dict['duplication'][i] = jsn['duplication'][i]

	return sparse_dict

def print_pruned(prn):
	lines2write = [ [x, y, prn[x][y]]  for x in prn.keys() for y in prn[x]]
	if args.tag:
		[ ell.insert(0, args.tag) for ell in lines2write ]
	phial_out = open(args.flat_out, 'w')
	for preline in lines2write:
		field_count = len(preline)
		line = ("%s" + "\t%s"*(field_count-1) + "\n") % tuple(preline)
		phial_out.write(line)
	phial_out.close()

jason_file=open(args.json_in)
jason_str = jason_file.read()
jason = json.loads(jason_str)

jason_condensed = prune_jason(jason)
print_pruned(jason_condensed)
