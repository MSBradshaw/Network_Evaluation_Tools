import pandas as pd

naming = pd.read_csv('Data/9606.protein.info.v11.0.txt',sep='\t')

name_map = {naming['protein_external_id'][i]:naming['preferred_name'][i] for i in range(naming.shape[0])}

skip_count = 0
good_count = 0
with open('Data/string_edge_list.tsv','r') as input:
    with open('Data/string_edge_list_common_names.tsv','w') as output:
        for line in input:
            line = line.strip()
            # it space delimited don't know why...
            edge = line.split('   ')
            if edge[0] in name_map and edge[1] in name_map:
                output.write(name_map[edge[0]] + '\t' + name_map[edge[1]] + '\n')
                good_count += 1
            else:
                print(edge[0] + ' ' + edge[1])
                skip_count += 1
print(skip_count)
print(good_count)
