if __name__ == "__main__":
    
    import argparse
    import pandas as pd
    import os
    from subprocess import run

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-input', '--input_path', help='in path', required=True)
    parser.add_argument('-output', '--output_path', help='out path', required=True)
    parser.add_argument('-coef', help='coef val', type=float, required=True)


    args = vars(parser.parse_args())

    i_path = args['input_path']
    o_path = args['output_path']
    coef = args['coef']

    def get_taxlist(subset):
    
        taxlist = []
        
        for i in list(subset['staxids'].values):
            
            taxlist.extend(map(float, str(i).split(';')))
            
        taxlist = list(map(int, taxlist))
        ancestors = str(taxlist)[1:-1].replace(', ', ';')
        
        return ancestors
    
    diam_res = pd.read_csv(i_path, sep='\t').dropna(subset=['staxids'])
    
    LCA_res = {}
    
    for cnt in diam_res['qseqid'].unique():
        
        cnt_subset = diam_res[diam_res['qseqid'] == cnt]
        cnt_subset = cnt_subset[cnt_subset['bitscore'] >= cnt_subset['bitscore'].max() * coef]
        ancestors = get_taxlist(cnt_subset)
        result = run(f"echo '{ancestors}'| taxonkit lca -s ';'| cut -f2 |taxonkit lineage -R", 
                    shell=True, capture_output=True, text=True)   
        lineage_res = result.stdout.replace('\n', '').split('\t')
        LCA_res[cnt] = dict(zip(lineage_res[2].split(';'), lineage_res[1].split(';')))
        
        top_hit_id = cnt_subset['staxids'].values[0]
        tophit_res = run(f"echo '{top_hit_id}'| taxonkit lca -s ';'| cut -f2 |taxonkit lineage -R", 
                    shell=True, capture_output=True, text=True).stdout.replace('\n', '').split('\t')[1]
        LCA_res[cnt]['Top hit'] = tophit_res
    
    LCA_res = pd.DataFrame(LCA_res).T
    
    LCA_res = pd.DataFrame(LCA_res).T
    LCA_res.to_csv(o_path, sep='\t')
