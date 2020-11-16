def parse_gene_set_gmt(filepath):
    gs_lookup = dict()
    with open(filepath) as f:
        for line in f:
            line = line.rstrip().split("\t")
            name = line[0]
            genes = set(line[2:])
            gs_lookup[name] = genes
    return gs_lookup
