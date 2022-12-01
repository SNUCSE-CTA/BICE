from itertools import product
import os

def generate_querys(path_data_graph):
    sizes = [ 50, 100, 150, 200 ]
    
    os.system('mkdir -p {}'.format('sparse'))
    os.system('mkdir -p {}'.format('nonsparse'))

    for size in sizes:
        print(size)
        #os.system('mkdir -p {}/{}'.format('sparse', size))
        #os.system('mkdir -p {}/{}'.format('nonsparse', size))
        cmd = './makequery_default {} 100 {} {} {}'.format(
            path_data_graph, 
            str(size), 
            'sparse/{}.igraph'.format(size), 
            'nonsparse/{}.igraph'.format(size)
        )
        os.system(cmd)

if __name__=='__main__':
    os.system('sh compile.sh')
    path_data_graph = 'yeast.igraph'
    generate_querys(path_data_graph)