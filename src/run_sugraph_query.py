import os
import queue
from itertools import product
import json

def get_job_grapes(binary, data_path, query_path, output_path, log_path):
    re = "timeout 610s ./" + binary + " " +\
    "1 -f -gfu" + " " + data_path + " " +\
    query_path + " " + "-file" + " " +\
    ">>" + " " + log_path
    #"-o" + " " + output_path 
    return re

def go_job_grapes(binary, data_path, query_path, output_path, log_path):
    re = "timeout 610s ./" + binary + " " +\
    "1 -f -gfu" + " " + data_path + " " +\
    query_path + " " + "-file" + " " +\
    ">>" + " " + log_path
    return os.system(re)

def get_job_daf(binary, data_path, query_path, output_path, log_path):
    re = "timeout 610s ./" + binary + " " +\
    "-dg" + " " + data_path + " " +\
    "-qg" + " " + query_path + " " +\
    ">>" + " " + log_path
    return re

def go_job_daf(binary, data_path, query_path, log_path = "log.txt"):
    re = "timeout 610s ./" + binary + " " +\
    "-dg" + " " + data_path + " " +\
    "-qg" + " " + query_path + " " +\
    ">>" + " " + log_path
    return os.system(re)

def get_job_cfql(binary, data_path, query_path, output_path, log_path):
    re = "timeout 610s ./" + binary + " " +\
    "-d" + " " + data_path + " " + "-df FA" + " " +\
    "-q" + " " + query_path + " " + " -qf FA" + " " +\
    "-o" + " " + output_path + " " + ">>" + " " + log_path
    return re

def go_job_cfql(binary, data_path, query_path, output_path, log_path):
    re = "timeout 610s ./" + binary + " " +\
    "-d" + " " + data_path + " " + "-df FA" + " " +\
    "-q" + " " + query_path + " " + " -qf FA" + " " +\
    "-o" + " " + output_path + " " + ">>" + " " + log_path
    return os.system(re)

def get_data_path(name):
    return "../data/" + name + ".gfu"

def get_query_path(dir, name, density, sz, idx):
    return "../" + dir + "/" + name + "/" + density + "/" + str(sz) + "/q" + str(idx) + ".gfu"

def get_output_path(binary, name, density, sz, idx):
    dir_path = "../out/" + binary + "/" + name + "/" + density + "/" + str(sz) 
    os.system("mkdir -p" + " " + dir_path)
    return dir_path + "/" + str(idx) + ".txt"

def exp(
    binary, name, density, size, query_indces, log_json_path, 
    dir = 'finalQuery20200428'
):
    log_dict_queries = {}

    num_solved_quries = 0
    average_processing_time = 0.0

    for query_index in query_indces:
        #temporary buffer for each query file
        buffer = "buffer.txt"
        #buffer = 'buffer_' + log_json_path.split('/')[-1]
        os.system("touch" + " " + buffer)

        #generate
        log_dict_query = {}

        #execute
        qp = get_query_path(dir, name, density, size, query_index)
        dp = get_data_path(name)
        
        exit_code = go_job_daf(binary = binary, data_path = dp, query_path = qp, log_path = buffer) 

        lines = []
        with open(buffer, 'r') as f:
            while True:
                line = f.readline()
                if line:
                    lines.append(line)
                else:
                    break
    
        if exit_code == 0:
            log_dict_query['complete'] = True  
            log_dict_query['rec_num'] = int(lines[-6].split()[-1])
            log_dict_query['cad_num'] = int(lines[-5].split()[-1])
            log_dict_query['ans_num'] = int(lines[-4].split()[-1])
            log_dict_query['fil_time'] = float(lines[-3].split()[-1])
            log_dict_query['ver_time'] = float(lines[-2].split()[-1])
            log_dict_query['pro_time'] = float(lines[-1].split()[-1])

            num_solved_quries = num_solved_quries + 1
            average_processing_time = average_processing_time + log_dict_query['pro_time']
        else:
            log_dict_query['complete'] = False
            log_dict_query['last_id'] = int(lines[-1].split()[-1])
            average_processing_time = average_processing_time + 6e5
        
        #store
        log_dict_queries[query_index] = log_dict_query

        #delete buffer
        os.system('rm' + ' ' + buffer)
    
    log_dict_queries['num_solved'] = num_solved_quries
    if num_solved_quries == 0:
        log_dict_queries['avg_pro_time'] = 0
    else:
        log_dict_queries['avg_pro_time'] = average_processing_time / len(query_indces)
    with open(log_json_path, "w") as jf:
        json.dump(log_dict_queries, jf, indent = 4)
    return log_dict_queries

def exp_with_log(
    binary, name, density, size, query_indces, log_json_path, buffer_dir,
    dir = 'finalQuery20200428'
):
    log_dict_queries = {}

    num_solved_quries = 0
    average_processing_time = 0.0

    for query_index in query_indces:
        #temporary buffer for each query file
        buffer = buffer_dir + "/" + "_".join(
            [ binary, name, density, str(size) ]
        ) + ".txt"
        os.system("touch" + " " + buffer)

        #generate
        log_dict_query = {}

        #execute
        qp = get_query_path(dir, name, density, size, query_index)
        dp = get_data_path(name)
        
        exit_code = go_job_daf(binary = binary, data_path = dp, query_path = qp, log_path = buffer) 

        lines = []
        with open(buffer, 'r') as f:
            while True:
                line = f.readline()
                if line:
                    lines.append(line)
                else:
                    break
    
        if exit_code == 0:
            log_dict_query['complete'] = True  
            log_dict_query['rec_num'] = int(lines[-6].split()[-1])
            log_dict_query['cad_num'] = int(lines[-5].split()[-1])
            log_dict_query['ans_num'] = int(lines[-4].split()[-1])
            log_dict_query['fil_time'] = float(lines[-3].split()[-1])
            log_dict_query['ver_time'] = float(lines[-2].split()[-1])
            log_dict_query['pro_time'] = float(lines[-1].split()[-1])

            num_solved_quries = num_solved_quries + 1
            average_processing_time = average_processing_time + log_dict_query['pro_time']
        else:
            log_dict_query['complete'] = False
            log_dict_query['last_id'] = int(lines[-1].split()[-1])

        #store
        log_dict_queries[query_index] = log_dict_query
    
    log_dict_queries['num_solved'] = num_solved_quries
    if num_solved_quries == 0:
        log_dict_queries['avg_pro_time'] = 0
    else:
        log_dict_queries['avg_pro_time'] = average_processing_time / num_solved_quries

    with open(log_json_path, "w") as jf:
        json.dump(log_dict_queries, jf, indent = 4)
    return log_dict_queries

def experiment(
    ours = 'ours',
    debugging = 'debugging', buffer = 'buffer',
    num_queries = 20, names = [ "pcms", "REDDIT-MULTI-5K", "IMDB-MULTI" ],
    sizes = [ 8, 16 ], densities = [ 'randomwalk', 'bfs' ]
):
    query_indces = [ i for i in range(num_queries) ]
  
    os.system('rm -r {}'.format(debugging))
    os.system('mkdir -p {}'.format(debugging))

    os.system('rm -r {}'.format(buffer))
    os.system('mkdir -p {}'.format(buffer))

    os.system('make clean')
    os.system('make {}'.format(ours))
    
    num_error_quries = 0

    for name, density, size in product(names, densities, sizes):
        print(name, density, size)
        log_ours_path = debugging + "/" + "_".join([ ours, name, density, str(size) ]) + ".json"
        
        results = exp_with_log(binary = ours, name = name, 
        density = density, size = size,
        query_indces = query_indces, buffer_dir = buffer, log_json_path = log_ours_path)        
             
def comparison(
    original = 'simp', ours = 'ours',
    debugging = 'debugging', buffer = 'buffer',
    query_indces = [ 0 ], names = [ "pcms", "REDDIT-MULTI-5K", "IMDB-MULTI" ],
    sizes = [ 8, 16 ], densities = [ 'randomwalk', 'bfs' ], erase = True
):
    if erase == True:
        os.system('rm -r {}'.format(debugging))
        os.system('mkdir -p {}'.format(debugging))

        os.system('rm -r {}'.format(buffer))
        os.system('mkdir -p {}'.format(buffer))

    os.system('make clean')
    os.system('make {}'.format(ours))
    os.system('make {}'.format(original))
    
    num_error_quries = 0

    for name, density, size in product(names, densities, sizes):
        print(name, density, size)
        log_original_path = debugging + "/" + "_".join([ original, name, density, str(size) ]) + ".json" 
        log_ours_path = debugging + "/" + "_".join([ ours, name, density, str(size) ]) + ".json"

        results1 = exp_with_log(binary = original, name = name, 
        density = density, size = size,
		query_indces = query_indces, buffer_dir = buffer, log_json_path = log_original_path)
        
        results2 = exp_with_log(binary = ours, name = name, 
        density = density, size = size,
        query_indces = query_indces, buffer_dir = buffer, log_json_path = log_ours_path)        

        for k in query_indces:
            if results1[k]['complete'] == True and results2[k]['complete'] == True:
                try:
                    assert results1[k]['ans_num'] == results2[k]['ans_num']
                    #assert results1[k]['rec_num'] == results2[k]['rec_num']
                except BaseException:
                    num_error_quries = num_error_quries + 1
                    print("Solutions of {}th query aren't coincide".format(k))
    
    print('num_error_quries', num_error_quries)

def light_test():
    name_exp = 'light_test'
    comparison(original = 'veq', ours = 'adv_van',
    debugging = 'log/{}/debugging'.format(name_exp), 
    buffer = 'log/{}/buffer'.format(name_exp),
    query_indces = [ i for i in range(100) ],
    names = [ "pcms", "REDDIT-MULTI-5K", "IMDB-MULTI", "pdbs" ],
    sizes = [ 8, 16, 32, 64 ], densities = [ 'randomwalk', 'bfs' ])

def specific_test():
    name_exp = 'specific_test'
    comparison(original = 'veq', ours = 'adv_van',
    debugging = 'log/{}/debugging'.format(name_exp), 
    buffer = 'log/{}/buffer'.format(name_exp),
    query_indces = [ 16 ],
    names = [ "IMDB-MULTI" ],
    sizes = [ 32 ], densities = [ 'randomwalk' ], erase = True)

def ppigo_comparison():
    name_exp = 'ppigo_comparison'
    comparison(original = 'eq', ours = 'mm',
    debugging = 'log/{}/debugging'.format(name_exp), 
    buffer = 'log/{}/buffer'.format(name_exp),
    query_indces = [ i for i in range(100) ],
    names = [ "ppigo" ],
    sizes = [ 32, 64 ], densities = [ 'randomwalk' ])

def COLLAB_comparison():
    name_exp = 'COLLAB_comparison'
    comparison(original = 'eq', ours = 'mm',
    debugging = 'log/{}/debugging'.format(name_exp), 
    buffer = 'log/{}/buffer'.format(name_exp),
    query_indces = [ i for i in range(100) ],
    names = [ "COLLAB" ],
    sizes = [ 32 ], densities = [ 'randomwalk' ])

if __name__=="__main__":
    #vertification()
    #experiment()
    
    light_test()
    #specific_test()
    #ppigo_comparison()
    #COLLAB_comparison()