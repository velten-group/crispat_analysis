import os
import sys, getopt
import json
import time
import pandas as pd
import crispat

def main():
    t_start = time.time()
    
    # Parse command line arguments
    argv = sys.argv[1:]
    try:
        options, args = getopt.getopt(argv, "c:", ["config="])
    except:
        print('Incorrect arguments!')
        
    for name, value in options:
        if name in ('-c', '--config'):
            config_filename = value
    
    # Parse config file
    config = json.load(open(config_filename))
    
    # run crispat
    crispat.ga_binomial(config['input_file'], 
                        config['out_dir'] + "binomial/",
                        start_gRNA = 0, 
                        parallelize = True, 
                        n_processes = 32)
        
    # Save run time
    t_end = time.time()
    print(t_end - t_start)
    run_time = pd.DataFrame({'method': ['binomial'], 'time': [t_end - t_start]})
    run_time.to_csv(config['out_dir'] + 'binomial/run_time.csv', index = False)
    
    print('Done')

if __name__ == "__main__":
    main()