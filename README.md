usage: bigwig_hmm.py [-h] [-i INPUTFILE] [-g GENOME] [-n NUM_STATES]
                     [-o OUTPUTFILE]

Create bedfile of HMM states from bigwig

optional arguments:
  -h, --help            show this help message and exit
  
  -i INPUTFILE, --inputfile INPUTFILE
                        path to the bigwig file(s)
                        
  -g GENOME, --genome GENOME
                        genome (i.e. hg38)
                        
  -n NUM_STATES, --num_states NUM_STATES
                        number of hidden states
                        
  -o OUTPUTFILE, --outputfile OUTPUTFILE
                        output filename (endswith:#_state_HMM.bed)
                        
                        
Wrapper based on https://github.com/jmschrei/pomegranate                          
