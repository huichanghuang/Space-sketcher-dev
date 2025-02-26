import os
import sys
import json
import time
import logging
import sys
import io
import shutil
import base64
import subprocess
from datetime import datetime
# from dnbc4tools.__init__ import __root_dir__
from __init__ import __root_dir__

# Extracted common path construction logic into a separate function
def get_common_path_part():  
    return '/'.join(str(__root_dir__).split('/')[0:-4])

# Use os.makedirs with exist_ok=True instead of os.system
# Construct target path using os.path.join for better portability
# Handle exceptions and provide meaningful error messages
def str_mkdir(arg):
    """
    Exceptions:
    - If the directory already exists, no action is taken.
    - If there is no permission to create the directory, raises a PermissionError.
    - If any other OS error occurs, it is raised as-is.
    """
    try:
        os.makedirs(arg)  # Try creating the directory
    except FileExistsError:
        pass  # If the directory already exists, take no action
    except PermissionError:
        raise PermissionError("Permission denied to create the directory")  # If no permission to create, raise an exception
    except Exception as e:
        raise (f"Failed to create directory {arg}: {e}")  

# Use string formatting instead of concatenation for readability
# Handle exceptions and provide meaningful error messages    
def change_path():
    common_path = get_common_path_part()
    try:
        os.environ['PATH'] += f":{common_path}/bin" 
        os.environ['LD_LIBRARY_PATH'] = f"{common_path}/lib"
    except Exception as e:
        raise (f"Failed to update environment variables: {e}")  

# Construct bin path using os.path.join for better portability
def bin_path():
    return os.path.join(get_common_path_part(), 'bin')  
    
def rm_temp(*args):
    """
    Remove specified files or directories, including their contents if they are directories.
    
    Parameters:
    *args (str): A variable-length argument list containing the file/directory paths to be removed.

    Note:
    - Symbolic links pointing to directories are skipped without attempting removal.
    - Any encountered exceptions during the removal process are caught and printed to stdout.
    """
    for filename in args:
        try:
            if os.path.exists(filename):
                if os.path.isdir(filename):
                    if os.path.islink(filename) or os.path.realpath(filename) != filename:
                        print(f"Skipped symbolic link: {filename}")
                        continue
                    shutil.rmtree(filename)
                else:
                    os.remove(filename)
            else:
                pass
        except Exception as e:
            print(f"Error removing {filename}: {e}")

def get_formatted_time():
    current_time = datetime.now()
    formatted_time = current_time.strftime('%Y-%m-%d %H:%M:%S')
    return formatted_time

def construct_logfile_path(log_dir):
    """Constructs the logfile path."""
    today = time.strftime('%Y%m%d', time.localtime(time.time()))
    return f'{log_dir}/log/{today}.txt'

def validate_log_dir(log_dir):
    """Validates the log directory."""
    if not os.path.exists(log_dir):
        os.makedirs(log_dir, exist_ok=True)
    if not os.access(log_dir, os.W_OK):
        raise ValueError(f"Log directory '{log_dir}' is not writable.")

def setup_logging(name, log_dir):
    """Sets up logging configuration."""
    validate_log_dir(log_dir)
    logfile = construct_logfile_path(log_dir)
    
    logger = logging.getLogger(name)
    
    if not logger.handlers:
        logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        
        file_handler = logging.FileHandler(logfile, encoding="utf8")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.ERROR)
        console_handler.setFormatter(formatter)
        
        logger.addHandler(file_handler)
        logger.addHandler(console_handler)
    
    return logger

def logging_call(popenargs, name, log_dir):
    logger = setup_logging(name, log_dir)
    logger.info('Executing command: %s', ''.join(popenargs))
    try:
        output = subprocess.check_output(popenargs, shell=True, stderr=subprocess.STDOUT, universal_newlines=True)
        logger.info('%s', output)
    except subprocess.CalledProcessError as e:
        logger.error('Command failed with exit code %d', e.returncode)
        logger.error('%s', e.output)
        raise e
    except Exception as e:
        raise e

def start_print_cmd(arg, log_dir):
    validate_log_dir(log_dir)
    logfile = construct_logfile_path(log_dir)
    logging.basicConfig(filename=logfile,level=logging.INFO, format='%(message)s')
    logger = logging.getLogger()
    logger.info(arg)
    subprocess.check_call(arg, shell=True)


class StdoutAdapter:
    def __init__(self, handler):
        self.handler = handler

    def write(self, message):
        record = logging.LogRecord(
            name='stdout',
            level=logging.INFO,
            pathname=None,
            lineno=None,
            msg=message.rstrip('\n'),
            args=None,
            exc_info=None
        )
        self.handler.emit(record)

    def flush(self):
        self.handler.flush()

def logfunc(func):
    """A decorator to log function execution and output."""
    from functools import wraps
    default_outdir = '.'

    @wraps(func)
    def wrapper(*args, **kwargs):
        logger = logging.getLogger(func.__name__)
        outdir = kwargs.pop('logdir', default_outdir)
        validate_log_dir(outdir)
        logfile = construct_logfile_path(outdir)
        date_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        stream_handler = logging.StreamHandler(open(logfile, 'a', encoding='utf-8'))
        stream_handler.setLevel(logging.INFO)
        if not logger.handlers:
            logger.addHandler(stream_handler)

        adapter = StdoutAdapter(stream_handler)
        original_stdout = sys.stdout
        sys.stdout = adapter

        try:
            stream_handler.setFormatter(date_formatter)
            result = func(*args, **kwargs)
        finally:
            sys.stdout = original_stdout

        return result
    return wrapper


def judgeFilexits(*args):
    # Flatten args and filter out empty strings
    files_to_check = [file for arg in args for file in arg.split(',') if file]
    # Check for empty file names
    if any(not file for file in files_to_check):
        print("Error: Received empty file name(s).")
        return
    # Use a set comprehension to check file existence and collect missing files
    missing_files = {file for file in files_to_check if not os.path.exists(file)}
    # If any files are missing, print error messages and raise a custom exception
    if missing_files:
        error_msgs = [" ------------------------------------------------"]
        for file in missing_files:
            error_msgs.append("Error: Cannot find input file or directory '{}'".format(file))
        error_msgs.append(" ------------------------------------------------")
        print("\n".join(error_msgs), end="\n\n")
        raise FileNotFoundError("One or more input files do not exist.")


def hamming_distance(chain1, chain2):
    """
    Compute the Hamming distance between two DNA sequences.
    Args:
        chain1 (str): The first DNA sequence.
        chain2 (str): The second DNA sequence.
    Raises:
        TypeError: If either input is not a string.
        ValueError: If the lengths of the input sequences differ.
    Returns:
        int: The Hamming distance between the two sequences.
    """
    if not (isinstance(chain1, str) and isinstance(chain2, str)):
        raise TypeError("Both inputs must be strings")
    if len(chain1) != len(chain2):
        raise ValueError("Both strings must have the same length")
    
    return len(list(filter(lambda x: ord(x[0]) ^ ord(x[1]), zip(chain1, chain2))))

def read_json(file):
    """
    Read and parse a JSON file.
    Args:
        file (str): The path to the JSON file.
    Returns:
        dict or list: The parsed JSON data. If an error occurs during parsing,
                      None is returned, and an error message is printed.
    """
    try:
        with open(file, 'r', encoding='utf8') as fp:
            json_data = json.load(fp)
        return json_data
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON from file {file}: {e}")
        return None

def seq_comp(seq):
    """
    Compute the numerical representation of a DNA sequence.
    Args:
        seq (str): The DNA sequence.
    Raises:
        ValueError: If the input is not a non-empty string or contains invalid nucleotides.
    Returns:
        str: The numerical representation of the input sequence.
    """
    NT_COMP = {'A': '0', 'C': '1', 'G': '2', 'T': '3'}
    if not isinstance(seq, str) or len(seq) == 0:
        raise ValueError("Input sequence must be a non-empty string")
    if not all(n in NT_COMP for n in seq.upper()):
        raise ValueError("Input sequence contains invalid nucleotides")
    
    length = len(seq) - 1
    sum = 0
    for k, v in enumerate(seq.upper()):
        sum += int(NT_COMP[v]) * (4 ** (length - k))
    return str('%010x' % sum).upper()

def png_to_base64(file, base64_path):
    """
    Convert a PNG image file to a Base64-encoded string and write it to an HTML file.
    Args:
        file (str): The path to the PNG image file.
        base64_path (str): The path to the output HTML file.
    Returns:
        None
    Raises:
        FileNotFoundError: If the input PNG file does not exist.
    """
    if not os.path.isfile(file):
        print(f"File {file} does not exist")
        return
    with open(file, "rb") as f:
        base64_data = base64.b64encode(f.read())
        s = base64_data.decode()
        with open(base64_path, 'w') as base64_path_f:
            base64_path_f.write(f'<img src=data:image/png;base64,{s}>')


def csv_datatable(file,outfile):
    import pandas as pd
    if not os.path.exists(file):
        print(f"File {file} does not exist.")
        return
    try:
        df= pd.read_csv(open(file),encoding="utf-8",dtype=str,)
        fw = open(outfile,'w')
        for index, row in df.iterrows():
            fw.write('<tr><td>'+row['gene']+'</td>'\
                +'<td>'+row['cluster']+'</td>'\
                +'<td>'+row['p_val_adj']+'</td>'\
                +'<td>'+row['p_val']+'</td>'\
                +'<td>'+row['avg_log2FC']+'</td>'\
                +'<td>'+row['pct.1']+'</td>'\
                +'<td>'+row['pct.2']+'</td>'\
            )
        fw.close()
    except Exception as e:
        print(f"An error occurred: {e}")

# atac fragments gz and index
def bgzip_index(fragments, threads):
    bgzip_cmd = [os.path.join(f"{bin_path()}", 'bgzip'), "--force", "--threads", threads, fragments]
    tabix_cmd = [os.path.join(f"{bin_path()}", 'tabix'),'--force' ,'-p', 'bed', f'{fragments}.gz']
    try:
        subprocess.run(bgzip_cmd, check=True)
        subprocess.run(tabix_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during compression or indexing: {e}")
        sys.exit(1)

### generate index for bam
def create_index(threads,bam,outdir):
    try:
        bam_index_cmd = '%s/software/samtools index -@ %s %s'%(__root_dir__,threads,bam)
        logging_call(bam_index_cmd,'count',outdir)
    except Exception as e:
        print('build csi index for bam')
        bam_index_cmd = '%s/software/samtools index -c -@ %s %s'%(__root_dir__,threads,bam)
        logging_call(bam_index_cmd,'count',outdir)

from datetime import timedelta
from functools import wraps
def add_log(func):
    """
    logging start and done.
    """
    logFormatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    module = func.__module__
    name = func.__name__
    logger_name = f"{module}.{name}"
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)

    consoleHandler = logging.StreamHandler(sys.stderr) ##标准输出
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)

    @wraps(func) #在装饰器前加@wraps(func)能帮助保留原有函数的名称和文档字符串
    def wrapper(*args, **kwargs):
        logger.info(f"{name} start...")
        start = time.time()
        result = func(*args, **kwargs) ###执行函数
        end = time.time()
        used = timedelta(seconds=end - start)
        logger.info(f"{name} done. time used: {used}")
        return result 
    return wrapper