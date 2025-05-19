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
from pathlib import Path
from datetime import timedelta
from functools import wraps
from typing import Union, List
from space_sketcher.__init__ import __root_dir__

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

def csv2dict(filename: Path | str, sep=","):
    import pandas as pd
    
    df = pd.read_csv(filename, sep=sep, header=None)
    return dict(zip(df.iloc[:, 0], df.iloc[:, 1]))

# Extracted common path construction logic into a separate function
def get_common_path_part():
    return '/'.join(str(__root_dir__).split('/')[0:-1])

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
    return os.path.join(get_common_path_part(), '.venv/bin')  
    
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


def execute_and_log(
    command: Union[str, List[str]],  # 支持字符串或列表形式的命令
    name: str,                      # 日志标识（如模块名）
    log_dir: str,                   # 日志目录
    shell: bool = True,             # 是否使用shell执行
    log_level: str = "INFO"         # 日志级别（INFO/ERROR）
) -> None:
    """
    执行命令并将结果记录到日志（不返回输出）
    
    Args:
        command: 要执行的命令
        name:    日志名称（用于区分来源）
        log_dir: 日志存储目录
        shell:   是否使用shell模式执行
        log_level: 日志级别（INFO/ERROR）
    """
    # 确保日志目录存在
    os.makedirs(log_dir, exist_ok=True)
    
    # 设置日志（所有命令记录到同一文件）
    log_file = os.path.join(log_dir, "command_execution.log")
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    
    # 避免重复添加handler
    if not logger.handlers:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)
        
    # 判断命令是否需要执行
    if not os.path.exists(os.path.join(log_dir, f".{name}.done")):
        # 记录执行的命令
        cmd_str = command if isinstance(command, str) else " ".join(command)
        logger.info("[COMMAND] %s", cmd_str)
        # 执行命令并记录结果
        try:
            result = subprocess.run(
                command,
                shell=shell,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True
            )
            if result.stdout:
                logger.info("[OUTPUT]\n%s", result.stdout)
            if result.stderr:
                logger.info("[ERROR OUTPUT]\n%s", result.stderr)
                
        except subprocess.CalledProcessError as e:
            logger.error("[FAILED] Exit code: %d\nCommand: %s\nError Output:\n%s", 
                        e.returncode, cmd_str, e.stderr)
            raise
        except Exception as e:
            logger.error("[UNEXPECTED ERROR] %s", str(e))
            raise
    else:
        logger.info(f".{name}.done already exits, skip running {name}")

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

def gunzip(in_file):
    """
    命令行 gunzip 命令
    """
    import subprocess
    import tempfile

    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".txt") as temp_file:
            subprocess.check_call(f"gunzip -c {in_file} > {temp_file.name}", shell=True)
            return temp_file.name
    except Exception as e:
        raise ValueError(f"解压缩文件出错: {e}")

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

