"""
Helper functions for cnr pipeline

- fx_name
- is_read1
- is_paired
- list_fq_dir
- list_dir
- list_file
- list_fx
- list_fq_dir
"""

import os
import re
import yaml
import warnings
from pathlib import Path


def fx_name_parts(x: str) -> list:
    """
    Return the sample_name, rep_name, read1/2 
    
    Args:
        x (str): Path to the fastq file
    
    Returns:
        list: return each part of fx name {sample}_{unit}_{read12}{ext}
        
    Examples:
        >>> x = "/path/to/file/RNAseq_mESC_DMSO_2h_rep1_1.fq.gz"
        >>> fx_name_parts(x)
        ['RNAseq_mESC_DMSO_2h', 'rep1', '1', '.fq.gz']
    """
    if isinstance(x, str):
        # name
        name = fx_name(x, fix_pe=True, fix_rep=True)
        # unit|replicates
        p1 = re.compile("_(r(ep)?[0-9]+)_(R)?[12].f(ast)?q", flags=re.IGNORECASE)
        s1 = p1.search(Path(x).name)
        unit = s1.group(1) if bool(s1) else None
        # read12
        p2 = re.compile('(R?[12]).f(ast)?q', flags=re.IGNORECASE)
        s2 = p2.search(Path(x).name)
        read12 = s2.group(1) if bool(s2) else None
        # extension
        ext = Path(x).suffix
        if x.endswith('.gz'):
            xname = Path(x).stem
            ext = Path(xname).suffix + '.gz' # .fq.gz
        out = [name, unit, read12, ext]
    else:
        msg = f"Expect str in fx_name(), got {type(x)}"
        warnings.warn(msg, UserWarning)
        out = None
    return out


def fx_name(
    x: str,
    return_ext: bool = False,
    fix_pe: bool = False,
    fix_rep: bool = False,
) -> str:
    """
    Return the name or file_extension of the fastq file

    Args:
        x (str): Path to the fastq file
        return_ext (bool): Only return the extension of the extension
        fix_pe (bool): remove suffix _R1 or _1 in filename, indicate read1 or read2
        fix_rep (bool): remove replicates suffix _rep1 in filename

    Returns:
        str: filename or extension

    Examples:
        >>> x = "/data/biodata/tmp_1.fq.gz"
        >>> fx_name(x)
        tmp_1
        >>> fx_name(x, return_ext=True)
        .fq
    """
    if isinstance(x, str):
        if x.endswith(".gz"):
            x = Path(x).stem  # remove .gz
        if return_ext:
            out = Path(x).suffix  # ext
        else:
            out = Path(x).stem  # name
            # remove suffix
            if fix_pe:
                out = re.sub("_(R)?[12]$", "", out, flags=re.IGNORECASE)
            if fix_rep:
                out = re.sub("_rep[0-9]+", "", out, flags=re.IGNORECASE)
    else:
        msg = f"Expect str in fx_name(), got {type(x)}"
        warnings.warn(msg, UserWarning)
        out = None
    return out


def is_read1(x: str, read2: bool = False):
    xname = fx_name(x)
    if isinstance(xname, str):
        return xname.endswith("2") if read2 else xname.endswith("1")


def is_paired(x: str, y: str):
    """
    Check if x and y file are paired-end fastq
    Require the filenames should be consistent, except the suffix: 1/2

    Args:
        x (str): filename of fq file
        y (str): filename of fq file

    Returns:
        bool: True or False
    """
    xname = fx_name(x)
    yname = fx_name(y)
    return xname[:-1] == yname[:-1]


def list_dir(
    path,
    full_names: bool = True,
    recursive: bool = False,
    include_dirs: bool = True,
) -> list:
    """
    List all the files within the directory
    see: list.dirs() in R
    see answers on: https://stackoverflow.com/a/3207973

    Args:
        path (str): Path to directory
        full_names (bool): return the full_name of the files
        recrusive (bool): List dir recursively
        include_dirs (bool): inlcude directories

    Returns:
        list: A list of files
    """
    out = []
    if isinstance(path, str):
        if Path(path).is_dir():
            n = 0  # level=0
            for root, d, f in os.walk(path):
                dirs = [str(Path(root) / i) for i in d] if full_names else d
                files = [str(Path(root) / i) for i in f] if full_names else f
                out.extend(files)
                if include_dirs:
                    out.extend(dirs)
                if not recursive:
                    break  # first level
        else:
            msg = f"list_dir() failed, path not a directory: {path}"
            warnings.warn(msg, UserWarning)
    else:
        msg = f"list_dir() failed, expect str, got {type(path)}"
        warnings.warn(msg, UserWarning)
    return out


def list_file(
    path: str = ".",
    pattern: str = "*",
    full_names: bool = True,
    recursive: bool = False,
    include_dirs: bool = False,
) -> list:
    """
    Search files by the pattern, within directory fnmatch.fnmatch()
    see base::list.files() function in R

    Args:
        path (str): Path to a directory
        pattern (str): Pattern of the files, see fnmatch
            \*      matches everything
            \?      matches any single character
            [seq]   matches any character in seq
            [!seq]  matches any char not in seq
            An initial period in FILENAME is not special.
            Both FILENAME and PATTERN are first case-normalized
            if the operating system requires it.
            If you don't want this, use fnmatchcase(FILENAME, PATTERN).
        full_names (bool): Return the full_name of files
        recursive (bool): List the files recursively
        include_dirs (bool): Inlcude directories

    Returns:
        list: A list of files

    Examples:
        >>> list_file("./", "*.fq")
        >>> list_file("/data/biosoft", "*txt", recursive=True)
    """
    files = list_dir(path, full_names, recursive, include_dirs)
    files = [f for f in files if fnmatch.fnmatch(Path(f).name, pattern)]
    return sorted(files)


def list_fx(path: str = ".", fx_type: str = "fq", recursive: bool = False):
    """
    List all fastq/a files within path

    Args:
        path (str): Path to a directory
        fx_type (str): fx type, ['fa', 'fq', 'fx'], default 'fq'
        recursive (bool): recrusively

    Returns:
        list: A list of fastx files
    """
    if fx_type not in ["fa", "fq", "fx"]:
        msg = f"unknown fx_type={fx_type}, expect: [fa, fq, fx]"
        warnings.warn(msg, UserWarning)
        return []
    fx_list = list_dir(
        path, full_names=True, include_dirs=False, recursive=recursive
    )
    if len(fx_list) > 0:
        ext = "a|q" if fx_type.endswith("x") else fx_type[-1]
        p = re.compile(f"\.f(ast)?({ext})(\.gz)?$", flags=re.IGNORECASE)
        fx_list = [i for i in fx_list if p.search(i)]
    return sorted(list(set(fx_list)))


def list_fq_dir(path, paired_only: bool = True, recursive: bool = False) -> list:
    """
    List fastq files within directory and its subdirectory

    Args:
        path (str): Path to the fastq directory
        paired_only (bool): Only list paired end fastq files, _1.fq + _2.fq

    Returns:
        list: A list of fastq files, paired end files were grouped together

    Examples:
        >>> list_fq_dir("./")
        [["f1_1.fq.gz", "f1_2.fq.gz"], ...]
    """
    fq_list = list_fx(path, fx_type="fq", recursive=recursive)
    if paired_only:
        out = []
        fq2_list = fq_list.copy()  # a copy of the list
        for a in fq_list:
            if is_read1(a):
                fq2_list.pop(fq2_list.index(a))  # remove fq1 from fq2_list
                b_list = [i for i in fq2_list if is_paired(a, i)]
                if len(b_list) == 1:
                    b = b_list[0]
                    fq2_list.pop(fq2_list.index(b))  # remove fq2 from fq2_list
                    out.append([a, b])
    else:
        out = fq_list
    return out


def save_to_yaml(data: dict, file_path: str) -> str:
    """
    Save a list of dictionaries to a YAML file.

    Args:
        data (list): List of dictionaries.
        file_path (str): Path to the output YAML file.
    """
    try:
        with open(file_path, "w") as file:
            yaml.dump(data, file, default_flow_style=False)
    except Exception as e:
        msg = f"failed write to YAML file: {file_path}, with error: {e}"
        warnings.warn(msg, UserWarning)
