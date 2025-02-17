"""Collection of basic functions.
"""
from functools import partial
from collections import defaultdict
import re
import os


def add_path(directory, filenames):
    """
    Prepend a directory path to a filename or an arbitrarily nested list
    of filenames.
    Parameters
    ----------
    directory : str
        Base path to prepend.
    filenames : str or list of str
        A single filename or a nested list of filenames.
    Returns
    -------
    str or list of str
        The same structure as filenames, but with directory prepended
        to each filename.
    """
    if isinstance(filenames, str):
        # If a single filename, simply join
        return os.path.join(directory, filenames)
    # Otherwise assume it's a (possibly nested) list and recurse
    return [add_path(directory, fn) for fn in filenames]


def list_fastqs(directory, assay_name=None, use_absolute_paths=False):
    """
    List all valid FASTQ files (ending with .fastq.gz, not 'undetermined')
    from a directory. Optionally filter by an assay name. If requested,
    return absolute paths.
    Parameters
    ----------
    directory : str
        Directory to search.
    assay_name : str or None, default None
        If provided, only files whose first two underscore-delimited segments
        match assay_name are returned (case-insensitive).
    use_absolute_paths : bool, default False
        Whether to return full absolute paths or just filenames.

    Returns
    -------
    list of str
        List of filenames or absolute paths (depending on use_absolute_paths).
    """
    def is_valid_fastq(fname_fastq, assay=None):
        """
        Checks if fname_fastq has .fastq.gz suffix, excludes 'undetermined',
        and optionally checks if the assay is in the filename.
        """
        if not fname_fastq.endswith('.fastq.gz'):
            return False
        if 'undetermined' in fname_fastq.lower():
            return False
        # If assay provided, check if it matches the beginning of the file
        if assay is not None:
            # For example, if assay is "ABC_XYZ", the fastq should start with "ABC_XYZ" (case-insensitive)
            # This check below is consistent with the original code's logic, but can be adjusted as needed.
            normalized_fname = fname_fastq.replace('-', '_').split('_')[:2]
            if assay.upper() not in "_".join(normalized_fname).upper():
                return False
        return True
    # Gather valid FASTQ files
    fastq_files = [
        fn for fn in os.listdir(directory)
        if is_valid_fastq(fn, assay_name)
    ]

    # Convert to absolute paths if requested
    if use_absolute_paths:
        fastq_files = add_path(directory, fastq_files)
    return fastq_files


def list_fastqs_multilane(directory, assay_name=None, use_absolute_paths=False):
    """
    List valid FASTQ files in a directory, optionally filtering by assay name,
    and group them according to lane by removing '_Lxxx' from filenames.
    Parameters
    ----------
    directory : str
        Directory to search.
    assay_name : str or None, default None
        Filter by assay name if provided.
    use_absolute_paths : bool, default False
        Whether to return absolute paths or just filenames.
    Returns
    -------
    dict of {str : list of str}
        A dictionary whose keys are the filename prefixes (with _Lxxx removed),
        and values are lists of matching FASTQ files.
    """
    fastq_files = list_fastqs(directory, assay_name, use_absolute_paths)
    grouped_fastqs = {}
    for fname in fastq_files:
        # Remove pattern like "_L001", "_L002", etc.
        key = re.sub(r'_L\d\d\d', '', fname)
        grouped_fastqs.setdefault(key, []).append(fname)
    return grouped_fastqs


def list_element_is_substring(flag_list, target_str):
    """
    Check if any element of flag_list appears (as a substring) in target_str.

    Case-insensitive check.

    Parameters
    ----------
    flag_list : list of str
        Strings to test for presence in target_str.
    target_str : str
        The string in which we search for substrings.

    Returns
    -------
    bool
        True if at least one element from flag_list appears in target_str,
        otherwise False.
    """
    target_upper = target_str.upper()
    return any(flag.upper() in target_upper for flag in flag_list)


def nested_defaultdict(num_layers, default_factory):
    """
    Create a nested defaultdict with the specified number of layers.
    The final (innermost) layer is initialized with `default_factory`.

    Example
    -------
    d = nested_defaultdict(2, list)
    d['a']['b'].append(1)  # works without intermediate dictionary creation

    Parameters
    ----------
    num_layers : int
        Depth of nesting.
    default_factory : callable
        Factory for the innermost layer (e.g., list, dict, int, etc.).

    Returns
    -------
    defaultdict
        A defaultdict that automatically creates nested dictionaries
        up to `num_layers` deep, with the final layer using `default_factory`.
    """
    if num_layers < 1:
        raise ValueError("Number of layers must be >= 1.")
    if num_layers == 1:
        return defaultdict(default_factory)
    return defaultdict(lambda: nested_defaultdict(num_layers - 1, default_factory))




