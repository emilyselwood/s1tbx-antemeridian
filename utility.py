import logging
import os
import platform
import re
import sys
import subprocess
import zipfile as zf


def unpack_files(scene, out_path):
    """
    unpack selected files from zip file
    """

    # open archive and iterate on contents
    out_files = []
    with zf.ZipFile(scene) as z:

        files = z.namelist()
        for f in files:
            z.extract(f, out_path)
            out_files.append(os.path.join(out_path, f))

    return out_files


def match_file(files, exp):
    """
    match list of pathname strings with regular expression
    """

    # match single file or throw
    out = match_files(files, exp)
    if len(out) != 1:
        RaiseRuntimeException(
            'Expected single file matching expression {} in dataset: {} - found {} files'.format(exp, scene,
                                                                                                 len(out_files)))

    return out[0]


def match_files(files, exp):
    """
    match list of pathname strings with regular expression
    """

    # for each entry in argument
    out = []
    for f in files:

        # apply regexp on pathname
        m = re.match(exp, f)
        if m:
            # update match list
            out.append(f)

    return out


def find_items(obj, field):
    """
    recursively extract key values from dictionary
    """

    # for all key value pairs
    values = []
    for key, value in obj.items():

        # record value of key match
        if key == field:
            values.append(value)

        # recursive call on nested dict
        elif isinstance(value, dict):
            results = find_items(value, field)
            for result in results:
                values.append(result)

        # loop through contents in array
        elif isinstance(value, list):
            for item in value:

                # recursive call on nested dict
                if isinstance(item, dict):
                    results = find_items(item, field)
                    for result in results:
                        values.append(result)

    return values


def execute(name, arguments):
    """
    create and execute sub-process
    """

    base_env = os.environ.copy()
    if "LD_LIBRARY_PATH" not in base_env and platform.system() != "Windows":
        base_env["LD_LIBRARY_PATH"] = "."
    # create sub-process with argument list
    p = subprocess.Popen([name] + arguments, env=base_env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    code = p.poll()
    logging.debug(out.decode("utf-8"))
    logging.debug(err.decode("utf-8"))
    logging.info(f"return code: {code}")
    return out, err, code
