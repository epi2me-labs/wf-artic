"""Scrape versions of conda packages."""

from collections import namedtuple
import os
import subprocess


try:
    import pandas as pd
except ImportError:
    pass


PackageInfo = namedtuple(
    'PackageInfo', ('Name', 'Version', 'Build', 'Channel'))


def scrape_data(as_dataframe=False, include=None, version_dir=None):
    """Return versions of conda packages in base environment."""
    cmd = """
. ~/conda/etc/profile.d/mamba.sh;
micromamba activate;
micromamba list;
    """
    proc = subprocess.run(cmd, shell=True, check=True, capture_output=True)
    versions = dict()
    for line in proc.stdout.splitlines()[3:]:
        items = line.decode().strip().split()
        if len(items) == 3:
            # sometimes channel isn't listed :/
            items.append("")
        if include is None or items[0] in include:
            versions[items[0]] = PackageInfo(*items)
    if version_dir is not None:
        for fname in os.listdir(version_dir):
            print("Reading versions from file:", fname)
            try:
                with open(os.path.join(version_dir, fname), 'r') as fh:
                    for line in fh.readlines():
                        name, version = line.strip().split(',')
                        versions[name] = PackageInfo(name, version, '', '')
            except Exception as e:
                print(e)
                pass
    if as_dataframe:
        versions = pd.DataFrame.from_records(
            list(versions.values()),
            columns=PackageInfo._fields)
    return versions
