"""Utility to get versions and decide if we need to update."""
# TODO Move this functionality to epi2melabs, and import in labslauncher

import argparse
import functools
import json
import re
import subprocess
import sys

from cachetools import cached, TTLCache
from cachetools.keys import hashkey
import docker
from github import Github
import requests
import semver


def process_tags(tags_data, prefix):
    """Parse tags and check with semver."""
    tags = list()
    for t in tags_data:
        name = t['name']
        if prefix is not None:
            if name[0] != prefix:
                continue
            try:
                semver.parse(name[1:])
            except ValueError:
                continue
            else:
                tags.append(name[1:])
        else:
            try:
                semver.parse(name)
            except ValueError:
                continue
            else:
                tags.append(name)

    return tags


def proxieskey(*args, proxies=None, **kwargs):
    """Key function to allow hashing function below."""
    key = hashkey(*args, **kwargs)
    if proxies is not None:
        key += tuple(sorted(proxies.items()))
    return key


@cached(cache=TTLCache(maxsize=1, ttl=300), key=proxieskey)
def _get_image_meta(image, proxies=None):
    """Retrieve meta data from docker hub for tags of an image.

    :param image: image name.
    """
    if proxies is None:
        proxies = dict()
    tags = list()
    addr = 'https://hub.docker.com/v2/repositories/{}/tags'.format(image)
    while True:
        response = requests.get(addr, proxies=proxies)
        tags_data = json.loads(response.content.decode())
        tags.extend(tags_data['results'])
        if tags_data['next'] is not None:
            addr = tags_data['next']
        else:
            break
    return tags


def get_image_tags(image, prefix=None, proxies=None):
    """Retrieve tags from dockerhub of an image.

    :param image: image name, organisation/repository.
    :param prefix: prefix by which to filter images.

    :returns: sorted list of tags, newest first, ordered by semver.
        Or the list [None] if an error occurs fetching tag meta information.
    """
    try:
        tags_data = _get_image_meta(image, proxies=proxies)
    except Exception:
        sys.stderr.write("Failed to fetch image information from dockerhub\n")
        return [None]

    tags = process_tags(tags_data, prefix)

    ordered_tags = [
        '{}{}'.format('' if prefix is None else prefix, x) for x in
        sorted(
            tags, reverse=True,
            key=functools.cmp_to_key(semver.compare))]
    return ordered_tags


@functools.lru_cache(5)
def get_image_meta(image, tag, proxies=None):
    """Retrieve meta data from docker hub for a tag.

    :param image: image name.
    :param tag: image tag.

    """
    tags_data = _get_image_meta(image, proxies=proxies)
    for t in tags_data:
        name = t['name']
        if name == tag:
            return t
    raise IndexError("Tag was not found: \"{}\"".format(tag))


def docker_newest_tag(image, prefix=None, client=None, proxies=None):
    """Find the newest available local tag of an image.

    :param image: a docker repo image.
    :param client: a docker client.
    """
    if client is None:
        client = docker.from_env()

    tags = get_image_tags(image, proxies=proxies)
    sorted_tags = sorted(
        tags, reverse=True, key=functools.cmp_to_key(semver.compare)
    )
    all_tags = [
        '{}{}'.format('' if prefix is None else prefix, x) for x in sorted_tags
    ]

    if len(all_tags) > 0:
        return all_tags[0]
    else:
        return None


def github_newest_tag(repository, prefix=None, token=None):
    """Find the newest available tag of a repository.

    :param repository: GitHub org/repository.
    :param token: GitHub personal access token.
    """
    g = Github(token)
    tags_data = [
            {'name': tag.name} for tag in g.get_repo(repository).get_tags()
        ]
    tags = process_tags(tags_data, prefix)
    sorted_t = sorted(
        tags, reverse=True, key=functools.cmp_to_key(semver.compare)
    )
    all_tags = [
        '{}{}'.format('' if prefix is None else prefix, x) for x in sorted_t
    ]

    if len(all_tags) > 0:
        return all_tags[0]
    else:
        return None


def conda_newest_tag(repository, prefix=None):
    """Find the newest available tag of a repository.

    :param repository: conda repository.
    """
    conda_output = subprocess.run(
        ['conda', 'search', '-c', 'bioconda', '-f', repository],
        stdout=subprocess.PIPE
    )
    tags_data = list()

    for line in conda_output.stdout.decode('utf-8').split("\n"):
        if line.startswith(repository):
            package, name, repo_hash, repo = re.split(" +", line.rstrip())
            tags_data.append({'name': name})

    tags = process_tags(tags_data, prefix)
    sorted_t = sorted(
        tags, reverse=True, key=functools.cmp_to_key(semver.compare))
    all_tags = [
        '{}{}'.format('' if prefix is None else prefix, x) for x in sorted_t
    ]

    if len(all_tags) > 0:
        return all_tags[0]
    else:
        return None


def main():
    """Get Versions of aux software."""
    parser = argparse.ArgumentParser(
        description='''Get latest tags''')

    parser.add_argument(
        '-t', '--token', required=True, dest="token",
        help="GitHub Access Token"
    )

    parser.add_argument(
        '-d', '--docker_registry', required=True,
        dest="docker_registry", help="Docker registry"
    )

    parser.add_argument(
        '-g', '--github_repository', required=True,
        dest="github_repository", help="GitHub repository"
    )

    parser.add_argument(
        '-o', '--tool', required=True,
        dest="tool", help="Tool"
    )

    parser.add_argument(
        '-p', '--prefix', required=False,
        dest="prefix", default=None
    )

    args = parser.parse_args()

    docker_newest = docker_newest_tag(args.docker_registry)

    sys.stdout.write(f"DOCKER:{docker_newest}\n")

    github_newest = github_newest_tag(
            args.github_repository,
            prefix=args.prefix,
            token=args.token
        )

    sys.stdout.write(f"GITHUB:{github_newest}\n")

    conda_newest = conda_newest_tag(args.tool)

    sys.stdout.write(f"CONDA:{conda_newest}\n")

    try:
        if semver.compare(
                docker_newest,
                conda_newest
                ) == -1:
            sys.stdout.write(f"ACTION:{conda_newest}\n")
        else:
            sys.stdout.write("ACTION:NO_UPDATE\n")
    except ValueError:
        sys.stdout.write(f"ACTION:{conda_newest}\n")


if __name__ == "__main__":
    main()
