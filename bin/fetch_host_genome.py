#!/usr/bin/env python3
"""Download the pinned human reference for KneadData and check it byte for byte.

`kneaddata_database --download human_genome bowtie2` fetches whatever the
installed KneadData release points at, and nothing records which file arrived:
KneadData 0.12 gets the T2T-CHM13 index, older releases get hg37, and the
server has since added a v0.2 tarball alongside. A host database that changes
under a run changes which reads count as human, so this pins one tarball by
SHA-256 and the six Bowtie2 index files inside it, and refuses anything else.

    bin/fetch_host_genome.py /path/to/kneaddata_databases
    bin/fetch_host_genome.py --check /path/to/Homo_sapiens_hg39_T2T_Bowtie2_v0.1

The first form downloads, verifies and unpacks into
<dir>/Homo_sapiens_hg39_T2T_Bowtie2_v0.1/, which is what --host_genome takes.
The second verifies a copy that already exists, such as a site's shared one.
Standard library only.
"""

import argparse
import hashlib
import os
import sys
import tarfile
import urllib.request

NAME = "Homo_sapiens_hg39_T2T_Bowtie2_v0.1"
URL = f"https://huttenhower.sph.harvard.edu/kneadData_databases/{NAME}.tar.gz"
TARBALL_SHA256 = "01bda80e810c56e26092a8e7ab2a4a0bb32c517cf77628cbcddedfaf0738a966"
TARBALL_BYTES = 3824531235
FILES = {
    # Bowtie2 index, prefix hg_39: T2T-CHM13 v1.1 (GenBank GCA_009914755.3),
    # chr1-22, X and chrM. No Y; CHM13v2.0 added it.
    "hg_39.1.bt2": "8be88fcedea8d3e8b9ecc361006d26dee4d4779049d86b73a648453ccfcef8d3",
    "hg_39.2.bt2": "b8b3696bae4470182deeef12b6781a8e4fb92c59b77105ca3e2606eb6543d803",
    "hg_39.3.bt2": "53d8aa8c6b3445832d0c9f70e6204e0fc9b4141f1aee1de701baf2278d4e48e9",
    "hg_39.4.bt2": "051c5dafaf1b20c4e83ce3963454ba7ef4e405e0e7b5735099693309251f3fcd",
    "hg_39.rev.1.bt2": "a4fca628d60a2365dadce138f702fecb1c3b3e6fac9c5293d04fc140043febb3",
    "hg_39.rev.2.bt2": "8d6153f66908300f309b062b04b386d2841199c5edb07fcbdc17e59911e71c9b",
}
CHUNK = 8 << 20


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def check(index_dir):
    """Every pinned file present with the pinned hash; return the problems."""
    problems = []
    for name, want in FILES.items():
        path = os.path.join(index_dir, name)
        if not os.path.isfile(path):
            problems.append(f"missing {name}")
        elif sha256(path) != want:
            problems.append(f"{name} does not match the pinned SHA-256")
    return problems


def download(dest):
    """Stream the tarball to dest, hashing as it arrives."""
    h = hashlib.sha256()
    part = dest + ".part"
    with urllib.request.urlopen(URL) as r, open(part, "wb") as out:
        for block in iter(lambda: r.read(CHUNK), b""):
            h.update(block)
            out.write(block)
    if h.hexdigest() != TARBALL_SHA256:
        os.remove(part)
        sys.exit(f"ERROR: {URL} does not match the pinned SHA-256 "
                 f"(got {h.hexdigest()}). The server copy has changed; not unpacking.")
    os.replace(part, dest)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("dir", help="where to unpack (or, with --check, the index directory)")
    p.add_argument("--check", action="store_true", help="only verify an existing index directory")
    p.add_argument("--keep-tarball", action="store_true", help="keep the downloaded .tar.gz")
    args = p.parse_args()

    if args.check:
        index_dir = args.dir
    else:
        os.makedirs(args.dir, exist_ok=True)
        index_dir = os.path.join(args.dir, NAME)
        tarball = os.path.join(args.dir, NAME + ".tar.gz")
        if not os.path.isfile(tarball) or sha256(tarball) != TARBALL_SHA256:
            print(f"Downloading {URL} ({TARBALL_BYTES / 1e9:.1f} GB)", file=sys.stderr)
            download(tarball)
        os.makedirs(index_dir, exist_ok=True)
        with tarfile.open(tarball) as tar:
            members = [m for m in tar.getmembers() if m.isfile() and os.path.basename(m.name) in FILES]
            for m in members:
                m.name = os.path.basename(m.name)
                tar.extract(m, index_dir, **({"filter": "data"} if hasattr(tarfile, "data_filter") else {}))
        if not args.keep_tarball:
            os.remove(tarball)

    problems = check(index_dir)
    if problems:
        sys.exit(f"ERROR: {index_dir} is not {NAME}:\n  " + "\n  ".join(problems))
    print(f"OK {index_dir} is {NAME} (tarball sha256 {TARBALL_SHA256})")
    print(f"Use: --host_genome {os.path.abspath(index_dir)}")


if __name__ == "__main__":
    main()
