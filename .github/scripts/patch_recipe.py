#!/usr/bin/env python
"""Patch the virmake bioconda recipe for a new release.

Stdlib-only, no YAML dependency: the recipe structure is ours, so targeted
regex replacement with loud failure on any missing pattern is simpler and
safer than a round-trip that could reformat the whole file.
"""

import argparse
import re
import sys

URL_TEMPLATE = (
    "https://github.com/Rounge-lab/VirMake/releases/download/"
    "v{version}/virmake-{version}.tar.gz"
)


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("recipe", help="path to meta.yaml")
    parser.add_argument("--version", required=True)
    parser.add_argument("--sha256", required=True)
    args = parser.parse_args()

    if not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", args.version):
        fail(f"version {args.version!r} is not a plain PEP 440 X.Y.Z triplet")
    if not re.fullmatch(r"[0-9a-f]{64}", args.sha256):
        fail("sha256 must be a 64-character hex digest")

    try:
        content = open(args.recipe).read()
    except OSError as e:
        fail(f"cannot read {args.recipe}: {e}")

    replacements = {
        "version": (
            re.compile(r"(?m)^(\s*version:\s*).*$"),
            rf'\g<1>"{args.version}"',
        ),
        "url": (
            re.compile(r"(?m)^(\s*url:\s*).*$"),
            rf"\g<1>{URL_TEMPLATE.format(version=args.version)}",
        ),
        "sha256": (
            re.compile(r"(?m)^(\s*sha256:\s*).*$"),
            rf'\g<1>"{args.sha256}"',
        ),
    }
    for key, (pattern, repl) in replacements.items():
        content, n = pattern.subn(repl, content)
        if n == 0:
            fail(f"no {key}: line found in {args.recipe} - recipe structure changed?")

    with open(args.recipe, "w") as fh:
        fh.write(content)
    print(f"patched {args.recipe}: version={args.version} sha256={args.sha256}")


if __name__ == "__main__":
    main()
