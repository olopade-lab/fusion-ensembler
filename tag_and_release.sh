#!/bin/bash
# to run: ./tag_and_release.sh <VERSION_FOR_TAG>
set -euf -o pipefail

if [ -z "$1" ]
  then
    VERSION="unknown"
  else
    VERSION=$1
fi

POLYFUSE_VERSION=$(python3 -c "import polyfuse; print(polyfuse.__version__)")

if [[ $POLYFUSE_VERSION == "$VERSION" ]]
then
    echo "Version requested matches package version: $VERSION"
else
    echo "[ERROR] Version mismatch. User request: '$VERSION' while package version is: '$POLYFUSE_VERSION'"
    exit 1
fi

echo "Creating tag"
git tag -a "$VERSION" -m "Polyfuse $VERSION"

echo "Pushing tag"
git push origin --tags

rm -f dist/*

echo "======================================================================="
echo "Starting clean builds"
echo "======================================================================="
python3 setup.py sdist
python3 setup.py bdist_wheel

echo "======================================================================="
echo "Done with builds"
echo "======================================================================="


echo "======================================================================="
echo "Push to PyPi. This will require your username and password"
echo "======================================================================="
twine upload dist/*
