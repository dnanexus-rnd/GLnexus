#!/bin/bash -e

# Prepare for pushing documentation
mkdir html
cd html
git clone https://github.com/dnanexus-rnd/GLnexus.git .
git checkout gh-pages
git config --global user.name $GIT_NAME
git config --global user.email $GIT_EMAIL

# Make doxygen documentation and publish to gh-pages
cd ..
timeNow=`date`; echo  $'\n\n'Documentation Last Updated:$timeNow >> README.md
doxygen Doxyfile
cd html
git add --all .
git commit -m "Auto-updating Doxygen developer documentation"
git push https://$GH_TOKEN@github.com/dnanexus-rnd/GLnexus gh-pages
