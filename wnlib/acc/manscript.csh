#!/bin/csh -f

set nonomatch

foreach manfile (*/*.man side/*/*.man)
  if ( ! -e ../doc/man/`basename $manfile` ) then
    echo "linking $manfile"
    rm -f ../doc/man/`basename $manfile`
    ln -s ../../acc/$manfile ../doc/man/
  endif
end

pushd ../doc/man

foreach manfile (*.man)
  ln -s ../man/$manfile ../mantext/`basename $manfile .man`.txt
end

popd
