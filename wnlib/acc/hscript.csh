#!/bin/csh -f

set nonomatch

foreach hfile (*/*.h)
  if ( ! -e h/`basename $hfile` ) then
    echo "linking $hfile"
    rm -f h/`basename $hfile`
    ln -s ../$hfile h
  endif
end

