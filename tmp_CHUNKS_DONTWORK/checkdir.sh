#!/usr/bin/env bash

 ## The target and source can contain spaces as 
 ## long as they are quoted. 
 target="/home/jean/Codes/DustRemission/4Andre/Dust"
 source="/home/jean/Codes/Pynoptic/Interpolation";

 echo target==$target
 echo source==$source
 
 while true; do 
   ## Watch for new files, the grep will return true if a file has
   ## been copied, modified or created.
   echo "inotifywatch -e modify -e create -e moved_to -t 1 "$source" 2>/dev/null | grep total &&" 
   inotifywatch -e modify -e create -e moved_to -t 1 "$source" 2>/dev/null | grep total && 

   ## The -u option to cp causes it to only copy files 
   ## that are newer in $source than in $target. Any files
   ## not present in $target will be copied.
   cp -vu "$source"/* "$target"/
 done
