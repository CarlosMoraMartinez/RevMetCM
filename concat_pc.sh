  #!/bin/bash

  pcs2=$(tr -d "\[\]," <<<$pcs )
  for pc in $pcs2
  do 
    cat $pc >> all_$ont_id'.pc'
  done