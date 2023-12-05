
#Create cmd files for all
for nanop in FAR73740 FAR73782 FAR73806 FAR74611 FAR76967 FAR77015 FAR77731 FAR77735 FAR78499 FAR79983; 
do echo $nanop;  
    for f in $(ls template*.cmd);  
        do cat $f | sed "s/######/$nanop/" > $nanop'_'$f;  
    done; 
done

#create output results 
for nanop in FAR73740 FAR73782 FAR73806 FAR74611 FAR76967 FAR77015 FAR77731 FAR77735 FAR78499 FAR79983; 
do 
echo $nanop;
mkdir $nanop;

done

#create .config
for nanop in FAR73740 FAR73782 FAR73806 FAR74611 FAR76967 FAR77015 FAR77731 FAR77735 FAR78499 FAR79983; 
do echo $nanop;  
    for f in $(ls template*.config);  
        do cat $f | sed "s/######/$nanop/" > $nanop'_'$f;  
    done; 
done