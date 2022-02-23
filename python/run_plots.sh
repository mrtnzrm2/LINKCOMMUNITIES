D=("h" "v")
LC=34
model=zz_model

for k in {1..24}
do
    for d in "${D[@]}"
    do
        echo ${k}
        python3 plot_flatmap_linkcommunities.py true ${model} true tracto2016 0 ${LC} ${d} ${k}
    done
done
