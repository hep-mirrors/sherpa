sherpa="$1"
setup="$2"

# integrate first, such that all runs below start with the same rng state
"$sherpa" -e 0 "$setup"

# do OTF variation run
"$sherpa" -A OTF "$setup"

run_var() {
    "$sherpa" SCALE_VARIATIONS:None RSF:=$1 FSF:=$2 \
        -A Explicit__ME.MUR=$3__ME.MUF=$4__LHAPDF=93300 \
        "$setup"
}

# do explicit variation runs
run_var 1.00 1.00 1 1
run_var 4.00 1.00 2 1
run_var 4.00 4.00 2 2
run_var 1.00 4.00 1 2
run_var 0.25 1.00 0.5 1
run_var 0.25 0.25 0.5 0.5
run_var 1.00 0.25 1 0.5
