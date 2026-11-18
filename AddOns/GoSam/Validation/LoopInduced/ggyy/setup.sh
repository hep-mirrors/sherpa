SHERPA=${1:-"Sherpa"}
$SHERPA Setup.yaml && \
./makelibs &&
echo "Now run Sherpa as normal."
