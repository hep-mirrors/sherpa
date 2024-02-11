import yoda
import subprocess
import sys

# histogram patterns
nominal_patterns = [r"^/MC_XS/XS$"]
var_patterns = [r"^/MC_XS/XS\[.*\]"]
unpatterns = [r"^/MC_XS/XS\[EXTRA.*\]"]

# read all variation histograms from the on-the-fly run results
dfs = yoda.readYODA("OTF.yoda.gz",
                    patterns=var_patterns, unpatterns=unpatterns)

# run yodadiff for each variation
n_diffs_found = 0
for name, df in dfs.items():
    print(name)

    # write on-the-fly variation result into file "a"
    df.setPath("/MC_XS/XS")
    yoda.writeYODA([df], "a.yoda")

    # write explicit variation result into file "b"
    variation_name = name.split('[')[-1][:-1]
    df = yoda.readYODA("Explicit__{}.yoda.gz".format(variation_name),
                       patterns=nominal_patterns,
                       unpatterns=unpatterns,
                       asdict=False)[0]
    yoda.writeYODA([df], "b.yoda")

    # now we can simply use yodadiff
    result = subprocess.run(["yodadiff", "a.yoda", "b.yoda"])
    if result.returncode != 0:
        n_diffs_found += 1
        print("ERROR: diff detected (yodadiff return code:{})".format(
            result.returncode))

# report and exit
if n_diffs_found == 0:
    print("SUCCESS: no diffs detected")
else:
    print("ERROR: diff detected for {} of {} variations".format(
        n_diffs_found, len(dfs)))
    sys.exit(1)
