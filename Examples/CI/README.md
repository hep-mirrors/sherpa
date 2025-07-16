This directory is used to run tests via the GitLab CI. The commands to run
these tests is defined in `.gitlab-ci.yml` (at the project root).

Each subdirectory corresponds to a single test (= Sherpa run), and must at
least contain a `Sherpa.yaml` file.

Some subdirectories might be used by more than one test. The additional tests
then run in a separate work directory on the test runner, but use the same
`Sherpa.yaml` file, only with additional options provided via the command line.
