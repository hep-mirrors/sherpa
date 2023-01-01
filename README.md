# The Sherpa event generator framework for high-energy particle collisions

Sherpa is an event-generation framework for high-energy
particle collisions.

For more information, visit [our homepage](http://sherpa-team.gitlab.io)
or our [project site on GitLab](https://gitlab.com/sherpa-team/sherpa).
You can file a bug report on our [issue tracker](https://gitlab.com/sherpa-team/sherpa/issues)
or send us an [e-mail](sherpa@projects.hepforge.org).

To install use the commands
```
cmake -S <sherpadir> -B <builddir> [+ optional configuration options, see cmake -LAH]
cmake --build <builddir> [other build options, e.g. -j 8]
cmake --install <builddir>
```
Required dependencies are LHAPDF and libzip, they can be automaticall installed by adding
the options
```
-DSHERPA_ENABLE_INSTALL_LHAPDF=ON -DSHERPA_ENABLE_INSTALL_LIBZIP=ON
```
to the first `cmake` command. There are some optional features which can be compiled 
into Sherpa, please read the info page (`info Manual/Sherpa.info`) 
or run `cmake -LAH` to find out more about the available options.
