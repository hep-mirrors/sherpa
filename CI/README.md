The CentOS images are stable over many years, so one does not need to update them.
The automated weekly rebuilds of these images are just tiny updates. There is no need to touch them.
The Arch is updated frequently, so each week the image is de-facto new. But because Arch builds don't use all the dependencies,
the updates are smooth and there is no need to touch the Arch configuration.

The Fedora is in a different position: the Fedora versions are bumped each 6 months.
Means that at some point the automated builds with older Fedoras will fail because the will be no required dependency packages.
(The reason why the packages wuld be absent: HEPrpms as of now supports only ~4 last Fedora versions.)
So one will have to bump the versions of Fedora about every 6 months.

To do it:
1) Replace "fedoraX" in the .gitlab-ci.yaml with "fedoraY", Y=X+1. And create a MR.
2) To trigger the image building, the commit should contain "FedoraDockerFileY" string.
   Alternatively one can wait till the next weekly build.
3) One can optionally delete the older images from the registry.
