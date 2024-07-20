.. _Getting help:

############
Getting help
############


If Sherpa exits abnormally, first check the Sherpa output
for hints on the reason of program abort, and try to figure
out what has gone wrong with the help of the Manual. Note
that Sherpa throwing a ``normal_exit`` exception does not
imply any abnormal program termination! When using AMEGIC++
Sherpa will exit with the message:

.. code-block:: console

   New libraries created. Please compile.

In this case, follow the instructions given in
:ref:`Running Sherpa with AMEGIC++`.

If this does not help, contact the Sherpa team (see the
Sherpa Team section of the website
`sherpa.hepforge.org <http://sherpa.hepforge.org>`_), providing
all information on your setup. Please include

#. A complete tarred and gzipped set of the :option:`Sherpa.yaml` config file
   leading to the crash. Use the status recovery directory
   ``Status__<date of crash>`` produced before the program abort.
#. The command line (including possible parameters) you used to start Sherpa.
#. The installation log file, if available.


.. _Authors:

#######
Authors
#######


Sherpa was written by the Sherpa Team, see
`<https://sherpa-team.gitlab.io/team.html>`_.


.. _Copying:

#######
Copying
#######


Sherpa is free software.
You can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the
Free Software Foundation. You should have received a copy
of the GNU General Public License along with the source
for Sherpa; see the file COPYING. If not, write
to the Free Software Foundation, 59 Temple Place, Suite 330,
Boston, MA  02111-1307, USA.

Sherpa is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

Sherpa was created during the Marie Curie RTN's HEPTOOLS, MCnet and LHCphenonet.
The MCnet Guidelines apply, see the file GUIDELINES and
`<https://www.montecarlonet.org/guidelines>`_.
