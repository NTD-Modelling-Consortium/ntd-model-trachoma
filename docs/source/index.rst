.. NTDMC Trachoma model documentation master file, created by
   sphinx-quickstart on Sun Sep  8 14:31:04 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to NTDMC Trachoma model's documentation!
================================================

On this page:

.. contents::
   :backlinks: none

.. _installation:

Installation
~~~~~~~~~~~~

We recommend that you install the model within a dedicated Python
virtual environment. See `venv — Creation of virtual environments
<https://docs.python.org/3/library/venv.html>`_ for more information
about working with virtual environments.

The model can be installed with the standard `pip` package management
utility, specifying the URL of the project's GitHub repository
followed by the string ``@v<x.y.z>`` where ``<x.y.z>`` is the version
number.  The latest release's version number can be found `here <https://github.com/NTD-Modelling-Consortium/ntd-model-trachoma/tags>`_.

.. code:: shell

   pip install git+https://github.com/NTD-Modelling-Consortium/ntd-model-trachoma.git@v1.0.1

Running the model
~~~~~~~~~~~~~~~~~

Currently, the model is expected to be used as a library through the
unique entry point function `run_single_simulation`.  The calling code
is expected to build the objects required to be passed as parameter of
this model entry point function.

.. code:: python

   from ntd_trachoma import run_single_simulation

   vals, results = run_single_simulation(
       params, vals, timesim, burnin, demog, beta, MDA_times, MDAData,
       vacc_times, VaccData, outputTimes, doSurvey, doIHMEOutput,
       numpy_state,
   )

The arguments expected by the ``run_single_simulation`` function are
described below.

.. _expected-arguments:

Expected arguments
------------------

NumPy arrays are one-dimensional arrays with one element per
individual in the population.

* ``params`` (``dict``) Dictionary made of the following keys:

  * ``N`` (``int``)  The population size.
  * ``av_I_duration`` (``float``) Average duration of the I state
  * ``av_ID_duration`` (``float``) Average duration of the ID state
  * ``av_D_duration`` (``float``) Average duration of the D state
  * ``inf_red`` (``float``) Missing description
  * ``dis_red`` (``float``) Missing description
  * ``min_ID`` (``int``) Minimum duration for the ID state
  * ``min_D`` (``int``) Minimum duration of the D state
  * ``v_1`` (``float``) Missing description
  * ``v_2`` (``float``) Missing description
  * ``phi`` (``float``) Missing description
  * ``epsilon`` (``float``) Missing description
  * ``MDA_Eff`` (``float``) Efficacy of MDA treatment
  * ``rho`` Correlation parameter for systematic non-compliance
  * ``n_inf_sev`` Missing description

* ``vals`` is a dictionary made of the following keys:

  * ``IndI`` (``numpy.ndarray[int]``) Infected status for each
    individual (`0` or `1`)
  * ``IndD`` (``numpy.ndarray[int]``) Diseased status for each
    individual (`0` or `1`)
  * ``T_latent`` (``numpy.ndarray[int]``)  Latent period (weeks)
  * ``T_ID`` (``numpy.ndarray[int]``)  ID period (weeks)
  * ``T_D`` (``numpy.ndarray[int]``)  Diseased period (weeks)
  * ``Ind_latent`` (``numpy.ndarray[float]``) Latent base period
    (weeks)
  * ``Ind_ID_period_base`` (``numpy.ndarray[float]``) ID base period
    (weeks)
  * ``bact_load`` ``numpy.ndarray[float]``  Bacterial load
  * ``Age`` (``numpy.ndarray[int]``)  Age
  * ``N_MDA`` (``int``)  Total number of MDA rounds.
  * ``No_Inf`` ``numpy.ndarray[int]``  Total number of infections.
  * ``True_Prev_Disease_children`` (``list[float]``) Diseased
    prevalence among children aged 9 to 15 years old.  One list
    element per simulated week.
  * ``vaccinated`` (``numpy.ndarray[bool]``) Individuals vaccination
    status
  * ``treatProbability`` (``numpy.ndarray[float]``) MDA treatment
    probability

* ``MDAData`` (``list[list[float, int, int, float, int, int]]``) One
  element per MDA round.  Each element is a list of 6 elements:

  * The MDA round time (``float``), e.g. ``2022.53`` for a round
    sometimes mid-june. Round times are expressed as iteration indices
    within the ``MDA_times`` argument, see below.
  * The population minimum age (``int``)
  * The population maximum age (``int``)
  * The MDA coverage value  (``float``)
  * The MDA campaign index (``int``)
  * The total number of MDA campaigns (``int``)

* ``MDA_times`` (``list[int]``) One element per MDA round.  MDA times
  as a number of iterations, inclusive of burnin time.

* ``VaccData`` (``list[list[float, int, int, float, int, int]]``) One
  element per vaccination round.  Each element is a list of 6
  elements:

  * The vaccination round time (``float``), e.g. ``2022.53`` for a
    round sometimes mid-june. Round times are expressed as iteration
    indices within the ``vaccination_times`` argument, see below.
  * The population minimum age (``int``)
  * The population maximum age (``int``)
  * The vaccination coverage value  (``float``)
  * The vaccination campaign index (``int``)
  * The total number of vaccination campaigns (``int``)

* ``vacc_times`` (``list[int]``) One element per vaccination round.
  Vaccination times as a number of iterations, inclusive of burnin
  time.

* ``outputTimes`` (``list[int]``) Iteration indices at which to record
  output.  See <simulation outputs>.

* ``doSurvey`` (``bool``) Whether or not to carry surveys.

* ``doIHMEOutput`` (``bool``) Whether or not to build IHME outputs.

* ``numpy_state`` NumPy random generator state as returned by ``numpyp.random.get_state()``. See `numpy.random.get_state <https://numpy.org/doc/1.26/reference/random/generated/numpy.random.get_state.html>`_.

.. _outputs:
  
Outputs
-------

The ``run_single_simulation`` function returns a 2-tuple ``(vals,
results)``, which elements are described below.

* ``vals`` (``dict``) The input dictionary ``vals`` with values
  mutated to reflect the final state of the population.  See
  :ref:`expected-arguments`.
* ``results`` (``list[Result]``) Elements are instances of the
  ``ntd-trachoma.Result`` class:

  .. code:: python

     @dataclass
     class Result:
         time: float
         IndI: ndarray
         IndD: ndarray
         Age:ndarray
         NoInf: ndarray
         nMDA:Optional[ndarray] = None
         nMDADoses: Optional[ndarray] = None
         nVacc:Optional[ndarray] = None
         nVaccDoses: Optional[ndarray] = None
         nSurvey: Optional[int] = None
         surveyPass: Optional[int] = None
         elimination: Optional[int] = None
         propMDA: Optional[ndarray] = None    
         propVacc: Optional[ndarray] = None


  The ``results`` list contains one element per output time,
  defined by the ``outputTimes`` list argument to
  ``run_single_simulation``.  See :ref:`expected-arguments`.


